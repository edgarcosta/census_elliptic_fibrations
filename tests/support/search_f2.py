import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import os
import shlex
import shutil
import tempfile

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = json.loads((ROOT / 'tests/fixtures/f2-recovery/case.json').read_text())


def call(*args, ok=True):
    result = subprocess.run([sys.executable, str(ROOT / 'production/search_f2.py'),
                             *map(str, args)], capture_output=True, text=True)
    assert (result.returncode == 0) == ok, result.stdout + result.stderr
    return result


def setup(root, candidates=None):
    source = FIXTURE['source']
    row = dict(q=2, key=source['key'], codes=source['codes'], L=source['L'],
               rank=1, torsion=[], sha_reg='224', candidates=candidates or [[16, '14']],
               stratum=[1, False, []], source=dict(input_ordinal=0,
               stage1_receipt_sha256='a'*64, receipt='fixture', compressed_sha256='b'*64,
               keys_sha256=None))
    inputs = root / 'queue.jsonl.gz'
    with gzip.open(inputs, 'wt') as stream:
        stream.write(json.dumps(row)+'\n')
    manifest = root / 'manifest.json'
    manifest.write_text(json.dumps(dict(input_sha256=hashlib.sha256(inputs.read_bytes()).hexdigest(),
        input_count=1, source_count=1, existing_certified=0, source_identity='fixture',
        counts=dict(open_splits=1 if candidates and len(candidates)>1 else 0),
        source_manifest=json.loads((ROOT / 'data/mw/f2-recovery/source-manifest.json').read_text()))))
    run = root / 'run'
    call('prepare', '--inputs', inputs, '--manifest', manifest,
         '--cohort', 'pilot', '--sample-size', 1, '--output', run)
    return run


def retained_batch(run):
    import importlib.util
    spec = importlib.util.spec_from_file_location('recovery_test', ROOT / 'tests/support/recover_f2.py')
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    directory = run / 'batches/00000'
    directory.mkdir(parents=True)
    stage = helper.setup(directory)
    stage.rename(directory / 'b5')
    stage = directory / 'b5'
    prepared = json.loads((stage / 'prepared.json').read_text())
    manifest = json.loads((run / 'prepared.json').read_text())['manifest']['source_manifest']
    prepared['source_manifest'] = {k: manifest[k] for k in ('calibration_head', 'all_kernel_dependencies_compared')}
    (stage / 'prepared.json').write_text(json.dumps(prepared))
    producer = json.loads((stage / 'replay/producer.json').read_text())
    producer['prepared_sha256'] = helper.digest(stage / 'prepared.json')
    (stage / 'replay/producer.json').write_text(json.dumps(producer))
    receipt = json.loads((stage / 'replay/0000.json').read_text())
    receipt['producer_sha256'] = helper.digest(stage / 'replay/producer.json')
    (stage / 'replay/0000.json').write_text(json.dumps(receipt))
    helper.audit(stage)
    (directory / 'result.json').write_text(json.dumps(dict(
        input_sha256=helper.digest(run / 'requests/00000.jsonl'),
        prepared_sha256=helper.digest(run / 'prepared.json'), outcomes=[dict(key=4999168,
        reason='independent height and division audit', attempts=[dict(bound=5, status='certified')])])) )


def test_unique_split_and_resume(root):
    run = setup(root)
    kernel = os.environ.get('SEARCH_F2_KERNEL')
    if kernel:
        call('run', '--output', run, '--kernel', kernel, '--jobs', 1)
    else:
        retained_batch(run)
    call('verify', '--output', run)
    records = [json.loads(x) for x in (run / 'refinements.jsonl').read_text().splitlines()]
    assert records[0]['claims']['mw_basis']['status'] == 'exact'
    assert records[0]['claims']['mw_basis']['value']['gram']['upper'] == ['14']
    assert set(records[0]['claims']) == {'mw_basis'}
    assert json.loads((run / 'verification.json').read_text())['new_bases'] == 1
    before = {p: p.read_bytes() for p in (run / 'batches').rglob('*') if p.is_file()}
    if kernel:
        call('run', '--output', run, '--kernel', kernel, '--jobs', 1)
    call('verify', '--output', run)
    assert all(p.read_bytes() == data for p, data in before.items())
    result = next((run / 'batches').glob('*/result.json'))
    data = json.loads(result.read_text())
    assert [x['bound'] for x in data['outcomes'][0]['attempts']] == [5]
    raw = next((run / 'batches').glob('*/b5/replay/0000.jsonl'))
    native = json.loads(raw.read_text())
    native['key'] += 1
    raw.write_text(json.dumps(native)+'\n')
    assert 'digest' in call('verify', '--output', run, ok=False).stderr


def test_candidate_conflict(root):
    run = setup(root, [[16, '14'], [4, '56']])
    retained_batch(run)
    call('verify', '--output', run)
    record = json.loads((run / 'refinements.jsonl').read_text())
    assert record['claims']['sha_order']['value'] == '16'
    assert record['claims']['regulator']['value'] == '14'
    assert json.loads((run / 'verification.json').read_text())['new_splits'] == 1
    rest = root / 'remaining'
    call('prepare', '--inputs', root / 'queue.jsonl.gz', '--manifest', root / 'manifest.json',
         '--cohort', 'remaining', '--sample-size', 1, '--output', rest)
    call('verify', '--output', rest)
    combined = root / 'combined'
    call('summarize', '--runs', run, rest, '--output', combined)
    summary = json.loads((combined / 'summary.json').read_text())
    assert summary['selected'] == summary['new_bases'] == summary['new_splits'] == 1
    assert summary['remaining_splits'] == 0
    assert 'duplicate' in call('summarize', '--runs', run, run, '--output', root / 'duplicate', ok=False).stderr




def test_failed_resume(root):
    kernel = os.environ.get('SEARCH_F2_KERNEL')
    if not kernel:
        return
    run = setup(root)
    marker = root / 'fail-audit'
    marker.touch()
    wrapper = root / 'magma-gate'
    wrapper.write_text('#!/bin/sh\nif [ "$1" = "-bn" ] && [ -e ' +
        shlex.quote(str(marker)) + ' ]; then exit 7; fi\nexec ' +
        shlex.quote(shutil.which('magma')) + ' "$@"\n')
    wrapper.chmod(0o755)
    args = ('run', '--output', run, '--kernel', kernel, '--magma', wrapper, '--jobs', 1)
    call(*args, ok=False)
    marker.unlink()
    call(*args)
    call('verify', '--output', run)
    result = json.loads((run / 'verification.json').read_text())
    assert result['new_bases'] == 1
    assert list((run / 'batches/00000/failed').glob('*-result.json'))
    assert len((run / 'jobs.tsv').read_text().splitlines()) == 3


def test_input_and_status_checks(root):
    run = setup(root)
    inputs = root / 'queue.jsonl.gz'
    manifest = root / 'manifest.json'
    row = json.loads(gzip.decompress(inputs.read_bytes()))
    for field, bad, message in [('key', row['key']+1, 'key'),
                                ('candidates', [], 'empty'),
                                ('rank', 0, 'rank')]:
        altered = dict(row, **{field: bad})
        inputs.write_bytes(gzip.compress((json.dumps(altered)+'\n').encode()))
        data = json.loads(manifest.read_text())
        data['input_sha256'] = hashlib.sha256(inputs.read_bytes()).hexdigest()
        manifest.write_text(json.dumps(data))
        result = call('prepare', '--inputs', inputs, '--manifest', manifest,
                      '--cohort', 'pilot', '--sample-size', 1, '--output', root / ('bad-'+field), ok=False)
        assert message in result.stderr, result.stderr


if __name__ == '__main__':
    for test in (test_unique_split_and_resume, test_candidate_conflict, test_input_and_status_checks, test_failed_resume):
        with tempfile.TemporaryDirectory(dir=ROOT / 'artifacts') as directory:
            test(Path(directory))
