import copy
import gzip
import hashlib
import json
import shlex
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = json.loads((ROOT / 'tests/fixtures/f2-recovery/case.json').read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def call(*args, ok=True):
    result = subprocess.run([sys.executable, str(ROOT / 'production/recover_f2.py'),
                             *map(str, args)], capture_output=True, text=True)
    assert (result.returncode == 0) == ok, result.stdout + result.stderr
    return result


def setup(root):
    inputs = root / 'inputs.jsonl'
    inputs.write_text(json.dumps(FIXTURE['request']) + '\n')
    shard = root / 'shard.jsonl.gz'
    with gzip.open(shard, 'wt') as stream:
        stream.write(json.dumps(FIXTURE['source']) + '\n')
    receipt = root / 'shard.json'
    receipt.write_text(json.dumps(dict(q=2, count=1, source='fixture-source',
                                      sha256=digest(shard))))
    summary = root / 'summary.json'
    summary.write_text(json.dumps(dict(complete=True, count=1, generators_proved=0,
                                      source='fixture-source')))
    manifest = root / 'sources.json'
    kernel = json.loads((ROOT / 'data/mw/f2-recovery/source-manifest.json').read_text())
    manifest.write_text(json.dumps(dict(input_sha256=digest(inputs), input_count=1,
        calibration_head=kernel['calibration_head'],
        all_kernel_dependencies_compared=kernel['all_kernel_dependencies_compared'],
        summary='summary.json', shards=[dict(receipt='shard.json', count=1,
        receipt_sha256=digest(receipt), compressed_sha256=digest(shard))])))
    run = root / 'run'
    call('prepare', '--inputs', inputs, '--manifest', manifest,
         '--evidence-root', root, '--output', run)
    (run / 'replay').mkdir()
    raw = run / 'replay/0000.jsonl'
    raw.write_text(json.dumps(FIXTURE['native']) + '\n')
    (run / 'replay/producer.json').write_text(json.dumps(dict(
        method='retained-native-bundle-v1', executable_sha256=FIXTURE['executable_sha256'],
        kernel_revision=kernel['calibration_head'],
        kernel_dependencies=kernel['all_kernel_dependencies_compared'],
        prepared_sha256=digest(run / 'prepared.json'))))
    (run / 'replay/0000.json').write_text(json.dumps(dict(
        request_sha256=digest(run / 'requests/0000.jsonl'),
        output_sha256=digest(raw), producer_sha256=digest(run / 'replay/producer.json'))))
    return run


def audit(run):
    call('audit', '--output', run, '--jobs', 1)


def imported(run):
    call('import', '--output', run)
    return [json.loads(line) for line in (run / 'refinements.jsonl').read_text().splitlines()]


# A deeper certificate must survive import while existing exact claims stay intact.
def test_deeper_certificate_and_resume(root):
    run = setup(root)
    audit(run)
    records = imported(run)
    assert len(records) == 1
    record = records[0]
    assert record['key'] == '4999168'
    assert record['claims']['mw_basis']['status'] == 'exact'
    assert record['claims']['mw_basis']['value']['gram'] == {'dimension': 1, 'upper': ['14']}
    assert set(record['claims']) == {'mw_basis'}
    before = {p: p.read_bytes() for p in (run / 'audit').rglob('*') if p.is_file()}
    audit(run)
    assert all(p.read_bytes() == value for p, value in before.items())
    assert imported(run) == records
    summary = json.loads((run / 'import.json').read_text())
    assert summary['certified_total'] == 1 and summary['new_splits_credited'] == 0


# Failed work cannot promote a basis, and a wrong-key replay is an error.
def test_incomplete_and_wrong_key(root):
    run = setup(root)
    raw = run / 'replay/0000.jsonl'
    original = raw.read_bytes()
    receipt = run / 'replay/0000.json'
    receipt_data = json.loads(receipt.read_text())
    receipt.unlink()
    raw.unlink()
    assert imported(run)[0]['claims']['mw_basis']['status'] == 'unresolved'
    bad = copy.deepcopy(FIXTURE['native'])
    bad['key'] += 1
    raw.write_text(json.dumps(bad) + '\n')
    receipt_data['output_sha256'] = digest(raw)
    receipt.write_text(json.dumps(receipt_data))
    result = call('import', '--output', run, ok=False)
    assert 'key' in result.stderr.lower(), result.stderr
    raw.write_bytes(original)
    receipt_data['output_sha256'] = digest(raw)
    receipt.write_text(json.dumps(receipt_data))
    audit(run)
    raw.write_bytes(original + original)
    receipt_data['output_sha256'] = digest(raw)
    receipt.write_text(json.dumps(receipt_data))
    result = call('import', '--output', run, ok=False)
    assert 'count' in result.stderr.lower() or 'duplicate' in result.stderr.lower()


# An incomplete retry must keep a previously imported exact basis.
def test_exact_import_cannot_downgrade(root):
    run = setup(root)
    audit(run)
    imported(run)
    before = (run / 'refinements.jsonl').read_bytes()
    (run / 'audit/0000.json').unlink()
    result = call('import', '--output', run, ok=False)
    assert 'exact' in result.stderr.lower(), result.stderr
    assert (run / 'refinements.jsonl').read_bytes() == before


# A capped search is evidence of an attempt, not a basis certificate.
def test_capped_search(root):
    run = setup(root)
    raw = run / 'replay/0000.jsonl'
    result = copy.deepcopy(FIXTURE['native'])
    result['certified'] = False
    raw.write_text(json.dumps(result) + '\n')
    receipt = run / 'replay/0000.json'
    data = json.loads(receipt.read_text())
    data['output_sha256'] = digest(raw)
    receipt.write_text(json.dumps(data))
    assert imported(run)[0]['claims']['mw_basis']['status'] == 'unresolved'
    summary = json.loads((run / 'import.json').read_text())
    assert summary['recovered'] == 0 and summary['unresolved'] == 1


# Torsion coordinates are part of the published basis evidence.
def test_wrong_torsion_point(root):
    run = setup(root)
    raw = run / 'replay/0000.jsonl'
    data = copy.deepcopy(FIXTURE['native'])
    data['torsion_points'].append([0, 1, 0, 1])
    raw.write_text(json.dumps(data) + '\n')
    receipt = run / 'replay/0000.json'
    data = json.loads(receipt.read_text())
    data['output_sha256'] = digest(raw)
    receipt.write_text(json.dumps(data))
    call('audit', '--output', run, '--jobs', 1, ok=False)
    assert imported(run)[0]['claims']['mw_basis']['status'] == 'unresolved'


# GNU parallel retries a failed process without crediting its incomplete output.
def test_failed_audit_resume(root):
    run = setup(root)
    wrapper = root / 'magma-once'
    marker = root / 'failed-once'
    actual = shutil.which('magma')
    wrapper.write_text('#!/bin/sh\n'
        + 'if [ "$1" = "-bn" ] && [ ! -e ' + shlex.quote(str(marker)) + ' ]; then\n'
        + 'touch ' + shlex.quote(str(marker)) + '; exit 7; fi\n'
        + 'exec ' + shlex.quote(actual) + ' "$@"\n')
    wrapper.chmod(0o755)
    args = ('audit', '--output', run, '--magma', wrapper, '--jobs', 1)
    call(*args, ok=False)
    assert marker.exists() and not (run / 'audit/0000.json').exists()
    assert imported(run)[0]['claims']['mw_basis']['status'] == 'unresolved'
    call(*args)
    assert imported(run)[0]['claims']['mw_basis']['status'] == 'exact'
    log = (run / 'audit/jobs.tsv').read_text().splitlines()
    assert len(log) == 3
    assert int(log[1].split()[6]) != 0 and int(log[2].split()[6]) == 0


# Producer identity must agree with the prepared kernel, even before auditing.
def test_wrong_kernel_identity(root):
    run = setup(root)
    producer = run / 'replay/producer.json'
    data = json.loads(producer.read_text())
    data['kernel_revision'] = '0' * 40
    producer.write_text(json.dumps(data))
    receipt = run / 'replay/0000.json'
    data = json.loads(receipt.read_text())
    data['producer_sha256'] = digest(producer)
    receipt.write_text(json.dumps(data))
    result = call('import', '--output', run, ok=False)
    assert 'kernel' in result.stderr.lower(), result.stderr


if __name__ == '__main__':
    (ROOT / 'artifacts').mkdir(exist_ok=True)
    for test in (test_deeper_certificate_and_resume, test_incomplete_and_wrong_key,
                 test_exact_import_cannot_downgrade, test_capped_search, test_wrong_torsion_point, test_failed_audit_resume,
                 test_wrong_kernel_identity):
        with tempfile.TemporaryDirectory(dir=ROOT / 'artifacts') as directory:
            test(Path(directory))
