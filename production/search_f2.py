#!/usr/bin/env python3
"""Search accessible open F_2 bases through bound 7 and audit each certificate."""

import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import random
import shlex
import shutil
import subprocess
import sys
import tarfile
from types import SimpleNamespace

import recover_f2 as recovery
from recover_f2 import (AUDITOR, ROOT, certified_value, digest, encoded, immutable,
                        lines, locked, matrix_determinant, read, require, save)

DATA = ROOT / 'data/mw/f2-accessible'
RECOVERY = ROOT / 'data/mw/f2-recovery'
SEED = 20261119
BOUNDS = (5, 6, 7)


def valid_row(row):
    a1, a2, a3, a4, a6 = row['codes']
    require(row['q'] == 2 and a2 == 0 and row['rank'] in (1, 2, 3), 'queue field/rank')
    require(all(isinstance(c, int) and 0 <= c < 2**d for c, d in zip(row['codes'], (3, 5, 7, 9, 13))), 'coefficient range')
    require(row['key'] == (a1 << 29) + (a3 << 22) + (a4 << 13) + a6, 'queue key')
    require(row['candidates'], 'empty candidate intersection')
    require(len({s for s, _ in row['candidates']}) == len(row['candidates']), 'duplicate candidates')
    require(all(s > 0 and Fraction(r) > 0 and s*Fraction(r) == Fraction(row['sha_reg'])
                for s, r in row['candidates']), 'candidate product')
    require(row['stratum'] == [row['rank'], bool(a1), row['torsion']], 'queue stratum')


def bounds(row):
    return [b for b in BOUNDS if row['rank'] != 1 or
            min(Fraction(r) for _, r in row['candidates']) <= 4 + 2*b]


def export(args):
    manifest = read(RECOVERY / 'source-manifest.json')
    summary = read(args.evidence_root / manifest['summary'])
    recovered_summary = read(RECOVERY / 'recovery-summary.json')
    archive = RECOVERY / recovered_summary['archive']
    require(digest(archive) == recovered_summary['archive_sha256'], 'recovery archive digest')
    with tarfile.open(archive) as stream:
        recovered = [json.loads(x) for x in stream.extractfile('refinements.jsonl')]
    require(all(r['claims']['mw_basis']['status'] == 'exact' for r in recovered), 'recovery status')
    recovered_keys = {int(r['key']) for r in recovered}
    require(len(recovered_keys) == len(recovered) == recovered_summary['import_result']['recovered'], 'recovery key count')
    require(recovered_keys == {r['key'] for r in lines(RECOVERY / 'replay-candidates.jsonl')}, 'recovery keys')
    with args.refined.open() as stream:
        refinements = list(csv.DictReader(stream))
    refined = {int(r['key']): r for r in refinements}
    require(len(refined) == len(refinements), 'duplicate refined key')
    rows, seen, matched, recovered_seen = [], set(), set(), set()
    counts = Counter()
    for entry in manifest['shards']:
        path = args.evidence_root / entry['receipt']
        require(digest(path) == entry['receipt_sha256'], 'source receipt digest')
        receipt = read(path)
        require(receipt['q'] == 2 and receipt['source'] == summary['source'], 'source identity')
        shard = path.with_suffix('.jsonl.gz')
        require(digest(shard) == receipt['sha256'] == entry['compressed_sha256'], 'source shard digest')
        n = 0
        with gzip.open(shard, 'rt') as stream:
            for ordinal, line in enumerate(stream):
                source = json.loads(line)
                key = source['key']
                require(source['q'] == 2 and key not in seen, 'source field/duplicate key')
                seen.add(key)
                n += 1
                counts['source_count'] += 1
                pairs = source['split']['candidates']
                if key in refined:
                    change = refined[key]
                    require(json.loads(change['old_candidate_list']) == [s for s, _ in pairs], 'refined source candidates')
                    keep = json.loads(change['new_candidate_list'])
                    require(keep and set(keep) <= {s for s, _ in pairs}, 'refined candidates')
                    pairs = [[s, r] for s, r in pairs if s in keep]
                    matched.add(key)
                counts['open_splits'] += len(pairs) > 1
                if source['generators']['certified']:
                    counts['source_certified'] += 1
                    continue
                if key in recovered_keys:
                    counts['recovered_excluded'] += 1
                    recovered_seen.add(key)
                    continue
                row = dict(q=2, key=key, codes=source['codes'], L=source['L'],
                    rank=source['derived']['rank'], torsion=source['generators']['torsion'],
                    sha_reg=source['generators']['sha_reg'], candidates=pairs,
                    stratum=[source['derived']['rank'], bool(source['codes'][0]), source['generators']['torsion']],
                    source=dict(input_ordinal=ordinal, receipt=entry['receipt'],
                        stage1_receipt_sha256=entry['receipt_sha256'], compressed_sha256=entry['compressed_sha256'],
                        keys_sha256=receipt.get('keys_sha256')))
                valid_row(row)
                if bounds(row):
                    rows.append(row)
                    counts[f'eligible_rank_{row["rank"]}'] += 1
                else:
                    counts['height_excluded'] += 1
        require(n == receipt['count'] == entry['count'], 'source count')
    require(summary['complete'] and counts['source_count'] == summary['count'] and
            counts['source_certified'] == summary['generators_proved'], 'source summary')
    require(matched == set(refined) and recovered_seen == recovered_keys, 'missing overlay key')
    rows.sort(key=lambda r: r['key'])
    args.output.mkdir(parents=True, exist_ok=True)
    path = args.output / 'queue.jsonl.gz'
    immutable(path, gzip.compress(b''.join(encoded(r) for r in rows), mtime=0))
    result = dict(input_sha256=digest(path), input_count=len(rows), source_manifest=manifest,
        source_manifest_sha256=digest(RECOVERY / 'source-manifest.json'),
        refined_sha256=digest(args.refined), refined_count=len(refined),
        recovery_archive_sha256=digest(archive), source_identity=summary['source'],
        source_summary_sha256=digest(args.evidence_root / manifest['summary']),
        source_count=counts['source_count'], existing_certified=counts['source_certified'] + len(recovered_keys),
        counts=dict(counts), seed=SEED, pilot_size=4096)
    immutable(args.output / 'manifest.json', encoded(result))
    print(json.dumps({k: v for k, v in result.items() if k != 'source_manifest'}, sort_keys=True))


def sample(rows, size):
    groups = defaultdict(list)
    for row in rows:
        groups[json.dumps(row['stratum'])].append(row['key'])
    require(size >= len(groups), 'sample too small for joint strata')
    rng = random.Random(SEED)
    chosen = {rng.choice(sorted(groups[g])) for g in sorted(groups)}
    chosen.update(rng.sample(sorted({r['key'] for r in rows} - chosen), min(size, len(rows)) - len(chosen)))
    return chosen


def prepare(args):
    manifest = read(args.manifest)
    require(digest(args.inputs) == manifest['input_sha256'], 'queue digest')
    with gzip.open(args.inputs, 'rt') as stream:
        rows = [json.loads(line) for line in stream]
    require(len(rows) == manifest['input_count'] == len({r['key'] for r in rows}), 'queue key count')
    for row in rows:
        valid_row(row)
        require(bounds(row), 'inaccessible queue row')
    pilot = sample(rows, min(args.sample_size, len(rows)))
    selected = [r for r in rows if (r['key'] in pilot) == (args.cohort != 'remaining')]
    files = {}
    for i in range(0, len(selected), 8):
        name = f'requests/{i//8:05d}.jsonl'
        immutable(args.output / name, b''.join(encoded(r) for r in selected[i:i+8]))
        files[name] = digest(args.output / name)
    require(args.sample_size > 0, 'sample size must be positive')
    identity = dict(method='accessible-f2-bases-v1', script_sha256=digest(Path(__file__)),
        recovery_sha256=digest(ROOT / 'production/recover_f2.py'), audit_sha256=digest(AUDITOR))
    immutable(args.output / 'prepared.json', encoded(dict(files=files, manifest=manifest,
        manifest_sha256=digest(args.manifest), cohort=args.cohort, sample_size=args.sample_size,
        seed=SEED, selected=len(selected), identity=identity)))
    print(json.dumps(dict(selected=len(selected), batches=len(files), cohort=args.cohort)))


def prepared(output, batch=None):
    metadata = read(output / 'prepared.json')
    require(metadata['identity'] == dict(method='accessible-f2-bases-v1', script_sha256=digest(Path(__file__)),
        recovery_sha256=digest(ROOT / 'production/recover_f2.py'), audit_sha256=digest(AUDITOR)), 'changed search producer')
    names = metadata['files'] if batch is None else [f'requests/{batch}.jsonl']
    for name in names:
        require(name in metadata['files'] and digest(output / name) == metadata['files'][name], 'prepared input digest')
    return metadata


def kernel_manifest(metadata):
    manifest = metadata['manifest']['source_manifest']
    return {k: manifest[k] for k in ('calibration_head', 'all_kernel_dependencies_compared')}


def make_stage(path, rows, b, metadata):
    requests = [dict(q=2, key=r['key'], codes=r['codes'], L=r['L'], b=b) for r in rows]
    immutable(path / 'requests/0000.jsonl', b''.join(encoded(r) for r in requests))
    immutable(path / 'prepared.json', encoded(dict(files={'requests/0000.jsonl': digest(path / 'requests/0000.jsonl')},
        source_manifest=kernel_manifest(metadata))))
    return requests


def value(row, request, native, checked):
    determinant = matrix_determinant(native['gram'])
    pairs = [[s, r] for s, r in row['candidates'] if Fraction(r) == determinant]
    require(len(pairs) == 1, 'certified determinant outside refined candidates')
    sha, regulator = pairs[0]
    source = dict(derived=dict(rank=row['rank']), generators=dict(sha_reg=row['sha_reg'], torsion=row['torsion']),
        split=dict(candidates=pairs, sha=sha, regulator=regulator))
    return certified_value(source, request, native, checked)


def worker(args, metadata):
    batch = args.batch
    rows = lines(args.output / f'requests/{batch}.jsonl')
    directory = args.output / f'batches/{batch}'
    directory.mkdir(parents=True, exist_ok=True)
    result = directory / 'result.json'
    if result.exists():
        _, old, _ = check_batch(args.output, batch, metadata)
        if not any('retry pending' in o['reason'] for o in old):
            return
        failed = directory / 'failed'
        failed.mkdir(exist_ok=True)
        number = len(list(failed.iterdir()))
        result.rename(failed / f'{number}-result.json')
        for b in BOUNDS:
            stage = directory / f'b{b}'
            if stage.exists():
                stage.rename(failed / f'{number}-b{b}')
    outcomes = {r['key']: dict(key=r['key'], reason='bounds exhausted without certified basis', attempts=[]) for r in rows}
    pending = list(rows)
    for b in BOUNDS:
        eligible = [r for r in pending if b in bounds(r)]
        if not eligible:
            continue
        stage = directory / f'b{b}'
        requests = make_stage(stage, eligible, b, metadata)
        options = SimpleNamespace(output=stage, batch='0000', jobs=1, kernel=args.kernel, magma=args.magma)
        try:
            recovery.parallel(options, 'replay')
        except subprocess.SubprocessError as error:
            for row in eligible:
                outcomes[row['key']]['attempts'].append(dict(bound=b, status='failed', reason=str(error)))
            continue
        _, natives = recovery.native_batch(stage, '0000')
        for row, request, native in zip(eligible, requests, natives):
            if native['certified']:
                value(row, request, native, None)
        try:
            recovery.parallel(options, 'audit')
        except subprocess.SubprocessError as error:
            for row in eligible:
                outcomes[row['key']]['attempts'].append(dict(bound=b, status='audit_failed', reason=str(error)))
            continue
        audit = recovery.audit_result(stage, '0000')
        for row, request, native in zip(eligible, requests, natives):
            key = row['key']
            attempt = dict(bound=b, status='certified' if native['certified'] else 'uncertified',
                found_rank=len(native['basis']), saturated=native['saturated'],
                unresolved_primes=native['unresolved_primes'])
            outcomes[key]['attempts'].append(attempt)
            if native['certified']:
                value(row, request, native, audit['results'].get(str(key)))
                require(str(key) in audit['results'], 'missing independent audit key')
                outcomes[key]['reason'] = 'independent height and division audit'
                pending.remove(row)
    for row in pending:
        attempts = outcomes[row['key']]['attempts']
        if any(a['status'] in ('failed', 'audit_failed') for a in attempts):
            outcomes[row['key']]['reason'] = 'failed arithmetic or audit; retry pending'
        elif attempts and attempts[-1]['found_rank'] < row['rank']:
            outcomes[row['key']]['reason'] = 'rank short after bound 7'
        elif attempts:
            outcomes[row['key']]['reason'] = 'saturation or covering unresolved after bound 7'
    immutable(result, encoded(dict(input_sha256=digest(args.output / f'requests/{batch}.jsonl'),
        prepared_sha256=digest(args.output / 'prepared.json'), outcomes=list(outcomes.values()))))
    check_batch(args.output, batch, metadata)
    require(not any('retry pending' in o['reason'] for o in outcomes.values()),
            'failed arithmetic attempt retained; repeat run to retry')


def check_batch(output, batch, metadata):
    directory = output / f'batches/{batch}'
    rows = lines(output / f'requests/{batch}.jsonl')
    result = read(directory / 'result.json')
    require(result['input_sha256'] == digest(output / f'requests/{batch}.jsonl') and
            result['prepared_sha256'] == digest(output / 'prepared.json'), 'batch identity')
    require([x['key'] for x in result['outcomes']] == [r['key'] for r in rows], 'outcome keys')
    values, statuses = {}, {}
    for b in BOUNDS:
        stage = directory / f'b{b}'
        if not stage.exists():
            continue
        require(read(stage / 'prepared.json')['source_manifest'] == kernel_manifest(metadata),
                'bound source identity')
        requests, natives = recovery.native_batch(stage, '0000')
        if natives is None:
            for request in requests:
                statuses[request['key'], b] = 'failed'
            continue
        audit = recovery.audit_result(stage, '0000')
        run_producer = output / 'run-producer.json'
        if run_producer.exists():
            identity = read(run_producer)
            native_producer = read(stage / 'replay/producer.json')
            require(identity['prepared_sha256'] == digest(output / 'prepared.json') and
                    identity['kernel_sha256'] == native_producer['executable_sha256'], 'run native producer')
            if audit is not None:
                require(identity['magma_sha256'] == read(stage / 'audit/producer.json')['executable_sha256'],
                        'run audit producer')
        expected = [r for r in rows if b in bounds(r) and r['key'] not in values]
        require(requests == [dict(q=2, key=r['key'], codes=r['codes'], L=r['L'], b=b) for r in expected],
                'bound request association')
        for row, request, native in zip(expected, requests, natives):
            statuses[row['key'], b] = ('audit_failed' if audit is None else
                                       'certified' if native['certified'] else 'uncertified')
            if native['certified']:
                checked = audit['results'].get(str(row['key'])) if audit else None
                exact = value(row, request, native, checked)
                if exact is not None:
                    values[row['key']] = exact
    for row, outcome in zip(rows, result['outcomes']):
        require([a['bound'] for a in outcome['attempts']] ==
                [b for b in bounds(row) if row['key'] not in values or b <= values[row['key']]['bound']], 'attempt bounds')
        for attempt in outcome['attempts']:
            require(attempt['status'] == statuses.get((row['key'], attempt['bound'])), 'attempt status')
        require((any(a['status'] == 'certified' for a in outcome['attempts'])) == (row['key'] in values), 'certificate status')
    return rows, result['outcomes'], values


def run(args):
    metadata = prepared(args.output, args.batch)
    if args.batch is not None:
        with locked(args.output, f'.batch-{args.batch}.lock'):
            worker(args, metadata)
        return
    require(args.jobs > 0, 'jobs must be positive')
    kernel = args.kernel.resolve()
    manifest = metadata['manifest']['source_manifest']
    src = kernel.parent.parent / 'src'
    for name, expected in manifest['all_kernel_dependencies_compared'].items():
        require(digest(src / name) == expected, f'kernel dependency: {name}')
    require(subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=src, text=True).strip() ==
            manifest['calibration_head'], 'kernel revision')
    magma = Path(shutil.which(args.magma)).resolve()
    immutable(args.output / 'run-producer.json', encoded(dict(kernel_sha256=digest(kernel),
        magma_sha256=digest(magma), prepared_sha256=digest(args.output / 'prepared.json'))))
    command = shlex.join([sys.executable, str(Path(__file__).resolve()), 'run', '--output', str(args.output),
        '--kernel', str(args.kernel.resolve()), '--magma', args.magma, '--batch']) + ' {}'
    batches = [Path(name).stem for name in metadata['files']]
    args.output.joinpath('batches').mkdir(exist_ok=True)
    subprocess.run(['parallel', '-j', str(args.jobs), '--halt', 'soon,fail=1', '--resume-failed',
        '--joblog', str(args.output / 'jobs.tsv'), '--results', str(args.output / 'logs'), command],
        input='\n'.join(batches)+'\n', text=True, check=True)
    verify(args)


def verify(args):
    metadata = prepared(args.output)
    records, reasons, strata, by_bound = [], Counter(), Counter(), Counter()
    total, splits = 0, 0
    producer = dict(**metadata['identity'], prepared_sha256=digest(args.output / 'prepared.json'))
    producer_id = hashlib.sha256(encoded(producer)).hexdigest()
    for name in metadata['files']:
        batch = Path(name).stem
        rows, outcomes, values = check_batch(args.output, batch, metadata)
        for row, outcome in zip(rows, outcomes):
            key = row['key']
            exact = values.get(key)
            strata[json.dumps(row['stratum'])] += 1
            if exact:
                by_bound[exact['bound']] += 1
            total += exact is not None
            claims = dict(mw_basis=dict(status='exact' if exact else 'unresolved',
                method='native-bounded-magma-audit-v1', reason=outcome['reason'], value=exact))
            if exact and len(row['candidates']) > 1:
                sha = Fraction(row['sha_reg']) / Fraction(exact['regulator'])
                require(sha.denominator == 1, 'nonintegral Sha')
                claims.update(regulator=dict(status='exact', method='certified-basis-determinant-v1', value=exact['regulator']),
                    sha_order=dict(status='exact', method='Sha-Reg-product-v1', value=str(sha)))
                splits += 1
            reasons[outcome['reason']] += 1
            source = row['source']
            records.append(dict(schema='census-refinement', schema_version=1, q=2, key=str(key),
                key_encoding_version=1, input_ordinal=source['input_ordinal'], stage='bounded', producer_id=producer_id,
                stage1_receipt_sha256=source['stage1_receipt_sha256'],
                input_identity=dict(source_identity=metadata['manifest']['source_identity'], **source),
                status='complete' if exact else 'unresolved', claims=claims,
                evidence=dict(batch=f'batches/{batch}/result.json', attempts=outcome['attempts'])))
    require(len(records) == metadata['selected'], 'selected count')
    previous = args.output / 'refinements.jsonl'
    if previous.exists():
        current = {r['key']: r for r in records}
        for old in lines(previous):
            if old['claims']['mw_basis']['status'] == 'exact':
                require(old['claims'] == current[old['key']]['claims'], 'existing exact claim changed')
    immutable(args.output / 'import-producer.json', encoded(producer))
    save_records = b''.join(encoded(r) for r in records)
    temporary = previous.with_suffix('.new')
    temporary.write_bytes(save_records)
    temporary.replace(previous)
    stage_cpu = {}
    for stage in ('replay', 'audit'):
        user, system = 0.0, 0.0
        for path in (args.output / 'batches').glob(f'*/b*/{stage}/0000.stderr'):
            for line in path.read_text().splitlines():
                if 'User time (seconds):' in line:
                    user += float(line.split(':', 1)[1])
                if 'System time (seconds):' in line:
                    system += float(line.split(':', 1)[1])
        stage_cpu[stage] = dict(user_seconds=user, system_seconds=system)
    result = dict(selected=len(records), new_bases=total, new_splits=splits, unresolved=len(records)-total,
        certified_total=metadata['manifest']['existing_certified']+total,
        open_total=metadata['manifest']['source_count']-metadata['manifest']['existing_certified']-total,
        reasons=dict(reasons), joint_strata=dict(strata), bases_by_bound=dict(by_bound), stage_cpu=stage_cpu,
        records_sha256=digest(previous),
        prepared_sha256=digest(args.output / 'prepared.json'), evidence_bytes=sum(p.stat().st_size for p in args.output.rglob('*') if p.is_file() and p.name != 'verification.json'))
    save(args.output / 'verification.json', result)
    print(json.dumps(result, sort_keys=True))


def summarize(args):
    records, summaries, manifest = {}, [], None
    for directory in args.runs:
        directory = directory.resolve()
        verify(SimpleNamespace(output=directory))
        metadata = prepared(directory)
        if manifest is None:
            manifest = metadata['manifest']
        require(metadata['manifest'] == manifest, 'different summary populations')
        for row in lines(directory / 'refinements.jsonl'):
            require(row['key'] not in records, 'duplicate summary key; use only one pilot measurement')
            records[row['key']] = row
        summaries.append(read(directory / 'verification.json'))
    require(len(records) == manifest['input_count'], 'summary does not cover accessible queue')
    new_bases = sum(r['claims']['mw_basis']['status'] == 'exact' for r in records.values())
    new_splits = sum('sha_order' in r['claims'] for r in records.values())
    result = dict(selected=len(records), new_bases=new_bases, new_splits=new_splits,
        unresolved=len(records)-new_bases, certified_total=manifest['existing_certified']+new_bases,
        open_total=manifest['source_count']-manifest['existing_certified']-new_bases,
        remaining_splits=manifest['counts']['open_splits']-new_splits, runs=summaries)
    save(args.output / 'summary.json', result)
    print(json.dumps(result, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('export')
    p.add_argument('--evidence-root', type=Path, required=True)
    p.add_argument('--refined', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--inputs', type=Path, default=DATA / 'queue.jsonl.gz')
    p.add_argument('--manifest', type=Path, default=DATA / 'manifest.json')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--cohort', choices=['pilot', 'remaining'], required=True)
    p.add_argument('--sample-size', type=int, default=4096)
    p = sub.add_parser('run')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--kernel', type=Path, required=True)
    p.add_argument('--magma', default='magma')
    p.add_argument('--jobs', type=int, default=4)
    p.add_argument('--batch')
    p = sub.add_parser('verify')
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('summarize')
    p.add_argument('--runs', type=Path, nargs='+', required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output = args.output.resolve()
    if args.command == 'run' and args.batch:
        run(args)
    else:
        with locked(args.output):
            globals()[args.command](args)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError, subprocess.SubprocessError) as error:
        sys.exit(str(error))
