#!/usr/bin/env python3
"""Replay and independently audit retained F_2 bases before separate-record import."""

import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from stage1 import digest, environment, locked, save

ROOT = Path(__file__).resolve().parents[1]
AUDITOR = ROOT / 'production/recovery_audit.m'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_text())


def lines(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def immutable(path, data):
    if path.exists():
        require(path.read_bytes() == data, f'changed immutable file: {path}')
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + '.new')
        temporary.write_bytes(data)
        temporary.replace(path)


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n').encode()


def matrix_determinant(gram):
    a = [[Fraction(x) for x in row] for row in gram]
    r = len(a)
    require(all(len(row) == r for row in a), 'Gram shape')
    require(all(a[i][j] == a[j][i] for i in range(r) for j in range(r)), 'Gram symmetry')
    determinant = Fraction(1)
    for j in range(r):
        require(a[j][j] > 0, 'Gram positivity')
        determinant *= a[j][j]
        for i in range(j + 1, r):
            factor = a[i][j] / a[j][j]
            for k in range(j + 1, r):
                a[i][k] -= factor * a[j][k]
    return determinant


def prepared(output, batch=None):
    manifest = read(output / 'prepared.json')
    for name, expected in manifest['files'].items():
        if batch is not None and name != f'requests/{batch}.jsonl':
            continue
        require(digest(output / name) == expected, f'changed prepared input: {name}')
    return manifest


def prepare(args):
    inputs, manifest = lines(args.inputs), read(args.manifest)
    require(digest(args.inputs) == manifest['input_sha256'], 'input digest')
    require(len(inputs) == manifest['input_count'], 'input count')
    by_key = {r['key']: r for r in inputs}
    require(len(by_key) == len(inputs), 'duplicate input key')
    for r in inputs:
        a1, a2, a3, a4, a6 = r['codes']
        require(r['q'] == 2 and a2 == 0 and r['b'] in (5, 6, 7), 'request field/bound')
        require(r['key'] == (a1 << 29) + (a3 << 22) + (a4 << 13) + a6, 'request key')
    summary_path = args.evidence_root / manifest['summary']
    summary = read(summary_path)
    require(summary['complete'], 'incomplete source summary')
    sources, count, certified = {}, 0, 0
    for entry in manifest['shards']:
        path = args.evidence_root / entry['receipt']
        require(digest(path) == entry['receipt_sha256'], 'source receipt digest')
        receipt = read(path)
        shard = path.with_suffix('.jsonl.gz')
        require(receipt['q'] == 2 and receipt['source'] == summary['source'], 'source identity')
        require(digest(shard) == receipt['sha256'] == entry['compressed_sha256'], 'shard digest')
        n = 0
        with gzip.open(shard, 'rt') as stream:
            for ordinal, line in enumerate(stream):
                row = json.loads(line)
                n += 1
                certified += row['generators']['certified'] is True
                if row['key'] not in by_key:
                    continue
                request = by_key[row['key']]
                require(row['key'] not in sources, 'duplicate source key')
                require(row['q'] == 2 and row['codes'] == request['codes'] and row['L'] == request['L'], 'source key association')
                require(not row['generators']['certified'] and row['split']['proved'], 'source not eligible')
                sources[row['key']] = dict(row=row, input_ordinal=ordinal,
                    receipt=entry['receipt'], stage1_receipt_sha256=entry['receipt_sha256'],
                    compressed_sha256=entry['compressed_sha256'], keys_sha256=receipt.get('keys_sha256'))
        require(n == receipt['count'] == entry['count'], 'source count')
        count += n
    require(count == summary['count'] and certified == summary['generators_proved'], 'summary counts')
    require(set(sources) == set(by_key), 'missing source key')
    files = {}
    for index in range(0, len(inputs), 4):
        name = f'requests/{index // 4:04d}.jsonl'
        immutable(args.output / name, b''.join(encoded(r) for r in inputs[index:index + 4]))
        files[name] = digest(args.output / name)
    immutable(args.output / 'sources.json', encoded({str(k): v for k, v in sources.items()}))
    files['sources.json'] = digest(args.output / 'sources.json')
    result = dict(files=files, input_sha256=digest(args.inputs),
        source_manifest_sha256=digest(args.manifest), source_manifest=manifest,
        summary_sha256=digest(summary_path), source_identity=summary['source'],
        source_count=count, existing_certified=certified, recovery_count=len(inputs))
    immutable(args.output / 'prepared.json', encoded(result))


def native_batch(output, batch):
    metadata = prepared(output, batch)
    request_path = output / f'requests/{batch}.jsonl'
    requests = lines(request_path)
    path = output / f'replay/{batch}.jsonl'
    receipt_path = path.with_suffix('.json')
    if not receipt_path.exists():
        return requests, None
    receipt = read(receipt_path)
    require(receipt['request_sha256'] == digest(request_path), 'replay request digest')
    require(receipt['output_sha256'] == digest(path), 'replay output digest')
    producer = output / 'replay/producer.json'
    require(receipt['producer_sha256'] == digest(producer), 'replay producer digest')
    identity = read(producer)
    require(identity['prepared_sha256'] == digest(output / 'prepared.json'), 'replay input identity')
    require(identity['method'] == 'retained-native-bundle-v1', 'replay method')
    kernel = metadata['source_manifest']
    require(identity['kernel_revision'] == kernel['calibration_head']
            and identity['kernel_dependencies'] == kernel['all_kernel_dependencies_compared'],
            'replay kernel identity')
    executable_digest = identity['executable_sha256']
    require(len(executable_digest) == 64 and set(executable_digest) <= set('0123456789abcdef'),
            'replay executable digest')
    results = lines(path)
    require(len(results) == len(requests), 'replay count or duplicate key')
    for request, result in zip(requests, results):
        require((result['q'], result['key'], result['bound']) == (2, request['key'], request['b']), 'replay key or bound')
    return requests, results


def check_kernel(kernel, manifest):
    src = kernel.parent.parent / 'src'
    for name, expected in manifest['all_kernel_dependencies_compared'].items():
        require(digest(src / name) == expected, f'kernel dependency: {name}')
    require(subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=src, text=True).strip() == manifest['calibration_head'], 'kernel revision')


def parallel(args, stage):
    metadata = prepared(args.output, args.batch)
    if stage == 'replay':
        kernel = args.kernel.resolve()
        manifest = metadata['source_manifest']
        check_kernel(kernel, manifest)
        producer = dict(method='retained-native-bundle-v1', executable=str(kernel),
            executable_sha256=digest(kernel), kernel_revision=manifest['calibration_head'],
            kernel_dependencies=manifest['all_kernel_dependencies_compared'],
            prepared_sha256=digest(args.output / 'prepared.json'))
    else:
        executable = Path(shutil.which(args.magma)).resolve()
        producer = dict(method='magma-height-division-v1', script_sha256=digest(AUDITOR),
            executable=str(executable), executable_sha256=digest(executable),
            version=subprocess.check_output([executable, '-V'], text=True, env=environment(args.output)).strip(),
            replay_producer_sha256=digest(args.output / 'replay/producer.json'),
            prepared_sha256=digest(args.output / 'prepared.json'))
    directory = args.output / stage
    directory.mkdir(exist_ok=True)
    immutable(directory / 'producer.json', encoded(producer))
    batches = [Path(name).stem for name in metadata['files'] if name.startswith('requests/')]
    if args.batch is not None:
        require(args.batch in batches, 'unknown batch')
        worker(args, stage, producer)
        return
    for batch in batches:
        native_batch(args.output, batch)
    command = shlex.join([sys.executable, str(Path(__file__).resolve()), stage,
        '--output', str(args.output.resolve()), '--jobs', '1',
        *(['--kernel', str(args.kernel.resolve())] if stage == 'replay' else ['--magma', args.magma]), '--batch']) + ' {}'
    subprocess.run(['parallel', '-j', str(args.jobs), '--halt', 'soon,fail=1',
        '--resume-failed', '--joblog', str(directory / 'jobs.tsv'), '--results', str(directory / 'logs'), command],
        input='\n'.join(batches) + '\n', text=True, check=True)
    for batch in batches:
        requests, results = native_batch(args.output, batch)
        require(results is not None, f'missing replay: {batch}')
        if stage == 'audit':
            require(audit_result(args.output, batch) is not None, f'missing audit: {batch}')


def magma_literal(value):
    if isinstance(value, list):
        return '[' + ', '.join(map(magma_literal, value)) + ']'
    return str(value)


def worker(args, stage, producer):
    output, batch = args.output, args.batch
    request_path = output / f'requests/{batch}.jsonl'
    directory = output / stage
    receipt_path = directory / f'{batch}.json'
    if receipt_path.exists():
        if stage == 'replay':
            native_batch(output, batch)
        else:
            audit_result(output, batch)
        return
    if stage == 'replay':
        require(digest(args.kernel) == producer['executable_sha256'], 'changed executable')
        with request_path.open('rb') as stdin, (directory / f'{batch}.jsonl').open('wb') as stdout, (directory / f'{batch}.stderr').open('wb') as stderr:
            subprocess.run(['/usr/bin/time', '-v', str(args.kernel.resolve()), '--bundle'],
                stdin=stdin, stdout=stdout, stderr=stderr, check=True, env=environment(directory))
        save(receipt_path, dict(request_sha256=digest(request_path),
            output_sha256=digest(directory / f'{batch}.jsonl'), producer_sha256=digest(directory / 'producer.json')))
        native_batch(output, batch)
    else:
        requests, results = native_batch(output, batch)
        require(results is not None, 'missing replay')
        inputs = directory / f'{batch}.m'
        tuples = []
        for r, n in zip(requests, results):
            if not n['certified']:
                continue
            tuples.append('<' + ', '.join(magma_literal(v) for v in
                (r['key'], r['b'], r['codes'], n['basis'], n['gram'], n['rank'], n['torsion'], n['torsion_points'])) + '>')
        immutable(inputs, ('[*' + ',\n'.join(tuples) + '*]\n').encode())
        raw = directory / f'{batch}.stdout'
        with raw.open('wb') as stdout, (directory / f'{batch}.stderr').open('wb') as stderr:
            subprocess.run(['/usr/bin/time', '-v', producer['executable'], '-bn',
                'inputs:=' + str(inputs.resolve()), str(AUDITOR)], stdout=stdout, stderr=stderr,
                env=environment(directory), timeout=300, check=True)
        audited = {}
        for line in raw.read_text().splitlines():
            if line.startswith('AUDIT|'):
                _, key, checks, mu, cmax, d = line.split('|')
                require(key not in audited, 'duplicate audit key')
                audited[key] = dict(division_checks=int(checks), mu=mu, Cmax=cmax, D=int(d))
        require(set(audited) == {str(n['key']) for n in results if n['certified']}, 'audit key count')
        save(receipt_path, dict(results=audited, replay_sha256=digest(output / f'replay/{batch}.jsonl'),
            transcript_sha256=digest(raw), producer_sha256=digest(directory / 'producer.json')))


def audit_result(output, batch):
    path = output / f'audit/{batch}.json'
    if not path.exists():
        return None
    result = read(path)
    require(result['replay_sha256'] == digest(output / f'replay/{batch}.jsonl'), 'audit replay digest')
    require(result['transcript_sha256'] == digest(output / f'audit/{batch}.stdout'), 'audit transcript digest')
    producer = output / 'audit/producer.json'
    require(result['producer_sha256'] == digest(producer), 'audit producer digest')
    identity = read(producer)
    require(identity['script_sha256'] == digest(AUDITOR), 'changed audit script')
    require(identity['replay_producer_sha256'] == digest(output / 'replay/producer.json'), 'audit replay producer')
    require(identity['prepared_sha256'] == digest(output / 'prepared.json'), 'audit input identity')
    return result


def certified_value(row, request, n, checked=None):
    rank = row['derived']['rank']
    require(len(n['basis']) == len(n['gram']) == n['rank'] == rank, 'basis rank')
    det = matrix_determinant(n['gram'])
    require(row['split']['candidates'] == [[row['split']['sha'], row['split']['regulator']]], 'source split')
    require(det == Fraction(row['split']['regulator']) == Fraction(n['determinant']), 'determinant/split')
    require(Fraction(n['sha_reg']) == Fraction(row['generators']['sha_reg']) == det * row['split']['sha'], 'Sha*Reg preservation')
    require(n['torsion'] == row['generators']['torsion'], 'torsion factors')
    require(Fraction(n['mu']) < 6 + 2 * request['b'] - Fraction(n['Cmax']) == Fraction(n['missing_height_bound']), 'covering inequality')
    if checked is None:
        return None
    require(Fraction(checked['Cmax']) == Fraction(n['Cmax']) and checked['D'] == n['D'], 'local audit')
    return dict(
        basis=[[str(x) for x in point] for point in n['basis']],
        gram=dict(dimension=rank, upper=[str(Fraction(n['gram'][i][j])) for i in range(rank) for j in range(i, rank)]),
        torsion_points=[[str(x) for x in point] for point in n['torsion_points']],
        bound=request['b'], regulator=str(det), covering_mu=str(Fraction(n['mu'])),
        missing_height_bound=str(Fraction(n['missing_height_bound'])), division_audit=checked)


def import_records(args):
    metadata = prepared(args.output)
    sources = read(args.output / 'sources.json')
    records = []
    exact = 0
    producer = dict(method='retained-f2-basis-recovery-v1',
        prepared_sha256=digest(args.output / 'prepared.json'),
        importer_sha256=digest(Path(__file__)), auditor_sha256=digest(AUDITOR),
        replay_producer_sha256=digest(args.output / 'replay/producer.json'),
        audit_producer_sha256=digest(args.output / 'audit/producer.json') if (args.output / 'audit/producer.json').exists() else None)
    producer_id = hashlib.sha256(encoded(producer)).hexdigest()
    for name in metadata['files']:
        if not name.startswith('requests/'):
            continue
        batch = Path(name).stem
        requests, native = native_batch(args.output, batch)
        audit = audit_result(args.output, batch) if native else None
        for i, request in enumerate(requests):
            key = str(request['key'])
            source = sources[key]
            row = source['row']
            claim = dict(status='unresolved', method='retained-native-covering-and-magma-v1',
                         reason='missing replay or independent audit', value=None)
            n = native[i] if native else None
            if n and not n['certified']:
                claim['reason'] = 'native search did not certify a basis'
            if n and n['certified']:
                value = certified_value(row, request, n, audit['results'].get(key) if audit else None)
                if value is not None:
                    exact += 1
                    claim.update(status='exact', reason='full independent height and division audit', value=value)
            records.append(dict(schema='census-refinement', schema_version=1, q=2, key=key,
                key_encoding_version=1, input_ordinal=source['input_ordinal'], stage='bounded',
                producer_id=producer_id, stage1_receipt_sha256=source['stage1_receipt_sha256'],
                input_identity=dict(source_identity=metadata['source_identity'], receipt=source['receipt'],
                    compressed_sha256=source['compressed_sha256'], keys_sha256=source['keys_sha256']),
                status='complete' if claim['status'] == 'exact' else 'unresolved', claims=dict(mw_basis=claim),
                evidence=dict(request=name, replay=f'replay/{batch}.json', audit=f'audit/{batch}.json' if audit else None)))
    previous = args.output / 'refinements.jsonl'
    if previous.exists():
        current = {r['key']: r for r in records}
        for old in lines(previous):
            claim = old['claims']['mw_basis']
            if claim['status'] == 'exact':
                require(old['key'] in current and current[old['key']]['claims']['mw_basis'] == claim,
                        'import would replace an existing exact basis')
    save(args.output / 'import-producer.json', producer)
    data = b''.join(encoded(record) for record in records)
    path = args.output / 'refinements.jsonl'
    temporary = path.with_suffix('.new')
    temporary.write_bytes(data)
    temporary.replace(path)
    save(args.output / 'import.json', dict(producer_id=producer_id, records_sha256=digest(path),
        selected=len(records), recovered=exact, unresolved=len(records) - exact,
        certified_total=metadata['existing_certified'] + exact,
        open_total=metadata['source_count'] - metadata['existing_certified'] - exact, new_splits_credited=0))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--inputs', type=Path, default=ROOT / 'data/mw/f2-recovery/replay-candidates.jsonl')
    p.add_argument('--manifest', type=Path, default=ROOT / 'data/mw/f2-recovery/source-manifest.json')
    p.add_argument('--evidence-root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    for stage in ('replay', 'audit'):
        p = sub.add_parser(stage)
        p.add_argument('--output', type=Path, required=True)
        p.add_argument('--jobs', type=int, default=4)
        p.add_argument('--batch')
        if stage == 'replay':
            p.add_argument('--kernel', type=Path, required=True)
        else:
            p.add_argument('--magma', default='magma')
    p = sub.add_parser('import')
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output = args.output.resolve()
    if args.command in ('replay', 'audit'):
        require(args.jobs > 0, 'jobs must be positive')
        if args.batch is None:
            with locked(args.output, f'.{args.command}.lock'):
                parallel(args, args.command)
        else:
            with locked(args.output, f'.{args.command}-{args.batch}.lock'):
                parallel(args, args.command)
    else:
        with locked(args.output):
            (prepare if args.command == 'prepare' else import_records)(args)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError, subprocess.SubprocessError) as error:
        sys.exit(str(error))
