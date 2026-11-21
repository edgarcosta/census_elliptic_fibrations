#!/usr/bin/env python3
"""Resumable prefix jobs for the native stage 1 CSV v1 writer."""

import argparse
from contextlib import contextmanager
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
TOTALS = {2: 1397598, 3: 1307723786}
PREFIXES = {2: 88, 3: 1860}
HEADER = ('schema_version,q,key,kind,minimality,status,degree,sign,lower_L,'
          'fibres,torsion,sha_reg,admissible_sha,split_method,d2,sha2_dimensions,'
          'refinement_status,message').split(',')
SELECTION = {'generic_fibre_kinds': ['smooth'], 'minimalities': ['everywhere minimal']}
BUFFER = 131072


def digest(path):
    with open(path, 'rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def identity(binary, magma):
    sha = hashlib.sha256()
    files = [ROOT / 'spec']
    for directory, patterns in [('magma', ('*.m', 'spec')), ('native', ('*.cpp', '*.hpp')),
                                ('production', ('*.py', '*.m'))]:
        for pattern in patterns:
            files.extend((ROOT / directory).rglob(pattern))
    for path in sorted(files):
        sha.update(str(path.relative_to(ROOT)).encode('ascii') + b'\0')
        sha.update(path.read_bytes())
    return {'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT,
                                              text=True).strip(),
            'source_sha256': sha.hexdigest(), 'executable_sha256': digest(binary),
            'magma_sha256': digest(magma),
            'magma_version': subprocess.check_output([magma, '-V'], text=True,
                                                      env=environment(ROOT), timeout=30).strip()}


def sync_dir(path):
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def save(path, data):
    temporary = path.with_suffix(path.suffix + '.new')
    with temporary.open('w', encoding='ascii') as stream:
        json.dump(data, stream, sort_keys=True, indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    sync_dir(path.parent)


@contextmanager
def locked(output, name='.lock'):
    output.mkdir(parents=True, exist_ok=True)
    with (output / name).open('a') as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('output is locked by another command') from None
        yield


def record_lines(path):
    listing = subprocess.run(['zstd', '--list', path], check=False,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if listing.returncode:
        raise ValueError('zstd --list: ' + listing.stderr.strip())
    fields = listing.stdout.splitlines()[1].split()
    if fields[:2] != ['1', '0'] or 'XXH64' not in fields:
        raise ValueError('expected one zstd frame with a content checksum and no skips')
    with subprocess.Popen(['zstd', '-dc', path], stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE) as process:
        try:
            for line in iter(lambda: process.stdout.readline(BUFFER + 1), b''):
                if len(line) > BUFFER:
                    raise ValueError('oversized CSV record')
                if not line.endswith(b'\n'):
                    raise ValueError('missing final CSV newline')
                yield line.decode('ascii')
            errors = process.stderr.read().decode('ascii')
            if process.wait():
                raise ValueError(f'zstd ({process.returncode}): ' + errors.strip())
        finally:
            if process.poll() is None:
                process.terminate()


def key_stream(path):
    with path.open('rb') as stream:
        while block := stream.read(BUFFER):
            if len(block) % 8:
                raise ValueError('truncated binary key')
            for (key,) in struct.iter_unpack('<Q', block):
                yield key


def validate(directory, q, entry):
    index, count = entry['prefix'], entry['count']
    name = f'keys_{index:04d}'
    keys = directory / (name + '.bin')
    if keys.stat().st_size != 8 * count:
        raise ValueError('key count differs from package count')
    degrees = (2, 6) if q == 2 else (4, 8)
    expected = {'file_format_version': '1', 'key_encoding_version': '1',
                'characteristic': str(q), 'prefix_index': str(index),
                'prefix_form_1_degree': str(degrees[0]),
                'prefix_form_2_degree': str(degrees[1]),
                'prefix_form_1_key': str(entry['form_1_key']),
                'prefix_form_2_key': str(entry['form_2_key']),
                'generic_fibre_kinds': 'smooth', 'minimalities': 'everywhere minimal',
                'number_of_keys': str(count), 'cell': f'smooth|everywhere minimal|{count}',
                'complete': '1'}
    lines = (directory / (name + '.stats')).read_text(encoding='ascii').splitlines()
    actual = dict(line.split('=', 1) for line in lines)
    if actual != expected or len(lines) != len(expected):
        raise ValueError('key metadata differs from package count or prefix forms')
    with (directory / 'keys.csv').open(encoding='ascii', newline='') as stream:
        table = list(csv.DictReader(stream))
    wanted = {k: v for k, v in expected.items() if k in (
        'file_format_version', 'key_encoding_version', 'characteristic', 'prefix_index',
        'prefix_form_1_degree', 'prefix_form_1_key', 'prefix_form_2_degree', 'prefix_form_2_key',
        'number_of_keys')}
    wanted.update(generic_fibre_kind='smooth', minimality='everywhere minimal')
    if table != [wanted]:
        raise ValueError('key CSV metadata differs from prefix forms')
    records = csv.DictReader(record_lines(directory / 'rows.csv.zst'), strict=True)
    if records.fieldnames != HEADER:
        raise ValueError('wrong stage 1 CSV header')
    ordered_keys = key_stream(keys)
    seen = 0
    for row in records:
        key = next(ordered_keys, None)
        if (key is None or row.get('key') != str(key) or row.get('q') != str(q)
                or row.get('schema_version') != '1' or row.get('status') != 'ok'
                or row.get('kind') != 'smooth' or row.get('minimality') != 'everywhere minimal'
                or None in row or any(value is None for value in row.values())):
            raise ValueError('record status, classification or ordered key correspondence failed')
        if (not all(row[field] for field in ('degree', 'sign', 'lower_L', 'fibres', 'torsion',
                                            'sha_reg', 'admissible_sha', 'split_method'))
                or row['d2'] or row['sha2_dimensions'] or row['message']
                or row['split_method'] not in ('constraints', 'unresolved')):
            raise ValueError('incomplete stage 1 arithmetic')
        seen += 1
    if seen != count or next(ordered_keys, None) is not None:
        raise ValueError('record count differs from package count')
    return seen


def environment(scratch):
    env = os.environ.copy()
    env.pop('MAGMA_STARTUP_FILE', None)
    env.update(MAGMA_TEMP_DIR=str(scratch), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', BLIS_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    return env


def execute(command, log, scratch, output=None, keys=None, q=None):
    with tempfile.TemporaryFile(dir=scratch) as incoming:
        if keys:
            for key in key_stream(keys):
                incoming.write(f'{q} {key}\n'.encode('ascii'))
            incoming.seek(0)
        with log.open('wb') as errors:
            target_path = output or log.with_suffix('.stdout.log')
            with target_path.open('wb') as target:
                result = subprocess.run(list(map(str, command)), cwd=ROOT,
                    stdin=incoming if keys else subprocess.DEVNULL,
                    stdout=target, stderr=errors, env=environment(scratch))
        if result.returncode:
            raise ValueError(f'producer failed ({result.returncode}); see {log}')


def magma_call(args, scratch, **assignments):
    execute([args.magma, '-bn', f'field:={args.q}',
             *[f'{key}:={value}' for key, value in assignments.items()], ROOT / 'production/keys.m'],
            scratch / 'magma.log', scratch)


def counts(args, indices, scratch):
    path = scratch / 'counts.csv'
    magma_call(args, scratch, action='count', prefixes=','.join(map(str, indices)), destination=path)
    with path.open(encoding='ascii', newline='') as stream:
        table = [{key: int(value) for key, value in row.items()} for row in csv.DictReader(stream)]
    if [row['prefix'] for row in table] != indices or any(row['count'] < 0 for row in table):
        raise ValueError('package prefix counts are incomplete')
    if (indices == list(range(1, PREFIXES[args.q] + 1))
            and sum(row['count'] for row in table) != TOTALS[args.q]):
        raise ValueError('full census count differs from pinned expected count')
    return table


def file_names(index):
    return ['keys.csv', f'keys_{index:04d}.bin', f'keys_{index:04d}.stats', 'rows.csv.zst']


def receipt(directory, args, entry, producer):
    path = directory / 'receipt.json'
    data = json.loads(path.read_text(encoding='ascii'))
    wanted = {'receipt_version': 1, 'schema_version': 1, 'key_encoding_version': 1,
              'q': args.q, 'selections': SELECTION, 'prefix': entry, 'count': entry['count'],
              'producer': producer}
    if any(data.get(key) != value for key, value in wanted.items()):
        raise ValueError(f'incompatible or corrupt receipt: {path}')
    files = {name: digest(directory / name) for name in file_names(entry['prefix'])}
    if data.get('files') != files:
        raise ValueError(f'file digest differs from receipt: {path}')
    validate(directory, args.q, entry)
    return data


def process_prefix(args, entry, producer):
    index = entry['prefix']
    final = args.output / 'prefixes' / f'{index:04d}'
    if (final / 'receipt.json').exists():
        receipt(final, args, entry, producer)
        return
    scratch_root = args.output / 'scratch'
    scratch = scratch_root / f'prefix-{index:04d}'
    for incomplete in (scratch, final):
        if incomplete.exists():
            failed = args.output / 'failed'
            failed.mkdir(exist_ok=True)
            archived = failed / f'prefix-{index:04d}-{time.time_ns()}'
            os.replace(incomplete, archived)
            sync_dir(failed)
            sync_dir(incomplete.parent)
            print(f'stage1: retained incomplete attempt at {archived}', file=sys.stderr)
    scratch.mkdir()
    print(f'stage1: starting prefix {index:04d} ({entry["count"]} keys)', file=sys.stderr, flush=True)
    started = time.monotonic()
    magma_call(args, scratch, action='write', prefix=index, destination=scratch / 'keys')
    enumeration = time.monotonic() - started
    native_start = time.monotonic()
    execute([args.native, 'invariants'], scratch / 'native.log', scratch,
            output=scratch / 'rows.csv.zst', keys=scratch / f'keys_{index:04d}.bin', q=args.q)
    native_seconds = time.monotonic() - native_start
    validate(scratch, args.q, entry)
    data = {'receipt_version': 1, 'schema_version': 1, 'key_encoding_version': 1,
            'q': args.q, 'selections': SELECTION, 'prefix': entry, 'count': entry['count'],
            'producer': producer, 'enumeration_seconds': enumeration, 'native_seconds': native_seconds,
            'files': {name: digest(scratch / name) for name in file_names(index)}}
    for path in scratch.iterdir():
        if path.is_file():
            with path.open('rb') as stream:
                os.fsync(stream.fileno())
    sync_dir(scratch)
    os.replace(scratch, final)
    sync_dir(final.parent)
    save(final / 'receipt.json', data)
    print(f'stage1: completed prefix {index:04d}', file=sys.stderr, flush=True)


def selection(args):
    if args.prefixes is None:
        return list(range(1, PREFIXES[args.q] + 1))
    indices = [int(value) for value in args.prefixes.split(',')]
    if not indices or len(set(indices)) != len(indices) or any(i < 1 or i > PREFIXES[args.q] for i in indices):
        raise ValueError('prefixes must be distinct valid indices')
    return sorted(indices)


def summary(args, table, producer):
    count = sum(entry['count'] for entry in table)
    partial = args.prefixes is not None
    if not partial and (len(table) != PREFIXES[args.q] or count != TOTALS[args.q]):
        raise ValueError('full census count differs from pinned expected count')
    if args.expected_count is not None and count != args.expected_count:
        raise ValueError('expected count differs from verified prefix population')
    result = {'q': args.q, 'partial': partial, 'count': count, 'prefixes': [e['prefix'] for e in table],
              'producer': producer,
              'gaps': ['independent smooth/minimal label check is open', 'T1 angle rank is not included']}
    return result


def manifest_base(args, producer):
    return {'manifest_version': 1, 'q': args.q, 'producer': producer,
            'selections': SELECTION, 'schema_version': 1, 'key_encoding_version': 1,
            'package_prefix_count': PREFIXES[args.q]}


def read_manifest(args, producer):
    path = args.output / 'manifest.json'
    if not path.exists():
        raise ValueError('missing manifest; run prepare first')
    saved = json.loads(path.read_text(encoding='ascii'))
    if any(saved.get(key) != value for key, value in manifest_base(args, producer).items()):
        raise ValueError('incompatible manifest producer, selections or schema')
    return saved['prefixes']


def counted(args, indices):
    temporary = Path(tempfile.mkdtemp(prefix='count-', dir=args.output / 'scratch'))
    table = counts(args, indices, temporary)
    shutil.rmtree(temporary)
    return table


def prepare(args):
    indices = selection(args)
    producer = identity(args.native, args.magma)
    with locked(args.output):
        for name in ('scratch', 'prefixes'):
            (args.output / name).mkdir(exist_ok=True)
        table = counted(args, indices)
        result = summary(args, table, producer)
        manifest = args.output / 'manifest.json'
        if manifest.exists() and read_manifest(args, producer) != table:
            raise ValueError('package count or prefix form changed')
        if not manifest.exists():
            save(manifest, {**manifest_base(args, producer), 'prefixes': table})
        return result


def run(args):
    if args.prefixes is None or len(indices := selection(args)) != 1:
        raise ValueError('run requires exactly one prefix')
    producer = identity(args.native, args.magma)
    table = read_manifest(args, producer)
    entry = next((entry for entry in table if entry['prefix'] == indices[0]), None)
    if entry is None:
        raise ValueError('prefix is absent from prepared manifest')
    with locked(args.output / 'locks', f'{indices[0]:04d}.lock'):
        process_prefix(args, entry, producer)
    return summary(args, [entry], producer)


def verify(args):
    (args.output / 'verification.json').unlink(missing_ok=True)
    indices = selection(args)
    producer = identity(args.native, args.magma)
    with locked(args.output):
        (args.output / 'scratch').mkdir(exist_ok=True)
        table = counted(args, indices)
        known = {entry['prefix']: entry for entry in read_manifest(args, producer)}
        if any(known.get(entry['prefix']) != entry for entry in table):
            raise ValueError('package count or prefix form changed')
        if any(not (args.output / 'prefixes' / f'{entry["prefix"]:04d}' / 'receipt.json').exists()
               for entry in table):
            raise ValueError('missing completed prefix receipt')
        for directory in sorted((args.output / 'prefixes').iterdir()):
            if (directory / 'receipt.json').exists():
                index = int(directory.name)
                if index not in known:
                    raise ValueError('receipt has no counted prefix provenance')
                print(f'stage1: revalidating prefix {index:04d}', file=sys.stderr, flush=True)
                receipt(directory, args, known[index], producer)
                print(f'stage1: revalidated prefix {index:04d}', file=sys.stderr, flush=True)
        result = summary(args, table, producer)
        save(args.output / 'verification.json', result)
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'run', 'verify'))
    parser.add_argument('--q', type=int, required=True, choices=(2, 3))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--prefixes')
    parser.add_argument('--expected-count', type=int)
    parser.add_argument('--native', type=Path, default=ROOT / 'artifacts/native/census-native')
    parser.add_argument('--magma', default=os.environ.get('MAGMA', 'magma'))
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.native = args.native.resolve()
    args.magma = shutil.which(args.magma)
    if not args.magma:
        parser.error('magma executable is required')
    try:
        result = {'prepare': prepare, 'run': run, 'verify': verify}[args.command](args)
        print(json.dumps(result, sort_keys=True))
    except (ValueError, OSError, KeyError, csv.Error, subprocess.SubprocessError) as error:
        print(f'stage1: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
