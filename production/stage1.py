#!/usr/bin/env python3
"""Resumable prefix jobs for the native stage 1 CSV v1 writer."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import struct
import subprocess
import sys
import tempfile
import threading
import time

ROOT = Path(__file__).resolve().parents[1]
TOTALS = {2: 1397598, 3: 1307723786}
PREFIXES = {2: 88, 3: 1860}
HEADER = ('schema_version,q,key,kind,minimality,status,degree,sign,lower_L,'
          'fibres,torsion,sha_reg,admissible_sha,split_method,d2,sha2_dimensions,'
          'refinement_status,message').split(',')
SELECTION = {'generic_fibre_kinds': ['smooth'], 'minimalities': ['everywhere minimal']}
BUFFER = 131072
STOP = threading.Event()


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
def locked(output):
    output.mkdir(parents=True, exist_ok=True)
    with (output / '.lock').open('a') as stream:
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


class Monitor:
    def __init__(self, disk_root=None):
        self.disk_root = disk_root
        self.disk = 0
        self.aggregate = 0
        self.parent = 0
        self.worker = 0
        self.done = threading.Event()
        self.thread = threading.Thread(target=self.sample_loop)

    def sample(self):
        if self.disk_root:
            sizes = []
            for path in self.disk_root.rglob('*'):
                try:
                    if path.is_file():
                        sizes.append(path.stat().st_size)
                except FileNotFoundError:
                    pass
            self.disk = max(self.disk, sum(sizes))
        processes = {}
        for path in Path('/proc').glob('[0-9]*/status'):
            try:
                data = dict(line.split(':', 1) for line in path.read_text().splitlines())
                processes[int(path.parent.name)] = (int(data['PPid']),
                    int(data.get('VmRSS', '0 kB').split()[0]) * 1024,
                    int(data.get('VmHWM', '0 kB').split()[0]) * 1024)
            except (OSError, ValueError, KeyError):
                continue
        owned = {os.getpid()}
        while True:
            more = {pid for pid, (parent, _, _) in processes.items() if parent in owned}
            if more <= owned:
                break
            owned |= more
        self.aggregate = max(self.aggregate, sum(processes.get(pid, (0, 0, 0))[1] for pid in owned))
        self.parent = max(self.parent, processes.get(os.getpid(), (0, 0, 0))[2])
        self.worker = max([self.worker] + [processes[pid][2] for pid in owned - {os.getpid()}])

    def sample_loop(self):
        while not self.done.is_set():
            self.sample()
            self.done.wait(.05)

    def __enter__(self):
        self.thread.start()
        return self

    def __exit__(self, *args):
        self.done.set()
        self.thread.join()
        self.sample()


def environment(scratch):
    env = os.environ.copy()
    env.pop('MAGMA_STARTUP_FILE', None)
    env.update(MAGMA_TEMP_DIR=str(scratch), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', BLIS_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    return env


def execute(command, log, cpu, scratch, output=None, keys=None, q=None):
    with log.open('wb') as errors:
        with (output.open('wb') if output else log.with_suffix('.stdout.log').open('wb')) as target:
            process = subprocess.Popen(['taskset', '-c', str(cpu), *map(str, command)], cwd=ROOT,
                                       stdin=subprocess.PIPE if keys else subprocess.DEVNULL,
                                       stdout=target, stderr=errors, env=environment(scratch),
                                       start_new_session=True)
            failures = []

            def feed():
                try:
                    with process.stdin:
                        for key in key_stream(keys):
                            process.stdin.write(f'{q} {key}\n'.encode('ascii'))
                except (OSError, ValueError) as error:
                    failures.append(error)

            feeder = threading.Thread(target=feed) if keys else None
            if feeder:
                feeder.start()
            try:
                while process.poll() is None:
                    if STOP.wait(.02):
                        raise InterruptedError('stage 1 interrupted; incomplete prefixes will be recomputed')
                if process.returncode or failures:
                    raise ValueError(f'producer failed ({process.returncode}); see {log}')
            finally:
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                if feeder:
                    feeder.join()


def magma_call(args, scratch, cpu, **assignments):
    execute([args.magma, '-bn', f'field:={args.q}',
             *[f'{key}:={value}' for key, value in assignments.items()], ROOT / 'production/keys.m'],
            scratch / 'magma.log', cpu, scratch)


def counts(args, indices, scratch, cpu):
    path = scratch / 'counts.csv'
    magma_call(args, scratch, cpu, action='count', prefixes=','.join(map(str, indices)), destination=path)
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


def process_prefix(args, entry, producer, cpu):
    index = entry['prefix']
    final = args.output / 'prefixes' / f'{index:04d}'
    if (final / 'receipt.json').exists():
        return
    scratch_root = args.output / 'scratch'
    scratch = scratch_root / f'prefix-{index:04d}'
    if scratch.exists():
        failed = args.output / 'failed'
        failed.mkdir(exist_ok=True)
        archived = failed / f'prefix-{index:04d}-{time.time_ns()}'
        os.replace(scratch, archived)
        sync_dir(failed)
        sync_dir(scratch_root)
        print(f'stage1: retained incomplete attempt at {archived}', file=sys.stderr)
    scratch.mkdir()
    print(f'stage1: starting prefix {index:04d} ({entry["count"]} keys)', file=sys.stderr, flush=True)
    started = time.monotonic()
    magma_call(args, scratch, cpu, action='write', prefix=index, destination=scratch / 'keys')
    enumeration = time.monotonic() - started
    native_start = time.monotonic()
    execute([args.native, 'invariants'], scratch / 'native.log', cpu, scratch,
            output=scratch / 'rows.csv.zst', keys=scratch / f'keys_{index:04d}.bin', q=args.q)
    native_seconds = time.monotonic() - native_start
    validate(scratch, args.q, entry)
    if STOP.is_set():
        raise InterruptedError('stage 1 interrupted before receipt')
    data = {'receipt_version': 1, 'schema_version': 1, 'key_encoding_version': 1,
            'q': args.q, 'selections': SELECTION, 'prefix': entry, 'count': entry['count'],
            'producer': producer, 'enumeration_seconds': enumeration, 'native_seconds': native_seconds,
            'files': {name: digest(scratch / name) for name in file_names(index)}}
    for path in scratch.iterdir():
        if path.is_file():
            with path.open('rb') as stream:
                os.fsync(stream.fileno())
    sync_dir(scratch)
    if final.exists():
        shutil.rmtree(final)
    os.replace(scratch, final)
    sync_dir(final.parent)
    save(final / 'receipt.json', data)
    print(f'stage1: completed prefix {index:04d}', file=sys.stderr, flush=True)


def cpu_list(text, threads):
    allowed = os.sched_getaffinity(0)
    groups = {}
    for cpu in sorted(allowed):
        path = Path(f'/sys/devices/system/cpu/cpu{cpu}/topology')
        key = tuple((path / name).read_text().strip() for name in ('physical_package_id', 'core_id'))
        groups.setdefault(key, []).append(cpu)
    first = [group[0] for group in groups.values()]
    selected = [int(cpu) for cpu in text.split(',')] if text else first[:threads]
    if threads < 1 or len(selected) != threads or len(set(selected)) != threads or not set(selected) <= allowed:
        raise ValueError('cpus must name exactly one distinct allowed CPU per worker')
    cores = [key for cpu in selected for key, group in groups.items() if cpu in group]
    if len(set(cores)) != threads:
        raise ValueError('each worker requires its own physical core')
    return selected


def available_memory():
    values = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    available = int(values['MemAvailable'].split()[0]) * 1024
    # cgroup limits can be smaller than host MemAvailable in a container.
    for directory, limit, used in [('/sys/fs/cgroup', 'memory.max', 'memory.current')]:
        try:
            maximum = (Path(directory) / limit).read_text().strip()
            current = int((Path(directory) / used).read_text())
            if maximum != 'max':
                available = min(available, max(0, int(maximum) - current))
        except OSError:
            pass
    return available


def budget(args):
    worker = args.worker_memory_bytes or 1024**3
    required = (worker * args.threads + 256 * 1024**2) * 6 // 5
    if available_memory() < required:
        raise ValueError(f'insufficient available RAM: require {required} bytes')
    return required


def selection(args):
    if args.prefixes is None:
        return list(range(1, PREFIXES[args.q] + 1))
    indices = [int(value) for value in args.prefixes.split(',')]
    if not indices or len(set(indices)) != len(indices) or any(i < 1 or i > PREFIXES[args.q] for i in indices):
        raise ValueError('prefixes must be distinct valid indices')
    return sorted(indices)


def summary(args, table, producer, monitor=None, cpus=None):
    count = sum(entry['count'] for entry in table)
    partial = args.prefixes is not None
    if not partial and (len(table) != PREFIXES[args.q] or count != TOTALS[args.q]):
        raise ValueError('full census count differs from pinned expected count')
    if args.expected_count is not None and count != args.expected_count:
        raise ValueError('expected count differs from verified prefix population')
    result = {'q': args.q, 'partial': partial, 'count': count, 'prefixes': [e['prefix'] for e in table],
              'producer': producer, 'stream_buffer_bytes': BUFFER * 4,
              'gaps': ['T0 orbit/mass check is open', 'T1 angle rank is not included']}
    if monitor:
        result.update(cpus=cpus, aggregate_peak_rss_bytes=monitor.aggregate,
                      parent_peak_rss_bytes=monitor.parent, worker_peak_rss_bytes=monitor.worker)
    return result


def run(args):
    indices = selection(args)
    cpus = cpu_list(args.cpus, args.threads) if args.command != 'verify' else [min(os.sched_getaffinity(0))]
    producer = identity(args.native, args.magma)
    with locked(args.output), Monitor() as monitor:
        for name in ('scratch', 'prefixes'):
            (args.output / name).mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='count-', dir=args.output / 'scratch') as temporary:
            table = counts(args, indices, Path(temporary), cpus[0])
        if args.command != 'verify':
            budget(args)
        manifest = args.output / 'manifest.json'
        base = {'manifest_version': 1, 'q': args.q, 'producer': producer,
                'selections': SELECTION, 'schema_version': 1, 'key_encoding_version': 1,
                'package_prefix_count': PREFIXES[args.q]}
        if manifest.exists():
            saved = json.loads(manifest.read_text(encoding='ascii'))
            if any(saved.get(key) != value for key, value in base.items()):
                raise ValueError('incompatible manifest producer, selections or schema')
            known = {entry['prefix']: entry for entry in saved['prefixes']}
            if any(known.get(entry['prefix'], entry) != entry for entry in table):
                raise ValueError('package count or prefix form changed')
        else:
            if args.command == 'verify':
                raise ValueError('missing manifest')
            known = {}
        completed = set()
        for directory in sorted((args.output / 'prefixes').iterdir()):
            if STOP.is_set():
                raise InterruptedError('stage 1 interrupted during receipt revalidation')
            if (directory / 'receipt.json').exists():
                index = int(directory.name)
                if index not in known:
                    raise ValueError('receipt has no counted prefix provenance')
                print(f'stage1: revalidating prefix {index:04d}', file=sys.stderr, flush=True)
                receipt(directory, args, known[index], producer)
                if STOP.is_set():
                    raise InterruptedError('stage 1 interrupted during receipt revalidation')
                completed.add(index)
                print(f'stage1: revalidated prefix {index:04d}', file=sys.stderr, flush=True)
        for entry in table:
            if args.command == 'verify' and entry['prefix'] not in completed:
                raise ValueError('missing completed prefix receipt')
            known[entry['prefix']] = entry
        if args.command != 'verify':
            save(manifest, {**base, 'prefixes': [known[i] for i in sorted(known)]})
            queues = [table[i::args.threads] for i in range(args.threads)]

            def worker(queue, cpu):
                try:
                    for entry in queue:
                        if STOP.is_set():
                            raise InterruptedError('stage 1 interrupted')
                        process_prefix(args, entry, producer, cpu)
                except BaseException:
                    STOP.set()
                    raise

            with ThreadPoolExecutor(max_workers=args.threads) as pool:
                futures = [pool.submit(worker, queue, cpu) for queue, cpu in zip(queues, cpus)]
                errors = []
                for future in futures:
                    try:
                        future.result()
                    except Exception as error:
                        errors.append(error)
                if errors:
                    raise next((error for error in errors if not isinstance(error, InterruptedError)), errors[0])
        monitor.sample()
        result = summary(args, table, producer, monitor, cpus)
        if STOP.is_set():
            raise InterruptedError('stage 1 interrupted')
        save(args.output / ('verification.json' if args.command == 'verify' else 'summary.json'), result)
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('run', 'verify'))
    parser.add_argument('--q', type=int, required=True, choices=(2, 3))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--prefixes')
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--cpus')
    parser.add_argument('--expected-count', type=int)
    parser.add_argument('--native', type=Path, default=ROOT / 'artifacts/native/census-native')
    parser.add_argument('--magma', default=os.environ.get('MAGMA', 'magma'))
    parser.add_argument('--worker-memory-bytes', type=int, default=0)
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.native = args.native.resolve()
    args.magma = shutil.which(args.magma)
    if not args.magma:
        parser.error('magma executable is required')
    if args.worker_memory_bytes < 0:
        parser.error('worker memory must be nonnegative')
    for sig in (signal.SIGTERM, signal.SIGINT):
        signal.signal(sig, lambda *_: STOP.set())
    try:
        result = run(args)
        print(json.dumps(result, sort_keys=True))
    except (ValueError, OSError, KeyError, csv.Error, subprocess.SubprocessError) as error:
        print(f'stage1: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
