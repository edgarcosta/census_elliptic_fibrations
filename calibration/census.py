#!/usr/bin/env python3
import argparse
import concurrent.futures as cf
import gzip
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import shutil
import signal
import struct
import sys
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'scripts'))
from pipeline import Pipeline, cpu
from tables import publish_shard, read_fragment, full_row, expand

def atomic_json(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.partial')
    with tmp.open('w') as f:
        json.dump(value, f, indent=2, sort_keys=True)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    tmp.replace(path)

def source_id():
    h = hashlib.sha256()
    files = sorted([ROOT / 'census.py', ROOT / 'build.sh', *ROOT.glob('reference/*'), *ROOT.glob('scripts/*.py'), *ROOT.glob('scripts/*.m'),
                    *[p for p in (ROOT / 'src').rglob('*') if p.suffix in ('.py', '.m', '.cpp', '.hpp', '.inc') or p.name == 'spec']])
    for p in files:
        h.update(str(p.relative_to(ROOT)).encode())
        h.update(p.read_bytes())
    return h.hexdigest()

def initialize(cpus):
    os.chdir(ROOT)
    os.setsid()
    os.environ.pop("MAGMA_STARTUP_FILE", None)
    idx = multiprocessing.current_process()._identity[-1] - 1
    os.sched_setaffinity(0, {cpus[idx % len(cpus)]})
    for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[k] = '1'

def shard(task):
    q, keypath, offset, count, label, out, identity, scalar = task
    out = Path(out)
    dest = out / 'shards'
    dest.mkdir(parents=True, exist_ok=True)
    manifest = dest / (label + '.json')
    with Path(keypath).open('rb') as f:
        f.seek(offset * 8)
        data = f.read(count * 8)
    if len(data) != count * 8:
        raise RuntimeError('short key read')
    keyhash = hashlib.sha256(data).hexdigest()
    if manifest.exists():
        raise RuntimeError(f'output already exists: {manifest}; use a fresh directory')
    keys = [x[0] for x in struct.iter_unpack('<Q', data)]
    if len(set(keys)) != len(keys):
        raise RuntimeError('duplicate key within shard')
    work = out / 'work' / label
    pipe = Pipeline(work, scalar)
    before, start = cpu(), time.monotonic()
    rows = pipe.run(q, keys)
    with gzip.open(ROOT / 'reference' / f'f{q}-expected.jsonl.gz', 'rt') as f:
        expected = {r['key']: r for r in map(json.loads, f)}
    for r in rows:
        compare_records(expected[r['key']], r, f'F_{q}/{r["key"]}')
    tables = publish_shard(q, rows, dest, label, scalar)
    total = cpu() - before
    pipe.timings['writer_and_dispatch'] = max(0, total - sum(pipe.timings.values()))
    m = dict(format=2, q=q, label=label, count=count, worker_pid=os.getpid(),
             worker_cpu=next(iter(os.sched_getaffinity(0))), source=identity, keys_sha256=keyhash,
             input=dict(file=os.path.relpath(keypath, out), offset=offset),
             tables=tables, bytes=sum(t['bytes'] for t in tables.values()), cpu=total, wall=time.monotonic()-start,
             regression_compared=True, work_bytes=sum(p.stat().st_size for p in work.rglob('*') if p.is_file()),
             stages=pipe.timings, max_rss_kib=__import__('resource').getrusage(__import__('resource').RUSAGE_CHILDREN).ru_maxrss,
             split_proved=sum(r['split']['proved'] for r in rows),
             generators_proved=sum(r['generators']['certified'] for r in rows))
    atomic_json(manifest, m)
    shutil.rmtree(work)
    return m

def parallel(fn, tasks, cpus, label):
    results = []
    start = last = time.monotonic()
    with cf.ProcessPoolExecutor(max_workers=len(cpus), initializer=initialize, initargs=(cpus,)) as pool:
        pending = {pool.submit(fn, t) for t in tasks}
        try:
            while pending:
                done, pending = cf.wait(pending, timeout=1, return_when=cf.FIRST_COMPLETED)
                for f in done:
                    results.append(f.result())
                now = time.monotonic()
                if now - last >= 30 or not pending:
                    count = sum(r['count'] for r in results)
                    print(f'{label}: {len(results)}/{len(tasks)} shards, {count:,} surfaces, {now-start:.1f}s', flush=True)
                    last = now
        except BaseException:
            for p in pool._processes.values():
                try:
                    os.killpg(p.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
            raise
    return results

def compare_records(expected, actual, path=''):
    if type(expected) is not type(actual):
        raise RuntimeError(f'MISMATCH {path}: type {type(actual).__name__}, expected {type(expected).__name__}')
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            raise RuntimeError(f'MISMATCH {path}: fields differ')
        for key in expected:
            compare_records(expected[key], actual[key], path + '.' + key)
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise RuntimeError(f'MISMATCH {path}: lengths differ')
        for i, (x, y) in enumerate(zip(expected, actual)):
            compare_records(x, y, path + f'[{i}]')
    elif expected != actual:
        raise RuntimeError(f'MISMATCH {path}: {actual!r}, expected {expected!r}')

def summary(manifests):
    n = sum(m['count'] for m in manifests)
    total = sum(m['cpu'] for m in manifests)
    stages = {}
    for m in manifests:
        for k, v in m['stages'].items():
            stages[k] = stages.get(k, 0) + v / n
    return dict(count=n, cpu=total, cpu_per_surface=total/n,
                shard_wall_seconds=sorted(m['wall'] for m in manifests),
                worker_cpus=sorted({m['worker_cpu'] for m in manifests}),
                stage_cpu_per_surface=stages, compressed_bytes_per_surface=sum(m['bytes'] for m in manifests)/n,
                max_child_rss_kib=max(m['max_rss_kib'] for m in manifests),
                generators_proved=sum(m['generators_proved'] for m in manifests),
                split_proved=sum(m['split_proved'] for m in manifests),
                table_bytes_per_surface={k:sum(m['tables'][k]['bytes'] for m in manifests)/n for k in manifests[0]['tables']},
                max_work_bytes=max(m['work_bytes'] for m in manifests))

def main():
    p = argparse.ArgumentParser(description='Measure weak scaling and one fixed-total tail run.')
    p.add_argument('--workers', default='1,16,64,256')
    p.add_argument('--cpus', help='ordered, distinct CPU IDs; default: allowed physical cores first')
    p.add_argument('--tail-workers', type=int, default=256)
    p.add_argument('--limit', type=int, help='rows per fixture for a local trial')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--scalar', action='store_true')
    a = p.parse_args()
    workers = list(map(int, a.workers.split(',')))
    if not workers or min(workers) < 1 or a.tail_workers < 1 or (a.limit is not None and a.limit < a.tail_workers):
        p.error('positive workers required; --limit must cover every tail worker')
    allowed = os.sched_getaffinity(0)
    if a.cpus:
        cpus = list(map(int, a.cpus.split(',')))
    else:
        groups = {}
        for c in sorted(allowed):
            topology = Path(f'/sys/devices/system/cpu/cpu{c}/topology')
            core = tuple((topology / name).read_text().strip() for name in ('physical_package_id', 'core_id'))
            groups.setdefault(core, []).append(c)
        cpus = [c for layer in range(max(map(len, groups.values()))) for group in groups.values() for c in group[layer:layer+1]]
    if len(cpus) != len(set(cpus)) or not set(cpus) <= allowed or len(cpus) < max(*workers, a.tail_workers):
        p.error('not enough distinct allowed CPUs')
    if a.output.is_absolute() or '..' in a.output.parts:
        p.error('--output must be relative to the calibration directory')
    os.chdir(ROOT)
    a.output.mkdir(parents=True, exist_ok=False)
    os.environ.pop('MAGMA_STARTUP_FILE', None)
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[name] = '1'
    identity = source_id()
    import subprocess
    atomic_json(a.output / 'host.json', dict(source=identity, cpus=cpus, workers=workers,
                tail_workers=a.tail_workers, limit=a.limit, python=sys.version,
                cpu_model=next(s.split(':', 1)[1].strip() for s in Path('/proc/cpuinfo').read_text().splitlines() if s.startswith('model name')),
                topology=subprocess.check_output(['lscpu', '-e=CPU,CORE,SOCKET,NODE,ONLINE'], text=True),
                build=json.loads((ROOT / 'build/manifest.json').read_text())))
    cert = Pipeline(a.output / 'certificate', a.scalar)
    result = cert.call('lfun', ['--certify-tables'], stage='certify_tables')
    (a.output / 'table-certificate.txt').write_text(result.stdout)
    if result.stdout.count('TABLE_OK|') != 24:
        raise RuntimeError('table certification incomplete')
    reports = []
    for q in (2, 3):
        with gzip.open(ROOT / 'reference' / f'f{q}-expected.jsonl.gz', 'rt') as f:
            expected = list(map(json.loads, f))
        if a.limit:
            expected = expected[:a.limit]
        if len(expected) < a.tail_workers:
            p.error('fixture too small for tail worker count')
        field = a.output / f'f{q}'
        field.mkdir()
        decode_path = (field / 'decode.txt').resolve()
        decode_path.write_text(''.join(' '.join(map(str, [r['key'], *r['codes']])) + '\n' for r in expected))
        decode = Pipeline(field / 'decode', a.scalar)
        result = decode.command(['taskset', '-c', str(cpus[0]), 'magma', '-bn', f'q:={q}',
                                 f'input:={decode_path}', 'scripts/decode_check.m'], None, 'decode_check')
        (field / 'decode.log').write_text(result.stdout + result.stderr)
        if f'DECODE_OK|{len(expected)}' not in result.stdout:
            raise RuntimeError('Magma decoding check incomplete')
        keypath = (field / 'keys.bin').resolve()
        keypath.write_bytes(b''.join(struct.pack('<Q', r['key']) for r in expected))
        for mode, n in [('weak', n) for n in workers] + [('fixed', a.tail_workers)]:
            folder = field / f'{mode}-{n}'
            folder.mkdir()
            if mode == 'weak':
                slices = [(0, len(expected))] * n
            else:
                size, extra = divmod(len(expected), n)
                slices = [(i*size + min(i, extra), size + (i < extra)) for i in range(n)]
            tasks = [(q, str(keypath), offset, count, f'w{i:04d}', str(folder), identity, a.scalar)
                     for i, (offset, count) in enumerate(slices)]
            start = time.monotonic()
            ms = parallel(shard, tasks, cpus[:n], f'F_{q} {mode} {n}')
            wall = time.monotonic() - start
            for m in ms:
                records, _ = read_fragment(folder / 'shards', m)
                if q == 3:
                    records = expand(q, records, a.scalar)
                expected_rows = expected[slices[int(m['label'][1:])][0]:][:m['count']]
                for want, got in zip(sorted(expected_rows, key=lambda r: r['key']), records):
                    compare_records(full_row(want), got, f'F_{q}/{want["key"]}')
            report = summary(ms)
            report.update(q=q, mode=mode, workers=n, wall=wall,
                          rows_per_second=report['count']/wall,
                          shard_rows=[count for _, count in slices],
                          decode_check_cpu=decode.timings['decode_check'])
            atomic_json(folder / 'summary.json', report)
            reports.append(report)
            atomic_json(a.output / 'progress.json', reports)
    atomic_json(a.output / 'summary.json', dict(status='ok', source=identity, runs=reports,
                certificate_cpu=cert.timings['certify_tables']))

if __name__ == '__main__':
    main()
