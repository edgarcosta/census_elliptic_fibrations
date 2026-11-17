import csv
import fcntl
import hashlib
import io
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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests/compare"))
from records import decompress

CLI = ROOT / "production/stage1.py"
BINARY = ROOT / "artifacts/native/census-native"


def call(command, q, output, *extra, ok=True):
    result = subprocess.run([sys.executable, CLI, command, '--q', str(q),
                             '--output', output, *map(str, extra)],
                            cwd=ROOT, capture_output=True, text=True)
    assert (result.returncode == 0) == ok, result.stdout + result.stderr
    return result


def compress(path, text, checksum=True):
    with path.open('wb') as output:
        subprocess.run(['zstd', '-q', '--check' if checksum else '--no-check'],
                       input=text, stdout=output, check=True)


def rows(output, prefix):
    return list(csv.DictReader(io.StringIO(decompress(
        (output / 'prefixes' / f'{prefix:04d}' / 'rows.csv.zst').read_bytes()))))


def wait_for(predicate, process):
    deadline = time.monotonic() + 120
    while not predicate():
        assert process.poll() is None, process.communicate()
        assert time.monotonic() < deadline, 'worker did not reach interruption point'
        time.sleep(.02)


def worker_args():
    physical = {}
    for cpu in sorted(os.sched_getaffinity(0)):
        topology = Path(f'/sys/devices/system/cpu/cpu{cpu}/topology')
        key = tuple((topology / name).read_text().strip()
                    for name in ('physical_package_id', 'core_id'))
        physical.setdefault(key, cpu)
    cpus = list(physical.values())[:2]
    assert len(cpus) == 2
    return ['--threads', '2', '--cpus', ','.join(map(str, cpus))]


def test_revalidation_stop(base):
    output = base / 'revalidated'
    call('run', 3, output, '--prefixes', '2,3')
    summary = output / 'summary.json'
    summary.unlink()
    path = output / 'prefixes/0003/receipt.json'
    data = json.loads(path.read_text())
    data['count'] += 1
    path.write_text(json.dumps(data))
    path = output / 'prefixes/0002/keys.csv'
    retained = path.read_bytes()
    path.unlink()
    os.mkfifo(path)
    process = subprocess.Popen([sys.executable, CLI, 'run', '--q', '3', '--output', output,
                                '--prefixes', '2,3'], cwd=ROOT,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    failures = []

    def slow_storage():
        try:
            # Separate FIFOs give the digest and validation readers their own EOF.
            # Send the stop request while the digest reader is waiting for data.
            for attempt in range(2):
                with path.open('wb', buffering=0) as stream:
                    if attempt == 0:
                        process.send_signal(signal.SIGTERM)
                    stream.write(retained)
                    if attempt == 0:
                        path.unlink()
                        os.mkfifo(path)
        except OSError as error:
            failures.append(error)

    feeder = threading.Thread(target=slow_storage, daemon=True)
    feeder.start()
    try:
        stdout, stderr = process.communicate(timeout=30)
        feeder.join(timeout=1)
        assert not feeder.is_alive() and not failures, failures
        assert process.returncode != 0, stdout
        assert b'stage 1 interrupted' in stderr, stderr
        assert b'revalidating prefix 0002' in stderr, stderr
        assert b'revalidating prefix 0003' not in stderr, stderr
        assert not summary.exists()
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()


def test_changed_full_count(base):
    wrapper = base / 'changed-count-magma'
    marker = base / 'enumeration-started'
    magma = shutil.which(os.environ.get('MAGMA', 'magma'))
    wrapper.write_text('#!/usr/bin/env python3\nimport csv,pathlib,subprocess,sys\n'
                       f'if "-V" in sys.argv: sys.exit(subprocess.run([{magma!r}, "-V"]).returncode)\n'
                       'if "action:=write" in sys.argv:\n'
                       f' pathlib.Path({str(marker)!r}).touch(); sys.exit(7)\n'
                       f'r=subprocess.run([{magma!r}, *sys.argv[1:]])\n'
                       'if r.returncode == 0 and "action:=count" in sys.argv:\n'
                       ' p=pathlib.Path(next(a.split(":=",1)[1] for a in sys.argv if a.startswith("destination:=")))\n'
                       ' with p.open() as f: reader=csv.DictReader(f); fields=reader.fieldnames; rows=list(reader)\n'
                       ' rows[1]["count"]=str(int(rows[1]["count"])+1)\n'
                       ' with p.open("w") as f: writer=csv.DictWriter(f,fields); writer.writeheader(); writer.writerows(rows)\n'
                       'sys.exit(r.returncode)\n')
    wrapper.chmod(0o755)
    for command in ('run', 'pilot'):
        output = base / command
        result = call(command, 3, output, '--magma', wrapper, ok=False)
        assert 'full census count differs from pinned expected count' in result.stderr, result.stderr
        assert not marker.exists(), 'enumeration started with a changed full count'
        assert not (output / 'manifest.json').exists()
        assert not (output / 'pilot.json').exists()


def test_worker_error(base):
    reference = base / 'reference'
    call('run', 3, reference, *worker_args(), '--prefixes', '3')
    first = int(rows(reference, 3)[0]['key'])
    held = base / 'held'
    wrapper = base / 'fail-second'
    wrapper.write_text('#!/usr/bin/env python3\nimport pathlib,signal,sys,time\n'
                       'data=sys.stdin.buffer.read()\n'
                       f'held=pathlib.Path({str(held)!r})\n'
                       f'if data.startswith(b"3 {first}\\n"):\n'
                       ' while not held.exists(): time.sleep(.01)\n'
                       ' print("prefix three failure", file=sys.stderr)\n sys.exit(7)\n'
                       'held.touch()\nsignal.pause()\n')
    wrapper.chmod(0o755)
    output = base / 'failed'
    result = call('run', 3, output, *worker_args(), '--prefixes', '2,3',
                  '--native', wrapper, ok=False)
    assert 'producer failed (7)' in result.stderr, result.stderr
    assert 'prefix-0003/native.log' in result.stderr, result.stderr
    assert not (output / 'summary.json').exists()


def test_failed_logs(base):
    marker = base / 'failed-once'
    wrapper = base / 'fail-once'
    wrapper.write_text('#!/usr/bin/env python3\nimport pathlib,subprocess,sys\n'
                       f'marker=pathlib.Path({str(marker)!r})\n'
                       'if not marker.exists():\n marker.touch()\n'
                       ' print("retained failure evidence", file=sys.stderr)\n sys.exit(7)\n'
                       f'sys.exit(subprocess.run([{str(BINARY)!r}, "invariants"]).returncode)\n')
    wrapper.chmod(0o755)
    output = base / 'retried'
    args = ['--prefixes', '2', '--native', wrapper]
    call('run', 3, output, *args, ok=False)
    original = (output / 'scratch/prefix-0002/native.log').read_bytes()
    assert b'retained failure evidence' in original
    call('run', 3, output, *args)
    assert any(p.read_bytes() == original for p in output.rglob('native.log')), 'failed log was deleted'
    call('verify', 3, output, *args)


def test_publication_recovery(base):
    output = base / 'published'
    call('run', 3, output, '--prefixes', '1,2')
    path = output / 'prefixes/0002/receipt.json'
    original = json.loads(path.read_text())
    path.unlink()
    call('verify', 3, output, '--prefixes', '1,2', ok=False)
    call('run', 3, output, '--prefixes', '1,2')
    assert path.exists(), 'resume skipped a published directory without its receipt'
    assert json.loads(path.read_text())['files'] == original['files']
    call('verify', 3, output, '--prefixes', '1,2')


def test_magma_version(base):
    version = base / 'version'
    version.write_text('review-version-one')
    wrapper = base / 'magma-launcher'
    magma = shutil.which(os.environ.get('MAGMA', 'magma'))
    wrapper.write_text('#!/usr/bin/env python3\nimport os,pathlib,sys\n'
                       f'if "-V" in sys.argv:\n print(pathlib.Path({str(version)!r}).read_text()); sys.exit()\n'
                       f'os.execv({magma!r}, [{magma!r}, *sys.argv[1:]])\n')
    wrapper.chmod(0o755)
    output = base / 'versioned'
    call('run', 3, output, '--prefixes', '1', '--magma', wrapper)
    version.write_text('review-version-two')
    result = call('run', 3, output, '--prefixes', '1', '--magma', wrapper, ok=False)
    assert 'incompatible manifest' in result.stderr


def test_magma_spec(base):
    shutil.copy2(ROOT / 'spec', base / 'spec')
    for directory in ('magma', 'production'):
        shutil.copytree(ROOT / directory, base / directory)
    command = [sys.executable, base / 'production/stage1.py', 'run', '--q', '3',
               '--output', base / 'results', '--prefixes', '1', '--native', BINARY]
    result = subprocess.run(command, cwd=base, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    spec = base / 'magma/spec'
    spec.write_text(spec.read_text() + '\n')
    result = subprocess.run(command, cwd=base, capture_output=True, text=True)
    assert result.returncode != 0, result.stdout + result.stderr
    assert 'incompatible manifest' in result.stderr, result.stdout + result.stderr


def test_legacy_pilot_budget(base):
    result = call('run', 3, base / 'identity', '--prefixes', '1')
    producer = json.loads(result.stdout)['producer']
    path = base / 'old-pilot.json'
    path.write_text(json.dumps({'passed': True, 'q': 3, 'producer': producer,
                               'worker_peak_rss_bytes': 10**18, 'parent_peak_rss_bytes': 1}))
    result = call('run', 3, base / 'full', '--pilot-data', path, ok=False)
    assert 'largest-prefix memory measurement' in result.stderr, result.stderr


def test_repilot(base):
    output = base / 'pilot'
    args = worker_args()
    call('pilot', 3, output, *args)
    wrapper = base / 'native-launcher'
    wrapper.write_text('#!/usr/bin/env python3\nimport os,sys\n'
                       f'os.execv({str(BINARY)!r}, [{str(BINARY)!r}, *sys.argv[1:]])\n')
    wrapper.chmod(0o755)
    call('pilot', 3, output, *args, '--native', wrapper)
    assert json.loads((output / 'pilot.json').read_text())['passed']


def test_commands(base):
    allowed = sorted(os.sched_getaffinity(0))
    physical = {}
    for cpu in allowed:
        topology = Path(f'/sys/devices/system/cpu/cpu{cpu}/topology')
        key = tuple((topology / name).read_text().strip()
                    for name in ('physical_package_id', 'core_id'))
        physical.setdefault(key, cpu)
    cpus = list(physical.values())[:2]
    args = ['--threads', str(len(cpus)), '--cpus', ','.join(map(str, cpus))]
    fresh = base / 'fresh'
    result = call('run', 3, fresh, *args, '--prefixes', '1,2,3')
    summary = json.loads(result.stdout)
    assert summary['partial'] and summary['count'] == 4494
    assert summary['cpus'] == cpus
    assert summary['aggregate_peak_rss_bytes'] > 0
    assert [len(rows(fresh, i)) for i in (1, 2, 3)] == [0, 1122, 3372]
    call('verify', 3, fresh, '--prefixes', '1,2,3', '--expected-count', '4494')
    call('verify', 3, fresh, '--expected-count', '4494', ok=False)
    call('verify', 3, fresh, '--prefixes', '1,2,3', '--expected-count', '1307723786', ok=False)

    # Hold the real computation before its wrapper exits, so no receipt can race the stop.
    resumed = base / 'resumed'
    paused = base / 'pause-native'
    marker = base / 'native-held'
    first = int(rows(fresh, 3)[0]['key'])
    paused.write_text('#!/usr/bin/env python3\nimport pathlib,signal,subprocess,sys\n'
                      'data=sys.stdin.buffer.read()\n'
                      f'r=subprocess.run([{str(BINARY)!r}, "invariants"], input=data, stdout=sys.stdout.buffer)\n'
                      f'm=pathlib.Path({str(marker)!r})\n'
                      f'if data.startswith(b"3 {first}\\n") and not m.exists():\n'
                      ' m.touch()\n signal.pause()\n'
                      'sys.exit(r.returncode)\n')
    paused.chmod(0o755)
    recovery_args = [*args, '--native', str(paused)]
    call('run', 3, resumed, *recovery_args, '--prefixes', '1,2')
    receipt = resumed / 'prefixes/0002/receipt.json'
    original = receipt.read_bytes()
    process = subprocess.Popen([sys.executable, CLI, 'run', '--q', '3', '--output', resumed,
                                *recovery_args, '--prefixes', '1,2,3'], cwd=ROOT,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        wait_for(marker.exists, process)
        process.send_signal(signal.SIGTERM)
        process.communicate(timeout=30)
        assert process.returncode != 0
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()
    assert not (resumed / 'prefixes/0003/receipt.json').exists()
    call('run', 3, resumed, *recovery_args, '--prefixes', '1,2,3')
    assert receipt.read_bytes() == original
    for i in (1, 2, 3):
        assert rows(resumed, i) == rows(fresh, i)
        keys = [int(r['key']) for r in rows(resumed, i)]
        assert len(keys) == len(set(keys))
    single = base / 'single'
    call('run', 3, single, '--threads', '1', '--cpus', cpus[0], '--prefixes', '1,2,3')
    assert all(rows(single, i) == rows(fresh, i) for i in (1, 2, 3))

    for fault in ('receipt', 'truncated', 'unchecked', 'key', 'identity', 'row'):
        damaged = base / fault
        shutil.copytree(fresh, damaged)
        directory = damaged / 'prefixes/0002'
        path = directory / 'receipt.json'
        data = json.loads(path.read_text())
        if fault == 'receipt':
            data['count'] += 1
        elif fault == 'identity':
            data['producer']['executable_sha256'] = '0' * 64
        elif fault == 'truncated':
            p = directory / 'rows.csv.zst'
            p.write_bytes(p.read_bytes()[:-4])
            data['files'][p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
        elif fault == 'unchecked':
            p = directory / 'rows.csv.zst'
            compress(p, decompress(p.read_bytes()).encode('ascii'), checksum=False)
            data['files'][p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
        else:
            p = directory / 'keys_0002.bin'
            blob = bytearray(p.read_bytes())
            if fault == 'key':
                blob[:8], blob[8:16] = blob[8:16], blob[:8]
                p.write_bytes(blob)
                data['files'][p.name] = hashlib.sha256(blob).hexdigest()
            else:
                text = decompress((directory / 'rows.csv.zst').read_bytes())
                text = text.replace(',ok,', ',error,', 1)
                # A valid frame must not make a failed arithmetic row complete.
                compress(directory / 'rows.csv.zst', text.encode('ascii'))
                p = directory / 'rows.csv.zst'
                data['files'][p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
        path.write_text(json.dumps(data))
        result = call('verify', 3, damaged, '--prefixes', '1,2,3', ok=False)
        if fault == 'truncated':
            assert 'zstd --list: ' in result.stderr, result.stderr
        call('run', 3, damaged, *args, '--prefixes', '1,2,3', ok=False)
        assert json.loads(path.read_text()) == data

    with open(fresh / '.lock', 'a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        call('run', 3, fresh, *args, '--prefixes', '1,2,3', ok=False)
        call('verify', 3, fresh, '--prefixes', '1,2,3', ok=False)
    call('run', 3, base / 'memory', *args, '--prefixes', '2',
         '--worker-memory-bytes', str(10**18), ok=False)
    assert not list((base / 'memory').glob('prefixes/*/receipt.json'))
    call('run', 3, fresh, *args, '--prefixes', '1,1', ok=False)
    call('run', 3, fresh, '--threads', '2', '--cpus', f'{cpus[0]},{cpus[0]}',
         '--prefixes', '1', ok=False)

    # Native failures with valid compressed output must never publish a receipt.
    for fail in ('exit', 'unsupported'):
        script = base / f'native-{fail}'
        script.write_text('#!/usr/bin/env python3\nimport subprocess,sys\n'
                          f'p=subprocess.run([{str(BINARY)!r}, "invariants"], '
                          'input=b"3 0\\n" if ' + repr(fail == 'unsupported') + ' else sys.stdin.buffer.read(), '
                          'stdout=sys.stdout.buffer)\nsys.exit(1 if ' + repr(fail == 'exit') + ' else p.returncode)\n')
        script.chmod(0o755)
        failed = base / f'failed-{fail}'
        call('run', 3, failed, *args, '--prefixes', '2', '--native', script, ok=False)
        assert not list(failed.glob('prefixes/*/receipt.json'))

    f2 = base / 'f2'
    call('run', 2, f2, *args, '--prefixes', '1,2')
    assert len(rows(f2, 2)) == 4208
    call('verify', 2, f2, '--prefixes', '1,2', '--expected-count', '4208')
    reference = base / 'reference.m'
    out = base / 'reference.csv'
    listed = base / 'listed.bin'
    reference.write_text('SetQuitOnError(true);\nSetNthreads(1);\nSetColumns(0);\n'
                         'AttachSpec("spec");\n'
                         f'Write("{out}", FibrationInvariantCSV());\n'
                         'for q in [2, 3] do\nR<s,t> := PolynomialRing(GF(q), 2);\n'
                         'keys := FibrationRepresentativeKeys(R, 2 : GenericFibreKinds := ["smooth"], '
                         'Minimalities := ["everywhere minimal"]);\n'
                         f'F := Open("{listed}" cat Sprint(q), "wb");\n'
                         'WriteBytes(F, &cat[[(k div 256^j) mod 256 : j in [0..7]] : k in keys]);\n'
                         'delete F;\nfor k in keys[1..4] do\n'
                         f'Write("{out}", FibrationInvariantCSV(FibrationInvariantRow(k, aInvariantsFromKey(R,k))));\n'
                         'end for;\nend for;\nquit;\n')
    check = subprocess.run([os.environ.get('MAGMA', 'magma'), '-bn', reference],
                           cwd=ROOT, capture_output=True, text=True)
    assert check.returncode == 0, check.stdout + check.stderr
    expected = list(csv.DictReader(out.open()))
    for q, directory in [(2, f2), (3, fresh)]:
        assert (Path(str(listed) + str(q)).read_bytes() ==
                (directory / 'prefixes/0002/keys_0002.bin').read_bytes())
        assert rows(directory, 2)[:4] == [r for r in expected if r['q'] == str(q)]
    pilot = base / 'stage1-pilot'
    call('pilot', 3, pilot, *args)
    measurement = json.loads((pilot / 'pilot.json').read_text())
    assert measurement['passed'] and measurement['worker_peak_rss_bytes'] > 0
    probe = measurement['memory_probe']
    assert probe['prefix']['count'] == 1594323
    assert probe['key_bytes'] == 1594323 * 8
    assert probe['worker_peak_rss_bytes'] > 0
    assert measurement['worker_peak_rss_bytes'] >= probe['worker_peak_rss_bytes']
    controls = list(csv.DictReader((pilot / 'controls/reference.csv').open()))
    assert any(r['torsion'] not in ('', '[]') for r in controls)
    assert any(r['sha_reg'] == '3/4' for r in controls)
    inflated = base / 'inflated-pilot.json'
    measurement['worker_peak_rss_bytes'] = 10**18
    inflated.write_text(json.dumps(measurement))
    call('run', 3, base / 'measured-memory', *args, '--prefixes', '2',
         '--pilot-data', inflated, ok=False)
    call('run', 3, base / 'missing-pilot', *args,
         '--pilot-data', base / 'absent.json', ok=False)
    call('pilot', 2, base / 'pilot-f2', *args)
    print('STAGE1_OK: both fields, exact streams, recovery, failures, resources and pilot')


if __name__ == '__main__':
    with tempfile.TemporaryDirectory(prefix='stage1-', dir=ROOT / 'artifacts') as directory:
        base = Path(directory)
        for test in (test_revalidation_stop, test_changed_full_count, test_worker_error,
                     test_failed_logs, test_publication_recovery, test_magma_version, test_magma_spec,
                     test_legacy_pilot_budget, test_repilot, test_commands):
            output = base / test.__name__
            output.mkdir()
            test(output)
