import csv
import fcntl
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import shlex
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests/compare"))
from records import decompress

CLI = ROOT / "production/stage1.py"
BINARY = ROOT / "build/native/census-native"


def call(command, q, output, *extra, ok=True):
    result = subprocess.run([sys.executable, CLI, command, '--q', str(q),
                             '--output', output, *map(str, extra)],
                            cwd=ROOT, capture_output=True, text=True)
    assert (result.returncode == 0) == ok, result.stdout + result.stderr
    return result


def run_prefixes(q, output, indices, *extra):
    call('prepare', q, output, '--prefixes', indices, *extra)
    for index in indices.split(','):
        call('run', q, output, '--prefixes', index, *extra)


def compress(path, text, checksum=True):
    with path.open('wb') as output:
        subprocess.run(['zstd', '-q', '--check' if checksum else '--no-check'],
                       input=text, stdout=output, check=True)


def rows(output, prefix):
    return list(csv.DictReader(io.StringIO(decompress(
        (output / 'prefixes' / f'{prefix:04d}' / 'rows.csv.zst').read_bytes()))))


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
    output = base / 'prepare'
    result = call('prepare', 3, output, '--magma', wrapper, ok=False)
    assert 'full census count differs from pinned expected count' in result.stderr, result.stderr
    assert not marker.exists(), 'enumeration started with a changed full count'
    assert not (output / 'manifest.json').exists()


def test_failed_count_logs(base):
    wrapper = base / 'fail-count-magma'
    magma = shutil.which(os.environ.get('MAGMA', 'magma'))
    wrapper.write_text('#!/usr/bin/env python3\nimport subprocess,sys\n'
                       f'if "-V" in sys.argv: sys.exit(subprocess.run([{magma!r}, "-V"]).returncode)\n'
                       'print("count stdout diagnostic")\n'
                       'print("count stderr diagnostic", file=sys.stderr)\n'
                       'sys.exit(7)\n')
    wrapper.chmod(0o755)
    for command in ('prepare', 'verify'):
        output = base / command
        result = call(command, 3, output, '--prefixes', '1', '--magma', wrapper, ok=False)
        assert 'producer failed (7); see ' in result.stderr, result.stderr
        log = Path(result.stderr.rsplit('; see ', 1)[1].strip())
        assert log.is_relative_to(output), result.stderr
        assert log.exists(), f'failed count log was deleted: {log}'
        assert 'count stderr diagnostic' in log.read_text()
        assert 'count stdout diagnostic' in log.with_suffix('.stdout.log').read_text()


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
    call('prepare', 3, output, *args)
    call('run', 3, output, *args, ok=False)
    original = (output / 'scratch/prefix-0002/native.log').read_bytes()
    assert b'retained failure evidence' in original
    call('run', 3, output, *args)
    assert any(p.read_bytes() == original for p in output.rglob('native.log')), 'failed log was deleted'
    call('verify', 3, output, *args)


def test_publication_recovery(base):
    output = base / 'published'
    run_prefixes(3, output, '1,2')
    path = output / 'prefixes/0002/receipt.json'
    original = json.loads(path.read_text())
    log = path.with_name('native.log')
    with log.open('ab') as stream:
        stream.write(b'published attempt diagnostic\n')
    original_log = log.read_bytes()
    path.unlink()
    result = call('verify', 3, output, '--prefixes', '1,2', ok=False)
    assert 'missing completed prefix receipt' in result.stderr, result.stderr
    assert 'revalidating prefix' not in result.stderr, result.stderr
    run_prefixes(3, output, '1,2')
    assert path.exists(), 'resume skipped a published directory without its receipt'
    assert json.loads(path.read_text())['files'] == original['files']
    assert any(p.read_bytes() == original_log for p in (output / 'failed').glob('*/native.log')), \
        'resume deleted the unreceipted publication log'
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
    run_prefixes(3, output, '1', '--magma', wrapper)
    version.write_text('review-version-two')
    for command in ('prepare', 'run', 'verify'):
        result = call(command, 3, output, '--prefixes', '1', '--magma', wrapper, ok=False)
        assert 'incompatible manifest' in result.stderr


def test_magma_spec(base):
    shutil.copy2(ROOT / 'spec', base / 'spec')
    for directory in ('magma', 'production'):
        shutil.copytree(ROOT / directory, base / directory)
    command = [sys.executable, base / 'production/stage1.py', 'prepare', '--q', '3',
               '--output', base / 'results', '--prefixes', '1', '--native', BINARY]
    result = subprocess.run(command, cwd=base, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    spec = base / 'magma/spec'
    spec.write_text(spec.read_text() + '\n')
    result = subprocess.run(command, cwd=base, capture_output=True, text=True)
    assert result.returncode != 0, result.stdout + result.stderr
    assert 'incompatible manifest' in result.stderr, result.stdout + result.stderr


def test_commands(base):
    fresh = base / 'fresh'
    run_prefixes(3, fresh, '1,2,3')
    result = call('verify', 3, fresh, '--prefixes', '1,2,3')
    summary = json.loads(result.stdout)
    assert summary['partial'] and summary['count'] == 4494
    assert not any('mass' in gap.lower() for gap in summary['gaps']), summary['gaps']
    assert [len(rows(fresh, i)) for i in (1, 2, 3)] == [0, 1122, 3372]
    call('verify', 3, fresh, '--prefixes', '1,2,3', '--expected-count', '4494')
    assert (fresh / 'verification.json').exists()
    call('verify', 3, fresh, '--expected-count', '4494', ok=False)
    assert not (fresh / 'verification.json').exists(), 'failed verify retained a successful result'
    call('verify', 3, fresh, '--prefixes', '1,2,3', '--expected-count', '1307723786', ok=False)

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
        call('run', 3, damaged, '--prefixes', '2', ok=False)
        assert json.loads(path.read_text()) == data

    with open(fresh / 'locks/0002.lock', 'a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        call('run', 3, fresh, '--prefixes', '2', ok=False)
        call('run', 3, fresh, '--prefixes', '1')
    for selection in ('1,2', '1,1'):
        call('run', 3, fresh, '--prefixes', selection, ok=False)
    call('run', 3, fresh, ok=False)
    call('run', 3, base / 'unprepared', '--prefixes', '2', ok=False)

    # Native failures with valid compressed output must never publish a receipt.
    for fail in ('exit', 'unsupported'):
        script = base / f'native-{fail}'
        script.write_text('#!/usr/bin/env python3\nimport subprocess,sys\n'
                          f'p=subprocess.run([{str(BINARY)!r}, "invariants"], '
                          'input=b"3 0\\n" if ' + repr(fail == 'unsupported') + ' else sys.stdin.buffer.read(), '
                          'stdout=sys.stdout.buffer)\nsys.exit(1 if ' + repr(fail == 'exit') + ' else p.returncode)\n')
        script.chmod(0o755)
        failed = base / f'failed-{fail}'
        call('prepare', 3, failed, '--prefixes', '2', '--native', script)
        call('run', 3, failed, '--prefixes', '2', '--native', script, ok=False)
        assert not list(failed.glob('prefixes/*/receipt.json'))

    f2 = base / 'f2'
    run_prefixes(2, f2, '1,2')
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


def test_parallel(base):
    output = base / 'parallel'
    marker = base / 'failed-once'
    wrapper = base / 'controlled-native'
    wrapper.write_text('#!/usr/bin/env python3\nimport os,pathlib,subprocess,sys\n'
                       'scratch=os.environ["MAGMA_TEMP_DIR"]\n'
                       f'marker=pathlib.Path({str(marker)!r})\n'
                       'if scratch.endswith("prefix-0003") and not marker.exists():\n'
                       ' marker.touch(); sys.exit(7)\n'
                       f'sys.exit(subprocess.run([{str(BINARY)!r}, "invariants"]).returncode)\n')
    wrapper.chmod(0o755)
    call('prepare', 3, output, '--prefixes', '1,2,3', '--native', wrapper)
    command = shlex.join([sys.executable, str(CLI), 'run', '--q', '3', '--output', str(output),
                          '--native', str(wrapper), '--prefixes']) + ' {}'
    options = ['parallel', '-j2', '--halt', 'soon,fail=1', '--resume-failed',
               '--joblog', str(output / 'jobs.tsv'), '--results', str(output / 'logs'), command]
    result = subprocess.run(options, input='1\n2\n3\n', text=True, capture_output=True)
    assert result.returncode != 0, result.stdout + result.stderr
    assert (output / 'prefixes/0001/receipt.json').exists()
    retained = (output / 'prefixes/0001/receipt.json').read_bytes()
    result = subprocess.run(options, input='1\n2\n3\n', text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (output / 'prefixes/0001/receipt.json').read_bytes() == retained
    call('verify', 3, output, '--prefixes', '1,2,3', '--native', wrapper, '--expected-count', '4494')
    serial = base.parent / 'test_commands' / 'fresh'
    assert all(rows(output, i) == rows(serial, i) for i in (1, 2, 3))
    path = output / 'prefixes/0002/rows.csv.zst'
    path.write_bytes(path.read_bytes()[:-4])
    result = subprocess.run(options, input='1\n2\n3\n', text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    call('verify', 3, output, '--prefixes', '1,2,3', '--native', wrapper, ok=False)
    wrapper.write_text(wrapper.read_text() + '\n')
    call('prepare', 3, output, '--prefixes', '1,2,3', '--native', wrapper, ok=False)


if __name__ == '__main__':
    with tempfile.TemporaryDirectory(prefix='stage1-', dir=ROOT / 'build') as directory:
        base = Path(directory)
        for test in (test_failed_count_logs, test_changed_full_count,
                     test_failed_logs, test_publication_recovery, test_magma_version, test_magma_spec,
                     test_commands, test_parallel):
            output = base / test.__name__
            output.mkdir()
            test(output)
