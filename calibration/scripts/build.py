import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
p = argparse.ArgumentParser()
p.add_argument('--jobs', type=int, default=4)
p.add_argument('--scalar', action='store_true', help='omit native executables')
a = p.parse_args()
if a.jobs < 1:
    p.error('--jobs must be positive')
Path('build').mkdir(exist_ok=True)
Path('bin').mkdir(exist_ok=True)
os.environ['TMPDIR'] = str(ROOT / 'build')
cxx = os.environ.get('CXX', 'c++')
missing = []
if sys.version_info < (3, 10):
    missing.append('Python 3.10+: sudo apt-get install python3; or sudo dnf install python3')
if not shutil.which(cxx):
    missing.append(f'C++17 compiler ({cxx}): sudo apt-get install g++; or sudo dnf install gcc-c++')
if not shutil.which('magma'):
    missing.append('Magma: install your licensed Magma distribution and put magma on PATH (no distribution package)')
if not shutil.which('taskset'):
    missing.append('taskset: sudo apt-get install util-linux; or sudo dnf install util-linux')
checks = [
    ('GMP C++ development files', '#include <gmpxx.h>\nint main(){mpz_class n=3; return n==3?0:1;}', ['-lgmpxx', '-lgmp'], 'libgmp-dev', 'gmp-devel'),
    ('NTL development files', '#include <NTL/lzz_pXFactoring.h>\nint main(){NTL::zz_p::init(3);}', ['-lntl', '-lgmp'], 'libntl-dev', 'ntl-devel'),
    ('FLINT development files', '#include <flint/fmpz_mat.h>\nint main(){fmpz_mat_t m;fmpz_mat_init(m,1,1);fmpz_mat_clear(m);}', ['-lflint', '-lgmp'], 'libflint-dev', 'flint-devel'),
    ('Boost JSON headers (Boost 1.75+)', '#include <boost/json/src.hpp>\n#include <boost/property_tree/json_parser.hpp>\nint main(){boost::json::value v=1;}', [], 'libboost-json-dev', 'boost-devel'),
]
if shutil.which(cxx):
    for i, (name, source, libs, apt, dnf) in enumerate(checks):
        probe = Path(f'build/probe-{i}.cpp')
        probe.write_text(source + '\n')
        r = subprocess.run([cxx, '-std=c++17', str(probe), *libs, '-o', f'build/probe-{i}'], capture_output=True, text=True)
        Path(f'build/probe-{i}.log').write_text(r.stdout + r.stderr)
        if r.returncode:
            missing.append(f'{name}: sudo apt-get install {apt}; or sudo dnf install {dnf} (details: build/probe-{i}.log)')
try:
    sys.path.insert(0, str(ROOT / 'build/python'))
    import pyarrow
except ImportError:
    missing.append('PyArrow: python3 -m pip install --target build/python pyarrow')
if missing:
    print('\n'.join('Missing or unusable: ' + s for s in missing), file=sys.stderr)
    raise SystemExit(1)

common = ['-O3', '-std=c++17', '-Wall', '-Wextra']
avx2 = 'avx2' in Path('/proc/cpuinfo').read_text().split()
jobs = []
for mode in (['scalar'] if a.scalar else ['native']):
    flags = ['-fno-tree-vectorize'] if mode == 'scalar' else ['-march=native']
    for name, files, libs in [
        ('lfun', ['lfun-csign/' + s + '.cpp' for s in ['counter', 'tate', 'field', 'evaluate']], ['-lntl', '-lgmp']),
        ('stages', ['c-stages/core.cpp', 'c-stages/stages.cpp', 'lfun-tate/tate.cpp', 'search-c/search2.cpp', 'search-c/search3.cpp'], ['-lntl', '-lflint', '-lgmpxx', '-lgmp']),
        ('d2', ['d2-native/d2.cpp'], ['-lntl', '-lgmp']),
        ('derive', ['census-invariants/derive.cpp'], []),
        ('extras', ['census-invariants/extras.cpp'], ['-lntl', '-lgmp']),
    ]:
        extra = ['-DUSE_AVX2'] if name == 'lfun' and mode == 'native' and avx2 else []
        jobs.append((name + '-' + mode, [cxx, *common, *flags, *extra, *['src/' + f for f in files], *libs, '-o', f'bin/{name}-{mode}']))

def build(job):
    name, cmd = job
    print('Building ' + name, flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    Path('build', name + '.log').write_text(r.stdout + r.stderr)
    if r.returncode:
        raise RuntimeError(f'Build failed: build/{name}.log')
    return name, hashlib.sha256(Path('bin', name).read_bytes()).hexdigest()

with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
    hashes = dict(pool.map(build, jobs))
manifest = dict(compiler=subprocess.check_output([cxx, '--version'], text=True).splitlines()[0],
                machine=platform.machine(), flags=common, native=not a.scalar, avx2=avx2,
                binaries=hashes)
for mode in (['scalar'] if a.scalar else ['native']):
    result = subprocess.run([str(ROOT / 'bin' / ('lfun-'+mode)), '--certify-tables'], capture_output=True, text=True)
    Path('build', 'table-certificate-'+mode+'.txt').write_text(result.stdout+result.stderr)
    if result.returncode or result.stdout.count('TABLE_OK|') != 24:
        raise RuntimeError('Exact table certification failed: '+mode)
manifest['tables_certified'] = True
manifest['pyarrow'] = pyarrow.__version__
Path('build/manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print('Build complete. Run: python3 census.py --help')
