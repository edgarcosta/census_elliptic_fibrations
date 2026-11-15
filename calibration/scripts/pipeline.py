import csv
from fractions import Fraction as F
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src/hard/code'))
from analyze import analyze_row
from lattice import _det
spec = importlib.util.spec_from_file_location('split_policy', ROOT / 'src/sha-reg-split/analyze.py')
policy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(policy)

def packed(q, key):
    sizes = [3, 7, 9, 13] if q == 2 else [5, 9, 13]
    values = []
    for size in reversed(sizes):
        values.insert(0, key % q**size)
        key //= q**size
    if key:
        raise ValueError('key outside encoding range')
    return [values[0], 0, *values[1:]] if q == 2 else [0, values[0], 0, *values[1:]]

def coefficients(n, q):
    a = []
    while n:
        a.append(n % q)
        n //= q
    return a or [0]

def cpu():
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time() + r.ru_utime + r.ru_stime

def stable(value):
    if isinstance(value, dict):
        return {k: stable(v) for k, v in value.items() if k not in ('ms', 'search_ms', 'saturation_ms', 'filter_ms')}
    if isinstance(value, list):
        return [stable(v) for v in value]
    return value

class Pipeline:
    def __init__(self, work, scalar=False):
        self.work = Path(work)
        self.work.mkdir(parents=True, exist_ok=True)
        self.mode = 'scalar' if scalar else 'native'
        self.timings = {}

    def call(self, name, args=(), data=None, stage=None, allow_error=False):
        cmd = [str(ROOT / 'bin' / (name + '-' + self.mode)), *map(str, args)]
        return self.command(cmd, data, stage or name, allow_error)

    def command(self, cmd, data, stage, allow_error=False, timeout=None):
        start = cpu()
        env = os.environ.copy()
        env.pop('MAGMA_STARTUP_FILE', None)
        try:
            r = subprocess.run(cmd, env=env, input=data, capture_output=True, text=True, cwd=ROOT, timeout=timeout)
        except subprocess.TimeoutExpired as err:
            self.timings[stage] = self.timings.get(stage, 0) + cpu() - start
            raise RuntimeError(f'{stage} timed out after {timeout}s') from err
        self.timings[stage] = self.timings.get(stage, 0) + cpu() - start
        with (self.work / (stage + '.log')).open('a') as f:
            f.write(r.stderr)
            if r.returncode:
                f.write(r.stdout)
        if r.returncode and not allow_error:
            raise RuntimeError(f'{stage} exited {r.returncode}; see {self.work / (stage + ".log")}')
        return r

    def stages(self, rows, bundle=True, stage='generators'):
        if not rows:
            return []
        r = self.call('stages', ['--bundle'] if bundle else [], ''.join(json.dumps(x) + '\n' for x in rows), stage)
        out = [json.loads(s) for s in r.stdout.splitlines()]
        if [(x['q'], x['key']) for x in out] != [(x['q'], x['key']) for x in rows]:
            raise RuntimeError('stages output does not match the input keys')
        return out

    def fallback_d2(self, row):
        pars = [coefficients(row['codes'][i], 3) for i in (1, 3, 4)]
        path = self.work / 'fallback-input.m'
        path.write_text('<' + str(row['key']) + ',' + str(pars) + '>\n')
        r = self.command(['magma', '-bn', 'IN:=' + str(path), 'src/d2-native/fallback.m'], None, 'd2_fallback')
        lines = [s.split('|') for s in r.stdout.splitlines() if s.startswith('FALLBACK|')]
        if len(lines) != 1 or 'COMPLETE' not in r.stdout:
            (self.work / 'd2_fallback.log').write_text(r.stdout + r.stderr)
            raise RuntimeError('Magma d2 fallback did not complete')
        return dict(genus=int(lines[0][2]), d2=int(lines[0][3]), method='proved_class_group')

    def descent(self, row):
        codes = row['codes']
        script = f'''SetQuitOnError(true);SetSeed(20261114);SetNthreads(1);SetColumns(0);
K<t> := FunctionField(GF(2));
function Decode(n)
 p:=K!0;i:=0;
 while n gt 0 do p+:=(n mod 2)*t^i;n div:=2;i+:=1;end while;
 return p;
end function;
E:=EllipticCurve([Decode(n):n in {codes}]);
try
 SV,SF:=TwoIsogenySelmerGroups(E);
 printf "DESCENT|%o|%o\\n",Valuation(#SV,2),Valuation(#SF,2);
catch err
 printf "ERROR|%o\\n",err;
end try;
quit;
'''
        path = self.work / 'descent.m'
        path.write_text(script)
        try:
            r = self.command(['magma', '-bn', str(path)], None, 'descent', allow_error=True, timeout=25)
        except RuntimeError:
            return dict(status='timeout', sha2_upper=None)
        lines = [s.split('|') for s in r.stdout.splitlines() if s.startswith('DESCENT|')]
        if r.returncode or len(lines) != 1:
            with (self.work / 'descent.log').open('a') as f:
                f.write(str(row['key']) + '\n' + r.stdout)
            return dict(status='unresolved', sha2_upper=None)
        v, f = map(int, lines[0][1:])
        upper = v + f - row['rank'] - 1
        if upper < 0:
            raise RuntimeError('negative Sha[2] upper bound')
        return dict(status='ok', v=v, f=f, sha2_upper=upper)

    def run(self, q, keys):
        if not keys:
            return []
        start = cpu()
        rows = [dict(q=q, key=k, codes=packed(q, k), b=4 if q == 2 else 3) for k in keys]
        source = ''.join(' '.join(map(str, [q, r['key'], 0, *r['codes']])) + '\n' for r in rows)
        self.timings['decode'] = cpu() - start
        inp = self.work / 'input.txt'
        inp.write_text(source)
        locpath = self.work / 'local.csv'
        lr = self.call('lfun', ['--local-output', locpath], source)
        ls = list(csv.DictReader(io.StringIO(lr.stdout)))
        if len(ls) != len(rows):
            raise RuntimeError('L-function output count mismatch')
        local = {k: [] for k in keys}
        for l in csv.DictReader(locpath.open()):
            if int(l['q']) != q or int(l['key']) not in local:
                raise RuntimeError('unexpected local-data key')
            local[int(l['key'])].append({k: v if k in ('place', 'kodaira') else int(v) for k, v in l.items() if k not in ('q', 'key')})
        for row, l in zip(rows, ls):
            if (int(l['q']), int(l['key']), l['status'], l['other_candidate']) != (q, row['key'], 'ok', ''):
                raise RuntimeError('L-function failure: ' + str(l))
            row['L'] = list(map(int, l['L_polynomial'].split()))
            row['sign'] = int(l['sign'])
            row['local'] = sorted(local[row['key']], key=lambda x: (x['degree'], x['place']))
        natives = self.stages([{k: r[k] for k in ('q', 'key', 'b', 'codes', 'L')} for r in rows])
        derive = []
        types = {'II': 1, 'III': 2, 'IV': 3, 'IV*': 4, 'III*': 5, 'II*': 6}
        for row, n in zip(rows, natives):
            row['native'] = n
            row['rank'] = n['rank']
            row['sha_reg'] = n['sha_reg']
            row['torsion_order'] = math.prod(n['torsion'])
            ff = []
            for l in row['local']:
                ks = l['kodaira']
                code = types.get(ks, 8 if ks.endswith('*') else 7)
                number = int(ks[1:-1] if code == 8 else ks[1:]) if code >= 7 else 0
                ff.extend([l['degree'], code, number, l['tamagawa']])
            native_local = sorted(n['local'], key=lambda x: (x[1], x[0]))
            if native_local != [[l[k] for k in ('place', 'degree', 'kodaira', 'tamagawa', 'vd', 'conductor')] for l in row['local']]:
                raise RuntimeError('local engines disagree')
            derive.append([q, row['key'], 1, row['torsion_order'], 1, 0, len(row['L']), *row['L'], len(row['local']), *ff])
        dp = self.work / 'derive.txt'
        dp.write_text(''.join(' '.join(map(str, x)) + '\n' for x in derive))
        derived = [json.loads(x) for x in self.call('derive', [dp, 0]).stdout.splitlines()]
        extras = list(csv.reader(io.StringIO(self.call('extras', [q, inp]).stdout)))
        if len(derived) != len(rows) or len(extras) != len(rows):
            raise RuntimeError('derived invariant output count mismatch')
        start = cpu()
        for row, d, ex in zip(rows, derived, extras):
            if d['key'] != row['key'] or int(ex[0]) != row['key'] or F(d['sha_reg']) != F(row['sha_reg']) or d['rank'] != row['rank']:
                raise RuntimeError('derived invariants disagree')
            row['derived'] = {k: v for k, v in d.items() if k not in ('sha_candidate', 'disc_candidate', 'reg_candidate', 'split_known', 'q', 'key')}
            row['derived'].update(j_degree=int(ex[1]), j_insep_exponent=int(ex[2]))
            row['fibres'] = ';'.join(':'.join(str(l[k]) for k in ('place', 'degree', 'kodaira', 'tamagawa', 'vd', 'conductor')) for l in row['local'])
            ar, issues = analyze_row(q, [row['key'], 0, row['rank'], d['geom_rank'], row['sha_reg'],
                'x'.join(map(str, row['native']['torsion'])) or '1', row['torsion_order'], d['br_disc'], d['rho'],
                sum(l['degree'] * l['vd'] for l in row['local']), row['fibres']])
            if issues:
                raise RuntimeError(f'admissible-list inconsistency {row["key"]}: {issues}')
            row['sha'] = [m*m for m in ar['adm']]
            row['initial_sha'] = row['sha'][:]
            row['d2'] = dict(status='not_needed' if q == 3 else 'not_applicable')
        self.timings['admissible'] = cpu() - start
        todo = [r for r in rows if q == 3 and len(r['sha']) > 1]
        if todo:
            data = ''
            for r in todo:
                pars = [coefficients(r['codes'][i], 3) for i in (1, 3, 4)]
                data += ' '.join(map(str, [r['key'], *[v for a in pars for v in [len(a), *a]]])) + '\n'
            values = list(csv.DictReader(io.StringIO(self.call('d2', ['--stream'], data).stdout)))
            if len(values) != len(todo):
                raise RuntimeError('d2 output count mismatch')
            for row, v in zip(todo, values):
                if int(v['key']) != row['key']:
                    raise RuntimeError('d2 key mismatch')
                d = self.fallback_d2(row) if v['status'] == 'fallback' else dict(genus=int(v['genus']), d2=int(v['d2']), method='native')
                bound = sum(2 if l['kodaira'].startswith('I') and l['kodaira'].endswith('*') and l['kodaira'][1:-1].isdigit() and int(l['kodaira'][1:-1]) % 2 == 0 and l['tamagawa'] == 4 else 1 for l in row['local'] if l['tamagawa'] % 2 == 0)
                dims = [n for n in range(max(0, d['d2'] - bound - row['rank']), d['d2'] + bound - row['rank'] + 1) if n % 2 == 0]
                row['sha'] = [n for n in row['sha'] if any(n % 2 == 1 if dim == 0 else (n & -n).bit_length()-1 >= dim for dim in dims)]
                if not row['sha']:
                    raise RuntimeError('d2 eliminated every Sha candidate')
                row['d2'] = dict(status='ok', **d, component_bound=bound, sha2_dimensions=dims, sha=row['sha'][:])
        start = cpu()
        for row in rows:
            row['split_evidence'] = []
            row['descent'] = dict(status='not_needed')
            if row['rank'] > 0:
                f = policy.filters(row, row['native'])
                row['sha'] = f['combined_sha']
                row['split_evidence'].append({k: stable(v) for k, v in f.items() if k not in ('native', 'filter_ms', 'q', 'key')})
        self.timings['split_filters'] = cpu() - start
        for b in (range(5, 8) if q == 2 else [4]):
            todo = [r for r in rows if len(r['sha']) > 1 and (r['rank'] != 1 or min(F(r['sha_reg']) / n for n in r['sha']) <= 4 + 2*b)]
            ns = self.stages([dict(q=q, key=r['key'], b=b, codes=r['codes'], L=r['L']) for r in todo], bundle=False, stage='split_search')
            for row, n in zip(todo, ns):
                f = policy.filters(dict(row, b=b), n)
                row['sha'] = f['combined_sha']
                row['split_evidence'].append({k: stable(v) for k, v in f.items() if k not in ('filter_ms', 'q', 'key')})
        for row in rows:
            if q == 2 and len(row['sha']) > 1 and row['codes'][0] and any(n % 2 == 0 for n in row['sha']):
                row['descent'] = self.descent(row)
                upper = row['descent'].get('sha2_upper')
                if upper is not None and upper < 2:
                    row['sha'] = [n for n in row['sha'] if n % 2]
            if not row['sha']:
                raise RuntimeError('split eliminated every Sha candidate')
        out = []
        for row in rows:
            pairs = [[n, str(F(row['sha_reg']) / n)] for n in row['sha']]
            n = stable(row['native'])
            n.pop('q'); n.pop('key'); n.pop('local')
            n['verdict'] = 'proved' if n['certified'] else 'rank_short' if len(n['basis']) < n['rank'] else 'saturation_unresolved'
            out.append(dict(q=q, key=row['key'], codes=row['codes'], L=row['L'], sign=row['sign'], local=row['local'],
                derived=row['derived'], initial_sha=row['initial_sha'], d2=row['d2'],
                split=dict(proved=len(pairs) == 1, sha=pairs[0][0] if len(pairs) == 1 else None,
                    regulator=pairs[0][1] if len(pairs) == 1 else None, candidates=pairs,
                    evidence=row['split_evidence'], descent=row['descent']), generators=n))
        return out
