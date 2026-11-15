from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'build/python'))
sys.path.append(str(ROOT))
import pyarrow as pa
import pyarrow.parquet as pq
from pipeline import packed
pa.set_cpu_count(1)
pa.set_io_thread_count(1)

COUNTS = {2: 1397598, 3: 1307723786}
META = {'schema': 'census-v2', 'package_commit': '0bdf62a', 'key_encoding': '1',
        'chart': 'homogeneous (s,t) evaluated at (t,1); inf is s/t infinity',
        'rational_encoding': 'reduced numerator/denominator strings; integers omit /1',
        'rho_F2': 'legacy column name; Picard number over the base field q',
        'fibres': 'place:degree:Kodaira:Tamagawa:discriminant-valuation:conductor:trace',
        'd2': 'dimension over F_2 of Pic(C)[2] for the cubic 2-division cover; null unless computed',
        'certificate': 'exhaustive bounded search and covering radius; replay with the shipped source',
        'polynomial_codes': 'base-q coefficients, constant term first',
        'generator_points': '[x_num,x_den,y_num,y_den]; identity []; codes are decimal strings',
        'gram_upper': 'row-major upper triangle; exact rational strings'}


def digest(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()


def schema(q, kind):
    seq = lambda ty: pa.list_(pa.field('element', ty))
    u8, s, nums = pa.uint8(), pa.string(), seq(pa.int64())
    key = pa.uint32() if q == 2 else pa.uint64()
    if kind == 'lean':
        fields = [('key', key), ('D', u8), ('eps', pa.int8()), ('c1_to_cm', nums),
                  ('torsion', seq(u8)), ('fibres', s), ('sha_reg', s), ('d2', pa.int16()),
                  ('d2_status', s), ('d2_method', s), ('admissible_sha', nums), ('generator_status', s)]
    elif kind == 'generators':
        fields = [('key', key), ('basis', seq(seq(s))), ('gram_upper', seq(s)),
                  ('torsion_points', seq(seq(s))), ('bound', u8), ('certificate_verdict', s)]
    else:
        fields = [('key', key), ('a1', s), ('mw_torsion', s), ('singular_fibres', s),
                  ('L_function', nums), ('L_sign', pa.int8()), ('analytic_rank', u8), ('P2', nums),
                  ('zeta_sign', pa.int8()), ('rho_F2', u8), ('rho_geom', u8), ('artin_height', s),
                  ('br_disc', pa.int64())] + [(f'pts_F{q**i}', pa.int64()) for i in range(1, 6)]
        fields += [('places', s), ('sha_reg', s), ('d2', pa.int16()), ('d2_status', s), ('d2_method', s),
                   ('admissible_sha', nums), ('admissible_reg', seq(s)), ('sha_order', s),
                   ('regulator', s), ('split_status', s), ('generator_status', s)]
    return pa.schema(fields, metadata={k.encode(): str(v).encode() for k, v in dict(META, q=q, table=kind).items()})


def a1_label(q, key):
    n = packed(q, key)[0]
    terms = []
    for i in range(2, -1, -1):
        c = n // q**i % q
        if c:
            factors = (['s' if i == 1 else f's^{i}'] if i else [])
            j = 2-i
            factors += (['t' if j == 1 else f't^{j}'] if j else [])
            if c != 1: factors.insert(0, str(c))
            terms.append('*'.join(factors))
    return ' + '.join(terms) or '0'


def fibre_string(local):
    return ';'.join(':'.join(str(l[k]) for k in ('place', 'degree', 'kodaira', 'tamagawa', 'vd', 'conductor', 'trace')) for l in local)


def local_from_string(text):
    local = []
    for part in text.split(';') if text else []:
        p, d, k, c, v, f, t = part.split(':')
        local.append(dict(place=p, degree=int(d), kodaira=k, tamagawa=int(c), vd=int(v), conductor=int(f), trace=int(t)))
    return local


def full_row(r):
    q, k, d, g, sp = r['q'], r['key'], r['derived'], r['generators'], r['split']
    local = sorted(r['local'], key=lambda l: (l['degree'], l['kodaira'], l['conductor'], l['tamagawa']))
    result = dict(key=k, a1=a1_label(q, k), mw_torsion=' x '.join('Z/'+str(n) for n in g['torsion']) or 'trivial',
                  singular_fibres=';'.join(f'{l["degree"]}:{l["conductor"]}:{l["tamagawa"]}:{l["kodaira"]}' for l in local),
                  L_function=r['L'], L_sign=r['sign'], analytic_rank=d['rank'], P2=d['P2'],
                  zeta_sign=d['zeta_sign'], rho_F2=d['rho'], rho_geom=d['rho_geom'],
                  artin_height=str(d['height']) if d['height'] else 'inf', br_disc=d['br_disc'])
    result.update({f'pts_F{q**i}': v for i, v in enumerate(d['points'], 1)})
    result.update(places=fibre_string(r['local']), sha_reg=d['sha_reg'], d2=r['d2'].get('d2'),
                  d2_status=r['d2']['status'], d2_method=r['d2'].get('method'),
                  admissible_sha=[v[0] for v in sp['candidates']], admissible_reg=[v[1] for v in sp['candidates']],
                  sha_order=str(sp['sha']) if sp['proved'] else None, regulator=sp['regulator'],
                  split_status='proved' if sp['proved'] else 'unresolved', generator_status=g['verdict'])
    return result


def lean_row(r):
    D = len(r['L'])-1
    return dict(key=r['key'], D=D, eps=r['sign'], c1_to_cm=r['L'][1:D//2+1],
                torsion=r['generators']['torsion'], fibres=fibre_string(r['local']), sha_reg=r['derived']['sha_reg'],
                d2=r['d2'].get('d2'), d2_status=r['d2']['status'], d2_method=r['d2'].get('method'),
                admissible_sha=[v[0] for v in r['split']['candidates']], generator_status=r['generators']['verdict'])


def generator_row(r):
    g = r['generators']
    if not g['certified']: return None
    rank = len(g['basis'])
    return dict(key=r['key'], basis=g['basis'], gram_upper=[g['gram'][i][j] for i in range(rank) for j in range(i, rank)],
                torsion_points=g['torsion_points'], bound=g['bound'], certificate_verdict=g['verdict'])


def expand(q, rows, scalar=False):
    data, partial = [], []
    types = {'II': 1, 'III': 2, 'IV': 3, 'IV*': 4, 'III*': 5, 'II*': 6}
    for r in rows:
        D, eps = r['D'], r['eps']
        if eps not in (-1, 1) or len(r['c1_to_cm']) != D//2: raise ValueError('invalid half-L')
        L = [1] + r['c1_to_cm'] + [0] * (D-D//2)
        for i in range(D//2+1):
            c = eps*q**(D-2*i)*L[i]
            if i == D-i and L[i] != c: raise ValueError('inconsistent central L coefficient')
            L[D-i] = c
        local = local_from_string(r['fibres'])
        ff = []
        for l in local:
            ks = l['kodaira']; code = types.get(ks, 8 if ks.endswith('*') else 7)
            n = int(ks[1:-1] if code == 8 else ks[1:]) if code >= 7 else 0
            ff.extend([l['degree'], code, n, l['tamagawa']])
        tor = 1
        for n in r['torsion']: tor *= n
        data.append([q, r['key'], 1, tor, 1, 0, len(L), *L, len(local), *ff])
        partial.append(dict(q=q, key=r['key'], L=L, sign=eps, local=local,
                            generators=dict(torsion=r['torsion'], verdict=r['generator_status']),
                            d2=dict(d2=r['d2'], status=r['d2_status'], method=r['d2_method'])))
    if not rows: return []
    mode = 'scalar' if scalar else 'native'
    proc = subprocess.run([str(ROOT/'bin'/('derive-'+mode)), '/dev/stdin', '0'],
                          input=''.join(' '.join(map(str, d))+'\n' for d in data), text=True, capture_output=True)
    if proc.returncode: raise RuntimeError('lean reconstruction failed: '+proc.stderr)
    derived = [json.loads(s) for s in proc.stdout.splitlines()]
    if len(derived) != len(rows): raise RuntimeError('lean reader count mismatch')
    result = []
    for source, r, d in zip(rows, partial, derived):
        if d['key'] != r['key'] or d['sha_reg'] != source['sha_reg']: raise ValueError('lean BSD product mismatch')
        sha = source['admissible_sha']
        if not sha or sha != sorted(set(sha)) or any(s < 1 for s in sha): raise ValueError('invalid Sha candidates')
        pairs = [[s, str(F(source['sha_reg'])/s)] for s in sha]
        r['derived'] = d
        r['split'] = dict(candidates=pairs, proved=len(pairs)==1, sha=sha[0] if len(sha)==1 else None,
                          regulator=pairs[0][1] if len(sha)==1 else None)
        result.append(full_row(r))
    return result


def write_parquet(path, rows, q, kind):
    path = Path(path)
    table = pa.Table.from_pylist(rows, schema=schema(q, kind))
    temp = path.with_suffix('.partial')
    pq.write_table(table, temp, compression='zstd', compression_level=9, row_group_size=4096, write_statistics=True)
    if not pq.read_table(temp).equals(table): raise RuntimeError('Parquet round trip mismatch')
    if temp.stat().st_size > 64000000: raise RuntimeError('Parquet shard exceeds 64 MB; lower --chunk-size')
    with temp.open('rb') as f: os.fsync(f.fileno())
    temp.replace(path)
    return dict(file=path.name, bytes=path.stat().st_size, sha256=digest(path), count=len(rows),
                min_key=rows[0]['key'] if rows else None, max_key=rows[-1]['key'] if rows else None)


def publish_shard(q, rows, dest, label, scalar=False):
    dest = Path(dest); dest.mkdir(parents=True, exist_ok=True)
    rows = sorted(rows, key=lambda r: r['key'])
    if len({r['key'] for r in rows}) != len(rows): raise ValueError('duplicate table key')
    full = [full_row(r) for r in rows]
    fullpath = dest/(label+'.surfaces.parquet')
    fm = write_parquet(fullpath, full, q, 'surfaces')
    gens = [g for r in rows if (g := generator_row(r)) is not None]
    gm = write_parquet(dest/(label+'.generators.parquet'), gens, q, 'generators')
    if q == 3:
        lean = [lean_row(r) for r in rows]
        lm = write_parquet(dest/(label+'.lean.parquet'), lean, q, 'lean')
        if expand(q, pq.read_table(dest/lm['file']).to_pylist(), scalar) != pq.read_table(fullpath).to_pylist():
            raise RuntimeError('full/lean exact reconstruction mismatch')
        # Both durable tables and their exact round trips precede deletion.
        fullpath.unlink()
        result = dict(lean=lm, generators=gm)
    else:
        result = dict(surfaces=fm, generators=gm)
    return result


def read_fragment(folder, manifest):
    q = manifest['q']; tables = manifest['tables']; kind = 'lean' if q == 3 else 'surfaces'
    for m in tables.values():
        path = Path(folder)/m['file']
        if digest(path) != m['sha256']: raise RuntimeError('table checksum mismatch: '+str(path))
    records = pq.read_table(Path(folder)/tables[kind]['file']).to_pylist()
    if len(records) != manifest['count'] or any(a['key'] >= b['key'] for a,b in zip(records,records[1:])):
        raise RuntimeError('table key order/count mismatch')
    gens = pq.read_table(Path(folder)/tables['generators']['file']).to_pylist()
    proved = {r['key'] for r in records if r['generator_status']=='proved'}
    if {r['key'] for r in gens} != proved or len(gens) != len(proved): raise RuntimeError('certified generator key mismatch')
    return records, gens
