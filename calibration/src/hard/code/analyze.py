"""Lower bounds on Mordell-Weil generator heights from #Sha * Reg.

python3 -I code/analyze.py --q 2 --out f2/derived.csv --summary f2/summary.json f2/out/*.csv

Input rows (code/entry.m): key, tag, rank, geometric rank, C = #Sha*Reg, torsion
invariants, torsion order, N0 = #Br * |disc NS| (Artin-Tate), arithmetic Picard
number, degree of the discriminant, bad fibres place:deg:type:c:ordDisc:ordCond,
L-polynomial coefficients, CPU seconds.
"""
import argparse, json, math, sys, os
from fractions import Fraction as Fr
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lattice import fibre_data

HEADER = ('q,key,tag,rank,geometric_rank,sha_times_reg,torsion,torsion_order,'
          'br_times_disc,picard_number,reducible_fibres,admissible_sha,'
          'admissible_reg,reg_determined,min_reg,height_lower_bound,'
          'po_lower_bound,po_values,po0_section_possible')


def isqrt_exact(n):
    r = math.isqrt(n)
    return r if r * r == n else None


def parity_ok(D, rho):
    """An even lattice of rank rho: rho odd -> D even; rho even -> 4 | D or
    (-1)^(rho/2) det = (-1)^(rho/2+1) D = 1 mod 4 (2-adic Jordan decomposition)."""
    if rho % 2 == 1:
        return D % 2 == 0
    if D % 4 == 0:
        return True
    if D % 2 == 0:
        return False
    return ((-1) ** (rho // 2 + 1) * D) % 4 == 1


def corrections(fibres):
    """All values of the sum of corrections over rational components."""
    s = {Fr(0)}
    for f in fibres:
        s = {a + b for a in s for b in f[2]}
    return s


def analyze_row(q, f):
    key, tag = int(f[0]), int(f[1])
    if f[2] == 'ERROR':
        return None, ['error:' + ','.join(f[3:])]
    r, rg = int(f[2]), int(f[3])
    C = Fr(f[4])
    tors, T = f[5], int(f[6])
    N0, rho, degdisc = Fr(f[7]), int(f[8]), int(f[9])
    issues = []
    if degdisc != 24:
        issues.append('degdisc')
    fibres, red = [], []
    for item in f[10].split(';') if f[10] else []:
        place, e, typ, c, od, oc = item.split(':')
        e, c = int(e), int(c)
        det, k, corr, _ = fibre_data(typ, c, e)
        if k:
            fibres.append((det, k, corr))
            red.append(f'{typ}@{place}' + (f'/c{c}' if True else ''))
    detT = 1
    for d, _, _ in fibres:
        detT *= d
    if N0.denominator != 1:
        issues.append('N0-not-integer')
    N0 = int(N0)
    if C * detT / (T * T) != N0:
        issues.append('artin-tate')
    if rho != 2 + sum(k for _, k, _ in fibres) + r:
        issues.append('rho')
    # admissible #Sha = m^2: m^2 | N0 and the parity of disc NS(X) = N0/m^2
    adm = []
    m = 1
    while m * m <= N0:
        if N0 % (m * m) == 0 and parity_ok(N0 // (m * m), rho):
            adm.append(m)
        m += 1
    po_values = {}
    if r == 0:
        # C = #Sha
        if C.denominator != 1 or isqrt_exact(int(C)) is None:
            issues.append('rank0-sha-not-square')
        adm = [m for m in adm if m * m == C]
        if not adm:
            issues.append('rank0-sha-inadmissible')
    elif r == 1:
        corr = corrections(fibres)
        keep = []
        for m in adm:
            reg = C / (m * m)
            pos = sorted({int((reg + c - 4) / 2) for c in corr
                          if (reg + c - 4) >= 0 and ((reg + c - 4) / 2).denominator == 1})
            if pos:
                keep.append(m)
                po_values[m] = pos
        adm = keep
        # A non-torsion section kP+U with (P.O) = 0 needs k^2 Reg + c = 4.
        po0 = any(k * k * C / (m * m) + c == 4 for m in adm for c in corr
                  for k in range(1, math.isqrt(int(4 * m * m / C)) + 2))
    if r >= 1 and not adm:
        issues.append('no-admissible-sha')
    out = dict(q=q, key=key, tag=tag, rank=r, geometric_rank=rg, C=C, torsion=tors,
               T=T, N0=N0, rho=rho, red=red, adm=adm, issues=issues)
    if r >= 1 and adm:
        regs = [C / (m * m) for m in adm]
        out['regs'] = regs
        out['min_reg'] = min(regs)
        if r == 1:
            out['lb_h'] = float(out['min_reg'])
            out['lb_po'] = min(min(v) for v in po_values.values())
        else:
            lb = float(out['min_reg']) ** (1.0 / r)
            out['lb_h'] = lb
            out['lb_po'] = max(0, math.ceil((lb - 4) / 2 - 1e-12))
        out['po_values'] = po_values
        if r == 1:
            out['po0'] = int(po0)
    return out, issues


def row_csv(o):
    adm = o['adm']
    regs = o.get('regs', [])
    po = o.get('po_values', {})
    return ','.join(map(str, [
        o['q'], o['key'], o['tag'], o['rank'], o['geometric_rank'], o['C'], o['torsion'],
        o['T'], o['N0'], o['rho'], ';'.join(o['red']),
        ' '.join(str(m * m) for m in adm),
        ' '.join(str(x) for x in regs),
        int(o['rank'] >= 1 and len(adm) == 1),
        o.get('min_reg', ''),
        '' if 'lb_h' not in o else f"{o['lb_h']:.6g}",
        o.get('lb_po', ''),
        ' '.join(f"{m * m}:" + '|'.join(map(str, po[m])) for m in adm if m in po),
        o.get('po0', '')]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--q', type=int, required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--summary', required=True)
    ap.add_argument('--top', type=int, default=40)
    ap.add_argument('files', nargs='+')
    a = ap.parse_args()
    by_rank = Counter()
    det_by_rank = Counter()
    po_hist = defaultdict(Counter)
    lbh_hist = defaultdict(Counter)
    sha_amb = defaultdict(Counter)
    issues = Counter()
    issue_keys = defaultdict(list)
    top = []
    rank1_no_po0 = 0
    nrows = 0
    cpu = 0.0
    seen = set()
    with open(a.out, 'w') as out:
        out.write(HEADER + '\n')
        for fn in a.files:
            for line in open(fn):
                f = line.rstrip('\n').split(',')
                if f[0] in seen:
                    issues['duplicate'] += 1
                    continue
                seen.add(f[0])
                nrows += 1
                o, iss = analyze_row(a.q, f)
                for i in iss:
                    issues[i] += 1
                    if len(issue_keys[i]) < 20:
                        issue_keys[i].append(f[0])
                if o is None:
                    continue
                cpu += float(f[-1])
                out.write(row_csv(o) + '\n')
                r = o['rank']
                by_rank[r] += 1
                if r >= 1 and 'lb_po' in o:
                    det_by_rank[r] += len(o['adm']) == 1
                    po_hist[r][o['lb_po']] += 1
                    lbh_hist[r][math.floor(o['lb_h'])] += 1
                    sha_amb[r][len(o['adm'])] += 1
                    if r == 1 and not o['po0']:
                        rank1_no_po0 += 1
                    top.append((o['lb_h'], o['lb_po'], r, o['key'], str(o['C']), o['torsion'],
                                ' '.join(str(m * m) for m in o['adm']), ';'.join(o['red'])))
    pos = sum(v for k, v in by_rank.items() if k >= 1)
    hard = sum(c for r in po_hist for po, c in po_hist[r].items() if po >= 1)
    top_h = sorted(top, key=lambda t: (-t[0], -t[1]))[:a.top]
    top_po = sorted(top, key=lambda t: (-t[1], -t[0]))[:a.top]
    summary = dict(
        q=a.q, rows=nrows, cpu_seconds=cpu, by_rank=dict(sorted(by_rank.items())),
        reg_determined_by_rank=dict(sorted(det_by_rank.items())),
        positive_rank=pos, po_lower_bound_at_least_1=hard,
        fraction_po0_search_cannot_solve=(hard / pos if pos else None),
        rank1_without_any_po0_section=rank1_no_po0,
        search_bound_cannot_solve={b: sum(c for r in po_hist for po, c in po_hist[r].items() if po > b)
                                   for b in range(0, 11)},
        po_lower_bound_hist={r: dict(sorted(h.items())) for r, h in sorted(po_hist.items())},
        height_lower_bound_floor_hist={r: dict(sorted(h.items())) for r, h in sorted(lbh_hist.items())},
        number_of_admissible_sha_hist={r: dict(sorted(h.items())) for r, h in sorted(sha_amb.items())},
        issues=dict(issues), issue_examples=dict(issue_keys),
        top_by_height=top_h, top_by_po=top_po)
    json.dump(summary, open(a.summary, 'w'), indent=1, default=str)
    print(json.dumps({k: summary[k] for k in ['rows', 'by_rank', 'reg_determined_by_rank',
          'positive_rank', 'po_lower_bound_at_least_1', 'fraction_po0_search_cannot_solve',
          'issues']}, default=str))


if __name__ == '__main__':
    main()
