from fractions import Fraction


def require(condition, message):
    if not condition:
        raise ValueError(message)


def source_identity(row):
    require(row['q'] == 2, 'wrong field for ordinary F_2 Sha[2]')
    codes = row['codes']
    require(len(codes) == 5 and all(type(x) is int for x in codes), 'coefficient codes')
    a1, a2, a3, a4, a6 = codes
    require(a2 == 0 and all(0 <= x < 2**w for x, w in zip(
        [a1, a3, a4, a6], [3, 7, 9, 13])), 'coefficient code range')
    require(row['key'] == (a1 << 29) + (a3 << 22) + (a4 << 13) + a6, 'source key identity')
    rank = row['derived']['rank']
    require(type(rank) is int and rank >= 0, 'Mordell-Weil rank')
    product = Fraction(row['derived']['sha_reg'])
    pairs = row['split']['candidates']
    require(pairs and len({n for n, _ in pairs}) == len(pairs), 'empty or duplicate candidate list')
    require(all(type(n) is int and n > 0 and Fraction(reg) * n == product
                for n, reg in pairs), 'candidate Sha*Reg disagreement')
    if row['split']['proved']:
        require(pairs == [[row['split']['sha'], row['split']['regulator']]], 'conflicting exact split')
    return dict(q=2, key=str(row['key']), codes=list(map(str, codes)), rank=rank,
                sha_reg=str(product))


def refine_sha2(row, defect=None):
    identity = source_identity(row)
    original = row['split']['candidates']
    pairs = [[str(n), str(Fraction(reg))] for n, reg in original]
    method = 'ordinary-f2-selmer-pairing-v1'
    dimension = dict(status='not_computed', value=None, method=method,
                     reason='saved_descent_not_available')
    feasible = {}
    if not row['codes'][0]:
        require(defect is None, 'supersingular pairing is not applicable')
        dimension.update(status='not_applicable', reason='supersingular_generic_fibre')
    elif row['split']['descent']['status'] == 'ok':
        descent = row['split']['descent']
        f, v, r = descent['f'], descent['v'], identity['rank']
        require(all(type(x) is int and x >= 0 for x in (f, v)), 'Selmer dimensions')
        upper = f + v - r - 1
        require(upper >= 0 and upper == descent['sha2_upper'], 'saved Selmer upper bound')
        for c in range(0, min(v, upper) + 1, 2):
            if defect is not None and defect != c:
                continue
            d = upper - c
            survivors = [[str(n), str(Fraction(reg))] for n, reg in original
                         if (n % 2 == 1 if d == 0 else n % 2**d == 0)]
            if survivors:
                feasible[c] = (d, survivors)
        require(feasible, 'empty feasible defect intersection')
        dimensions = sorted({d for d, _ in feasible.values()})
        orders = {n for _, survivors in feasible.values() for n, _ in survivors}
        pairs = [p for p in pairs if p[0] in orders]
        dimension.update(status='exact' if len(dimensions) == 1 else 'bound',
                         value=dimensions[0] if len(dimensions) == 1 else None,
                         candidates=dimensions, reason='residue_pairing' if defect is not None
                         else 'feasible_alternating_defects',
                         f=f, v=v, rank=r, upper_bound=upper,
                         defect_status='exact' if defect is not None else 'bound',
                         defect_candidates=sorted(feasible))
    else:
        require(defect is None, 'pairing without successful saved descent')
    split = dict(status='exact' if len(pairs) == 1 else 'bound', candidates=pairs,
                 method=method, sha=pairs[0][0] if len(pairs) == 1 else None,
                 regulator=pairs[0][1] if len(pairs) == 1 else None)
    return dict(sha2=dimension, split=split)
