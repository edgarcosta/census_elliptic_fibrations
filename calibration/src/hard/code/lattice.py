"""Fibre lattices and section corrections for the generator-height bounds.

A reducible fibre is given by (Kodaira type, c, e): c the number of components
over the residue field (Magma's Tamagawa number), e the degree of the place.
Frobenius of the residue field acts on the components of one geometric fibre by a
diagram automorphism fixing the identity component, determined by (type, c).
"""
from fractions import Fraction as Fr
from functools import lru_cache
import re


def parse_type(s):
    """'I5' -> ('I', 5), 'I2*' -> ('I*', 2), 'IV*' -> ('IV*', None)."""
    m = re.fullmatch(r'I(\d+)(\*?)', s)
    if m:
        return ('I*' if m.group(2) else 'I', int(m.group(1)))
    assert s in ('II', 'III', 'IV', 'II*', 'III*', 'IV*'), s
    return (s, None)


def _root_lattice(kind, n):
    """Gram matrix (negative definite) of the non-identity components and the
    labels needed for the Frobenius action."""
    edges = []
    if kind == 'I':
        # cycle T0..T(n-1); non-identity T1..T(n-1), labelled 1..n-1 -> index i-1
        labels = list(range(1, n))
        edges = [(i, i + 1) for i in range(1, n - 1)]
    elif kind == 'I*':
        # T1 near leaf, chain C0..Cn, far leaves T2, T3; T0 (identity) on C0.
        labels = ['T1'] + [('C', j) for j in range(n + 1)] + ['T2', 'T3']
        edges = [('T1', ('C', 0))] + [(('C', j), ('C', j + 1)) for j in range(n)]
        edges += [('T2', ('C', n)), ('T3', ('C', n))]
    elif kind == 'III':
        labels = ['A']
    elif kind == 'IV':
        labels = ['T1', 'T2']
        edges = [('T1', 'T2')]
    elif kind == 'IV*':
        # E6: centre Z, arms Z-A1-(T0), Z-B1-B2, Z-D1-D2.
        labels = ['Z', 'A1', 'B1', 'B2', 'D1', 'D2']
        edges = [('Z', 'A1'), ('Z', 'B1'), ('B1', 'B2'), ('Z', 'D1'), ('D1', 'D2')]
    elif kind == 'III*':
        # E7: chain of 7 with a node on the 4th; T0 at one end of the long chain.
        labels = list(range(7))
        edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (2, 6)]
    elif kind == 'II*':
        labels = list(range(8))
        edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (4, 7)]
    else:
        labels = []
    idx = {l: i for i, l in enumerate(labels)}
    g = [[-2 if i == j else 0 for j in range(len(labels))] for i in range(len(labels))]
    for a, b in edges:
        g[idx[a]][idx[b]] += 1
        g[idx[b]][idx[a]] += 1
    return labels, idx, g


def _frobenius(kind, n, c, labels):
    """The permutation of the non-identity components as a dict label -> label."""
    perm = {l: l for l in labels}
    if kind == 'I':
        if c == n or n <= 2:
            assert c == n
        else:
            assert c == (2 if n % 2 == 0 else 1), (kind, n, c)
            perm = {i: n - i for i in labels}
    elif kind == 'I*':
        if c == 4:
            pass
        elif c == 2:
            perm['T2'], perm['T3'] = 'T3', 'T2'
        else:
            assert n == 0 and c == 1, (kind, n, c)
            perm['T1'], perm['T2'], perm['T3'] = 'T2', 'T3', 'T1'
    elif kind in ('IV', 'IV*'):
        if c == 1:
            if kind == 'IV':
                perm = {'T1': 'T2', 'T2': 'T1'}
            else:
                perm.update({'B1': 'D1', 'D1': 'B1', 'B2': 'D2', 'D2': 'B2'})
        else:
            assert c == 3
    elif kind in ('III', 'III*'):
        assert c == 2
    elif kind in ('II', 'II*'):
        assert c == 1
    return perm


def _orbits(perm):
    seen, orbs = set(), []
    for l in perm:
        if l in seen:
            continue
        o, x = [], l
        while x not in o:
            o.append(x)
            x = perm[x]
        seen.update(o)
        orbs.append(o)
    return orbs


def _det(m):
    m = [[Fr(x) for x in row] for row in m]
    n, d = len(m), Fr(1)
    for i in range(n):
        p = next((r for r in range(i, n) if m[r][i] != 0), None)
        if p is None:
            return Fr(0)
        if p != i:
            m[i], m[p] = m[p], m[i]
            d = -d
        d *= m[i][i]
        for r in range(i + 1, n):
            f = m[r][i] / m[i][i]
            for k in range(i, n):
                m[r][k] -= f * m[i][k]
    return d


@lru_cache(maxsize=None)
def fibre_data(typ, c, e):
    """(|det| of the Frobenius-invariant part of the fibre lattice, its rank,
    the set of correction values over rational components, scaled by e,
    and the product of the Frobenius orbit sizes on all geometric components)."""
    kind, n = parse_type(typ)
    labels, idx, g = _root_lattice(kind, n if n is not None else 0)
    perm = _frobenius(kind, n, c, labels)
    orbs = _orbits(perm)
    k = len(orbs)
    gram = [[sum(g[idx[a]][idx[b]] for a in oa for b in ob) for ob in orbs] for oa in orbs]
    det = abs(_det(gram)) * e ** k if k else Fr(1)
    orbit_product = _prod(e * len(o) for o in orbs)
    corr = {Fr(0)}
    if kind == 'I':
        rat = [i for i in range(n) if (n - i) % n == i % n] if c != n else list(range(n))
        corr = {Fr(i * (n - i), n) for i in rat}
    elif kind == 'I*':
        corr = {Fr(0), Fr(1)} if c in (2, 4) else {Fr(0)}
        if c == 4:
            corr.add(1 + Fr(n, 4))
    elif kind == 'III':
        corr = {Fr(0), Fr(1, 2)}
    elif kind == 'III*':
        corr = {Fr(0), Fr(3, 2)}
    elif kind == 'IV' and c == 3:
        corr = {Fr(0), Fr(2, 3)}
    elif kind == 'IV*' and c == 3:
        corr = {Fr(0), Fr(4, 3)}
    return det, k, frozenset(e * x for x in corr), orbit_product


def _prod(xs):
    r = 1
    for x in xs:
        r *= x
    return r
