# Sha and regulator candidates

`FibrationCandidates(E)` returns every pair allowed by the necessary constraints
below for an everywhere minimal elliptic K3 polynomial model over GF(2)(t)
or GF(3)(t). Attach `spec` first.

```c++
F<t> := FunctionField(GF(2));
E := EllipticCurve([F | 0, 0, 1, t^7, 0]);
sha, regulators := FibrationCandidates(E);
assert sha eq [1, 4];
assert regulators eq [Rationals() | 14, 7/2];
```

The first list contains distinct positive square integers in increasing order.
The second contains the corresponding positive rational regulators.
Their products equal `ShaRegulatorProduct` exactly. A singleton proves the
split, conditional on the supplied complete arithmetic. An ambiguous list
retains every survivor; it does not assert that every survivor occurs.
In rank zero the regulator is 1 and Sha equals `ShaRegulatorProduct`.

`FibrationCandidates(data, fibres)` accepts a
[FibrationInvariants record](native-invariants.md) and complete
[FibrationLocalData records](native-local.md) for the same elliptic K3.
This overload avoids repeating those computations. It reads `Rank`,
`PicardNumber`, `ArtinTateProduct` and `ShaRegulatorProduct`, and the local
fields `Type`, `Index`, `Tamagawa` and `Degree`. Correct supplied arithmetic
and the K3 hypothesis are caller preconditions. Basic ranges and positivity
are checked; inconsistent supplied data can give two empty lists.

## Necessary constraints

Write S = #Sha * Reg, N = #Br * abs(disc NS), r = rational MW rank and
rho = base-field Picard number. The BSD and Artin-Tate normalizations are
those of [FibrationInvariants](native-invariants.md#exact-formulas).
Here #Br = #Sha, and the finite Brauer group of a K3 over a finite field
has square order, including in characteristic 2. Thus a Sha order s is a
positive square divisor of N, and D = N/s is an integer.
The even NS lattice gives these further conditions:

```text
rho odd:  D is even;
rho even: D = 0 mod 4, or (-1)^(rho/2+1)*D = 1 mod 4.
```

Square divisors are enumerated from the prime factorization of N. The code
does not scan all integers up to sqrt(N). All comparisons use integers or
rational numbers, including for inputs larger than 2^53.

For positive rank, the height of a section satisfies

```text
h(P) = 4 + 2*(P.O) - sum_v degree(v)*contr_v(P),  (P.O) >= 0 integral.
```

A rational section meets a component fixed by residue-field Frobenius.
The possible local corrections are:

| Fibre | Tamagawa number | Corrections |
| --- | --- | --- |
| I_n | n | i*(n-i)/n, 0 <= i < n |
| I_n, n even | 2, nonsplit | 0, n/4 |
| I_n, n odd | 1, nonsplit | 0 |
| I_n* | 4 | 0, 1, 1+n/4 |
| I_n* | 2 | 0, 1 |
| I_0* | 1 | 0 |
| III, III* | 2 | 0, 1/2; 0, 3/2 respectively |
| IV, IV* | 3 | 0, 2/3; 0, 4/3 respectively |
| IV, IV* | 1 | 0 |
| II, II* | 1 | 0 |

The split row also covers I_1 and I_2. Each closed-place correction is
multiplied by its residue degree, including at infinity. Cartesian sums
give a superset C of possible total corrections; correlations between
different fibres are not assumed. The table applies to wild fibres too.

In rank one, Reg = S/s must equal 4+2*n-c for some integer n >= 0 and
c in C. In higher rank, the least positive possible height is

```text
h_min = min_{c in C} (4 + 2*max(0, floor((c-4)/2)+1) - c).
```

Hermite's inequality gives Reg >= h_min^r / gamma_r^r. For ranks 2 through
5 the exact gamma_r^r values are 4/3, 2, 4 and 8. Above rank 5, the
classical bound gamma_r <= (4/3)^((r-1)/2) gives the conservative exact test

```text
Reg^2 * (4/3)^(r*(r-1)/2) >= h_min^(2*r).
```

This may retain more candidates than sharper Hermite constants would.
No found-section, d2, descent, quotient-lattice or saturation filter is used.
No geometric NS discriminant or integral lattice is inferred.

## Checks and provenance

Run `make check/native_candidates_magma.m`. The existing runner, `make check`
and CI matrix discover this test; no extra discovery entry is required.
The fixed sample contains eight equations in both characteristics, including
rank zero, nontrivial torsion, non-rational places, singleton splits and
genuine ambiguity. Expected lists come from the Python prototypes below.
The minimum-height Hermite filter removes additional higher-rank candidates.
Supplied arithmetic controls cover the correction table, parity, ranks
through 20 and square orders above 2^53. Every returned pair multiplies
to the exact supplied product.

Fixtures are equations and arithmetic controls in test source; no census
dataset is added. Prototype paths are relative to `artifacts/mw/`.
These sources are provenance, not runtime dependencies.

| Prototype | SHA-256 |
| --- | --- |
| `hard/code/analyze.py` | `5b10dd729a1134696557a4249158f1763dc12c525c2b2dcd8410f5a697800243` |
| `hard/code/lattice.py` | `19c2093361887e45419b6f4ff55acca4268dc7e4291dce8c235c17a0cc1b619e` |
| `sha-reg-split/analyze.py`, minimum-height filter | `76ff3bac68625483b7d14961fa1352022165a9e42398fa79f63148195bc102d4` |
