# Exact derived K3 invariants

`FibrationInvariants(E)` computes the invariants below for an everywhere
minimal elliptic K3 polynomial model over GF(2)(t) or GF(3)(t).
Attach `spec` from the repository root first. It obtains the signed
L-polynomial, bad-place records and complete torsion through
[`FibrationLPolynomial`](native-lfun.md),
[`FibrationLocalData`](native-local.md) and
[`FibrationTorsion`](native-torsion.md).

`FibrationInvariants(E, L, fibres, torsion)` performs the derivation from
already computed, complete data. This overload accepts any finite constant
field F_q. The caller supplies the correct L-polynomial, component cycles
and residue-field Tamagawa numbers for E, and its torsion invariant factors.
It checks the integer polynomial's constant coefficient, torsion factor
order and resulting cohomology degree. It does not certify the supplied
arithmetic or the K3 hypothesis. Fibre records need `ComponentCycles` and
`Tamagawa`; a cycle already includes conjugate fibres at a closed place.

## Result record

| Field | Meaning |
| --- | --- |
| `Rank`, `GeometricRank` | Rational and geometric Mordell-Weil ranks |
| `PicardNumber`, `GeometricPicardNumber` | Base-field and geometric Picard numbers |
| `P2` | det(1 - Frobenius*T on H^2), in Z[T], of degree 22 |
| `ZetaSign` | Leading coefficient of P2 divided by q^22, either -1 or 1 |
| `PointCounts` | Counts on the smooth surface over F_q through F_(q^5) |
| `BrauerHeight` | Formal Brauer group height 1 through 10, or `Infinity()` |
| `ArtinTateProduct` | #Br(X) times the absolute base-field NS discriminant |
| `NeronSeveriSplittingDegree` | Least constant-field extension degree over which every geometric NS class is defined |
| `JDegree` | Degree of the reduced j-map; zero for a constant map |
| `JInseparabilityExponent` | e such that the nonconstant j-map has inseparable degree p^e; zero for a constant map |
| `IsIsotrivial` | Whether j is constant |
| `ShaRegulatorProduct` | Exact rational product #Sha(E) * Reg(E) |

All arithmetic uses Magma integers and rationals. Infinite formal height
is distinct from every integer. The record determines neither a geometric
NS discriminant nor an integral NS lattice. The
[candidate interface](native-candidates.md) lists the admissible Sha orders
and regulators and identifies splits forced by these constraints.

## Exact formulas

Write the component cycle lengths as n_i. The rational rank r is the
multiplicity of 1-q*T in L. Geometric rank is the total degree, with
multiplicities, of cyclotomic factors in L(T/q). Shioda-Tate gives

```text
rho = 2 + r + number of component cycles,
rho_geom = 2 + geometric rank + sum n_i,
P2(T) = L(T) * (1-q*T)^2 * product_i (1-(q*T)^n_i).
```

The NS splitting degree is the LCM of the component cycle lengths and
the orders of the cyclotomic factors. An empty collection has LCM 1.
Newton identities recover the first five Frobenius traces from P2;
the count over F_(q^n) is 1 + q^(2*n) + trace(Frobenius^n on H^2).

For q=p^f and L(T)=sum c_i*T^i, the smallest Newton slope is

```text
lambda = min({1} union {v_p(c_i)/(f*i) : i > 0 and c_i != 0}).
```

The formal height is 1/(1-lambda) when lambda < 1, and `Infinity()` when
lambda = 1. Component factors have slope 1, so L suffices here.
The j-map is reduced by `jInvariant(E)`; repeated p-th roots until its
derivative is nonzero give the inseparability exponent.

Remove exactly the rational root multiplicities:

```text
Q = L/(1-q*T)^r,
G = P2/(1-q*T)^rho,
#Sha * Reg = q * Q(1/q) * (#torsion)^2 / product_v c_v,
#Br * abs(disc NS) = q * G(1/q).
```

These are the chi(O_X)=2 normalizations. Every closed place contributes
one Tamagawa factor c_v, including infinity. There is no residue-degree
exponent on c_v. The [BSD justification](native-mathematics.md#frobenius-twist-functional-equation-sign)
fixes the conventions. The Artin-Tate product is a positive integer;
the Sha-regulator product remains an exact rational number.

## Independent checks

Run `make check/native_derived_magma.m`. The existing runner, `make check`
and CI matrix discover this entry point. Each CI test has a 30-minute budget.

Eight fixed equations include rank zero with torsion of order 6 and
non-rational bad places. Both characteristics have nontrivial torsion and
supersingular examples. A torsion group with factors [2,4] has
Sha-regulator product 3/4. The finite formal heights are 1, 3, 4 and 8.
Independent Magma `LFunction`, `AnalyticInformation`,
`FrobeniusActionOnTrivialLattice`, `ArtinTateFormula`, Newton polygons and
the zeta power series supply comparisons. The tests also check both
Shioda-Tate identities and BSD/Artin-Tate consistency.

For the supersingular equation y^2+y=x^3+t^6*x over GF(2),
L=(1-2*T)*(1-(2*T)^6) and all 13 non-identity component classes are fixed.
Over GF(2^12), all 22 eigenvalues are q. The test supplies L=(1-q*T)^7
and checks P2=(1-q*T)^22 and point counts 1+q^(2*n)+22*q^n.
The last count is odd and exceeds 2^53; the leading P2 coefficient is
2^264. A separate Frobenius-matrix control over GF(4) checks finite
height 3 with v_q(q)=1.

Fixtures are equations in test source; no census dataset is added.
The prototypes below were read for conventions and are not runtime
dependencies. Paths are relative to the prototype root `artifacts/mw/`.

| Source | SHA-256 |
| --- | --- |
| `census-invariants/oracle.m` | `2f761af4ea9d9a6630bbccd3b5323434552e1e7dcddbc8432e1f22482128ae5a` |
| `census-invariants/ns_degree.m` | `55fb52aefaabdff7df195c58be4c643b1e69a8b0ac84c617f11d110d0331c40c` |
| `census-invariants/derive.cpp` | `cff4e2da73d90f81cd0bc720ff511c2f6e929cf8b6bc9f0b14b2b44ed16550b5` |
