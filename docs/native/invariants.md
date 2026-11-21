# Exact derived K3 invariants

`FibrationInvariants(E)` computes the invariants below for an everywhere
minimal elliptic K3 polynomial model over GF(2)(t) or GF(3)(t).
It obtains the signed
L-polynomial, bad-place records, and complete torsion through
[`FibrationLPolynomial`](lfunction.md),
[`FibrationLocalData`](local.md) and
[`FibrationTorsion`](torsion.md).

`FibrationInvariants(E, L, fibres)` accepts a complete L-polynomial
and fibre component cycles over any finite constant field F_q.
The ranks and Sha-regulator product come from `AnalyticInformation(E, L)`,
which computes torsion and Tamagawa numbers from E.
The overload also checks the integer polynomial's constant coefficient and
resulting cohomology degree. It does not certify the supplied arithmetic or
the K3 hypothesis. Fibre records need only `ComponentCycles`; a cycle already
includes conjugate fibres at a closed place.

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
NS discriminant nor an integral NS lattice.

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

These are the [chi(O_X)=2 normalizations](proofs/frobenius-sign.md). Every closed place contributes
one Tamagawa factor c_v, including infinity. There is no residue-degree
exponent on c_v. The Artin-Tate product is a positive integer;
the Sha-regulator product remains an exact rational number.

## Independent checks

Magma's `LFunction`, `FrobeniusActionOnTrivialLattice`, `ArtinTateFormula`,
Newton polygons and the zeta power series supply independent comparisons.

For the supersingular equation y^2+y=x^3+t^6*x over GF(2),
L=(1-2*T)*(1-(2*T)^6) and all 13 non-identity component classes are fixed.
Over GF(2^12), all 22 eigenvalues are q. The test supplies L=(1-q*T)^7
and checks P2=(1-q*T)^22 and point counts 1+q^(2*n)+22*q^n.
The last count is odd and exceeds 2^53; the leading P2 coefficient is
2^264. A separate Weil-polynomial extension control over GF(4) checks finite
height 3 with v_q(q)=1.

## Native command

`census-native invariants --diagnostic` derives these fields from the same in-memory
L-polynomial, local records, and complete torsion used by its diagnostic CSV.
All polynomial coefficients, products, and counts use arbitrary-precision
integers. The Sha-regulator product uses a reduced rational with a positive
denominator. Infinite formal height is written `infinity`.

| CSV columns | Magma record fields |
| --- | --- |
| `rank`, `geometric_rank` | `Rank`, `GeometricRank` |
| `picard_number`, `geometric_picard_number` | `PicardNumber`, `GeometricPicardNumber` |
| `P2`, `zeta_sign`, `point_counts` | `P2`, `ZetaSign`, `PointCounts` |
| `brauer_height`, `artin_tate_product` | `BrauerHeight`, `ArtinTateProduct` |
| `ns_splitting_degree` | `NeronSeveriSplittingDegree` |
| `j_degree`, `j_inseparability_exponent` | `JDegree`, `JInseparabilityExponent` |
| `sha_reg` | `ShaRegulatorProduct` |

Lists use ascending polynomial coefficients or increasing extension degree.

Native arithmetic accepts prime-field models. Cyclotomic factors are generated
exactly by FLINT: degree at most 20 and phi(n)>=sqrt(n/2) bound orders by 800.
