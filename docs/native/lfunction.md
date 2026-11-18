# Signed L-polynomials

`FibrationLPolynomial(E)` computes the signed L-polynomial of an elliptic K3
model over `GF(2)(t)` or `GF(3)(t)`.
The input uses polynomial coefficients a1, a2, a3, a4, a6 of degrees at most
2, 4, 6, 8 and 12. The equation must define an everywhere minimal elliptic K3;
the intrinsic checks this through `FibrationProperties`. Coefficients need
not be normalized. Other constant fields, coefficients with poles, and the
other census cells are rejected.

```c++
F<t> := FunctionField(GF(2));
E := EllipticCurve([F | 0, 0, 1, t^6, 0]);
L := FibrationLPolynomial(E);
assert L`Degree eq 7 and L`Sign eq 1;
assert L`Lower eq [-2, 0, 0];
assert Coefficients(L`Polynomial) eq [1, -2, 0, 0, 0, 0, -64, 128];
```

## Result record

| Field | Meaning |
| --- | --- |
| `Polynomial` | L in Z[T], with constant coefficient 1 |
| `Degree` | D = sum of degree-weighted conductor exponents minus 4 |
| `Sign` | Functional-equation sign, either -1 or 1 |
| `Lower` | [c1,...,c_floor(D/2)]; c0 = 1 is implicit |
| `TraceSums` | Euler trace sums S1,...,S_floor(D/2) |
| `Height`, `TwistHeight` | Minimal discriminant degrees divided by 12 for E and E^(p) |
| `TamagawaExponent` | v_p(C(E^(p))/C(E)), with one Tamagawa factor per closed place |

For L(T) = sum c_i*T^i, recover every omitted coefficient by

```text
c_(D-i) = Sign*p^(D-2*i)*c_i.
```

For even D and negative sign the middle coefficient must vanish. All
coefficient arithmetic is integral, and each Newton division is checked.
`Lower` has exactly floor(D/2) entries, including any terminal zeros.

## Local data, sign and counts

The [local records](local.md) supply the full conductor, including
wild terms, and the minimized discriminants and Tamagawa numbers of E
and its Frobenius twist. The twist has coefficients a_i^p. Its height need
not be 2. The [sign](proofs/frobenius-sign.md) is

```text
(-1)^(TwistHeight-Height+TamagawaExponent).
```

The conductor determines D, while the height
difference and Tamagawa valuations determine the sign.

For each degree through floor(D/2), the computation enumerates one point
from every full Frobenius orbit on the affine base. Ordinary good-fibre
traces use the [certified exact tables](tables.md), with their
quadratic or Artin-Schreier twist. Supersingular good fibres use Magma's
finite-field elliptic-curve cardinality. Bad fibres use the minimized
local traces. Infinity is included once, using its local record or the
weight-2i chart coefficients for a good fibre.

For a good closed place of degree d, the trace recurrence is

```text
a_0 = 2, a_1 = a, a_r = a*a_(r-1)-p^d*a_(r-2).
```

Its contribution to S_(d*r) is d*a_r. At a bad place use a^r, where a
is 1, -1 or 0. Newton's identity gives
n*c_n = sum_(j=1)^n S_j*c_(n-j). The known degree and sign complete L
from these coefficients.

## Independent Magma comparisons

The Magma test in `tests/native_lfun_magma.m` compares each fixed
equation with Magma's built-in `LFunction`.

The controls cover both signs, wild reduction, reciprocal charts and
unnormalized equations. Literal signs and Newton sums supplement polynomial
equality. `LFunction` shares Magma's local routines with the package.
Some F_3 inputs exceed the independent oracle's time limit and are excluded.

## Native L-polynomials

`census-native invariants --diagnostic` computes the signed L-polynomial for
each supported key. The input and cell labels follow the [model format](models.md),
and the `degree` and `sign` columns come from the [local records](local.md).
The command appends:

| Column | Encoding |
| --- | --- |
| `lower_L` | c1,...,c_floor(degree/2), separated by spaces, including terminal zeros |
| `L_polynomial` | All coefficients c0,...,c_degree in ascending powers, separated by spaces |
| `trace_sums` | S1,...,S_floor(degree/2), separated by spaces |

An empty list is `[]`. Other census cells retain their labels and have empty
arithmetic columns. Supported rows have status `partial`, which records their
eligibility for arithmetic.

The default `--backend auto` checks CPU support at runtime and selects AVX2
on supported x86 hosts. `--backend scalar` selects scalar counting;
`--backend avx2` requires AVX2 and reports an error when it is unavailable.
AVX2 functions have their own compiler targets; the executable is built
without a global AVX2 requirement. Vector evaluation works in both fields;
the F_3 ordinary trace sum also uses vector gathers and exact integer residues.
On other architectures the executable uses scalar counting.
