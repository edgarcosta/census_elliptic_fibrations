# Complete torsion for elliptic K3 models

`FibrationTorsion(E)` returns the torsion invariant factors and every torsion
point of a supported elliptic K3 equation.
The constant field must be `GF(2)` or `GF(3)`. The five coefficients must be
polynomials of degrees at most 2, 4, 6, 8 and 12, respectively.
Their homogeneous equation must be everywhere minimal. Coefficients need
not be normalized. Other models are rejected before computing torsion.

```c++
F<t> := FunctionField(GF(2));
E := EllipticCurve([F | t^2, 0, 1, t^8, t^6]);
factors, points := FibrationTorsion(E);
assert factors eq [2];
assert E![1/t^2, 1/t^3] in points;
```

The invariant factors are positive integers greater than 1 in divisibility
order. Trivial torsion has factors `[]` and the single point `E!0`.
The second result is a sequence of points on the input curve, with the
identity first. No order is promised for the remaining points.

The Magma routine calls `TorsionSubgroup(E)` and returns its invariant
factors and the images of all group elements on the input curve.

The [complete search bound](proofs/torsion-bound.md) uses the
geometric component maxima and the degree weights from the local records:

```text
Cmax = sum_v degree(v)*maximum_component_correction(v) <= 6,
m = floor((Cmax-4)/2) <= 1.
```

If m is negative, the identity is the whole torsion group. Otherwise, the
search enumerates all monic denominators g of degree at most m and all numerators A
of degree at most 4+2m with gcd(A, g)=1. For each x=A/g^2, exact quadratic
roots give all y=B/g^3 on the input Weierstrass equation. The search retains
only polynomial B of degree at most 6+3m, sections with (P.O)<=m, and exact
height zero. The intersection includes the reciprocal chart at infinity.
Wild fibres use their Kodaira component corrections; their discriminant
cost and the proof's global bound include wild ramification.

The identity is inserted separately. Repeated point addition gives the
orders of all enumerated points and their least common multiple e.
Elliptic torsion has at most two invariant factors, so they are
`[#points div e, e]` after removing factors equal to 1. This includes
characteristic-primary torsion. The native implementation does not call
Magma, search for free generators, or certify a basis.

## Native torsion

`census-native invariants --diagnostic` computes complete torsion for each [supported key](models.md).
The arithmetic columns include:

| Column | Encoding |
| --- | --- |
| `torsion` | Invariant factors greater than 1 in divisibility order, separated by spaces; `[]` for trivial torsion |
| `torsion_points` | All torsion points, separated by semicolons; `O` first, then `x_numerator:x_denominator:y_numerator:y_denominator` |

Each coordinate polynomial is encoded by the nonnegative integer
`sum coefficient_i*q^i`. Rational coordinates
are reduced with monic denominators. Other census cells retain their labels
and have empty torsion columns. Supported rows have status `partial`, which
records their eligibility for arithmetic.

Exact rational local heights filter out every positive-height section.
The local heights use the psi2, c4 and psi3 cases of Theorem 5.2 in
Silverman, "Computing heights on elliptic curves", Math. Comp. 51 (1988),
rescaled to Shioda corrections; the completeness of the search rests on
this filter as well as on the proved bound.
The computation includes infinity on the weight-4 and weight-6 coordinate
charts. It uses the original minimal discriminant valuations, including wild
ramification. Heights use arbitrary precision.
