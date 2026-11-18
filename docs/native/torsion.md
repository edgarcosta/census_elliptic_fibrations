# Complete torsion for elliptic K3 models

Thanks for helping check the torsion computation. Build and comparison commands
live in [native/README.md](../../native/README.md#complete-torsion).

`FibrationTorsion(E)` returns the torsion invariant factors and every torsion
point of a supported elliptic K3 equation. Attach `spec` first.
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

The [torsion proof](proofs/torsion-bound.md)
gives a complete search bound for the native implementation. It uses the
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

The fixed Magma test in `tests/native_torsion_magma.m` checks expected
groups and coordinates, distinct points, the identity first, and rejected models.
It covers both characteristics, trivial torsion, orders 2, 3, 4, 5 and 6,
and invariant factors `[2,2]`, `[2,4]` and `[2,6]`.
Other cases include rank zero, unnormalized coefficients and wild infinity fibres.
The degree-one denominator in the example and its reciprocal require m=1;
the reciprocal torsion point has coordinates `(t^6,t^9)`.

Fixed equations are test source, with no installed prototype dependency
or census dataset. The rank-zero keys 1346994688 over GF(2) and 75547326549
over GF(3) come from the first rank-zero entries in the reference samples
below. The test confirms rank zero by evaluating a fresh L-polynomial at
1/p. The remaining controls are explicit equations checked with Magma.

| Reference under the prototype `artifacts/mw/` root | SHA-256 |
| --- | --- |
| `c-stages/reference.m` | `fe11663e7886648f2464969558cc9949822a608434ce01d11f23a01f14c7d27c` |
| `fable/mwk3.m` | `df303efe87f03a8e27b1434e798d97d1311ef1cbfdd7f09397c8b32416468ccf` |
| `census-cost/sample_f2.m` | `daec4ba812200ddfc1df2a5f51cca54e95e9a0218ba1cd505909da231d8e0391` |
| `census-cost/sample_f3.m` | `948f2ce7d02c729cfa7a08bf3d28c25077229ef199c0402810d94e9b82fee0e3` |

## Native torsion

`census-native invariants --diagnostic` computes complete torsion for each supported key.
The input and cell labels follow the [model format](models.md).
The arithmetic columns include:

| Column | Encoding |
| --- | --- |
| `torsion` | Invariant factors greater than 1 in divisibility order, separated by spaces; `[]` for trivial torsion |
| `torsion_points` | All torsion points, separated by semicolons; `O` first, then `x_numerator:x_denominator:y_numerator:y_denominator` |
| `torsion_orders` | Point orders in the same order, separated by spaces |

Each coordinate polynomial is encoded by the nonnegative integer
`sum coefficient_i*q^i`, using arbitrary precision. Rational coordinates
are reduced with monic denominators. Other census cells retain their labels
and have empty torsion columns. Supported rows have status `partial`, which
records their eligibility for arithmetic.

`complete_torsion` in `native/invariants/invariants.cpp` consumes the model and its
shared local records. It repeats neither Tate reduction nor discriminant
factorization. Degree-weighted Kodaira maxima give the bound above.
The search enumerates every permitted A and monic g. In characteristic 2,
Gaussian elimination solves the linear map `B -> B^2+linear*B`.
In characteristic 3, exact polynomial square roots solve the quadratic.
Both roots are retained when distinct.

Exact rational local heights filter out every positive-height section.
The local heights use the psi2, c4 and psi3 cases of Theorem 5.2 in
Silverman, "Computing heights on elliptic curves", Math. Comp. 51 (1988),
rescaled to Shioda corrections; the completeness of the search rests on
this filter as well as on the proved bound.
The computation includes infinity on the weight-4 and weight-6 coordinate
charts. It uses the original minimal discriminant valuations, including wild
ramification. Heights and encoded coordinate integers use arbitrary precision.
The search has at most 17,497 points, so group orders fit in native integers.
Repeated addition checks point orders and membership in the enumerated set;
the exponent and group order determine the invariant factors. No Magma
process, generator search, or basis certification is used.

The native comparison in `tests/compare/native_torsion.m` uses the
same 17 equations as the Magma torsion test, each on both charts.
The characteristic-3 equation with nonzero a1 and a3 is normalized before
key encoding. The reference uses `FibrationTorsion`, which calls
`TorsionSubgroup`; the native command must reproduce all points, orders,
and invariant factors. The F_3 samples include C4, C2 x C2 and C2 x C6.
The test includes both sections requiring m=1 and runs the command with
a failing Magma executable as the only PATH entry.
The shared test list discovers both entries for the runner, Make, and CI.
The fixed equations live in test source; no census dataset is committed.

The native port uses the following prototype sources for point operations,
local heights, group factors and exact quadratic solutions. They are
provenance references, with no runtime dependency.

| Source under the prototype `artifacts/mw/` root | SHA-256 |
| --- | --- |
| `c-stages/core.cpp` | `f853ab9907861ca22547af31ef21680fa1a83618bb39594ab234beef56872c28` |
| `c-stages/stages.cpp` | `02ca9e62fc0067bd7c467455e837e207aaa6d5c20e22bf3acc2c4506b4b1ef97` |
| `search-c/search2.cpp` | `5b67992c70f4b4bcb538a7edc2f0b65a716bfd41366d83f049c6750f2df3c2dc` |
| `search-c/search3.cpp` | `05d1ceaa38e0638ff5edaeea50669e4f1be455e97d866fdc95466ee58dafb3f2` |
