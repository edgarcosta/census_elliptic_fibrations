# Complete Magma torsion for elliptic K3 models

`FibrationTorsion(E)` returns the torsion invariant factors and every torsion
point of a supported elliptic K3 equation. Attach `spec` first.
The constant field must be `GF(2)` or `GF(3)`. The five coefficients must be
polynomials of degrees at most 2, 4, 6, 8 and 12, respectively.
Their homogeneous equation must be everywhere minimal. Coefficients need
not be normalized. Other models are rejected before the bounded search.

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

The [torsion proof](native-mathematics.md#complete-torsion-search-from-the-k3-fibre-corrections)
gives a complete search bound. The implementation uses the geometric
component maxima and the degree weights from `FibrationLocalData(E)`:

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
characteristic-primary torsion. The implementation does not call
`TorsionSubgroup`, search for free generators, or certify a basis.

Run `make check/native_torsion_magma.m` for the fixed sample. The test
compares complete point sets, invariant factors, point orders, inverses, and
addition with Magma's `TorsionSubgroup`. It covers both characteristics,
trivial torsion, orders 2, 3 and 6, invariant factors `[2,4]`, rank zero,
unnormalized coefficients, wild infinity fibres, and rejected models.
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
