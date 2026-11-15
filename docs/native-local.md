# Magma local fibre records

`FibrationLocalData(E : FrobeniusTwist := false)` returns the bad-place records
of an elliptic curve over `GF(2)(t)` or `GF(3)(t)`. Attach `spec` first.
The curve may have nonminimal coefficients. Each row uses local minimal data.
Constant curves with good reduction everywhere return an empty sequence.

```c++
F<t> := FunctionField(GF(2));
E := EllipticCurve([F | 0, 0, 1, t^6, 0]);
fibres := FibrationLocalData(E);
twisted := FibrationLocalData(E : FrobeniusTwist := true);
assert fibres[1]`DiscriminantValuation eq 24;
assert twisted[1]`DiscriminantValuation eq 12;
```

| Field | Meaning |
| --- | --- |
| `Place` | Integer `sum c_i*p^i` for the monic irreducible polynomial `sum c_i*t^i`, or the string `"inf"` |
| `Degree` | Residue-field degree of the closed place; 1 at infinity |
| `Type` | Kodaira family `I`, `I*`, `II`, `III`, `IV`, `IV*`, `III*` or `II*` |
| `Index` | The index n in `I_n` or `I_n*`; zero for the other families |
| `DiscriminantValuation` | Valuation of the local minimal discriminant |
| `Conductor` | Full conductor exponent, including wild ramification |
| `Tamagawa` | Tamagawa number over the residue field |
| `Trace` | 1 for split multiplicative reduction, -1 for nonsplit multiplicative reduction, 0 for additive reduction |
| `ComponentCycles` | Sorted cycle lengths of base-field Frobenius on geometric non-identity fibre components |

Finite rows are sorted by `Place`, followed by infinity. For a closed place
of degree d with m geometric components, the sum of `ComponentCycles` is
`d*(m-1)`. These cycles already include the d conjugate fibres. The omitted
identity component contributes to the fibre class in the trivial lattice.
An irreducible bad fibre has an empty cycle list.

The package uses Magma's [LocalInformation](https://magma.maths.usyd.edu.au/magma/handbook/text/1575)
for Tate reduction and [FrobeniusActionOnReducibleFiber](https://magma.maths.usyd.edu.au/magma/handbook/text/1582)
for component action. Ogg's identity is `f = v(Delta_min)+1-m`.
For example, the characteristic-2 curve above has type `I_9*`, m = 14,
and conductor 11. Replacing its conductor by the tame additive value 2
would lose the wild contribution 9.

With `FrobeniusTwist := true`, the input coefficients become `a_i^p`, where
p is the characteristic. Tate reduction then minimizes at every bad place,
including infinity. This is the Frobenius twist of the curve over the function
field: powers of t also change. Its discriminant divisor need not have degree
24, even when the original curve is a smooth minimal K3 model. This intrinsic
does not compute a functional-equation sign.

Run `make check/native_local_magma.m` for the fixed controls. The test compares
reciprocal charts, pins split and nonsplit multiplicative traces and cycles,
checks wild fibres and minimized twists, and verifies Ogg's identity.
Five smooth minimal K3 fixtures also have discriminant divisor of degree 24.
The infinity twist controls use `y^2+y=x^3+t^d*x` for d = 6 and 7,
whose twists have minimal discriminant degrees 12 and 48 respectively.

The reference prototypes were read when packaging these routines; they are
not runtime dependencies. Paths below are relative to their `artifacts/mw/`
root. The fixed equations are test source; no census data file is added.

| Reference | SHA-256 |
| --- | --- |
| `lfun-csign/oracle.m` | `79ba87a2ee1f72aa0bdde54cddd4643796ddd0cf6dc6ac6c804b57058ea996de` |
| `lfun-csign/infinity-check.m` | `0ac0e8e6073131f6662210135c6858dd60f34e8ef6c0c872d45daaf01e3f6487` |
| `census-invariants/oracle.m` | `2f761af4ea9d9a6630bbccd3b5323434552e1e7dcddbc8432e1f22482128ae5a` |
| `lfun-csign/test_csign.py` | `a70cae73dcb6564438564787e25e1b305d3b8661074e183cde5b1a10ea6faed9` |
