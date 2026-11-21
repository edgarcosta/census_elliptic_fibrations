# Local fibre records

`FibrationLocalData(E : FrobeniusTwist := false)` returns the bad-place records
of an elliptic curve over `GF(2)(t)` or `GF(3)(t)`.
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

## Native local data and sign

Each input row is a prime q and census key in the `a_i(t,1)` chart.
Smooth, everywhere minimal K3 rows include `fibres`, `twist_fibres`, `degree`,
`sign`, `height`, `twist_height` and `tamagawa_exponent`.
The `degree` column is `sum(Degree*Conductor)-4`, the L-polynomial degree.
Their status is `partial` to record their eligibility for arithmetic.
Other census cells have `unsupported` status and empty arithmetic columns.

Each fibre is encoded as
`place:degree:type:index:v_delta:conductor:c:trace:cycles`.
Semicolons separate fibres; spaces separate component-cycle lengths.
A known empty fibre or cycle list is `[]`.
The fields and sorting follow the Magma record above.
Tamagawa numbers are over the residue field of each closed place.

The sign is (-1)^(twist_height-height+tamagawa_exponent), as
[proved here](proofs/frobenius-sign.md). Each Tamagawa number contributes
once, without a residue-degree exponent. Wild conductors enter the L-degree.
Local failures report the input line and place and exit with nonzero status.
