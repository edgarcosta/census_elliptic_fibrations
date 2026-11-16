# Local fibre records

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

## Native local data and sign

`make census-native` builds the native command. From the repository root, run:

```sh
printf '2 4884480\n3 372699133912\n' | artifacts/native/census-native invariants --diagnostic
```

Each input row is a prime q and census key in the `a_i(t,1)` chart.
Smooth, everywhere minimal K3 rows include `fibres`, `twist_fibres`, `degree`,
`sign`, `height`, `twist_height` and `tamagawa_exponent`.
Their status is `partial` while the remaining invariant computations are pending.
Other census cells have `unsupported` status and empty arithmetic columns.

Each fibre is encoded as
`place:degree:type:index:v_delta:conductor:c:trace:cycles`.
Semicolons separate fibres; spaces separate component-cycle lengths.
A known empty fibre or cycle list is `[]`.
The fields and sorting follow the Magma record above.
Tamagawa numbers are over the residue field of each closed place.

`local_models` in `native/model.cpp` factors the discriminant once, then reduces
both the original and Frobenius-twisted equation at each bad place.
It uses the infinity chart before taking p-th powers and minimizes the twist
independently. Component cycles include the closed-place degree.
The shared Weierstrass invariants also serve model classification.

`functional_equation_sign` in `native/lfun.cpp` consumes these records without
repeating local reduction. It checks the original discriminant degree is 24,
the twist degree is divisible by 12, and the conductor gives an L-degree
between 0 and 20. The sign uses the [proved exponent](native-lfun.md#local-data-sign-and-counts)
`twist_height-height+tamagawa_exponent`.
Each closed-place Tamagawa number contributes once, without a degree exponent.
The native executable uses NTL and GMP and starts no Magma process.
Local failures report the input line and place and exit with nonzero status.

Run `make check/native_local.m` for comparisons through `census-native invariants --diagnostic`.
The test generates 422 Magma reference rows from fixed equations and keys,
including reciprocal charts and 100 distinct recurrence keys per characteristic.
It compares every local and twist field, component cycles, conductor degree,
height, Tamagawa exponent and sign. The eleven fixed equations also compare
with the signed package L-polynomial computation.
Each characteristic must have at least 100 distinct keys and every Kodaira
family among the original and twist fibres, including `I_n*` with n > 0.
The sample also covers split and nonsplit fibres, nonrational closed places,
wild terms, both signs, and twist heights 1, 2 and 4.
The test checks Ogg identities and both discriminant sums, and runs the native
command with a failing Magma executable on its search path.
The existing runner, Makefile and CI matrix discover the named test.

The Tate engine and component-cycle formulas were adapted from these prototypes.
They are reference sources, with no installed prototype dependency.
No census dataset is added; reference rows are generated during testing.

| Source | SHA-256 |
| --- | --- |
| `lfun-csign/tate.cpp` | `a2bd160ae5b1b5c8847a50361bedbfc8f55adc7c5cd8ae73bec304bebc3c8161` |
| `lfun-csign/tate.hpp` | `b556141722ecf04a5003011eb06ac1bc601ab983bb9fce8e1b235172dc8f527a` |
| `census-invariants/derive.cpp` | `cff4e2da73d90f81cd0bc720ff511c2f6e929cf8b6bc9f0b14b2b44ed16550b5` |
