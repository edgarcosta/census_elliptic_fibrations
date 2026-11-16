# Signed L-polynomials

`FibrationLPolynomial(E)` computes the unique signed L-polynomial of an
elliptic K3 model over GF(2)(t) or GF(3)(t). Attach `spec` from the repository
root before calling it.

The input uses polynomial coefficients a1, a2, a3, a4, a6 of degrees at most
2, 4, 6, 8, 12. The equation must define an everywhere minimal elliptic K3.
The intrinsic checks this through `FibrationProperties`. It accepts
unnormalized equations within these bounds. Other constant fields,
rational coefficients with poles, and the other census cells are rejected.

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

The [local records](native-local.md) supply the full conductor, including
wild terms, and the minimized discriminants and Tamagawa numbers of E
and its Frobenius twist. The twist has coefficients a_i^p. Its height need
not be 2. The sign is

```text
(-1)^(TwistHeight-Height+TamagawaExponent).
```

The [sign proof](native-mathematics.md#frobenius-twist-functional-equation-sign)
fixes the BSD normalization and square-order hypotheses. The conductor
determines D, while the height difference and Tamagawa valuations determine
the sign.

For each degree through floor(D/2), the computation enumerates one point
from every full Frobenius orbit on the affine base. Ordinary good-fibre
traces use the [certified exact tables](native-tables.md), with their
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
from these coefficients. The Magma implementation calls neither `LFunction`
nor a native executable.

## Independent comparisons

Run `make check/native_lfun_magma.m`. The test is also discovered by
`run_tests.m` and the CI matrix, whose job budget is 30 minutes.
Each independent `LFunction` call runs in a separate Magma process,
without attaching this package, with a 120 CPU-second cap. A required
reference that times out fails the test and is identified as a timeout.
A consistency check cannot satisfy the reference comparison.

The fixed sample contains four formerly unresolved-sign examples:
GF(2) keys 4884480 and 11821184, and GF(3) keys 372699133912 and
320009441124. These include both signs and rank zero in both fields.
The test checks that their coefficients around the middle vanish.
Their fresh `LFunction` results supply independent polynomial equality.

The remaining equations cover the three wild/twist controls in the sign
proof, ordinary GF(2) fibres and unnormalized equations in both fields.
Every equation is also checked in its reciprocal chart against the same
independent reference, since that chart change is an isomorphism.
The test checks the complete integral polynomial, sign, degree, lower-half
encoding, functional equation and Newton sums.

Fixtures are equations and keys in test source. No census dataset or
prototype path is needed at runtime. The earlier verification's ten
GF(3) timeouts remain timeouts; they are outside this fixed sample and
are not counted as independent matches here.

The prototype interfaces and fixture provenance are recorded below.
Paths are relative to the existing prototype root `artifacts/mw/`.
These hashes identify reference inputs, not installed dependencies.

| Source | SHA-256 |
| --- | --- |
| `verify-lfun/verify.m` | `17fd26aea244d8db16fdcfafc728157cd0620414c11061d2dd180019eae4326e` |
| `lfun-csign/oracle.m` | `79ba87a2ee1f72aa0bdde54cddd4643796ddd0cf6dc6ac6c804b57058ea996de` |
| `lfun-counter/open_signs.csv` | `d43d5370d8dd50b69c878c93a2ca7c730eff99037c7c4ec358c9b2a6d4ab4615` |
| `lfun-counter/f2_reference.csv` | `4f9fad2460bfe347dfea17e82e95cf71bb286d3f1921015133347d2cc828cd27` |
| `lfun-counter/f3_reference.csv` | `2794f18443bbaf08f69422b65fecff9a4acbb839a672348be6402e6b71c1d7ce` |

## Native command

Build with `make census-native` and run `artifacts/native/census-native invariants --diagnostic`.
Each input line is a field size and census key, for example `2 4884480`.
The diagnostic CSV includes `degree`, `sign`, `lower_L`, `L_polynomial`
and `trace_sums`, alongside the model, local data, and torsion columns.
The coefficient lists use ascending powers, separated by spaces.
`lower_L` omits c0 and has floor(degree/2) entries, including terminal zeros.
An empty list is `[]`; unsupported census cells have blank arithmetic fields.
Supported rows remain `partial` while the remaining invariants are assembled.
The release row encoding is described in [native-format.md](native-format.md).

`--coefficients` reads `q a1 a2 a3 a4 a6` instead of `q key`.
Each coefficient is its nonnegative base-q polynomial code, with the constant
term as the least significant digit and the same degree bounds as above.
This input uses the census normal forms: a2 = 0 over F_2, and a1 = a3 = 0
over F_3. The command derives the key, local data, conductor degree, and sign
from the coefficients. Each row contains only the field and coefficient codes.
Other polynomial equations can be passed directly to `FibrationLPolynomial`
in Magma, whose input supports unnormalized coefficients.

The native counter reuses the original and minimized-twist local records.
It counts through floor(degree/2), using certified tables and one
representative per full Frobenius orbit. Bad places contribute their local
Euler series once, and good infinity is counted in its weight-2i chart.
Supersingular traces use a quadratic character sum over F_3 and elimination
of the quadratic form Tr((x^3+a2*x^2+a4*x+a6)/a3^2) over F_2.
Newton divisions and functional-equation completion use arbitrary precision
integers, with a check on the middle coefficient when the degree is even.
The native process invokes no Magma executable and reads no prototype files.

The default `--backend auto` checks CPU support at runtime and selects AVX2
on supported x86 hosts. `--backend scalar` selects scalar counting;
`--backend avx2` requires AVX2 and reports an error when it is unavailable.
AVX2 functions have their own compiler targets; the executable is built
without a global AVX2 requirement. Vector evaluation works in both fields;
the F_3 ordinary trace sum also uses vector gathers and exact integer residues.
On other architectures the executable uses scalar counting.

`make census-native-scalar` builds `artifacts/native/census-native-scalar`
with the AVX2 code excluded. It accepts the same inputs and options,
with an unavailable-backend error for `--backend avx2`.

Run `make check/native_lfun.m` for the native comparison. It generates 54
package reference records, covering the four former-open-sign keys, the
seven local controls and eight recurrence keys per field, in both charts.
The sample includes both signs, even and odd degrees, wild conductors,
good and bad infinity, ordinary curves and everywhere supersingular curves.
Every full polynomial, lower coefficient and counted Euler sum equals Magma.
The test compares automatic, scalar, available AVX2, and scalar-only builds
on identical mixed-field inputs. Coefficient input is checked on the same
sample, and unsupported cells retain blank arithmetic fields.
The existing Magma comparison supplies independent completed `LFunction`
references for the former-open-sign and wild controls.

The native port uses these prototype sources; none is an installed dependency.

| Source under `artifacts/mw/` | SHA-256 |
| --- | --- |
| `lfun-csign/counter.cpp` | `611e767bb08b77ead608ef29ea157b548b8c65732ca0233b96fea8cc55874c31` |
| `lfun-csign/evaluate.cpp` | `32dabcb98fdf3234d769d6fe421c539f06c5221ea9939e43e42d59f4bd5772c0` |
| `lfun-csign/field.cpp` | `43cacd63b6f23bb2402aacaf02c65cfe31c7933ef5ecf3560a990860ed28198c` |
