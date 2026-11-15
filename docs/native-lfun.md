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
from these coefficients. The implementation calls neither `LFunction`
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
