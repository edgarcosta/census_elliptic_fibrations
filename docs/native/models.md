# Native census models

The native command decodes census keys into equations on both charts and
classifies their generic fibres and minimality. Thanks for helping check
these computations against the Magma package.
See [native/README.md](../../native/README.md) for build and run commands.

## Input and output

Pass one decimal `q key` pair per line to `census-native invariants`.
Only the prime fields `q=2` and `q=3` are supported. Keys satisfy
`0 <= key < 2^32` and `0 <= key < 3^27`, respectively.
The [key guide](../keys.md) specifies the normalized coefficient blocks.

The command rejects malformed or oversized integers before decoding.
Blank lines are ignored. An invalid row stops the command with a nonzero
exit status and a line-numbered error on standard error.
Earlier valid rows may already have been written.

The input chart is `a_i(t,1)`. You can state it explicitly with
`invariants --chart 'a_i(t,1)'`; any other chart argument is rejected.
The output is uncompressed CSV. `--diagnostic` selects the same CSV. Its model columns are:

| Column | Meaning |
| --- | --- |
| `q`, `key` | Prime field and key re-encoded from the decoded coefficients |
| `kind` | `smooth`, `regular cusp`, or `nonregular` |
| `minimality` | `everywhere minimal` or `nonminimal somewhere` |
| `status` | `partial` for smooth everywhere minimal equations; `unsupported` for the other five cells |
| `chart` | `a_i(t,1)`, quoted as a CSV field |
| `a1`, `a2`, `a3`, `a4`, `a6` | Coefficients of the five polynomials in the input chart |
| `infinity_a1`, `infinity_a2`, `infinity_a3`, `infinity_a4`, `infinity_a6` | Coefficients in the chart `a_i(1,t)` |

Polynomial lists contain space-separated prime-field integers in ascending
powers, with trailing zeros removed and the zero polynomial written as `0`.
To change charts, reverse a list padded to its homogeneous degree plus one.
The five degrees are `2,4,6,8,12`. The two charts satisfy
`a_i(1,u) = u^(2i)*a_i(1/u,1)`, with `u=1/t`.
The fibre coordinates change by `x_infinity=u^4*x`, `y_infinity=u^6*y`.
The `partial` status records eligibility for arithmetic. These rows contain
the decoded coefficients, classification, [local fibre and sign data](local.md),
[complete torsion](torsion.md), [signed L-polynomials](lfunction.md),
[derived invariants](invariants.md), and [Sha/regulator candidates](candidates.md).

## Classification and minimality

A nonzero discriminant gives a smooth generic fibre.
For a singular generic fibre over `F_2(t)`, a regular cusp has
`a1=a3=0` and at least one of `a4,a6` is not a square.
Over `F_3(t)`, a regular cusp has `a2=a4=0` and `a6` is not a cube.
These power tests take place over the function field before algebraic closure.
The remaining singular equations are nonregular.

Minimality uses integral coordinate translations and weighted divisibility
at closed places, including infinity. Outside the cuspidal prefixes,
the coefficient degrees restrict possible nonminimal places to degree at
most two. In the cuspidal prefixes, rational places suffice unless every
place is bad. The [mathematics guide](../mathematics.md#local-minimality)
gives the local congruences and the degree bounds.

The native computation follows the Magma implementations in
`magma/models/models.m`, `classification.m`, `minimality.m` and `forms.m`.
The package exposes `aInvariantsFromKey`, `FibrationKey`,
`FibrationProperties` and `IsEverywhereMinimal` through the root `spec`.

## Comparisons

`tests/compare/native_model.m` compares the public command with those
package intrinsics on fixed formulas, chart exchanges, every coefficient
position, quadratic nonminimal places and a fixed integer recurrence.
The recurrence supplies 100 distinct keys in each characteristic.
The twelve six-cell formulas are shared with `tests/native_model_magma.m`
in `tests/support/_native_model_samples.m`.
The Magma checks also verify chart gluing and rejected keys and models.
The tests generate their comparison file at runtime; no census dataset is
required for these model comparisons.
