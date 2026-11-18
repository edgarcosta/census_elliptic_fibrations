# Native census models

The native command decodes [census keys](../keys.md) into equations on both charts and
classifies their [generic fibres](../mathematics.md#generic-fibre-kind) and [minimality](../mathematics.md#local-minimality).

## Input and output

Pass one decimal `q key` pair per line to `census-native invariants`.
Only the prime fields `q=2` and `q=3` are supported.

The command rejects malformed or oversized integers before decoding.
Blank lines are ignored. An invalid row stops the command with a nonzero
exit status and a line-numbered error on standard error.
Earlier valid rows may already have been written.

The input chart is `a_i(t,1)`.
The output is uncompressed CSV. `--diagnostic` selects the same model
columns:

| Column | Meaning |
| --- | --- |
| `q`, `key` | Prime field and input key |
| `kind` | `smooth`, `regular cusp`, or `nonregular` |
| `minimality` | `everywhere minimal` or `nonminimal somewhere` |
| `status` | `partial` for smooth everywhere minimal equations; `unsupported` for the other five cells |
| `a1`, `a2`, `a3`, `a4`, `a6` | Coefficients of the five polynomials in the input chart |
| `infinity_a1`, `infinity_a2`, `infinity_a3`, `infinity_a4`, `infinity_a6` | Coefficients in the chart `a_i(1,t)` |

Polynomial lists contain space-separated prime-field integers in ascending
powers, with trailing zeros removed and the zero polynomial written as `0`.
To change charts, reverse a list padded to its homogeneous degree plus one.
The five degrees are `2,4,6,8,12`. The two charts satisfy
`a_i(1,u) = u^(2i)*a_i(1/u,1)`, with `u=1/t`.
The fibre coordinates change by `x_infinity=u^4*x`, `y_infinity=u^6*y`.
The `partial` status records eligibility for arithmetic; these rows contain
the decoded coefficients, classification, and [local fibre and sign data](local.md).

## Classification and minimality

The native computation follows the Magma implementations in
`magma/models/models.m`, `classification.m`, `minimality.m` and `forms.m`.
The package exposes `aInvariantsFromKey`, `FibrationKey`,
`FibrationProperties` and `IsEverywhereMinimal` through the root `spec`.
