# Native model records

The native model command lets you inspect census keys and the two chart
equations. Thanks for helping check these records against the Magma package.

Build with `make census-native`, then pass one decimal `q key` pair per line:

```sh
printf '2 1346371584\n3 951893735796\n' | artifacts/native/census-native invariants
```

Only the prime fields `q=2` and `q=3` are supported. Keys satisfy
`0 <= key < 2^32` and `0 <= key < 3^27`, respectively. The command rejects
malformed or oversized integers before decoding. Blank lines are ignored.
An invalid row stops the command with a nonzero exit status and a line-numbered
error on standard error. Earlier valid rows may already have been written.

The input chart is `a_i(t,1)`. You can state it explicitly with
`invariants --chart 'a_i(t,1)'`; any other chart argument is rejected.
The output is uncompressed CSV with these columns:

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
`a_i(1,u) = u^(2i)*a_i(1/u,1)`.

This command currently computes model data only. All arithmetic fields are
absent, and neither partial nor unsupported rows enter BSD computations.
The compressed invariant record format will be added with the arithmetic
writer.

The computations follow `src/models.m`, `src/classification.m`,
`src/minimality.m` and `src/forms.m`. `make check/native_model.m` compares
the public command with those package intrinsics on fixed formulas, chart
exchanges, every coefficient position, quadratic nonminimal places and a
fixed integer recurrence. The twelve six-cell formulas are shared with
`make check/native_model_magma.m` in `tests/_native_model_samples.m`.
The tests generate their comparison file at runtime; no census dataset is
required.
