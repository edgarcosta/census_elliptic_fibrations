# Native records

## Native model command

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
The output is uncompressed CSV beginning with these model columns:

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

Supported rows also contain [local records and signs](native-local.md#native-local-data-and-sign)
and [complete torsion](native-torsion.md#native-torsion).
Their status remains `partial` while the remaining invariants are pending.
Unsupported rows have empty arithmetic columns. The stage 1 Magma reference
format is described below; connecting it to the native writer is a separate step.

The computations follow `src/models.m`, `src/classification.m`,
`src/minimality.m` and `src/forms.m`. `make check/native_model.m` compares
the public command with those package intrinsics on fixed formulas, chart
exchanges, every coefficient position, quadratic nonminimal places and a
fixed integer recurrence. The twelve six-cell formulas are shared with
`make check/native_model_magma.m` in `tests/_native_model_samples.m`.
The tests generate their comparison file at runtime; no census dataset is
required.

## Stage 1 invariant rows

The Magma reference codec constructs and encodes the lean inputs to the native
invariants writer. Attach `spec` from the repository root to use these intrinsics:

```text
row := FibrationInvariantRow(key, W);
header := FibrationInvariantCSV();
text := FibrationInvariantCSV(row);
copy := FibrationInvariantRow(text);
```

`W` is a homogeneous five-coefficient census model over GF(2) or GF(3).
The caller supplies the full coefficient key `FibrationKey(W)`.
The constructor classifies
`W` with `FibrationProperties`. For a smooth, everywhere minimal model it uses
`FibrationLPolynomial`, `FibrationLocalData`, `FibrationTorsion`,
`FibrationInvariants` and `FibrationCandidates`. The other five census cells
receive `unsupported` arithmetic, with their classifications intact.

For partial results or failures, begin with
`FibrationInvariantRow(key, q, kind, minimality)` and assign known fields.
Its supported-cell status starts as `partial`; other cells start as `unsupported`.
Unknown optional fields contain the empty string `""`. A known empty list is an
empty Magma sequence or list. The encoder checks the wire types and stage 1
status rules. It does not certify caller-supplied arithmetic or key provenance.

The string overload is a reference codec for one row. Release-file reading,
manifest validation, sidecar validation, zstd compression, sharding, and resumption belong
to stage 3. No prototype path is an installed dependency.

### Wire format

Version 1 uses ASCII CSV with the header returned above. Each record has 18
columns, in the order below, and carries its schema version and field size.
Package keys include the prefix coefficients. The enclosing job must retain
the prefix forms and key-encoding version for release validation. This row alone
is not a release manifest.

`Z` is a signed decimal integer of arbitrary precision; `N` is a nonnegative `Z`.
Canonical integers have no leading plus or leading zero, except for zero itself.
`Q` is a reduced `Z/N` with positive denominator, including `/1` for integers.
No number passes through floating point. Scalar bounds below are inclusive.

Lists contain space-separated numbers. Empty optional cells mean unknown;
`[]` means a known empty list. Nonempty lists have no brackets. Fibre tuples
use colons between entries and semicolons between tuples. CSV quotes surround
cells containing commas, double quotes, CR or LF; double quotes are doubled.
The encoder adds no record separator. Set `SetColumns(0)` before using Magma
`Write` to prevent line wrapping. The decoder accepts exactly one record,
without its trailing separator; CR and LF inside quoted messages are retained.

| CSV column | Magma field | Type and meaning |
| --- | --- | --- |
| `schema_version` | `Version` | `1` |
| `q` | `Q` | `2` or `3` |
| `key` | `Key` | `N`, at most `2^64-1` |
| `kind` | `Kind` | `smooth`, `regular cusp`, `nonregular` |
| `minimality` | `Minimality` | `everywhere minimal`, `nonminimal somewhere` |
| `status` | `Status` | `ok`, `partial`, `unsupported`, `error` |
| `degree` | `Degree` | Optional `N`, 0..20, the degree D of L |
| `sign` | `Sign` | Optional `-1` or `1`, the functional-equation sign |
| `lower_L` | `LowerL` | Optional `Z` list, c1 through c_floor(D/2); c0=1 is implicit |
| `fibres` | `Fibres` | Optional list of eight-entry tuples described below |
| `torsion` | `Torsion` | Optional invariant factors, each at least 2, in divisibility order; `[]` is trivial torsion |
| `sha_reg` | `ShaReg` | Optional positive `Q`, the exact product #Sha * Reg |
| `admissible_sha` | `AdmissibleSha` | Optional sorted distinct positive square `N` list |
| `split_method` | `SplitMethod` | Optional `constraints` or `unresolved` in stage 1 |
| `d2` | `D2` | Unknown in stage 1; reserved for the optional stage 2 uint8 value |
| `sha2_dimensions` | `Sha2Dimensions` | Unknown in computed stage 1 rows; reserved for a stage 2 uint8 list |
| `refinement_status` | `RefinementStatus` | Unknown, `not_needed`, `unsupported` or `error` in stage 1 |
| `message` | `Message` | ASCII diagnostic text; empty when there is no diagnostic |

A fibre is `(place,degree,type,n,v_delta,conductor,c,trace)`. The place is `inf`
or its nonnegative polynomial code from `FibrationLocalData`. Degree is 1..255;
infinity has degree 1. The Kodaira family is `I`, `II`, `III`, `IV`, `I*`, `IV*`,
`III*` or `II*`; n carries the index of `I_n` or `I_n*`. The index, discriminant
valuation, conductor and Tamagawa number c fit uint16; c is positive. Trace is
-1, 0 or 1. For example:

```text
81:4:I:6:6:1:2:-1;inf:1:II*:0:12:4:1:0
```

The lower half reconstructs L by
`c_(D-i) = sign*q^(D-2*i)*c_i`. For even D, the central coefficient is stored
once. The lean row omits derived views: local component cycles, ranks, Picard
numbers, P2, point counts, formal Brauer height, Artin-Tate product, NS splitting
degree, and j-data. See
[native-invariants.md](native-invariants.md) for their packaged computation.
A singleton Sha list determines Reg from `sha_reg/#Sha`. Ambiguous lists remain
complete; the encoder never chooses a candidate or infers a geometric lattice.

### Status meanings

`ok` means the stage 1 arithmetic fields are complete and the candidate list is
nonempty. An ambiguous split can still have status `ok`. `partial` retains known
fields while leaving the rest unknown. `unsupported` records a computation
outside the supported domain. `error` records a failed computation and can retain
known inputs and partial results. A known empty candidate list is distinct from
an unknown list; the computed constructor marks that inconsistent result `error`.
Arithmetic fields remain unknown for the five unsupported census cells.

A singleton uses `split_method=constraints` and `refinement_status=not_needed`.
An ambiguous list uses `split_method=unresolved`. Its refinement status is
`unsupported` over GF(2), where the planned stage 2 method does not apply.
Over GF(3) the refinement status stays unknown: eligibility and refinement have
not been computed. Unknown does not assert eligibility or exclusion.

The stage 2 vocabulary also includes `exact` and `bound`, and later split methods
include `d2` and `basis`. Version 1 of this stage 1 codec rejects these claims
and nonempty refinement values. It accepts `[]` for a supplied known empty
`Sha2Dimensions` list and preserves it separately from unknown. The computed
constructor never supplies that list. A later schema must define computed
refinements and their hypotheses before these claims can be emitted.

### Checks

`make check/native_records_magma.m` tests literal CSV bytes, escaped quotes and
newlines, exact integers beyond 2^53, rational reduction, and unknown and empty
lists. It covers all six cells in both characteristics, nontrivial torsion,
rank-zero singletons, and complete ambiguous Sha lists. Fixed curves also compare
reconstructed L with Magma's independent `LFunction`. The test is discovered by `run_tests.m`,
`make check` and the CI matrix, whose jobs have a 30-minute limit.
Test inputs are source-level equations and arithmetic controls, not census data
files.
