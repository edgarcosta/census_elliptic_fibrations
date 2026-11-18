# Exact compressed invariant records

The native invariants writer and Magma reference codec encode the same exact
arithmetic records.
Attach `spec` from the repository root to use these intrinsics:

```c++
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

The string overload is a reference codec for one row. It does not read
release files or validate manifests.

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
| `degree` | `Degree` | Optional `N`, 0 to 20, the degree D of L |
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
or its nonnegative polynomial code from `FibrationLocalData`. Degree is 1 to 255;
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
degree, and j-data.
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

### Compressed native output

`census-native invariants` reads one `q key` pair per line and writes a single
zstd frame at level 9 to standard output. Both fields can occur in one stream.
`--diagnostic` selects the expanded, uncompressed comparison CSV described
in the model and arithmetic guides.

The frame has a content checksum and begins with the column header, followed
by one record per input. Arithmetic uses shared in-memory local data, torsion,
and L-polynomials.

A computation failure emits an `error` record containing known fields and an
ASCII diagnostic, continues with later inputs, and gives a nonzero exit status.
Malformed input stops at that line and reports its line number on standard
error. The preceding records form a complete compressed frame. Output or input
I/O failures also give a nonzero exit status; failed output can be incomplete.
