# Complete F_2 comparison

Thanks for helping check the census. The full comparison needs Magma, the
[native build dependencies](../../native/README.md), and Python 3.11 or later.
See the [build and test commands](../../native/README.md#complete-f_2-comparison).
The check also runs through `tests/run_tests.m`, `make check`, and
`make check-native-invariants`.

The package enumerates all six F_2 cells using `FibrationRepresentativeKeys`.
The enumeration source under `magma/` matches package base `c1560f2` apart from
import paths and the coordinate conversions added to `magma/models/models.m`.
The smooth, everywhere minimal population contains 1,397,598 distinct keys.
`tests/fixtures/full_f2/sources.json` pins the population SHA-256 over sorted keys,
each encoded as a four-byte little-endian unsigned integer without a header.
The other five cell counts are 4,134, 680, 38, 7 and 11.

The published `census.csv.zip` contains 1,249,569 rows. It is stored through
Git LFS, with the source repository, source commit and both archive and CSV
hashes in `sources.json`. The key population is generated during the check.

## Coordinate changes

`tests/native/f2_orbit_match.c` applies all six matrices in
GL(2, F_2), together with `x -> x + r`, `y -> y + p*x + w`, where
`r = p^2 + a1*p`, `deg(p) = 2` and `deg(w) = 6` as homogeneous forms.
There are 6,144 changes preserving `a2 = 0`. Every published model is
matched against the generated representatives; raw key equality is not
required. Matches must be one to one and land in the smooth minimal cell.

The matcher takes three file arguments: representative CSV, published keys
and match CSV. The representative CSV has `key,cell,prefix`
columns, with cell zero denoting smooth minimal models. The match CSV
records the matrix index and the masks for p and w, as well as missing,
multiple or wrong-cell matches.

## Arithmetic checks and evidence

The public native command computes expanded invariants once for every
smooth minimal key. Up to four worker processes run batches of 10,000 keys,
with each worker checking identities and comparing published values through
a read-only SQLite connection. Results are written in key order. Every row
must have the requested key, smooth minimal kind and complete arithmetic.
`tests/compare/native_records.m` checks compressed records against Magma, including
the fixed F_3 models.

Every row is checked for its functional equation, conductor degree,
discriminant degree, Ogg relation and component cycles. Further checks
cover rank parity, Shioda-Tate ranks, Newton point counts, P2 reconstruction,
the BSD leading term, the BSD/Artin-Tate relation and candidate products.
These checks reuse the local, derived and candidate
test helpers. They test identities; the published table provides the
independent arithmetic comparison on matched rows.

All 17 shared columns are compared on every matched row. These include
a1 orbit class, torsion, the bad-fibre multiset, L-polynomial, sign,
analytic rank, P2, and zeta sign. The remaining columns are rational and
geometric Picard numbers, Brauer height, Artin-Tate product and five point
counts. The a1 class matches by construction because coordinate changes
preserve it; this column does not provide an independent arithmetic check.
The published a1 text is also checked against its packed key. Fibre
positions can move under the base transformation, so comparison uses degree,
conductor, Tamagawa number and Kodaira symbol. Native runs use a
failing Magma executable to detect a subprocess fallback.

Evidence is written under `artifacts/full-f2/`. `summary.json` records
coverage, agreement counts, worker count and setup and native-loop timings.
`matches.csv` records every orbit search;
`unmatched-published.csv` and `unmatched-native.csv` record every unmatched
row, including an empty file with a header when appropriate.
`disagreements.csv` records every differing field. The expected outcome is
zero unmatched published rows and 148,029 unmatched native rows. The published
archive omits those native rows, so they receive only the identity checks.

CI runs one job per test file. The `tests/compare/native_full_f2.m` job fetches only
the census archive through LFS. It retains the summary, correspondence,
unmatched rows, disagreements, and logs.
The job has a 30-minute limit. The complete F_3 census, Sha refinement,
generators, and release manifests are outside this check.
