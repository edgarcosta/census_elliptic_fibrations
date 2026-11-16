# Complete F_2 comparison

Thanks for helping check the census. Run `make check-native-invariants-f2`
from the repository root after `git lfs pull`. The check needs Magma, the
native build dependencies and Python 3.11 or later. It also runs through
`run_tests.m`, `make check` and `make check-native-invariants`.

The package enumerates all six F_2 cells using `FibrationRepresentativeKeys`.
The enumeration source in `src/` agrees with package base `c1560f2`.
The smooth, everywhere minimal population contains 1,397,598 distinct keys.
`tests/native/data/full_f2/sources.json` pins its SHA-256 over sorted keys,
each encoded as a four-byte little-endian unsigned integer without a header.
The other five cell counts are 4,134, 680, 38, 7 and 11.

The published `census.csv.zip` contains 1,249,569 rows. It is stored through
Git LFS, with the source repository, source commit and both archive and CSV
hashes in `sources.json`. The key population and Magma references are
generated during the check; no prototype directory is needed at runtime.

## Coordinate changes

`tests/native/f2_orbit_match.c` uses the orbit search from
`verify-f2run/scripts/orbit_match.c`. It applies all six matrices in
GL(2, F_2), together with `x -> x + r`, `y -> y + p*x + w`, where
`r = p^2 + a1*p`, `deg(p) = 2` and `deg(w) = 6` as homogeneous forms.
There are 6,144 changes preserving `a2 = 0`. Every published model is
matched against the generated representatives; raw key equality is not
required. Matches must be one to one and land in the smooth minimal cell.

The matcher takes four file arguments: representative CSV, published keys,
match CSV and hit-count binary. The representative CSV has `key,cell,prefix`
columns, with cell zero denoting smooth minimal models. The match CSV
records the matrix index and the masks for p and w, as well as missing,
multiple or wrong-cell matches. The binary contains one uint32 hit count
per smooth minimal representative in enumeration order.

## Arithmetic checks and evidence

The public native command computes both compressed records and expanded
invariants for every smooth minimal key. Batches bound the test's memory;
they do not add release sharding. Every compressed row must be complete,
have the requested key and agree with its expanded arithmetic. Unknown
refinement fields and unresolved Sha lists retain their format meanings.

Every row is checked for its functional equation, conductor degree,
discriminant degree, Ogg relation and component cycles. Further checks
cover rank parity, Shioda-Tate ranks, Newton point counts, P2 reconstruction,
the BSD leading term, the BSD/Artin-Tate relation and candidate products.
These checks reuse the local, derived and candidate
test helpers. They test identities; the published table provides the
independent arithmetic comparison on matched rows.

All 17 shared columns are compared on every matched row. These include
a1 orbit class, torsion, the bad-fibre multiset, L-polynomial, sign,
analytic rank, P2 and zeta sign. The remaining columns are rational and
geometric Picard numbers, Brauer height, Artin-Tate product and five point
counts. Fibre positions can move under the base transformation, so comparison
uses degree, conductor, Tamagawa number and Kodaira symbol. Fixed F_3 models from the shared fixture are
compared with `FibrationInvariantRow` byte for byte. Native runs use a
failing Magma executable to detect a subprocess fallback.

Evidence is written under `artifacts/full-f2/`. `summary.json` records
coverage and agreement counts. `matches.csv` records every orbit search;
`unmatched-published.csv` and `unmatched-native.csv` record every unmatched
row, including an empty file with a header when appropriate.
`disagreements.csv` records every differing field. The expected outcome is
zero unmatched published rows and 148,029 unmatched native rows.

CI discovers `tests/native_full_f2.m`, obtains the LFS input and retains
the summary, correspondence, unmatched rows, disagreements and logs.
The job has a 30-minute limit. The complete F_3 census, Sha refinement,
generators and release manifests are outside this check.
