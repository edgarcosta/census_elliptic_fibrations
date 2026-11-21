# Complete F_2 comparison

The full comparison needs Magma, the [native build
dependencies](../../native/README.md), Python 3.11 or later, and GNU parallel.

The package enumerates all six F_2 cells using `FibrationRepresentativeKeys`.
The smooth, everywhere minimal population contains 1,397,598 distinct keys.
`tests/fixtures/published-f2/sources.json` pins the population SHA-256 over sorted keys,
each encoded as a four-byte little-endian unsigned integer without a header.
The other five cell counts are 4,134, 680, 38, 7 and 11.

The published `census.csv.zip` contains 1,249,569 rows. The test downloads it
from the public revision pinned in `sources.json` and verifies both archive
and CSV hashes before reading rows. The key population is generated during the check.

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

The public native command computes expanded invariants for every smooth minimal key.

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
conductor, Tamagawa number and Kodaira symbol.

Evidence is written under `build/full-f2/`. `summary.json` records
coverage and agreement counts.
`matches.csv` records every orbit search;
`unmatched-published.csv` and `unmatched-native.csv` record every unmatched
row, including an empty file with a header when appropriate.
`disagreements.csv` records every differing field. The expected outcome is
zero unmatched published rows and 148,029 unmatched native rows. The published
archive omits those native rows, so they receive only the identity checks.

CI uploads the summary, correspondence, unmatched rows, disagreements and logs
when the full-comparison job fails.
The complete F_3 census, Sha refinement,
generators, and release manifests are outside this check.
