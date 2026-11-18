# Published F_2 census fixtures

The sample fixture checks 16 rows from the published F_2 census table, which is incomplete.
It compares coordinate-equivalent models and saved invariants, with no network
request during the test.

Source: [Census of elliptic and quasi-elliptic K3 surfaces over F2](https://github.com/Census-of-elliptic-fibrations/Census-of-elliptic-and-quasi-elliptic-K3-surfaces-over-F2),
revision `9d1b3b85b702811da14d0e53a3d768356e998f27`, file `census.csv.zip`,
Git blob `0f1839c0e0f8746ae8a93b9fda456c250aebcb16`, member `census.csv`.
The table has 1,249,569 of the 1,397,598 smooth everywhere-minimal census rows.
This sample makes no claim of complete-census coverage. Its observed a1 classes
are sampling strata, not full-census classification buckets.

Within each a1 class, preserve CSV order and select zero-based positions
0, n div 3, (2*n) div 3 and n-1. The observed class sizes n are:

| a1 | Rows |
| --- | ---: |
| 0 | 173977 |
| s^2 | 524493 |
| s*t | 510215 |
| s^2 + s*t + t^2 | 40884 |

`source_row` records the one-based data-row position, excluding the header.
The other columns are copied unchanged from `census.csv`: key, a1, torsion
label, bad-fibre multiset, L-polynomial coefficients in ascending degree,
L-sign and analytic rank. Each fibre record gives its degree, conductor
exponent, Tamagawa number and Kodaira symbol. Derived zeta, Picard and height
columns are omitted.

The `census.csv` key encodes a1, a3, a4, a6 in blocks of 3, 7, 9, 13 bits,
with a2 = 0. Bit j in a degree-d block multiplies s^j*t^(d-j).
The test decodes each key through `aInvariantsFromKey` and matches it against
`FibrationRepresentativeKeys`. It searches the allowed changes with q of
degree 2, r = q^2+a1*q, v of degree 6 and all of GL(2, F_2).
It compares models after the substitutions, allowing different representative
keys. The coordinate formulas follow Silverman's *The Arithmetic of Elliptic
Curves*, Table 3.1.

For both the table's model and the package's representative, Magma independently computes
`TorsionSubgroup`, `LocalInformation`, `LFunction` and `AnalyticInformation`.
Bad-fibre names are omitted because base changes move the places; multiplicity
is retained. The sign is checked from the leading L-polynomial coefficient.
The test reads only this fixture and imports no enumeration code from the F_2 census table repository.

To recreate the sample, use a local checkout of the F_2 census table repository
containing the pinned revision. Run from this repository's root:

```sh
python3 tests/support/extract_published_f2.py /path/to/f2-census-table artifacts/sample.csv
cmp artifacts/sample.csv tests/fixtures/published-f2/sample.csv
```

The extractor reads the pinned Git object, scans the archive once, checks its
observed class sizes and writes exactly 16 rows. It uses Python's standard
library and does not fetch anything. Store its output under an ignored scratch
directory when comparing it locally.

Run `make check/published_f2.m` to compare it with the package.
