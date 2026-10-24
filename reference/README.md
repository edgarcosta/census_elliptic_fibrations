# Independent reference calculations

Thanks for checking the census against an independent calculation.
These plain Magma files retain the scan-based orbit calculations and
local tests. They are outside the package spec, and production code
does not import them. Package comparisons call the public intrinsics.

## Files

| Path | Calculation |
| --- | --- |
| f2/elliptic_fibration_orbits_f2.m | F2 prefix stabilizers, scanned model orbits and coefficient keys |
| f2/k3_minimality_f2.m | Original smooth F2 minimality filter |
| f2/census_driver.m | Complete smooth F2 scan and minimality counts |
| f3/ef_common.m | Binary forms and coefficient encodings |
| f3/ef3_models.m | F3 triples [C,A,B], coordinate changes and discriminant |
| f3/ef3_orbits.m | F3 prefixes, fibre scans and staged Burnside counts |
| f3/k3_minimality_f3.m | Original smooth F3 minimality filter |
| f3/ef3_stream.m | F3 scan driver and historical output format |
| f3/ef3_burnside.m | Full-group F3 Burnside calculation |
| f3/bench_f3.m | Manual scan and writing benchmark |
| support/classification.m | Direct square/cube tests and local congruences in both charts |
| support/singular.m | Independent singular-orbit representatives and their cell counts |
| checks/ | Comparisons with the package, classification, Burnside, and Tate checks |
| f3/checks/ | Standalone scan, key, group-action and output checks |

The F3 plain files use triples [C,A,B]; the package uses five a-invariants.
Reference and package prefix indices and representatives can differ.
Comparisons match prefix orbits and model orbits; literal representative
keys can differ. The original smooth minimality
filters do not supply the package's all-model certificates.

## Running checks

Run from the repository root. For the complete reference suite use the
Makefile target `check-reference`; for one check use
`reference-check/<basename>.m`. A bounded classification check is:

```sh
make reference-check/check_classification.m
```

The launcher starts a separate Magma process, enables error exits before
loading the check, and verifies a result marker written after the assertions.
Its local files are under the ignored `artifacts/` directory.
The Makefile target is silent on success and prints the log on failure.
This avoids relying on the exit status of an unguarded Magma `load`.

The `check-slow` target runs the wrappers in `examples/slow/`, which launch
their reference checks in child processes. These longer checks are outside
the fast suite and CI matrix. The scan drivers and F3 benchmark are manual;
they are not additional test targets.

## What is independently checked

The singular classification enumerates every discriminant-zero orbit using
the reference transformations. It tests squares and cubes directly and
checks local congruences in both affine charts, including infinity.
It compares both properties and bad-place certificates with the package,
and checks invariance under prefix stabilizers and overlapping bad places.

The total-count check compares independent Burnside counts and singular
cell counts with the public census. The F2 scan comparison covers all
prefixes, including the smooth minimal/nonminimal split. Over F3, all
prefixes are checked for total and smooth counts, while the smooth
minimal/nonminimal split is sampled. Selected representative sets are
compared as orbits, and smooth local tests are compared with Magma's
Tate algorithm. This is not an exhaustive independent recount of the
F3 smooth minimal/nonminimal split.

For the definitions of the six cells, see [the mathematics](../docs/mathematics.md).
For package output, use [the public file format](../docs/keys.md).
