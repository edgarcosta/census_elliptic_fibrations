# Heavy reference checks

Thanks for helping check the census exhaustively. Run these checks from
the repository root with Magma attached through the existing launcher:

```sh
make check-heavy REFERENCE_PROCESSES=8
```

The target runs every `reference/heavy/check_*.m` through the same launcher
as `check-reference`. A disagreement prints `MISMATCH` and exits nonzero.
Every child must exit zero and write the completion marker. Checks are
silent on success. Temporary files stay under `artifacts/` and listing
files are removed after validation. Both `make check-reference` and
`make check-heavy` need Magma, a POSIX shell, standard shell utilities, and
`python3` for `scripts/run_reference_check.sh` to write its Magma launcher.
The fast `Tests/reference_runner.m` test uses that launcher and has the same
requirements. No C++ compiler is needed.

| Check | Coverage |
| --- | --- |
| `check_recount.m` | Independent scan against all six package cells on all 88 F2 and 1,860 F3 prefixes; bijective prefix matching and complete scans |
| `check_listing.m` | Every representative in every prefix and all six cells, including empty files; counts, distinct integer keys, key round trips, properties, and completion metadata |
| `check_f2_global.m` | Full-group Burnside on long models with free a2; prefix totals and singular cells from permutation-group orbits |
| `check_f3_global.m` | Completed-square full-group Burnside, prefix totals and singular cells from permutation-group orbits |
| `check_prefix_f2.m`, `check_prefix_f3.m` | Permutation-group orbits on complete coefficient fibres for ten F2 and fifteen F3 target prefixes, with an independent long-model classifier |
| `check_equation_rows.m` | Six complete quotient-equation enumerations against public equation rows, with independently enumerated prefix and translation weights |
| `check_coefficients_f2.m` | Complete nonminimal coefficient union and discriminant-solved nonregular locus; exhaustive cusp coefficients |

`REFERENCE_PROCESSES` controls the recount, listing, and equation-row child
processes. Other checks run in one process. The target runs checks in
sequence, so these process counts do not accumulate. Permutation-group
checks can use several GB of memory per process.

The original full listing of 1,309,216,705 keys took 27,533 seconds
(about 7.6 hours) on 24 workers with Magma 2.30-1. The independent recount
took about 35 seconds for F2 on eight workers and 734 seconds for F3 on
24 workers. These are measurements of the original scripts, before this
port. The new Magma F2 coefficient check took 97 seconds. The additional
permutation-group checks have no measured full-target run time; some F3
targets can take hours.

For a short smoke run of the two exhaustive checks, select a few public
prefix indices and the two check files:

```sh
make check-heavy REFERENCE_PROCESSES=2 HEAVY_PREFIX_INDICES=1,2 \
  HEAVY_FILES='reference/heavy/check_recount.m reference/heavy/check_listing.m'
```

`HEAVY_PREFIX_INDICES` is a comma-separated list valid in both fields.
It restricts only the recount and listing. `HEAVY_TARGET_LIMIT=1` restricts
each permutation-group prefix check to its first target. Unset these
variables and leave `HEAVY_FILES` at its default for the complete target.

The shorter coefficient check over F3, the equation-weight checks on every
prefix, and the engine identity-locus checks belong to `check-reference`.
The engine check inspects its eight Boolean locus patterns; the public
equation-row check independently tests the corresponding six cells.
The formulas themselves are already checked by `Tests/equation_counts.m`.

Distinct integer keys do not independently establish that two models lie
in different coordinate-change orbits. `check_prefix_f2.m` and
`check_prefix_f3.m` compare per-cell permutation-group orbit counts with
`NumberOfFibrations`; they do not inspect listed representatives. All checks
still share Magma arithmetic and the stated Weierstrass equivalence and
minimality criteria.
