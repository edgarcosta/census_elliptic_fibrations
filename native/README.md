# Native command

Thanks for helping develop and check the native command. Its entry point is
`main.cpp`. Key decoding, classification, and local reduction are in
`models/`; exact trace tables, the twist sign and signed L-polynomials are in
`lfunction/`.
Complete bounded torsion, derived invariants, and Sha/regulator candidates
are in `invariants/`. The exact CSV codec, zstd writer, and shared output
helpers are in `records/`.
The [model guide](../docs/native/models.md) describes its input and output.

From the repository root, build and run it with:

```sh
make census-native
artifacts/native/census-native --help
```

The build needs GNU Make 4.2 or later, a C++17 compiler, NTL, FLINT, GMP, and zstd.
On Debian or Ubuntu, install `g++`, `make`, `python3`, `libntl-dev`,
`libflint-dev`, `libgmp-dev`, `libzstd-dev`, and `zstd` for the build
and its checks. The `zstd` package also provides `zstdcat` for inspecting output.

`help`, `--help`, and `-h` print help to standard output and exit with status
0. Missing commands, unknown commands, and extra help arguments print usage
to standard error and exit with status 2.

Decode keys and compare the native and Magma model computations with:

```sh
printf '2 1346371584\n3 951893735796\n' | artifacts/native/census-native invariants --diagnostic
make check/native_model.m check/native_model_magma.m
```

Run the build and command checks through the Magma runner:

```sh
make check/native_build.m
magma -bn filename:=tests/native/native_build.m tests/run_tests.m
```

These checks also need Magma and Python 3. `make check` runs them with the
[package tests and examples](../README.md#tests-and-examples).

## Trace tables and point counts

The [table guide](../docs/native/tables.md) describes the equations, field
conventions and independent checks. From the repository root:

```sh
printf '3 2 0\n3 2 zero\n' | artifacts/native/census-native count
make check/native_tables_magma.m check/native_tables.m
```

## Local fibres and the twist sign

The [local guide](../docs/native/local.md) describes the fibre records,
wild conductors and minimized Frobenius twists. From the repository root:

```sh
printf '2 4884480\n3 372699133912\n' | artifacts/native/census-native invariants --diagnostic
make check/native_local_magma.m check/native_local.m
```

The CSV appends local and twist fibres, L-polynomial degree, sign, heights
and the Tamagawa valuation difference to the model columns.
Unsupported census cells have empty local and sign fields.

## Complete torsion

The [torsion guide](../docs/native/torsion.md) describes the group factors,
point encodings, and complete K3 search bound. From the repository root:

```sh
printf '2 2151677952\n' | artifacts/native/census-native invariants --diagnostic
make check/native_torsion_magma.m check/native_torsion.m
```

The CSV appends `torsion`, `torsion_points`, and `torsion_orders` to the
local and sign columns. Unsupported census cells have empty torsion fields.

## Signed L-polynomials

The [L-polynomial guide](../docs/native/lfunction.md) describes the
coefficient encoding, the counting backends and the Magma comparisons.
From the repository root:

```sh
make census-native census-native-scalar
printf '2 4884480\n3 372699133912\n' | artifacts/native/census-native invariants --diagnostic
printf '2 0 0 1 64 0\n' | artifacts/native/census-native invariants --diagnostic --coefficients --backend scalar
make check/native_lfun_magma.m check/native_lfun.m
```

The CSV appends `lower_L`, `L_polynomial`, and `trace_sums` to the torsion
columns. Unsupported census cells have empty L-polynomial fields.
`--backend auto|scalar|avx2` selects the point-counting backend, and
`artifacts/native/census-native-scalar` is built without the AVX2 code.
`--coefficients` reads base-q coefficient codes `q a1 a2 a3 a4 a6` instead
of `q key`.

## Derived invariants and candidates

The [invariant guide](../docs/native/invariants.md) gives the exact K3
formulas, field meanings, and supplied-data interfaces. The
[candidate guide](../docs/native/candidates.md) gives the necessary
constraints on Sha orders and regulators. From the repository root:

```sh
make census-native
printf '2 4884480\n3 372699133912\n' | artifacts/native/census-native invariants --diagnostic
make check/native_derived_magma.m check/native_derived.m
make check/native_candidates_magma.m check/native_candidates.m
```

The CSV appends derived fields and aligned candidate lists to the
L-polynomial columns. A singleton supplies the individual Sha order,
regulator, Brauer order, and absolute base-field NS discriminant.
An unresolved split leaves those four individual fields blank.
Unsupported census cells have empty derived and candidate fields.

## Compressed invariant records

The [record guide](../docs/native/records.md) defines the exact CSV schema,
status meanings, and Magma codec. The default output is one zstd frame at
level 9 with a checksum. From the repository root:

```sh
printf '2 1077936128\n3 320009441124\n' | artifacts/native/census-native invariants > artifacts/rows.csv.zst
zstdcat artifacts/rows.csv.zst
make check/native_records_magma.m check/native_records.m
make check-native-invariants
```

The frame begins with the 18-column header. Unknown arithmetic stays empty;
known empty lists use `[]`. Every admissible Sha order is retained.
`--diagnostic` selects the expanded, uncompressed CSV used above.
Malformed input ends the stream after its preceding records, with a
line-numbered error. A computation failure retains completed arithmetic
in an `error` row and continues; either failure gives a nonzero exit status.
The native computation starts no Magma process.

## Compiler and library paths

The Makefile accepts `CXX`, `CPPFLAGS`, `CXXFLAGS`, `LDFLAGS`, `LDLIBS`, and
`NATIVE_LDLIBS` as command-line overrides. `CXXFLAGS` defaults to
`-O2 -Wall -Wextra -Wpedantic`; the build adds `-std=c++17`.
`NATIVE_LDLIBS` defaults to `-lntl -lflint -lgmpxx -lgmp -lzstd`.

For headers or libraries installed outside the compiler's default paths,
pass include and library search paths explicitly:

```sh
make census-native CXX=g++ CPPFLAGS=-I/opt/local/include LDFLAGS=-L/opt/local/lib
```

With GCC, `CPATH` and `CPLUS_INCLUDE_PATH` add header search directories,
and `LIBRARY_PATH` adds library search directories at link time. On Linux,
`LD_LIBRARY_PATH` adds shared-library search directories when the command
runs. These variables take colon-separated paths. Set them in the calling
shell when using a local library installation.

After changing build flags or paths, remove the binaries in `artifacts/native/`
or use `make -B census-native census-native-scalar` to force a rebuild.

## Complete F_2 comparison

The [verification guide](../docs/native/verification.md) describes the complete
F_2 population, coordinate-orbit matching, identities, and published columns.
This check also needs Python 3.11 or later. From the repository root:

```sh
git lfs pull --include='tests/fixtures/full_f2/census.csv.zip' --exclude=''
make check-native-invariants-f2
```

The comparison writes its summary, matches, missing rows, disagreements, and
logs under `artifacts/full-f2/`. CI retains those files with a 30-minute limit.
