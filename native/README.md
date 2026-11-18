# Native command

Thanks for helping develop and check the native command. Its entry point is
`main.cpp`. Key decoding, classification, and local reduction are in
`models/`; exact trace tables and the twist sign are in `lfunction/`.
The [model guide](../docs/native/models.md) describes its input and output.

From the repository root, build and run it with:

```sh
make census-native
artifacts/native/census-native --help
```

The build needs GNU Make 4.2 or later, a C++17 compiler, NTL and GMP.
On Debian or Ubuntu, install `g++`, `make`, `python3`, `libntl-dev` and
`libgmp-dev` for the build and its checks.

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

## Compiler and library paths

The Makefile accepts `CXX`, `CPPFLAGS`, `CXXFLAGS`, `LDFLAGS`, `LDLIBS`, and
`NATIVE_LDLIBS` as command-line overrides. `CXXFLAGS` defaults to
`-O2 -Wall -Wextra -Wpedantic`; the build adds `-std=c++17`.
`NATIVE_LDLIBS` defaults to `-lntl -lgmpxx -lgmp`.

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

After changing build flags or paths, remove `artifacts/native/census-native`
or use `make -B census-native` to force a rebuild.
