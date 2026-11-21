# Native command

From the repository root, build and run the native command with:

```sh
make census-native
build/native/census-native --help
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
printf '2 1346371584\n3 951893735796\n' | build/native/census-native invariants --diagnostic
make check/native_model.m
```

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

Executables are built under `build/`: commands in `build/native/` and test
programs in `build/tests/native/`. Test scratch and comparison output also
use `build/`. After changing build flags or paths, run `make clean` to remove
`build/`, or `make -B census-native` to force a rebuild.
