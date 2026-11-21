# Native command

From the repository root, build and run the native command with:

```sh
make census-native
artifacts/native/census-native --help
```

The build needs GNU Make 4.2 or later and a C++17 compiler. It uses only the
C++ standard library. On Debian or Ubuntu, install `g++`, `make`, and
`python3` for the build and its smoke test.

`help`, `--help`, and `-h` print help to standard output and exit with status
0. Missing commands, unknown commands, and extra arguments print usage to
standard error and exit with status 2.

## Compiler and library paths

The Makefile accepts `CXX`, `CPPFLAGS`, `CXXFLAGS`, `LDFLAGS`, `LDLIBS`, and
`NATIVE_LDLIBS` as command-line overrides. `CXXFLAGS` defaults to
`-O2 -Wall -Wextra -Wpedantic`; the build adds `-std=c++17`.
`NATIVE_LDLIBS` is empty for this command.

For headers or libraries installed outside the compiler's default paths,
pass include and library search paths explicitly:

```sh
make census-native CXX=g++ CPPFLAGS=-I/opt/local/include LDFLAGS=-L/opt/local/lib
```

After changing build flags or paths, remove `artifacts/native/census-native`
or use `make -B census-native` to force a rebuild.
