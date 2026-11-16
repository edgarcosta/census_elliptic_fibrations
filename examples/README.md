# Package examples

This directory holds runnable examples of counting and listing fibrations
over GF(2) and GF(3).

- `counts.m`: counts all model orbits over both fields and checks the totals.
- `representatives.m`: lists quasi-elliptic K3 representatives, decodes a key and checks its properties.

From the repository root, run the two examples separately:

```sh
magma -bn filename:=examples/counts.m run_tests.m
magma -bn filename:=examples/representatives.m run_tests.m
```

Both examples run in CI and in `make check`. That target also runs the
tests, using the processor count reported by `nproc` (one if unavailable).
See the [interface](../docs/interface.md) for selections and return values,
and the [parent README](../README.md) for installation and printed examples.
