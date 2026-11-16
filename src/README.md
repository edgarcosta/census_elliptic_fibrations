# Package source

Thanks for helping develop the census. This directory implements the Magma
intrinsics and their private routines.

## Enumeration

- `declarations.m`: the `FibrationCensus` verbosity levels.
- `forms.m`: binary-form bases, coefficient vectors, symmetric-power actions and power tests; native power tests are in `native/model.cpp`.
- `models.m`: input checks, coordinate changes, normalization, discriminants, and key conversion through `FibrationKey` and `aInvariantsFromKey`; native counterpart: `native/model.cpp`.
- `engine.m`: affine translation quotients, invariant locus counts and representative enumeration.
- `f2.m`, `f3.m`: characteristic-specific prefix tables and affine fibres with their three loci.
- `minimality.m`: local congruences, nonminimal loci, bad-place certificates and `IsEverywhereMinimal`; native counterpart: `native/model.cpp`.
- `classification.m`: generic-fibre kinds, singular/nonregular loci and `FibrationProperties`; native per-model classification is in `native/model.cpp`.
- `census.m`: public prefix, orbit-count, equation-count and representative intrinsics with selections.
- `files.m`: binary representative files, CSV summaries and completion metadata through `WriteFibrationRepresentativeKeys`.

## Per-surface arithmetic

- `trace_tables.m`: finite fields and exact ordinary trace tables; native counterparts: `native/field.cpp` and `native/lfun.cpp`.
- `local.m`: Tate fibre records, component Frobenius cycles and minimized twists; native counterpart: `native/model.cpp`.
- `lfun.m`: signed L-polynomials from fibre traces and the functional equation; native counterpart: `native/lfun.cpp`.
- `torsion.m`: complete bounded zero-height sections and torsion invariant factors; native counterpart: `native/invariants.cpp`.
- `invariants.m`: ranks, cohomology, point counts, formal height and BSD/Artin-Tate products; native counterpart: `native/invariants.cpp`.
- `candidates.m`: admissible Sha orders and matching regulators; native counterpart: `native/invariants.cpp`.
- `records.m`: invariant rows and CSV encoding; native counterparts: `native/records.cpp` and `native/format.cpp`.

From Magma in the repository root, attach the package:

```c++
AttachSpec("spec");
```

The root `spec` lists these modules. Run the tests from the repository root:

```sh
magma -bn run_tests.m
```

CI runs the tests and examples. See the
[interface](../docs/interface.md) for calls and return values, the
[mathematics](../docs/mathematics.md) for the loci and quotient construction,
and the [parent README](../README.md) for installation.
