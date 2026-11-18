# Package source

Thanks for helping develop the census. This directory implements the Magma
intrinsics and their private orbit, classification, and file routines.

- `enumeration/declarations.m`: the `FibrationCensus` verbosity levels.
- `models/forms.m`: binary-form bases, coefficient vectors, symmetric-power actions and power tests.
- `models/models.m`: input checks, coordinate changes, normalization, discriminants, and key conversion through `FibrationKey` and `aInvariantsFromKey`.
- `enumeration/engine.m`: affine translation quotients, invariant locus counts and representative enumeration.
- `enumeration/f2.m`, `enumeration/f3.m`: characteristic-specific prefix tables and affine fibres with their three loci.
- `models/minimality.m`: local congruences, nonminimal loci, bad-place certificates and `IsEverywhereMinimal`.
- `models/classification.m`: generic-fibre kinds, singular/nonregular loci and `FibrationProperties`.
- `enumeration/census.m`: public prefix, orbit-count, equation-count and representative intrinsics with selections.
- `models/local.m`: bad-place records, full conductors and minimized Frobenius twists through `FibrationLocalData`.
- `lfunction/trace_tables.m`: exact ordinary trace tables through `FibrationTraceTable`.
- `invariants/torsion.m`: complete torsion factors and point sets through `FibrationTorsion`.
- `records/files.m`: binary representative files, CSV summaries and completion metadata through `WriteFibrationRepresentativeKeys`.

From Magma in the repository root, attach the package:

```c++
AttachSpec("spec");
```

The root `spec` includes `magma/spec`, whose entries are relative to this
directory. Run the tests from the repository root:

```sh
magma -bn tests/run_tests.m
```

CI runs the tests and examples. See the
[interface](../docs/interface.md) for calls and return values, the
[mathematics](../docs/mathematics.md) for the loci and quotient construction,
and the [parent README](../README.md) for installation.
