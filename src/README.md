# Package source

This directory implements the Magma intrinsics and their private orbit,
classification, and file routines.

- `declarations.m`: the `FibrationCensus` verbosity levels.
- `forms.m`: binary-form bases, coefficient vectors, symmetric-power actions and power tests.
- `models.m`: input checks, coordinate changes, normalization, discriminants, and key conversion through `FibrationKey` and `aInvariantsFromKey`.
- `engine.m`: affine translation quotients, invariant locus counts and representative enumeration.
- `f2.m`, `f3.m`: characteristic-specific prefix tables and affine fibres with their three loci.
- `minimality.m`: local congruences, nonminimal loci, bad-place certificates and `IsEverywhereMinimal`.
- `classification.m`: generic-fibre kinds, singular/nonregular loci and `FibrationProperties`.
- `census.m`: public prefix, orbit-count, equation-count and representative intrinsics with selections.
- `files.m`: binary representative files, CSV summaries and completion metadata through `WriteFibrationRepresentativeKeys`.

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
