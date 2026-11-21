# Package source

This directory implements the Magma intrinsics and their private orbit,
classification, and file routines.

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
- `lfunction/lfun.m`: signed L-polynomials from fibre traces and the functional equation through `FibrationLPolynomial`.
- `invariants/torsion.m`: complete torsion factors and point sets through `FibrationTorsion`.
- `invariants/invariants.m`: exact ranks, cohomology, formal height and BSD/Artin-Tate products through `FibrationInvariants`.
- `invariants/candidates.m`: admissible Sha orders and corresponding regulators through `FibrationCandidates`.
- `records/records.m`: exact invariant rows and versioned CSV through `FibrationInvariantRow` and `FibrationInvariantCSV`.
- `records/files.m`: binary representative files, CSV summaries and completion metadata through `WriteFibrationRepresentativeKeys`.

See the [root README](../README.md) for installation and test commands,
the [interface](../docs/interface.md) for calls and return values, and the
[mathematics](../docs/mathematics.md) for the loci and quotient construction.
