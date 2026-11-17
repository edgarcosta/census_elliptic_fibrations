# Independent census checks

Thanks for helping check the census. These tests compare the public interface
with independent calculations and data:

- `burnside.m` derives the whole-space affine action by direct substitution
  over both fields. It checks the total orbit count, including singular and
  nonminimal equations, but does not separate the six cells. See the
  [group action and Burnside argument](../docs/mathematics.md#whole-space-burnside-oracle).
- `equation_counts.m` checks all six normalized-equation counts against closed
  formulas. `equation_weights.m` checks prefix stabilizers, translation images
  and per-prefix weights. Divisibility by the translation-image size does not
  rule out quotienting by a larger shift subspace.
- `archive_fixtures.m` and `archive_invariants.m` compare coordinate-orbit
  multiset hashes for all 88 F_2 prefixes and two F_3 prefixes. The F_3 sample
  omits smooth nonminimal models and some stabilizer types; archived singular
  flags do not distinguish the four singular cells. See
  [archive provenance and coverage](data/oldcode-archive/README.md).
- `published_f2.m` checks coordinate-orbit matches and independent invariants
  for 16 rows from the incomplete published F_2 table. It does not establish
  complete-census coverage. See [fixture provenance](data/published-f2/README.md).
- `minimality_tate_f2.m` and `minimality_tate_f3.m` compare complete certificates
  with Magma's Tate algorithm in both integral charts on 256 reproducible
  smooth models per field. Explicit fixtures cover translated nonminimal
  equations, multiple bad places and wild reduction; this oracle excludes
  singular models.
- `invariance.m` checks public properties, transported certificates and decoded
  normalization under coordinate changes. `orbit_listing.m` checks disjoint,
  exhaustive representative orbits under independent substitutions on six
  selected fibres. These do not cover all prefixes.

The remaining tests check classification fixtures, counts and selections,
coefficient encodings, coordinate substitutions, affine loci, minimality,
and saved files. `published_counts.m` checks the printed census numbers;
`runner.m` checks shell exit status for successful and failing assertions.
Underscore-prefixed files provide shared helpers. See the
[root README](../README.md#tests-and-examples) for test commands and prerequisites.

CI runs each test file as its own job with a 30-minute limit.

## Native tests

Run the native checks from the repository root:

```sh
make -j2 check-native-invariants
```

The full suite also runs these checks. The top-level `native_*.m` files
call the test programs in `tests/native/` through the shared runner.
See the root README for native build dependencies and commands for
running one test file.

- `native_build.m` checks native build dependencies, command help and
  invalid invocations.
- `native_model.m` compares native coefficients, both charts, six-cell
  labels and invalid inputs with Magma.
- `native_model_magma.m` checks fixed model keys, both charts, six-cell
  classification and encoding bounds.
- `native_local_magma.m` checks Tate records on both charts, split and
  nonsplit fibres, wild conductors and minimized Frobenius twists.
- `native_torsion_magma.m` compares complete bounded torsion groups,
  point orders, rank zero and infinity with Magma.
- `native_torsion.m` compares native torsion points, orders and invariant
  factors with the same Magma sample on both charts.
- `native_derived.m` compares exact native K3 fields and identities with
  Magma, including constant-field integer controls.
- `native_candidates.m` compares complete native Sha and regulator lists
  with Magma, including singleton-only values and exact supplied-data controls.
- `native_tables_magma.m` checks exact ordinary tables in both
  characteristics through degree 12, direct small-field counts and stable hashes.
