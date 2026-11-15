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
`native_build.m` checks native build dependencies, command help and invalid invocations.
Underscore-prefixed files provide shared helpers. See the
[root README](../README.md#tests-and-examples) for test commands and prerequisites.
