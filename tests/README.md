# Independent census checks

These tests compare the public interface with independent calculations and data:

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
  [archive provenance and coverage](fixtures/oldcode-archive/README.md).
- `published_f2.m` checks coordinate-orbit matches and independent invariants
  for 16 rows from the incomplete published F_2 table. It does not establish
  complete-census coverage. See [fixture provenance](fixtures/published-f2/README.md).
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
`native/native_build.m` checks the native build, help aliases, rejected
invocations, build ordering before native tests, and the guard against duplicate
test basenames.
`compare/native_model.m` compares native coefficients, both charts, six-cell
labels and rejected input with Magma. The key tests pin monomials independently.
`trace_tables.m` checks all table digests and independent low-degree traces;
`compare/native_tables.m` checks native digests and interleaved requests.
`native_local_magma.m` checks reciprocal charts, wild fibres and twist controls.
The model comparison also checks local records and signs.
`native_lfun_magma.m` compares with an unattached `LFunction` oracle;
`compare/native_lfun.m` compares package results across counting backends.
`native_derived_magma.m` checks nonprime-field controls; derived and candidate
comparisons check the command and supplied candidate engine.
The record tests compare compressed bytes and exercise codec validation.
Shared helpers live in `support/` and fixtures in `fixtures/`.
The runner, Make, and CI use `python3 tests/list_tests.py` to list test
entries in `tests/`, `tests/native/`, and `tests/compare/`,
excluding the runner and helpers. See the [root README](../README.md#tests) for test
commands and prerequisites.

Native tests also need the [native build prerequisites](../native/README.md).

CI runs each test file as its own job with a 30-minute limit.
