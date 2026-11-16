# AGENTS.md

An untracked `AGENTS.local.md` in this directory, when present, adds setup
for one machine; read both files. This file points agents using any harness
to the repository conventions.

## Where the conventions live

Read the relevant document before changing the behavior it describes. Keep
details in their owning documents and update them with the implementation.

- [README.md](README.md) - installation, scope, counts, and test commands.
- [docs/interface.md](docs/interface.md) - public calls, selections, and certificates.
- [docs/mathematics.md](docs/mathematics.md) - equivalence, counting, and classification.
- [docs/keys.md](docs/keys.md) - coefficient encodings and representative files.
- [src/README.md](src/README.md) and [spec](spec) - package modules and attachment.
- [tests/README.md](tests/README.md) - independent checks and their limits;
  fixture READMEs under `tests/data/` record provenance and reproduction.
- [.github/workflows/tests.yml](.github/workflows/tests.yml) and
  [Makefile](Makefile) - CI discovery and commands runnable locally.

## Working in this checkout

- Implement the mathematics in Magma first. A native implementation also
  needs a Magma implementation for comparison.
- Keep the equivalence relation and the separate generic-fibre and minimality
  classifications in `docs/mathematics.md` consistent with the code.
  Review mathematical preconditions, normalizations and termination when
  changing an algorithm; passing tests alone do not establish correctness.
- Expose public operations as intrinsics named for their results, using
  `UpperCamelCase` and established mathematical spelling such as
  `aInvariantsFromKey`. Use `NumberOf<Thing>` for counts and constructor
  overloads or `<Result>From<Source>` for conversions. Private helpers are
  `snake_case` functions or procedures.
- Put essential arguments before the colon and named options after it.
  An intrinsic docstring starts with a noun phrase for the result, has no
  period at the end, and names the default of every option.
- For Magma whitespace edits, use one space after commas, none immediately
  inside brackets, and spaces around `:=`. Preserve line breaks, indentation,
  and comment placement. Fence Magma examples in Markdown as `c++`.
- Attach this checkout's `spec` and run tests from its root using the README
  commands. Standalone Magma scripts begin with `SetQuitOnError(true);`
  and finish with `quit;`; run them with `magma -bn`.
  Files loaded by `run_tests.m`, including tests and examples, leave error
  handling and process exit to that runner.
- Put tests in lowercase `tests/`; native tests belong in `tests/native/`.
  The current runner discovers top-level `.m` files only. Add the runner and
  CI integration when introducing native tests. Shared test helpers start
  with `_`, which excludes them from test discovery.
- Tests check census mathematics or real failures. Prefer independent formulas,
  substitutions, algorithms or published data to expected values produced by the code under test.
  Use named procedures that identify each scenario and assertions on its
  results. Successful test scenarios are silent.
- Put census data and data fixtures in Git LFS, adding the matching path to
  `.gitattributes`. Record the source revision, extraction procedure and
  comparison limits beside a fixture. Distinguish a file checksum from a
  mathematical comparison of its contents.
- Treat the census as one joint project. Name the data or code being
  discussed, such as the published F_2 census table.
