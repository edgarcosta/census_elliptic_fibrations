# Contributing to the elliptic K3 census

Thank you for helping develop and check the census. These conventions keep
the mathematics, its implementations and the published data consistent.

Read the README for installation and test commands, `magma/README.md` for the
package layout, and the relevant guide under `docs/` before changing behavior.
Keep details in their owning documents and update them with the implementation.
`tests/README.md` describes the independent checks and their limits.

## Implementations

Implement the mathematics in Magma first. Native implementations need a
Magma implementation for comparison. Keep the equivalence relation and the
separate generic-fibre and minimality classifications in `docs/mathematics.md`
consistent with the code. When changing an algorithm, establish its correctness
by testing it and by reviewing its mathematical preconditions, normalizations
and termination.

Public Magma operations are intrinsics named for their results, using
`UpperCamelCase` and established mathematical spelling such as
`aInvariantsFromKey`. Use `NumberOf<Thing>` for counts and constructor
overloads or `<Result>From<Source>` for conversions. Private functions and
procedures use `snake_case`.

Put essential arguments before the colon and named options after it.
An intrinsic docstring starts with a noun phrase for the result, has no
period at the end, and names the default of every option.

Magma code has one space after commas and around `:=`, and none
immediately inside brackets. Spacing changes keep line breaks, indentation
and comment placement. Fence Magma examples in Markdown as `c++`.

Standalone Magma scripts begin with `SetQuitOnError(true);`, end with `quit;`,
and run with `magma -bn`. Tests and examples loaded by `run_tests.m` leave
error handling and process exit to that runner.

## Tests and data

Tests live in `tests/`; native tests belong in `tests/native/`.
The runner discovers only top-level `.m` files, so integrate native tests
with the runner and CI. Shared helpers have names starting with `_`,
which excludes them from test discovery.

Test census mathematics or real failures. Use independent formulas,
substitutions, algorithms or published data as oracles. Name each scenario with
a procedure, assert on its results and keep successful scenarios silent. Run the
test commands in `README.md` from the root of the checkout under test, with
that checkout's `spec` as the only copy of the package attached.

Store census data and fixtures with Git LFS, including the matching paths in
`.gitattributes`. Record the source revision, extraction procedure and
comparison limits beside each fixture. A file checksum establishes byte
integrity; a mathematical comparison establishes what the data represents.

## Project wording

The census is one joint project: name data and code by what they contain
or where they were published, such as the published F_2 census table.
