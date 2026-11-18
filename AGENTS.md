# AGENTS.md

An untracked `AGENTS.local.md` in this directory, when present, adds setup
for one machine; read both files. This file points agents using any harness
to the repository conventions.

## Where the conventions live

- [CONTRIBUTING.md](CONTRIBUTING.md) - implementation, test, data and writing rules.
- [README.md](README.md) - installation, scope, counts, and test commands.
- [docs/interface.md](docs/interface.md) - public calls, selections, and certificates.
- [docs/mathematics.md](docs/mathematics.md) - equivalence, counting, and classification.
- [docs/keys.md](docs/keys.md) - coefficient encodings and representative files.
- [magma/README.md](magma/README.md) and [spec](spec) - package modules and attachment.
- [tests/README.md](tests/README.md) - independent checks and their limits;
  fixture READMEs under `tests/data/` record provenance and reproduction.
- [.github/workflows/tests.yml](.github/workflows/tests.yml) and
  [Makefile](Makefile) - CI discovery and commands runnable locally.

## Working in this checkout

Read and follow [CONTRIBUTING.md](CONTRIBUTING.md) before changing code,
tests, data or documentation; it holds the conventions shared with human
contributors.
