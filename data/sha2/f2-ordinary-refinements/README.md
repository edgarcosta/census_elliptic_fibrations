# Ordinary-F_2 Sha[2] refinement evidence

Thanks for preserving these completed computations. `evidence.tar.gz` contains
the 2026-11-17 exact-refinement inputs and results, with a byte inventory in
`manifest.json`. Both files use Git LFS. The archive preserves original bytes;
its timestamps and ownership are normalized for reproducible packaging.

The source revisions recorded in `inputs/effect_summary.json` are calibration
`45a0aa6f2fa7a006095da0201f77ced69cc70dfe` and note `3a9579d`.
The archived pipeline differs from that calibration revision; `sources.json`
and the inventory pin its exact bytes.

The inventory includes the original unresolved-row snapshot, public-row
comparison, source note and pipeline, and all three result CSVs.
It preserves every saved checkpoint and transcript, the computation scripts
and the mathematical review. It also includes the independent bound check
and the independent residue sample at two precisions. The original manifests retain their historical paths as provenance.
The importer resolves these files inside the archive and never executes them.
No installed private code or external arithmetic program is required.

The source population has 1,397,598 rows and 6,008 unresolved splits.
The import retains 1,191 bound settlements and 2,227 residue settlements,
leaving 2,590 unresolved splits: 1,004 ordinary and 1,586 supersingular.
All 2,265 residue targets completed; 38 retain two candidates.
Five additional bound refinements retain two candidates, making 43 shortened
unresolved lists. There is no credit for new arithmetic or new residue work.

Reproduction of the original archive uses the saved `sha2-exact` directory
and the five source files identified by its `sources.json`.
Select every file in its `manifest.json`. Also select the original manifest,
sources, rows, selection, summary and conflicts. Keep the review, residue-check
review, crosscheck review, validation, validation sources, release validation,
diversity and method-check selection. Include every file under `checkpoints`, `logs`,
`residue-check` and `crosscheck`, excluding Python bytecode caches.
Store the five sources under `inputs/<basename>`. Sort archive member paths,
use regular files with mode 0644, zero owner/group identifiers and timestamps,
and gzip with an empty original filename and timestamp zero.
The JSON manifest records every member's SHA-256 and byte size, the archive
SHA-256, and the population counts used by the importer.

The saved reviews establish the arithmetic scope. Import checks byte integrity,
associations, matrix ranks and local sums; it does not reprove local-image
completeness or the pairing theorem. Both residue computations share Magma's
Selmer-group implementation. The independent-order controls do not cover a
nonzero defect with quadratic a1. Preserve `review.md` and `residue-check.md`
when using the imported claims; precision agreement alone is insufficient.
