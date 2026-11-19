# Separate refinement records

Thanks for helping preserve the census evidence. The F_2 recovery command
writes `census-refinement` version 1 records beside the unchanged source rows.
It implements the retained-basis part of the bounded stage.

Each record identifies `q`, a decimal-string `key`, `key_encoding_version`,
`stage`, `input_ordinal`, and `producer_id`. The ordinal belongs to the source
shard named in `input_identity.receipt`. The source receipt digest, compressed
shard digest, key-file digest and retained package source identity identify the
original population. These are retained-run receipts, not production stage 1
CSV receipts.

`claims.mw_basis` has its own status, method, reason and value. An exact value
contains four-code rational-function coordinates, the row-major upper triangle
of the rational Gram matrix, torsion points, regulator, search bound and covering bound. It also records independent division results. Integer codes and reduced rationals
are decimal strings; an empty point denotes the identity. An unresolved value
is null, with the partial native result retained in its evidence file.

`producer_id` is the SHA-256 of canonical JSON for `import-producer.json`:
sorted keys, compact separators and a final newline. That manifest identifies
the prepared inputs, importer, auditor, replay producer and audit producer.
Replay pins the executable and all nine kernel dependencies. Audit pins its
Magma executable, version and script. Evidence paths refer to the run directory.
Receipts bind the request, raw result, transcript and selected producers.

Import rejects changed dependencies, duplicate or wrong-key results, invalid
Gram matrices, determinant/split disagreement and contradictory exact bases.
Missing or uncertified results remain unresolved. Basis promotion requires
substitution, `HeightPairingMatrix`, torsion, local-bound, covering and division
checks by Magma. A doubled point supplies a positive division control.

Recovery records contain only the MW basis claim. They add no split or Sha[2]
claims and never overwrite the retained source rows. Existing exact basis
values cannot be replaced by an incomplete retry. `import.json` reports counts,
the records digest and zero new splits credited. Later bounded, structure and
d2 producers require their own implementations and contracts.

The accessible-basis producer uses the same basis value and independent audit.
Its producer manifest pins the search scripts and prepared cohort; each batch
pins replay and Magma producers through the retained recovery receipts.
For a previously ambiguous split, a certified determinant selects a remaining
candidate pair. The record then also carries `claims.regulator` and
`claims.sha_order`, each with exact status, method and decimal-string value.
The original split and Sha[2] claims remain in their source overlay.
`verification.json` counts these new splits separately from basis gains.
