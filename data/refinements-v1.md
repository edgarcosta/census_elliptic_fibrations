# Separate refinement records

Thanks for helping keep mathematical claims and their evidence together.
The saved Sha[2] importer writes `census-refinement` version 1 records beside
the original census data. Stage 1 CSV v1 is unchanged.

Each record has field `q = 2`, a canonical decimal-string `key`,
`key_encoding_version = 1`, `stage = "bounded"`, an `input_ordinal`, and a
`producer_id`. Integers in candidate orders, coefficient codes and keys are
strings. Regulators are reduced rational strings with positive denominators;
integers omit `/1`. Dimensions and ordinals are JSON integers.

This adapter imports a retained unresolved-row snapshot. Its ordinal belongs
to `inputs/unresolved_raw.json` inside the evidence archive, rather than a
production prefix. `input_identity` names that snapshot, its byte digest,
the canonical source-row digest and the surface's arithmetic identity.
`stage1_receipt_sha256` is null because the saved snapshot is the dependency.
These identities do not claim a new enumeration or a production prefix receipt.

`producer_id` is the SHA-256 of `import-producer.json`, encoded as sorted,
compact JSON with a final newline. The producer identifies the importer and
composition sources, the evidence archive and its manifest by content digest.
Each evidence reference names an archive member and its original byte digest.
Evidence paths are relative to the import directory. Keep the archive,
manifest, producer, records and `import.json` together.

`claims.sha2` describes the dimension of Sha[2], separately from the split.
For a successful ordinary descent it uses `d = f + v - rank - 1 - c_E`.
The defect is even and lies between zero and `min(v, f + v - rank - 1)`.
Without a residue value, all feasible defects remain constraints; a dimension
is exact only if every compatible defect gives the same dimension.
A saved residue pairing supplies its exact defect. Generic supersingular
fibres are `not_applicable`; missing saved descent is `not_computed`.
Neither status asserts dimension zero.

`claims.split` contains the surviving `[Sha order, regulator]` pairs.
Dimension zero permits only odd orders; positive dimension d requires
`2^d` to divide the order. Dimension 2 permits both orders 4 and 16.
The filter preserves Sha*Reg and intersects the supplied candidates.
A singleton split is `exact`; a longer list is `bound`.
A processed record can therefore have stage status `unresolved` while its
Sha[2] dimension is exact. Completion of processing does not certify a basis
or a full Brauer group structure.

The public `pipeline.sha2.refine_sha2` applies these constraints to a retained
row, with an optional proved defect. `pipeline.records.compose_sha2` joins a
verified imported record with later retained evidence. It preserves source
quantities and prior exact dimensions, intersects bounds, and fails on wrong
identities, empty intersections or conflicting exact claims. Its optional
`prior_sha2` argument supplies a dimension claim from another verified producer.
The caller remains responsible for verifying that producer's evidence.

An intersection creates a fresh Sha[2] claim with `method = "intersection"`,
the surviving `candidates`, and a `prior_evidence` list of contributing claims.
Each contributor keeps its original method, candidates and evidence.
Repeated composition extends this list with the original contributors of
each intersection. A single applicable claim is copied unchanged;
`not_computed` and `not_applicable` claims do not constrain the intersection.
