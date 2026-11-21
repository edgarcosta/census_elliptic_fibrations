# Saved F_2 Sha[2] import

The public command imports the completed ordinary-F_2 refinement and verifies
its saved evidence:

```sh
python3 production/refine.py import-sha2 --output artifacts/f2-sha2 &&
python3 production/refine.py verify-sha2 --output artifacts/f2-sha2
```

Run from this checkout with Python 3 and the Git LFS data fetched.
The commands use only the standard library and perform no residue searches.
They validate the original input and artifact hashes, all successful pairing
checkpoints, local matrix sums, and the independent 64-row residue sample.
Failed and timed-out attempts remain evidence; they supply no exact claims.
The commands reject changed inputs, conflicting pairings and empty candidate
intersections. Repeating an import accepts identical output; changed output
requires investigation and a separate directory.

The output contains 6,008 separate refinement records, the complete evidence
archive, its manifest, an importer identity and a summary. It preserves 3,418
settlements and 43 shortened unresolved lists. The saved population then has
1,395,008 proved splits and 2,590 open splits. These are imported results;
`new_arithmetic` is zero. The published table and retained source rows are
separate inputs to later release work.

`verify-sha2` reconstructs every record from the retained evidence and compares
all output bytes, including producer identities. Verify before consuming
`refinements.jsonl`. `pipeline.records.compose_sha2` composes a verified record
with a retained row by field, key, coefficient codes, rank and Sha*Reg.
It preserves generators and other source quantities, intersects paired Sha
orders and regulators, and rejects conflicting exact split or dimension claims.
An inapplicable or missing descent cannot erase a prior exact dimension.

See [the record contract](../data/refinements-v1.md) for claim statuses and
[the evidence inventory](../data/sha2/f2-ordinary-refinements/README.md) for arithmetic limits.
`make check/sha2_import.m` exercises the public commands and composition.
