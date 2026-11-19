# Retained F_2 basis replay inputs

Thanks for helping recover the retained basis certificates. These inputs
request 3,085 previously successful native searches over F_2: 1,522 at bound
5, 467 at bound 6 and 1,096 at bound 7. They contain no replay results.

`replay-candidates.jsonl` contains one `(q, key, codes, L, b)` request per
line for `stages-native --bundle`. The input SHA-256 is
`87545c1b886acf925e15ead063373acfdbf40674c145248536c629f3187e39ab`.
Keys and coefficients use the existing kernel request encoding.

`source-manifest.json` identifies the 371 retained shard receipts, their
compressed data digests and the nine kernel dependencies compared with the
reviewed calibration revision
`45a0aa6f2fa7a006095da0201f77ced69cc70dfe`. Rebuild that revision on the target
machine. The manifest's relative paths refer to the original evidence
directory `artifacts/mw/mw-finish-f2/` in the worktree that retained the run.
Resolve these provenance paths against that evidence directory.

The requests were extracted by scanning the saved F_2 arithmetic shards.
A request selects a row whose bound-4 generator record is uncertified but
whose deeper split evidence contains a certified full-rank native basis.
Its exact Gram determinant equals the row's proved regulator. The requests
retain that row's field, key, coefficient codes, L-polynomial and successful
search bound. The input and source manifest are unchanged copies of the
planning evidence.

Source digests establish byte identity. Retained certificate and determinant
checks establish which saved results these requests select. Neither check
replays the searches or supplies an independent arithmetic audit.

The separate refinement records must preserve the existing certified
bases, split claims and exact Sha[2] claims. If replay and the required
independent audit agree, the expected certified basis count is 1,045,076,
up from 1,041,991, with zero new splits credited. The recovery command replays and independently audits these requests before
writing basis-only refinement records.

## Recovered evidence

`recovered-evidence.tar.gz` retains the 3,085 replayed and independently audited
bases, the separate `refinements.jsonl` records, producer manifests, source
receipts, native output and Magma transcripts. Extract it into a new directory;
its internal paths are relative to that directory. `recovery-summary.json`
records its digest, counts and local measurements. These records add 3,085
bases to the retained 1,041,991, giving 1,045,076 certified bases and 352,522
open bases. They credit zero new splits.

[The recovery commands](../../../production/README.md#retained-f_2-basis-recovery)
and [record contract](../../refinements-v1.md) describe verification and resume.
This is a retained-run refinement overlay; composing a full release remains
separate work.
