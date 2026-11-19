# Retained F_2 recovery control

Thanks for helping check the recovery path. `case.json` retains key 4999168,
whose bound-4 generator record is rank-short and whose bound-5 basis has
Gram matrix `[14]`. Its proved split is Sha order 16 and regulator 14.

The source row is ordinal 366 of `p0002-s00000.jsonl.gz` in the saved F_2
run. The compressed shard digest is
`7c9b79eabe731785f339efd479500c29861468bd8f3a773b93d368162587441c`;
its retained producer identity is
`332efa40ac87b5539a9c08b7f777f4d093203559d8fc30d99cb7de7ea756937b`.
The request is the first line of `data/mw/f2-recovery/replay-candidates.jsonl`.
The native result was replayed with `stages-native --bundle` built from
calibration revision `45a0aa6f2fa7a006095da0201f77ced69cc70dfe`.

The fixture was extracted by matching the key in that compressed shard and
retaining the request, complete source row and replay output. The test builds
a one-row source receipt and manifest from these bytes. Magma independently
checks the equation, height matrix, torsion, local bounds and division before
import. The fixture alone establishes no arithmetic claim.

The command tests also cover missing and capped work, wrong and duplicate keys,
invalid torsion evidence, and incompatible kernel identities. They check failed
GNU parallel job recovery and preserve an exact basis after an incomplete retry. This control does not cover the full
recovery population or implement the later refinement stages.
