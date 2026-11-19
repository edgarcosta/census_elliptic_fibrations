# Accessible F_2 basis inputs

Thanks for helping finish the F_2 bases. This queue contains 131,678 open
bases that can benefit from searches through denominator bound 7. It excludes
all 3,085 independently recovered bases and 220,844 rank-one rows whose
remaining regulator candidates exceed 18.

The queue contains 84,520 rank-one, 46,592 rank-two and 566 rank-three rows.
The starting certified count is 1,045,076, with 352,522 bases still open.
The refined split count remains 2,590 open. Input export earns no new basis
or split credit.

`queue.jsonl.gz` is a deterministic gzip stream of canonical JSON lines,
sorted by key. Each line retains coefficient codes, L-polynomial, rank,
torsion factors, Sha*Reg, refined (Sha, regulator) candidates, sampling stratum
and source-shard association. It contains no new exact Sha[2] claim.
`manifest.json` pins its digest, all 371 source receipts and compressed shards,
the retained summary, the completed recovery archive and the refined list.
Its counts partition the entire 1,397,598-row source population.

The producer is `production/search_f2.py export`. It verifies source digests and counts. It intersects the saved candidates
with the 3,461 rows in the exact Sha[2] refinement CSV, then excludes the
3,085 recovered basis keys.
It rejects missing overlay keys, empty intersections and duplicate source keys.
Source hashes establish provenance; they do not independently certify arithmetic.

The reproducible pilot selects 4,096 keys using Python's random generator
with seed 20261119. It chooses one uniformly from each of nine joint strata
(rank, ordinary/supersingular generic fibre, torsion factors), then samples
uniformly without replacement from all other keys. This deliberately includes
rare strata. The remaining cohort is its exact complement. A smaller local
sample uses the same rule and supplies compatibility evidence only.

[The search commands](../../../production/README.md#accessible-f_2-basis-search)
use the unchanged pinned calibration kernel and independently audit every
promoted basis with Magma. Local evidence and measurements are retained in
`local-evidence.tar.gz` and `local-summary.json`. Target-machine pilot and
full-queue measurements remain pending.
