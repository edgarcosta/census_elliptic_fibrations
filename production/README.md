# Native stage 1 runs

The driver writes representative keys, exact native CSV v1 records and prefix
metadata. Run from a pinned Linux checkout with Python 3.11 or later, Magma,
GNU parallel, zstd and the native build dependencies.

```sh
make -j4 census-native &&
python3 production/stage1.py prepare --q 3 --output results/f3 &&
seq 1 1860 |
parallel -j512 --halt soon,fail=1 --resume-failed --joblog results/f3/jobs.tsv --results results/f3/logs 'python3 production/stage1.py run --q 3 --output results/f3 --prefixes {}' &&
python3 production/stage1.py verify --q 3 --output results/f3 --expected-count 1307723786
```

`prepare` checks counts, prefix forms and producer identity before writing
the immutable manifest.
`--magma` selects the executable; its default is `$MAGMA` or `magma` from PATH.

Each `run` handles one prepared prefix. To resume, repeat preparation, the
same GNU parallel invocation and verification. Keep the input order, command,
output directory and job log unchanged. If the job log is lost, submit every
prefix; the driver validates completed receipts before skipping them.

Both fields use these commands: F_2 has 88 prefixes and F_3 has 1,860.
Use separate output directories for each field. For a bounded check,
pass `--prefixes 1,2,3` to `prepare` and `verify`, then run each selected prefix.
This reports `partial: true`; the prepared selection is fixed.

## Output and verification

Each directory under `prefixes/` contains keys in the
[binary format](../docs/keys.md#binary-bytes), CSV and stats metadata,
compressed rows and `receipt.json`. Receipts pin selections, prefix forms,
encoding and schema versions, counts, file hashes and producer identity.
Keep the checkout pinned while a run is in progress; any commit changes that identity.
Jobs lock their prefix under `locks/`. Incomplete output without a receipt
moves under `failed/` before recomputation, retaining logs. Corrupt or
incompatible completed receipts stop the command for diagnosis.

The serial `verify` command recomputes package counts and prefix forms,
validates every required receipt and writes `verification.json`.
Stream validation checks the zstd checksum, CSV header, smooth minimal labels
and exact ordered key correspondence. Empty prefixes also need receipts.
Verification does not recompute arithmetic rows or certify Mordell-Weil bases.
A failed command can leave an earlier verification file; use its exit status
to judge the result.
