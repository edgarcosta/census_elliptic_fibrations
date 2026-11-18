# Native stage 1 runs

Thanks for helping run and check the census. This driver retains exact native
CSV v1 records, representative keys and prefix metadata for later computations.
Its output is intermediate data; it does not certify Mordell-Weil bases or
produce a release. The F_3 coefficient-orbit mass check is certified; the
independent check of smooth and minimal labels, and angle ranks, remain open.

Run from a pinned checkout on Linux, with Python 3.11 or later, Magma, GNU
parallel, and the native build dependencies in the root README. Python runs
the `zstd` command-line tool to read compressed rows; no Python packages are
required. Build with `make -j4 census-native`. When using local zstd development
files, match the headers to the installed library and test compilation and
linking first.

Magma runs with `-bn` and one thread; numerical-library thread counts are also
one. The driver removes `MAGMA_STARTUP_FILE` from child environments and
supplies its own scratch directory. `--magma` selects the executable; its
default is `$MAGMA` or `magma` from PATH. Receipts contain no licence credentials.

## Prepare, run, verify

The serial pilot compares fixed native rows with Magma and independent
`LFunction` values, then checks enumeration of a largest prefix. It retains
1,594,323 keys for the largest F_3 prefix. Its arithmetic controls include the
native Sha lists `[1,4,9]` and `[1]` for keys 2729765910770 and 66249649205.
Each pilot attempt uses new directories; a failed attempt removes the earlier
`pilot.json`. It records timings, not resource measurements.

`prepare` checks package counts, prefix forms and producer identity, then
writes an immutable manifest. Full preparation requires a successful pilot
with matching identity and at least 150,000,000,000 free bytes, less the bytes
already retained under `prefixes/`. Scratch and failed attempts do not count
towards this reserve. `--pilot-data PATH` selects the pilot file; the default
is `stage1-pilot/pilot.json` beside the output directory.

`run --prefixes N` handles exactly one prepared prefix. It checks identity,
validates a completed receipt before skipping it, and otherwise computes and
publishes the prefix. GNU parallel schedules the independent jobs:

```sh
python3 production/stage1.py pilot --q 3 --output results/f3-pilot &&
python3 production/stage1.py prepare --q 3 --output results/f3 --pilot-data results/f3-pilot/pilot.json &&
seq 1 1860 |
parallel -j512 --halt soon,fail=1 --resume-failed --joblog results/f3/jobs.tsv --results results/f3/logs 'python3 production/stage1.py run --q 3 --output results/f3 --prefixes {}' &&
python3 production/stage1.py verify --q 3 --output results/f3 --expected-count 1307723786
```

`-j512` permits up to 512 jobs; choose the job count from host measurements.
The driver does not pin CPUs or check available RAM. `--halt soon,fail=1`
stops new launches after one failure and waits for running jobs.
`--resume-failed` retries failed jobs and starts unfinished jobs, using
`--joblog` to identify successful sequence numbers. `--results` saves each
job's stdout, stderr and sequence number under the specified directory.

To resume, repeat `prepare`, the same parallel invocation, and `verify`.
Keep the input order, job command, output directory and job log unchanged.
Preparation must run even when the job log records every job as successful:
it checks identity before GNU parallel can skip those jobs. Verification
then checks that every required file and receipt still exists and is valid.
If the job log is lost, submit all prefixes; each job validates its receipt
before skipping completed work. Job output files can be overwritten on retry;
producer logs retained under the output directory preserve earlier attempts.

Both fields use the same commands. F_2 has 88 prefixes and 1,397,598 rows;
F_3 has 1,860 prefixes and 1,307,723,786 rows. Empty prefixes also need receipts.
For a bounded check, use `--prefixes 1,2,3` on `prepare` and `verify`, then call
`run` separately for each selected prefix. Such verification reports
`partial: true`, even when the selected count matches `--expected-count`.
The prepared selection is fixed; use another output for a different selection.
Keep pilot and run outputs for different fields separate.

## Retained data and recovery

Each prefix has a directory under `prefixes/` containing eight-byte
little-endian keys, CSV and stats metadata, compressed rows and `receipt.json`.
Receipts pin selections, prefix forms, encoding/schema versions, counts and
file SHA-256 digests. They also pin the Git commit, arithmetic and orchestration
source digest, native executable digest, and Magma launcher digest and version.
The source digest includes `magma/*.m`, both specs and the production scripts.

Stream validation checks the zstd content checksum, CSV header, successful
smooth/minimal classification and exact ordered key correspondence. It uses
bounded buffers. Decimal native input uses a temporary disk file of up to 25,509,168 bytes
per largest F_3 prefix. Allow about 13.1 GB for 512 such jobs. Allow this space as well as the preparation reserve.
Counting and listing share package mathematics; their agreement does not
prove enumeration independently.

Each job locks its prefix under `locks/`; different prefixes run independently.
Jobs leave the manifest and verification file alone. Run preparation before jobs
and verification after jobs. The final serial verifier recomputes package
counts and prefix forms, validates all retained receipts, and alone writes
`verification.json`. It does not recompute every arithmetic row.
Local validation measured about 301,000 rows per second, implying about
72 minutes for all F_3 rows; measure this on the EPYC.

Send stop requests to GNU parallel. On the checked version 20240222, SIGHUP
stops new launches and waits for running jobs; SIGTERM or Ctrl-C terminates
running jobs and their descendants. A SIGHUP drain can exit 0 with unstarted
prefixes, so the recipe then attempts verification. That rejects missing
receipts before revalidating retained prefixes. Repeat preparation, the parallel
invocation and verification to finish the run. The driver uses default signal
handling and starts no separate child sessions. Enumeration and arithmetic
calls have no automatic timeout.

Incomplete scratch or published directories without receipts move under
`failed/` before recomputation, retaining logs. Files and parent directories
are flushed before atomic publication. Corrupt or incompatible completed
receipts stop the command and remain available for diagnosis.
Keep a running checkout pinned: any commit, including documentation changes,
changes receipt identity. Retain keys, metadata, rows, logs, pilot data and
receipts. A failed command can leave an earlier verification file;
use the command's exit status to judge its result.

Run `make check/run_stage1.m` for the command checks and the native suite for
arithmetic controls. The command checks cover GNU parallel failure, cancellation
and resume when it is installed. They do not simulate low free disk space or
measure EPYC resource peaks. Monitor free space, including failed attempts,
and use ordinary host tools to measure memory and choose the job count.
