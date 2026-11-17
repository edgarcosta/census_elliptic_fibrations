# Native stage 1 runs

Thanks for helping run and check the census. This driver retains exact native
CSV v1 records, representative keys and prefix metadata for later computations.
Its output is intermediate data; it does not certify Mordell-Weil bases or
produce a release. The orbit/mass check and angle rank remain open.

Run from a pinned checkout on Linux, with Python 3.11 or later, Magma, `taskset`,
and the native build dependencies in the root README. Python runs the `zstd`
command-line tool to read compressed rows; no Python packages are required.
Build with `make -j4 census-native`. Check NTL, FLINT, GMP, zstd headers and the
linker path before detaching a remote job. When using a local zstd development
setup, match its headers to the installed library and test compilation and
linking before building. A failed setup directory should be moved aside before
retrying.

Each worker needs its own physical core. Select one allowed logical CPU from
each core. The driver rejects duplicate cores and restricts each worker to its
assigned CPU. Without `--cpus` it selects allowed physical cores. Magma runs
with `-bn` and one thread; numerical-library thread counts are also one. The
driver removes
`MAGMA_STARTUP_FILE` from child environments and supplies its own scratch
directory. `--magma` selects the executable; its default is `$MAGMA`
or `magma` from PATH. No licence credentials are written to receipts.

```sh
python3 production/stage1.py pilot --q 3 --threads 4 --output results/stage1-pilot
python3 production/stage1.py run --q 3 --threads 256 --output results/stage1-f3
python3 production/stage1.py verify --q 3 --output results/stage1-f3 --expected-count 1307723786
```

Add `--cpus 0,2,4,6` with the actual allowed CPU list for a four-worker pilot.
The full run reads `stage1-pilot/pilot.json` beside its output directory;
`--pilot-data PATH` selects another measurement file. A successful pilot at the
same producer identity is required for a full launch. It compares fixed native
rows with Magma and independent `LFunction` values, interrupts and resumes a
prefix, compares recovered files, and rejects corrupt and stale receipts.
Each pilot attempt uses new measurement and recovery directories, so the same
pilot output can be reused after rebuilding a producer.
The F_3 controls retain the native Sha lists `[1,4,9]` and `[1]` for keys
2729765910770 and 66249649205.

The pilot records enumeration time, native throughput, sampled aggregate RSS,
parent RSS, worker high-water RSS, and sampled peak disk bytes. RAM is sampled every
50 ms across the driver's process tree. Short-lived peaks can fall between
samples. Before launch the driver reserves 120% of measured worker peak times
worker count plus at least 256 MiB for the parent. A full launch also requires
150,000,000,000 free bytes on the first launch. On resume, bytes retained
under `prefixes/` count towards this reserve; scratch and failed attempts do not.
The pilot enumerates a prefix with the largest package count (1,594,323 keys
for F_3) and includes its sampled memory peak.
Other prefixes can have different peaks; monitor the host during production.
A bounded run without pilot measurements reserves 1 GiB per worker.
`--worker-memory-bytes` raises or sets this estimate for a separately measured
worker workload.

Both fields use the same commands. For a bounded check, select distinct prefix
indices with `--prefixes 1,2,3` on `run` and `verify`. Such verification reports
`partial: true` and the selected population, even if an expected count is given.
Without this option verification requires all 88 F_2 or 1,860 F_3 prefixes,
including empty ones, and totals 1,397,598 or 1,307,723,786 respectively.
Pilot and full commands reject a changed total before enumerating keys.
A prefix selection may be extended when resuming the same output directory.

One lock protects each output directory, including verification. Each prefix
has a directory under `prefixes/` containing eight-byte little-endian keys,
the writer's CSV and stats metadata, compressed rows and `receipt.json`.
Receipts pin selections, prefix forms, encoding/schema versions, counts, and
file SHA-256 digests. They also pin the Git commit, arithmetic and orchestration
source digest, native executable digest, and Magma launcher digest and version.
Stream validation checks the zstd content checksum,
CSV header, successful smooth/minimal classification and exact ordered key
correspondence. Validation uses bounded buffers; it keeps no whole-census key
set in memory. Counting and listing share package mathematics, so agreement
checks orchestration; it does not prove enumeration independently.

Send SIGTERM or SIGINT to stop. The driver stops its live children; completed
receipts stay intact. Repeat the same `run` command in the same pinned checkout
to resume. Only validated completed prefixes are skipped. Incomplete work in
owned `scratch/` directories is moved under `failed/` before recomputation,
retaining logs from every attempt. The driver prints each prefix start and
completion to stderr; a start without a completion identifies stalled work.
Startup revalidation also prints each prefix start and completion. Stop requests
are checked between these prefixes; the current prefix finishes validation first.
Enumeration and arithmetic calls have no automatic timeout; use SIGTERM and
resume if a prefix stalls. Files and their parent directories are flushed before
an atomic receipt publication. A corrupt or incompatible completed receipt stops
the command and is retained for diagnosis; it is never silently replaced.
Retain keys, metadata, rows, logs, measurements, and receipts.

Keep a running checkout unchanged. Develop later refinements in another
checkout so source identities remain compatible. Verification recomputes public
package counts and prefix forms before checking receipts and rows. It does not
recompute every row's arithmetic. Run `make -j2 check/run_stage1.m` for the driver
checks in `tests/compare/run_stage1.m`, and the native suite for arithmetic controls.

Validation and key feeding share one Python process. Local validation measured
about 301,000 rows per second, or about 72 minutes for all F_3 rows. Resume
and verify revalidate retained rows serially. The small-prefix pilot throughput
includes process startup costs and is not a full-run wall-time estimate.
Measure the EPYC run before projecting its completion time.

Limitations: Command tests do not simulate low free disk space. The test
directories have no free-space injection; the resumed disk reserve calculation
is checked by inspection. Monitor free space during production, especially
after failed attempts.
