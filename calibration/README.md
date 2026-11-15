# Calibrate the census arithmetic

Thank you for measuring the pipeline on the target machine. This directory
builds and runs the arithmetic and table writer on frozen F_2 and F_3 fixtures.

The source comes from the local `epyc-bundle` calibration, with its existing
arithmetic checks. The included Magma census package is revision `0bdf62a`.
The `reference` files use Git LFS; fetch them before building. Build the native
executables on each target machine.

## Build

Requirements: Python 3.10+, PyArrow 26.0.1, C++17, Boost JSON headers (1.75+),
NTL, FLINT, GMP, util-linux (`taskset` and `lscpu`), Git LFS and licensed Magma.
The build checks the libraries and reports missing dependencies. Its temporary
files and executables stay in `build` and `bin`.

From this directory, use `./build.sh --jobs 4`. The default builds five native
executables and certifies all 24 point-count tables through degree 12.
`--scalar` selects the scalar build and must also be passed to `census.py`.
For a local Python dependency, install PyArrow with
`python3 -m pip install --target build/python pyarrow==26.0.1`.

Magma calls use `magma -bn`, with `MAGMA_STARTUP_FILE` unset. Each script starts
with `SetQuitOnError(true);` and ends with `quit;`. The target version is
2.29-10; the local trial uses 2.30-1.

## Measure

Run `python3 census.py --output results/epyc` from this directory. Use a fresh
output directory for each run. Runs are sequential across worker counts and
fields; each child process is pinned to one allowed CPU. By default, distinct
physical cores precede sibling threads. `--cpus` supplies an explicit ordered
list. Numerical libraries and PyArrow use one thread per process.

At 1, 16, 64 and 256 workers, every worker receives all 2,000 F_2 rows or all
5,000 F_3 rows in a separate shard. Repeated keys live in separate files.
This measures weak scaling, with 512,000 F_2 and 1,280,000 F_3 row evaluations
at 256 workers. Each shard starts fresh executables, so startup CPU is included
at the same batch size across worker counts.

One further run per field partitions the original fixture into exactly 256
nonempty shards. This fixed-total run records the tail and startup overhead
of a small sample. It includes only 2,000 or 5,000 row evaluations.

For a small local trial, use
`python3 census.py --workers 1,2,4 --tail-workers 4 --limit 13 --output results/trial`.
This runs the same arithmetic, comparison, and reporting paths.

## Read the measurements

`results/epyc/summary.json` is written only after every run succeeds.
`progress.json` records completed runs. Each run has its own `summary.json`,
Parquet fragments and JSON receipts under `f2` or `f3`. `host.json` records CPU
placement, topology, compiler, executable digests and source identity.

Reports give total CPU seconds, CPU seconds per surface, wall seconds,
surfaces per second, shard sizes, compressed bytes and maximum child RSS.
The RSS value is a process high-water mark, not the aggregate resident memory
of concurrent workers. Shard CPU and wall times exclude fixture loading and
comparison. Run wall time includes worker startup, fixture checks and
collection; serial final readback is outside the timed interval.

The Linux process pool uses fork to start every worker before dispatching
shards, including short fixed-total shards. Each run requires a distinct
worker PID for every shard.

`stage_cpu_per_surface` separates decode, L-function/local data (`lfun`),
bounded generators, derived invariants (`derive`, `extras`), admissible lists,
d2, split filters/search, descent, and writer/dispatch overhead. Absent stages
were not reached. The L-function and generator executables each bundle some
local setup; these timings do not split their internal algorithms. Exact table
certification and Magma key-decoding checks are recorded separately.

Every arithmetic record is compared with its frozen fixture. Parquet is read
back exactly, and F_3 lean rows must reconstruct all full fields. The fixture
provenance and independent coverage are recorded in `reference`. Regression
agreement does not independently prove every certificate.

This measures the existing smooth, everywhere-minimal arithmetic pipeline.
It excludes enumeration, final packing, unbounded generator searches and
lattice assembly. It does not calibrate the other fibration classes.

After building, run `python3 -m unittest discover -s . -p 'test_census.py'`
from this directory for the timing and worker-count regression checks.
