#!/usr/bin/env python3
"""Record process-tree CPU, wall time and sampled aggregate RSS on Linux."""

import json
from pathlib import Path
import resource
import subprocess
import sys
import time

out = Path(sys.argv[1])
start = time.monotonic()
p = subprocess.Popen(sys.argv[2:])
peak = samples = 0
while p.poll() is None:
    parents = {}
    rss = {}
    for path in Path('/proc').glob('[0-9]*/stat'):
        try:
            fields = path.read_text().rsplit(')', 1)[1].split()
            pid = int(path.parent.name)
            parents[pid] = int(fields[1])
            rss[pid] = int(fields[21]) * resource.getpagesize()
        except (FileNotFoundError, ProcessLookupError, PermissionError, IndexError):
            continue
    descendants = {p.pid}
    while True:
        new = {pid for pid, parent in parents.items() if parent in descendants}
        if new <= descendants:
            break
        descendants |= new
    peak = max(peak, sum(rss.get(pid, 0) for pid in descendants))
    samples += 1
    time.sleep(0.05)
code = p.wait()
r = resource.getrusage(resource.RUSAGE_CHILDREN)
out.write_text(json.dumps(dict(command=sys.argv[2:], exit_status=code,
    wall_seconds=time.monotonic()-start, user_seconds=r.ru_utime,
    system_seconds=r.ru_stime, max_process_rss_kib=r.ru_maxrss,
    sampled_peak_aggregate_rss_bytes=peak, samples=samples,
    sampling_interval_seconds=0.05,
    memory_limit='Sampled sum of descendant RSS; shared pages counted per process, peaks between samples may be missed.'), indent=2)+'\n')
sys.exit(code)
