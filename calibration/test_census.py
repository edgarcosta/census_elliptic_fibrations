import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parent
LAUNCH = r'''
import multiprocessing
from multiprocessing.context import ForkServerProcess
from pathlib import Path
import sys
import time

sys.path.insert(0, sys.argv[1])
import census

mode, output = sys.argv[2:]
if mode == 'timing':
    multiprocessing.set_start_method('fork', force=True)
    offset = [0.0]
    real_cpu, real_wall = census.cpu, time.monotonic
    census.cpu = lambda: real_cpu() + offset[0]
    time.monotonic = lambda: real_wall() + offset[0]
    real_open, real_compare = census.gzip.open, census.compare_records

    def fixture_open(*args, **kwargs):
        offset[0] += 10000
        return real_open(*args, **kwargs)

    def compare(*args, **kwargs):
        offset[0] += 10000
        return real_compare(*args, **kwargs)

    census.gzip.open = fixture_open
    census.compare_records = compare
    workers, limit = '1', '1'
else:
    multiprocessing.set_start_method('forkserver', force=True)
    real_start = ForkServerProcess.start

    def slow_start(self):
        real_start(self)
        time.sleep(2)

    ForkServerProcess.start = slow_start
    workers, limit = '4', '4'

sys.argv = ['census.py', '--workers', workers, '--tail-workers', workers,
            '--limit', limit, '--output', output]
census.main()
'''


class CalibrationTests(unittest.TestCase):
    def run_calibration(self, mode):
        results = ROOT / 'results'
        results.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='regression-', dir=results) as temp:
            output = Path(temp) / 'run'
            env = os.environ.copy()
            env.pop('MAGMA_STARTUP_FILE', None)
            env['PYTHONDONTWRITEBYTECODE'] = '1'
            proc = subprocess.run(
                [sys.executable, '-c', LAUNCH, str(ROOT), mode,
                 str(output.relative_to(ROOT))],
                cwd=ROOT, env=env, capture_output=True, text=True, timeout=120)
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            summary = json.loads((output / 'summary.json').read_text())
            receipts = [json.loads(p.read_text())
                        for p in output.glob('f*/*/shards/*.json')]
            return summary, receipts

    def test_fixture_checks_are_outside_shard_timings(self):
        summary, receipts = self.run_calibration('timing')
        self.assertEqual(len(summary['runs']), 4)
        for receipt in receipts:
            with self.subTest(q=receipt['q'], label=receipt['label']):
                self.assertTrue(receipt['regression_compared'])
                self.assertLess(receipt['cpu'], 1000)
                self.assertLess(receipt['wall'], 1000)

    @unittest.skipUnless(len(os.sched_getaffinity(0)) >= 4, 'requires four CPUs')
    def test_short_shards_use_every_requested_worker(self):
        summary, receipts = self.run_calibration('workers')
        self.assertEqual(len(summary['runs']), 4)
        for run in summary['runs']:
            with self.subTest(q=run['q'], mode=run['mode']):
                self.assertEqual(len(run['worker_cpus']), 4)
        for q in (2, 3):
            field = [m for m in receipts if m['q'] == q]
            self.assertEqual(len(field), 8)
            self.assertEqual(len({m['worker_pid'] for m in field}), 8)


if __name__ == '__main__':
    unittest.main()
