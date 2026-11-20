import copy
import csv
import importlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest

ROOT = Path.cwd()
sys.path.insert(0, str(ROOT))


def command(*args):
    return subprocess.run([sys.executable, 'production/refine.py', *map(str, args)],
                          text=True, capture_output=True)


class Sha2EvidenceImport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        scratch = ROOT / 'artifacts/test-tmp'
        scratch.mkdir(parents=True, exist_ok=True)
        temporary = tempfile.TemporaryDirectory(prefix='sha2-import-', dir=scratch)
        cls.addClassCleanup(temporary.cleanup)
        cls.output = Path(temporary.name)
        result = command('import-sha2', '--output', cls.output)
        if result.returncode:
            raise AssertionError('The public Sha[2] importer must accept the saved evidence:\n'
                                 + result.stdout + result.stderr)
        cls.summary = json.loads(result.stdout)
        cls.records = {int(r['key']): r for r in map(json.loads,
            (cls.output / 'refinements.jsonl').read_text().splitlines())}
        with tarfile.open(cls.output / 'evidence.tar.gz') as archive:
            cls.raw = {r['key']: r for r in json.load(archive.extractfile('inputs/unresolved_raw.json'))}
            cls.outcomes = list(csv.DictReader(io.StringIO(
                archive.extractfile('outcomes.csv').read().decode())))
        cls.api = importlib.import_module('pipeline.sha2')
        cls.compose = staticmethod(importlib.import_module('pipeline.records').compose_sha2)

    def test_all_saved_rows_round_trip(self):
        self.assertEqual(self.summary['rows'], 6008)
        self.assertEqual(self.summary['settlements'], 3418)
        self.assertEqual(self.summary['shortened_unresolved'], 43)
        self.assertEqual(self.summary['splits_proved'], 1395008)
        self.assertEqual(self.summary['splits_open'], 2590)
        self.assertEqual(self.summary['new_arithmetic'], 0)
        self.assertEqual(set(self.records), set(self.raw))
        for out in self.outcomes:
            record = self.records[int(out['key'])]
            claim = record['claims']['sha2']
            self.assertEqual([int(p[0]) for p in record['claims']['split']['candidates']],
                             json.loads(out['new_candidate_list']))
            if out['new_sha2_dimension']:
                self.assertEqual((claim['status'], claim['value']),
                                 ('exact', int(out['new_sha2_dimension'])))
            elif out['ordinary'] == 'False':
                self.assertEqual(claim['status'], 'not_applicable')
            else:
                self.assertNotEqual(claim['status'], 'exact')
            source = self.raw[int(out['key'])]
            composed = self.compose(source, record)
            self.assertEqual(composed['generators'], source['generators'])
            for sha, reg in composed['split']['candidates']:
                from fractions import Fraction
                self.assertEqual(sha * Fraction(reg), Fraction(source['derived']['sha_reg']))
        result = command('verify-sha2', '--output', self.output)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout), self.summary)
        result = command('import-sha2', '--output', self.output)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout), self.summary)

    def test_independent_formula_controls(self):
        controls = [(536889474, [4], 2), (536938506, [4], 2),
                    (543690369, [1, 4], None), (808473738, [1, 4], None),
                    (1357551361, [4], 2), (842431119, [4, 16], 2)]
        for key, orders, dimension in controls:
            result = self.api.refine_sha2(self.raw[key])
            self.assertEqual([int(x[0]) for x in result['split']['candidates']], orders)
            self.assertEqual(result['sha2']['value'], dimension)
            self.assertEqual(result['sha2']['status'], 'bound' if dimension is None else 'exact')
        source = next(r for r in self.raw.values() if not r['codes'][0])
        result = self.api.refine_sha2(source)
        self.assertEqual(result['sha2']['status'], 'not_applicable')
        self.assertEqual(result['split']['candidates'],
                         [[str(n), reg] for n, reg in source['split']['candidates']])
        with self.assertRaisesRegex(ValueError, 'empty'):
            self.api.refine_sha2(self.raw[536889474], defect=2)

    def test_all_saved_rows_against_magma(self):
        # The private routine has no public Magma caller; compare it with the
        # Python entry point on every saved row, with and without saved defects.
        cases = []
        for key, source in sorted(self.raw.items()):
            claim = self.records[key]['claims']['sha2']
            defects = [None]
            if claim.get('defect_status') == 'exact':
                defects += claim['defect_candidates']
            descent = source['split']['descent']
            for defect in defects:
                expected = self.api.refine_sha2(source, defect=defect)
                dimension = expected['sha2']
                args = [str(bool(source['codes'][0])).lower(),
                        str(descent['status'] == 'ok').lower()]
                args += [str(descent.get(name, 0)) for name in ('f', 'v')]
                args += [str(source['derived']['rank']), str(descent.get('sha2_upper', 0)),
                         str([n for n, _ in source['split']['candidates']])]
                want = [json.dumps(dimension['status']),
                        str(dimension.get('candidates', [])),
                        str(dimension.get('defect_candidates', [])),
                        str([int(n) for n, _ in expected['split']['candidates']])]
                call = 'sha2_refinement(' + ', '.join(args)
                call += ' : Defect := ' + str(-1 if defect is None else defect) + ')'
                cases.append('    error if ' + call + ' ne <' + ', '.join(want)
                             + '>, "Sha[2] mismatch for key ' + str(key)
                             + ', defect ' + str(defect) + '";')
        script = self.output / 'compare-sha2.m'
        script.write_text('\n'.join([
            'SetQuitOnError(true);', 'AttachSpec("spec");',
            'import "magma/arithmetic/sha2.m": sha2_refinement;',
            'procedure test_saved_sha2_formula()', *cases, 'end procedure;',
            'test_saved_sha2_formula();', 'quit;', '']))
        result = subprocess.run([os.environ.get('MAGMA', 'magma'), '-bn', str(script)],
                                text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_composition_rejects_conflicts_and_keeps_exact_claims(self):
        key = 536889474
        source, record = self.raw[key], self.records[key]
        with self.assertRaisesRegex(ValueError, 'conflict|empty'):
            self.compose(source, record, prior_sha2={'status': 'exact', 'value': 0})
        bad = copy.deepcopy(record)
        bad['q'] = 3
        with self.assertRaisesRegex(ValueError, 'field|identity'):
            self.compose(source, bad)
        bad = copy.deepcopy(record)
        bad['key'] = '536938506'
        with self.assertRaisesRegex(ValueError, 'key|identity'):
            self.compose(source, bad)
        proved = copy.deepcopy(source)
        proved['split'].update(proved=True, sha=1, regulator=source['derived']['sha_reg'],
                               candidates=[source['split']['candidates'][0]])
        with self.assertRaisesRegex(ValueError, 'conflict|empty'):
            self.compose(proved, record)
        skipped = next(r for r in self.records.values()
                       if r['claims']['sha2']['status'] == 'not_applicable')
        old = {'status': 'exact', 'value': 0, 'method': 'brauer_invariant_factors'}
        composed = self.compose(self.raw[int(skipped['key'])], skipped, prior_sha2=old)
        self.assertEqual(composed['sha2'], old)

    def test_repeated_composition_intersects_all_dimension_claims(self):
        source, record = self.raw[536889474], self.records[536889474]
        composed = self.compose(source, record)
        self.assertEqual(composed['sha2'], record['claims']['sha2'])
        self.assertEqual(self.compose(composed, record), composed)

        source, record = self.raw[5515777], self.records[5515777]
        exact = {'status': 'exact', 'value': 2}
        composed = self.compose(source, record, prior_sha2=exact)
        self.assertEqual(composed['split']['candidates'], [[4, '128'], [16, '32']])
        for prior in [{'status': 'exact', 'value': 4},
                      {'status': 'bound', 'candidates': [4, 6]}]:
            with self.subTest(prior=prior):
                with self.assertRaisesRegex(ValueError, 'conflicting Sha\\[2\\] dimension'):
                    self.compose(composed, record, prior_sha2=prior)
        for prior in [{'status': 'not_computed', 'value': None},
                      {'status': 'bound', 'candidates': [0, 2, 4]}]:
            with self.subTest(prior=prior):
                repeated = self.compose(composed, record, prior_sha2=prior)
                self.assertEqual(repeated['sha2']['status'], 'exact')
                self.assertEqual(repeated['sha2']['value'], 2)
                self.assertEqual(repeated['split']['candidates'], [[4, '128'], [16, '32']])
                self.assertFalse(repeated['split']['proved'])
        bounded = self.compose(source, record,
                               prior_sha2={'status': 'bound', 'candidates': [0, 2]})
        incoming = copy.deepcopy(record)
        incoming['claims']['sha2'] = {'status': 'bound', 'candidates': [0, 4]}
        intersected = self.compose(bounded, incoming,
                                   prior_sha2={'status': 'bound', 'candidates': [0, 2, 4]})
        self.assertEqual((intersected['sha2']['status'], intersected['sha2']['value']),
                         ('exact', 0))
        self.assertEqual(intersected['split']['candidates'], [[1, '512']])

        source, record = self.raw[540934786], self.records[540934786]
        priors = [dict(status='bound', value=None, candidates=values,
                       method='producer-' + name)
                  for name, values in [('A', [0, 2]), ('B', [0, 4]), ('C', [0, 2])]]
        original_priors = copy.deepcopy(priors)
        for status in ['not_computed', 'not_applicable']:
            with self.subTest(status=status):
                incoming = copy.deepcopy(record)
                incoming['claims']['sha2'] = dict(status=status, value=None)
                composed = self.compose(source, incoming, prior_sha2=priors[0])
                self.assertEqual(composed['sha2'], priors[0])
                for count, prior in enumerate(priors[1:], start=2):
                    previous = copy.deepcopy(composed)
                    repeated = self.compose(composed, incoming, prior_sha2=prior)
                    self.assertEqual(repeated['sha2'], dict(
                        status='exact', value=0, candidates=[0], method='intersection',
                        prior_evidence=priors[:count]))
                    self.assertEqual(repeated['split']['candidates'], [[1, '234'], [9, '26']])
                    self.assertFalse(repeated['split']['proved'])
                    self.assertEqual(composed, previous)
                    self.assertEqual(priors, original_priors)
                    composed = repeated

        composed, record = self.raw[547761537], self.records[547761537]
        prior = dict(status='bound', candidates=[0, 2], method='producer-A')
        for iteration in range(2):
            composed = self.compose(composed, record, prior_sha2=prior)
            with self.subTest(iteration=iteration, evidence='prior_evidence'):
                self.assertEqual(len(composed['sha2']['prior_evidence']), 2)
            with self.subTest(iteration=iteration, evidence='refinement_evidence'):
                self.assertEqual(len(composed['split']['refinement_evidence']), 1)

    def test_changed_output_or_evidence_is_rejected(self):
        path = self.output / 'refinements.jsonl'
        saved = path.read_bytes()
        try:
            path.write_bytes(saved + saved.splitlines(keepends=True)[0])
            result = command('verify-sha2', '--output', self.output)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('changed', result.stderr)
        finally:
            path.write_bytes(saved)
        path = self.output / 'evidence.tar.gz'
        saved = path.read_bytes()
        try:
            path.write_bytes(saved + b'changed')
            result = command('verify-sha2', '--output', self.output)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('digest', result.stderr)
        finally:
            path.write_bytes(saved)


if __name__ == '__main__':
    unittest.main()
