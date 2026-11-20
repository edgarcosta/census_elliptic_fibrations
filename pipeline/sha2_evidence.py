import collections
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile

from .sha2 import refine_sha2, require, source_identity

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / 'data/sha2/f2-ordinary-refinements'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n').encode()


def matrix_rank(values, size):
    basis = {}
    for i in range(size):
        x = sum(values[size*i+j] << j for j in range(size))
        while x:
            pivot = x.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = x
                break
            x ^= basis[pivot]
    return len(basis)


def pairing(values, size, rank):
    require(len(values) == size**2 and all(x in (0, 1) for x in values), 'pairing entries')
    require(all(values[size*i+i] == 0 for i in range(size)) and all(
        values[size*i+j] == values[size*j+i] for i in range(size) for j in range(size)),
        'non-alternating pairing')
    require(matrix_rank(values, size) == rank and rank % 2 == 0, 'pairing rank')


class Sha2Evidence:
    def __init__(self, archive, manifest):
        self.archive = archive
        self.manifest = json.loads(manifest)
        require(digest(archive) == self.manifest['archive_sha256'], 'evidence archive digest')
        with tarfile.open(fileobj=io.BytesIO(archive), mode='r:gz') as tar:
            entries = tar.getmembers()
            require(all(e.isfile() for e in entries) and len({e.name for e in entries}) == len(entries),
                    'duplicate or unsupported evidence member')
            self.files = {e.name: tar.extractfile(e).read() for e in entries}
        require(self.files.keys() == self.manifest['files'].keys(), 'evidence inventory')
        for name, expected in self.manifest['files'].items():
            data = self.files[name]
            require(digest(data) == expected['sha256'] and len(data) == expected['bytes'],
                    f'evidence digest: {name}')
        for entry in self.json('manifest.json'):
            require(digest(self.files[entry['path']]) == entry['sha256'], 'original artifact digest')
        for entry in self.json('sources.json'):
            require(digest(self.files['inputs/' + Path(entry['path']).name]) == entry['sha256'],
                    'original source digest')
        require(self.json('conflicts.json') == [], 'saved pairing conflicts')
        self.raw = self.json('inputs/unresolved_raw.json')
        self.by_key = {r['key']: r for r in self.raw}
        require(len(self.by_key) == len(self.raw) == self.manifest['outcomes'], 'duplicate source key or row count')
        release = {r['key']: r for r in self.json('inputs/unresolved.json')}
        require(release.keys() == self.by_key.keys(), 'original release keys')
        for r in self.raw:
            source_identity(r)
            old = release[r['key']]
            require(old['analytic_rank'] == r['derived']['rank']
                    and old['admissible_sha'] == [n for n, _ in r['split']['candidates']],
                    'original release association')
        self.outcomes = self.csv('outcomes.csv')
        require({int(r['key']) for r in self.outcomes} == self.by_key.keys()
                and len(self.outcomes) == len(self.raw), 'outcome keys')
        for name, expected in [('settled.csv', [r for r in self.outcomes
                if len(json.loads(r['new_candidate_list'])) == 1]),
                ('refined.csv', [r for r in self.outcomes
                if r['new_candidate_list'] != r['old_candidate_list']])]:
            actual = self.csv(name)
            require(len(actual) == len(expected) and all(
                all(r[k] == e[k] for k in r) for r, e in zip(actual, expected)), 'saved CSV disagreement')
        self.results = self.audit_pairings()
        self.producer = dict(method='saved-ordinary-f2-sha2-v1', evidence_sha256=digest(archive),
            evidence_manifest_sha256=digest(manifest), implementation={
                name: digest((ROOT / name).read_bytes()) for name in (
                    'pipeline/sha2.py', 'pipeline/sha2_evidence.py', 'pipeline/records.py', 'production/refine.py')})
        self.producer_id = digest(encoded(self.producer))

    def json(self, name):
        return json.loads(self.files[name])

    def csv(self, name):
        return list(csv.DictReader(io.StringIO(self.files[name].decode())))

    def audit_pairings(self):
        results = collections.defaultdict(list)
        for name in sorted(self.files):
            if not name.startswith('checkpoints/'):
                continue
            checkpoint = self.json(name)
            if checkpoint['status'] != 'ok':
                continue
            key, v, f, c = (checkpoint[x] for x in ('key', 'v', 'f', 'c_E'))
            text = self.files[checkpoint['log']].decode()
            matches = re.findall(r'^RESULT\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|([^\n]*)', text, re.M)
            require(len(matches) == 1 and f'COMPLETE|{key}' in text and 'error' not in text.lower(),
                    'incomplete pairing transcript')
            k0, v0, f0, c0, matrix = matches[0]
            values = json.loads(matrix)
            require(tuple(map(int, (k0, v0, f0, c0))) == (key, v, f, c)
                    and values == checkpoint['matrix'], 'pairing transcript association')
            pairing(values, v, c)
            local = re.findall(r'^PLACE key=(\d+) at=.* MATRIX=(\[[^\n]*\])$', text, re.M)
            total = [0] * (v*v)
            for local_key, entries in local:
                entries = json.loads(entries)
                require(int(local_key) == key and len(entries) == len(total)
                        and all(x in (0, 1) for x in entries), 'local pairing association')
                total = [x ^ y for x, y in zip(total, entries)]
            if checkpoint['column']:
                j = checkpoint['column'] - 1
                require(v == 2 and j in (0, 1), 'one-column pairing domain')
                total[2*j+1-j] = total[2*(1-j)+j]
            require(total == values, 'local pairing sum')
            if checkpoint['label'] == 'residue':
                row = self.by_key[key]
                descent = row['split']['descent']
                require((descent['f'], descent['v']) == (f, v), 'pairing Selmer dimensions')
                claims = refine_sha2(row, defect=c)
                require(claims['sha2']['value'] == checkpoint['sha2_dim']
                        and [int(n) for n, _ in claims['split']['candidates']] == checkpoint['after'],
                        'pairing formula disagreement')
                results[key].append((name, checkpoint))
        for attempts in results.values():
            require(len({(r['c_E'], r['sha2_dim'], tuple(r['after'])) for _, r in attempts}) == 1,
                    'conflicting exact pairing claims')
        require(set(self.json('selection.json')['targets']) == results.keys(), 'unfinished residue targets')
        sample_keys = {int(r['key']) for r in self.json('residue-check/sample.json')}
        for name in ('residue-check/results-p200-s1.jsonl', 'residue-check/results-p320-s2.jsonl'):
            sample = [json.loads(line) for line in self.files[name].splitlines()]
            require(len(sample) == len(sample_keys) and {r['key'] for r in sample} == sample_keys,
                    'independent sample coverage')
            for r in sample:
                expected = results[r['key']][0][1]
                require(r['status'] == 'ok' and (r['f'], r['v'], r['c_E']) ==
                        (expected['f'], expected['v'], expected['c_E']), 'independent residue disagreement')
                pairing([x for row in r['matrix'] for x in row], r['v'], r['c_E'])
        return results

    def records(self):
        outcomes = {int(r['key']): r for r in self.outcomes}
        records = []
        for ordinal, row in enumerate(self.raw):
            key = row['key']
            attempts = self.results.get(key, [])
            defect = attempts[0][1]['c_E'] if attempts else None
            claims = refine_sha2(row, defect)
            out = outcomes[key]
            descent = row['split']['descent']
            require(out['ordinary'] == str(bool(row['codes'][0]))
                    and int(out['r']) == row['derived']['rank']
                    and json.loads(out['old_candidate_list']) == [n for n, _ in row['split']['candidates']]
                    and [int(n) for n, _ in claims['split']['candidates']] == json.loads(out['new_candidate_list']),
                    'outcome source or survivor disagreement')
            require((out['method'] == 'residue') == bool(attempts), 'outcome pairing association')
            if descent['status'] == 'ok':
                require((int(out['f']), int(out['v'])) == (descent['f'], descent['v']), 'outcome Selmer association')
            else:
                require(not out['f'] and not out['v'], 'outcome without saved Selmer dimensions')
            if out['new_sha2_dimension']:
                require(claims['sha2']['status'] == 'exact'
                        and claims['sha2']['value'] == int(out['new_sha2_dimension']), 'outcome dimension disagreement')
            else:
                require(claims['sha2']['status'] != 'exact', 'unproved exact dimension')
            if attempts:
                require(int(out['c_E']) == defect, 'outcome defect disagreement')
            refs = ['outcomes.csv', 'review.md', 'residue-check.md', 'crosscheck.md']
            refs += [p for name, r in attempts for p in (name, r['log'])]
            evidence = [dict(archive='evidence.tar.gz', member=name,
                             sha256=self.manifest['files'][name]['sha256']) for name in refs]
            records.append(dict(schema='census-refinement', schema_version=1, q=2, key=str(key),
                key_encoding_version=1, stage='bounded', input_ordinal=ordinal,
                producer_id=self.producer_id, stage1_receipt_sha256=None,
                input_identity=dict(kind='retained-unresolved-snapshot', surface=source_identity(row),
                    snapshot='inputs/unresolved_raw.json',
                    snapshot_sha256=self.manifest['files']['inputs/unresolved_raw.json']['sha256'],
                    row_sha256=digest(encoded(row))),
                status='complete' if claims['split']['status'] == 'exact' else 'unresolved',
                method='saved-ordinary-f2-sha2-v1', reason=out['reason'], claims=claims, evidence=evidence))
        return records

    def summary(self, records):
        settled = sum(r['claims']['split']['status'] == 'exact' for r in records)
        shortened = sum(len(r['claims']['split']['candidates']) > 1 and
            len(r['claims']['split']['candidates']) < len(self.by_key[int(r['key'])]['split']['candidates'])
            for r in records)
        m = self.manifest
        require(settled == m['settlements'] and shortened == m['shortened_unresolved']
                and len(records) - settled == m['splits_open']
                and m['previous_splits_proved'] + settled == m['splits_proved']
                and m['splits_proved'] + m['splits_open'] == m['population'], 'saved population accounting')
        return dict(rows=len(records), settlements=settled, shortened_unresolved=shortened,
                    splits_proved=m['splits_proved'], splits_open=m['splits_open'],
                    new_arithmetic=0, producer_id=self.producer_id,
                    sha2_statuses=dict(collections.Counter(r['claims']['sha2']['status'] for r in records)))
