import copy
from fractions import Fraction

from .sha2 import require, source_identity


def compose_sha2(source, record, prior_sha2=None):
    identity = source_identity(source)
    require(record['schema'] == 'census-refinement' and record['schema_version'] == 1
            and record['key_encoding_version'] == 1 and record['stage'] == 'bounded',
            'incompatible refinement schema')
    require(record['q'] == identity['q'] and record['key'] == identity['key'],
            'refinement field/key identity')
    require(record['input_identity']['surface'] == identity, 'refinement source identity')
    result = copy.deepcopy(source)
    dimension = copy.deepcopy(record['claims']['sha2'])
    applicable = {'exact', 'bound'}
    for prior in (source.get('sha2'), prior_sha2):
        if prior is not None and prior['status'] in applicable:
            if dimension['status'] not in applicable:
                dimension = copy.deepcopy(prior)
            else:
                left = {prior['value']} if prior['status'] == 'exact' else set(prior['candidates'])
                right = {dimension['value']} if dimension['status'] == 'exact' else set(dimension['candidates'])
                values = sorted(left & right)
                require(values, 'conflicting Sha[2] dimension claims')
                evidence = []
                for claim in (dimension, prior):
                    for contributor in (claim['prior_evidence']
                                        if claim.get('method') == 'intersection' else [claim]):
                        if contributor not in evidence:
                            evidence.append(contributor)
                if len(evidence) == 1:
                    dimension = copy.deepcopy(evidence[0])
                else:
                    dimension = dict(status='exact' if len(values) == 1 else 'bound',
                                     value=values[0] if len(values) == 1 else None, candidates=values,
                                     method='intersection', prior_evidence=copy.deepcopy(evidence))
    result['sha2'] = dimension
    incoming_pairs = record['claims']['split']['candidates']
    require(all(Fraction(reg) * int(n) == Fraction(identity['sha_reg'])
                for n, reg in incoming_pairs), 'refinement Sha*Reg conflict')
    allowed = {(int(n), Fraction(reg)) for n, reg in incoming_pairs}
    pairs = [[n, str(Fraction(reg))] for n, reg in source['split']['candidates']
             if (n, Fraction(reg)) in allowed]
    if dimension['status'] in applicable:
        dims = [dimension['value']] if dimension['status'] == 'exact' else dimension['candidates']
        pairs = [[n, reg] for n, reg in pairs
                 if any(n % 2 == 1 if d == 0 else n % 2**d == 0 for d in dims)]
    require(pairs, 'empty split intersection or conflicting exact split')
    result['split'].update(candidates=pairs, proved=len(pairs) == 1,
                           sha=pairs[0][0] if len(pairs) == 1 else None,
                           regulator=pairs[0][1] if len(pairs) == 1 else None)
    evidence = dict(producer_id=record['producer_id'], evidence=record['evidence'])
    previous = result['split'].setdefault('refinement_evidence', [])
    if evidence not in previous:
        previous.append(evidence)
    return result
