import collections
import json
import math
import sys
import time
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
MW = HERE.parent
sys.path.insert(0,str(MW/'hard/code'))
from lattice import fibre_data, _det

GAMMA_POWERS = {1:F(1),2:F(4,3),3:F(2),4:F(4),5:F(8)}

def square_integer(x):
    return x.denominator==1 and x>=0 and math.isqrt(x.numerator)**2==x.numerator

from functools import lru_cache

@lru_cache(maxsize=4096)
def cached_corrections(fibres):
    cs={F(0)}
    for item in fibres.split(';'):
        place,e,typ,c,vd,cond=item.split(':')
        cs={x+y for x in cs for y in fibre_data(typ,int(c),int(e))[2]}
    return cs

def corrections(row):
    return cached_corrections(row['fibres'])

def min_height(cs, nmin=0):
    return min(4+2*max(nmin,math.floor((c-4)/2)+1)-c for c in cs)

def invariant_sha(row):
    cs=corrections(row);r=row['rank'];hmin=min_height(cs);S=F(row['sha_reg'])
    return [n for n in row['sha'] if S/n>=hmin**r/GAMMA_POWERS.get(r, F(r**r))]

def pos(height,cs):
    return sorted({int((height+c-4)/2) for c in cs if height+c>=4 and ((height+c-4)/2).denominator==1})

def filters(row,native):
    start=time.process_time_ns()
    r=row['rank']; G=[[F(x) for x in rr] for rr in native['gram']];k=len(G)
    S=F(row['sha_reg']); cs=corrections(row);det=_det(G);b=row['b']
    assert r==native['rank'] and S==F(native['sha_reg'])
    assert k<=r
    H=min_height(cs,b+1);mu=F(native['mu'])
    sha=row['sha']
    inv=invariant_sha(row)
    height=[n for n in sha if k!=r or square_integer(det*n/S)]
    cert=[n for n in sha if not native['certified'] or S/n==det]
    negative=[];evidence={}
    for n in sha:
        R=S/n;why=[]
        if r==1:
            ps=pos(R,cs)
            assert ps
            if max(ps)<=b and (k==0 or R<G[0][0]):
                why.append('rank_one_generator_absent')
        if k<r and mu<H:
            lower=det*(H-mu)**(r-k)/GAMMA_POWERS.get(r-k, F((r-k)**(r-k)))
            if R<lower:
                why.append('quotient_Hermite_bound')
        if why:evidence[str(n)]=why
        else:negative.append(n)
    combined=[n for n in height if n in cert and n in negative and n in inv]
    assert combined,(row,native,evidence)
    return dict(q=row['q'],key=row['key'],rank=r,found_rank=k,b=b,initial_sha=sha,initial_pairs=[[n,str(S/n)] for n in sha],invariant_sha=inv,height_sha=height,certificate_sha=cert,negative_sha=negative,combined_sha=combined,combined_pairs=[[n,str(S/n)] for n in combined],det=str(det),mu=str(mu),missing_height_lower_bound=str(H),corrections=list(map(str,sorted(cs))),negative_evidence=evidence,native=native,filter_ms=(time.process_time_ns()-start)/1e6)

def main():
    rows=[json.loads(s) for s in (HERE/'population.jsonl').read_text().splitlines()]
    out=[]
    for q in (2,3):
        native={r['key']:r for r in map(json.loads,(HERE/f'bounded-f{q}.jsonl').read_text().splitlines())}
        rs=[r for r in rows if r['q']==q]
        assert len(rs)==len(native)
        out.extend(filters(r,native[r['key']]) for r in rs)
    (HERE/'bounded-analysis.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in out))
    for q in (2,3):
        rs=[r for r in out if r['q']==q]
        print(q,'initial',len(rs),'found ranks',collections.Counter((r['rank'],r['found_rank']) for r in rs))
        for method in ['invariant','height','certificate','negative','combined']:
            rem=[r for r in rs if len(r[method+'_sha'])>1]
            print(method,'remaining',len(rem),'ranks',collections.Counter(r['rank'] for r in rem))
        print('native ms',sum(sum(r['native']['ms'])+r['native']['search_ms'] for r in rs)/len(rs))
if __name__=='__main__':main()
