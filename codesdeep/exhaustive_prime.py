#!/usr/bin/env python3
"""All normalized prime-field maps GF(3), GF(5), with exact Walsh signs.

Normalize F(0)=F(1)=0 by subtracting unique affine map; every planar
F is equivalent under that operation to exactly one listed normalized F.
"""
import itertools
import json
from collections import Counter
from planar_classes import base_gauss, abs_sq_is_q, equiv

def check(p):
    gauss=base_gauss(p,1)
    visited=0;planar=0;weak=0;failures=0;epsilon=Counter()
    for ys in itertools.product(range(p),repeat=p-2):
        vals=[0,0]+list(ys)
        visited+=1
        if any(len({(vals[(x+a)%p]-vals[x])%p for x in range(p)})!=p
               for a in range(1,p)):continue
        planar+=1
        for c in range(1,p):
            common=None
            for b in range(p):
                counts=[0]*p
                for x in range(p):counts[(c*vals[x]+b*x)%p]+=1
                assert abs_sq_is_q(counts,p,p)
                matches=[(sgn,j) for sgn in (-1,1) for j in range(p)
                         if equiv(counts,[sgn*gauss[(k-j)%p] for k in range(p)])]
                if len(matches)!=1 or (common is not None and common!=matches[0][0]):
                    failures+=1
                    raise AssertionError(('weakness',p,vals,c,b,counts,matches))
                common=matches[0][0]
            epsilon[common]+=1;weak+=1
    print(json.dumps({'p':p,'normalized_maps':visited,'planar_maps':planar,
                      'components':weak,'walsh_coefficients':weak*p,
                      'signs':dict(epsilon),'failures':failures},sort_keys=True))
    return planar,weak

def main():
    assert check(3)==(2,4)
    assert check(5)==(4,16)
    print(json.dumps({'summary':'PASS','coverage':'all maps GF(3) and GF(5) via affine normalization'},sort_keys=True))
if __name__=='__main__':main()
