#!/usr/bin/env python3
"""Reconstruct minimal GF(9) errata and both printed example enumerators.

Independent of verify/ and attack/ repository scripts. Uses the direct field
implementation in codesdeep/planar_classes.py for larger examples.
"""
import json
from collections import Counter
from math import ceil
from planar_classes import Field

def dist(E,a,level):
    f=[E.trace[E.mul[a][E.power(x,2)]] for x in range(E.q)]
    targets=([0] if level==0 else [a for a in range(1,E.p)
             if (pow(a,(E.p-1)//2,E.p)==1)==(level=='square')])
    D=[x for x in range(1,E.q) if f[x] in targets]
    h=Counter(sum(E.pair[b][x]!=0 for x in D) for b in range(E.q))
    z=h[0]
    assert z>0 and E.q%z==0
    ratio=E.q//z;dim=0
    while ratio>1:
        assert ratio%E.p==0
        ratio//=E.p;dim+=1
    return {'defining_set_indices':D,'n':len(D),'k':dim,'enumerator':{str(k):h[k] for k in sorted(h)},
            'total_messages':sum(h.values()),'zero_messages':h[0]}

def griesmer(p,k,d):return sum((d+p**i-1)//p**i for i in range(k))

def main():
    for p,m in ((3,2),(3,5),(5,3)):
        E=Field(p,m)
        print(json.dumps({'field':f'GF({p}^{m})','poly':E.poly,'primitive':E.primitive},sort_keys=True))
        a_values=[1,E.primitive] if m==2 else [E.primitive]
        levels=[0,'square','nonsquare'] if m==2 else ['square','nonsquare']
        for a in a_values:
            for level in levels:
                r=dist(E,a,level)
                print(json.dumps({'p':p,'m':m,'coefficient':a,'level':level,**r},sort_keys=True))
                if m==2 and a==1 and level==0:assert r['n']==4
                if m==2 and a==E.primitive and level==0:assert (r['n'],r['k'])==(0,0)
                if m==2 and a==1 and level in ('square','nonsquare'):assert (r['n'],r['k'])==(2,1)
                if (p,m,level)==(3,5,'nonsquare'):assert (r['n'],r['k'])==(90,5)
                if (p,m,level)==(3,5,'square'):assert (r['n'],r['k'])==(72,5)
                if (p,m,level)==(5,3,'nonsquare'):assert (r['n'],r['k'])==(60,3)
                if (p,m,level)==(5,3,'square'):assert (r['n'],r['k'])==(40,3)
    for printed in [{40:24,48:40,52:60},{28:40,32:60,40:24}]:
        print(json.dumps({'printed_weights':printed,'sum':1+sum(printed.values()),
                          'required_q_3_5':3**5,'matches_q_5_3':5**3},sort_keys=True))
        assert 1+sum(printed.values())==125!=243
    g=griesmer(3,4,7)
    print(json.dumps({'Griesmer_p3_k4_d7':g,'claimed_optimal_n':12,
                      'bound_excludes_d7':g>12},sort_keys=True))
    assert g==12
    print(json.dumps({'summary':'PASS','GF9_minimal_witness':'explicit integer basis and enumerated sets',
                      'printed_enumerators':'contradict 3^5 and match 5^3',
                      'griesmer':'insufficient for claimed optimality'},sort_keys=True))
if __name__=='__main__':main()
