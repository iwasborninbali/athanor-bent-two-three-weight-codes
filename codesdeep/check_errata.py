#!/usr/bin/env python3
"""Independent tuple-field checker of minimal errata witnesses and corrected examples.

Does NOT import planar_classes.py or errata.py; no third-party packages.
"""
import json
from collections import Counter

FIELDS=[(3,2,[1,0,1],(1,1)),
        (3,5,[1,0,0,0,2,1],(0,1,0,0,0)),
        (5,3,[1,0,1,1],(2,1,0))]

def reproduce(p,m,mod,w):
    q=p**m
    zero=(0,)*m
    one=(1,)+(0,)*(m-1)
    es=[tuple((i//p**j)%p for j in range(m)) for i in range(q)]
    def add(a,b):return tuple((x+y)%p for x,y in zip(a,b))
    def mul(a,b):
        acc=[0]*(2*m-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):acc[i+j]+=x*y
        for k in range(2*m-2,m-1,-1):
            c=acc[k]%p
            for j in range(m):acc[k-m+j]-=c*mod[j]
        return tuple(v%p for v in acc[:m])
    def power(a,k):
        out=one
        while k:
            if k&1:out=mul(out,a)
            a=mul(a,a);k>>=1
        return out
    def tr(a):
        out=zero
        for k in range(m):out=add(out,power(a,p**k))
        assert all(v==0 for v in out[1:])
        return out[0]
    assert all(power(a,q-1)==one for a in es if a!=zero)
    classes={0:{0},'square':{i for i in range(1,p) if pow(i,(p-1)//2,p)==1},
             'nonsquare':{i for i in range(1,p) if pow(i,(p-1)//2,p)!=1}}
    cs=[one,w] if m==2 else [w]
    levels=(0,'square','nonsquare') if m==2 else ('square','nonsquare')
    for c in cs:
        trace_values={x:tr(mul(c,mul(x,x))) for x in es}
        for lev in levels:
            d=[x for x in es if x!=zero and trace_values[x] in classes[lev]]
            count=Counter(sum(tr(mul(beta,x))!=0 for x in d) for beta in es)
            z=count[0];size=q//z;rank=0
            while size>1:
                assert size%p==0
                size//=p;rank+=1
            assert sum(count.values())==q
            obj={'p':p,'m':m,'coefficient':list(c),'level':lev,'n':len(d),
                 'k':rank,'enumerator':{str(k):count[k] for k in sorted(count)}}
            print(json.dumps(obj,sort_keys=True))
            if (p,m,c,lev)==(3,2,one,0):assert (len(d),rank)==(4,2)
            if p==3 and m==2 and c==w and lev==0:assert (len(d),rank)==(0,0)
            if p==3 and m==2 and c==one and lev!='0' and lev!=0:assert (len(d),rank)==(2,1)
            if p==3 and m==5 and lev=='square':assert count==Counter({0:1,42:72,48:90,54:80})
            if p==3 and m==5 and lev=='nonsquare':assert count==Counter({0:1,54:80,60:72,66:90})
            if p==5 and m==3 and lev=='square':assert count==Counter({0:1,28:40,32:60,40:24})
            if p==5 and m==3 and lev=='nonsquare':assert count==Counter({0:1,40:24,48:40,52:60})

def main():
    for p,m,mod,w in FIELDS:reproduce(p,m,mod,w)
    assert sum((7+3**i-1)//3**i for i in range(4))==12
    assert 1+24+40+60==125!=3**5
    assert 1+40+60+24==125!=3**5
    print(json.dumps({'independent_tuple_arithmetic':'PASS',
                      'checked':'GF9 endpoint; GF(3^5), GF(5^3) complete code enumerators; sums; Griesmer'},sort_keys=True))
if __name__=='__main__':main()
