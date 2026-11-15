#!/usr/bin/env python3
"""Complete sparse-family scan, not a scan over all planar maps.
Enumerate F(x)=x^2+a*x^d on GF(3^m), all a != 0, 3 <= d <= q-2,
for m=3,4. Planarity exact; test EVERY nonzero component coefficient via exact cyclotomic Walsh data.
This is only an exhaustive scan of the stated BINOMIAL FAMILY, not of
all planar maps in these fields.
"""
from collections import Counter
from coulter_matthews import Field,verify_walsh

def run(m):
    E=Field(m);q=E.q
    add=[[E.add(x,y) for y in range(q)] for x in range(q)]
    sub=[[E.sub(x,y) for y in range(q)] for x in range(q)]
    f2=[E.power(x,2) for x in range(q)]
    nr=0;ps=0;tested=0;weakpass=0;failing=[];dist=Counter()
    sample=[]
    for d in range(3,q-1):
        fd=[E.power(x,d) for x in range(q)]
        num=0
        for a in range(1,q):
            nr+=1
            vals=[add[f2[x]][E.mul(a,fd[x])] for x in range(q)]
            planar=True
            for delta in range(1,q):
                seen=bytearray(q)
                plus=add[delta]
                for x in range(q):
                    z=sub[vals[plus[x]]][vals[x]]
                    if seen[z]:planar=False;break
                    seen[z]=1
                if not planar:break
            if not planar:continue
            ps+=1;num+=1
            if len(sample)<8:sample.append((d,a))
            for c in range(1,q):
                fx=[E.trace[E.mul(c,v)] for v in vals]
                s,b,p=verify_walsh(E,fx)
                tested+=1
                if s is None:failing.append((d,a,c,b,p))
                else:weakpass+=1
        if num:dist[d]=num
    print(f'GF(3^{m}) q={q} poly={E.poly} primitive={E.generator} '
          f'sparse_cases={nr} planar_cases={ps} exact_Walsh_components_tested={tested} '
          f'weak_pass={weakpass} weak_fail={len(failing)}')
    print(f'  planar_counts_by_exponent={dict(dist)}')
    print(f'  first_planar_(exponent,a)={sample} failing_witnesses_first={failing[:3]}')
    assert not failing

if __name__=='__main__':
    run(3)
    run(4)
