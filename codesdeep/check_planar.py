#!/usr/bin/env python3
"""Independent tuple-polynomial, direct complex-Walsh cross-check at c=1.

Uses none of planar_classes.py's field code: separate tuple arithmetic,
independent exponentiation, full derivatives and every Walsh index for c=1.
Numerical phase verification is supplementary to exact main checker.
"""
import cmath
import json
from math import pi

# class labels, coefficients and irreducible moduli as independently documented
# by the primary log. Recompute arithmetic, NOT imported from primary code.
INPUTS = [
(3,1,[1,1], [('prime_x2',[(1,2)])]),
(3,2,[1,0,1], [('2.1_FF',[(1,2)])]),
(3,3,[1,0,2,1], [('3.1_FF',[(1,2)]),('3.2_Albert',[(1,4)])]),
(3,4,[1,0,1,1,1], [('4.1_FF',[(1,2)]),('4.2_CM',[(1,14)]),('4.3_Dickson_BH_ZP',[(1,36),(2,10),(2,4)])]),
(3,5,[1,0,0,0,2,1], [('5.1_FF',[(1,2)]),('5.2_Albert',[(1,4)]),('5.3_Albert',[(1,10)]),('5.4_CM',[(1,14)]),('5.5_TST_DingYuan',[(1,10),(1,6),(2,2)]),('5.6_TST_DingYuan',[(1,10),(2,6),(2,2)]),('5.7_ACW',[(1,90),(1,2)]),('5.8_CK',[(1,162),(1,108),(2,84),(1,2)])]),
(5,1,[1,1], [('prime_x2',[(1,2)])]),
(5,2,[1,1,1], [('p2_FF',[(1,2)])]),
(5,3,[1,0,1,1], [('p3_FF',[(1,2)]),('p3_Albert',[(1,6)])]),
]

def crosscheck(p,m,mod,classes):
    q=p**m
    elements=[tuple((v//p**i)%p for i in range(m)) for v in range(q)]
    index={v:k for k,v in enumerate(elements)}
    zero=(0,)*m;one=(1,)+(0,)*(m-1)
    def plus(a,b):return tuple((x+y)%p for x,y in zip(a,b))
    def minus(a,b):return tuple((x-y)%p for x,y in zip(a,b))
    def times(a,b):
        out=[0]*(2*m-1)
        for i,u in enumerate(a):
            if not u:continue
            for j,v in enumerate(b):out[i+j]+=u*v
        for i in range(2*m-2,m-1,-1):
            t=out[i]%p
            if t:
                for j in range(m):out[i-m+j]-=t*mod[j]
        return tuple(out[i]%p for i in range(m))
    def power(a,e):
        v=one
        while e:
            if e&1:v=times(v,a)
            a=times(a,a);e//=2
        return v
    def trace(a):
        s=zero;v=a
        for i in range(m):
            s=plus(s,v);v=power(v,p)
        assert not any(s[1:])
        return s[0]
    assert all(power(a,q-1)==one for a in elements if a!=zero)
    # The next three arrays arise from tuple-based arithmetic, not from
    # planar_classes Field/logarithm tables.
    ADD=[[index[plus(a,b)] for b in elements] for a in elements]
    SUB=[[index[minus(a,b)] for b in elements] for a in elements]
    MUL=[[index[times(a,b)] for b in elements] for a in elements]
    TR=[trace(x) for x in elements]
    PAIR=[[TR[MUL[b][x]] for x in range(q)] for b in range(q)]
    roots=[cmath.exp(2j*pi*j/p) for j in range(p)]
    for name,terms in classes:
        F=[]
        for x in elements:
            y=zero
            for coefficient,degree in terms:
                term=power(x,degree)
                for _ in range(coefficient):y=plus(y,term)
            F.append(index[y])
        for a in range(1,q):
            ds={SUB[F[ADD[a][x]]][F[x]] for x in range(q)}
            assert len(ds)==q,(p,m,name,'planarity',a)
        # Evaluate every Walsh point for scalar component c=1 by complex sum,
        # independently of primary integer cyclotomic character checks.
        phases=[TR[F[x]] for x in range(q)]
        W=[sum(roots[(phases[x]+PAIR[b][x])%p] for x in range(q))
           for b in range(q)]
        base=W[0]
        assert abs(abs(base)**2-q)<q*1e-8,(p,m,name,'W0 magnitude')
        err=max(min(abs(w-base*root) for root in roots) for w in W)
        assert err<1e-7,(p,m,name,'phase',err)
        print(json.dumps({'p':p,'m':m,'class':name,'derivatives_checked':q-1,
                          'checked_component':1,'walsh_indices':q,
                          'maximum_complex_phase_residual':round(err,12),
                          'status':'PASS'},sort_keys=True))

def main():
    num=0;wcount=0
    for p,m,mod,classes in INPUTS:
        crosscheck(p,m,mod,classes)
        num+=len(classes);wcount+=len(classes)*p**m
    print(json.dumps({'independent_representatives':num,'c1_walsh_coefficients':wcount,
                      'result':'PASS','scope':'one component per class, not all c'},sort_keys=True))
if __name__=='__main__':main()
