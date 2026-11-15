#!/usr/bin/env python3
"""Pointwise check of Lemmas 2.7, 2.8, 2.11, 2.14 for Tr(a x^2).
Direct trace-zero intersections for every beta in the listed small fields.
"""
import json, math
import numpy as np
from verify_codes import Field, sign_from_counts

PARAMS=[(3,m) for m in range(2,8)]+[(5,m) for m in range(2,6)]+[(7,m) for m in (2,3,4)]+[(11,m) for m in (2,3)]+[(13,m) for m in (2,3)]

def check(g,a):
    p,m,q=g.p,g.m,g.q
    f=g.fvalues(a,2)
    eps=sign_from_counts(p,m,f)
    chi=-1 if p%4==3 else 1
    sq={i*i%p for i in range(1,p)}
    # Completing the square: f*(beta) = Tr(-beta^2/(4a)).
    inv=g.power(g.mul(4%p,a),q-2)
    c=g.mul((p-1),inv)
    dual=np.array([int(g.trace[g.mul(c,g.power(b,2))]) for b in range(q)],dtype=int)
    assert dual[0]==0
    dual_sign=sign_from_counts(p,m,dual)
    assert dual_sign==eps*chi**m
    if m%2==0:
        k=eps*chi**(m//2)*p**((m-2)//2)
        theory=[p**(m-1)+(p-1)*k if j==0 else p**(m-1)-k for j in range(p)]
    else:
        s=eps*chi**((m-1)//2)*p**((m-1)//2)
        theory=[p**(m-1)]+[p**(m-1)+(s if j in sq else -s) for j in range(1,p)]
    assert np.bincount(f,minlength=p).tolist()==theory,('lemma 2.7',p,m,a)
    # Note p=odd, beta=0 excluded in Lemmas 2.11/2.14.
    jf=g.log[np.flatnonzero(f==0)]
    js=g.log[np.flatnonzero(np.isin(f,list(sq)))]
    jn=g.log[np.flatnonzero((f!=0)&~np.isin(f,list(sq)))]
    # zero element must not occur among any of these, except f^-1(0).
    jf=jf[jf>=0]
    L=q-1
    tested=0
    for i in range(L):
        b=int(g.exp[i]); d=int(dual[b])
        n0=1+int(np.count_nonzero(g.tshift[(jf+i)%L]==0))
        ns=int(np.count_nonzero(g.tshift[(js+i)%L]==0))
        nn=int(np.count_nonzero(g.tshift[(jn+i)%L]==0))
        if m%2==0:
            zero_expected=p**(m-2)+(p-1)*k if d==0 else p**(m-2)
            common=(p-1)*(p**(m-2)-k)//2
            if d==0:es=en=common
            elif d in sq:es=(p-1)*(p**(m-2)+k)//2;en=common
            else:es=common;en=(p-1)*(p**(m-2)+k)//2
        else:
            v=eps*chi**((m+1)//2)*(p-1)*p**((m-3)//2)
            zero_expected=p**(m-2) if d==0 else p**(m-2)+(v if d in sq else -v)
            S=eps*chi**((m-1)//2)*p**((m-1)//2)
            T=eps*chi**((m-3)//2)*p**((m-3)//2)
            if d==0: es=(p-1)*(p**(m-2)+S)//2;en=(p-1)*(p**(m-2)-S)//2
            elif d in sq: es=en=(p-1)*(p**(m-2)-T)//2
            else: es=en=(p-1)*(p**(m-2)+T)//2
        assert (n0,ns,nn)==(zero_expected,es,en),('lemmas 2.11/2.14',p,m,a,b,d,(n0,ns,nn),(zero_expected,es,en))
        assert n0+ns+nn==p**(m-1)
        tested+=1
    return {'p':p,'m':m,'a':a,'epsilon':eps,'dual_sign':dual_sign,
            'beta_checked':tested,'status':'CONFIRMED'}

def main():
    for p,m in PARAMS:
        g=Field(p,m)
        for a in (1,g.g):
            print(json.dumps(check(g,a),sort_keys=True),flush=True)
    print(json.dumps({'status':'COMPLETE','type':'done'}),flush=True)

if __name__=='__main__':main()
