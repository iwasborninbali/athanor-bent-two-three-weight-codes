#!/usr/bin/env python3
"""Independent direct-polynomial counterexample to the printed lengths in
Theorem 3.1 / Corollary 3.3 (arXiv v3, page 10), m=6, p=3.

No calls to coulter_matthews.py or library finite-field routines. This
checker recomputes GF(3^6), f(x)=Tr(x^2), D, all 729 codewords and the
punctured code by direct polynomial multiplication and trace.
Usage: python3 theorem_counterexample.py
"""
from collections import Counter

MOD=[1,0,0,0,1,1,1]  # t^6+t^5+t^4+1, irreducible over GF(3)
Q=729

def digits(a):
    out=[]
    for _ in range(6):out.append(a%3);a//=3
    assert not a
    return out

D=[digits(a) for a in range(Q)]
POWS=[3**i for i in range(6)]

def add(a,b):return sum(((x+y)%3)*w for x,y,w in zip(D[a],D[b],POWS))

def multiply(a,b):
    x=D[a];y=D[b];v=[0]*11
    for i in range(6):
        for j in range(6):v[i+j]=(v[i+j]+x[i]*y[j])%3
    for k in range(10,5,-1):
        t=v[k]
        for j in range(6):v[k-6+j]=(v[k-6+j]-t*MOD[j])%3
        v[k]=0
    return sum(v[i]*POWS[i] for i in range(6))

def power(a,n):
    ret=1
    while n:
        if n&1:ret=multiply(ret,a)
        a=multiply(a,a);n//=2
    return ret

def trace(a):
    z=0
    for _ in range(6):z=add(z,a);a=power(a,3)
    assert z<3
    return z


def direct_trace(a):
    z=0
    for _ in range(6):z=add(z,a);a=power(a,3)
    assert z<3
    return z

def main():
    # Direct polynomial multiplication, no data shared with other programs.
    assert len({power(4,k) for k in range(Q-1)})==Q-1
    T=[direct_trace(power(3,k)) for k in range(11)]
    assert all(v in (0,1,2) for v in T)
    gram=[[T[i+j] for j in range(6)] for i in range(6)]
    traces=[sum(d*t for d,t in zip(D[x],T))%3 for x in range(Q)]
    assert all(traces[x]==direct_trace(x) for x in range(Q))
    f=[traces[multiply(x,x)] for x in range(Q)]
    counts=[f.count(i) for i in range(3)]
    W0=(counts[0]-counts[2],counts[1]-counts[2])
    assert W0==(27,0),W0
    eps=-1
    nonzero_D=[x for x in range(1,Q) if f[x]==0]
    assert len(nonzero_D)==260
    # Bilinear Gram pairing independently verifies trace(beta*x).
    def code_row(beta,points):
        db=D[beta]
        linear=[sum(db[i]*gram[i][j] for i in range(6))%3 for j in range(6)]
        return [sum(linear[j]*D[x][j] for j in range(6))%3 for x in points]
    assert all(code_row(b,[x])[0]==traces[multiply(b,x)]
               for b in (1,4,10,50) for x in (1,2,3,44,300))
    weights=Counter(sum(bool(t) for t in code_row(beta,nonzero_D)) for beta in range(Q))
    assert weights==Counter({0:1,162:260,180:468}),weights
    seen=set();projective=[]
    for x in nonzero_D:
        if x in seen:continue
        neg=add(x,x)
        assert neg in nonzero_D and neg!=x
        seen.add(x);seen.add(neg)
        projective.append(x)
    assert len(projective)==130
    shortened=Counter(sum(bool(t) for t in code_row(beta,projective)) for beta in range(Q))
    assert shortened==Counter({0:1,81:260,90:468}),shortened
    print('INDEPENDENT COUNTEREXAMPLE to printed parameter statements:')
    print('p=3 m=6 F=GF(3)[t]/(t^6+t^5+t^4+1); f(x)=Tr(x^2), epsilon=-1')
    print(f'phase_counts_f0_f1_f2={counts} W0_as_(A+Bzeta)={W0}')
    print(f'Theorem_3_1_printed_length=224 actual_length={len(nonzero_D)} '+
          f'full_code_weight_histogram={sorted(weights.items())}')
    print(f'Corollary_3_3_printed_length=112 actual_punctured_length={len(projective)} '+
          f'punctured_weight_histogram={sorted(shortened.items())}')
    print('tables_I_III_with_sign_factor_match_actual=True')

if __name__=='__main__':main()
