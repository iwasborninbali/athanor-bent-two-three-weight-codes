#!/usr/bin/env python3
"""Independent GF(81) minimality-containment witness checker.
Uses direct polynomial reduction for every operation, with no use of the
logarithm/trace tables from coulter_matthews.py.
Usage: python3 verify_witness.py attack/minimality_witness.json
"""
import json
import sys


def digits(x):
    a=[]
    for _ in range(4):a.append(x%3);x//=3
    assert x==0
    return a

def num(a):return sum(d*3**i for i,d in enumerate(a))

def add(x,y):return num([(a+b)%3 for a,b in zip(digits(x),digits(y))])

def multiply(x,y,poly):
    a=digits(x);b=digits(y)
    t=[0]*7
    for i in range(4):
        for j in range(4):t[i+j]=(t[i+j]+a[i]*b[j])%3
    for k in range(6,3,-1):
        top=t[k]
        for j in range(4):t[k-4+j]=(t[k-4+j]-top*poly[j])%3
        t[k]=0
    return num(t[:4])

def power(x,n,poly):
    acc=1
    while n:
        if n&1:acc=multiply(acc,x,poly)
        x=multiply(x,x,poly);n//=2
    return acc

def trace(x,poly):
    s=0
    for _ in range(4):
        s=add(s,x);x=power(x,3,poly)
    assert s<3
    return s

def compute(poly,c,betas):
    D=[x for x in range(1,81) if trace(multiply(c,power(x,14,poly),poly),poly)==0]
    words=[[trace(multiply(beta,x,poly),poly) for x in D] for beta in betas]
    return D,words

def main(path):
    obj=json.load(open(path,encoding='utf-8'))
    poly=obj['field_modulus_coefficients_low_first']
    assert poly==[1,0,1,1,1]
    c=obj['c'];betas=obj['betas']
    assert c==1 and betas==[1,10]
    D,words=compute(poly,c,betas)
    assert D==obj['defining_set_indices']
    assert words==obj['codewords']
    supports=[{i for i,v in enumerate(w) if v} for w in words]
    assert supports[0]<supports[1]
    assert len(supports[0])==12 and len(supports[1])==18
    assert betas[0]!=betas[1] and betas[0]!=81-betas[1]
    # GF(3) scalar multiples are 1*beta or 2*beta (2 encoded as constant).
    assert all(betas[0]!=multiply(a,betas[1],poly) for a in (1,2))
    print('WITNESS VERIFIED GF(81), modulus=t^4+t^3+t^2+1, c=1, exponent=14')
    print(f'D_size={len(D)} D={D}')
    print(f'betas={betas} weights={[len(s) for s in supports]}')
    print(f'support_positions_0based={[sorted(s) for s in supports]}')
    print('proper_support_inclusion=True beta_not_proportional=True nonminimal=True')

if __name__=='__main__':
    assert len(sys.argv)==2
    main(sys.argv[1])
