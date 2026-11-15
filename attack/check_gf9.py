#!/usr/bin/env python3
"""Independent checker/enumerator for exhaustive_planar.c.

Represents GF(9) as pairs (a,b) with u^2=-1, NOT as C's encoded
integer tables. Backtracking prunes impossible partial derivative rows.
The Walsh test uses complex exponentials and tests a *single common*
phase sign, independent of the C program's integer pair templates.
Usage: python3 check_gf9.py
"""
import cmath
from collections import Counter
from itertools import product

P = 3
ELS = list(product(range(P), repeat=2))
ZERO = (0, 0)
ONE = (1, 0)

def add(a,b): return ((a[0]+b[0])%3, (a[1]+b[1])%3)
def sub(a,b): return ((a[0]-b[0])%3, (a[1]-b[1])%3)
def mul(a,b): return ((a[0]*b[0]-a[1]*b[1])%3, (a[0]*b[1]+a[1]*b[0])%3)
def tr(a): return (2*a[0])%3

# Convert between pair representation, and the digit index in C, only for
# the final digest and the published witness.
def n(a): return a[0]+3*a[1]
R = sorted(ELS,key=n)
A = {x:i for i,x in enumerate(R)}
ZETA = cmath.exp(2j*cmath.pi/3)


def is_planar(F):
    for a in R[1:]:
        differences = {sub(F[add(x,a)],F[x]) for x in R}
        if len(differences)!=9: return False
    return True


def weak_sign(F,c):
    # The p*= -3, m=2 transform has reference amplitude -3.
    signs=set()
    for beta in R:
        z=sum(ZETA ** ((tr(mul(c,F[x]))+tr(mul(beta,x)))%3) for x in R)
        z/= -3
        matches=[sgn for sgn in (-1,1) for j in range(3) if abs(z-sgn*ZETA**j)<1.e-8]
        if len(matches)!=1: return None
        signs.add(matches[0])
    return next(iter(signs)) if len(signs)==1 else None


def partial_ok(F):
    # A derivative D_aF(x) is assigned iff both endpoints are known.
    for a in R[1:]:
        vals=set()
        for x in R:
            y=add(x,a)
            if x in F and y in F:
                d=sub(F[y],F[x])
                if d in vals:return False
                vals.add(d)
    return True


def exhaustive():
    F={ZERO:ZERO,ONE:ZERO}
    global checked,planars,digest,positive,negative,failures
    checked=9**7;planars=positive=negative=failures=0
    digest=14695981039346656037
    def dfs(i):
        global planars,digest,positive,negative,failures
        if i==9:
            assert is_planar(F)
            planars+=1
            for x in R:
                digest ^= n(F[x])
                digest=(digest*1099511628211)&((1<<64)-1)
            for c in R[1:]:
                s=weak_sign(F,c)
                if s is None:failures+=1
                elif s==1:positive+=1
                else:negative+=1
            return
        x=R[i]
        for v in R:
            F[x]=v
            if partial_ok(F):dfs(i+1)
            del F[x]
    dfs(2)
    print(f'INDEPENDENT exhaustive GF(9): checked={checked} planar={planars} '
          f'components={positive+negative+failures} weak_failures={failures} '
          f'positive_epsilon={positive} negative_epsilon={negative} fnv64={digest:016x}')


def check_witness():
    witness=[0,0,1,0,3,7,2,8,3]
    F={x:R[v] for x,v in zip(R,witness)}
    assert is_planar(F)
    print('EXPLICIT first_planar witness:',witness)
    for c in R[1:]:
        print(f'  c={n(c)} weak_sign_epsilon={weak_sign(F,c)}')

if __name__=='__main__':
    check_witness()
    exhaustive()
