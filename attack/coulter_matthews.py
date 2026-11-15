#!/usr/bin/env python3
"""Exact arithmetic: Coulter--Matthews planar monomial x^14 over GF(3^m).

No Sage, Magma, sympy, numpy or other third-party dependencies.
- Find an irreducible polynomial, choose a primitive element, build log tables.
- Verify planarity (all nonzero derivatives are permutations).
- Verify all Walsh coefficients in Z[zeta_3] EXACTLY against
    epsilon*(sqrt(-3))^m*zeta_3^j with one common epsilon per component.
- Compute trace-defined code weights directly, and its projective minimality.
- Compare Table I (m even) / Table II (m odd) from Tang et al.

m=4,5: all c != 0. m=7: two representatives c=1, primitive.
For exponent 14, gcd(14,3^m-1)=2 in these three degrees, so the two
representatives suffice for every c by x -> t*x (see REPORT.md proof).
Usage: python3 coulter_matthews.py
"""
import math
from collections import Counter
from itertools import product

P=3

def digits(v,m):
    out=[]
    for _ in range(m): out.append(v%3);v//=3
    return out

def polynomial_rem(poly,divisor):
    v=list(poly)
    while len(v)>=len(divisor):
        if v[-1]:
            factor=v[-1]
            offset=len(v)-len(divisor)
            for i,coef in enumerate(divisor):
                v[offset+i]=(v[offset+i]-factor*coef)%3
        while v and v[-1]==0:v.pop()
    return v

def irreducible(m):
    for coeff in product(range(3),repeat=m):
        if coeff[0]==0:continue
        f=list(coeff)+[1]
        for d in range(1,m//2+1):
            if any(not polynomial_rem(f,list(g)+[1])
                   for g in product(range(3),repeat=d)):
                break
        else:return f
    raise ValueError('No polynomial')

def prime_factors(n):
    out=[]
    f=2
    while f*f<=n:
        if n%f==0:
            out.append(f)
            while n%f==0:n//=f
        f+=1
    if n>1:out.append(n)
    return out

class Field:
    def __init__(self,m):
        self.m=m;self.q=3**m;self.poly=irreducible(m)
        self.D=[digits(i,m) for i in range(self.q)]
        self.pow3=[3**i for i in range(m)]
        factors=prime_factors(self.q-1)
        self.generator=next(a for a in range(2,self.q)
                            if all(self.slowpow(a,(self.q-1)//d)!=1 for d in factors))
        self.exp=[1]*(self.q-1)
        self.log=[-1]*self.q
        a=1
        for k in range(self.q-1):
            self.exp[k]=a
            assert self.log[a]==-1
            self.log[a]=k
            a=self.slowmul(a,self.generator)
        assert a==1 and self.log[0]==-1
        self.trace=[self.calc_trace(a) for a in range(self.q)]
        assert set(self.trace)=={0,1,2}
        # Bilinear pairing table: each row a bytearray indexed by x.
        self.pair=[bytearray(self.q) for _ in range(self.q)]
        for beta in range(self.q):
            row=self.pair[beta]
            for x in range(self.q):row[x]=self.trace[self.mul(beta,x)]

    def add(self,a,b):
        return sum(((x+y)%3)*t for x,y,t in zip(self.D[a],self.D[b],self.pow3))
    def sub(self,a,b):
        return sum(((x-y)%3)*t for x,y,t in zip(self.D[a],self.D[b],self.pow3))
    def slowmul(self,a,b):
        m=self.m;acc=[0]*(2*m-1);da=self.D[a];db=self.D[b]
        for i,ai in enumerate(da):
            if ai:
                for j,bj in enumerate(db):acc[i+j]=(acc[i+j]+ai*bj)%3
        for k in range(2*m-2,m-1,-1):
            x=acc[k]
            if x:
                for j in range(m):acc[k-m+j]=(acc[k-m+j]-x*self.poly[j])%3
        return sum(acc[i]*self.pow3[i] for i in range(m))
    def slowpow(self,a,k):
        v=1
        while k:
            if k&1:v=self.slowmul(v,a)
            a=self.slowmul(a,a);k//=2
        return v
    def mul(self,a,b):
        if not(a and b):return 0
        return self.exp[(self.log[a]+self.log[b])%(self.q-1)]
    def power(self,a,k):
        return 0 if a==0 else self.exp[(self.log[a]*k)%(self.q-1)]
    def calc_trace(self,a):
        z=0;b=a
        for _ in range(self.m):z=self.add(z,b);b=self.slowpow(b,3)
        assert z<3
        return z


def is_planar(E,F):
    q=E.q
    # Equivalent exact derivative permutation definition.
    for a in range(1,q):
        seen=bytearray(q)
        for x in range(q):
            d=E.sub(F[E.add(x,a)],F[x])
            if seen[d]:return False,a,x,d
            seen[d]=1
    return True,None,None,None


def sqrt_pstar_coeff(m):
    # (sqrt(-3))^m = A+B*zeta_3 since sqrt(-3)=1+2*zeta_3.
    k=(-3)**(m//2)
    return (k,0) if m%2==0 else (k,2*k)

def rot(A,B):
    # Multiply by zeta: zeta^2=-1-zeta
    return (-B,A-B)

def walsh_expected(m):
    A,B=sqrt_pstar_coeff(m)
    pos=[]
    for j in range(3):
        pos.append((A,B));A,B=rot(A,B)
    return {v:(sgn,j) for sgn in (1,-1) for j,(a,b) in enumerate(pos)
            for v in [(sgn*a,sgn*b)]}


def verify_walsh(E,fx):
    q=E.q;table=walsh_expected(E.m)
    sign=None
    for beta in range(q):
        cnt=[0,0,0]
        brow=E.pair[beta]
        for x in range(q):cnt[(fx[x]+brow[x])%3]+=1
        pair=(cnt[0]-cnt[2],cnt[1]-cnt[2])
        if pair not in table:return None,beta,pair
        sg,j=table[pair]
        if sign is not None and sign!=sg:return None,beta,pair
        sign=sg
    return sign,None,None


def code_distribution(E,fx):
    D=[x for x in range(1,E.q) if fx[x]==0]
    hist=Counter(sum(E.pair[beta][x]!=0 for x in D) for beta in range(E.q))
    return D,hist

def expected_distribution(m,eps):
    p=3
    if m%2==0:
        k=(-1)**(m//2)
        a=(p-1)*p**(m-2)
        b=(p-1)*(p**(m-2)+eps*k*p**((m-2)//2))
        return Counter({0:1,a:p**(m-1)-1+eps*k*(p-1)*p**((m-2)//2),
                        b:(p-1)*(p**(m-1)-eps*k*p**((m-2)//2))})
    a=(p-1)*p**(m-2)
    return Counter({0:1,a:p**(m-1)-1,
                    (p-1)*(p**(m-2)-p**((m-3)//2)):(p-1)//2*(p**(m-1)+p**((m-1)//2)),
                    (p-1)*(p**(m-2)+p**((m-3)//2)):(p-1)//2*(p**(m-1)-p**((m-1)//2))})

def check_minimal(E,D):
    # One codeword for each projective beta in GF(3)^m (beta and -beta).
    supp={}
    for beta in range(1,E.q):
        neg=E.sub(0,beta)
        if beta>neg:continue
        mask=sum(1<<i for i,x in enumerate(D) if E.pair[beta][x])
        supp[beta]=mask
    for a,A in supp.items():
        for b,B in supp.items():
            if a!=b and (A&~B)==0:return False,(a,b,A.bit_count(),B.bit_count()),len(supp)
    return True,None,len(supp)


def main():
    for m in (4,5,7):
        E=Field(m);q=E.q
        F=[E.power(x,14) for x in range(q)]
        planar,*details=is_planar(E,F)
        assert planar,details
        print(f'GF(3^{m}) q={q} poly={E.poly} primitive={E.generator} x^14 planar={planar} all_derivatives={q-1}')
        cset=list(range(1,q)) if m in (4,5) else [1,E.generator]
        signcounts=Counter();distcounts=Counter();mincounts=Counter()
        for c in cset:
            fx=[E.trace[E.mul(c,v)] for v in F]
            eps,beta,pair=verify_walsh(E,fx)
            assert eps is not None,(m,c,beta,pair)
            # Coulter-Matthews sign: (-1)^(m-1)*eta(c).
            eta=1 if E.power(c,(q-1)//2)==1 else -1
            assert eps==((-1)**(m-1))*eta,(m,c,eps,eta)
            signcounts[eps]+=1
            D,observed=code_distribution(E,fx)
            expected=expected_distribution(m,eps)
            assert observed==expected,(m,c,observed,expected)
            distcounts[tuple(sorted(observed.items()))]+=1
            minimal,witness,num=check_minimal(E,D)
            mincounts[(eps,minimal)]+=1
            if c in (1,E.generator):
                print(f'  c={c} eta={eta} epsilon={eps} D_length={len(D)} '+
                      f'weights={sorted(observed.items())} projective_classes={num} '
                      f'minimal={minimal} containment_witness={witness}')
        print(f'  tested_coefficients={len(cset)} epsilon_counts={dict(sorted(signcounts.items()))} '+
              f'weight_enumerator_types={len(distcounts)} minimality_counts={dict(sorted(mincounts.items()))} '+
              f'all_exact_Walsh_and_tables=PASS')
        if m==7:
            assert math.gcd(14,q-1)==2
            print('  q=2187: two square-class representatives only; gcd(14,q-1)=2')

if __name__=='__main__':main()
