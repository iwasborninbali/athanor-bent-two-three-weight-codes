#!/usr/bin/env python3
"""Exact exhaustive scalar-component spectra for published small-field EA/CCZ representatives.

No external libraries; deterministic polynomial-basis finite fields and integer
cyclotomic count vectors. Checks ALL nonzero components and ALL Walsh indices
for each listed representative; not all planar maps/equivalence classes.
"""
import itertools
import json
from collections import Counter

# Andreoli--Budaghyan--Coulter--Haukenes--Kaleyski--Piccione (2025),
# Appendix A Table 1, n=2,...,5. p=5 cases from Menichetti classification
# of commutative semifields of order p^3, not of arbitrary planar maps.
CLASSES = {
    (3,1): [('prime_x2', [(1,2)])],
    (3,2): [('2.1_FF', [(1,2)])],
    (3,3): [('3.1_FF', [(1,2)]), ('3.2_Albert', [(1,4)])],
    (3,4): [('4.1_FF', [(1,2)]),('4.2_CM',[(1,14)]),
            ('4.3_Dickson_BH_ZP',[(1,36),(2,10),(2,4)])],
    (3,5): [('5.1_FF',[(1,2)]),('5.2_Albert',[(1,4)]),
            ('5.3_Albert',[(1,10)]),('5.4_CM',[(1,14)]),
            ('5.5_TST_DingYuan',[(1,10),(1,6),(2,2)]),
            ('5.6_TST_DingYuan',[(1,10),(2,6),(2,2)]),
            ('5.7_ACW',[(1,90),(1,2)]),
            ('5.8_CK',[(1,162),(1,108),(2,84),(1,2)])],
    (5,1): [('prime_x2',[(1,2)])],
    (5,2): [('p2_FF',[(1,2)])],
    (5,3): [('p3_FF',[(1,2)]),('p3_Albert',[(1,6)])]
}

def digits(v,p,m):
    out=[]
    for _ in range(m):
        out.append(v%p);v//=p
    return out

def remainder(poly,div,p):
    r=list(poly)
    while r and r[-1]==0:r.pop()
    while len(r)>=len(div):
        c=r[-1];j=len(r)-len(div)
        if c:
            for k,v in enumerate(div):r[j+k]=(r[j+k]-c*v)%p
        while r and r[-1]==0:r.pop()
    return r

def irred(p,m):
    for cs in itertools.product(range(p), repeat=m):
        if cs[0]==0:continue
        f=list(cs)+[1]
        if all(remainder(f,list(g)+[1],p) for d in range(1,m//2+1)
               for g in itertools.product(range(p),repeat=d)):
            return f
    raise AssertionError('missing irreducible polynomial')

def convolution(a,b):
    p=len(a);out=[0]*p
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[(i+j)%p]+=x*y
    return out

def base_gauss(p,m):
    residues={x*x%p for x in range(1,p)}
    g=[0]+[1 if i in residues else -1 for i in range(1,p)]
    out=[1]+[0]*(p-1)
    for _ in range(m):out=convolution(out,g)
    return out

def equiv(a,b):
    # eval at primitive pth root equal iff a-b constant (multiple of Phi_p)
    delta=[u-v for u,v in zip(a,b)]
    return len(set(delta))==1

def abs_sq_is_q(n,p,q):
    # |sum_j n_j zeta_p^j|^2=q iff R0-q=R1=...=R_{p-1}
    r=[sum(n[j]*n[(j+t)%p] for j in range(p)) for t in range(p)]
    return all(r[t]==r[0]-q for t in range(1,p))

class Field:
    def __init__(self,p,m):
        self.p,self.m,self.q=p,m,p**m
        self.poly=irred(p,m)
        q=self.q
        self.d=[digits(x,p,m) for x in range(q)]
        self.base=[p**j for j in range(m)]
        self.add=[bytearray(q) for _ in range(q)]
        self.sub=[bytearray(q) for _ in range(q)]
        self.mul=[bytearray(q) for _ in range(q)]
        for a in range(q):
            for b in range(q):
                self.add[a][b]=sum((u+v)%p*k for u,v,k in zip(self.d[a],self.d[b],self.base))
                self.sub[a][b]=sum((u-v)%p*k for u,v,k in zip(self.d[a],self.d[b],self.base))
                product=[0]*(2*m-1)
                for i,ai in enumerate(self.d[a]):
                    if ai:
                        for j,bj in enumerate(self.d[b]):product[i+j]=(product[i+j]+ai*bj)%p
                for i in range(2*m-2,m-1,-1):
                    if product[i]:
                        c=product[i]
                        for j in range(m):product[i-m+j]=(product[i-m+j]-c*self.poly[j])%p
                self.mul[a][b]=sum(product[i]*self.base[i] for i in range(m))
        self.exp=[1]*(q-1); self.log=[-1]*q
        pf=[d for d in range(2,q) if (q-1)%d==0 and all(d%k for k in range(2,int(d**.5)+1))]
        def pow_slow(a,n):
            o=1
            while n:
                if n&1:o=self.mul[o][a]
                a=self.mul[a][a];n//=2
            return o
        self.primitive=next((a for a in range(1,q)
                             if all(pow_slow(a,(q-1)//d)!=1 for d in pf)),1)
        a=1
        for i in range(q-1):
            assert self.log[a]<0
            self.exp[i]=a;self.log[a]=i;a=self.mul[a][self.primitive]
        assert a==1
        self.trace=[]
        for a in range(q):
            s=0;v=a
            for _ in range(m):s=self.add[s][v];v=pow_slow(v,p)
            assert s<p
            self.trace.append(s)
        self.pair=[bytearray(self.trace[self.mul[b][x]] for x in range(q)) for b in range(q)]
    def power(self,a,n):return self.exp[self.log[a]*n%(self.q-1)] if a else 0
    def evaluate(self,terms):
        vals=[]
        for x in range(self.q):
            a=0
            for coef,degree in terms:
                z=self.power(x,degree)
                # coefficients are prime-subfield ints
                for _ in range(coef):a=self.add[a][z]
            vals.append(a)
        return vals
    def planar(self,F):
        q=self.q
        for a in range(1,q):
            seen=bytearray(q)
            ar=self.add[a]
            for x in range(q):
                v=self.sub[F[ar[x]]][F[x]]
                if seen[v]:return False,(a,x,v)
                seen[v]=1
        return True,None
    def check_component(self,F,c,gauss):
        p,q=self.p,self.q
        f=[self.trace[self.mul[c][F[x]]] for x in range(q)]
        sign=None
        sign_count=0
        # b=0 calibrates sign; compare rotations of one shared gauss^m
        for beta in range(q):
            counts=[0]*p
            row=self.pair[beta]
            for x in range(q):counts[(f[x]+row[x])%p]+=1
            assert abs_sq_is_q(counts,p,q),('not bent',c,beta,counts)
            candidates=[]
            for e in (1,-1):
                for j in range(p):
                    template=[e*gauss[(k-j)%p] for k in range(p)]
                    if equiv(counts,template):candidates.append((e,j))
            assert len(candidates)==1,('not weakly regular',c,beta,counts,candidates)
            if sign is None:sign=candidates[0][0]
            assert sign==candidates[0][0],('sign varies',c,beta,sign,candidates)
            sign_count+=1
        return sign,sign_count

def main():
    print(json.dumps({'method':'exact cyclotomic count vectors, Phi_p divisibility, exact derivative permutations',
                      'scope':'published small-field representatives ONLY; all nonzero c and all Walsh b',
                      'classes':sum(len(v) for v in CLASSES.values())},sort_keys=True))
    grand_classes=grand_components=grand_walsh=0
    for p,m in CLASSES:
        E=Field(p,m)
        g=base_gauss(p,m)
        for name,terms in CLASSES[(p,m)]:
            F=E.evaluate(terms)
            good,details=E.planar(F)
            assert good,(p,m,name,details)
            signs=Counter();walsh=0
            for c in range(1,E.q):
                sign,n=E.check_component(F,c,g)
                signs[sign]+=1;walsh+=n
            grand_classes+=1;grand_components+=E.q-1;grand_walsh+=walsh
            print(json.dumps({'class':name,'p':p,'m':m,'q':E.q,'modulus':E.poly,
                              'primitive':E.primitive,'terms':terms,'planar_derivatives':E.q-1,
                              'components':E.q-1,'walsh_coefficients':walsh,
                              'epsilon_counts':{str(k):signs[k] for k in sorted(signs)},
                              'result':'PASS'},sort_keys=True))
    print(json.dumps({'summary':'PASS','classes':grand_classes,'components':grand_components,
                      'walsh_coefficients':grand_walsh,'failures':0},sort_keys=True))
if __name__=='__main__':main()
