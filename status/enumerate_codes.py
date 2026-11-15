#!/usr/bin/env python3
"""Enumerate Tang et al. §III/§IV trace codes over GF(3^3), GF(3^4).
No dependencies beyond the Python standard library. Deterministic output.
Output witness.json with every word, support-subset witnesses and dimensions.
"""
from collections import Counter
from itertools import product
import json
from pathlib import Path

P = 3
POLYS = {3: (1, 2, 0, 1), 4: (2, 1, 0, 0, 1)}  # monic irreducible

class Field:
    def __init__(self, m):
        self.m = m
        self.mod = POLYS[m]
        self.q = P ** m
        self.elems = list(product(range(P), repeat=m))
        self.zero = (0,) * m
        self.one = (1,) + (0,) * (m-1)

    def add(self, x, y):
        return tuple((a+b) % P for a,b in zip(x,y))

    def multiply(self, x, y):
        m = self.m
        accum = [0] * (2*m-1)
        for i,a in enumerate(x):
            for j,b in enumerate(y):
                accum[i+j] = (accum[i+j] + a*b) % P
        for k in range(2*m-2, m-1, -1):
            c = accum[k]
            for j in range(m):
                accum[k-m+j] = (accum[k-m+j] - c*self.mod[j]) % P
            accum[k] = 0
        return tuple(accum[:m])

    def power(self, x, n):
        out = self.one
        while n:
            if n & 1:
                out = self.multiply(out,x)
            x = self.multiply(x,x)
            n >>= 1
        return out

    def trace(self, x):
        accum = self.zero
        y = x
        for _ in range(self.m):
            accum = self.add(accum,y)
            y = self.power(y,3)
        assert accum[1:] == (0,)*(self.m-1), ('trace in base field', x, accum)
        return accum[0]


def support(a):
    return {i for i,z in enumerate(a) if z != 0}

def scalar_multiple(a, b):
    """Vectors proportional over F_3, including zero where applicable."""
    return any(all(x == t*y % 3 for x,y in zip(a,b)) for t in (1,2))

def minimality(words):
    for i,c in enumerate(words):
        if not any(c):
            continue
        sa = support(c)
        for j,d in enumerate(words):
            if i == j or not any(d):
                continue
            sb = support(d)
            if sb.issubset(sa) and not scalar_multiple(c,d):
                return {"minimal": False, "containing_index": i, "contained_index": j,
                        "containing_word": c, "contained_word": d}
    return {"minimal": True}


def bent_test(field):
    q=field.q
    f={x:field.trace(field.multiply(x,x)) for x in field.elems}
    for beta in field.elems:
        counts=[0,0,0]
        for x in field.elems:
            counts[(f[x]+field.trace(field.multiply(beta,x)))%3] += 1
        a,b,c=counts
        norm2=a*a+b*b+c*c-a*b-b*c-a*c
        assert norm2==q, (field.m,beta,counts,norm2)
    return f


def main():
    result={"method":"explicit polynomial-basis finite fields and full trace-code enumeration",
            "cases":[]}
    for m in (3,4):
        K=Field(m)
        f=bent_test(K)
        for lev in ((0,1,2) if m==3 else (0,)):
            D=[x for x in K.elems if x!=K.zero and f[x]==lev]
            wordlist=[ [K.trace(K.multiply(beta,d)) for d in D] for beta in K.elems]
            assert len({tuple(w) for w in wordlist}) == K.q
            counts=dict(sorted(Counter(sum(t != 0 for t in w) for w in wordlist).items()))
            if (m,lev)==(3,0):
                assert counts=={0:1,4:12,6:8,8:6}, counts  # Table II
            if (m,lev)==(4,0):
                assert counts=={0:1,12:60,18:20}, counts  # Table I
            basis=[K.elems.index(tuple(1 if j==i else 0 for j in range(m))) for i in range(m)]
            entry={"m":m,"p":P,"level":lev,"modulus_coefficients_low_to_high":list(K.mod),
                   "D":D,"beta_order":K.elems,"words":wordlist,
                   "basis_word_indices":basis,
                   "length":len(D),"dimension":m,
                   "weight_multiplicities":{str(k):v for k,v in counts.items()},
                   "minimality":minimality(wordlist)}
            result['cases'].append(entry)
            print(f"GF(3^{m}) F(x)=x^2 f(x)=Tr(x^2) level={lev} "
                  f"n={len(D)} dimension={m} histogram={counts} "
                  f"minimal={entry['minimality']['minimal']}")
            if not entry['minimality']['minimal']:
                z=entry['minimality']
                print(f"  nonminimal witness beta indices {z['containing_index']} "
                      f"and {z['contained_index']}; containing={z['containing_word']} "
                      f"contained={z['contained_word']}")
    dst=Path(__file__).with_name('witness.json')
    dst.write_text(json.dumps(result,indent=2,separators=(',',': '),sort_keys=True)+'\n',encoding='utf-8')
    print(f"Wrote {dst.name}; all Walsh squared magnitudes checked exactly; no floating point.")

if __name__=='__main__':
    main()
