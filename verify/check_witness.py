#!/usr/bin/env python3
"""Small-field counterexamples verified without numpy or verify_codes.py.
F_9 = F_3[t]/(t^2+1); g=1+t, f_a(x)=Tr(a*x^2).
All additions, products, traces, weights and dimensions use plain tuples.
"""
import collections,json
p=3
zero=(0,0);one=(1,0);g=(1,1)
def add(a,b):return ((a[0]+b[0])%3,(a[1]+b[1])%3)
def mul(a,b):return ((a[0]*b[0]-a[1]*b[1])%3,(a[0]*b[1]+a[1]*b[0])%3)
def trace(a):return (2*a[0])%3
E=[(a,b) for a in range(3) for b in range(3)]
assert len(set(E))==9

def code(a):
    f={x:trace(mul(a,mul(x,x))) for x in E}
    count=collections.Counter(f.values())
    eps=(count[0]-3)//(-2)
    assert eps in (1,-1)
    D=[x for x in E if x!=zero and f[x]==0]
    weights=[sum(trace(mul(b,x))!=0 for x in D) for b in E]
    dist=collections.Counter(weights)
    dim=2-(0 if dist[0]==1 else 1 if dist[0]==3 else 2)
    printed_length=3-1+eps*(3-1)
    corrected_length=3-1+eps*(-1)*(3-1)
    assert corrected_length==len(D)
    return {'coefficient':a,'epsilon':eps,'f_value_counts':dict(sorted(count.items())),
            'D':D,'actual_n':len(D),'actual_dimension':dim,'weight_distribution_of_messages':dict(sorted(dist.items())),
            'theorem_3_1_printed_n':printed_length,'theorem_3_1_corrected_n':corrected_length,
            'corollary_3_3_printed_n':(3-1)//(3-1)+eps,
            'corollary_3_3_corrected_n':len(D)//2}

for a in (one,g):
    print(json.dumps(code(a),sort_keys=True))
assert code(one)['actual_n']==4 and code(one)['theorem_3_1_printed_n']==0
assert code(g)['actual_n']==0 and code(g)['actual_dimension']==0
print(json.dumps({'type':'done','status':'COMPLETE'}))
