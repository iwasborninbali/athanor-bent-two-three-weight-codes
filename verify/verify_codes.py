#!/usr/bin/env python3
"""Independent reproducible checks for Tang et al., arXiv:1507.06148v3.
No CAS: irreducible polynomial field, direct word enumeration, cyclic FFT,
additive FFT for bentness, and explicit projective-support containment.
Python 3 with numpy. Output is deterministic JSON lines.
"""
import argparse, collections, hashlib, itertools, json, math, sys
import numpy as np


def irreducible(p, m, a):
    for degree in range(1, m // 2 + 1):
        for z in itertools.product(range(p), repeat=degree):
            d = list(z) + [1]
            r = list(a)
            for k in range(m - degree, -1, -1):
                v = r[k + degree]
                if v:
                    for j in range(degree + 1):
                        r[k + j] = (r[k + j] - v * d[j]) % p
            if not any(r[:degree]):
                return False
    return True


def factors(n):
    out = []
    k = 2
    while k * k <= n:
        if n % k == 0:
            out.append(k)
            while n % k == 0:
                n //= k
        k += 1
    if n > 1:
        out.append(n)
    return out


class Field:
    def __init__(self, p, m):
        self.p, self.m, self.q = p, m, p**m
        self.poly = next(list(z) + [1] for z in itertools.product(range(p), repeat=m)
                         if z[0] and irreducible(p, m, list(z) + [1]))
        self.digits = np.array([[x // (p**j) % p for j in range(m)]
                                for x in range(self.q)], dtype=np.int16)
        def rawmul(a, b):
            x = [int(v) for v in self.digits[a]]
            y = [int(v) for v in self.digits[b]]
            r = [0] * (2*m-1)
            for i in range(m):
                for j in range(m):
                    r[i+j] += x[i]*y[j]
            for k in range(2*m-2, m-1, -1):
                t = r[k] % p
                for j in range(m):
                    r[k-m+j] -= t * self.poly[j]
            return sum((r[j] % p) * p**j for j in range(m))
        def rawpow(a, n):
            r = 1
            while n:
                if n&1:
                    r = rawmul(r,a)
                a = rawmul(a,a)
                n >>= 1
            return r
        self.g = next(g for g in range(2, self.q)
                      if all(rawpow(g, (self.q-1)//l)!=1 for l in factors(self.q-1)))
        self.exp = np.empty(self.q-1,dtype=np.int32)
        self.log = np.full(self.q,-1,dtype=np.int32)
        v = 1
        for i in range(self.q-1):
            self.exp[i] = v
            self.log[v] = i
            v = rawmul(v,self.g)
        assert v == 1 and np.all(self.log[1:] >= 0)
        self.trace = np.zeros(self.q,dtype=np.int16)
        for x in range(1,self.q):
            s,y = 0,x
            for _ in range(m):
                s = self.add(s,y)
                y = self.power(y,p)
            assert s < p, ('trace not in prime subfield',p,m,x,s)
            self.trace[x] = s
        self.tshift = self.trace[self.exp]

    def add(self,x,y):
        return sum(((int(self.digits[x,j])+int(self.digits[y,j]))%self.p)*self.p**j
                   for j in range(self.m))

    def power(self,x,n):
        if n==0: return 1
        if x==0: return 0
        return int(self.exp[(int(self.log[x])*n)%(self.q-1)])

    def mul(self,x,y):
        if x==0 or y==0: return 0
        return int(self.exp[(int(self.log[x])+int(self.log[y]))%(self.q-1)])

    def fvalues(self, a, d):
        val = np.zeros(self.q,dtype=np.int16)
        la = int(self.log[a])
        val[self.exp] = self.trace[self.exp[(la+d*np.arange(self.q-1,dtype=np.int64))%(self.q-1)]]
        return val


def sign_from_counts(p,m,f):
    ct = np.bincount(f,minlength=p)
    chi = (-1 if p%4 == 3 else 1)
    if m%2==0:
        unit = (p-1)*chi**(m//2)*p**((m-2)//2)
        n = int(ct[0]-p**(m-1))
    else:
        unit = chi**((m-1)//2)*p**((m-1)//2)
        n = int(ct[1]-p**(m-1))
    assert unit and n%unit==0 and abs(n//unit)==1,(p,m,ct.tolist(),unit)
    return n//unit


def expected(p,m,eps,kind):
    punct = kind.endswith('_punct')
    kind = kind.replace('_punct','')
    u = p**(m-2)
    if m%2==0:
        k = eps*(-1 if p%4==3 else 1)**(m//2)*p**((m-2)//2)
        if kind=='zero':
            n = p**(m-1)-1+(p-1)*k
            out = {0:1,(p-1)*u:n,(p-1)*(u+k):(p-1)*(p**(m-1)-k)}
        else:
            n=(p-1)*(p**(m-1)-k)//2
            out={0:1,(p-1)**2*u//2:(p+1)*p**(m-1)//2+(p-1)*k//2-1,
                 (p-1)*((p-1)*u-2*k)//2:(p-1)*(p**(m-1)-k)//2}
    else:
        ch = -1 if p%4==3 else 1
        S = eps*ch**((m-1)//2)*p**((m-1)//2)
        T = eps*ch**((m-3)//2)*p**((m-3)//2)
        ps = ch*p
        if kind=='zero':
            n=p**(m-1)-1
            h=p**((m-3)//2)
            out={0:1,(p-1)*u:p**(m-1)-1,
                 (p-1)*(u-h):(p-1)*(p**(m-1)+p*h)//2,
                 (p-1)*(u+h):(p-1)*(p**(m-1)-p*h)//2}
        elif kind=='nsq':
            n=(p-1)*(p**(m-1)-S)//2
            out={0:1,(p-1)**2*u//2:p**(m-1)-1,
                 (p-1)*((p-1)*u+(1-ps)*T)//2:(p-1)*(p**(m-1)+ch*S)//2,
                 (p-1)*((p-1)*u-(1+ps)*T)//2:(p-1)*(p**(m-1)-ch*S)//2}
        else:
            n=(p-1)*(p**(m-1)+S)//2
            out={0:1,(p-1)**2*u//2:p**(m-1)-1,
                 (p-1)*((p-1)*u+(1+ps)*T)//2:(p-1)*(p**(m-1)+ch*S)//2,
                 (p-1)*((p-1)*u+(ps-1)*T)//2:(p-1)*(p**(m-1)-ch*S)//2}
    if punct:
        assert n%(p-1)==0
        n //= p-1
        out={w//(p-1):a for w,a in out.items()}
    return n,dict(sorted(out.items()))


def representatives(g, D):
    return np.array([x for x in D if x==min(g.mul(x,a) for a in range(1,g.p))],dtype=np.int32)


def distributions(g,D,direct=True):
    L=g.q-1
    mask=np.zeros(L,dtype=np.float64)
    mask[g.log[D]]=1.0
    t=(g.tshift==0).astype(np.float64)
    raw=np.fft.ifft(np.fft.fft(t)*np.conj(np.fft.fft(mask)))
    residual=float(max(np.max(np.abs(raw.imag)),np.max(np.abs(raw.real-np.rint(raw.real)))))
    assert residual<1e-6,(g.p,g.m,residual)
    weightsfft=(len(D)-np.rint(raw.real).astype(int))
    if direct:
        ix=g.log[D].astype(np.int64)
        offsets=(range(L) if L<=15625 else sorted(set([0,1,L//2,L-1]+list(np.linspace(0,L-1,257,dtype=int)))))
        for i in offsets:
            wd=len(D)-np.count_nonzero(g.tshift[(ix+i)%L]==0)
            assert wd==weightsfft[i],('DIRECT_FFT_MISMATCH',g.p,g.m,len(D),i,wd,weightsfft[i])
    weights=np.r_[0,weightsfft]
    return dict(sorted((int(w),int(c)) for w,c in zip(*np.unique(weights,return_counts=True)))),residual


def rank_from_distribution(dist,p,m):
    zero_messages=dist.get(0,0)
    assert p**m%zero_messages==0
    return m-round(math.log(zero_messages,p))


def minimality(g,D,maxq=3125):
    if g.q>maxq: return None
    beta = representatives(g,g.exp)
    indices=g.log[D].astype(np.int64)
    masks=[]
    for b in beta:
        t=(g.tshift[(int(g.log[b])+indices)%(g.q-1)]!=0)
        mask=0
        for j in np.flatnonzero(t): mask |= 1<<int(j)
        masks.append(mask)
    if any(v==0 for v in masks): return {'minimal':False,'reason':'nonzero message maps to zero'}
    for i,a in enumerate(masks):
        for j,b in enumerate(masks):
            if i!=j and a & ~b == 0:
                return {'minimal':False,'projective_beta_pair':[int(beta[i]),int(beta[j])],
                        'weights':[a.bit_count(),b.bit_count()]}
    return {'minimal':True,'projective_codewords':len(masks)}


def verify_field(g, a, d, kinds, do_minimal, extra_d=None):
    p,m,q=g.p,g.m,g.q
    f=g.fvalues(a,d)
    if extra_d is not None: f=(f+g.fvalues(1,extra_d))%p
    eps=sign_from_counts(p,m,f)
    z=np.exp((2j*np.pi/p)*f).reshape([p]*m)
    spectrum=np.fft.fftn(z)
    bentres=float(np.max(np.abs(np.abs(spectrum)-math.sqrt(q))))
    assert bentres<1e-6,(p,m,a,d,bentres)
    ct=np.bincount(f,minlength=p)
    assert int(ct.sum())==q
    sq={x*x%p for x in range(1,p)}
    D0=np.flatnonzero(f==0); D0=D0[D0!=0]
    Ds=np.flatnonzero(np.isin(f,list(sq)))
    Dn=np.flatnonzero((f!=0)&~np.isin(f,list(sq)))
    choices={'zero':D0,'sq':Ds,'nsq':Dn}
    for kind in kinds:
        base=kind.replace('_punct','')
        D=choices[base]
        if kind.endswith('punct'): D=representatives(g,D)
        actual,fftres=distributions(g,D)
        n,pred=expected(p,m,eps,kind)
        discrepancies=[]
        if n!=len(D): discrepancies.append('length')
        if pred!=actual: discrepancies.append('weight_distribution')
        row={'type':'code','p':p,'m':m,'q':q,'poly':g.poly,'primitive':g.g,
             'a':a,'exponent':d,'extra_exponent':extra_d,'epsilon':eps,'kind':kind,'length':len(D),
             'dimension':rank_from_distribution(actual,p,m),
             'weight_distribution':actual,'table_prediction':pred,
             'table_length':n,'status':'CONFIRMED' if not discrepancies else 'REFUTED',
             'differences':discrepancies,'fft_rounding_residual':round(fftres,9),
             'walsh_bent_residual':round(bentres,9),
             'direct_checks':(q-1 if q-1<=15625 else len(set([0,1,(q-1)//2,q-2]+list(np.linspace(0,q-2,257,dtype=int)))))}
        if do_minimal: row['minimality']=minimality(g,D)
        print(json.dumps(row,sort_keys=True),flush=True)
    return eps,ct


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--quick',action='store_true',help='reduced cases to smoke test')
    ap.add_argument('--no-minimality',action='store_true')
    args=ap.parse_args()
    params=[(3,m) for m in range(2,11)]+[(5,m) for m in range(2,8)]+[(7,m) for m in (2,3,4,5)]+[(11,m) for m in (2,3,4)]+[(13,m) for m in (2,3,4)]+[(17,4),(19,3),(23,3)]
    if args.quick:params=[(3,m) for m in (2,3,4,5,6)]+[(5,3),(5,5)]
    families=['zero','sq','nsq','zero_punct','sq_punct','nsq_punct']
    print(json.dumps({'type':'metadata','version':'1','numpy':np.__version__,'params':params,
                      'methods':['direct log-index codeword count','multiplicative cyclic FFT',
                                 'additive finite-field Walsh FFT','projective support containment'],
                      'arxiv':'1507.06148v3'}),flush=True)
    for p,m in params:
        g=Field(p,m)
        for a in (1,g.g):
            eps,ct=verify_field(g,a,2,families,not args.no_minimality)
            print(json.dumps({'type':'value_counts','p':p,'m':m,'a':a,'d':2,
                              'epsilon':eps,'counts':ct.tolist()},sort_keys=True),flush=True)
        if (p,m) in ((3,6),(5,6)) and not args.quick:
            cases = [(1,10),(g.power(g.g,7),210)] if p==3 else [(1,26),(g.g,26)]
            for a,d in cases:
                eps,ct=verify_field(g,a,d,families,False)
                print(json.dumps({'type':'value_counts','p':p,'m':m,'a':a,'d':d,
                                  'epsilon':eps,'counts':ct.tolist()},sort_keys=True),flush=True)
        if not args.quick and p==3 and m in (4,5,6,7,8):
            i=5 if m==6 else 3
            if math.gcd(i,m)==1:
                d=(3**i+1)//2
                for a in (1,g.g):
                    eps,ct=verify_field(g,a,d,families,False)
                    assert eps==((-1)**(m-1) if a==1 else (-1)**m),('Coulter-Matthews sign',m,a,eps)
        if not args.quick and m==4 and p in (3,5,7):
            d=p**3+p**2-p+1
            eps,ct=verify_field(g,1,d,families,False,extra_d=2)
            assert eps==-1,('Helleseth-Kholosha binomial sign',p,eps)
    print(json.dumps({'type':'done','status':'complete'}),flush=True)

if __name__=='__main__':main()
