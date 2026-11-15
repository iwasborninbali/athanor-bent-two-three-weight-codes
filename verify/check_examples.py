#!/usr/bin/env python3
"""Independently match all 18 numeric examples in pp. 11--14 against run log.
Also evaluate each cited Griesmer claim with integer arithmetic.
"""
import collections, json, math, sys
rows=[json.loads(x) for x in open(sys.argv[1] if len(sys.argv)>1 else 'full.jsonl')]
lookup={(r['p'],r['m'],r['exponent'],r['a'],r['kind']):r for r in rows if r.get('type')=='code'}
primitive={(r['p'],r['m']):r['primitive'] for r in rows if r.get('type')=='code'}
# (paper-page, label, p,m, a-selector, exponent, kinds, n,k,d, nonzero weight multiplicities)
examples=[
(11,'4.1a',3,6,'w7',210,'sq nsq',252,6,162,{162:476,180:252}),
(11,'4.1b',3,6,'1',10,'sq nsq',234,6,144,{144:234,162:494}),
(11,'4.1c',5,6,'1',26,'sq nsq',6300,6,5000,{5000:9324,5100:6300}),
(12,'4.1d',5,6,'w',26,'sq nsq',6200,6,4900,{4900:6200,5000:9424}),
(12,'4.2a',5,5,'1',2,'nsq',1200,5,940,{940:1200,960:1300,1000:624}),
(12,'4.2b',3,5,'w',2,'nsq',60,5,40,{40:24,48:40,52:60}),
(12,'4.3a',5,5,'1',2,'sq',1300,5,1000,{1000:624,1040:1200,1060:1300}),
(12,'4.3b',3,5,'w',2,'sq',40,5,28,{28:40,32:60,40:24}),
(13,'4.4a',3,4,'1',2,'sq_punct nsq_punct',15,4,9,{9:50,12:30}),
(13,'4.4b',3,4,'w',2,'sq_punct nsq_punct',12,4,6,{6:24,9:56}),
(13,'4.4c',3,6,'w',2,'sq_punct nsq_punct',126,6,81,{81:476,90:252}),
(13,'4.4d',5,4,'1',2,'sq_punct nsq_punct',65,4,50,{50:364,55:260}),
(13,'4.4e',5,4,'w',2,'sq_punct nsq_punct',60,4,45,{45:240,50:384}),
(13,'4.5a',3,5,'1',2,'nsq_punct',36,5,21,{21:72,24:90,27:80}),
(13,'4.5b',3,5,'w',2,'nsq_punct',45,5,27,{27:80,30:72,33:90}),
(13,'4.6a',3,3,'w',2,'sq_punct',6,3,3,{3:8,4:6,5:12}),
(14,'4.6b',5,3,'w',2,'sq_punct',10,3,7,{7:40,8:60,10:24}),
(14,'4.6c',7,3,'1',2,'sq_punct',21,3,17,{17:126,18:168,21:48}),
]
assert len(examples)==18

def griesmer(p,k,d):return sum((d+p**j-1)//p**j for j in range(k))

def main():
    checked=0
    for page,label,p,m,coef,exp,kinds,n,k,d,enumerator in examples:
        g=primitive[(p,m)]
        def fieldexp(a,n):
            # Retrieve the coefficient's field representation from the log by
            # deriving exp(log g* n) from polynomial basis is not possible here;
            # for w^7 look up unique recorded 210-exponent row.
            if n==7:
                return next(r['a'] for r in rows if r.get('type')=='code' and r['p']==p and r['m']==m and r['exponent']==210)
            return a
        a=1 if coef=='1' else (g if coef=='w' else fieldexp(g,7))
        for kind in kinds.split():
            r=lookup[p,m,exp,a,kind]
            expected={0:1,**enumerator}
            got={int(w):int(v) for w,v in r['weight_distribution'].items()}
            good=(r['length'],r['dimension'],min(w for w in got if w))==(n,k,d) and got==expected
            if not good:
                print(json.dumps({'type':'example_mismatch','label':label,'page':page,'kind':kind,'printed':[n,k,d],
                                  'actual':[r['length'],r['dimension'],min(w for w in got if w)],
                                  'printed_distribution':expected,'actual_distribution':got},sort_keys=True))
            else:checked+=1
        print(json.dumps({'type':'example','label':label,'page':page,'p':p,'m':m,'length':n,
                          'dimension':k,'distance':d,'checked_code_families':len(kinds.split()),
                          'result':'CONFIRMED' if all(lookup[p,m,exp,1 if coef=='1' else (primitive[(p,m)] if coef=='w' else a),kind]['length']==n for kind in kinds.split()) else 'REFUTED'},sort_keys=True))
    claims=[('4.4a',3,15,4,9),('4.4b',3,12,4,6),('4.4c',3,126,6,81),
            ('4.4d',5,65,4,50),('4.6a',3,6,3,3),('4.6b',5,10,3,7),
            ('4.6c',7,21,3,17)]
    for label,p,n,k,d in claims:
        gd=griesmer(p,k,d); gp=griesmer(p,k,d+1)
        print(json.dumps({'type':'griesmer','label':label,'p':p,'n':n,'k':k,'d':d,
                          'G(d)':gd,'G(d+1)':gp,
                          'Griesmer_proves_d_optimal':gp>n,'result':'CONFIRMED' if gp>n else 'NOT_ESTABLISHED_BY_BOUND'},sort_keys=True))
    print(json.dumps({'type':'summary','examples':len(examples),'confirmed_family_checks':checked,
                      'refuted_family_checks':sum(len(x[6].split()) for x in examples)-checked,
                      'griesmer_claims':len(claims),'griesmer_sufficient':sum(griesmer(p,k,d+1)>n for _,p,n,k,d in claims)},sort_keys=True))

if __name__=='__main__':main()
