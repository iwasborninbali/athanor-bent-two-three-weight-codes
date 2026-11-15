#!/usr/bin/env python3
"""Exact small-degree replay of historical Helleseth--Kholosha
ternary monomial conjecture, equations (6)-(7), Tang et al. pp.14-15.
Checks bentness/one common Walsh phase/sign and Table I weights, NOT
an independent proof of the full rational dual trace expression (7).
Usage: python3 helleseth_kholosha.py
"""
from coulter_matthews import Field,verify_walsh,code_distribution,expected_distribution,check_minimal

def main():
    for k in (1,3):
        m=2*k; E=Field(m);q=E.q
        exponent=(q-1)//4+3**k+1
        c=E.power(E.generator,(3**k+1)//4)
        f=[E.trace[E.mul(c,E.power(x,exponent))] for x in range(q)]
        eps,beta,pair=verify_walsh(E,f)
        assert eps==(-1)**(k+1),(k,m,eps,beta,pair)
        # Paper predicts W(0)=-3^k; verify exact Gaussian period coordinates.
        n=[sum(t==r for t in f) for r in range(3)]
        W0=(n[0]-n[2],n[1]-n[2])
        assert W0==(-3**k,0),(m,W0)
        D,weights=code_distribution(E,f)
        if m>=4:
            assert weights==expected_distribution(m,eps),(m,weights)
        if D:
            minimum,witness,projective=check_minimal(E,D)
        else:
            minimum,witness,projective="N/A (zero-length code)",None,0
        print(f'GF(3^{m}) q={q} poly={E.poly} primitive={E.generator} '
              f'HK_k={k} exponent={exponent} c={c} W0_pair={W0} '
              f'epsilon={eps} global_Walsh=PASS '+
              f'Theorem_3_1_printed_length={3**(m-1)-1+eps*2*3**((m-2)//2)}')
        print(f'  defining_set_length={len(D)} weight_histogram={sorted(weights.items())} '
              f'Table_I={"PASS" if m>=4 else "DEGENERATE_m2"} '+
              f'minimal={minimum} projective_classes={projective} '
              f'containment_witness={witness}')
if __name__=='__main__':main()
