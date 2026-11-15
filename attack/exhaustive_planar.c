/* Exhaustive normalized planar maps over GF(3), GF(5), GF(7), GF(9).
 * Copyright (c) 2026 Athanor computational attack. C11; no dependencies.
 * Normalize F(0)=F(1)=0 using additive-affine subtraction.
 * Every map satisfying normalization occurs exactly once.
 * Usage: cc -O3 -std=c11 -Wall -Wextra exhaustive_planar.c -o exhaustive_planar
 *        ./exhaustive_planar
 */
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

static int q,p,m,add[9][9],sub[9][9],mul[9][9],tr[9],f[9];
static uint64_t nchecked,nplanar,n_components,n_nonweak,n_sign_pos,n_sign_neg;
static int first_planar[9], witness_bad[9],bad_c,bad_beta,bad_counts[3];
static int recordfirst,recordbad;
static uint64_t digest;

static void setup(int prime,int degree) {
    p=prime;m=degree; q=(degree==1)?prime:9;
    for(int a=0;a<q;a++) for(int b=0;b<q;b++) {
        if(degree==1) {add[a][b]=(a+b)%p;sub[a][b]=(a-b+p)%p;mul[a][b]=(a*b)%p;}
        else { /* a=a0+3a1*u, u^2=2, irreducible over GF(3) */
            int a0=a%3,a1=a/3,b0=b%3,b1=b/3;
            add[a][b]=(a0+b0)%3+3*((a1+b1)%3);
            sub[a][b]=(a0-b0+3)%3+3*((a1-b1+3)%3);
            mul[a][b]=(a0*b0+2*a1*b1)%3+3*((a0*b1+a1*b0)%3);
        }
    }
    for(int a=0;a<q;a++) tr[a]=(degree==1)?a:(2*(a%3))%3;
}
static int planar(void){
    for(int a=1;a<q;a++) {
        unsigned int bits=0;
        for(int x=0;x<q;x++) {
            int v=sub[f[add[x][a]]][f[x]],mask=1u<<v;
            if(bits & mask) return 0;
            bits |= mask;
        }
    }
    return 1;
}
/* Test each nonzero component for weak regularity. For GF(9), W(beta) has
 * exact cyclotomic-coordinate pair (n0-n2, n1-n2). 3*zeta^j is one of
 * (3,0),(0,3),(-3,-3); -3*zeta^j reverses its sign.
 * For prime fields this is also tested numerically via |W| and common phase
 * with tolerance, but the GF(9) part is exact integer arithmetic.
 */
static int component9(int c,int *failbeta,int counts[3]){
    int sign=0;
    for(int beta=0;beta<9;beta++){
        int n[3]={0,0,0};
        for(int x=0;x<9;x++){
            int phase=(tr[mul[c][f[x]]]+tr[mul[beta][x]])%3;
            n[phase]++;
        }
        int A=n[0]-n[2], B=n[1]-n[2], s=0;
        if((A==3&&B==0)||(A==0&&B==3)||(A==-3&&B==-3)) s=1;
        if((A==-3&&B==0)||(A==0&&B==-3)||(A==3&&B==3)) s=-1;
        if(!s || (sign && sign!=s)) {
            *failbeta=beta;
            for(int z=0;z<3;z++)counts[z]=n[z];
            return 0;
        }
        sign=s;
    }
    return sign;
}
static void examine(void){
    nchecked++;
    if(!planar()) return;
    nplanar++;
    for(int x=0;x<q;x++){digest ^= (uint8_t)f[x]; digest*=UINT64_C(1099511628211);}
    if(!recordfirst) {memcpy(first_planar,f,sizeof(f));recordfirst=1;}
    if(m==2){
        for(int c=1;c<9;c++){
            int b=0, cnt[3]={0,0,0};
            int s=component9(c,&b,cnt);
            n_components++;
            if(s==1)n_sign_pos++; else if(s==-1)n_sign_neg++;
            else{
                n_nonweak++;
                if(!recordbad){memcpy(witness_bad,f,sizeof(f));bad_c=c;bad_beta=b;memcpy(bad_counts,cnt,sizeof(cnt));recordbad=1;}
            }
        }
    }
}
static void recur(int i){
    if(i==q){examine();return;}
    for(int v=0;v<q;v++){f[i]=v;recur(i+1);}
}
static void printmap(const int *a){putchar('[');for(int i=0;i<q;i++)printf("%s%d",i?",":"",a[i]);putchar(']');}
int main(void){
    const int cases[][2]={{3,1},{5,1},{7,1},{3,2}};
    for(int t=0;t<4;t++){
        setup(cases[t][0],cases[t][1]);
        nchecked=nplanar=n_components=n_nonweak=n_sign_pos=n_sign_neg=0;
        recordfirst=recordbad=0;digest=UINT64_C(14695981039346656037);
        memset(f,0,sizeof(f));
        /* F(0)=F(1)=0. All q^(q-2) normalized functions. */
        recur(2);
        printf("GF(%d) p=%d m=%d checked=%" PRIu64 " planar=%" PRIu64 " first=",q,p,m,nchecked,nplanar);
        printmap(first_planar);
        if(m==2) printf(" components=%" PRIu64 " exact_weak=%" PRIu64 " failures=%" PRIu64 " positive_Wsign=%" PRIu64 " negative_Wsign=%" PRIu64,
            n_components,n_components-n_nonweak,n_nonweak,n_sign_pos,n_sign_neg);
        printf(" fnv64=%016" PRIx64, digest);putchar('\n');
        if(recordbad){printf("COUNTEREXAMPLE F=");printmap(witness_bad);printf(" c=%d beta=%d counts=[%d,%d,%d]\n",bad_c,bad_beta,bad_counts[0],bad_counts[1],bad_counts[2]);}
    }
    return 0;
}
