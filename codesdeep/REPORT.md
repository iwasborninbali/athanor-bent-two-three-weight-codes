# codesdeep — consolidated errata and planar-component audit

**Date:** 2026-11-15. **Laboratory:** Athanor (Aleksei Kudriashov). **Target:** C. Tang, N. Li, Y. Qi, Z. Zhou, T. Helleseth, *Linear codes with two or three weights from weakly regular bent functions*, *IEEE Trans. Inf. Theory* **62** (2016), 1166–1176, DOI [10.1109/TIT.2016.2518678](https://doi.org/10.1109/TIT.2016.2518678), [arXiv:1507.06148v3](https://arxiv.org/abs/1507.06148). Page numbers refer to the **16 numbered pages of v3**, as in the upstream audits; they are not journal pagination.

**Evidence convention:** `PROVED` means a complete written argument appears here; `REPLAYED-EXACT` means exhaustive within the *stated finite set* with an included deterministic program, complete log, SHA-256, and separate checker; `PRIOR ART` denotes an identified external result; `BOUNDED-SEARCH` is **not** a proof of absence or novelty. No paid services.

## Question

1. Assemble **every issue identified so far by this laboratory** in Tang et al., distinguishing demonstrated incorrect mathematical assertions, printing discrepancies, and a proof/optimality justification gap. Determine whether each has been *verified* corrected between arXiv versions, rather than assuming a newer version repairs it.
2. Section II-B, p. 2: for a planar map `F:GF(p^m)->GF(p^m)` and every nonzero `c`, must `f_c(x)=Tr(c F(x))` be **weakly regular bent**, not merely bent? List and test **every currently documented small-field class** through `GF(3^5)` and `GF(5^3)`; decide what is and is not settled by published classification.

## Results with evidence status

### A. Errata: consolidated identified issues

| ID / location (arXiv v3) | Printed assertion / issue | Correct statement and smallest explicit witness | Evidence / revision status |
|---|---|---|---|
| **E1. Theorem 3.1 and Corollary 3.3**, p. 10 | For **even** `m`, displayed zero-set and punctured lengths omit `chi=((-1)/p)^(m/2)`. | `n=p^(m-1)-1 + epsilon*(p-1)*chi*p^((m-2)/2)`; `n_punct=(p^(m-1)-1)/(p-1) + epsilon*chi*p^((m-2)/2)`. In `GF(9)=F_3[t]/(t^2+1)`, `f=Tr(x^2)`, `epsilon=-1`, defining set `{4,5,7,8}` using base-3 index `a+3b`; **actual** `n=4`, punctured `n=2`, **printed** `n=0`, punctured `n=0`. | **PROVED + REPLAYED-EXACT** (`errata.py`, `check_errata.py`). **Still printed incorrectly in v3**, despite its own **Table I and Table III** already containing the correct sign. These tables are *not* wrong. Whether v1/v2 differed is **unverified**. |
| **E2. Dimension/nondegeneracy at `m=2`**, Theorem 3.1/Corollary 3.3 and square/nonsquare assertions in §IV (notably Theorem 4.1/Corollary 4.4), pp. 10–13 | Unqualified dimension `m` (and, for empty codes, nontrivial two-weight status). | In the same `GF(9)`, with `w=1+t` primitive, `f=Tr(w x^2)` has **empty** zero defining set and code `[0,0]`, rather than dimension 2. With `f=Tr(x^2)`, `D_square={3,6}`, `D_nonsquare={1,2}` each has `n=2`, `k=1` rather than `k=2` (and puncturing has length 1 and dimension 1). **Safe repair for these claims: exclude `m=2`**, e.g. require `m>=4` for the even-degree formulas; no claim is made that this is the most general necessary condition. | **REPLAYED-EXACT + separate tuple checker**; upstream `verify/full.jsonl` found 14 dimension/table discrepancies, all at `m=2`, across its sampled primes. **Still unqualified in v3**; earlier-version history **not established**. |
| **E3. Two examples after Theorems 4.2 and 4.3**, p. 12 | Both say `p=3,m=5` and `k=5`, but print respectively `1+24z^40+40z^48+60z^52` with `[60,5,40]`, and `1+40z^28+60z^32+24z^40` with `[40,5,28]`. Each enumerator has **125** rather than `3^5=243` total words. | For the **stated** `p=3,m=5`, `f=Tr(w x^2)` with `w` primitive: **nonsquare:** `[90,5,54]`, `1+80z^54+72z^60+90z^66`; **square:** `[72,5,42]`, `1+72z^42+90z^48+80z^54`. The printed enumerators instead exactly describe `p=5,m=3`, `k=3`: `[60,3,40]` and `[40,3,28]`. A likely transposition is **not** evidence of author intent. | **REPLAYED-EXACT**, two separate implementations of full codeword counts; **still incorrect in v3**. No explicit later correction located by upstream priority audit. |
| **E4. Corollary 4.4 Griesmer justification**, p. 13 | Declaring the ternary `[12,4,6]` example optimal *because of the Griesmer bound*. | `G_3(4,7)=ceil(7)+ceil(7/3)+ceil(7/9)+ceil(7/27)=7+3+1+1=12`. Thus the bound **does not exclude `[12,4,7]`**; a separate nonexistence argument is needed to deduce optimality. **Optimality itself is not refuted.** | **PROVED (one-line exact arithmetic)** + separate checker, **justification gap still in v3**; not a demonstrated false code parameter or false optimality claim. |

**This is the complete list of issues *found in `verify/`, `attack/`, `status/`, and `priority/`***, not a guarantee that no other mistakes exist in the paper. E1 is one missing-factor defect affecting two printed statements; E3 comprises two independently false examples. Editorial wording such as “any a planar” is not counted as a substantive mathematical error.

**Full argument for E1.** With `p*=((-1)/p)*p`, the paper's `(2)` and `f*(0)=0` give `W_f(0)=epsilon*(sqrt(p*))^m=epsilon*chi*p^(m/2)` for even `m`. This is a rational integer, so every Galois conjugate `sum_x zeta_p^(a f(x))`, `a` nonzero in `F_p`, equals it. Character orthogonality yields `p * #{x:f(x)=0}=p^m+(p-1)W_f(0)`. Remove `x=0`; homogeneity means the remaining set is a union of `F_p^*`-orbits of size `p-1`, giving the punctured expression. In `GF(9)`, `W_f(0)=+3` and `(sqrt(-3))^2=-3`, so `epsilon=-1`; the sign omitted in the printed formula is decisive.

**Version distinction: no verified v1→v2→v3 repair.** The upstream [`priority/REPORT.md`](../priority/REPORT.md) explicitly records failed attempts to retrieve v1 and v2, and direct retrieval of those versioned PDFs was again unsuccessful in this audit. What **is** verified is that v3 internally contains the *correct* Legendre factor in Tables I/III and earlier Walsh identities, while its theorem and corollary **still omit** it. Thus “already correct elsewhere in v3” **must not** be relabeled “fixed in a later arXiv version.” No item is safely classed as a demonstrated earlier-to-later repair. Comparison with an authenticated IEEE print/PDF remains outstanding.

### B. Planar weak regularity: complete computation for known representatives, not all maps

**Published scope and precise families.** Andreoli, Budaghyan, Coulter, Haukenes, Kaleyski, Piccione, *On a classification of planar functions in characteristic three*, *Cryptography and Communications* **17** (2025), 823–853, [DOI](https://doi.org/10.1007/s12095-025-00781-y), **Appendix A, Table 1 (printed p. 840, PDF page index 17)**, enumerates the **currently known** CCZ-inequivalent representatives through degree five:

| Field | Representatives independently computed; 2025 table IDs | No. |
|---|---|---:|
| `GF(3)` | `x^2` (prime field) | 1 |
| `GF(3^2)` | `2.1 x^2` (FF) | 1 |
| `GF(3^3)` | `3.1 x^2` (FF), `3.2 x^4` (Albert twisted field) | 2 |
| `GF(3^4)` | `4.1 x^2` (FF), `4.2 x^14` (Coulter–Matthews), `4.3 x^36+2x^10+2x^4` (Dickson / Budaghyan–Helleseth / Zhou–Pott overlap) | 3 |
| `GF(3^5)` | `5.1 x^2`; `5.2 x^4`; `5.3 x^10`; `5.4 x^14` (CM); `5.5 x^10+x^6+2x^2`; `5.6 x^10+2x^6+2x^2` (Ten-Six-Two / Ding–Yuan); `5.7 x^90+x^2` (At–Cohen–Weng sporadic); `5.8 x^162+x^108+2x^84+x^2` (Coulter–Kosick sporadic) | 8 |
| `GF(5)` | `x^2` | 1 |
| `GF(5^2)` | `x^2` (known FF class) | 1 |
| `GF(5^3)` | `x^2`, `x^6` (Albert twisted field; commutative semifields of order `p^3` classified by Menichetti, as summarized in the 2025 paper §3) | 2 |
| **Total** | **19 distinct tested representatives** | **19** |

The Zhou–Pott family has degree `2m` and can coincide with Dickson/Budaghyan–Helleseth at small parameters; the 2025 classification places the `GF(3^4)` representative in exactly this overlap. **Coulter–Matthews is the only known non-Dembowski–Ostrom planar family** according to the 2025 survey (§3). The 2025 paper **does not claim that Table 1 exhausts arbitrary planar maps**; rather, it classifies the *known* examples. Similarly, Knuth's `p^2` and Menichetti's `p^3` **commutative-semifield** classifications concern **quadratic/DO** planar maps, not hypothetical unknown nonquadratic planar maps. The planar **monomial** classification for `m<=4` also is *not* a classification of all maps.

**Exact run (`planar_classes.py` / `planar_classes.log`):** for each of those 19 polynomials over the explicitly logged quotient-field modulus, checked all `q-1` nonzero derivative permutations, all `q-1` nonzero `c`, and all `q` Walsh frequencies `b`. **19/19 planar, 2,514 scalar components, 522,990 Walsh coefficients, no failures**. Exactly half of components at each field/class have `epsilon=+1` versus `-1`, as shown per line in the log. Cyclotomic arithmetic is exact for both `p=3` **and `p=5`**: for counts `n_j` of phases `zeta_p^j`, require `R_0-q=R_1=...=R_(p-1)` where `R_t=sum_j n_j n_(j+t)`, precisely equivalent to squared Walsh norm `q`; compare every count vector modulo the all-ones vector with a cyclic rotation of `epsilon * G_p^m`, where `G_p=sum_{a=1}^{p-1}(a|p)zeta_p^a`. A *single* epsilon must work for every b. These are integer checks, not numerical rounding.

**Independent check (`check_planar.py` / log):** different tuple-based quotient-field arithmetic, independent multiplication/exponentiation, all derivatives for all 19, and direct complex Walsh sums for **`c=1` and every b**: **2,533 independently checked Walsh coefficients, all pass**, maximum phase residual logged. This is genuinely independent field arithmetic but **does not independently replay all 2,514 components**. The full coverage is the exact primary run; its broader structural interpretation uses the following complete argument and published CM theorem.

**PROVED: all DO-planar scalar components are weakly regular, in *every odd-characteristic field*.** A DO polynomial `F(x)=sum_{i,j} a_ij x^(p^i+p^j)` is quadratic over `F_p`, so for `c!=0` the function `Q_c(x)=Tr(c F(x))` is an `F_p`-quadratic form. Its polar bilinear form `B_c(x,y)=Q_c(x+y)-Q_c(x)-Q_c(y)` is nondegenerate: if `a!=0` were in the radical, the additive derivative `Q_c(x+a)-Q_c(x)` would be constant; however, because `F` is planar, `x -> F(x+a)-F(x)` is a bijection of `GF(p^m)`, and its nonzero trace scalar cannot be constant. By diagonalization over `F_p` (odd p), `Q_c=sum_i d_i X_i^2` with all `d_i!=0`, hence `W_{Q_c}(0)=prod_i (d_i|p) * G_p^m`, with one sign `epsilon_c` independent of b. For any linear Walsh term `Tr(bx)`, nondegeneracy gives a unique `v` with `B_c(v,x)=Tr(bx)`; completing the square yields `W_{Q_c}(b)=zeta_p^(-Q_c(v)) W_{Q_c}(0)`. This is exactly weak regularity. **No finite computational result is required for this theorem.**

**EA invariance (proof).** For `H(x)=A_out(F(A_in(x)))+L(x)` with invertible affine `F_p`-linear input/output maps and affine `L`, `Tr(cH(x))` is a nonzero scalar component of F after an invertible output dual map, composed with an affine input permutation and plus an affine scalar function. Changing variables in Walsh sums only permutes frequencies and multiplies by roots of unity. Therefore bentness and a common weak-regularity sign are invariant under EA equivalence. For planar DO functions, the cited 2025 §2.3 records that CCZ equivalence coincides with EA equivalence. **Coulter–Matthews components** were already proved weakly regular by Helleseth–Hollmann–Kholosha–Wang–Xiang, *IEEE TIT* **55** (2009), 5272–5283, [Theorem 1.2 / open PDF](https://bpb-us-w2.wpmucdn.com/sites.udel.edu/dist/a/6387/files/2017/08/pbentJune09-fin-1eyo5xc.pdf). Thus *all currently documented families* above have weakly regular components by prior theorems, independently reinforced by the computations; neither statement is a universal theorem about unknown planar maps.

**Complete-map exception for the smallest fields:** `exhaustive_prime.py` visits **all** `3^(3-2)=3` normalized maps over `GF(3)` and **all** `5^(5-2)=125` over `GF(5)`: respectively **2 and 4** planar maps, all **4 and 16** nonzero components weakly regular, every frequency exactly checked. Every map can be uniquely normalized by subtracting `F(0)+(F(1)-F(0))*x`, which preserves planarity and scalar-component weak regularity. This is a genuine all-map result for `GF(3)` and `GF(5)`, not just representatives. The upstream [`attack/exhaustive_planar.c`](../attack/exhaustive_planar.c) and its independent Python backtracker [`attack/check_gf9.py`](../attack/check_gf9.py) likewise exhaust **4,782,969** normalized `GF(9)` maps; **1,296** planar; **10,368** components pass (matching ordered digest `4bedb2f3ee993315`). That `GF(9)` full-map result is **prior repository work**, not newly rerun or claimed as a fresh discovery here.

**Answer to “settled for these fields by classification plus computation?”** **NO for `GF(3^3)`, `GF(3^4)`, `GF(3^5)`, `GF(5^2)`, `GF(5^3)` as a theorem about *every possible planar mapping*.** The named representatives exhaust *currently known classes*, not all unknown classes; all-DO/CM theorems prove their positive cases but do not exclude an unknown non-DO, non-CM function. For `GF(3)`, `GF(5)`, and (by prior exhaustive repository work) `GF(9)`, the question **is settled affirmatively within those entire fields** by literal all-map enumeration and affine normalization. A later vectorial bent example with non-weakly-regular components does *not* settle Tang's question unless it is also planar. No counterexample or universal proof was verified by our bounded 2026 search.

## Priority check

**Searched within repository first:** [`verify/REPORT.md`](../verify/REPORT.md) (signed lengths, 14 `m=2` dimension cases, both examples, Griesmer); [`attack/REPORT.md`](../attack/REPORT.md) (complete sign derivation, `GF(9)` exhaustive 1,296 planar-map enumeration, sparse `GF(27)/GF(81)` and CM); [`status/REPORT.md`](../status/REPORT.md) (earlier literature and false inference from vectorial bent); [`priority/REPORT.md`](../priority/REPORT.md) (targeted IEEE errata/comment/title/DOI citations, Springer/IEEE/MathOverflow/Math StackExchange/OpenAlex queries; documented access limits). Read also the upstream implementation and witnesses cited there. **Nothing in E1–E4 is newly discovered by codesdeep.**

**Independent public literature check, 2026-11-15:** searched target DOI/arXiv ID with “correction”, “erratum”, “v1”, “v2”, length/Griesmer/example phrases; searched “planar weakly regular components” with “counterexample”/“universal theorem”, small-field classifications `p^2`, `p^3`, `p^4`, Ding–Yuan, Zhou–Pott, TST/CK/ACW, and 2025–2026 publications. Read the open-access [Andreoli et al. 2025 article and Appendix A Table 1](https://link.springer.com/article/10.1007/s12095-025-00781-y) (plus its PDF), the public [2009 CM theorem PDF](https://bpb-us-w2.wpmucdn.com/sites.udel.edu/dist/a/6387/files/2017/08/pbentJune09-fin-1eyo5xc.pdf) as cited by upstream, the [2019 workshop abstract](https://boolean.w.uib.no/files/2019/06/abstract15.pdf), the [2024 vectorial non-weak-regular paper's metadata](https://doi.org/10.1109/TIT.2024.3434481), and [Beierle–Beyne 2025 degree-bound work](https://escholarship.org/uc/item/6697p4nj). The latter explicitly concerns complete **monomial** classification in named ranges, not all planar maps. The **2019 abstract already says all then-known planar families yielded weakly regular components**: positive claims about known families are **PRIOR ART**, not novel. The 2025 survey explicitly lists the 8 known `GF(3^5)` classes including the two sporadic classes absent from the upstream binomial search. No accessible universal proof, planar counterexample, or exact independent publication of this particular 19-class computational log was identified. **Not-found ≠ new.**

**Access gaps in priority:** v1/v2 arXiv PDFs were inaccessible, no trusted v1-v2-v3 diff; journal correction tab and full IEEE print inaccessible; not all ~190 citing works individually read; full citation graph not obtained. These gaps preclude an “every version corrected” assertion and a proof of priority.

## Evidence and replay

From the root of this repository, with **CPython 3 standard library only**:

```sh
python3 codesdeep/errata.py > /tmp/codesdeep-errata.log
cmp codesdeep/errata.log /tmp/codesdeep-errata.log
python3 codesdeep/check_errata.py > /tmp/codesdeep-check-errata.log
cmp codesdeep/check_errata.log /tmp/codesdeep-check-errata.log
python3 codesdeep/planar_classes.py > /tmp/codesdeep-planar.log
cmp codesdeep/planar_classes.log /tmp/codesdeep-planar.log
python3 codesdeep/check_planar.py > /tmp/codesdeep-check-planar.log
cmp codesdeep/check_planar.log /tmp/codesdeep-check-planar.log
python3 codesdeep/exhaustive_prime.py > /tmp/codesdeep-prime.log
cmp codesdeep/exhaustive_prime.log /tmp/codesdeep-prime.log
(cd codesdeep && sha256sum -c SHA256SUMS)
```

`planar_classes.py` constructs an irreducible quotient for each field, independent additive/multiplicative tables, a primitive generator, the trace, all derivatives and all scalar components. `check_planar.py` deliberately uses tuple arithmetic rather than importing that program. `errata.py` reconstructs GF(9) and both corrected example enumerators, while `check_errata.py` independently reconstructs the same witnesses by a separate tuple implementation. The `GF(9)` all-map backtracking is **not rerun** in codesdeep; refer to the explicit upstream programs and hashes. `SHA256SUMS` covers this report, every new source and every new log; it does not include itself to avoid self-reference. All program outputs are deterministic. No proprietary numerical software, paid journals or commercial API used.

## Not covered

- **Not** all planar maps on `GF(27)`, `GF(81)`, `GF(243)`, `GF(25)` or `GF(125)`; only known inequivalent representatives and mathematical reductions, not an all-map classification. Nor every coefficient variant of a named family *as a separate finite map*, though the class-level proof above covers EA variants and every scalar coefficient for each representative.
- **Not** arbitrary field degrees, higher-prime/characteristic infinite-family enumeration, all future sporadic discoveries, or the full weak-regularity question in all fields.
- **Not** a fully sourced v1/v2-v3 comparison, direct IEEE printed page check, publisher erratum audit, or complete global citation/priority search.
- **Not** a proof that `[12,4,6]_3` is suboptimal (only a proof the cited Griesmer bound cannot establish optimality), nor a proof that all printed numerical examples beyond the two flagged are universally correct.
- **Not** a new independent reconstruction of Helleseth et al.'s full 2009 Coulter–Matthews theorem; the full DO proof is written above.

## Next steps

1. **Most promising:** search for a *genuinely nonquadratic, non-Coulter–Matthews planar map* (or prove a finite exclusion) in a field **not already completely enumerated**. Start with a symmetry-reduced, verifiably exhaustive `GF(27)` search or a certified degree-restricted `GF(243)` search; independent planarity/weak-regularity checkers and a stated orbit-completeness argument are essential. A found **planar** component with changing Walsh sign would settle the universal question negatively, but no such example is claimed.
2. Obtain authenticated arXiv **v1 and v2** and publisher print to produce a version-by-version errata diff; mark exact corrections only after direct comparison.
3. Resolve the `[12,4,6]_3` optimality argument by a standalone nonexistence certificate for `[12,4,7]_3` or an explicit counter-code. Do not infer it from Griesmer alone.
4. Extend exact cyclotomic checking to the new **2025 sporadic `GF(3^6)` classes** (seven classes) and their nontrivial isotopic splits, with class-specific witnesses; distinguish computations from the already-published DO theorem.
