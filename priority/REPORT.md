# Independent priority audit: Tang et al. (2016)

**Date:** 2026-11-15. **Client:** Athanor lab (Aleksei Kudriashov).

**Target:** C. Tang, N. Li, Y. Qi, Z. Zhou and T. Helleseth, "Linear codes with two or three weights from weakly regular bent functions," *IEEE Transactions on Information Theory* 62 (2016), 1166-1176; [DOI 10.1109/TIT.2016.2518678](https://doi.org/10.1109/TIT.2016.2518678); [arXiv:1507.06148v3 PDF](https://pdfs.assets.alphaxiv.org/1507.06148v3.pdf), 14 August 2015. All page numbers below refer to the **numbered 16-page arXiv v3**, not journal pagination.

## Audit standard and verdicts

"Already published" requires an accessible reference expressly stating the same result/error; "partly published" means an equivalent formula or genuinely narrower result exists in published material but not the complete finding; "not found" reports only a bounded negative literature search, **not proof of novelty**. Repository computations are primary evidence of the finding, not independently published prior art.

| Claim | Verdict | Reference and exact location |
|---|---|---|
| **1. Missing even-m factor in lengths** | **PARTLY PUBLISHED** | Tang et al., arXiv v3 **Theorem 3.1 and Corollary 3.3, p. 10**, omit the factor; **equation (2), p. 2; Lemma 2.7, p. 4; Tables I and III, pp. 10-11**, already incorporate the correct sign. No explicit later correction, erratum or comment located. |
| **2. False unqualified m=2 dimension assertions** | **NOT FOUND** | No publication correcting Tang et al. **Theorem 3.1/Corollary 3.3, pp. 10-11, or Theorem 4.1/Corollary 4.4, pp. 11-13**, using these GF(9) endpoint witnesses located. |
| **3. Two impossible example enumerators** | **NOT FOUND** | No published notice correcting the examples after **Theorems 4.2 and 4.3, p. 12**, located. |
| **4. Weakly regular components in exhaustive/sparse planar tests** | **PARTLY PUBLISHED** | Helleseth et al., *IEEE TIT* **55** (2009), 5272-5283, **Theorem 1.2** proves the Coulter-Matthews family, [DOI](https://doi.org/10.1109/TIT.2009.2030465). Cesmelioglu, Meidl and Pott, [2019 public workshop abstract](https://boolean.w.uib.no/files/2019/06/abstract15.pdf), paragraph beginning "The classical examples", states components of all **then-known** planar families are weakly regular. Neither gives the exact enumerations or a universal theorem. |

## Claim 1: published ingredients but no explicit correction found

Writing chi = (Legendre(-1,p))^(m/2), the correct even-m lengths are

    n = p^(m-1)-1 + epsilon*(p-1)*chi*p^((m-2)/2)
    n_punctured = (p^(m-1)-1)/(p-1) + epsilon*chi*p^((m-2)/2).

Theorem 3.1 and Corollary 3.3 print these without chi. For f in the authors' RF class, f*(0)=0 and equation (2) gives W_f(0)=epsilon*(sqrt(p*))^m=epsilon*chi*p^(m/2), since p*=Legendre(-1,p)*p. Orthogonality gives p*#{x:f(x)=0}=p^m+(p-1)*W_f(0), proving the first signed length upon removing x=0; projective puncturing divides by p-1. Hence the error is real whenever p=3 mod 4 and m=2 mod 4.

**Important correction to the claim wording:** The original **Table I (p. 10) and Table III (p. 11) already CONTAIN the Legendre sign in weights and multiplicities**. They do **not** inherit the omission; the displayed theorem/corollary parameters contradict the tables. The proof and signed table entries are already in the 2015/2016 source, which merits a **partly published** priority classification, but we found no published erratum flagging the missing factor.

Independent finite witness: [verify/REPORT.md](../verify/REPORT.md) and [verify/check_witness.py](../verify/check_witness.py) construct GF(9)=F_3[t]/(t^2+1), f(x)=Tr(x^2), epsilon=-1; actual zero-set length 4 (punctured length 2), against printed 0 (punctured 0). [attack/REPORT.md](../attack/REPORT.md) provides an independent argument and another larger-field witness.

Search for a statement explicitly calling out this typo: article title/DOI/author combination plus "correction", "erratum", "comment", "Theorem 3.1", "Legendre", "sign", "code length", and targeted citing-journal papers. **No exact published correction found**.

## Claim 2: missing dimension qualifications at m=2

The claim is about the *dimension* of trace evaluation codes, not about whether the example functions are bent. A reproducible GF(9) witness in [verify/witnesses.jsonl](../verify/witnesses.jsonl) / [verify/REPORT.md](../verify/REPORT.md):

- For primitive w=1+t, f(x)=Tr(w*x^2) is weakly regular bent and in RF, but its zero-level defining set is **empty**: code dimension **0**, not asserted m=2.
- For f(x)=Tr(x^2), each square- and nonsquare-level defining set has two elements yet yields code **dimension 1**, not 2.

This is an endpoint failure, not an attack on the tested m>=3 cases. Later papers using their own m>=3/m>=4 hypotheses are not automatically corrections to Tang's explicit m=2 statements. Searches joining title, theorem numbers, GF(9), "rank", "dimension", "degenerate", "counterexample", "correction", and "erratum" found no accessible paper stating these particular witnesses or a specific repair. **Verdict: NOT FOUND.**

## Claim 3: the two p=3, m=5 examples

In arXiv v3 **p. 12**, the example immediately after **Theorem 4.2** prints a ternary [60,5,40] code and enumerator

    1 + 24*z^40 + 40*z^48 + 60*z^52

while the example immediately after **Theorem 4.3** prints [40,5,28] and

    1 + 40*z^28 + 60*z^32 + 24*z^40.

Each sums to **125** codewords, while a ternary 5-dimensional linear code necessarily has **3^5 = 243**. The numerical verification in [verify/REPORT.md](../verify/REPORT.md) and [verify/check_examples.py](../verify/check_examples.py) shows the printed enumerators are consistent with the **p=5,m=3, dimension 3** families, not the printed p=3,m=5 settings. This strongly suggests a swapped field parameter and printed dimension but is not evidence of author intent. Title/DOI and example-specific term searches in journals and public citation listings found no independent published correction. **Verdict: NOT FOUND.**

## Claim 4: exact computational evidence versus earlier mathematical results

Target paper **Section II-B, p. 2** asks whether the scalar functions Tr(c F(x)) derived from **every planar F:GF(p^m)->GF(p^m)** are weakly regular bent for all c!=0. Planarity alone ensures bentness of the components; it does **not** directly settle weak regularity.

**Reproducible unpublished-in-the-literature-to-us counts:**

- [attack/exhaustive_planar.c](../attack/exhaustive_planar.c), [log](../attack/exhaustive_planar.log), and [independent backtracking check](../attack/check_gf9.py): all **4,782,969** maps GF(9)->GF(9) normalized by F(0)=F(1)=0, **1,296** planar maps, **10,368** nonzero scalar components, all weakly regular (5,184 of each common Walsh sign, zero failures). Normalization subtracts an affine map, preserving both planarity and weak regularity; this exhausts normalized GF(9) maps only.
- [attack/binomial_probe.py](../attack/binomial_probe.py) and [log](../attack/binomial_probe.log): **restricted** polynomial family F(x)=x^2+a*x^d, a!=0, 3<=d<=q-2. GF(27): 598 candidates, **78 planar**, **2,028/2,028** components weakly regular. GF(81): 6,160 candidates, **420 planar**, **33,600/33,600** components weakly regular. This does not classify all planar GF(27) or GF(81) maps. The qualifying exponents are additive/Dembowski-Ostrom-type.

**Exact earlier statements checked:**

1. T. Helleseth, H. D. L. Hollmann, A. Kholosha, Z. Wang and Q. Xiang, ["Proofs of two conjectures on ternary weakly regular bent functions"](https://bpb-us-w2.wpmucdn.com/sites.udel.edu/dist/a/6387/files/2017/08/pbentJune09-fin-1eyo5xc.pdf), *IEEE TIT* **55** (2009), 5272-5283: **Lemma 1.1** says planar iff all scalar components are bent; **Theorem 1.2** proves weak regularity of Coulter-Matthews trace components. Quadratic planar functions have weakly regular components by quadratic Gauss sums. This is **family-specific**, not universal.
2. A. Cesmelioglu, W. Meidl and A. Pott, ["Vectorial bent functions in odd characteristic and their components"](https://boolean.w.uib.no/files/2019/06/abstract15.pdf), **2019 workshop abstract**, says components of **all known planar functions** were weakly regular, noting those families are quadratic except Coulter-Matthews. This is a prior **observation about known constructions**, not a general proof or an enumeration of 1,296 maps.
3. L. Koelsch and A. Polujan, ["Value Distributions of Perfect Nonlinear Functions"](https://doi.org/10.1007/s00493-023-00067-y), *Combinatorica* **44** (2024), 231-268, **Section 6, Proposition 6.1, Theorem 6.4**, proves value-distribution/plateaued-planar results under stated assumptions; not universal weak regularity.
4. A. Cesmelioglu and W. Meidl, ["Vectorial Bent Functions With Non-Weakly Regular Components"](https://doi.org/10.1109/TIT.2024.3434481), *IEEE TIT* **70** (2024), 9214-9226, **abstract and Introduction**, produces vectorial bent functions with non-weakly-regular components. **The constructed vectorial functions are not shown planar** and are NOT counterexamples to Tang's question.
5. Y. Wu and Y. Pan, ["Linear codes from planar functions and related covering codes"](https://doi.org/10.1016/j.ffa.2024.102535), *Finite Fields Appl.* **101** (2025), article 102535: related planar-function code constructions; no identified general component-weak-regularity proof.
6. *Handbook of Finite Fields*, **Theorem 9.3.66**, states planarity iff every nonzero trace component is (generalized) **bent**; importantly this is weaker than "weakly regular bent."

Searches for the exact counts "1296" with GF(9), "2028" with GF(27), and "33600" with GF(81), and for planar counterexamples and universal theorems, found no exact publication. **Verdict: PARTLY PUBLISHED** (family-specific positive results); **no universal proof or genuine planar-function counterexample verified** as of this search. Do not interpret "not found" as a mathematical proof that the question remains unresolved worldwide.

## Search protocol, hits and important access gaps

| Venue / database | What was searched or read | Limitation |
|---|---|---|
| **IEEE Xplore and IEEE TIT index** | DOI 10.1109/TIT.2016.2518678, full title and authors with "erratum", "correction", "corrigendum", "comment on", theorem and example terms; target's bibliographic page | Direct IEEE Xplore corrections/comments tab not accessible; none found in indexed results. **Publisher database not exhaustively checked.** |
| **arXiv** | Read original **v3** PDF and visually inspected numbered **pp. 2, 10, 12**, including formula, tables and two example enumerators | Attempts to retrieve **v1 and v2** via arxiv.org, export.arxiv.org and alphaXiv failed. **No v1-v2-v3 textual diff completed**; version identity is not claimed. |
| **Google Scholar, Semantic Scholar, OpenAlex, citing graph** | Publicly indexable citation searches by DOI/title, search hits on Scispace and citing-paper bibliographies (the surfaced Scispace record displayed around **190 citations**) | Complete Google Scholar/Semantic Scholar/OpenAlex citing-paper export inaccessible. **The approximately 180-190 citing papers were not all read.** |
| **Coding-theory journals** | Indexed articles and accessible statements in *IEEE TIT*, *Designs, Codes and Cryptography*, *Finite Fields and Their Applications*, *Cryptography and Communications*, *Advances in Mathematics of Communications* | Targeted sample, **not every full text**. Relevant non-correction articles include S. Mesnager, *Cryptogr. Commun.* **9** (2017), 71-84, DOI 10.1007/s12095-016-0186-5; F. Ozbudak and R. M. Pelen, *IEEE TIT* **68** (2022), 3014-3027, DOI 10.1109/TIT.2022.3145337; Koelsch/Polujan (2024); Wu/Pan (2025). |
| **MathOverflow and Math StackExchange** | Targeted site-indexed searches with "planar", "weakly regular", "Tang", "two/three weights", "GF(9)", "counterexample", "erratum" | No germane explicit correction, proof or planar counterexample located; not an exhaustive crawl. |
| **Exact-phrase searches** | Factor / Legendre symbol and theorem numbers; two enumerators and p=3,m=5; GF(9) "1296", GF(27) "2028", GF(81) "33600" | No publicly indexed exact match outside the Athanor repository returned. |

**Priority and reproducibility limits.** Evidence in [verify/REPORT.md](../verify/REPORT.md), [attack/REPORT.md](../attack/REPORT.md), and [status/REPORT.md](../status/REPORT.md) is independently testable, but is neither third-party publication nor proof of literature novelty. A finite computation cannot prove the arbitrary-field planar question. Searches without publisher-fulltext and full citing-graph access must be characterized as **bounded**. The first three findings should be reported as **apparently uncorrected**, pending a publisher check and a complete citations audit.

**Bottom line:** (1) PARTLY PUBLISHED, original tables show the sign but explicit correction not found; (2) NOT FOUND; (3) NOT FOUND; (4) PARTLY PUBLISHED for specific planar families only, exact enumeration not found, universal question not verified solved.
