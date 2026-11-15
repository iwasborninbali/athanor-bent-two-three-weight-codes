# State of the art and reproducible ideas: Tang--Li--Qi--Zhou--Helleseth

Date: 2026-11-15. Status: literature-based research note with **limited independently
replayed finite computations**, not a claim of a new theorem.

## Scope

Primary source: Chunming Tang, Nian Li, Yanfeng Qi, Zhengchun Zhou, Tor Helleseth,
*Linear codes with two or three weights from weakly regular bent functions*,
IEEE Transactions on Information Theory **62** (2016), 1166--1176,
DOI [10.1109/TIT.2016.2518678](https://doi.org/10.1109/TIT.2016.2518678).
Read the 16-page arXiv v3 PDF, dated August 14, 2015,
[arXiv:1507.06148](https://arxiv.org/pdf/1507.06148).
**Page numbers below refer to the numbered arXiv PDF**, not the journal pagination.
The paper studies $p$ odd, $q=p^m$, $f:\mathbb F_q\to\mathbb F_p$ weakly
regular bent in the subclass $RF$ with $f(0)=0$ and
$f(ax)=a^hf(x)$ for $a\in\mathbb F_p^*$, $\gcd(h-1,p-1)=1$.
It considers the defining set $D_{f,0}=\{x\in\mathbb F_q^*:f(x)=0\}$
and trace code $C_D=\{(\operatorname{Tr}(\beta d))_{d\in D}:\beta\in\mathbb F_q\}$,
as well as nonzero square/nonsquare levels and puncturings.

Question inventory is deliberately separated into (A) an open question explicitly
stated by the authors, (B) an inherited open problem claimed solved by the paper,
(C) an older conjecture recalled as already proved, and (D) future work rather
than a formal conjecture. This prevents assigning spurious 2016-onward solutions.

Method: inspected source PDF pp. 1--2, 10--16, including the displayed
formulas/tables; searched public arXiv, journal pages, preprints, indexed citing
literature and targeted MathOverflow, Math StackExchange and OEIS searches.
OpenAlex public API access via the available web interface returned an error;
therefore no complete Semantic Scholar/OpenAlex citation-graph audit is claimed.
No paid services, Magma, SaaS computation or external executable code were used.

## Findings

### A. Open: universal weak regularity of components of planar functions

**Verbatim**, primary paper Section II.B, p. 2, last third of page:

> It is still an open problem whether the derived p-ary function from any a planar function is a weakly regular bent function.

Mathematically precise question: for *every* planar function
$F:\mathbb F_{p^m}\to\mathbb F_{p^m}$ and every $b\ne0$, is
$f_b(x)=\operatorname{Tr}_{\mathbb F_{p^m}/\mathbb F_p}(bF(x))$
weakly regular bent? The source has a harmless but confusing typo around the
definition of a planar mapping; the domain and codomain for $F$ are both
$\mathbb F_q$. Planarity already guarantees bentness of each $f_b$;
**weak** regularity is the issue.

**Status: OPEN / no universal proof or planar counterexample verified** as of this
search. The latest literature examined does **not** justify the assertion that
this was solved. In particular:

* A. Cesmelioglu and W. Meidl, *Vectorial Bent Functions With Non-Weakly
  Regular Components*, IEEE TIT **70** (2024), 9214--9226,
  [DOI](https://doi.org/10.1109/TIT.2024.3434481), proves existence of
  vectorial bent functions with non-weakly-regular components and studies
  Walsh distributions and codes. **Not** a counterexample to the quoted
  problem: a vectorial bent function $\mathbb F_p^n\to\mathbb F_p^r$ with
  arbitrary $r$ need not be planar ($r=n$).
* L. Koelsch and A. Polujan, *Value Distributions of Perfect Nonlinear
  Functions*, Combinatorica **44** (2024), 231--268,
  [DOI](https://doi.org/10.1007/s00493-023-00067-y), Section 6,
  Proposition 6.1 and Theorem 6.4, develops constraints and characterizations
  of planar functions (not a proof of universal component weak regularity).
* The earlier [2019 vectorial-bent abstract](https://boolean.w.uib.no/files/2019/06/abstract15.pdf)
  explicitly notes that the components of all then-known planar functions
  were weakly regular, while discussing non-weakly-regular components of
  more general vectorial functions. This is historical evidence, not an
  up-to-date exhaustive classification.

**Evidence type:** verified source statement plus identified later
papers/abstracts; *negative literature search only*, not a proof of openness.

### B. Already solved in Tang et al.: Ding--Ding weight distributions

Section I, p. 1, and Section III, p. 10, say that K. Ding and C. Ding asked
for the distribution of the trace code defined by
$f(x)=\operatorname{Tr}(F(x))$ for general planar $F$.
Tang et al. explicitly say they solve that problem. The precise theorems
are Section III **Theorem 3.1, Table I (even $m$)** and **Theorem 3.2,
Table II (odd $m$)**, pp. 10--11, under their $RF$ hypotheses; their
projective puncturing is **Corollaries 3.3--3.4, Tables III--IV**, pp. 10--11.
Section IV **Theorems 4.1--4.3, Tables V--VII**, pp. 11--12, treats
square/nonsquare defining sets; **Corollaries 4.4--4.6, Tables VIII--X**,
pp. 12--13, treats the corresponding puncturings.
Source of the inherited question: K. Ding and C. Ding,
[A class of two-weight and three-weight codes and their applications in
secret sharing](https://arxiv.org/abs/1503.06512) (2015).

**Status: SOLVED IN THE TARGET PAPER within the stated weak-regular/homogeneous
class, NOT a proof for all hypothetical planar functions.** The universal
planar-component question A is still logically separate. Thus the paper's
wording should not be read as an unconditional resolution of A.

**Evidence type:** primary source, exact theorem/table locations; finite
checks below cover two representative instances, not the entire theorem.

### C. Older Helleseth--Kholosha conjecture: solved before 2016

Section V.B, pp. 14--15, quotes a 2006 conjecture on the weak regularity and
specific Walsh transform of the ternary monomial
$f(x)=\operatorname{Tr}_1^{2k}(c x^{(3^{2k}-1)/4+3^k+1})$
for $k$ odd with a special coefficient $c$; its complete mathematical
statement is **equations (6) and (7)** on pp. 14--15.
Tang et al. state that Helleseth, Hollmann, Kholosha, Wang and Xiang proved
the weak regularity in 2009 using Stickelberger's theorem, while Gong,
Helleseth, Hu and Kholosha established the stated dual/Walsh formula in 2012.
Bibliography references **[13]** and **[12]**, respectively:

* T. Helleseth et al., *Proofs of two conjectures on ternary weakly regular
  bent functions*, IEEE TIT **55** (2009), 5272--5283,
  [DOI](https://doi.org/10.1109/TIT.2009.2030465).
* G. Gong et al., *On the Dual of Certain Ternary Weakly Regular Bent
  Functions*, IEEE TIT **58** (2012), 2237--2243,
  [DOI](https://doi.org/10.1109/TIT.2011.2176312).

**Status: SOLVED PRE-2016**, not a conjecture posed by Tang et al.
**Evidence type:** target paper's explicit retrospective attribution plus
bibliographic metadata; this project did not independently reconstruct the
Stickelberger/dual-function arguments.

### D. Generalizing to more general functions: partially addressed

Section VI, p. 16, ends with a future-work sentence about constructing
few-weight codes from more general $f(x)$. It is **not** a precise conjecture
with a unique yes/no resolution. Several different extensions now exist:

* S. Mesnager, *Linear codes with few weights from weakly regular bent
  functions based on a generic construction*, Cryptography and
  Communications **9** (2017), 71--84,
  [DOI](https://doi.org/10.1007/s12095-016-0186-5);
  open [preprint](https://eprint.iacr.org/2015/1103.pdf),
  Section 5, **Theorems 2--3, Tables 2--4**, pp. 15--16: a *different*
  first-generic-construction family with $m+1$ dimension and three weights.
  The preprint predates the 2016 journal article, so it must not be described
  as a post-2016 discovery, although the journal appeared in 2017.
* G. Xu, L. Qu, G. Luo, *Minimal linear codes from weakly regular bent
  functions*, SETA 2020, [open paper](https://seta-2020.etu.ru/assets/files/program/paper/paper-8.pdf),
  **Theorem 13 and Table 2, p. 9**: for even $m\ge4$, a different code
  $C_U$ of length $p^{2m}-1$ and dimension $2m+1$ is minimal even though it
  violates the Ashikhmin--Barg sufficient ratio test. It does **not** say
  that Tang's original $C_D$ is minimal.
* F. Ozbudak and R. M. Pelen, *Two or Three Weight Linear Codes from
  Non-Weakly Regular Bent Functions*, IEEE TIT **68** (2022), 3014--3027,
  [DOI](https://doi.org/10.1109/TIT.2022.3145337): builds codes from
  selected signs/dual-bent components, obtains two/three weights,
  investigates minimality, and gives distributions for a generalized
  Maiorana--McFarland subclass. This **partially answers the broader
  research direction** but not the planar-function question A.
* R. M. Pelen, *Three weight ternary linear codes from non-weakly regular
  bent functions*, Advances in Mathematics of Communications **18** (2024),
  711--737, [DOI](https://doi.org/10.3934/amc.2022020),
  published online in 2022, **Tables 1--4**: flexible defining sets from
  non-weakly-regular dual-bent functions.
* W. Jin, K. Li, L. Qu, *Several classes of minimal linear codes from
  weakly regular and non-weakly regular bent functions*,
  Discrete Applied Mathematics **373** (2025), 53--76,
  [DOI](https://doi.org/10.1016/j.dam.2025.04.044),
  **Sections 3--4**: four-to-six-weight minimal codes violating the
  Ashikhmin--Barg condition, from both classes of bent functions.
  This is an alternative construction, not a minimality theorem for the
  original Tang et al. trace codes.
* Additional nearby work: S. Yang, T. Zhang, P. Li,
  *Linear Codes from Two Weakly Regular Plateaued Balanced Functions*,
  Entropy **25** (2023), article 369,
  [DOI](https://doi.org/10.3390/e25020369), Sections 4--5;
  and the 2024 vectorial paper above.

**Status: PARTLY SOLVED AS A BROAD PROGRAM**, with identifiable advances and
clear differences in construction. There is no single theorem covering
all $f$ or resolving the planar question.

### E. Minimality of the original defining-set codes

The 2016 article gives weights, distances and secret-sharing/graph applications;
it **does not state a general minimality criterion** for its original
$C_{D_{f,0}}$ and its square/nonsquare variants. Being two-weight or
three-weight does not by itself imply minimality. The cited later minimal-code
papers mostly use different constructions; their statements must not be
transferred to the original codes without a proof.

**Status: finite instances decisively NONMINIMAL**, with complete enumerations
and independently checked explicit witnesses below. General classification
of the original families: **NOT COVERED / not asserted solved**.

## Evidence

Everything in this section can be replayed with CPython 3 and its standard
library; no installation or network access is needed:

```sh
cd status
python3 enumerate_codes.py > /tmp/enum.new.log
cmp /tmp/enum.new.log enumerate_codes.log
python3 check_witness.py > /tmp/check.new.log
cmp /tmp/check.new.log check_witness.log
sha256sum -c SHA256SUMS
```

Files: `enumerate_codes.py` (direct polynomial-basis field arithmetic),
`enumerate_codes.log`, `witness.json` (all codewords and first support
violations), `check_witness.py` (independent, **field-free** checker of
GF(3)-linear spans, closure, dimensions, weight counts, support inclusions),
`check_witness.log` and `SHA256SUMS`.

Fields are $\mathbb F_3[t]/(t^3+2t+1)$ and
$\mathbb F_3[t]/(t^4+t+2)$. The first modulus is irreducible because
it has no root in $\mathbb F_3$. The second has no $\mathbb F_3$ root
and remainders $t$, $t+1$, $t+1$ upon division by the three monic
irreducible quadratics $t^2+1$, $t^2+t+2$, $t^2+2t+2$, respectively;
thus it is irreducible. The enumerator independently checks the absolute
squared Walsh magnitudes using the exact integer identity
$|a+b\omega+c\omega^2|^2=a^2+b^2+c^2-ab-bc-ac$ for a primitive
cube root $\omega$ (all 3^m spectral points, no floating point).

Actual full-code enumerations for $f(x)=\operatorname{Tr}(x^2)$:

| $p,m$, level $f(x)$ | $[n,k]$ | complete weight enumerator | minimal? |
|---|---|---|---|
| $3,3$, 0 | $[8,3]$ | $1+12z^4+8z^6+6z^8$ | no |
| $3,3$, 1 (squares) | $[6,3]$ | $1+6z^2+12z^4+8z^6$ | no |
| $3,3$, 2 (nonsquares) | $[12,3]$ | $1+8z^6+6z^8+12z^{10}$ | no |
| $3,4$, 0 | $[20,4]$ | $1+60z^{12}+20z^{18}$ | no |

The level-zero results reproduce **Theorem 3.2/Table II** and
**Theorem 3.1/Table I**, respectively. Level-one/two cases are additional
finite checks of the defining sets in **Section IV**, not theorem-wide
claims. All enumerated trace maps are injective ($3^k$ distinct words).

Explicit witnesses of nonminimality are in `witness.json` and the logs.
For example, for $m=3$, level zero, the containing word
`[2,1,2,1,1,2,2,1]` and contained word `[2,1,0,0,1,2,0,0]` belong to
the enumerated code, are not scalar multiples, and the second's support
is contained in the first's. For $m=4$, level zero, the log gives a
20-coordinate pair and a separate checker validates the inclusion and
non-proportionality. This is a full **finite** certificate of nonminimality,
not a computer-assisted assertion about all $m$. `check_witness.py` does
not import or reuse field operations from `enumerate_codes.py`.

**Reasoning behind the finite minimality certificate:** by definition a
linear code is minimal only if any nonzero codeword whose support is
contained in that of another is a scalar multiple of the latter.
A single verified pair of nonproportional words with that support
inclusion disproves minimality for the corresponding specific code.
This is the complete argument used after codeword membership has been
validated by the independent linear-span checker.

## Not covered

* No global theorem about universal planar component weak regularity, and
  no non-weakly-regular **planar** function found or checked. "Open" means
  no verifiable resolution in this bounded search, not a proof of absence.
* No full formal proof of the 2016 general weight-distribution formulas,
  2009/2012 Walsh conjecture, or any later paper's claims; we audited
  their stated results and exact pointers. Their proofs remain external.
* No computation for large $m$, arbitrary odd primes, nonquadratic planar
  families, or all members of $RF$. No independent finite-field package
  implementation; the second checker verifies linear codeword structure,
  not the field multiplication table from first principles.
* No exhaustive 2016--2026 citation graph: the public OpenAlex endpoint
  returned an access error in this interface, and no exhaustive Semantic
  Scholar citation list was obtained. Targeted MathOverflow, Math
  StackExchange and OEIS queries gave no germane resolution, but this is
  negative search evidence only.
* No claim that the 2022/2024 vectorial non-weak-regular examples are planar;
  no claim that later minimal codes have the same defining sets as Tang's.
* No paid journal full-text or computational services were called. The
  2022 IEEE and 2024 IEEE status above is based on accessible public abstracts
  / metadata and later references where full text was not freely retrieved.

## Next steps

1. **Minimality classification for the original six defining-set families.**
   Derive necessary/sufficient support containment conditions from the
   dual Walsh orbit data (rather than the weight ratio only); seek a
   theorem in $m,p,\varepsilon$ and separate even/odd cases. **First step:**
   augment the independent checker to enumerate all projective codeword
   pairs for $p=3$, $m=3,4,5$ and all square/nonsquare puncturings;
   publish first failing pair for each negative case, plus full outputs.
2. **Targeted planar-function weak-regularity search.** Search small fields
   for genuinely planar $F$ outside the known quadratic and Coulter--Matthews
   templates, or prove a bounded classification. **First step:** implement
   two independent planarity tests (every derivative $F(x+a)-F(x)$ is a
   permutation for $a\ne0$; a direct difference table), then exact
   cyclotomic Walsh sign checks for every nonzero component. Export
   coefficients, multiplication polynomial and explicit witnesses.
3. **Remove Tang's homogeneity constraint.** For weakly regular bent $f$
   not in $RF$, track the character sums of the scalar orbits of the
   dual and determine when the zero-set trace code still has only two
   or three weights. **First step:** derive a general formula for
   $\operatorname{wt}(c_\beta)$ in terms of the $p-1$ Walsh transforms
   of $sf$ and evaluate it for one published bent example outside $RF$.
4. **Sharpen code bounds and projectivity.** For the projectively punctured
   codes, investigate exact dual distance and whether the two-weight
   constructions yield previously unknown strongly regular graph
   parameters. **First step:** compute exact dual weight enumerators
   from the full enumerator/MacWilliams transform in small cases and
   cross-check by explicit dual-kernel enumeration; list parameters
   and comparison with known tables before claiming novelty.
5. **Mechanically check a tractable theorem fragment.** Port the
   trace-code construction and finite Walsh orthogonality identity into
   Lean or another open proof assistant, separating field-independent
   characters from the case split $m$ odd/even. **First step:** formalize
   the weight-counting identity for a finite character group and the
   definition of minimality; commit a complete proof script and build log.

No direction above is yet claimed a new result. Every proposed extension
requires program+log+SHA-256 and an independent witness/checker where
possible, or a complete written argument for a mathematical theorem.
