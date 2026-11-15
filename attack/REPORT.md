# Computational attack: two-/three-weight codes from weakly regular bent functions

## Scope

**Source examined:** C. Tang, N. Li, Y. Qi, Z. Zhou, and T. Helleseth, *Linear codes with two or three weights from weakly regular bent functions*, arXiv:1507.06148v3 (14 August 2015), https://arxiv.org/pdf/1507.06148v3; published in *IEEE Transactions on Information Theory* **62** (2016), 1166–1176. Page numbers below are the **printed arXiv PDF page numbers**, not the journal pagination. Its PDF has 16 pages. No paid services were used.

**Inventory of problems and conjectures stated/discussed in that paper:**

1. **Explicitly open at writing, universal weak regularity of planar components.** Section II-B, **p. 2**, last paragraph above the definition of `RF` (arXiv PDF page index 1):
   > It is still an open problem whether the derived p-ary function from any a planar function is a weakly regular bent function.

   Precise interpretation: for every planar map `F: GF(p^m) -> GF(p^m)` in odd characteristic and every `c != 0`, must `Tr(c*F(x))` be weakly regular bent? The paper already states that the components are bent, but weak regularity is the unresolved part. **Status of this attack: no counterexample found within the explicit ranges below.** This must not be read as a resolution for all finite fields. We did not establish the present-day literature status of this universal question.

2. **Ding–Ding weight-distribution problem (reported as solved).** Introduction, **p. 1**, and Section III immediately before Theorem 3.1, **p. 10**: determine the weight distribution of `C_{D_f}` with `D_f={x != 0 : Tr(F(x))=0}` for general planar `F`, extending `F(x)=x^2`. The authors say they solve the question by working with their class `RF` of weakly regular homogeneous bent functions; see **Theorems 3.1–3.2, pp. 10–11, Tables I–II**. This is a *historical motivating problem*, **not an open conjecture left by the authors**. A universal, unconditional reading for all planar `F` requires the separate weak-regularity issue above and/or another argument. We independently tested the formulas in specific families and found the *printed length* inconsistency below.

3. **Helleseth–Kholosha historical conjecture (explicitly described as proved).** Section V-B, **pp. 14–15, equations (6)–(7)**, `m=2k`, `k` odd, `p=3`, primitive `alpha`; it asserts weak regularity and a prescribed Walsh transform for
   `f(x)=Tr_{3^m/3}(c*x^((3^m-1)/4+3^k+1))`, `c=alpha^((3^k+1)/4)`.
   The transform has leading coefficient `-3^k` and a specified `GF(3)` phase in equation (7). The paper explicitly records the 2009 proof of bentness and the **2012 proof of the conjecture by Gong et al.**, reference [12]. This is **not open in the paper**. We test exact Walsh amplitudes and a common phase sign for `k=1,3`; we **do not check the full dual-phase expression in (7)**.

4. **Future research direction, not a formal conjecture.** Section VI, **p. 16**, final sentence: extend constructions of few-weight codes to more general functions `f(x)`. There is no precise quantified conjecture or specified search range. Our sparse-binomial probes are a limited exploration of nearby examples, not a test of a falsifiable universal assertion in that sentence.

No other explicit unresolved conjecture was located in the paper. In particular, its minimality of codewords is *not asserted as a universal theorem or open conjecture*; we test it because it is relevant to secret-sharing applications and to the task's requested subject.

## Findings

### F1 — Counterexample to the parameter statements in Theorem 3.1 / Corollary 3.3 **(exact computation + complete argument + independent witness)**

**Location:** **p. 10**, **Theorem 3.1** states that, for even `m`, `C_{D_f}` has length

`n_printed = p^(m-1) - 1 + epsilon*(p-1)*p^((m-2)/2)`.

The **same omission** occurs in **Corollary 3.3**, p. 10, for the length of its one-representative-per-`GF(p)^*`-orbit puncturing. Both displayed length expressions omit the factor

`chi = ((-1)/p)^(m/2) = (-1)^(((p-1)/2)*(m/2))`.

**Correct lengths** (with the Walsh-sign convention of equation (2), p. 2) are

```
n_correct = p^(m-1) - 1 + epsilon*chi*(p-1)*p^((m-2)/2)
n_punctured_correct = (p^(m-1)-1)/(p-1) + epsilon*chi*p^((m-2)/2).
```

**Complete argument.** By the definition of `RF`, `f(0)=0`, and Proposition 2.4, p. 3, gives `f*(0)=0`. Equation (2), p. 2, therefore gives `W_f(0)=epsilon*(sqrt(p*))^m`, where `p*=((-1)/p)*p` and `epsilon` is `+1` or `-1`. For **even** `m`, this is the rational integer `epsilon*chi*p^(m/2)`. For every `a in GF(p)^*`, the cyclotomic automorphism `zeta_p -> zeta_p^a` fixes that rational number; consequently `sum_x zeta_p^(a*f(x)) = W_f(0)`. By root-of-unity orthogonality,

```
p * #{x : f(x)=0}
  = sum_x sum_{a in GF(p)} zeta_p^(a*f(x))
  = p^m + (p-1)*W_f(0).
```

Subtracting the zero element yields the first corrected length. For the punctured version, the defining set is a union of orbits of size `p-1` under multiplication by `GF(p)^*`, so divide by `p-1`. This proves the two corrected length expressions. They **differ** from the printed statements whenever `p = 3 (mod 4)` and `m = 2 (mod 4)`. The paper's **Table I, p. 10, and Table III, p. 11, already contain `chi`** in their multiplicities/weights, and agree with the correction. The issue is specifically the displayed *length parameters*, not the table entries.

**Independent explicit example (not relying on historical conjectures):** `GF(3^6)=GF(3)[t]/(t^6+t^5+t^4+1)`; represent elements with coefficients in `{0,1,2}`, and take `f(x)=Tr_{729/3}(x^2)`. The checker independently constructs the field by polynomial reduction, verifies its chosen generator has order 728, computes the trace values and codewords, and obtains:

```
counts of f(x)=0,1,2: 261, 234, 234
W_f(0) = 27, so epsilon = -1 because (sqrt(-3))^6 = -27
actual C_D length 260, dimension 6; weights 0:1, 162:260, 180:468
Theorem 3.1 printed length 224  [INCORRECT for this example]
actual punctured length 130; weights 0:1, 81:260, 90:468
Corollary 3.3 printed length 112 [INCORRECT for this example].
```

Replay `python3 attack/theorem_counterexample.py`, log `attack/theorem_counterexample.log`. The independent Helleseth–Kholosha sample in F4 also exhibits the mismatch with the *opposite* sign.

**Degenerate endpoint:** Theorem 3.1 is stated for all even `m` without an explicit lower bound. For the Helleseth–Kholosha function at `p=3,m=2,k=1` (field `GF(9)`, `c` a primitive element, exponent 6), the defining set is empty; the code is the zero-length zero code. The printed claim of a nondegenerate `[4,2]` two-weight code does not hold. The sign-corrected count is zero, but a lower-bound/nondegeneracy hypothesis is also required before claiming dimension `m` and *two distinct nonzero* weights.

### F2 — Exhaustive small-field check of the remaining planar-component question **(exact enumeration + independent enumeration)**

`attack/exhaustive_planar.c` enumerates **every normalized map** `F:GF(q)->GF(q)` with `F(0)=F(1)=0` for `q=3,5,7,9`, checking the defining condition that `x -> F(x+a)-F(x)` permutes the field for every nonzero `a`. For any planar map `F`, subtraction of the affine map `F(0)+x*(F(1)-F(0))` preserves planarity and weak regularity of all trace components (it merely shifts Walsh indices and multiplies by a common root of unity). Therefore **these normalized enumerations cover *all* planar maps** in those fields up to this affine subtraction.

```
field   normalized maps visited   planar maps   components tested exactly
GF(3)                   3                 2   --
GF(5)                 125                 4   --
GF(7)              16,807                 6   --
GF(9)           4,782,969             1,296   10,368 (all weak; 0 failures)
```

In the prime fields, the C search enumerates **planarity only**; it does *not* independently test weak regularity of their components. For `GF(9)`, the C test is exact in `Z[zeta_3]`: each Walsh coefficient must equal `s*3*zeta_3^j` for a **common** `s in {+1,-1}` across all frequencies. The program reports 5,184 occurrences of each `s` (equivalently, 5,184 of each epsilon in the paper's `m=2` convention). A **second Python implementation** (`attack/check_gf9.py`) represents field elements as pairs and uses *backtracking* rather than C's full base-9 enumeration, with an independent complex-Walsh check. It reproduces both the count and the **FNV-1a digest of the ordered planar-map list**: `4bedb2f3ee993315`. All eight components of the displayed planar witness are independently checked in `check_gf9.log`.

This **does not prove** that all planar-component functions are weakly regular in arbitrary degrees and characteristics; `GF(3^3)` and larger fields have vastly more maps than searched here.

### F3 — Targeted complete binomial-family checks in larger fields **(exact arithmetic)**

`attack/binomial_probe.py` exhaustively scans the **declared sparse family**

`F(x)=x^2+a*x^d`, `a in GF(3^m)^*`, `3<=d<=3^m-2`, for `m=3` and `m=4`, with all derivatives tested for permutation and **every component coefficient** `c in GF(3^m)^*` checked by exact three-phase cyclotomic Walsh comparison.

```
GF(27):    598 candidate polynomials, 78 planar,
           2,028 components passed; 0 weak-regularity failures.
           planar exponents: d=3:26, 6:13, 9:26, 18:13.
GF(81):  6,160 candidate polynomials, 420 planar,
          33,600 components passed; 0 weak-regularity failures.
          planar exponents: d=3:80, 6:40, 9:80, 10:30,
                            18:70, 27:80, 54:40.
```

All entries in this family are covered, not all degree-`m` maps. These exponents are additive and/or Dembowski–Ostrom-type; **no newly discovered nonquadratic planar family** is claimed.

### F4 — Replayed Coulter–Matthews nonquadratic functions and code tables **(exact Walsh + independent witness for nonminimality)**

**Source:** Section V-D, **p. 15**, equation (10), and **Theorem 3.1/Table I (p. 10); Theorem 3.2/Table II (pp. 10–11)**. Test `F(x)=x^14` (`(3^3+1)/2`), with `gcd(3,m)=1` and `f_c(x)=Tr(c*F(x))`, where `m=4,5,7`. For each degree, the code directly checks **every nonzero derivative** of `F` and exact Walsh phase counts for its tested components. For `m=4,5` it tests **all nonzero c**; for `m=7` it tests `c=1` and `c=g` (primitive). Their code weight distributions agree with the **table entries**, using the correct sign convention. All projective codeword supports are checked for proper containment.

| Field | c/sign | `C_D` length | Nonzero weights with multiplicities | Minimal? |
|---|---|---:|---|---|
| `GF(81)` | `c=1`, epsilon -1 | 20 | `12:60, 18:20` | **No** |
| `GF(81)` | `c=10` primitive, epsilon +1 | 32 | `18:32, 24:48` | **Yes** |
| `GF(243)` | `c=1` and `c=3` primitive, signs +1/-1 | 80 | `48:90, 54:80, 60:72` | **Yes** |
| `GF(2187)` | `c=1` and `c=3` primitive, signs +1/-1 | 728 | `468:756, 486:728, 504:702` | **Yes** |

In `GF(81)`, the direct minimality checker tested **all 80** nonzero `c`: 40 square coefficients gave nonminimal codes and all 40 nonsquares gave minimal codes. In `GF(243)`, **all 242** coefficients gave minimal codes. In `GF(2187)`, only the two named coefficients were explicitly calculated; the following elementary reduction covers the others in this **fixed monomial family**:

**Full square-class reduction argument.** The multiplicative group `GF(q)^*` is cyclic. For `q=81,243,2187`, `gcd(14,q-1)=2`. Thus the image of `u -> u^14` consists **exactly of squares**. If `c'/c = u^14`, then `Tr(c' x^14)=Tr(c (u x)^14)`. The change of variables `x -> ux` permutes the defining set, and the trace code is isometric after the corresponding reparameterization of `beta`. The Walsh spectra also differ only by a frequency permutation. Therefore every nonzero coefficient is equivalent to **one of the two representatives** `1` and a fixed primitive nonsquare `g`, and both weak regularity and code minimality are invariant within each class. This justifies the full-family classification for these specific degrees from the two representatives, **not** a general assertion about all planar functions.

**Explicit nonminimality witness:** Work over `GF(81)=GF(3)[t]/(t^4+t^3+t^2+1)`, encode `a0+a1*t+a2*t^2+a3*t^3` as `a0+3*a1+9*a2+27*a3`, let `f(x)=Tr(x^14)` and `D={x != 0:f(x)=0}`. The 20 defining-set field indices, in order, are

`[10,14,16,17,20,22,23,25,30,31,34,35,43,49,58,59,60,62,71,77]`.

Two nonproportional parameters `beta=1` and `beta=10=1+t^2` give codewords `c_beta=(Tr(beta*d))_{d in D}` of weights 12 and 18. Their **0-based** support-position sets are

```
S_1  = {2,3,5,6,8,9,10,11,14,15,16,17}
S_10 = {1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}
```

Hence `S_1` is a **proper subset** of `S_10`, while `1` is neither `10` nor `2*10` in the base field. This is a counterexample to the **extra** claim that every code in the family must be minimal; **the paper does not make that claim**. See `attack/minimality_witness.json` and the independently written direct-polynomial checker `attack/verify_witness.py`.

### F5 — Small-case check of the historical Helleseth–Kholosha conjecture **(exact Walsh; partial check of equation (7))**

`attack/helleseth_kholosha.py` implements **equation (6)** with `k=1,m=2` and `k=3,m=6`, primitive `alpha` found by finite-field enumeration and `c=alpha^((3^k+1)/4)`. It checks every Walsh frequency for a **common sign**, plus the claimed `W_f(0)=-3^k`. Results: `GF(9)` exponent 6, `c=4`, epsilon +1, `W0=-3`; `GF(729)` exponent 210, `c=703`, epsilon +1, `W0=-27`. The latter gives **length 224**, with weights `144:504, 162:224`, matching Table I and again contradicting the **printed length 260** in Theorem 3.1. At `GF(9)` the defining set is empty (degenerate endpoint). The **full dual-phase formula (7) is not checked**; the conjecture was already proved in the cited 2012 paper, not solved here.

## Evidence

All algorithms were written in this attack. Python programs use only the standard library; C program requires a C11 compiler. Field arithmetic uses explicitly generated irreducible polynomials and primitive elements, without Magma/Sage or paid services. All Walsh comparisons for `p=3` are **integer equalities in `Z[zeta_3]`**, except the explicitly independent `check_gf9.py` checker, which additionally uses numerical complex phases with tolerance `1e-8`. The enumerations and reported results in the logs have no randomized element.

**Replay from repository root:**

```sh
cc -O3 -std=c11 -Wall -Wextra -pedantic attack/exhaustive_planar.c -o /tmp/exhaustive_planar
/tmp/exhaustive_planar
python3 attack/check_gf9.py
python3 attack/binomial_probe.py
python3 attack/coulter_matthews.py
python3 attack/helleseth_kholosha.py
python3 attack/theorem_counterexample.py
python3 attack/verify_witness.py attack/minimality_witness.json
(cd attack && sha256sum -c SHA256SUMS)
```

**Output logs** are next to their corresponding programs, with the same stem and `.log` suffix. The `attack/SHA256SUMS` file gives the SHA-256 digest of **every program, log, witness and this report**, using filenames relative to `attack/`. In particular, the entire primary C enumeration is independently reproduced by a different Python backtracking implementation with matching result count and FNV-1a hash. The two substantive counterexamples have separate replay checkers or a full algebraic argument.

## Not covered

- No search enumerated all planar functions over `GF(27)`, `GF(81)`, `GF(243)`, or `GF(2187)`. The `GF(9)` enumeration is the largest **full planar-map enumeration** performed.
- No counterexample to the **universal** weak-regularity-of-planar-components question was found. The scope does not justify a universal positive claim or a claim about current open status.
- Only p-ary characteristic `3` components were Walsh-tested here; `GF(5)` and `GF(7)` were enumerated for planarity but not for their weak-regularity components. There is no search over `GF(25)` or characteristic `>=5` extensions.
- The Helleseth–Kholosha 2012 resolved conjecture's entire dual-phase formula (7) is not independently checked, nor are all families in the paper's Tables V–X.
- These computations do not prove novelty in the literature: the **parameter inconsistency is visible by comparing the paper's own theorem with Tables I/III**, and minimality results may have appeared elsewhere. No priority claim is made.
- The `m=2` zero-length endpoint is not a meaningful nondegenerate minimal-code classification. We report it solely to make the range of the printed theorem precise.

## Next steps

1. Search planar functions outside Dembowski–Ostrom/monomial classes, especially in larger `GF(3^m)` and odd-characteristic extension fields where a non-weakly-regular component might occur. State symmetry reductions before claiming exhaustivity.
2. Derive an exact canonical/isotopy classification for low-order planar maps; compare against the `GF(9)` full-function digest and expand to carefully restricted `GF(27)` classes.
3. Check every parameter expression in Sections III–IV against the signed Walsh convention and table formulas, with machine-checkable symbolic derivations and boundary conditions on `m`/`p`.
4. Investigate minimality of further code variants in Section IV, using explicit support-containment witnesses and coordinate-level independent checkers.
5. Independently verify the Helleseth–Kholosha full dual-phase equation (7) in `GF(729)`, beyond the constant phase and sign checked here.
