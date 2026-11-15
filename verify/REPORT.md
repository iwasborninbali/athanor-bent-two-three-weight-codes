# Independent verification: two-/three-weight codes from weakly regular bent functions

**Paper:** C. Tang, N. Li, Y. Qi, Z. Zhou and T. Helleseth, *Linear codes with two or three weights from weakly regular bent functions*, IEEE Trans. Inf. Theory **62** (2016), 1166–1176; [arXiv:1507.06148v3](https://arxiv.org/pdf/1507.06148), 16-page preprint. **Review date:** 2026-11-15. **Status:** independent computation, not a proof of the universal theorems.

## Scope

This audit treats an **individual statement as confirmed only on a specified finite test domain**, or where the complete proof is written below. A statement is **refuted** when a reproducible counterexample violates its literal wording. Numerical examples are checked exactly as *printed*, including their stated prime, extension degree, dimension and enumerator. `full.jsonl` records every tested defining set, code dimension, complete weight enumerator, the theoretical table prediction, a Walsh-spectrum magnitude check, and two counting methods. Six code constructions are tested for each `f`: zero, square, nonsquare and their projectively punctured counterparts.

`verify_codes.py` constructs its own `GF(p^m)` from a deterministically selected irreducible polynomial, enumerates all additive code messages, and counts trace-zero coordinates by (A) **direct shifted-index enumeration** and (B) **multiplicative cyclic FFT correlation**. An *additive* finite-field FFT separately checks bentness; coefficient signs come from exact output-value counts, not the paper's claimed sign. Direct counting is **complete for 366 constructions** with `q <= 15,625`; for **72 larger constructions**, it checks 259 predetermined beta values against the complete cyclic-FFT calculation. Largest tested field: `GF(17^4)` with `q = 83,521`. There are 28 distinct `(p,m)` pairs, 73 functions and **438** code constructions. Of the 438 table/dimension comparisons, **424 agree**; **14 fail**, all at `m = 2` due to dimension collapse / zero-weight rows. For `m >= 3`, all 378 tested constructions agree with Tables I–X. This is finite evidence, *not* an all-parameter verification.

A separate `check_lemmas.py` tests the pointwise value-count and trace-intersection formulas of Lemmas 2.7, 2.8, 2.11, 2.14 for `f(x)=Tr(a x^2)`: 34 functions, **27,540 nonzero beta** instances, complete direct intersections. An independent, standard-library-only `check_witness.py` uses `GF(9)=F_3[t]/(t^2+1)` and no code from the main verifier. `check_examples.py` checks all 18 printed numerical examples (27 family checks) and all seven examples invoking Griesmer optimality.

## Findings

### Checkable-claim inventory

| Location in arXiv v3 | Checkable statement | Status and exact coverage |
|---|---|---|
| §II.A, Lemma 2.1, p. 2; Lemmas 2.2–2.3, p. 3; Lemma 2.6, p. 4 | Cyclotomic/Gauss-sum and parity identities; sign of quadratic trace sums | **Partially confirmed** through quadratic Walsh counts/spectra at sampled primes; general number-theory identities not separately proved in this audit. |
| Definition of `RF`, pp. 2–3; Propositions 2.4–2.5, pp. 3–4 | Homogeneity, dual at zero, and dual stability | **Confirmed for tested trace-quadratic families** via completing the square and the lemma checker; **not checked** for arbitrary RF functions. |
| Lemmas 2.7–2.8, pp. 4–5 | Output-value counts of `f` and `f*` | **Confirmed** for all 34 tested trace quadratics with `3<=p<=13`, `2<=m<=7` as enumerated in `lemmas.jsonl`. |
| Lemmas 2.9–2.10, pp. 5–6; Lemmas 2.12–2.13, pp. 7–8 | Double exponential sums | **Not separately checked** by direct character-sum summation; their consequences for trace intersections are checked below. |
| Lemma 2.11, pp. 6–7; Lemma 2.14, pp. 8–10 | Pointwise counts `N_{f,beta}`, `N_{sq,beta}`, `N_{nsq,beta}` | **Confirmed** for 27,540 beta instances over 34 quadratic functions; complete direct trace-intersection test. |
| Theorem 3.1, p. 10, and Corollary 3.3, pp. 10–11 | **Printed even-`m` code lengths** | **REFUTED as written.** Missing factor `((-1)/p)^(m/2)`; proof and explicit `GF(9)` witnesses below. Weight rows in Tables I and III *do* contain this factor. |
| Tables I–IV, pp. 10–11; Theorems 3.1–3.2, Corollaries 3.3–3.4 | Zero-set code/punctured weight distributions and `[n,m]` assertions | **Confirmed for `m>=3`** in sampled fields; **refuted** for `m=2` whenever dimension drops (examples below). Correct signed even-`m` lengths match the tables. |
| Tables V–X, pp. 11–13; Theorems 4.1–4.3, Corollaries 4.4–4.6 | Square/nonsquare code and punctured weight distributions, dimensions | **Confirmed for `m>=3`** in sampled fields; at `m=2`, some length/dimension claims fail, e.g. `p=3,m=2,a=1` gives `[2,1,2]` rather than `[2,2,...]`. |
| Numerical examples after Theorems 4.1–4.3 and Corollaries 4.4–4.6, pp. 11–14 | 18 printed examples, their dimensions, distances and enumerators | **16 confirmed / 2 refuted** exactly as printed (25 of 27 family checks confirm); see next subsection. |
| Explicit Griesmer optimality assertions, pp. 13–14 | Seven `[n,k,d]` examples | **Six have a complete Griesmer exclusion of `d+1`**. For `[12,4,6]_3` it does **not**: `G_3(4,7)=12`. Optimality may still be true, but **cannot be inferred from this bound alone**. |
| §V, Table XI, p. 14; equations (5)–(10), pp. 14–15 | Weakly regular bent families and their Walsh signs | **Partial checks**: quadratics over all fields above; Helleseth–Kholosha `(6)` at `p=3,m=6`, `w^7 x^210`; Helleseth–Kholosha binomial `(8)` at `m=4`, `p=3,5,7`; Coulter–Matthews `(10)` for `p=3`, `(m,i)=(4,3),(5,3),(6,5),(7,3),(8,3)`, with both square/non-square coefficients. Dillon `(5)` and completeness of Table XI **not checked**. |
| §VI, p. 16; Introduction, p. 1 | Ding et al. planar-function open problem said to be solved; future general-function construction | **Not checkable by finite enumeration in full generality**; no new precisely formulated conjecture is posed. |
| No claim in paper | Minimal *codewords* (support minimality), separate from minimum Hamming distance | **Supplemental experiment**, not a paper theorem: exhaustive pairwise projective support tests for 204 constructions (`q<=3125`) yield 90 with all nonzero words minimal and 114 without; 234 larger/nonselected constructions not assessed. Witness pairs logged for failures. |

The 28 `(p,m)` pairs in the main program are `p=3, m=2..10`; `p=5, m=2..7`; `p=7, m=2..5`; `p=11, m=2..4`; `p=13, m=2..4`; `(17,4)`, `(19,3)`, `(23,3)`. Quadratic coefficients `1` and a primitive element are used for each. Extra nonquadratic functions are logged with their exponents. `p>=3` is always prime and tested, not arbitrary composite odd `p`.

### Literal counterexamples and typographical discrepancies

**(1) Omitted sign in Theorem 3.1 and Corollary 3.3 (p. 10).** Let `chi = (-1/p)`, the quadratic character of `-1` modulo `p`. The actual defining-set length for even `m` is

`|D_f| = p^(m-1) - 1 + epsilon * (p-1) * chi^(m/2) * p^((m-2)/2)`.

The theorem *prints* the same expression without `chi^(m/2)`. The punctured length is `|D_f|/(p-1)`, likewise missing `chi^(m/2)` in Corollary 3.3. **Complete argument:** for `f` in `RF`, `f*(0)=0`, so `W_f(0)=epsilon*(sqrt(p*))^m`, where `p*=chi*p`. For even `m`, each nonzero `y in F_p` sends `W_f(0)` to `sum_x zeta_p^(y f(x))`; Galois conjugation fixes the integer `W_f(0)`. By character orthogonality,

`#{x:f(x)=0} = (p^m + (p-1)*W_f(0))/p = p^(m-1) + epsilon*(p-1)*chi^(m/2)*p^((m-2)/2)`.

Removing `x=0` gives the corrected formula exactly; dividing by `p-1` gives the punctured version. **Standalone witness:** in `GF(9)=F_3[t]/(t^2+1)`, `f(x)=Tr(x^2)` has `epsilon=-1` and `D_f={(1,1),(1,2),(2,1),(2,2)}`, so `|D_f|=4` and punctured length 2. The *printed* expressions predict 0 and 0. `check_witness.py` independently exhausts all nine elements. At `p=3,m=6`, `f=Tr(w^7*x^210)` with `epsilon=+1`, direct calculation gives `|D_f|=224` while the printed Theorem 3.1 predicts 260; Table I correctly uses 224. These statements are **REFUTED as printed** but their evident corrected version is supported.

**(2) Unqualified `m=2` dimension claims (Theorems 3.1/4.1 and corollaries).** In the same `GF(9)`, for primitive `w=1+t`, `f(x)=Tr(w*x^2)` is weakly regular bent and homogeneous, with `epsilon=+1`, but `D_f` is **empty**, producing a length-zero, dimension-zero code, not an `[n,2]` two-weight code. For coefficient `1`, the square/non-square defining sets each have two elements but generate a code of **dimension 1**, not 2. These are structural counterexamples to the unqualified dimension assertion. For `m>=3`, no dimension failures were found in the tested domain. `full.jsonl` records 14 discrepant `m=2` constructions across the five tested primes.

**(3) Two inconsistent printed p. 12 examples:** The example after Theorem 4.2 specifies `p=3,m=5,f=Tr(wx^2)` with `[60,5,40]` and `1+24z^40+40z^48+60z^52`; the example after Theorem 4.3 repeats those field parameters with `[40,5,28]` and `1+40z^28+60z^32+24z^40`. A dimension-5 ternary code has `3^5=243` words, but **each printed enumerator sums to 125**, so both are impossible as written. Direct checks for the stated `GF(3^5)` yield `C_nsq: [90,5,54]`, `1+80z^54+72z^60+90z^66` and `C_sq: [72,5,42]`, `1+72z^42+90z^48+80z^54`. Both printed enumerators instead coincide with the corresponding `GF(5^3)` constructions (`C_nsq: [60,3,40]`; `C_sq: [40,3,28]`); thus the examples likely transpose `p` and `m` *and* misprint `k`. The other 16 distinct examples are confirmed, including the `GF(3^6)` monomial `w^7*x^210` and `GF(5^6)` examples.

**(4) Griesmer caveat:** Corollary 4.4, p. 13, says the ternary `[12,4,6]` example is optimal "due to the Griesmer bound". But `G_3(4,7)=7+3+1+1=12`, so the bound alone does **not** forbid distance 7. This audit **does not refute optimality**; a separate nonexistence theorem is required. For the other six cited optimality assertions, `G_p(k,d+1)>n` directly excludes greater distance. The claims about known optimal `[36,5,22]_3` and `[45,5,28]_3` codes are **not independently verified**.

## Evidence

The `verify/` branch includes this report, `SHA256SUMS`, and **eight individually readable evidence files**: four Python source files (`verify_codes.py`, `check_lemmas.py`, `check_examples.py`, `check_witness.py`) and their four full output logs (`full.jsonl`, `lemmas.jsonl`, `examples.jsonl`, `witnesses.jsonl`). `SHA256SUMS` contains the SHA-256 of every program, log and report. There is no dependency on Magma/Sage/paid services. The only external Python library is NumPy for main FFT/array counting; the `GF(9)` witness uses only the Python standard library.

To replay from the repository's `verify/` directory (Python 3 with NumPy installed):

```sh
cd verify
sha256sum -c SHA256SUMS
OPENBLAS_NUM_THREADS=1 python verify_codes.py > full.replay.jsonl
OPENBLAS_NUM_THREADS=1 python check_lemmas.py > lemmas.replay.jsonl
python check_examples.py full.replay.jsonl > examples.replay.jsonl
python check_witness.py > witnesses.replay.jsonl
cmp full.jsonl full.replay.jsonl
cmp lemmas.jsonl lemmas.replay.jsonl
cmp examples.jsonl examples.replay.jsonl
cmp witnesses.jsonl witnesses.replay.jsonl
```

All outputs are deterministic for the included scripts and use an explicit finite enumeration. Floating-point FFT quantities are *never* treated as exact evidence by themselves: nearest-integer outputs must have residual `<1e-6`, and a second exact trace count checks all beta (small `q`) or a fixed beta subset (large `q`). The paper was read from the publicly accessible arXiv v3 PDF. The exact PDF text was not vendored in the repository.

## Not covered

- No exhaustive proof for arbitrary weakly regular bent functions or all primes and degrees. Main tests use finitely many families, with full-beta direct checks only up to `q=15625`.
- No independent proof of every cyclotomic Gauss-sum lemma, arbitrary dual properties, the Dillon `(5)` family, or the **exhaustiveness** of the "known" families in Table XI.
- No general characterization of minimum-support/minimal codewords, secret-sharing access structures, strongly regular graph isomorphisms, association schemes, or a general-function solution to the future-work idea in §VI.
- No complete independent classification of best possible linear-code parameters. The six strict Griesmer checks are definitive for their tested examples; the `[12,4,6]_3` and two "almost optimal" claims need further coding-theory evidence.
- The connected repository's `main` branch is **never modified**. Only `verify/` contains new files. Read-back hashes on the branch are recorded in the delivery response.

## Next steps

1. Repair the three apparent typesetting defects in a formal erratum: two even-`m` lengths, and two p. 12 examples (including dimensions); explicitly state a safe lower bound on `m` for full dimension.
2. Extend the all-beta direct method beyond `q=15625`, add independent finite-field arithmetic implementations, and evaluate Dillon-type examples plus direct character-sum Lemmas 2.9–2.10/2.12–2.13.
3. Independently resolve the optimality of `[12,4,6]_3` with an explicit nonexistence certificate for `[12,4,7]_3`, and verify the two "almost optimal" reference codes.
