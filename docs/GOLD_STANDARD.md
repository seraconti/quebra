# Gold standard: published values the checks reproduce

This file is the reference for every externally published number the check battery is
validated against: the value, where it is printed, and how it was recomputed. The tests that
pin these values, and the docstrings, cite the papers directly.
`scripts/verify_gold_standard.py` recomputes every value in the recomputed columns of
sections 3 to 5 from the printed data, with numpy and scipy alone (section 7).

Tags:

- **MEASURED**: read on the printed page, by eye, or recomputed here from the printed data.
- **SIMULATED**: a seeded Monte Carlo value computed by `scripts/verify_gold_standard.py`,
  reported with its standard error.
- **UNVERIFIED**: taken from a secondary source; the primary was not opened.
- **NOT RE-CHECKED**: the source was not available when this file was checked, so its
  locators and quotations are not confirmed. The values in such rows are recomputed where
  they can be.

## 1. Sources

| key | citation | how checked |
|---|---|---|
| KL2020 | J. T. Kvaløy and B. H. Lindqvist, "A class of tests for trend in time censored recurrent event data", *Technometrics* 62(1):101-115, 2020, doi:10.1080/00401706.2019.1605936. Published online 21 Jun 2019. | journal PDF, MEASURED |
| KLv1 | Same authors and title, arXiv:1802.08339v1, 22 Feb 2018. | arXiv PDF, MEASURED |
| KLcode | github.com/jtkgithub/trendtests, files `trendtests.R`, `trendtests_multi.R`, `example_LHDdata.R`, `example_LHDdata_multi.R`. The same files are the Technometrics supplementary material, distributed as `8040713.zip`. | the zip, MEASURED; its identity with the github repository NOT RE-CHECKED |
| AH1991 | O. O. Aalen and E. Husebye, "Statistical analysis of repeated events forming renewal processes", *Statistics in Medicine* 10:1227-1240, 1991, doi:10.1002/sim.4780100806. The issue number is printed neither on the article nor in KL2020's reference list. | PDF, MEASURED |
| AD1952 | T. W. Anderson and D. A. Darling, "Asymptotic theory of certain 'goodness of fit' criteria based on stochastic processes", *Annals of Mathematical Statistics* 23:193-212, 1952. | PDF, NOT RE-CHECKED |
| AD1954 | T. W. Anderson and D. A. Darling, "A test of goodness of fit", *JASA* 49:765-769, 1954. | PDF, MEASURED |
| S1974 | M. A. Stephens, "EDF statistics for goodness of fit and some comparisons", *JASA* 69(347):730-737, 1974. | PDF, MEASURED |
| MM2004 | G. Marsaglia and J. Marsaglia, "Evaluating the Anderson-Darling distribution", *Journal of Statistical Software* 9(2):1-5, 2004, doi:10.18637/jss.v009.i02. Code archive `ADinf.c`, `AnDarl.c`. | PDF and C code, NOT RE-CHECKED |
| C2021 | S. Chatterjee, "A new coefficient of correlation", *JASA* 116(536):2009-2022, 2021, doi:10.1080/01621459.2020.1758115. | text extraction of the published PDF, not read by eye; NOT RE-CHECKED |
| GR2004 | C. Genest and B. Rémillard, "Tests of independence and randomness based on the empirical copula process", *Test* 13:335-369, 2004. | UNVERIFIED: from the references of `copula::serialIndepTest` (copula 1.1-7) |
| LEH2003 | B. H. Lindqvist, G. Elvebakk and K. Heggland, "The trend-renewal process for statistical analysis of repairable systems", *Technometrics* 45:31-44, 2003. | UNVERIFIED: from KL2020's reference list |
| KKG1989 | U. Kumar, B. Klefsjö and S. Granholm, "Reliability investigation for a fleet of load haul dump machines in a Swedish mine", *RESS* 26:341-361, 1989. | PDF, MEASURED (p. 341). KL2020's reference list (p. 115) and KLcode's `example_LHDdata.R` print volume 24. |
| KK1992 | U. Kumar and B. Klefsjö, "Reliability analysis of hydraulic systems of LHD machines using the power law process model", *RESS* 35:217-224, 1992. | PDF, MEASURED. Its Appendix, p. 224, prints the hydraulic gaps. |
| LCC2012 | J. F. Lawless, C. Çiğşar and R. J. Cook, "Testing for monotone trend in recurrent event processes", *Technometrics* 54:147-158, 2012. | UNVERIFIED: from the KLcode header and KL2020's reference list (p. 115) |

## 2. KL2020 against KLv1: the locator map

The two versions are not interchangeable. Section numbers, the second example and the
multi-process weights all differ. Cite KL2020 and use its locators; use KLv1 locators only
for material KLv1 alone contains.

| content | KLv1 | KL2020 | same? |
|---|---|---|---|
| eq (4) LR, single process | p. 4 | p. 103 | yes |
| eq (5) KS | p. 5 | p. 103 | yes |
| eq (6) CvM | p. 5 | p. 103 | yes, including the typo below |
| eq (7) AD | p. 5 | p. 104 | yes |
| eq (9) ELR | p. 6 | p. 104 | yes |
| Section 3.6, gamma_hat from the "sample standard deviation" of complete gaps | p. 6 | p. 104 | yes |
| eq (10) mu_tilde, sigma_tilde, gamma_tilde | p. 6 | p. 104 | yes |
| eq (11) sigma*, and the authors declining to use it | p. 7 | p. 104 | yes |
| eq (13) LR^m | p. 8 | p. 105 | yes (KL2020 writes T_ij/tau_j) |
| **eqs (14)-(16), multi-process weights** | **A_j ∝ gamma_j tau_j sqrt(N_j)**, p. 8 | **A_j ∝ sqrt(N_j)/gamma_j**, p. 105 | **NO: different tests** |
| Section 4.2 | "Other Tests for m Processes", pp. 8-9 | "Further Tests for m Processes", pp. 105-106 | text yes, title no |
| AD normal approximation "less well ... very skew distribution" | Section 4.2, p. 9 | Section 4.2, p. 106 | yes |
| linear rank and generalised Laplace tests | GL only, Section 4.1, p. 8 | Sections 5.1-5.2, p. 106 | no |
| simulation study, level properties | Section 5.1, p. 10 | Section 6.1, p. 107 | wording differs |
| Figure 1, level vs expected number of events, Weibull shapes 0.75 and 1.5 | p. 10 | p. 107 | |
| AD left out for m > 1 "as the Cramer-von Mises test had better level properties" | Section 5, p. 9 | Section 6, p. 106 | yes |
| TRP asymptotics, optimality of the eq (14) weights | none | Section 7, 7.1.2 | KL2020 only |
| single LHD machine | Section 6.1, pp. 14-15 | **Section 8.1**, pp. 111-113 | numbers identical; KL2020 also prints LinR |
| second example | Section 6.2, small bowel motility, m = 19, Table 4 p. 17 | **Section 8.2, hydraulic LHD, m = 6**, Table 4 p. 113 | **different data** |
| permutation remark (cites LCC2012 for its validity; "confirmed this in simulations not reported here") | Section 7, p. 17 | **Section 9**, p. 113 | yes |
| sigma_hat^2 written with 1/N(tau) | Appendix 1, p. 18 | Appendix A.2, p. 114 | yes |
| R code | "can be obtained from the authors" | github.com/jtkgithub/trendtests, p. 113 | |

**eq (6) typo, both versions.** The middle term prints `- i N(tau) (T_{i+1}^2 - T_i^2)/tau`;
integrating the definition gives `/tau^2`, and KLcode's `CvMtestobs` divides by `tau^2`.
`tests/test_checks_cvm.py::test_the_bracket_is_the_bridge_integral` settles it. MEASURED.

## 3. Single process: the load-haul-dump record

Data: KL2020 Table 1, p. 112 (KLv1 Table 1, p. 14). 36 failure times, the last at 1970 h.
Table 1 note: "The data are time censored at 2000 hr." Section 8.1 (p. 111) adds that this
is the authors' choice: "For the purpose of this example we considered the data to be time
censored at tau = 2000 hrs." KLv1 has the same sentence ending in "hours". The 36 gaps
differenced from Table 1 are the
list in `tests/fixtures/load_haul_dump.py`, sum 1970, censored gap 30. MEASURED.

| quantity | printed | locator (KL2020) | from | recomputed, 1/(N-1) | recomputed, 1/N |
|---|---|---|---|---|---|
| mu_hat | 54.72 | Table 2 row 1 | Section 3.6 | 54.722 | same |
| sigma_hat | 48.61 | Table 2 row 1 | Section 3.6 | **48.611** | 47.931 |
| gamma_hat | 0.888 | Table 2 row 1 | Section 3.6 | **0.8883** | 0.8759 |
| mu_tilde | 55.56 | Table 2 row 2 | eq (10) | 55.556 | n/a |
| sigma_tilde | 47.23 | Table 2 row 2 | eq (10) | 47.228 | n/a |
| gamma_tilde | 0.850 | Table 2 row 2 | eq (10) | 0.8501 | n/a |
| Weibull mu, sigma, gamma | 55.46, 47.22, 0.851 | Table 2 row 3 | MLE with the 30 h gap right-censored | 55.459, 47.219, 0.8514 (shape 1.179) | n/a |
| Laplace | 0.605 | Section 8.1 text | eq (4) without 1/gamma_hat (Section 3.1) | 0.6051 | same |
| LR | 0.681 (p 0.50) | Section 8.1 text; Table 3 | eq (4) | **0.6811** (p 0.4958) | 0.6908 |
| sigma* | 42.77 | Section 8.1 text | eq (11) | 42.771 | n/a |
| gamma* = sigma*/mu_hat | 0.782 | Section 8.1 text | eq (11) | 0.7816 | n/a |
| LR with gamma* | 0.774 (p 0.44) | Section 8.1 text | eq (4) | 0.7741 (p 0.4389) | n/a |

The Laplace, LR, sigma*, gamma* and LR-with-gamma* values are in the Section 8.1 TEXT, not in
Table 2. All MEASURED.

Table 3, p. 113, p-values. The table's note states no tail; the Section 8.1 text (p. 112)
calls them "Two-sided p-values for all tests". LR and ELR reproduce as two-sided normal
p-values, KS, CvM and AD as upper tails of their limits. KLv1 Table 3, p. 15, lacks the LinR
column.

| test | printed | recomputed, 1/(N-1) | recomputed, 1/N | statistic, 1/(N-1) |
|---|---|---|---|---|
| LR | 0.50 | 0.4958 | 0.4897 | 0.6811 |
| KS | 0.29 | 0.2864 | 0.2711 | 0.9850 |
| CvM | 0.13 | 0.1312 | 0.1240 | 0.30462 |
| AD | 0.086 | 0.0856 | 0.0795 | 2.05555 |
| ELR, a = 1/2 | 0.011 | 0.0115 | 0.0103 | 2.5283 |
| LinR | 0.76 | not recomputed | | |

Of the statistics in the last column only LR is printed (0.681, Section 8.1 text). The
others come from KLcode's closed forms, which agree with direct numerical integration of the
tied-down bridge to 1e-14, so they are cross-implementation values, not published ones.
p-values use the exact limits for CvM and AD, the Kolmogorov limit for KS and the normal for
LR and ELR. The AD p at 1/N is the exact limit's; `adinf` gives 0.0796.

### The divisor

- Section 3.6, p. 104: "obvious choices for estimators of mu and sigma are the sample mean
  mu_hat and sample standard deviation sigma_hat of the completely observed interevent times."
  No formula.
- Appendix A.2, p. 114 (KLv1 Appendix 1, p. 18): `sigma_hat^2 = (1/N(tau)) sum (X_i - mu_hat)^2`,
  inside a consistency argument where the divisor does not matter asymptotically.
- Every printed number (Table 2, Section 8.1, Table 3, Table 4, and KLv1 Section 6.2) lands
  on 1/(N-1), AD^m and CvM^m within the authors' Monte Carlo error, and misses on 1/N. KLcode's `findCV(sigma = "s")` computes `sqrt(var(xvec))`,
  and R's `var` divides by N-1.

So both divisors have a source, and only one matches the paper's numbers. The shipped
default is 1/(N-1) (`GAMMA_DEFAULT = GAMMA_COMPLETE_SAMPLE`). MEASURED.

## 4. Multiple processes

### 4.1 Hydraulic systems of six LHD machines, KL2020 Section 8.2

Data: KK1992 Appendix, p. 224, machines 1, 3, 9, 11, 17, 20, identical to KLcode
`example_LHDdata_multi.R`. Not printed in KL2020. Convention, p. 112: "for the purpose of this example we will follow their analysis
and assume that for each machine the data are time censored at the last recorded event
time", so each list's last gap ends at tau_j and is not an event.

| machine | N_j | tau_j | gamma_hat_j, 1/(N-1) | LR_j |
|---|---|---|---|---|
| 1 | 22 | 2496 | 1.1952 | 1.6538 |
| 3 | 24 | 3526 | 1.0864 | 0.6848 |
| 9 | 26 | 4743 | 1.1075 | 1.9106 |
| 11 | 27 | 2913 | 0.9277 | 0.2606 |
| 17 | 25 | 3230 | 1.0102 | 1.5297 |
| 20 | 22 | 3309 | 0.8672 | 0.1102 |

Table 4, p. 113, p-values (no tail stated; LR, ELR and GL reproduce as two-sided normal
p-values, AD and CvM are upper tails of simulated nulls):

| test | printed | recomputed | condition |
|---|---|---|---|
| LR^m | **0.019** | **0.0188** (LR^m 2.349) | KL2020 weights, 1/(N-1) |
| | | 0.0164 (2.399) | KL2020 weights, 1/N |
| | | 0.0073 (2.683) | KLv1 weights, 1/(N-1) |
| | | 0.0062 (2.739) | KLv1 weights, 1/N: `c1.statistic(..., weights=WEIGHTS_ARXIV_V1, gamma_estimator=GAMMA_COMPLETE)` |
| GL | 0.062 | 0.0623 (GL 1.864) | no gamma, no weights: checks the transcription, but one p to 3 decimals cannot catch a single 1 h error in a gap |
| AD^m | 0.076 | 0.0752 exact; 0.0750 ± 0.0004 at 4e5 draws (SIMULATED) | weights ∝ tau_j; exact by numerical inversion of the limit, draws from the limit |
| CvM^m | 0.064 | 0.0651 exact; 0.0650 ± 0.0004 at 4e5 draws (SIMULATED) | as AD^m; KLcode's example uses 1e4 draws, s.e. about 0.0025 |
| ELR^m | 0.003 | 0.0027 | `avec = 1 - 0.5 tau_j / max tau` as in KLcode |
| | | 0.0079 | `a_j tau_j = 4743/2` as the Section 8.2 text describes |

Only KL2020 weights with 1/(N-1) reproduce LR^m. The eq (10) and eq (11) estimators do not
either: 0.0140 with eq (10)'s gamma in both A_j and LR_j (KLcode's `sigma = "c"` keeps the
sample CV in A_j and gives 0.0129), and 0.0065 with eq (11)'s. The printed AD^m and CvM^m do not round from the exact values but
lie within one s.e. of a 1e4-draw simulation, the count KLcode's example uses; the paper
states none. The ELR^m value matches the code, not the text, so it is not used as a fixture.
MEASURED.

### 4.2 Small bowel motility, KLv1 Section 6.2 only

Data: AH1991 Table I, p. 1229. 19 subjects, 80 complete periods (minutes), one censored
period per subject; each record starts at a phase III. Transcribed in
`tests/fixtures/load_haul_dump.py` and checked against the page image. MEASURED.

The paper tests the stronger null that all 19 processes share one gap distribution and so
pools a single gamma_hat over the 80 complete periods (p. 15).

| quantity | printed (KLv1) | recomputed, 1/(N-1) | 1/N |
|---|---|---|---|
| mu_hat | 98.76 (p. 16) | 98.76 | same |
| sigma_hat | 52.62 | 52.62 | 52.29 |
| gamma_hat | 0.533 | 0.5328 | 0.5294 |
| Laplace^m, KLv1 eq (16) with gamma = 1 | 1.95, p = 0.051 | 1.9530, p 0.0508 | same |
| LR^m = Laplace^m / gamma_hat | 3.67, p = 0.00024 | 3.6658, p 0.000247 | 3.6890 |
| GL | p = 0.007 (Table 4, p. 17) | 2.685, p 0.0073 | same |
| Weibull mu, sigma, gamma | 104.49, 52.45, 0.502 | 104.50, 52.46, 0.502 (shape 2.092) | |

The printed p = 0.00024 equals 2 Phi(-3.67) of the rounded statistic. The unrounded 3.6658
gives 0.000247, which truncates, but does not round, to 0.00024, so the route is not
established. KLv1 Table 4 prints LR p < 0.0001, which contradicts its own text; the text value
is the fixture.

Shipped C1 cannot reproduce 3.67. It estimates gamma per segment and refuses a segment with
fewer than two complete gaps; subject 5 has one. Without subject 5, per-segment gamma gives
LR^m = 4.49 under the shipped defaults (KL2020 weights, 1/(N-1)) and 4.37 under KLv1 weights
with 1/N. A pooled-gamma estimator would be new functionality.

## 5. Limiting distributions

### 5.1 Anderson-Darling A^2

| level | S1974 Table 1A part 1.0, p. 732 | exact limit | upper tail at S1974 value |
|---|---|---|---|
| 15% | 1.610 | 1.6212385 | 0.1523 |
| 10% | 1.933 | 1.9329578 | 0.1000 |
| 5% | 2.492 | 2.4923672 | 0.0500 |
| 2.5% | 3.070 | **3.0774642** | 0.0252 |
| 1% | 3.857 | **3.8781250** | 0.0102 |

"Exact limit" is the limiting law evaluated by numerical inversion of its characteristic
function (AD1954 eq (7), p. 768), checked against the series of AD1954 eq (8), p. 768, to
2e-12. The same series is cited elsewhere as AD1952 eq (4.38), p. 204: NOT RE-CHECKED. MM2004
p. 2 is quoted as saying the published 99th percentile "3.857 is actually 3.878125..." and as
giving 1.9329578327415937304, 2.4923671600494096176 and 3.8781250216053948842; the digits agree
with the exact limit to 1e-11, but the page is NOT RE-CHECKED. No published source was located
for the exact 2.5% point. AD1954 p. 766 prints the asymptotic points 1.933, 2.492 and 3.857
at .10, .05 and .01, so the 3.857 that MM2004 is quoted as correcting is AD1954's. S1974
prints its row "For all n ≥ 5"; KL2020 p. 104 gives only the 5% value 2.492, citing AD1954.
MEASURED.

**`adinf` accuracy.** MM2004 p. 3 is quoted as stating |error| < 2e-6 for 0 < z < 2 and
< 8e-7 for z ≥ 2 (NOT RE-CHECKED). Against the exact limit, the shipped `adinf` errs by
-1.9e-5 at z = 1.0, -1.1e-5 at z = 1.933, +8.1e-6 at 2.492, -2.6e-6 at 3.878 and +8.4e-6 at
z = 5; the largest errors are 1.95e-5 at z = 0.97 below z = 2 and 9.1e-6 at z = 2.59 above
it. That is immaterial at alpha = 0.05, but about ten times the quoted accuracy. Whether the
shipped coefficients reproduce MM2004's own `ADinf.c`, or differ from it, needs the C file,
which is NOT RE-CHECKED. MEASURED against the exact limit.

### 5.2 Cramér-von Mises W^2

| level | AD1952 Table 1, p. 203 (NOT RE-CHECKED) | S1974 (modified W^2) | series (the Bessel-K form; cited as AD1952 eq (4.35), p. 202: NOT RE-CHECKED) |
|---|---|---|---|
| 10% | 0.34730 (AD1954 p. 766: .3473) | 0.347 | 0.3473049 |
| 5% | 0.46136 | 0.461 | 0.4613613 |
| 2.5% | **not printed** | 0.581 | 0.5806147 |
| 1% | 0.74346 | 0.743 | 0.7434593 |
| 0.1% | 1.16786 | | 1.1678583 |

AD1952 Table 1 is described as tabulating z at a1(z) = .01 to .99 in steps of .01, plus
.999, with no .975 row, so that 0.58061 has no published source among these (NOT
RE-CHECKED). The series values are MEASURED: every AD1952 value above is the series rounded to
five decimals, and 0.58061 is correct to five decimals.

## 6. Chatterjee's xi

C2021 was read by text extraction of the published PDF, not by eye, and is NOT RE-CHECKED.

- Tie-free form: eq (1), p. 2010, `1 - 3 sum|r_{i+1} - r_i| / (n^2 - 1)`.
- Form with ties: an unnumbered display on p. 2010,
  `1 - n sum|r_{i+1} - r_i| / (2 sum l_i (n - l_i))`, where r_i is the number of j with
  Y_(j) ≤ Y_(i) and l_i the number with Y_(j) ≥ Y_(i).
- X-ties are broken "uniformly at random" (p. 2010).
- Theorem 2.1, p. 2011: X, Y independent and Y continuous give sqrt(n) xi -> N(0, 2/5).
  The normal approximation is "roughly valid even for n as small as 20".
- Theorem 2.2: X, Y independent with Y arbitrary give N(0, tau^2), tau^2 from formula (3).

The repository calls the form with ties the tie-corrected form and cites it to C2021.
`.claude/qre_checks_reference.tex` numbers it eq (8), which is that document's own
numbering and is unrelated to KL2020's eq (8), the extended Lewis-Robinson integral.

## 7. Reproducing

```
python scripts/verify_gold_standard.py --klcode-zip <8040713.zip> --fixture tests/fixtures/load_haul_dump.py
```

Both arguments are optional. `--klcode-zip` is the KL2020 supplementary-material archive
(doi:10.1080/00401706.2019.1605936), distributed as `8040713.zip`; `--fixture` is the test
fixture that carries the three records.

The script transcribes its data from the printed pages (KL2020 Table 1, KK1992 p. 224, AH1991
Table I). It compares the hydraulic record with the authors' R code and all three records with
the fixture, recomputes every value in the recomputed columns of sections 3 to 5 under both
divisors, and prints each beside the value this file gives and the published one. Those values
are held in the script, and it checks that each one still occurs in this file, so an edit here
that the script does not share fails the run. It exits 1 if any value differs from its
recomputation by more than half a unit in its last printed digit.

Prose figures the script does not recompute: the `adinf` p 0.0796 in section 3, the eq (10)
and eq (11) values 0.0140, 0.0129 and 0.0065 in section 4.1, and the agreement tolerances
1e-14 (section 3) and 1e-11 (section 5.1).
