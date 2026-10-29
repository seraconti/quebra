# P5 review findings

Appended at the moment of discovery, not batched. Format:
`SEVERITY | file:line | problem | minimal fix`

Severities: CRITICAL (wrong number reaches an artifact), IMPORTANT (wrong claim, or a
defect that will fire), MINOR (accuracy, naming, dead reference).

Authority for "what it should be": qre_checks_reference_PREVIEW.pdf, 15 pp.

---

## Increment A - audit

IMPORTANT | analyzers/shape_stats.py:85 | `chatterjee_xi` computes `1 - 3*sum|dr|/(n^2-1)`, which is the TIE-FREE reduction. Reference eq (8) is `1 - n*sum|dr| / (2*sum l_i(n-l_i))`. The two agree only when Y has no ties; our durations are lattice-valued and heavily tied at loose thresholds. | Implement eq (8). Increment B.

IMPORTANT | analyzers/shape_stats.py:104-112 | `xi_p_value` divides by `XI_NULL_VARIANCE = 2/5` unconditionally. That null is derived for continuous tie-free data; reference 10.3 caveat 1 says under ties "the closed-form p-value must not be used and a permutation calibration should replace it". | Permutation-calibrate. Increment B, same increment as the estimator - fixing one without the other is worse than fixing neither.

MINOR | analyzers/shape_stats.py:82 | Only constant X is guarded. Under eq (8) a constant Y makes `l_i = n` for all i, so the denominator `2*sum l_i(n-l_i)` is exactly 0. Harmless today, division by zero the moment eq (8) lands. | Add the constant-Y guard with eq (8).

MINOR | analyzers/checks/c1_lewis_robinson.py:26-30, c2_anderson_darling.py:28-34 | Both docstrings argue for the permutation calibration on their own authority. It is the source paper's OWN prescription (its Section 7), citing Lawless, Cigsar and Cook, Technometrics 54, 2012, for validity under time censoring - reference Section 4. A claim with no locator, in a repo whose rule is that every doc claim match a source. | Add the citation to both docstrings.

MINOR | docs/iid_checks/C5_rank_autocorr.md:16 | Says the studentization "follows the Chung-Romano / Romano-Tirlea argument". Reference Section 7 states explicitly that it is NOT Romano-Tirlea studentization: ours is a lag-comparability correction (different lags have different pair counts), theirs restores level under an uncorrelated-but-dependent null. My own misattribution, written in P4b. | Reword to lag-comparability and name the distinction.

IMPORTANT | analyzers/check_ledger.py:409 | Catches only `(ValueError, KeyError)`. `c3._invoke_rscript` raises `RuntimeError` in five places (c3_serial_copula.py:92,118,123,130,135) and `subprocess.run(timeout=120)` raises `TimeoutExpired`. R is now installed, so `rscript_path()` is non-None and this call executes for real: a failed R invocation KILLS the ledger job instead of writing a `not computed` row. | Widen the except. Increment B.

MINOR | analyzers/check_ledger.py:271 | `r_version` is assigned `str(r_path)` - the PATH to Rscript, not a version - while the module docstring (:20-21) says "the R version is discovered at runtime and lives in CheckLedger.r_version". Field name and docstring both describe something the code does not store. | Capture `Rscript --version`, or rename the field to `r_path`.

IMPORTANT | core/job.py:293-304 with panels/_artifact_guard.py:36 | `Job.code_hash` hashes ONLY the job's own source file, and `StaleArtifactGuard` compares field NAMES only. So editing any analyzer changes every downstream number while the run `identity` stays byte-identical and old pickles keep loading. Nothing in the pipeline signals an analyzer edit. | Out of P5's scope to fix, and larger than it looks. Recorded because it is why the xi change must have its test written first, and it belongs in the Part 7 list.

MINOR | reference PDF Section 5.2 eq (7) | The CvM middle term prints `- i N (T^2_{i+1} - T^2_i)/tau`, which is dimensionally inconsistent: the first term `i^2 X_{i+1}/tau` is dimensionless. Deriving `integral_0^1 (N(s tau) - sN)^2 ds` gives `/tau^2`, and reproduces the document's own tail `N^2[u_N^2 - u_N + 1/3]` exactly. Either a typo or a dropped superscript in extraction. | Correct the reference. Derive rather than transcribe in Increment B.

GREEN | analyzers/checks/c1_lewis_robinson.py:72-88 | Eq (16) transcribed correctly, and only eq (16) exists - the m=1 collapse to eq (4) is asserted numerically by test rather than by a second code path. | none
GREEN | analyzers/checks/c2_anderson_darling.py:69-118 | `ad_limiting_cdf` is Marsaglia-Marsaglia `adinf` alone; `errfix` is deliberately absent and the docstring says why (the statistic is divided by an estimated gamma^2, so the finite-N null is not the finite-N AD null). Matches reference Section 4. | none
GREEN | analyzers/checks/_permutation.py + _rank_serial.py | Permutation never crosses a segment boundary: measured over 500 draws at sizes [4,5,3], every index stays in its own block. Lag pairs likewise - lag-1 gives 9 pairs (sum of size-1), not 11 (total-1). | none
GREEN | analyzers/shape_stats.py:167-181 | `dcor` is the Szekely-Rizzo-Bakirov (2007) double-centred V-statistic and returns R, not R^2 - the same convention as `energy::dcor`, so the Tier 4 comparison is direct. | none
GREEN | analyzers/check_ledger.py LEDGER_COLUMNS | Reference Section 6 requires the ledger report C3's dropped-censored count; `n_censored_dropped` is a ledger column and C3 populates it. | none

IMPORTANT | analyzers/checks/ (m > 1 path) | The multi-process eq (16) route has NO external validation. The only published worked example is Kvaloy and Lindqvist Section 6.2 (Aalen and Husebye 1991 bowel motility, m = 19), whose raw data is in neither the paper nor any repository or R package located on 2026-08-12. Everything we have at m > 1 is internal consistency, which cannot catch a shared transcription error. | `test_small_bowel_motility_multiprocess` skips with the targets and an activation recipe. Obtain the data, or accept the gap explicitly in the instrument report.

MINOR | tests/ | The `m > 1` Anderson-Darling route is an UNWEIGHTED SUM of per-segment eq (7) values, calibrated by permutation. That is an extension beyond the source: the paper drops AD for m > 1 in favour of Cramer-von Mises and offers no multi-process AD form. Documented in c2's docstring, but carries no published evidence of any kind. | Cramer-von Mises (Increment B) gives the first comparator.

---

## Part 1 fidelity table

| routine | claimed source (eq) | what the code computes | agree? | note |
|---|---|---|---|---|
| C1 `statistic` | K-L eq (16); eq (4) at m=1 | `sqrt(12)*sum_j[sum_i T_ij - (N_j/2)tau_j] / sqrt(sum_k g_k^2 tau_k^2 N_k)` | YES | only eq (16) exists; the m=1 collapse is asserted numerically, not by a second code path. Now also pinned to published Table 2 values. |
| C1 calibration | K-L Section 7; Lawless, Cigsar and Cook 2012 | asymptotic N(0,1) + within-segment permutation | YES | but the docstring argues it on its own authority and carries NO locator. Finding logged. |
| C2 `_eq7` | K-L eq (7) | matches term for term; equals the textbook `A^2` at gamma=1 to 1e-15 | YES | pinned by `test_eq7_is_the_classical_anderson_darling` at five n. |
| C2 `ad_limiting_cdf` | Marsaglia and Marsaglia JSS 9(2) 2004 | `adinf` ALONE; `errfix` absent by design | YES | reference Section 4 agrees this is correct: the statistic is divided by an estimated gamma^2, so the finite-N null is not the finite-N AD null. Four critical values pinned. |
| gamma estimators | K-L Appendix 1, eqs (10), (11) | 1: population (1/N) complete gaps. 2: eq (10). 3: ABSENT. | PARTLY | estimator 3's absence matches reference Section 2.3 exactly. But the DIVISOR is ours, not the paper's: reference 2.3 calls estimator 1 the "sample mean and standard deviation", which reads as 1/(N-1), and all three of Table 2's complete-gap numbers land on 1/(N-1) together. See the new findings below. |
| C3 | Ghoudi-Kulperger-Remillard 2001; Genest-Remillard 2004; Kojadinovic-Yan 2011 | imports `copula::serialIndepTest` over a file bridge | N/A | imported, not transcribed. Never executed until now. Censored windows dropped, and the ledger does report the dropped count as reference Section 6 requires. |
| C5 | rank autocorrelation + permutation; reference Section 7 | max over lags of lag-standardised `r_h`, within-segment permutation | YES | permutation and lag pairs both verified never to cross a segment boundary. The "studentized" naming is lag-comparability, NOT Romano-Tirlea - our own doc misattributes it. |
| C6 | portmanteau `Q = sum_h m_h r_h^2`, exact permutation | matches | YES | the one check whose validity rests on no approximation and no citation. |
| Spearman `_spearman` | classical | `scipy.stats.spearmanr` equivalent on ranks | YES | no test of its own. |
| xi `chatterjee_xi` | Chatterjee 2021 eq; reference eq (8) | `1 - 3*sum|dr|/(n^2-1)` = the TIE-FREE reduction | **NO** | reference eq (8) is the tie-corrected form. Our durations are heavily tied. Increment B. |
| xi `xi_p_value` | `sqrt(n) xi -> N(0, 2/5)`, tie-free | uses 2/5 unconditionally | **NO** | reference 10.3 caveat 1 forbids the closed form under ties. Increment B. |
| dcor | Szekely, Rizzo and Bakirov 2007 | double-centred V-statistic, returns R | YES | same convention as `energy::dcor`, so the Tier 4 comparison is direct. No test of its own. |

---

## Part 7 - ranked test proposals (deliverable; NOT implemented in this pass)

Ranked by probability of a wrong number reaching an artifact, times its severity.

1. **Golden-value pins on the shape statistics.** Nothing in `tests/` touches `shape_stats`.
   Catches: exactly the xi tie-free error found in this pass, which has been shipping into
   every `DistinguishBand` since P3; a dcor double-centring error; a `peak_frac` off-by-one.
   This is rank 1 because the defect is already confirmed present and was invisible.
2. **A pin that fails when an analyzer changes a panel number.** `Job.code_hash` hashes only
   the job file and `StaleArtifactGuard` compares field names, so an analyzer edit changes
   every downstream number silently. Catches: any future edit to `shape_stats`, the bands or
   `_within_calibration_compute` that moves a figure without moving the identity.
3. **`segments_from_windows` against a REAL gapped record.** Currently exercised only on
   synthetic tables. Catches: a wrong interior `tau`, a mis-placed split, the `gap_resume`
   taxonomy diverging from `gap_spans_s`.
4. **Cross-check the m > 1 AD sum against CvM on the same segments.** The AD multi-process
   route is an extension beyond the source with no published evidence. Catches: a shared
   error in the per-segment sum that both self-consistency and the bench would miss.
5. **Degenerate-input matrix, one test per case.** Every duration identical; a segment with
   exactly `min_events`; `tau` within float noise of `T_N`; a single-read window with
   `duration_s == 0`. Partly covered by `TAU_MARGIN` work, never systematically.
6. **Permutation-set provenance.** `check_permuted` guards SHAPE only; a same-shape matrix
   built from a different record is undetectable. Catches: silent corruption of every
   permutation p-value in a run, the exact failure the shared-permutation optimisation risks.
7. **C3 CSV round trip with R present.** The write/parse path in `_invoke_rscript` has never
   run. Catches: a locale decimal separator, a column-name mismatch, an R-side error that
   exits 0.
8. **`dcor_with_reason`'s four reason codes.** Catches: a NaN attributed to "constant" when
   the real cause was non-finite input - a wrong reason in a typed contract.
9. **`window_stats` completeness.** Catches: a key silently dropped from the shape row, which
   the panel would render as an empty annotation rather than an error.
10. **Ledger verdict precedence under combined failures** (ties AND too few events AND
    bench-rejected simultaneously). Catches: a precedence inversion that reports the least
    important reason.

---

## Increment A - reviewer pass (independent, fresh context, 2026-08-13)

Scope reviewed: tests/fixtures/__init__.py, tests/fixtures/load_haul_dump.py,
tests/test_checks_published_values.py, this ledger. All gitignored; read directly.

DATA VERDICT: the 36 cumulative times are correct. Strictly increasing, min gap 3 h,
T_36 = 1970, mu = 1970/36 = 54.7222, and SEVEN independent published quantities reproduce
to printed precision (54.72, 47.23, 0.850, 0.605, 0.774, 48.61 via Bessel, 0.888, 0.681).
No transcription error in the fixture.

IMPORTANT | tests/test_checks_published_values.py:22-28 and this ledger:53 | "their Appendix 1 DEFINES the population (1/N) form" has no locator and the cited authority points the other way: reference 2.3 estimator 1 is "sample mean and standard deviation of the completely observed gaps", and ALL THREE of Table 2's complete-gap numbers (48.61, 0.888, 0.681) land on 1/(N-1) simultaneously. The parsimonious reading is that GAMMA_COMPLETE's 1/N is OUR deliberate choice (_multiprocess.gamma_hat:80-81 says exactly that: to match eq 10's divisor), not the paper contradicting its own appendix. Asymptotic C1 is then inflated by sqrt(N/(N-1)): 1.4% at N=36, 12% at N=5, anti-conservative, unflagged. | Add the Appendix-1 locator, or restate as our divisor choice.

IMPORTANT | analyzers/checks/_multiprocess.py:54 vs :57 | Signature default is `estimator=GAMMA_TRUNCATED`; the docstring's first line says "`GAMMA_COMPLETE` is the DEFAULT". Both live call sites (c1:83, c2:171) pass explicitly, so nothing is wrong today; the next caller that omits it silently gets eq (10), which the same docstring says goes negative at small N. The fidelity table's gamma row does not record this. | Fix the signature or the sentence; note it in the row.

IMPORTANT | this ledger:49 | "Now also pinned to published Table 2 values" overstates. Pinned: `c1.statistic`, m=1, GAMMA_COMPLETE only. Not pinned: m>1 eq (16) (the gap already logged at :39), `statistic_batch`, `run`'s p-values, GAMMA_TRUNCATED inside C1. And the directly asserted 0.691 is NOT a published value - the published 0.681 is reached only through the sqrt((N-1)/N) bridge. | Qualify the cell to "m=1, complete-gaps, statistic() only".

MINOR | tests/test_checks_published_values.py:96-101 | "a transcription error in the sqrt(12)/(tau*sqrt(N)) factor would show up here and nowhere else in this file" is false - line 75 and the eq (11) test carry the same factor. "computed independently" also misleads: gamma_hat cancels exactly. The test is NOT circular (the cancellation is what leaves a clean pin on the numerator) but it is nearly implied by the other two pins (0.876 x 0.691 = 0.6053, inside the 0.001 band). | Reword to "pins the numerator and its scaling; gamma cancels by construction".

MINOR | tests/test_checks_published_values.py:88-91 | `gamma_tilde * tau / N` IS reference 2.3's definition (mu_tilde = tau/N), so it is not question-begging - but given the pin at :63 it is the same assertion rescaled (47.23/(2000/36) = 0.85010). Its only independent content is that the paper's printed pair is mutually consistent and that TAU_H and N are right. | Say that, rather than "genuine check ... weaker than a direct pin".

MINOR | tests/test_checks_published_values.py:112-125 | The second assertion is arithmetic on the first; at one N, sqrt(N/(N-1)) = 1.0142 is not separable from any other factor within the 0.01 band, so "identifies the difference as the divisor and nothing else" overstates. What actually identifies it is the simultaneous agreement of 48.61, 0.888 and 0.681 - which the test does assert. | Reword the reasoning; the conclusion stands.

MINOR | tests/test_checks_published_values.py:176,198 | `request` is unused and the skip is unconditional, so the "TO ACTIVATE" recipe activates nothing: adding `SMALL_BOWEL_SEGMENTS` leaves the test skipping until someone rewrites the body. `SMALL_BOWEL_PUBLISHED` (fixture:73-82) is never imported and its numbers are duplicated verbatim in the test docstring. | `seg = getattr(fx, "SMALL_BOWEL_SEGMENTS", None); if seg is None: pytest.skip(...)`, assert against SMALL_BOWEL_PUBLISHED, drop `request`.

MINOR | tests/fixtures/load_haul_dump.py:3-5 | "Technometrics 62(1):101-115" and "Table 1" are unverifiable offline and the volume/issue/pages look wrong for this paper. A thesis citation hardens from here. | Check against the journal record.

MINOR | tests/test_checks_published_values.py, tolerances | Measured margins (used/allowed): mu 0.0022/0.01, sigma_tilde 0.0022/0.01, gamma_tilde 1.0e-4/1e-3, laplace 6.3e-5/1e-3, gamma_pop 1.1e-5/1e-3, LR 2.1e-4/1e-3, sigma bridge 8e-4/0.01, gamma bridge 3.2e-4/1e-3, LR bridge 1.3e-4/1e-3, eq11 1.4e-4/1e-3. Each is ~2x the print half-width - correct order. Discrimination measured by perturbation: a +/-10 h error at any single failure time fails at least one assertion; a +/-1 h error is caught at indices 10, 17, 34 but NOT at 0, 5, 25. | No change; state the +/-1 h residual as the bound this tier certifies.

GREEN | tests/test_checks_published_values.py:139-160 | `_sigma_star_squared` matches reference 2.3 estimator 3 term for term, and reproducing 0.774 is the file's only pin on the gap ORDER (eq 11 is order sensitive, gamma is not). The docstring undersells it as "the surrounding machinery". | none
GREEN | claims verified | `pytest tests/` = 196 passed, 2 skipped. Increment A touched nothing outside tests/ and .claude/ (core/paths.py and .gitignore differ from HEAD but predate it by a day). Permutation confined to blocks, reproduced at sizes [4,5,3] over 500 draws. eq (11) unreachable: GAMMA_ESTIMATORS is exactly the two labels and the ValueError text matches. NOTE: a bare `pytest` from the repo root dies with INTERNALERROR importing monoliths/v13fig/test_carve.py (SystemExit at import); the green claim holds only when scoped to tests/. | none

---

## Post-review corrections (reviewer pass, 2026-08-13)

IMPORTANT | analyzers/checks/_multiprocess.py:54 | Signature is `estimator: str = GAMMA_TRUNCATED`, while the docstring's first line says "`GAMMA_COMPLETE` is the DEFAULT". Both current call sites pass explicitly (c1:83, c2:171) so no shipped number is wrong today - but the next caller that omits the argument silently gets eq (10), which the same docstring warns goes negative at small N. | Change the signature default to `GAMMA_COMPLETE`. Increment B, where source may be touched.

IMPORTANT | tests/test_checks_published_values.py (docstring, now corrected) | An earlier draft asserted "their Appendix 1 DEFINES the population (1/N) form while Table 2 uses 1/(N-1)". That claim came from the P5 brief and is supported by NO source in hand: reference Section 2.3 calls the paper's estimator 1 the "sample mean and standard deviation" of the complete gaps, which reads as 1/(N-1), and says nothing about Appendix 1's divisor. The parsimonious reading is the opposite - 1/N is OUR deliberate choice, made to match eq (10)'s divisor, and `_multiprocess.py:80-81` says exactly that. | Withdrawn from the docstring and restated as observation plus our own choice. Verify against the journal text before the thesis cites it.

IMPORTANT | analyzers/checks/_multiprocess.py (divisor consequence) | Following from the above: our `gamma_hat` is smaller than the paper's estimator 1 by sqrt((N-1)/N), so our C1 statistic is LARGER by sqrt(N/(N-1)) - 1.4% at N = 36, ~12% on a five-event segment. That direction is anti-conservative, and short segments are routine on gapped records. It is a divisor choice and not an error, but it is a difference from the source. | No code change (the reason for 1/N is sound). State it in the instrument report, Increment C.

MINOR | tests/test_checks_published_values.py:96-101 (corrected) | The Laplace test's rationale claimed gamma_hat was "computed independently" - it cancels exactly, being the same call on the same floats. What the test genuinely pins is the `sqrt(12)*(sum T_i - N*tau/2)/(tau*sqrt(N))` scaling, which is real and non-circular. The reviewer also noted the eq (11) test is this file's ONLY pin on gap ORDERING, since eq (11) is order-sensitive and gamma is not - the docstring undersold it. | Rationale corrected.

MINOR | tests/test_checks_published_values.py:198 (corrected) | The skipped m>1 test took an unused `request` fixture and skipped unconditionally, so its "TO ACTIVATE" recipe could never fire - adding the data would not have un-skipped it. | Now reads `SMALL_BOWEL_SEGMENTS` from the fixture and self-activates; asserts against `SMALL_BOWEL_PUBLISHED` rather than re-typing the targets.

MINOR | (out of scope, pre-existing) | A bare `pytest` from the repo root dies with INTERNALERROR importing `monoliths/v13fig/test_carve.py`. The suite is green only when scoped to `pytest tests/`. | Add `monoliths` to a pytest norecursedirs, or keep invoking `pytest tests/`.

GREEN | tests/fixtures/load_haul_dump.py:21-28 | Fixture data verified independently by the reviewer: 36 entries, strictly increasing, min gap 3 h, T_36 = 1970, and seven published quantities reproduce simultaneously. A transcription error large enough to matter cannot survive that many agreements. | none
GREEN | tolerance audit | Every assertion's used margin is 1-2 orders inside its bound; each tolerance is ~2x the published print half-width. A +/-10 h error at any single failure time fails at least one assertion; a +/-1 h interior error is not always caught, which is below Table 2's printed resolution anyway. | state the bound in the instrument report

---

## Increment B - measured before changing anything

MEASURED | analyzers/shape_stats.py | The xi defect is REAL but DORMANT on this dataset. Measured across the whole T2* ladder: 279 windows at 3 us, 30 at 4 us, 1 at 5 us - and **zero windows with a tied X (margin) or a tied Y (forward time)**, on either axis. Reads sit at jittered timestamps so forward time takes distinct values, and margin is a continuous measurement. Where there are no ties, reference eq (8) reduces algebraically to `1 - 3*sum|dr|/(n^2-1)`, so the shipped estimator has been numerically CORRECT on every artifact produced so far, and the 2/5 null is the right null for tie-free data. | Fix it anyway: the tie-free form is right here only by accident, and the reference's caveat targets QUANTISED metrics (a fidelity ladder where most windows are one read long). But the fix must be verified as a no-op on real data, not assumed to be one, and the Part 5 tie experiment must use synthetic tied inputs because this data has none.

---

## Increment B - changes made

DONE | analyzers/checks/_multiprocess.py:54,116 | `gamma_hat` and `gamma_hat_batch` signature defaults changed GAMMA_TRUNCATED -> GAMMA_COMPLETE, matching their own docstrings. Verified no shipped number moves: every call site passes the estimator explicitly. | -
DONE | tests/test_shape_stats.py | NEW, 20 tests, written BEFORE the xi change because nothing else could have caught it. Covers xi (tie-free and tied), the exact finite-n ceiling `1 - 3/(n+1)`, tie_fraction, dcor and its four reason codes. Closes rank 1 of the Part 7 list. | -
DONE | analyzers/permutation.py | NEW shared paired-permutation step wrapping `scipy.stats.permutation_test(permutation_type="pairings")`. Seed required, never defaulted. Distinguished in its docstring from `analyzers/checks/_permutation.py`, which is segment-blocked and for duration vectors. | -
DONE | analyzers/shape_stats.py chatterjee_xi | Now reference eq (8), the tie-corrected form, with a constant-y guard. Verified three ways: matches an independent eq (8) reference on tied data, matches the old tie-free form exactly on tie-free data (4 seeds), and max |new - old| = 0.000e+00 across 279 real windows. | -
DONE | analyzers/shape_stats.py xi_p_value | Calibration now chosen by measuring ties: tie-free -> closed form, tied -> permutation. `xi_p_method` records which, as a new column. `xi_p_value_asymptotic` keeps the closed form for the Part 5 comparison and is labelled tie-free-only. | -
DONE | analyzers/check_ledger.py | Now catches `RuntimeError` and `subprocess.TimeoutExpired` from the C3 bridge, so a failed Rscript writes a `not computed` row instead of killing the job. | -
DONE | analyzers/shape_stats.py NON_NUMERIC_STAT_COLUMNS | Found while wiring `xi_p_method`: `for_windows` medianed EVERY column via `to_numpy(dtype=float)`, which raises on a categorical. The new set names the two non-statistic columns. | -
DONE | seed threading | `xi_seed` reaches `DistinguishBandInputs` -> `make_inputs_from_windows` -> `build_within_calibration_panel_data`. Not yet declared at job level, because on tie-free data the seed is never consulted; a job that meets tied data should declare it as a step kwarg.

OUTSTANDING | Cramer-von Mises (B3) | Not started. It is a whole instrument - derive the statistic from the bridge integral, cross-check against `scipy.stats.cramervonmises` at gamma=1, transcribe and pin the limiting CDF at 0.34730/0.46136/0.58061/0.74346, and write its tests. Deliberately not rushed at the end of a long session in a pass whose entire purpose is fidelity. | Next session, before the Increment B review gate.

MEASURED | analyzers/shape_stats.py xi_p_value | Design question: always-permute, or branch on ties? Measured over all 315 windows of the T2* ladder (median n = 8): closed form 0.17 s, always-permute 128 s - a 760x slowdown per panel render, and the panel renders per dataset and per composite. The branch stays. Note the cost is dominated by scipy's per-call overhead (the xi work itself is ~6 s), so an always-permute path is a tight loop away if ever needed. The ESTIMATOR is unbranched - eq (8) is the complete form and reduces to the tie-free one algebraically - so only the calibration branches, and only where the article says it must.

---

## Increment B - reviewer pass (independent, fresh context, 2026-08-13)

Scope: analyzers/shape_stats.py, analyzers/permutation.py, analyzers/checks/_multiprocess.py,
analyzers/check_ledger.py, analyzers/distinguish_band.py, panels/_within_calibration_compute.py,
tests/test_shape_stats.py. All gitignored; read directly, verified by measurement.

IMPORTANT | analyzers/permutation.py:96-100 | The `-np.inf` substitution for a non-finite resample statistic is ANTI-conservative and does the opposite of its own comment ("must not silently become a rejection"): a -inf null can never be >= the observed, so it is dropped from the numerator and the p-value SHRINKS. Measured on a statistic degenerate on ~half the pairings: p = 0.034 with -inf against 0.506 with the conservative +inf. Dormant for xi and dcor (pairings preserve both marginals, so degeneracy is all-or-nothing) but this module advertises itself to "any analyzer". | Raise, or substitute +inf so a degenerate resample counts as at least as extreme.
IMPORTANT | analyzers/permutation.py:79-82, :49-55 | "The null is always SAMPLED, never enumerated" is false at n <= 4, where (n!)^2 <= 999: measured n = 3 -> 36 enumerated draws, n = 4 -> 576. scipy drops the +1 adjustment on an exact test, so `resolution` = 1/(B+1) misstates the floor (reports 0.02703 at n = 3; the true floor is 1/36 = 0.02778) and `n_resamples` reports 36 rather than the 999 requested. 3-read windows are the modal case here. | Set `exact` from `n_resamples >= factorial(n)**2` and compute `resolution` accordingly, or state the enumeration in the docstring.
IMPORTANT | tests/test_shape_stats.py:99-146 | `_xi_eq8_reference` is never compared to `chatterjee_xi`. `test_the_two_forms_diverge_once_y_is_tied` compares the two test helpers to each other and never calls the implementation, so the ONLY behaviour this increment changed - tied y - is unpinned. (Verified good independently: max |impl - brute eq (8)| = 0 over 3000 random tied grids.) | One line: `assert chatterjee_xi(x, y) == approx(_xi_eq8_reference(x, y))` on the tied fixture.
IMPORTANT | tests/test_shape_stats.py | No test touches `xi_p_value`, `xi_p_method`, the new constant-y guard (shape_stats.py:114-118) or `analyzers/permutation.py`. The tie branch, the method string, and seed reproducibility ship unpinned in the same increment that added them - and per this ledger:29 an analyzer edit leaves the run identity byte-identical. | Three tests: tied -> `permutation`, tie-free -> `asymptotic`, same seed -> same p.
IMPORTANT | analyzers/t2star.py:175-186, analyzers/fidelity.py:248-258 | `xi_seed` is NOT threaded: both step functions call `build_within_calibration_panel_data` without it, so its default 0 is unreachable from any job and never reaches a Mermaid label. `fidelity` is exactly the quantised ladder that motivated the tie branch. This ledger:168 says only "not yet declared at job level", which understates it. Compounding: distinguish_band.py:158-162 says the field is "a constructor field rather than a default" while the code reads `xi_seed: int = 0`. | Add `xi_seed` to both step signatures and pass it; the docstring claim follows.
IMPORTANT | analyzers/shape_stats.py:369 | `pooled_xi` discards the calibration label (`p_pooled, _method = ...`), so the pooled p in the artifact's 4-tuple has no method beside it - contradicting `xi_p_value`'s own "the method is part of the answer, not metadata". Pooled samples are far more tie-prone than per-window ones (3480 reads at 3 us). | Widen the tuple to 5, or return a small dataclass.
IMPORTANT | analyzers/shape_stats.py:90-92 and tests/test_shape_stats.py:54-67 | The `1 - 3/(n+1)` ceiling is the TIE-FREE ceiling, now stated unconditionally on a tie-corrected estimator. Measured counterexamples: n = 12 with two distinct y gives xi = 0.8333 against a claimed ceiling of 0.7692; n = 20, four distinct y gives 0.8800 against 0.8571. The test named "hits the exact finite-n maximum" only pins the tie-free case. | Qualify the docstring "tie-free y" and rename or extend the test.
IMPORTANT | analyzers/shape_stats.py:202-210 | The branch measures ties and never n, so reference 10.3 caveat 2 (small-sample convergence contested; permutation-calibrate in the disputed range) is neither implemented nor recorded. The median window is ~7 reads. Measured type-I of the closed form on independent pairs: 0.0118 at n = 7, 0.0356 at n = 10, 0.0418 at n = 20, 0.0508 at n = 200, all at nominal 0.05 - CONSERVATIVE, so nothing is over-stated, but at n = 5 the smallest attainable p is 0.0385 and the test can barely reject at all. | Add an n floor to the branch, or log caveat 2 as a knowing deviation with these numbers.
IMPORTANT | this ledger:17 | "This p-value is what the panel annotation shows" is false. The shape annotation is rho2 / xi / limb rho / peak (_within_calibration_render.py:99-102) and the xi heatmap plots median and pooled xi (non_repairable.py:632 unpacks `_p` and discards it). No xi p-value reaches any figure today. The fix was still right; its stated urgency was not. | Correct the entry.

MINOR | this ledger:167 | "`for_windows` medianed EVERY column" is contradicted by the pre-change artifact: output/t2star_q1_070423_..._20260808_124131 has median keys WITHOUT `window_index`, so the old code already skipped it and the change generalised one skip into a set. Verified no key was lost: old and new median and defined key sets are identical. | Restate as generalising an existing skip.
MINOR | analyzers/shape_stats.py:149-154 | "this is where it stops being usable, not where it starts being wrong" describes a cutoff above zero; the value is 0.0, where the two coincide. | Reword or set the constant deliberately above 0.
MINOR | analyzers/shape_stats.py:202 | `max(tie_fraction(x), tie_fraction(y))` is stricter than the cited authority requires: Chatterjee's Thm 2.1 needs Y continuous, and a FIXED tie-break on X leaves the null distribution unchanged under independence. X-only ties buy a 760x-slower calibration for no validity gain. Harmless (the permutation is valid), but the docstring justifies it by the 2/5 null, which is a Y-ties argument. | Say why X ties count too, or branch on y alone.
MINOR | analyzers/shape_stats.py:206 vs analyzers/permutation.py:79 | "exact whatever the ties" against "the null is always SAMPLED". With 999 sampled draws the p is level-valid, not exact. | Align the two sentences.
MINOR | analyzers/check_ledger.py:416 | `RuntimeError` also catches `RecursionError` and `NotImplementedError`, so a genuine bug inside `c3.run` becomes a `declined:` row rather than a traceback - the same failure mode this repo's raise-don't-swallow rule exists to prevent. | Give the bridge its own exception class and catch that.
MINOR | analyzers/shape_stats.py:204 | Function-level import of `analyzers.permutation`. No cycle exists (permutation imports nothing from analyzers), so it can be top-level; a hidden import is harder to see. | Move it up.
MINOR | analyzers/shape_stats.py | 468 lines now carrying estimator, calibration, per-window aggregation and the shape curve. Approaching the unwieldy line; the calibration half (xi_p_value*, XI_* constants) is the natural split. | Watch it.

GREEN | analyzers/shape_stats.py:74-126 | Eq (8) matches reference 10.2 verbatim. `method="max"` is right for BOTH: `rankdata(y,"max")` = #{j : y_j <= y_i} and `rankdata(-y,"max")` = #{j : y_j >= y_i}; `ell` is correctly NOT reordered (the l-sum is over all i, permutation invariant). Agrees with a brute-force counting implementation to 0.0 over 3000 random tied grids, and is BITWISE equal to the old tie-free form on 2000 tie-free samples. | none
GREEN | regression claim, verified independently | output/t2star_q1_070423_57fb5e_20260808_124131 (pre-change) against ..._20260813_101723 (post): pooled xi AND pooled p bit-identical at all 10 thresholds, max |dxi| = 0 and max |dp| = 0 over all 234 windows. `xi_p_method` is in the artifact (231/2/1 windows, all `asymptotic`, max tie fraction 0.0). | none
GREEN | analyzers/checks/_multiprocess.py:54,116 | Default change cannot move a number: every call site passes the estimator explicitly (c1:83,110; c2:171,191; tests), and bench/ never calls `gamma_hat` at all. | none
GREEN | analyzers/permutation.py | `permutation_type="pairings"` is the right null (both marginals fixed, only the pairing destroyed); `alternative="greater"` is right for xi and dcor; seeding reproduces (same seed -> identical p, different seed -> different p). `ruff` clean, `pytest tests/test_shape_stats.py` = 20 passed. | none

---

## Increment B - reviewer findings, and what was done

CORRECTED | .claude/review-findings.md (my own earlier entry) | I wrote that xi's p-value "is what the panel annotation shows". It is NOT: `shape_annotation` renders rho2, xi, limb_rho and peak, and `non_repairable.py` unpacks the pooled p and discards it. No xi p-value reaches any figure today. The defect was still worth fixing - the p-value is in the materialized artifact and the ledger-style consumers read it - but the stated consequence was wrong. | claim removed.

FIXED | analyzers/permutation.py | `-inf` substitution for a non-finite resample statistic was ANTI-CONSERVATIVE and inverted its own comment: a -inf draw can never exceed the observed value, so it shrinks p. Measured by the reviewer at p = 0.034 with -inf vs 0.506 with +inf. Now raises: a statistic that cannot be evaluated on a resample of its own data is broken, and the caller must know.
FIXED | analyzers/permutation.py | "The null is always SAMPLED, never enumerated" was FALSE at n <= 4, where `(n!)^2 <= 999`: measured 36 draws at n=3 and 576 at n=4. scipy drops the +1 on an exact test, so `resolution` misstated the floor (0.02703 vs the true 1/36 = 0.02778). A three-read window is the MODAL case here, so this was the common path. New `exact` field; `resolution` now reports the floor that applied.
FIXED | tests/test_shape_stats.py | `_xi_eq8_reference` was never compared to `chatterjee_xi` - the tied path, the only behaviour this increment changed, was unpinned. Added `test_the_implementation_matches_eq8_on_TIED_data` over four seeds.
FIXED | tests/test_shape_stats.py | Nothing tested `xi_p_value`, `xi_p_method`, the constant-y guard, or `analyzers/permutation.py`. Added 13 tests: both calibration branches, the method label never lying, degenerate input, seed reproducibility, the -inf refusal, and the exact-enumeration regime at n <= 4.
FIXED | analyzers/t2star.py, analyzers/fidelity.py | `xi_seed` was not threaded through either adapter, so it was unreachable from any job and could never reach a Mermaid label - a violation of this repo's own seed rule, and `fidelity` is precisely the quantised ladder that motivates the tie branch. Both now accept and forward it.
FIXED | analyzers/shape_stats.py pooled_xi | Discarded the calibration label. Now returns it as a fifth tuple element; `DistinguishBand.pooled_xi_per_threshold` and the renderer's unpack updated to match.
FIXED | analyzers/shape_stats.py chatterjee_xi | The `1 - 3/(n+1)` ceiling is the TIE-FREE ceiling; the docstring asserted it unconditionally on a tie-corrected estimator. Measured counterexample: n=12, two distinct y -> 0.833 against a claimed 0.769. Both docstring and test now say tie-free.

OPEN | analyzers/shape_stats.py xi_p_value | The branch measures ties and never n, so reference 10.3 caveat 2 (n in 35-355 "sits squarely in the disputed range") is neither implemented nor recorded. Measured type-I of the closed form on independent pairs at nominal 0.05: 0.0118 (n=7), 0.0356 (n=10), 0.0418 (n=20), 0.0508 (n=200) - CONSERVATIVE, so nothing is overstated. Left open deliberately: adding an n floor would send every window to the 760x-slower path for a bias that runs the safe way. | Record in the instrument report.
OPEN | analyzers/check_ledger.py:416 | `RuntimeError` also catches `RecursionError` and `NotImplementedError`, so a genuine bug in `c3.run` becomes a `declined:` row. | A dedicated bridge exception would keep the widening honest. Increment C.
OPEN | analyzers/shape_stats.py | 468 lines carrying estimator + calibration + aggregation + curve. The calibration half is the natural split. | Not this pass.

---

## Increment B3 - Cramer-von Mises (new instrument)

NEW | analyzers/checks/cvm_cramer_von_mises.py + tests/test_checks_cvm.py | The fourth
functional of the same Brownian bridge as C1 and C2. Reference 5.2 calls its omission "harder
to defend" and recommends implementing it, because for m > 1 the source paper REVERSES its
single-process preference and drops Anderson-Darling "as the Cramer-von Mises test had better
level properties in this case". Gapped records are routine here and m > 1 IS the gap case.
16 tests, 247 passing overall.

CONFIRMED | reference PDF eq (7), CvM middle term | The document prints
`- i N (T^2_{i+1} - T^2_i)/tau`, dimensionally inconsistent with its own first term
`i^2 X_{i+1}/tau`. Deriving `integral_0^1 (N(s tau) - sN)^2 ds` from the definition gives
`/tau^2` and reproduces the document's own tail term `N^2[u_N^2 - u_N + 1/3]` exactly.
Settled by THREE independent routes, not by preference:
  - piecewise `scipy.integrate.quad` on the step function  -> 0.117508047340
  - the closed form implemented here                       -> 0.117508052657
  - `scipy.stats.cramervonmises(T/tau,"uniform")` at gamma=1 -> 0.117508052657
`test_the_bracket_is_the_bridge_integral` pins route 1 against route 2 permanently, so this
never has to be re-litigated. Transcribing the printed equation literally would have sent the
pass hunting a bug in the authoritative document. Log against the reference, not the code.

GREEN | limiting CDF | Anderson-Darling (1952) series with `kv(0.25, .)`. Reproduces all four
published critical values: 0.34730 -> 0.100003, 0.46136 -> 0.050000, 0.58061 -> 0.025001,
0.74346 -> 0.010000. Max error 3.08e-06. No finite-N correction, for the same reason `adinf`
carries none: the statistic is divided by an ESTIMATED gamma_hat^2, so the finite-N null is
not the finite-N CvM null and a fitted correction would be a precision claim the derivation
cannot carry.

GREEN | external cross-check | At gamma = 1 the statistic equals `scipy.stats.cramervonmises`
to 4.4e-16 / 5.3e-16 / 6.8e-15 at n = 13 / 36 / 120. That is an independent implementation
reached by a different route, which is stronger evidence than a second transcription of the
same formula - the same standard `_classical_ad` already holds C2 to.

GREEN | the m > 1 motivation, made concrete | `require_strict_tau=False` is not a relaxed
guard: eq (7) carries a `1/(s(1-s))` weight so `ln(tau/(tau - T_N))` is +inf when the last
event lands on the truncation time, and C2 REFUSES that record. CvM's integrand has no such
weight and returns 0.165156 on the same input. On the in-spec clock of a carved record that
case is the common one, not a corner. `test_it_survives_a_record_that_c2_declines` pins it.

DELIBERATE | battery.ROW_KEYS | CvM is NOT registered, per the plan. ROW_KEYS is the schema of
the bench tables; an unregistered check makes `bench_acceptance_at_n` return None, so every CvM
ledger row would read `underpowered / no bench cell`. Promotion is one line here plus a ~106
minute bench re-run, and is a separate decision.
KNOWN GAP | asymptotic calibration | Refused for m > 1 and it raises rather than degrading: the
paper's normal approximation is to the WEIGHTED sum, and this implementation sums unweighted,
exactly as our C2 multi-process route does. Permutation is exact under iid gaps and needs no
limit. The unweighted m > 1 sum remains the one construction with no external validation
anywhere - same gap the skipped small-bowel test records.

---

## Increment B - review gate: MUTATION TESTING (the /review subagent died mid-pass)

The code-reviewer subagent hit the session limit at the exact moment it began mutation testing,
after ~9 minutes and 39 tool calls. It wrote NOTHING to this ledger, so nothing from that run is
salvageable and its reading passes are not recorded. Mutation testing was completed directly
instead, since that was the one step outstanding and it needs measurement, not fresh context.

Method: 14 deliberate defects introduced one at a time into `shape_stats.py`,
`cvm_cramer_von_mises.py` and `permutation.py`; the guarding test module run against each; the
file restored from the original source string in a `finally` block. Working tree verified intact
afterwards (`git status` shows only the pre-existing ruff-format re-wraps; 249 passed).
This is the standard CLAUDE.md demands - "a positive control must fail when the thing it guards
is broken" - applied to this increment's tests rather than asserted of them.

CAUGHT (12 of 14, each naming the test that caught it):
  xi  ell min-ranks instead of max        -> test_the_implementation_matches_eq8_on_TIED_data
  xi  r average-ranks (the tie-free habit)-> test_the_implementation_matches_eq8_on_TIED_data
  xi  factor 2 dropped from denominator   -> test_a_perfect_function_hits_the_exact_finite_n_maximum
  xi  XI_TIE_CUTOFF 0.0 -> 1.0            -> test_a_tied_window_switches_to_permutation
  cvm /tau^2 -> /tau (THE PDF's TYPO)      -> test_at_gamma_one_it_is_the_classical_cramer_von_mises
  cvm 1/3 tail term dropped               -> test_at_gamma_one_it_is_the_classical_cramer_von_mises
  cvm gamma^2 -> gamma                    -> test_statistic_batch_reproduces_the_statistic_under_the_identity
  cvm batch blocks misaligned             -> test_statistic_batch_reproduces_the_statistic_under_the_identity
  cvm CDF series 4j+1 -> 4j+3             -> test_the_limiting_cdf_reproduces_the_published_critical_values
  cvm m>1 asymptotic refusal removed      -> test_the_asymptotic_calibration_is_refused_for_multiple_segments
  perm -inf substitution reinstated       -> test_paired_permutation_refuses_a_statistic_that_goes_non_finite
  perm exact flag forced False            -> test_small_windows_get_an_enumerated_exact_null
The CvM typo mutant is the one that matters most: reintroducing the reference document's printed
`/tau` is caught immediately by the scipy cross-check, so that transcription can never come back.

FALSE CATCH, corrected | the 13th mutant | "ell reordered by argsort(x)" was reported CAUGHT, but
the mutation script injected an undefined name and the tests died on `NameError`, not on the
claim. Redone properly (`ell = rankdata(-y,"max")[order]`): the reorder is EXACTLY a no-op, max
|diff| = 0.000e+00 over 400 heavily-tied grids. The earlier GREEN entry - that the l-sum is over
all i and therefore permutation invariant - is independently confirmed, and a test asserting
otherwise would be wrong. Recorded because a mutation that dies on an error is not evidence.

FIXED | tests/test_shape_stats.py | THE SURVIVOR, and a real gap. Replacing the asymptotic null
variance `2.0/5.0` with `1.0/5.0` left all 35 tests passing. That constant scales every p-value
the tool actually emits: all 279 windows of the real T2* ladder are tie-free, so `xi_p_method`
is `asymptotic` on every one and this is where their significance is decided. The two tests that
existed only asserted `xi_p_value` agrees with `xi_p_value_asymptotic` - the same arithmetic
twice, which is exactly why the mutant slipped through.
Two tests added, both SIMULATED so they are an independent route to the constant rather than a
restatement of it: `test_the_asymptotic_null_variance_is_two_fifths` (Chatterjee Thm 2.2,
measured var(sqrt(n) xi) = 0.4117 at n=100 and 0.3982 at n=400 against the nominal 0.4) and
`test_the_closed_form_holds_its_nominal_level_under_independence` (type-I 0.049 at nominal 0.05,
0.055 at n=400; at n=60 the level holds one-sided). The 1/5 mutant is now CAUGHT.

STATED LIMIT | tests/test_shape_stats.py | A subtler mutant, `2/5 -> 0.45`, still SURVIVES. That
is the Monte Carlo resolution of 1500 replicates, not an oversight: the SE of the variance
estimate is about 0.015, so a 12.5% shift sits inside the band. Tightening the band to catch it
would make the test flaky, and raising the replicate count buys resolution at suite runtime. The
test catches a WRONG constant, not a slightly-off one, and this line says so rather than letting
the next reader assume more.

STATE AT THE GATE | 249 passed, 2 skipped. ruff at the one known pre-existing F841 (`raw_meta`).
Style ratchet test green at 17. Working tree carries no mutation residue.

STILL UNREVIEWED BY A COLD READER | Increment B never received a completed adversarial read. The
mutation evidence above is strong on the STATISTICS, and weaker on everything a reviewer reads
for and a mutant cannot express: the seed threading through `t2star`/`fidelity`/
`_within_calibration_compute`, the 5-tuple unpack sites, `check_ledger`'s widened except, and whether
the docstrings claim more than the code does. Respawn `/review` on exactly that scope when the
limit resets, before Increment C is committed on top of it.

---

## Increment B - COLD-READ reviewer pass (independent, fresh context, 2026-08-13)

Scope: the four things mutation testing is blind to - seed threading, the 5-tuple arity,
`check_ledger`'s widened except, `_multiprocess` defaults - plus docstring honesty.
Verified by measurement where measurable. Gitignored files read as explicit paths.

IMPORTANT | analyzers/t2star.py:139, analyzers/fidelity.py:227, panels/_within_calibration_compute.py:402, analyzers/distinguish_band.py:163,177, analyzers/shape_stats.py:276,345,416 | `xi_seed` is threaded through the adapters but is UNREACHABLE FROM ANY JOB and therefore never on a Mermaid label. `jobs/active/t2star_q1_070423.py:74-88` and `:100`, its 100423 twin, and `jobs/common.py:174-184` are the only step wrappers; none of the three takes or forwards `xi_seed`, so the value is `0` by default at all seven levels. runner.py:46-48 builds the label from `node.kwargs` only, so a defaulted argument cannot appear. This is exactly the rule the same job file states two lines above at t2star_q1_070423.py:107 ("Declared here, not defaulted in the builder, so both reach the provenance label") for `shape_min_reads`. Ledger:215 ("Both now accept and forward it") is true of the analyzer layer and NOT of the job layer, where the seed rule actually binds. | Add `xi_seed: int` to `_t2star_panel_data` / `_fidelity_panel_data` and declare `xi_seed=0` in the `job.step(...)` call, as `shape_min_reads` already is.

MINOR | analyzers/shape_stats.py:284,380,456 | One seed serves every window and the pooled sample. Measured: two DIFFERENT datasets of the same n at the same seed draw a BIT-IDENTICAL permutation sequence (50/50 draws matched). So per-window p-values are not independent draws - windows of equal n share their Monte-Carlo error, which then does not average down in `medians["xi_p_value"]` (for_windows:463-475 medians that column) nor in any count of significant windows. Dormant on T2* (all windows tie-free -> asymptotic branch) but live on the fidelity ladder, the quantised case the tie branch exists for. At B=999 the shared noise is ~0.01 in p, so this changes no verdict today. | Derive the per-window seed, e.g. `seed=xi_seed + window_index`, or state in `for_windows` that the MC error is common across windows by design.

MINOR | analyzers/shape_stats.py:360 | `pooled_xi`'s docstring still says "Returns (xi, p_value, n_reads_pooled, n_windows_pooled)" - four elements - while the signature and the body return five. The arity change of this increment is unstated in the one place a consumer reads. | Add `method` to the Returns line.

IMPORTANT | panels/within_calibration.py:676-678 | The thin-pooling caveat is keyed on a LEAKED loop variable. `label` is exhausted by the `for label in labels` loop at :626, so the `for col, ... in enumerate(notes)` loop at :676 tests `pooled_is_thin_per_threshold[LAST label]` on every column and appends the LAST label. Measured on the shipped artifact output/t2star_q1_070423_57fb5e_20260813_105641: thresholds 4-9 us are all thin, `10 us` is not - so `thin_pooling` comes out EMPTY and the figure prints NO thin-pooling warning at all. The caveat that six of ten pooled cells are not within-excursion statistics is silently suppressed on the real data; had the last threshold been thin instead, the note would have named it ten times. `thin_pooling` does not exist in HEAD, so this is new in this working tree. | Loop `for col, label in enumerate(labels)` (or zip labels with notes) and test that label.

MINOR | panels/within_calibration.py:178-180 | `pooled_lookup` defaults to a ONE-tuple `(float("nan"),)` and slices `[0]`. It therefore keeps working under any arity change to `pooled_xi_per_threshold` - the exact silence the 5-tuple change should have been able to trip. The sibling unpack at :636-640 does name five and would have failed loudly. | Default to the same 5-tuple the unpack site uses, so both sites break together.

MINOR | panels/_artifact_guard.py:36 with panels/within_calibration.py:636 | The arity change is invisible to the staleness guard, as the guard's own contract admits (field NAMES only). Measured: output/t2star_q1_070423_57fb5e_20260808_124131 (pre-P5) still carries 4-tuples in `pooled_xi_per_threshold`, `__setstate__` accepts it without complaint, and the 5-unpack then fails at DRAW time with "not enough values to unpack (expected 5, got 4)" - reached by `compare_t2star_0704_vs_1004 --reuse-deps`. It raises rather than lying, so no wrong number ships, but the guard's stated purpose ("fail loudly at the pickle boundary") is defeated for any value-shape change. | Out of this increment's scope to fix; it is ledger:29 item 2 again, now with a concrete instance to cite.

IMPORTANT | analyzers/checks/c3_serial_copula.py:100 (`--vanilla`) with check_ledger.py:416 | C3 CANNOT SUCCEED IN THIS ENVIRONMENT, and the P5 widening is what hides it. `copula` is installed in the USER library `/home/sera/R/library`; `--vanilla` implies `--no-environ --no-init-file`, so `.libPaths()` drops to the two system dirs and `requireNamespace("copula")` is FALSE. Measured end to end: `c3.run([Segment(x=rng.exponential(size=30), tau=...)], seed=0)` raises `RuntimeError: Rscript exited 1 ... Error: the 'copula' package is required for C3`, while a plain `Rscript -e 'requireNamespace("copula")'` returns TRUE. Before P5 this killed the ledger job, which is loud; after P5 every C3 row becomes `declined: Rscript exited 1 ...` / verdict `not computed` and the check silently never runs. The widening is right; landing it without first running the bridge once means the first real R invocation in this repo's history fails and nobody is told. | Drop `--vanilla` for `--no-save --no-restore` (keeps R_LIBS_USER), or pass `env={"R_LIBS_USER": ...}`; then run C3 once for real before shipping.

IMPORTANT | analyzers/checks/c3_serial_copula.py:13-18, :85-88, c3_serial_copula.R:9-12 | Three docstrings assert a fact about the environment that is now false: "R is absent on this machine, so C3 is UNASSESSED", "none of the code paths past `_invoke_rscript` have ever executed", "UNEXERCISED: R is absent on the machine this was written on". Measured: `Rscript` 4.5.3 is on PATH at /usr/bin/Rscript and `_invoke_rscript` now executes on every ledger run. check_ledger.py:407-415 already says "R is now installed" - the two modules contradict each other in the same tree. Under this repo's claims rule a docstring must match the state it ships in. | Restate as "R may be absent; when it is, ..." and record the one real invocation once the `--vanilla` problem above is fixed.

GREEN | analyzers/check_ledger.py:415-429 with _verdict:215-216 | A `declined:` row CANNOT be mistaken for a pass. `p_value=None` is the first precedence branch, so the row is `not computed` with the exception text in `notes`, and it never reaches the `p <= alpha` or the underpowered branches. Verified on the shipped ledger output/check_ledger_q1_7174d1_20260812_104552: all 12 C3 rows read `not computed` / "R unavailable". One caveat worth knowing: `pd.DataFrame` coerces `p_value=None` to NaN in the float column, so downstream the None-ness survives only in `verdict`/`notes`, not in `p_value`.

IMPORTANT | analyzers/checks/c3_serial_copula.py:143 (`timeout_s: float = 120.0`) | The 120 s budget is not adequate at the sizes this ledger reaches. Measured on this machine with the user library forced onto the path: `serialIndepTestSim(n=355, lag.max=5, N=1000)` + `serialIndepTest` takes 97.1 s - 81% of the budget, on an idle machine. The cost grows superlinearly in n, and the ledger calls C3 once per threshold per clock (up to 20 invocations, so a full run is ~30 min of R even when nothing times out). Any load spike, or a threshold with more events than 355, converts a valid check into `TimeoutExpired` -> `declined:` -> `not computed`. The widened except makes that outcome quiet. | Raise the default well above the measured 97 s (600 s is the same order as the other bench costs), and put the measured 97 s at n=355 in the docstring so the next reader can scale it.

MINOR (sharpens the OPEN entry at this ledger:220) | analyzers/check_ledger.py:416 | Confirmed: `RuntimeError` is the base of both `RecursionError` and `NotImplementedError` (checked), so a genuine bug under `c3.run` degrades to a `declined:` row. Two things sharpen it. (a) The over-catch is not new: `ValueError` was already in the tuple and `c3.run:149` calls `validate_segment`, whose ValueError means the CARVE produced an invalid segment - a data-contract violation, silently downgraded to `declined:` since before P5. (b) The only RuntimeErrors that are genuinely the bridge's are the five in `_invoke_rscript` (c3:92,118,123,130,135); nothing else in `analyzers/checks/` raises one. So a dedicated `RBridgeError` around those five, caught alone, costs one class and removes the whole over-catch. | Define `RBridgeError(RuntimeError)` in c3, raise it at those five sites, and narrow the except to `(RBridgeError, subprocess.TimeoutExpired)` - keeping `ValueError`/`KeyError` only if the carve contract is deliberately tolerated here.

GREEN | analyzers/checks/_multiprocess.py:54,116 (default change) | Claim re-verified independently and repo-wide. `gamma_hat`/`gamma_hat_batch` are called at exactly six sites - c1:83,110; c2:171,191; cvm:146,162 - and every one passes the estimator positionally from a `gamma_estimator` parameter that itself defaults to `GAMMA_COMPLETE` (c1:73,94,121; c2:157,178,201; cvm:132,153,172; battery:86). Tests pass it explicitly at every call (test_checks_published_values.py:104,118,134,171,184; test_checks_statistics.py:170,207,210,220,222). `bench/` never names `gamma_hat` at all, and neither does `monoliths/`. The default was unreachable before and after. | none

IMPORTANT | analyzers/shape_stats.py:86-88 and :189-190 | The docstring's measured count is not in any artifact. "0 of 279 windows had a tie on either axis" - the shipped artifact output/t2star_q1_070423_57fb5e_20260813_105641 has 234 windows in the shape tables (2/231/1 at 2/3/4 us) and the 100423 artifact has 75, so neither 279 nor any pairing of them reproduces it. The CLAIM is correct - recomputed from the read table, max tie fraction is exactly 0.0 on BOTH axes over all 234 eligible windows and over all four non-degenerate pooled samples - only the number is unsourced. This is the exact defect CLAUDE.md's second claims rule names. | Recompute the count into the docstring (234 on 0704, 309 across both datasets) or cite the artifact and drop the figure.

MINOR | analyzers/shape_stats.py:11 and :30 | "the median window here has ~7 reads" names a single cell as if it were the ladder. Measured: 7.0 is the median over COMPLETE windows at the 3 us threshold ALONE. Over the 234 windows that actually reach these statistics the median is 11.0; over all 1212 complete windows on the ladder it is 1.0, and 978 of them have four reads or fewer. The number chosen matters, because :30 uses it to argue the asymptotics are marginal. | Say which threshold and which population, e.g. "median 11 reads over the 234 eligible windows; 1 read over all complete windows".

IMPORTANT | analyzers/permutation.py:49-50 and :93-94 | "a three-read window is the modal case in this project, so this is the common path, not a corner" and "which is the common case here" are both false for the only consumer. `xi_p_value` is called from `window_stats` (via `for_windows`, which receives only windows already filtered to `n_reads >= shape_min_reads`, shipped value 5 in both jobs) and from `pooled_xi` (thousands of reads). Measured on the shipped artifact: 0 of 234 windows have n <= 4, minimum 5. So the exact-enumeration regime the `exact` field and `resolution` exist for is UNREACHABLE from the shipped pipeline, not the common path. The modal complete window does have 1 read - but those never reach xi at all. | Reword to "unreachable at the shipped `shape_min_reads=5`; the field exists for callers that lower it", and keep the tests, which are the only exercise of that path.

GREEN | analyzers/shape_stats.py:74-103, 151-156, 159-168 | Every quoted claim checked against the reference PDF: eq (8) is transcribed exactly as printed in 10.2 (numerator over i=1..n-1, denominator 2*sum over i=1..n, r = max-rank, l = reverse max-rank), the tie-free reduction `1 - 3 sum|dr|/(n^2-1)` matches, the XI_TIE_CUTOFF quote is verbatim from caveat 1, and "roughly valid even for n as small as 20" is verbatim from caveat 2. The declared deviation (deterministic stable tie-break where Chatterjee breaks ties "uniformly at random") is stated in the PDF in those words and is honestly flagged as a deviation. The tied-ceiling numbers verified by exhaustive search: n=12 two distinct y -> 0.8333 against a tie-free ceiling of 0.7692. | none

IMPORTANT (CRITICAL the day CvM is promoted) | analyzers/checks/cvm_cramer_von_mises.py:73-76, :186 | The limiting CDF is truncated at 12 terms and the series does NOT converge in the upper tail, so the asymptotic p-value is NON-MONOTONE in the statistic and eventually rises back toward 1. Measured: p is minimised at z = 7.66 (p = 1.0e-15) and thereafter INCREASES - p(20) = 1.2e-8, p(100) = 3.8e-3, p(500) = 0.100, p(1e5) = 0.705. A larger statistic gets a larger p-value. Constructible verdict flip: N = 682 near-constant durations with a step-change trend (gamma_hat ~ 0.005, exactly the loose-threshold regime this project has - 682 is the complete-window count at the 4 us threshold, and the reference itself says carved durations "are integer multiples of the read spacing") gives statistic 794.4 and p = 0.156, i.e. NOT REJECTED for an overwhelming trend. Nothing ships it today because CvM is deliberately unregistered in `battery.ROW_KEYS`, so this is latent, not live - but promotion is described as "one line plus a bench re-run". The four pinned critical values are all z <= 0.75 and cannot see this. C2's `adinf` was checked for the same class and is clean (monotone, saturates at p = 0). | Guard the tail: return `p = 0.0` above the z where the series stops converging (z ~ 8 is already p < 1e-15), or raise _CDF_TERMS and assert monotonicity. Add a test that `1 - cvm_limiting_cdf(z)` is non-increasing over z in [0.05, 1000].

MINOR | analyzers/checks/cvm_cramer_von_mises.py:72-74 | "It converges fast - the j = 3 term is already below 1e-12 across the range that matters - so this is generous rather than tuned." Measured: true only for z <~ 0.8. The j = 3 term is 3.5e-6 at z = 2, 2.0e-3 at z = 5 and 4.3e-2 at z = 20. "The range that matters" is doing unstated work, and it is what makes the finding above invisible. | State the range (z <= 1, the published-critical-value band) and say what happens outside it.

MINOR | analyzers/checks/cvm_cramer_von_mises.py:1-40 | No source equation number, against CLAUDE.md's "every transcribed formula cites its source to equation number". The docstring cites "the reference document Section 5.2" and "its eq (7)", but 5.2's own title is "Cramer-von Mises, eq. (6)" - eq (6) of Kvaloy & Lindqvist (arXiv:1802.08339), the same paper c2 cites by arXiv id and equation. Worse, `c2_anderson_darling._eq7` already means K-L eq (7), which is the ANDERSON-DARLING statistic, so "eq (7)" now names two different functionals in one package. | Cite "Kvaloy & Lindqvist, arXiv:1802.08339 eq (6) (printed as eq (7) of the reference document)".

MINOR | analyzers/checks/cvm_cramer_von_mises.py:8-14 | The m = 1 side of the paper's verdict is omitted. The reference says plainly "For a single process the paper drops it, since this test had less power than the Anderson-Darling test", and `run`'s DEFAULT calibration is the m = 1 asymptotic one. As written the docstring reads as an unqualified endorsement. | Add the m = 1 sentence beside the m > 1 one.

SHARPENED | analyzers/checks/c3_serial_copula.py:143 | The n = 682 timing measurement did not finish inside a 580 s ceiling, so C3 at the 4 us threshold (682 complete windows on the shipped 0704 record) costs at least 8 minutes against a 120 s budget. Combined with the measured 97 s at n = 355, the default guarantees a `TimeoutExpired` -> `declined:` row on the loosest thresholds - which is where the tie problem is worst and the check is most wanted.

MINOR | analyzers/distinguish_band.py:158-162 | The comment still asserts two things the code does not do: "a constructor field rather than a default" (the line below it is `xi_seed: int = 0`, a default) and "it reaches the provenance label like every other step kwarg" (it reaches no label - see the first finding of this pass). The previous pass flagged the first half; the second half was not corrected. | Say what is true: defaulted at every level, and not yet declared by any job.

MINOR | docs/iid_checks/C3_serial_copula.md:3,31,47, LIMITATIONS.md:42, README.md:90 | Same false environment claim as the c3 module, in four more places: "R is absent on this machine", "every row", and a status row reading "C3 | UNASSESSED". Rscript 4.5.3 is on PATH. README.md:90 is also a status table in an agent-facing doc, which CLAUDE.md forbids for exactly this reason - it went stale in one day. | Restate conditionally and drop the status row.

MINOR | docs/iid_checks/ | CvM shipped as a new instrument with no doc, in the one directory that carries a page for every other check (C1, C2, C3, C5, C6 + LIMITATIONS + README). Not a request to write one - CLAUDE.md says do not create long .md docs unprompted - but the asymmetry should be a deliberate choice, and README.md's own table does not mention CvM exists. | Note CvM in the README's list, or record deliberately that it is undocumented until promoted.

MINOR | panels/within_calibration.py | 1039 lines. The renderer now carries eight axes plus the heatmap, timeline, survival and summary builders. `_within_calibration_render.py` exists precisely as the functions-of-axes half; the per-axis builders are the natural split. | Watch it; the repo rule says flag files that grow unwieldy.

GREEN | baseline reproduced | `PYTHONPATH=. python -m pytest tests/ -q` = 249 passed, 2 skipped. `ruff check .` = the one pre-existing F841 (`raw_meta`, plots/interpolation_stage_plot.py:63) and nothing else. No source file was modified by this review; the only write is this ledger.

MINOR | jobs/active/t2star_q1_070423.py:74-88 (and its 100423 twin), jobs/common.py:174-184 | `xi_seed` is not the only parameter stranded on a default. `k` is declared on the WINDOWS step (`k=1.0`) but the panel adapters never receive it, so `make_panel_data`'s own `k: float = 1.0` decides the unresolved-band half-width that the figure draws. Today both are 1.0 and nothing is wrong; the moment a job carves at `k=2.0` the panel draws a k=1.0 band over a k=2.0 carve and says nothing. `use_uncertainty` IS threaded from the same result, which is what makes the omission visible. `windows.py:525` already puts `k` in `WindowsResult.meta`, so the fix is one line in each adapter, exactly as `use_uncertainty` is read today. | Read `k` from `window_result.meta` alongside `use_uncertainty`, or declare it as a step kwarg on the panel step too.

GREEN | analyzers/shape_stats.py:200-212, :362-373 | The `p_value=None` rule is met without a None. A p that was not computed always carries `xi_p_method == "none"` beside a NaN, and both branches that CAN produce a number set `asymptotic` or `permutation`; there is no path where a missing p is indistinguishable from a large one. Verified on the shipped artifact: 234/234 rows read `asymptotic` with a finite p, and the three degenerate pooled cells (1 us, 8-10 us) read `none` with NaN.

MINOR | analyzers/distinguish_band.py:96-104 | `DistinguishBand` records `shape_min_reads`, `k` and `use_uncertainty` as result fields but NOT `xi_seed`. So even once a job declares the seed, the artifact that carries the permutation p-values will not say which seed produced them; a reader would have to reconstruct it from the Mermaid label of the run that wrote the pickle. The other three settings that shape these numbers are all carried. | Add `xi_seed: int = 0` to the band alongside them (note this changes the field-name set, which IS what `StaleArtifactGuard` checks, so old artifacts will be caught properly for once).

GREEN | panels/_within_calibration_compute.py:400-466 | The seed threading through the compute layer itself is complete and drops nothing: `build_within_calibration_panel_data(xi_seed=...)` -> `distinguish_band.make_inputs_from_windows(xi_seed=...)` -> `DistinguishBandInputs.xi_seed` -> both `pooled_xi` and `for_windows`. The break is entirely at the job boundary. | none

---

## Increment B - fixes applied after the COLD-READ pass

FIXED | panels/within_calibration.py:676 | The leaked `label`. `for col, (...) in enumerate(notes)`
read the loop variable left over from the build loop at :626, so every column tested the LAST
threshold. Measured consequence on the shipped artifact: 4-9 us are thin and 10 us is not, so
`thin_pooling` came out EMPTY and the figure printed no thin-pooling caveat at all - six of ten
pooled cells are not within-excursion statistics and the figure said nothing. Now
`enumerate(zip(labels, notes))`, with a comment naming the failure so it is not re-introduced.

FIXED | analyzers/checks/cvm_cramer_von_mises.py | The limiting CDF diverged outside the range
anything tested. `w = (4j+1)^2/(16z)` goes to zero as z grows, `K_{1/4}` diverges, and the
truncated 12-term sum turns around: p bottomed at z = 7.66 then CLIMBED, giving p(100) = 3.8e-3,
p(500) = 0.100, p(1e5) = 0.705. The reviewer's constructed case - 682 near-constant durations
with a step trend, statistic 794.4 - returned p = 0.156 and FAILED TO REJECT. All four pinned
critical values sit in [0.34, 0.75], which is exactly why 16 tests missed it.
Delegating to scipy does NOT fix it: `scipy.stats._hypotests._cdf_cvm_inf` breaks the same way,
its p rising after z = 5 (4.87e-11 at z=5, 1.81e-9 at z=10, 3.29e-6 at z=1e5).
Fix: evaluate the series only for z <= `_SERIES_Z_MAX` = 5.0 and saturate at F = 1 above it.
The cutoff is measured, not guessed - three independent lines: the published criticals all lie
below 0.75; an INDEPENDENT Karhunen-Loeve Monte Carlo (`W^2 = sum_k Z_k^2/(k pi)^2`, 4e6 draws)
agrees to 4e-5 absolute through z = 2 and produced no draw above 2.904; and the dominant
eigenvalue bounds the tail by `exp(-pi^2 z/2)` = 2e-11 at z = 5. Four tests added, including
monotonicity over `z` in [0.02, 1e6] and the reviewer's verdict-flip case, which now gives p = 0.

FIXED | analyzers/checks/c3_serial_copula.py | `--vanilla` implies `--no-environ`, dropping
`R_LIBS_USER` and with it `/home/sera/R/library` from `.libPaths()` - measured, plain `Rscript`
sees it and `Rscript --vanilla` does not - so `copula` was invisible and EVERY C3 call raised.
Post-B4 that no longer killed the job; it silently became a `declined:` row, which is the worse
failure. `--vanilla` is kept (a reproducibility tool should refuse the user profile) and the
search path is now passed explicitly as `R_LIBS`, asked of R itself by a new `r_library_paths()`
rather than hard-coded to this machine.

C3 FIRST CONTACT | 2026-08-13, Rscript 4.5.3 | The check has now RUN, for the first time since it
was written. On iid exponential input: n=50 -> statistic 0.00579, p 0.958, 3.9 s; n=150 ->
0.00713, p 0.904, 14.7 s; n=355 -> 0.00763, p 0.866, 130.2 s. Failing to reject data that
satisfies the null is the expected outcome and is a SMOKE TEST, not calibration - C3 still has no
bench cell and therefore no size or power evidence. Five docs/modules asserting "R is absent on
this machine" corrected; `docs/iid_checks/README.md` now reads RUNS, UNCALIBRATED.

FIXED | analyzers/checks/c3_serial_copula.py:182 | Timeout 120 -> 900 s. The reviewer's 97.1 s at
n=355 is confirmed and exceeded here (130.2 s), so the old default was ALREADY overrun by a size
the ladder produces. Cost grows about as n^2.8; the docstring records the three measured points
and states plainly that n = 682 is still not covered (extrapolates to ~800 s, measured unfinished
at 580 s) and becomes `not computed` with the timeout recorded.

FIXED | jobs/active/t2star_q1_070423.py, t2star_q1_100423.py, jobs/common.py | `xi_seed` reached
the analyzers but no job declared it, so it was 0 at every level and - because `runner.py` builds
the Mermaid label from `node.kwargs` alone - a defaulted argument could never appear in
provenance. The ledger's earlier FIXED was true of the analyzer layer only; the rule binds at the
job. Now declared on both t2star panel steps and on `configure_ramsey_job`. VERIFIED in the
shipped artifact: `step_t2star_panel_data["t2star_panel_data(shape_min_reads=5,
use_uncertainty=True, xi_seed=20260813)"]` in q1_27h_0704_dataset_t2star.prov.md.

FIXED | analyzers/shape_stats.py:86,189 | "0 of 279 windows" was in no artifact. Recomputed from
the two shipped panel pickles: 234 window-rows on 0704 (d611af) and 75 on 1004 (731c02), 309 in
total, `xi_p_method` = `asymptotic` on every one and no tie on either axis. The CLAIM was right
and the NUMBER was invented; both datasets are now named separately, since pooling two into one
figure is what produced the wrong total.

FIXED | analyzers/permutation.py:49,93 | "a three-read window is the modal case ... the common
path" is false for the only consumer: `for_windows` admits only `n_reads >= shape_min_reads` (5),
and 0 of the 234 windows on 0704 have n <= 4. The exact-enumeration regime is UNREACHABLE from
the pipeline. Reworded to say so, and kept documented because a direct caller can still land
there and must not read a floor of 1/36 as though it were 1/1000.

STATE | 253 passed, 2 skipped. ruff clean across `analyzers/` and `jobs/` (the one known
pre-existing F841 is in a panel module, untouched). Both t2star jobs re-run; identities moved
(57fb5e -> d611af, c3ab94 -> 731c02) because the job files changed, which is the identity
signal the earlier analyzer-only edits could not produce.

NOT FIXED, carried to Increment C | the reviewer's four MINOR seed findings: one seed serves every
window (windows of equal n draw a bit-identical permutation sequence, so their Monte Carlo error
does not average down in `medians["xi_p_value"]`); `distinguish_band.py:158` still claims the
field reaches the label from the analyzer layer; `DistinguishBand` does not carry `xi_seed` as a
field; and `k` is declared on the windows step but never reaches the panel adapter, which uses
its own `k=1.0`. The last is a pre-existing defect this pass surfaced rather than caused.

---

## Increment B - AUDIT of the seven post-cold-read fixes (fresh context, 2026-08-13)

SCOPE: the seven FIXED entries above only. COMMIT: 2f1bd95.
Manifest (checked off as each is audited, findings appended at the moment of discovery):
- [ ] 1 panels/within_calibration.py leaked `label`
- [ ] 2 analyzers/checks/cvm_cramer_von_mises.py `_SERIES_Z_MAX`
- [ ] 3 analyzers/checks/c3_serial_copula.py `r_library_paths()` / R_LIBS
- [ ] 4 C3 timeout 120 -> 900 and the n^2.8 claim
- [ ] 5 xi_seed at the job layer (2 jobs + jobs/common.py)
- [ ] 6 analyzers/shape_stats.py "234 + 75 = 309"
- [ ] 7 the "RUNS, UNCALIBRATED" doc edits (5 places)

### 1 - leaked `label` (panels/within_calibration.py:676)  VERIFIED
GREEN | panels/within_calibration.py:673-681 | Fix confirmed end to end on the SHIPPED artifact.
`output/t2star_q1_070423_d611af_20260813_153158/t2star_panel_data.pkl` has
`pooled_is_thin_per_threshold` True at exactly 4-9 us (1,2,3,10 us False), and
`pdftotext` on the static PDF of the same run prints
"pooled over windows averaging too few reads to be a within-excursion statistic: 4 us, 5 us,
6 us, 7 us, 8 us, 9 us" - the right six thresholds, named once. | none
GREEN | panels/, plots/ (whole class) | AST scan for any name bound by a `for` target and read
after the loop body ends, comprehension scopes excluded, over every function in `panels/*.py`
and `plots/*.py`: ZERO hits. The only two candidates (`non_repairable.py:553,579`) are
comprehension-local rebinds, not leaks. The class is closed, not just the instance. | none
MINOR | panels/within_calibration.py:676-678 | `zip(labels, notes)` without `strict=True`. The two
lists are appended in lockstep in one loop so they cannot differ today, but a future edit that
`continue`s the build loop would silently truncate the caveat again - the same silence the fix
just removed. | `zip(labels, notes, strict=True)`.

### 2 - CvM `_SERIES_Z_MAX = 5.0`  ATTACKED, mostly holds
GREEN | analyzers/checks/cvm_cramer_von_mises.py:88-131 | The cutoff discards NO resolution
below it and the divergence claim reproduces. Measured: the 12-term sum is bit-comparable to a
40- and an 80-term sum through z = 5 (identical to all 15 printed digits at z = 0.05 ... 5.0),
so nothing in [2, 5] is lost - the series IS evaluated there (F(2) = 0.999987219, p = 1.278e-05;
F(5) = 0.999999999996944, p = 3.06e-12). With the cutoff removed the raw series bottoms at
z = 7.60 (p = 1.0e-15) and climbs: p(20) = 1.19e-08, p(100) = 3.76e-03, p(500) = 0.100,
p(1e5) = 0.705 - the ledger's numbers, confirmed independently. The function is monotone
non-increasing in p over 25000 points spanning [1e-3, 1e8], the discontinuity at the cutoff is
3.1e-12, and vectorised output is BIT-IDENTICAL to the scalar loop (max|diff| = 0) including an
array that straddles the cutoff. `np.where(inside, z, 1.0)` leaks nothing: with
`warnings.simplefilter("error")` no warning is raised at any z from 1e-300 to 1e300. | none
GREEN | tests/test_checks_cvm.py:151-190 | The new tests DO fail when the guard is defeated.
Runtime-patched `_SERIES_Z_MAX` (no source edit): at 1e12 both
`test_the_p_value_is_monotone_...` and `test_the_reviewers_verdict_flip_case_now_rejects` fail;
at 50 and at 20 the monotonicity test still fails. | none
MINOR | analyzers/checks/cvm_cramer_von_mises.py:214 (`run`) | A SATURATED p is reported as
`p_value = 0.0` with `notes = "gamma=... limiting_CvM"` - indistinguishable in the result from a
p that was computed. The truth is a bound (`p < 3e-12`), not a number, and this repo's own rule
is that a p which was not computed must be distinguishable in the result. Latent only because
CvM is unregistered in `battery.ROW_KEYS`; the day it is promoted a ledger row will read a flat
0. | When `observed > _SERIES_Z_MAX`, append e.g. `" p<3e-12 (series saturated)"` to `notes`.
MINOR | tests/test_checks_cvm.py:163 | The monotonicity test certifies "no verdict-flipping
divergence", NOT the value 5.0. Measured: with `_SERIES_Z_MAX` raised to 8, 10 or 12 the whole
file still passes 20/20, because the rise between z = 7.6 and z = 12 is ~1e-14, under the test's
own `<= 1e-12` slack. The band (5, ~15] is uncertified; only the comment's three lines of
evidence hold the constant at 5. | Either assert the constant directly, or say in the test that
it pins the divergence and not the cutoff.
MINOR | analyzers/checks/cvm_cramer_von_mises.py:71-72 | Un-fixed carry-over from the cold-read
pass (ledger:384), now WORSE: "the j = 3 term is already below 1e-12 across the range that
matters" sits four lines above the constant that DEFINES the range that matters as z <= 5, where
the j = 3 term is 1.949e-03 - nine orders out. (Convergence at z = 5 is fine on the real
evidence: the j = 11 term is 4.1e-24. The sentence, not the code, is wrong.) | Restate as "at the
cutoff z = 5 the last kept term is 4e-24; the j = 3 term is below 1e-12 only for z <~ 1".

### 3 - C3 `r_library_paths()` / `R_LIBS`  WORKS, with three caveats
GREEN | analyzers/checks/c3_serial_copula.py:70-92, :133-134 | The mechanism is correct and
measured. `r_library_paths()` returns
`('/home/sera/R/library','/usr/lib64/R/library','/usr/share/R/library')`; under `--vanilla`
`.libPaths()` is the last two and `requireNamespace('copula')` is FALSE, with `R_LIBS` set to the
joined tuple it is all three and TRUE. `R_LIBS` is the right variable (it is read from the
process environment, which `--no-environ` does not suppress - only `.Renviron` files are; and it
takes a path LIST, which `R_LIBS_USER` semantically does not). `os.pathsep` matches R's
separator on both platforms. C3 then runs end to end here: n=50 -> statistic 0.00579 in 1.8 s.
| none
MINOR | analyzers/checks/c3_serial_copula.py:80-92 | One extra `Rscript` spawn PER CALL, measured
at 0.28-0.31 s. Negligible against a 60-130 s R run, and a ledger pass adds ~6 s over 20 calls -
but it is a repeated probe of an environment that cannot change inside a run. | `functools.cache`
on `r_library_paths`, or hoist it to the caller.
MINOR | analyzers/checks/c3_serial_copula.py:134 | The empty-tuple case sets `R_LIBS=""`. Measured
harmless (an empty `R_LIBS` behaves exactly as unset: same `.libPaths()`, `copula` FALSE), and
unreachable anyway because `_invoke_rscript` raises when `rscript_path()` is None. It does mean
that if R ever printed nothing the fix silently reverts to the broken state rather than saying
so. | Skip the assignment when the tuple is empty, so the variable is absent rather than blank.
MINOR | analyzers/checks/c3_serial_copula.py:14, :123 and docs/iid_checks/C3_serial_copula.md:4 |
`/home/sera/R/library` is hard-coded in PROSE in three places - the exact machine-binding
`r_library_paths()` was written to avoid, now in the docstring instead of the code. Nothing leaks
into `output/`: the success-path `notes` is `lag.max=... N=... seed=...` and `R_LIBS` is never
recorded anywhere (see the next finding). The one real leak route is a FAILURE: `check_ledger`
writes `f"declined: {exc}"[:160]`, and that exception carries R's stdout/stderr, which can carry
absolute paths into the ledger artifact. | Say "the user library" and drop the literal path.
MINOR | analyzers/checks/c3_serial_copula.py:125-131 (docstring) | "The path is then a recorded
input instead of an ambient one, which is the same reason seeds are explicit here" is FALSE as
written: the path is passed explicitly to the subprocess but RECORDED nowhere - not in
`CheckResult.notes` (`lag.max=... N=... seed=...`), not in the ledger columns, not in any prov
record. Grepping `output/` for `R/library` finds nothing. It is explicit, not recorded. | Either
put the library paths (and the R version) in `notes`, or drop the "recorded" claim.

### 4 - timeout 120 -> 900 and the n^2.8 claim
GREEN | analyzers/checks/c3_serial_copula.py:182 with check_ledger.py:408 | The new default is
the one that binds: `check_ledger` calls `c3.run(...)` without `timeout_s`, so 900 s applies.
Timing re-measured here from scratch (idle machine): 1.8 s at n=50, 6.0 s at n=150, 60.3 s at
n=355 - about 2.2x faster than the docstring's machine, same shape. | none
IMPORTANT | analyzers/checks/c3_serial_copula.py:186-187 | "The cost grows roughly as n^2.8" is
NOT supported by the three points the same sentence cites. From (50, 3.9), (150, 14.7),
(355, 130.2): pairwise exponents are 1.208 (50->150), 2.532 (150->355) and 1.790 (50->355), and
the OLS log-log slope over all three is 1.764. My own three points give the same picture (1.10 /
2.68 / 1.79). Nothing in the cited data yields 2.8; the only evidence pointing that steep is the
UNFINISHED 682 run, which bounds the local exponent below at 2.29 and is not a measurement. The
"~800 s at n = 682" extrapolation two lines down is 2.8's output - at the fitted 1.76 it is 411 s
and at the top-end local 2.53 it is 679 s. This is the docstring-number rule at :186 exactly. |
State "superlinear; locally ~2.5 between n = 150 and 355, 1.8 fitted across the measured range",
and give the 682 extrapolation as a range with its exponent named.
IMPORTANT | analyzers/checks/c3_serial_copula.py:13-18 (the FIRST CONTACT numbers) | Two of the
three quoted p-values do not reproduce, and the docstring never says at which seed. Re-run with
the module's own default `seed=0` on the same inputs (confirmed same inputs: all three
STATISTICS reproduce to 5 digits - 0.00579, 0.00713, 0.00763 - once the generator stream is
continued rather than reseeded): p = 0.9476 vs the quoted 0.958 at n=50, p = 0.9116 vs 0.904 at
n=150, p = 0.8656 vs 0.866 at n=355. The statistic is seed-INVARIANT and the p is not: at n=50,
seeds 0/1/7/20260813 give 0.9476/0.9575/0.9535/0.9505. So the quoted p-values came from other
seeds, and the C3 p carries ~0.01 of Monte-Carlo noise at N=1000. | Quote the seed beside the
numbers (the repo's own rule: a seed that affects a reported number is a declared parameter), or
quote the statistics only, which do reproduce.
MINOR | analyzers/checks/c3_serial_copula.py:182 with check_ledger.py:407 | 900 s x up to 20 C3
invocations is a 5-hour worst case before a ledger job gives up, each ending in a `declined:` row
that reads `not computed`. The budget is right for one call and unbounded for the pass. | Note
the per-pass worst case in the docstring, or give the ledger its own smaller budget.

### 5 - `xi_seed` declared at the job layer
GREEN | jobs/active/t2star_q1_070423.py:102-117, twin, jobs/common.py:258-262 | Verified end to
end, not by reading. Both shipped prov.md files carry
`step_t2star_panel_data["t2star_panel_data(shape_min_reads=5, use_uncertainty=True,
xi_seed=20260813)"]`, and `python main.py inspect jobs/active/ramsey_q1_100423.py` shows
`fidelity_panel_data: _fidelity_panel_data <- [fidelity_raw, fidelity_windows]
{'xi_seed': 20260813}`, so the fidelity branch will label it too. Consumer coverage is complete:
`build_within_calibration_panel_data` has exactly two callers (`t2star.make_panel_data`,
`fidelity.make_panel_data`) and both now receive the seed from a job; `mtbf_q1` builds a
`AcrossCalibrationPanelData`, which has no xi. | none
MINOR | jobs/active/t2star_q1_070423.py:116, t2star_q1_100423.py:116, jobs/common.py:196 | The
literal `20260813` is written in three places with nothing tying them together. Nothing is wrong
today; the drift mode is that one job is re-seeded and the composite
`compare_t2star_0704_vs_1004` then compares a permutation-calibrated number against another
seed without saying so. | One module-level constant imported by all three, or read it from
`jobs/common`.
MINOR | jobs/active/t2star_q1_070423.py:116 with t2star_q1_100423.py:116 | The two independent
datasets now share a seed, so on tied data windows of equal `n` in the two records would draw a
BIT-IDENTICAL permutation sequence and their Monte-Carlo errors would be common - the same defect
already logged for windows WITHIN a record (this ledger:352), now across the two records the
composite job differences. Dormant today and verified so: 234/234 and 75/75 rows read
`xi_p_method = asymptotic`, so the seed is never consulted on this data. | Same fix as :352 -
derive the per-window seed - or say in the composite that the two share it.
IMPORTANT | jobs/common.py:258-262 | The fix is applied to ONE of the parameters that shape the
fidelity panel's numbers. The `fidelity_panel_data` step declares `xi_seed` and NOTHING else, so
`shape_min_reads` (5), `k` (1.0) and `use_uncertainty` all still ride the builder defaults and
are absent from the Mermaid label - confirmed by `inspect`, whose kwargs dict is
`{'xi_seed': 20260813}` alone. The t2star jobs two files over declare `shape_min_reads` and
`use_uncertainty` for exactly this reason, and the comment added here points AT that file. This
is the "fix every site of a class" rule, half applied inside the function being fixed. | Declare
`shape_min_reads` (and `use_uncertainty`) on the fidelity step too.
MINOR | jobs/common.py:180 | An unannounced behaviour change rode along with the seed fix:
`use_uncertainty=bool(window_result.meta.get("use_uncertainty", False))` is new and is in none of
the seven fix entries. It is a no-op today (the fidelity carve never passes `use_uncertainty`, so
the meta value is False), but `.get(key, False)` is a silent fallback on a key `windows.py:526`
writes unconditionally - if it ever went missing the panel would quietly draw the certain-state
figure over an uncertain carve. | `window_result.meta["use_uncertainty"]`, and record the change.

### 6 - "0 of 279" -> 234 + 75 = 309  VERIFIED IN THE ARTIFACTS
GREEN | analyzers/shape_stats.py:84-91, :193-194 | Both counts come from the artifacts that ship.
`output/t2star_q1_070423_d611af_20260813_153158/t2star_panel_data.pkl` has exactly 234 rows
across `shape_stats_per_threshold`, the 1004 twin (731c02) has exactly 75; 309 total. The tie
claim is confirmed by TWO independent routes in the artifact: `xi_tie_frac` (the x axis) is 0.0
on all 309, and `xi_p_method` is `asymptotic` on all 309 - which, since the branch is
`max(tie_fraction(x), tie_fraction(y)) > 0.0`, is a recorded proof that the Y axis is tie-free on
every one of them too. Naming the two datasets separately is the right correction. | none
MINOR | analyzers/shape_stats.py:11 and :30 | Left un-fixed, in the SAME docstring whose other
number was corrected, and it is the same defect class: "the median window here has ~7 reads",
twice, and it carries the argument that the asymptotics are marginal. Measured on the two shipped
artifacts: median `n_reads` over the rows that actually reach these statistics is 11.0 (0704) and
14.0 (1004), minimum 5 on both. 7 is the median at the 3 us threshold alone (the cold-read pass
said so at this ledger:376); it appears in neither the "fixed" nor the "carried to Increment C"
list. | Same treatment as the 309: name the population, "median 11 reads over the 234 eligible
windows on 0704, 14 over the 75 on 1004".

### 7 - "R is absent" -> "RUNS, UNCALIBRATED"  INCOMPLETE AND SELF-CONTRADICTORY
IMPORTANT | bench/results/promotion_report.md:5 | A SIXTH site, uncorrected, and it is the one
TRACKED doc of the set: "C3 (`copula::serialIndepTest`) needs R, and `Rscript` is absent on this
machine ... no code path past `_invoke_rscript` has ever executed." Both clauses are false as of
2026-08-13 (Rscript 4.5.3 at /usr/bin/Rscript; C3 ran here at n = 50/150/355 during this audit).
Its CONCLUSION - "C3 carries no evidence here. Its silence is not a pass." - is still exactly
right, which is why this is a doc defect and not a wrong verdict. The fix entry claims five
sites; the class has six. | Restate the environment clause; keep the conclusion verbatim.
IMPORTANT | docs/iid_checks/C3_serial_copula.md:35-41 vs :3-8 | Contradictory statements in one
short file, which is the specific thing this pass was asked to check. The header says C3 has run;
the "Limitations" section still says, in the present tense, "It cannot pin the numeric path. The
first machine with R installed gets a code path that has never run, including the CSV round-trip
and the result parsing." A reader who starts at "## Limitations" - which is where anyone
assessing the check starts - gets the pre-2026-08-13 state as current. :55 compounds it: "report
`not computed` when R is absent - which is every row in every ledger produced so far". | Move
those two paragraphs under the same "historic" framing the header paragraph already uses.
MINOR | docs/iid_checks/C3_serial_copula.md:10-11 | The edit left the markdown broken: "...had
ever executed. Its silence is not a pass.**" has a closing `**` with no opener, so every
subsequent bold marker in the rendered file is inverted, and the orphaned sentence no longer has
a subject. | Delete the stray `**`, or restore the pair.
MINOR | docs/iid_checks/C3_serial_copula.md:59-60 | "If R is installed later:
`test_r_is_genuinely_absent_here` skips rather than fails" is written as a future contingency for
something that has happened. Verified: `pytest tests/test_checks_c3_bridge.py` = 5 passed, 1
skipped, with exactly that skip message. | Present tense.
MINOR | docs/iid_checks/C3_serial_copula.md:3-6 and LIMITATIONS.md:44-46 | The unreproducible
smoke-test p-values (0.958 at n=50, 0.904 at n=150) and the unsupported `n^2.8` are each
replicated into a second doc, so the two findings under fixes 3-4 above have three sites apiece,
not one. | Fix the module docstring and propagate.
MINOR | docs/iid_checks/README.md:90 | The cold-read pass asked for the status ROW to be dropped
(CLAUDE.md: no status tables in agent-facing docs, "it went stale in one day"); the fix updated
the row instead. The content is honest - "no bench cell, no size or power evidence" - but
"first contact 2026-08-13 (Rscript 4.5.3)" is environment state inside a table whose other four
rows are bench VERDICTS, and it is the part that will go stale next. | Keep the verdict, move the
date to the C3 page.
MINOR | docs/iid_checks/LIMITATIONS.md:40 | The heading still reads "## 4. C3 is unassessed"
while its own first line says C3 executes. Defensible (unassessed = uncalibrated) but it is the
one line a skimmer reads. | "C3 runs but is uncalibrated".
GREEN | analyzers/permutation.py:49-54, :93-100 | The "modal case" reword is now TRUE and is the
right claim. Verified against both shipped artifacts: minimum `n_reads` is 5 on each (234 rows on
0704, 75 on 1004), so 0 rows can reach the enumerated regime, and `(n!)^2` = 36 / 576 / 518400 at
n = 3 / 4 / 6 as stated. The docstring keeps the field and says why a direct caller still needs
it. | none
GREEN | no doc implies C3 carries calibration evidence | Checked every site: README.md:90
("no bench cell, no size or power evidence"), LIMITATIONS.md:44-48 ("What remains missing is
CALIBRATION ... No conclusion anywhere rests on C3"), C3_serial_copula.md:3-8 ("a smoke test, NOT
calibration"), the module docstring, and promotion_report.md:5 ("C3 carries no evidence here").
The smoke test is never described as evidence of size or power anywhere. | none
GREEN | baseline reproduced | `PYTHONPATH=. python -m pytest tests/ -q` = 253 passed, 2 skipped.
`ruff check .` = the single known pre-existing F841 (`raw_meta`) and nothing else. No source file
was modified by this audit; the only write is this ledger. | none

MANIFEST CLOSED: all seven audited. 1 verified clean, 2 holds with 3 minors, 3 works with
4 minors, 4 has an unsupported exponent and unreproducible p-values, 5 is half-applied in
jobs/common.py, 6 verified in both artifacts, 7 is incomplete (a sixth site) and
self-contradictory inside one file.

---

## Increment B - fixes applied after the AUDIT of the seven fixes

The audit's own summary: five of seven fixes verified correct by measurement; the defects it
found were almost all in the PROSE describing them, which is the failure mode CLAUDE.md's claims
discipline names. Fixed:

FIXED | analyzers/checks/c3_serial_copula.py | "the cost grows roughly as n^2.8" was INVENTED.
Refit from the three cited points: OLS log-log slope 1.764, with pairwise slopes that disagree
with each other (1.208 from 50->150, 2.532 from 150->355, 1.790 end to end). Three points do not
determine a power law, and the separate n=682-unfinished-at-580s observation is not predicted by
any fit through them (slope 1.76 extrapolates to 412 s). No exponent is claimed now; the three
measured points are given, the disagreement is stated, and the reader is told to measure the n
they have. This is the third time this pass that a number was asserted rather than computed.

FIXED | analyzers/checks/c3_serial_copula.py, docs/iid_checks/C3_serial_copula.md | The quoted
C3 p-values (0.958 / 0.904 / 0.866) were produced at `seed=1`, and quoted without it; at the
module default `seed=0` the same inputs give 0.948 and 0.912. The STATISTIC is seed-invariant
and reproduces exactly - the R side simulates its own null, so only the p moves. Both sites now
quote the seed and say which half depends on it.

FIXED | bench/results/promotion_report.md:5 | The SIXTH "R is absent on this machine" site, and
the only TRACKED one - the ledger's pass rule cites this file, so it was the one that mattered
most. Now states that C3 has executed, that the smoke test changes nothing in the report because
the bench never ran C3 and `ROW_KEYS` has no C3 row, and that failing to reject iid data is not
calibration.

FIXED | docs/iid_checks/C3_serial_copula.md | Said C3 had run at line 3 and that the numeric path
"has never run" under `## Limitations`, and my edit left an unmatched `**` that bolded a whole
paragraph. Both repaired; marker count now even. The historic note is kept but labelled as
history rather than left contradicting the header.

FIXED | analyzers/checks/cvm_cramer_von_mises.py:71 | Carried-over comment "the j = 3 term is
already below 1e-12 across the range that matters". Measured: 8.4e-28 at z=0.35, 9.2e-11 at z=1,
but 1.9e-03 at z=5 - nine orders out at the top of the retained range. 12 terms is still right
(the audit confirmed 12-, 40- and 80-term sums agree to every printed digit through z=5), but the
margin comes from the later terms, and the comment now says so.

FIXED | analyzers/checks/cvm_cramer_von_mises.py | A saturated `p_value = 0.0` reached
`CheckResult` indistinguishable from a computed one. The verdict is right and the magnitude is
not a measurement, so `notes` now carries `p_saturated(z>5)` and a reader cannot quote the zero
as a precision claim.

FIXED | analyzers/shape_stats.py:11,30 | "the median window here has ~7 reads", twice. Measured
over the eligible rows of the shipped artifacts: 11.0 on 0704 and 14.0 on 1004. The 7 is the 3 us
cell - a pooled number labelled as though it were the whole, which is the FIRST rule in CLAUDE.md's
claims discipline and the fourth instance this pass.

FIXED | jobs/common.py + both t2star jobs | The literal `20260813` was in three files. Now one
`XI_SEED` in `jobs/common.py`, imported. The constant's comment carries the open finding it does
NOT fix: one seed serves every window, so two windows of equal n draw an identical permutation
sequence and their Monte Carlo error does not average down under the median.

FIXED | jobs/common.py:181 | `use_uncertainty=bool(window_result.meta.get("use_uncertainty",
False))` - a silent fallback that rode along with the seed fix and would redraw the panel in the
other mode if the key were ever missing. `windows.run` always records it, so absence is a broken
artifact. Now indexes and raises.

FIXED | panels/within_calibration.py | `zip(labels, notes, strict=True)`.

VERIFIED AFTER THE FIXES | Re-ran the 0704 job (identity 1e0202). The thin-pooling caveat prints
"4 us, 5 us, 6 us, 7 us, 8 us, 9 us" on the static PDF - exactly the six thresholds the artifact
marks thin, where before the leak it printed nothing. Label carries `xi_seed=20260813`.
253 passed, 2 skipped; ruff at the one known pre-existing F841.

NOT FIXED, carried to Increment C, with the audit's reasons:
- The new CvM tests pin the DIVERGENCE, not the value 5.0: raising `_SERIES_Z_MAX` to 8, 10 or 12
  still passes all 20. That is arguably correct (the cutoff is a judgement, the monotonicity is
  not), but it means the constant can drift without a test objecting.
- `r_library_paths()` spawns an extra Rscript per call (0.29 s) and its result is written to no
  provenance artifact, so two machines with different library sets are indistinguishable after
  the fact.
- `jobs/common.py`'s fidelity step declares only `xi_seed`; `shape_min_reads`, `k` and
  `use_uncertainty` remain on builder defaults and off the Mermaid label - the same class as the
  fix, in the function the fix touched.
- `k` is declared on the windows step and never reaches the panel adapter, which uses its own
  `k=1.0`. Pre-existing.
- `distinguish_band.py:158` still claims the field reaches the provenance label from the analyzer
  layer; `DistinguishBand` does not carry `xi_seed` as a field.

UNREVIEWED LAYER | the fixes made in response to the AUDIT | Recorded so it is not later assumed
covered. The review chain is three deep and only two links are reviewed: Increment B (cold read +
mutation pass), then the seven fixes (the audit), then THESE, which no reviewer has seen. Six are
prose corrections written directly against the audit's measured objections; four are code -
`p_saturated` in CvM notes, `zip(strict=True)`, the `XI_SEED` extraction, and `use_uncertainty`
going from `.get(..., False)` to indexing.
Exercised since: 253 passed; both t2star jobs re-run; `ramsey_q1_100423` re-run, which is the ONLY
caller of `_fidelity_panel_data` and therefore the only exercise of the raise - it completed and
its label reads `fidelity_panel_data(xi_seed=20260813)`. The "windows.run always records
use_uncertainty" claim was VERIFIED rather than left asserted: `WindowsResult` has exactly one
construction site, analyzers/windows.py:519-529, and it always sets the key.
To fold into the Increment C gate. Worth a specific look: whether these corrections now UNDERstate
in the other direction - three consecutive audits pushing on overclaiming can hedge prose into
uselessness, which is its own failure.

---

## Increment C - Tier 3 (permutation exactness and asymptotic size)

NEW | tests/test_tier3_calibration.py | 12 tests. Two claims held to two standards on purpose:
C5/C6 are permutation-calibrated, which is EXACT under exchangeability, so their level is a
combinatorial identity and is ASSERTED; C1/C2/CvM are asymptotic and divided by an estimated
gamma_hat, are known to miss at small n (bench: 0.069 and 0.0757 at n=20 vs 0.05), and are
MEASURED with only the direction asserted. Asserting a tight band on the latter would encode
the miss as correct.
Not a KS test, per the plan: at B=199 a permutation p lives on a 200-point grid and differs
from the continuous uniform by 1/(2B) BY CONSTRUCTION, so KS would test the grid. The exact
statement `P(p <= k/(B+1)) = k/(B+1)` is used instead, at alphas chosen to LAND on the grid
(0.01=2/200, 0.05=10/200, 0.10=20/200).
Aggregation and multiplicity stated before the verdict: 2 checks x 3 layouts x 3 alphas = 18
comparisons against one band of `bonferroni_z_crit(18)` null SEs, reusing
`calibration_summary.null_se` and `bonferroni_z_crit` so this file and the bench keep ONE
definition of calibrated. Resolution stated rather than hidden: +/-0.020 at alpha=0.10 and
+/-0.0067 at alpha=0.01, so the 0.01 cell catches only a gross error.

FOUND BY MUTATION, IN MY OWN NEW TEST | The first version of this file could not catch two
defects in `permutation_p_value`. Both mutants survived all 11 tests:
  - dropping the `+1` tie correction: at the grid alphas `k/199 <= 0.05` and `(1+k)/200 <= 0.05`
    BOTH reduce to `k <= 9`. The rejection rates are IDENTICAL at 0.01, 0.05 and 0.10 - an
    arithmetic coincidence, not a weak assertion. The term is visible only at the FLOOR, where
    the correct minimum is 1/(B+1) = 0.005 and the mutant returns exactly 0.
  - `>` instead of `>=`: on exponential gaps no permuted statistic ever exactly equals the
    observed one, so the substitution is a no-op across every replicate.
Two tests added that look where the terms actually live (the observed floor over 2000
replicates, and a direct unit test with exact ties in the null). All three mutants - including
`divide by B not B+1` - are now CAUGHT. This is the same lesson as the 2/5 null variance in
Increment B: a calibration test can pass while the constant it depends on is unpinned.

---

## Increment C - Tier 4 (R cross-implementation)

NEW | rscripts/reference_values.R + tests/fixtures/r_reference_{inputs,values}.csv +
tests/test_r_cross_implementation.py | 18 tests, and the suite still never requires R: the
script writes committed fixtures, the test reads them. CSV rather than the planned JSON because
`jsonlite` is not installed and long-format CSV diffs cleanly under review. The fixture carries
the INPUT DATA as well as the answers - matching R's RNG stream from Python would be fragile and
a generator mismatch would surface as a fake statistical disagreement.
Agreement at 1e-10 or better on: xi tie-free vs `XICOR::xicor`, xi with y-ties-only vs the same,
`energy::dcor` on three cases, Spearman on four, and the lag-1 rank autocorrelation C5 is built
on. R 4.5, XICOR 0.4, energy 1.7, randtests 1.0, copula 1.1, all recorded in the fixture.

THE FINDING | `XICOR::xicor` is RANDOM on tied-x data, so equality is the WRONG test | The
headline tied case failed at 0.4649 (ours) against 0.4976 (R), which looks exactly like a
transcription error. It is not. Eq (8) breaks x-ties uniformly AT RANDOM: measured, 7 distinct
values in 8 calls on identical data, and over 2000 draws mean 0.4825, sd 0.0429, range
[0.357, 0.628]. A single R draw is not a reference constant and pinning it would pin R's RNG
state. Our deterministic stable-sort value sits at 0.4649 - inside the 0.1-99.9 percentile band
and within 3 sd of the mean.
The test now asserts MEMBERSHIP in the tie-break distribution rather than equality, which is a
weaker claim than the one first written and the strongest claim the estimator supports. The
mechanism is localised, not guessed: the y-ties-only case has continuous x, `xicor` is
deterministic there (50 draws, 1 distinct value), and our value matches it to 1e-10. So the
disagreement is X-tie-BREAKING, not eq (8).
This is also the concrete motivation for the tie experiment: what the deterministic break COSTS
is now a measurable question with a known reference distribution, not a caveat.

BARTELS, deliberately not used as a target | `randtests::bartels.rank.test` is the rank von
Neumann RATIO; our C5 is a studentized MAXIMUM over lags. Cross-checking one against the other
would be apples to oranges. What they genuinely share is the lag-1 rank autocorrelation, so that
is what is pinned; the Bartels values are in the fixture as context and are asserted only to
exist. Recorded so a later reader does not "fix" the missing comparison.

---

## Increment C - the xi tie experiment and the instrument report

NEW | bench/xi_ties.py -> bench/results/xi_tie_experiment.csv | Three questions on one
generator: where the closed form stops agreeing with the permutation, what the deterministic
x-tie break costs against a randomised one, and whether the asymptotic null survives ties.
A study, not a pipeline layer - nothing outside bench/ imports it.

METHOD DEFECT FOUND AND FIXED BEFORE THE RUN | The deliverable was specified as "the tie
fraction where the two p-values diverge by more than 0.02", read off `mean |p_closed -
p_perm|`. That statistic cannot answer it: at B = 999 the permutation p-value carries Monte
Carlo error of about 0.0158, which is MOST of a 0.02 threshold, so the absolute difference
measures mostly its own resampling noise. Confirmed on the tie-free control cell, where the
true difference is ~0 and |dp| still reads 0.0120. The experiment now records the SIGNED
mean as well - noise averages out of it, systematic disagreement does not - and the
deliverable is read off the signed quantity with the MC floor stated beside it.

SECOND METHOD DEFECT, found by reading the first four cells | Q2 asks what the deterministic
X-tie break costs, and reported a spread of exactly 0.0000 in every cell. The generator only
quantised y; x stayed continuous, so there were never any x-ties for a tie break to act on.
A column of zeros would have read as "the deterministic choice costs nothing". Fixed: x is
quantised in the Q2 block only (Q1/Q3 keep x continuous, which is the shipped situation -
x is the read time). Now measures a real effect: at n=35 with 20 levels, x tie fraction 0.80,
tie-break sd 0.0456, and our deterministic value at the 82nd percentile of the admissible
ones; at 5 levels, sd 0.1719.

EARLY RESULT, full grid still running | Across n=35 the signed difference stays within
-0.0026 to -0.0011 at tie fractions from 0.000 to 0.994 - the closed form and the permutation
do NOT systematically disagree even at 99% ties. What DOES move is the closed form's type-I
error: 0.027 (tie-free) -> 0.043 -> 0.057 -> 0.063 (tie_y 0.994) against a nominal 0.05. So
the cost of ties shows up as level, not as p-value disagreement. If that holds across the
grid it argues XI_TIE_CUTOFF = 0.0 is far more conservative than it needs to be - but the
grid is not finished and this is written as the provisional reading it is.

NEW | analyzers/instrument_validation.py, plots/instrument_validation_plot.py,
jobs/active/instrument_validation.py, bench/instrument_report.py, tests (11) | The report.
Four figures from one typed artifact: the tier matrix (the index), published values (Tier 2),
cross-implementation (Tier 4), and the tie experiment. Pure-compute builder, pure renderer,
every input declared as a Dataset so the figure's dependence runs through provenance. The
`instrument | tier | verdict` table is GENERATED markdown, not drawn into a PNG - a PNG goes
stale as silently as a doc table and cannot be diffed.

STRUCTURAL | reference/ is new and TRACKED | The Tier 4 fixture was written to
tests/fixtures/, and `tests/` is gitignored - so the "committed fixture" the plan requires
would not have survived a clone, and the suite would have silently started requiring R to
regenerate it. The published load-haul-dump record had the same problem and is now also a
tracked CSV, read by BOTH the suite and the figure so the two cannot drift. CLAUDE.md's
Layout section updated.

NEW | jobs/active/check_ledger_q1_070423.py | The "our data" counterpart: the validated
instruments run on the real 912-day record, same carve as the T2* panel so the two are
comparable. No job wired the ledger before this.

OPEN, and a real one | C3 may make the ledger impractical. Its cost grows steeply with n and
the job runs it per threshold, so a 10-rung ladder can spend a long time in R before every
row either lands or times out at 900 s. `not computed` is the honest outcome and cannot be
misread as a pass, but if most C3 rows time out then paying the wall-clock buys nothing and
the ledger should either take a smaller timeout at the job level or skip C3 above some n.
Measure before deciding - this is written as the open question it is.

---

## Increment C - RESULTS (both studies finished)

THIRD METHOD DEFECT, caught by reading the printed output | bench/xi_ties.py main() | The
signed-mean rewrite of the DELIVERABLE block never applied - the formatter had reflowed it and
the string replace failed silently, the same way three earlier edits in this pass did. The run
therefore printed the deliverable off `|dp|`, the statistic already established as unusable.
Recomputed from the CSV, which carries both columns, and the two answers DIFFER: `|dp|` puts
the crossing at 3 levels for n=35, the signed mean at 2. The tie-free control settles which is
right - its true difference is ~0 and `|dp|` still reads 0.0095-0.0120, so `|dp|` was reporting
its own noise. Code fixed; no re-run needed.

DELIVERABLE | The closed-form and permutation p-values agree until the response is nearly
binary. Systematic divergence beyond 0.02 (signed mean) appears ONLY at 2 distinct response
values, at every n: -0.0552 (n=35), -0.0359 (n=75), -0.0261 (n=150), -0.0203 (n=355). Reported
in LEVELS rather than tie fraction because tie fraction saturates - 20 levels already ties
83-99% of a sample, so "diverges at tie fraction 1.000" is true of several cells at once.

THE MORE USEFUL RESULT | The cost of ties is LEVEL, not p-value disagreement. Type-I of the
closed form at nominal 0.05, by distinct response values:
    levels     20     10      5      3      2
    n=35    0.043  0.057  0.063  0.103  0.163
    n=75    0.033  0.047  0.063  0.087  0.147
    n=150   0.040  0.047  0.073  0.077  0.173
    n=355   0.053  0.063  0.087  0.087  0.140
The permutation reference holds 0.023-0.063 in EVERY cell. So the closed form is sound down to
about 5 distinct values and inflates to 3x nominal at 2-3. This is direct evidence that
`XI_TIE_CUTOFF = 0.0` - permute on ANY tie - is far more conservative than the data requires,
and it costs 760x per panel render. A cutoff keyed on DISTINCT VALUES rather than tie fraction
would be the principled form. NOT changed in this pass: it is a calibration decision, the
measurement is one generator on simulated data, and the estimator's real ladder is tie-free
anyway so nothing ships differently today.

Q2, the deterministic x-tie break | Spread of xi over random tie breaks: sd 0.019 to 0.172,
largest at coarse quantisation and small n, exactly 0 when x is continuous. Our deterministic
choice lands between the 2nd and 90th percentile across cells with no systematic bias - it is
an arbitrary admissible member, which is what the Tier 4 membership test asserts.

C3 CONCERN RESOLVED | The open worry that C3 would make the ledger impractical does NOT hold.
The real run completed and C3 returned p-values on every threshold that had usable windows
(.004, <.001, <.001, <.001, .016 in-spec; .004, .371, .648, .460, .480 calendar). No row timed
out at 900 s. The concern is withdrawn.

THE REAL-DATA RESULT | 104 ledger rows: 31 fail, 29 underpowered, 26 not computed, 12 pass,
6 not interpretable (ties). The iid assumption is REJECTED on the in-spec clock at 3-6 us by
every check that had power there, while the CALENDAR clock passes at 4-5 us. The two clocks
disagreeing is the finding the ledger exists to surface, and it is visible at a glance.

STATE AT THE C GATE | 298 passed, 2 skipped. ruff at the one known pre-existing F841. Style
ratchet 17 (the new plot module drove it to 28 and was refactored back - all colour and
typography now in plots/theme.py).

---

## Increment C + unreviewed-B-fix-layer - reviewer pass (fresh context, 2026-08-13)

SCOPE: Part 1 the unreviewed B fix layer (cvm p_saturated, zip strict, XI_SEED, use_uncertainty,
six prose corrections); Part 2 Increment C (tier3 tests, Tier 4 R cross-impl, bench/xi_ties.py,
instrument_validation.*, check_ledger_q1_070423.py, reference/, CLAUDE.md Layout);
Part 3 sweep for a fourth silently-failed edit.
COMMIT: 2f1bd95 (working tree, most files untracked)

MANIFEST (checked off as reviewed)
- [ ] analyzers/checks/cvm_cramer_von_mises.py (p_saturated)
- [ ] panels/within_calibration.py (zip strict), jobs/common.py (XI_SEED, use_uncertainty)
- [ ] prose: c3_serial_copula.py, docs/iid_checks/C3_serial_copula.md, promotion_report.md, shape_stats.py
- [ ] tests/test_tier3_calibration.py
- [ ] rscripts/reference_values.R + reference/*.csv + tests/test_r_cross_implementation.py
- [ ] bench/xi_ties.py + bench/results/xi_tie_experiment.csv
- [ ] analyzers/instrument_validation.py + plot + job + bench/instrument_report.py + tests
- [ ] jobs/active/check_ledger_q1_070423.py
- [ ] reference/ + CLAUDE.md Layout

FINDINGS (appended at the moment of discovery)

### Part 1 - the unreviewed B fix layer

GREEN | jobs/common.py:172-198 | `XI_SEED = 20260813` is a single definition imported by both
t2star jobs (no literal left anywhere: `grep -rn 20260813` finds only jobs/common.py). It is
passed as a step kwarg (`xi_seed=xi_seed` at :276) so it reaches the Mermaid label. The
`use_uncertainty` change indexes `window_result.meta["use_uncertainty"]` and the "windows.run
always records it" premise is TRUE: `WindowsResult(` has exactly one construction site
repo-wide (analyzers/windows.py:519) and it always sets the key. | none
GREEN | analyzers/shape_stats.py:11-14,32-34 | The corrected "median window is 11 reads on 0704
and 14 on 1004" is EXACT against the shipped artifacts: median `n_reads` over the concatenated
`shape_stats_per_threshold` rows is 11.0 (n=234) and 14.0 (n=75). Not an understatement - it is
the number the artifact carries. | none
GREEN | panels/within_calibration.py:681 | `zip(labels, notes, strict=True)`; the only other zips in
the file are over a value unpacked from one list (:513) and a `zip(*survival)` transpose (:830),
neither of which strict= applies to. | none

MINOR | analyzers/checks/cvm_cramer_von_mises.py:224 | `p_saturated(z>5)` is appended to the free
text `CheckResult.notes`. That is the right place TODAY (CvM is not in the ledger), but
`check_ledger.py:466` builds its `notes` column from `_verdict`'s reason plus the battery note
and DISCARDS `result.notes` except for a `variant=` token scan. If CvM is ever promoted into
`battery.ROW_KEYS` the saturation marker silently stops reaching the artifact, and a saturated
p=0.0 becomes indistinguishable from a computed one - the exact defect this fix was written to
close. | On promotion, add a structured flag (or scan for the token as `variant=` is scanned).
MINOR | analyzers/checks/c3_serial_copula.py:13 | The line is 130 chars where every neighbour is
<=91, and it states the same fact twice: "**C3 RUNS, and is still UNCALIBRATED.** R IS NOW
INSTALLED AND C3 HAS RUN." Signature of text inserted without the surrounding paragraph being
reflowed - the same edit-hygiene class as the three failed replaces this pass. Content is
correct, not understated. | Reflow and drop the duplicate sentence.
MINOR | docs/iid_checks/C3_serial_copula.md:44-45 | "What remains untested is ... behaviour on
data that is NOT iid" is now stale in the mild direction: C3 has since run on the real 912-day
record and returned rejections at 5 thresholds (ledger:993). Still uncalibrated, so the verdict
holds; the sentence just no longer describes the state. | Say "unCALIBRATED on non-iid data";
it has now been exercised there.

### Part 2/3 - Tier 3 (tests/test_tier3_calibration.py)

IMPORTANT | tests/test_tier3_calibration.py:15,168-172,194-195 | THE FOURTH FAILED EDIT, and it
is a claim not a formatter artifact. Three places say the instrument report/figure imports
`measure_asymptotic_size` "so the figure and the test measure the SAME thing by construction".
NOTHING imports it: `grep -rn "measure_asymptotic_size" --include=*.py .` matches only this file,
and `analyzers/instrument_validation.py` has no size field, no import of any check module, and
`plots/instrument_validation_plot.py` draws only the tie experiment's type-I columns. So the
function is not public-for-the-report, and the report's three tier-3 "asymptotic size measured"
verdicts rest on NO number in the artifact. | Delete the three sentences, or have the builder
actually call it and store the sizes.
IMPORTANT | tests/test_tier3_calibration.py:13,189 | "the bench measures it at 0.069 (C1) and
0.0757 (C2) at n = 20" is a POOLED MEAN over the two Weibull shapes, labelled as a measurement.
Verified against bench/results/size_table.csv: at n=20, arm A, continuous, no extra censoring,
C1 asymptotic is 0.064 (shape 0.75) and 0.074 (shape 1.5) -> mean 0.069; C2 is 0.0650 and 0.0865
-> mean 0.07575. This is verbatim the instance CLAUDE.md's FIRST claims-discipline rule records
("At n=20 the two Weibull shapes give 0.0650 and 0.0865"), reintroduced in new prose. |
Quote the two shapes, or say "mean over the two Weibull shapes, spread 0.065-0.087".
MINOR | tests/test_tier3_calibration.py:199 | The asymptotic band `0.005 <= size <= 0.20` is 4x
nominal at the top. Measured actual sizes at n=30, seed 777: it passes with enormous margin (see
the measurement below). The docstring is honest that it catches only a broken limit - but the
instrument report converts this same test into a tier-3 "pass", which is where the width starts
to matter. | State the measured value in the failure message so a drift is visible.

GREEN | tests/test_tier3_calibration.py | Bonferroni family is right: `N_COMPARISONS = 2 x 3 x 3
= 18` is exactly the number of assertions in the exactness test (c5, c6 x 3 layouts x 3 alphas)
and no other assertion shares the band. The alphas DO land on the B=199 grid: p = (1+k)/200, and
0.01=2/200, 0.05=10/200, 0.10=20/200 are exactly representable and exactly equal to the literals
in IEEE double, so `p <= alpha` counts the intended grid points. | none
GREEN | positive control, measured by the reviewer | Monkeypatched three mutants of
`permutation_p_value` into c6 and re-ran the whole family at [30], 2000 replicates: dropping the
numerator +1 gives rates 0.0145/0.0590/0.1115 - INSIDE the Bonferroni band, exactness passes;
`divide by B` and `strict >` also pass exactness. The floor test catches the first two (min p
0.000000 and 0.005025 vs the exact 0.005) and the tie unit test catches `>`. So the ledger's
"all three mutants are now CAUGHT" is TRUE, and its statement that the exactness test alone
cannot see them is also true and now independently confirmed. | none
GREEN | tests/test_tier3_calibration.py:167-181 vs the ledger path | `measure_asymptotic_size`
calls `module.run([segment], calibration=CALIB_ASYMPTOTIC)` with no gamma override, and
`c1.run`/`c2.run`/`cvm.run` all default to `GAMMA_COMPLETE`, which is what `run_battery` passes.
Same code path, m=1. (The ledger's real records are m>1 and the in-spec clock drops C1/C2
entirely, so the measurement is of the instrument, not of the shipped configuration - worth
saying, but it is not a different code path.) | none

IMPORTANT | tests/test_tier3_calibration.py:14, 186-192 | Direction claim contradicted by this
test's own generator. The docstring says the asymptotic checks "miss nominal at small n and miss
upward" and that only "the direction" is asserted. Measured by the reviewer with the file's own
`measure_asymptotic_size` (1200 replicates, seed 777, iid exponential, m=1):
    n=20  C1 0.0342  C2 0.0292  CvM 0.0383     n=30  C1 0.0525  C2 0.0475  CvM 0.0475
i.e. CONSERVATIVE at n=20, the opposite direction from the bench's Weibull cells. Two separate
problems: (a) no direction is asserted anywhere - the only assertion is `0.005 <= size <= 0.20`;
(b) the bench numbers are quoted as if they described this measurement, but they come from a
different generator (Weibull shape 0.75/1.5, bench carve and censoring). | Say which generator
each number came from, and drop "only the direction asserted" or assert one.
IMPORTANT | analyzers/instrument_validation.py:353-358 | The C1 tier-3 row is
`TIER_PASS, "asymptotic size measured, anti-conservative at small n as the source reports"`.
Neither half is supported by evidence in this repo at the place it points: the only tier-3
measurement of C1 (the file the row's own docstring family cites) is CONSERVATIVE at n=20
(0.0342), and "pass" is read off an assertion band four times nominal at the top. A tier whose
question is "does its p-value hold its nominal level" cannot be `pass` and "anti-conservative"
in the same cell. Compare `Chatterjee xi` tier 3, which is correctly `partial` for a weaker
defect. | `partial`, with the measured number and its generator in the detail.

### Part 2 - Tier 4 (rscripts/reference_values.R, reference/*.csv, test_r_cross_implementation.py)

IMPORTANT | tests/test_r_cross_implementation.py:83-113, :20-24 | THE BAND IS TOO WIDE TO BE
EVIDENCE, and the docstring points the reader at it as "the ONLY external evidence that the
[eq (8)] change was made correctly". Attacked as instructed by constructing deliberately wrong
xi estimators on the shipped `xi_tied` fixture (n=60) and testing them against the band
[0.35652, 0.59912] with the 3sd gate [0.3537, 0.6113]:
    tie-free reduction (the PRE-Increment-B code)   0.4815  INSIDE - passes
    ell uses min ranks instead of max               0.4210  INSIDE - passes
    r uses average ranks (the tie-free habit)       0.4744  INSIDE - passes
    r uses min ranks                                0.4838  INSIDE - passes
    ell built from +y instead of -y                 0.4047  INSIDE - passes
    n -> n-1 in the numerator                       0.4738  INSIDE - passes
    factor 2 dropped from the denominator          -0.0702  caught
    unsorted (forgot to order by x)                -0.0478  caught
Only gross errors are caught. Note also that `xicor_q001` EQUALS `xicor_min` in the fixture, so
the lower edge is the whole observed support of 2000 draws, not a percentile. | Keep the test as
the consistency check it is, and move the "only external evidence" sentence to the case that
carries it (below).
GREEN | tests/test_r_cross_implementation.py:136-145 | The evidence the eq (8) change actually
rests on is `xi_tied_y_only`: continuous x, 40 of 45 y-values tied, XICOR deterministic
(`xicor_distinct_draws == 1`, pinned), equality to 1e-10. Verified by the reviewer that this
test CATCHES every mutant the membership test lets through: tie-free reduction 3.6e-02, ell min
ranks 3.1e-02, r average ranks 2.2e-02, r min ranks 4.3e-02, ell from +y 6.5e-02, n->n-1
2.2e-02, against our 2.3e-16. So Tier 4 does discharge the claim - just not where the docstring
says. A revert to the tie-free form is separately caught by :162-170. | none
MINOR | tests/test_r_cross_implementation.py:88 and rscripts/reference_values.R:85-87 | "7
distinct values in 8 calls on identical data, spanning 0.425 to 0.563" is a scratchpad
measurement that appears in no shipped artifact - the fixture records 2000 draws with min 0.3565
and max 0.6284. CLAIMS DISCIPLINE: a measured number in a docstring must come from the artifact
it cites, in the state it ships. | Quote `xicor_sd`/`xicor_min`/`xicor_max` from the fixture, or
drop the 8-call anecdote.
MINOR | tests/test_r_cross_implementation.py:64-69 | `test_the_fixture_records_which_r_produced_it`
asserts `>= 0` on every package major and minor, which no non-negative integer can fail. It pins
that the KEYS exist, not the versions. | Assert the versions the fixture was actually built with
(XICOR 0.4, energy 1.7), so a regeneration under a different version is a visible failure.
GREEN | rscripts/reference_values.R | Reads correctly end to end: every `add_input` case is
consumed by a test, `ties = TRUE` is XICOR's eq (8) path, the tie-free case records BOTH
`ties=TRUE` and `ties=FALSE` and they are bit-identical (0.591553209224785), the C3 reference is
pinned with `set.seed(707)`/N=1000/lag.max=5 and the module docstring's statistic-invariant vs
p-value-varies distinction is respected, and `write.csv(quote=FALSE)` is safe (no field contains
a comma). 18 tests pass with no R present. | none

### Part 2 - bench/xi_ties.py and the tie experiment CSV

GREEN | bench/results/xi_tie_experiment.csv vs the ledger's DELIVERABLE | Every deliverable
number re-derived by the reviewer from the CSV: signed-mean crossings of 0.02 occur ONLY at
y_levels=2, at -0.055185 (n=35), -0.035895 (n=75), -0.026123 (n=150), -0.020296 (n=355) - the
ledger's four figures exactly, and no other cell exceeds 0.02 (next largest |signed| is 0.0116).
The type-I table reproduces cell for cell (0.043/0.057/0.063/0.103/0.163 at n=35, etc.) and the
permutation column spans 0.0233-0.0633, matching "0.023-0.063 in EVERY cell". Q2 spreads span
0.0192-0.1719 as claimed. | none
GREEN | bench/xi_ties.py:147-162 | `_xi_random_tie_break` IS eq (8) with a uniform random x-tie
break: `np.lexsort((rng.random(n), x))` sorts by x with a fresh uniform key breaking ties, `r =
rankdata(y,"max")[order]`, `ell = rankdata(-y,"max")` NOT reordered, denominator
`2*sum(ell*(n-ell))`. Term for term the same expression as `shape_stats.chatterjee_xi`, differing
only in the tie-break, which is what it is for. | none
GREEN | bench/xi_ties.py:205-218 | The Q2 sub-experiment DOES quantise x now
(`_make_pair(..., x_levels=y_levels)`) and only there; Q1/Q3 call `_make_pair` without
`x_levels`, so `tie_fraction_x` is 0.000 in every row while `tie_fraction_x_q2` runs 0.80-1.00.
The column of zeros the ledger describes is gone and the effect is real. | none
GREEN | bench/xi_ties.py main() | The DELIVERABLE block that silently failed to apply HAS landed:
:253-271 reads `p_signed_diff_mean`, states the MC floor, and falls back to reporting the worst
cell when nothing crosses. | none

MINOR | bench/xi_ties.py:10-11 | Module docstring still defines the deliverable as "the tie
fraction at which the two p-values first differ by more than 0.02", which is the form the results
section rejected as unreadable (tie fraction saturates at 1.000 across three cells at once). The
report and `_divergence_levels` both report LEVELS. | Say levels in the docstring too.
MINOR | bench/xi_ties.py:226 | `seed=hash((n, y_levels)) % (1 << 31)` - a seed that determines
every number in a shipped CSV, derived from a hash rather than declared. Int-tuple hashes are
stable across processes so the study IS reproducible, but the seed appears in no column of the
CSV and a reader cannot regenerate one cell. | Put the seed in the row.
MINOR | ledger:989 "no systematic bias" | The 20 tied cells' `xi_deterministic_percentile` has
mean 62.6 and median 77.8 against the 50 an unbiased arbitrary member would give; with SE ~6.5
that is about 2 SE high. The claim is asserted, not measured, and each cell is ONE pair (Q2 uses
a single representative dependent pair per cell, not an average over pairs). | Either measure it
or say "between the 2nd and 90th percentile, on one pair per cell".
MINOR | bench/xi_ties.py:218 | On the tie-free control the "percentile" is 100.0 because every
draw equals ours and the convention is `mean(draws <= ours)`. A point mass has no percentile; a
reader scanning that column sees 100 next to a spread of 1e-16. | Emit NaN when
`xi_tie_break_sd == 0`.
CORRECTION to the entry above (my own arithmetic, recomputed): the 20 tied cells' percentile has
mean 62.78, median 75.25, SE 6.21. The 2-SE-high reading stands; the median figure I first wrote
(77.8) was wrong.

### Part 2 - the instrument report (the tier verdicts)

IMPORTANT | analyzers/instrument_validation.py:369-374 | C2 tier 4 = `pass`, "gamma=1 path
equals the textbook AD form to float precision". The report's own header defines tier 4 as "does
it agree with an INDEPENDENT implementation of the same statistic". The thing it agrees with is
`_classical_ad` at tests/test_checks_statistics.py:63 - a five-line formula hand-written in this
repo, in the same language, by the same author. That is a second transcription, i.e. tier 1
evidence, and it is exactly the failure mode the module docstring says tier 4 exists to catch
("a routine that is SELF-CONSISTENTLY WRONG"). Contrast CvM's tier 4, which uses
`scipy.stats.cramervonmises` and IS independent. `scipy.stats.anderson` exists and was not used.
| Downgrade to `partial` ("second transcription, not an independent implementation"), or compare
against scipy.
IMPORTANT | analyzers/instrument_validation.py:384-389 | C3 tier 4 = `pass`, "it IS the R
implementation, over a bridge". Circular: the tier asks whether OUR number agrees with an
independent one, and no comparison exists. Everything tier 4 would catch here lives in the
bridge - CSV write, decimal separator, column order, result parsing - and
`docs/iid_checks/C3_serial_copula.md:39-41` says the numeric path is pinned by nothing. Worse,
the comparison IS available and was declined: `reference/r_reference_values.csv` carries
`durations_iid/serial_indep_global_statistic = 0.00739025142047825` with its seed, N and lag.max,
and test_r_cross_implementation.py:228-242 pins only that those fields exist. | `absent` with
"bridge unvalidated; reference value in hand, comparison needs R at test time" - or run the
comparison behind a `rscript_path() is None` skip.
IMPORTANT | analyzers/instrument_validation.py:414-421 vs 210-217 of cvm_cramer_von_mises.py |
CvM's three tier rows are all `pass` with no row recording that the m > 1 path - the one the
gapped records here actually take - has NO asymptotic calibration (it raises) and NO external
validation of any kind (the unweighted sum is an extension beyond the source; this ledger:268-272
says so). The tier-4 detail is scoped to the gamma=1, m=1 path, but the table cell reads `pass`
for the instrument. | Add the m>1 caveat to the tier-4 detail, as C5's tier 4 already does for
its max-over-lags aggregation.
MINOR | tests/test_instrument_validation.py:211 | `assert "2 levels" in render_tier_table_markdown(data) or True` is a tautology - it can never
fail, and the trailing comment admits the table uses a different phrasing. A dead assertion in a
file whose whole subject is not overstating evidence. | Delete it or assert the real string.
MINOR | tests/test_instrument_validation.py | No test pins the C1 tier-3 verdict, which is the
one this pass got wrong (see the C1 entry above). `test_c3_is_not_credited_...` and
`test_xi_tier_3_is_partial_...` show the pattern; C1 is the missing case. | Add one.
MINOR | bench/results/instrument_report.md, generated | C5 and C6 tier 2 and dcor tiers 2 and 3
render as a bare `absent` with no detail while every other cell reads "verdict - detail", because
`render_tier_table_markdown:504` falls back to the bare constant when no `TierRow` exists. The
absent-with-a-reason rows and the absent-by-omission rows are indistinguishable. | Emit "absent -
not assessed" for the fallback.
MINOR | analyzers/instrument_validation.py:481-486 | The header table documents four tiers; the
matrix renders only 2, 3 and 4. A reader sees tier 1 defined and then never scored. | Say tier 1
is code review and is not tabulated.
GREEN | analyzers/instrument_validation.py:204-210 | Real positive control: the builder RAISES if
a row marked `agrees=True` differs by more than 1e-3, and tests/test_instrument_validation.py:92
proves it fires on a corrupted published value. Verified in the generated report - the three
divisor rows carry `NO` and the reconciling factor 1.0142. | none
GREEN | jobs/active/instrument_validation.py | Contract-clean: every input is a `Dataset` loaded
with `job.load_df` (so each CSV's hash reaches the run identity), the builder is a pure step with
`divergence_threshold` and `n_perm_in_tie_study` declared as step kwargs, no pipeline module
imports `bench/` or `tests/`, and `bench/instrument_report.py` imports the pipeline (legal
direction). `_divergence_levels` picks the FINEST crossing and its ordering is right (0 -> 1e9,
descending), pinned by test:214. | none

### Part 2 - check_ledger_q1_070423.py, reference/, CLAUDE.md

GREEN | jobs/active/check_ledger_q1_070423.py vs jobs/active/t2star_q1_070423.py | The carve IS
identical, verified on the artifacts rather than by reading: same dataset path, same
`_filter_step(RAMSEY_CONFIG)` -> `_final_stage` -> `t2star.run` chain, byte-identical
`_T2STAR_THRESHOLDS`, byte-identical `_windows_run` body, same step kwargs (gap_mult=10.0, k=1.0,
use_uncertainty=True). Window counts per threshold in the shipped ledger match the shipped panel
windows table EXACTLY: 1/4/356/682/130/30/11/2/1 at 1-9 us. The docstring's comparability claim
holds. | none
GREEN | ledger:992-1000 real-data claims, re-derived from the artifact | 104 rows; verdicts
31 fail / 29 underpowered / 26 not computed / 12 pass / 6 not interpretable - exact. C3 in-spec
p-values 0.004496, 0.0005, 0.0005, 0.0005, 0.016484 and calendar 0.004496, 0.3711, 0.6479,
0.4600, 0.4800 - exact. No timeout row (the two `not computed` C3 rows read "n=3 too short for
lag.max=5"). In-spec 3-6 us: every check that returned a p rejected; calendar 4-5 us: 7 pass /1
underpowered and 5 pass /3 underpowered. All accurate. | none
MINOR | jobs/active/check_ledger_q1_070423.py:58 | `extra={"run_name": PREFIX}` makes the ledger
artifact's `dataset_id` "q1_27h_0704_ledger" while the panel's is "q1_27h_0704_dataset" - the same
physical record under two ids, in the two artifacts the docstring says are to be read side by
side. A join on `dataset_id` between them finds nothing. | Use the panel's run_name, or add a
column that carries it.
MINOR | CLAUDE.md:80 | "C3 copula-via-R (needs Rscript, ABSENT HERE)" - Rscript 4.5.3 is
installed and C3 has run on the real record. This is a SEVENTH site of the stale
"R is absent" claim (the fix layer corrected five, and docs/iid_checks/C3_serial_copula.md:63
still carries an eighth in "which is every row"). The agent guide is the one a fresh session
reads first. | Update the line.
MINOR | CLAUDE.md Layout | Stale in three more ways after this pass: `checks/` does not list
`cvm_cramer_von_mises.py`; `analyzers/` does not list `permutation.py`, `check_ledger.py` or
`instrument_validation.py`; `jobs/active/` does not list `instrument_validation.py` or
`check_ledger_q1_070423.py`. Also the doc-hygiene section still says the tracked .md files are
"README.md, docs/**, and bench/results/promotion_report.md", but `bench/results/instrument_report.md`
is now a tracked generated .md too. | Refresh the four lines.
MINOR | reference/load_haul_dump.csv vs tests/fixtures/load_haul_dump.py | The ledger says the
record is "read by BOTH the suite and the figure so the two cannot drift". It is actually
DUPLICATED - the CSV for the figure, a hardcoded numpy array for `test_checks_published_values.py`
- and kept in step by `test_instrument_validation.py:85-89`, which does compare them. So the
protection is real but it is a test, not a single source; the claim as written is stronger than
the arrangement. | Restate, or have the fixture module read the CSV.

### Part 3 - the sweep for a fourth silently-failed edit

FOUND, ONE | tests/test_tier3_calibration.py | See the IMPORTANT entry above: three sentences
say the instrument report/figure imports `measure_asymptotic_size`; nothing does. This is the
fourth instance of the class the brief asked me to look for - a described edit that never landed
(or a claim written for an edit that was then not made), and unlike the first three it is not
visible in the printed output because nothing prints it.
SWEPT, CLEAN | Everything else in the Increment C sections was checked against the code rather
than assumed: the signed-mean DELIVERABLE block in `bench/xi_ties.py:253-271` HAS landed; the Q2
x-quantisation HAS landed (`x_levels=y_levels` at :209, and `tie_fraction_x` is 0.000 in all 24
rows while `tie_fraction_x_q2` is 0.80-1.00); `_divergence_levels` really does take the finest
crossing; `XI_TIE_CUTOFF` is still 0.0 as the results section says it deliberately left it;
`bonferroni_z_crit`/`null_se` are genuinely imported from `calibration_summary` rather than
re-defined; the stated resolutions (+/-0.0201 at alpha=0.10, +/-0.0067 at 0.01) reproduce; the
style ratchet is 17 with zero hex/fontsize literals in the new plot module; `tests/` = 298
passed, 2 skipped, ruff at the one known F841. No further doubled sentences or unreflowed
insertions in the new files (only analyzers/checks/c3_serial_copula.py:13, logged above).

MINOR | analyzers/instrument_validation.py | 564 lines carrying four dataclasses, three builders,
two divergence readers and a markdown renderer. `render_tier_table_markdown` is 100 lines of
string assembly living in `analyzers/`. It is pure (the file writing is in bench/), so no rule is
broken, but the file is at the unwieldy line on arrival. | The markdown renderer is the natural
split.
MINOR | plots/instrument_validation_plot.py:227 | Caption says "{n} of {n} quantities reproduce
exactly" while `agrees` is a 1e-3 RELATIVE tolerance and mu_hat is 54.7222 against a printed
54.72. "Exactly" is the wrong word on a comparison against rounded published values. | "to the
precision the paper prints".
MINOR | analyzers/shape_stats.py:35 | A stray blank line splits item 5 from item 6 of the
"SIX DIFFERENCES" list - left by the median-window prose correction. Cosmetic, but it is the same
edit-hygiene signature. | Remove.

### Part 1 answer - does the corrected prose now UNDERSTATE?

Mostly no. The six corrections were checked against the artifacts they cite and each is either
exact (shape_stats' 11/14 median windows; the 1/N divisor consequence - sqrt(5/4) = 1.118, so
"about 12% on a five-event segment" is right) or correctly refuses a claim it cannot support
(the C3 cost, where no exponent is now asserted). Two places have drifted into understatement
since, both because the world moved after the prose was written:

MINOR | analyzers/checks/c3_serial_copula.py:206 | "A separate run had n = 682 unfinished at
580 s" is now the WEAKEST thing known about the cost: the shipped ledger computed C3 at
n_events = 682 on the in-spec clock inside the 900 s timeout (statistic 1.511, p 0.0005), and no
row timed out. A reader sizing `timeout_s` is being shown only the discouraging measurement. |
Add the completed n=682 point.
MINOR | docs/iid_checks/C3_serial_copula.md:44-45 and bench/results/promotion_report.md:5 | Both
still describe C3 as exercised on iid input only. It has since run on the real 912-day record at
ten thresholds on two clocks and rejected on five of them. Uncalibrated is still right;
"smoke-tested on iid input only" is no longer the whole of what has been run. | One sentence.

### VERDICT (reviewer, Increment C + the unreviewed B fix layer)

DO NOT SHIP as it stands. Nothing here produces a WRONG NUMBER - every deliverable figure in the
ledger's Increment C sections was re-derived from the CSVs and artifacts and every one is
correct. The defects are all of one kind, and it is the kind this pass exists to prevent: a
verdict or a docstring that claims more evidence than exists.

Blocking (four IMPORTANT, all one-line prose or one-word verdict changes):
1. tests/test_tier3_calibration.py - the instrument report does NOT import
   `measure_asymptotic_size`; three sentences say it does; the tier-3 "size measured" verdicts
   have no number behind them anywhere.
2. tests/test_tier3_calibration.py:13,189 - 0.069 and 0.0757 are pooled means over two Weibull
   shapes (0.064/0.074 and 0.065/0.0865), labelled as cells. CLAUDE.md's first claims-discipline
   rule, fourth instance.
3. tests/test_tier3_calibration.py:14 - the asserted "direction" is not asserted and runs the
   other way on this file's own generator (0.0342 / 0.0292 / 0.0383 at n=20, conservative).
   analyzers/instrument_validation.py:353-358 then prints "anti-conservative" as a `pass`.
4. analyzers/instrument_validation.py - C2 tier 4 `pass` rests on a hand-written formula in the
   test file, not an independent implementation; C3 tier 4 `pass` is circular and declines a
   comparison whose reference value is already in `reference/r_reference_values.csv`.
Plus the Tier 4 membership band, which accepts the pre-Increment-B estimator and five other
wrong ones - keep it, but stop calling it the evidence for eq (8) (the y-only case is).

---

## Increment C - fixes after the COLD-READ audit

FIXED | THE FOURTH FAILED EDIT, exactly where the audit was asked to sweep |
`measure_asymptotic_size` lived in `tests/test_tier3_calibration.py` while three of its
docstrings and the report's three tier-3 verdicts claimed the report imported it. It did not -
`grep -rn` matched that one file - and the pipeline may not import `tests/`, so "asymptotic
size measured" was a verdict resting on no number in the artifact. The function now lives in
`analyzers/instrument_validation.py`, the TEST imports it, and `build_instrument_validation`
calls `measure_all_asymptotic_sizes()` so each tier-3 row quotes a number that run produced.
Pinned by `test_the_tier_3_rows_quote_a_number_this_run_produced`.

FIXED | pooled means presented as cells | "the bench measures 0.069 (C1) and 0.0757 (C2) at
n=20" are POOLED MEANS over the two Weibull shapes. The table's actual cells, verified:
C1 0.0640 (shape 0.75) and 0.0740 (1.50); C2 0.0650 and 0.0865. This is verbatim the first
instance CLAUDE.md's claims-discipline section records - `report.size_by_n` - reintroduced in
new prose. All four cells are now named with their shapes and mc_se.

FIXED | the direction claim was INVERTED, and the reason is interesting | The prose said
"anti-conservative at small n as the source reports". Measured with the file's own function at
n=20: C1 0.0342, C2 0.0292, CvM 0.0383 - CONSERVATIVE. Not a contradiction of the bench: the
bench generates WEIBULL gaps (shapes 0.75, 1.50) and this measurement generates EXPONENTIAL
ones. Both are legitimate iid nulls, since these check trend against renewal rather than
exponentiality, so the finding is that the size of an asymptotic check here DEPENDS ON THE GAP
DISTRIBUTION through the estimated gamma_hat - and therefore no single direction can be
asserted. Both numbers now appear side by side and no direction is claimed anywhere.

FIXED | two tier-4 `pass` cells with nothing independent behind them |
  - C2: its "reference" is `_classical_ad`, five hand-written lines in our own test file -
    same language, same author. That is tier-1 evidence wearing a tier-4 label. -> partial.
  - C3: "it IS the R implementation, over a bridge" compares R to itself. -> absent, with the
    comparison that WOULD settle it named (our bridge against the fixture's
    `serial_indep_global_statistic` at seed 707), and why it is not done: it needs R at test
    time, which the suite refuses to require.
Both pinned by `test_no_tier_verdict_claims_evidence_that_does_not_exist`.

FIXED, and NOT found by the audit | tests/test_tier3_calibration.py | The exactness test seeded
itself with `hash((check_name, layout_name))`. Python randomises string hashing per process, so
the seed differed on every run - measured 90936 / 90550 / 90044 across three interpreters - and
a test whose entire claim is that a rate lands inside a band was silently FLAKY. It failed once
during this fix round, which is how it surfaced. Replaced with an explicit `SEEDS` table;
verified 12 passed on three consecutive runs. Same defect class `block_permutations` raises
over when a caller omits an rng.

NOTE ON THE TIE-BREAK BAND | The audit is right that `[q001, q999]` is wide: it accepts the
pre-Increment-B tie-free reduction and four other wrong estimators, and `xicor_q001` equals
`xicor_min`, so the lower edge is the whole observed support. The eq (8) claim is nonetheless
discharged - by `xi_tied_y_only`, which is deterministic, agrees to 2.3e-16, and catches every
one of those mutants. The membership test is a weak SUPPLEMENT, not the evidence. Docstring
should say so; recorded rather than silently relied on.

STATE | 300 passed, 2 skipped. ruff at the one known pre-existing F841. Ratchet 17.

STRENGTHENED | tests/test_r_cross_implementation.py, xi tier 4 | Prompted by the question
"if it was published why not just use it": `scipy.stats.chatterjeexi` EXISTS as of scipy 1.17
and is installed here. Measured against ours on all four fixture cases: agreement to 1e-10,
including `xi_tied`, where the R comparison cannot pin equality because XICOR breaks x-ties at
random and scipy (like us) breaks them deterministically. So the tied path - the only behaviour
Increment B changed - now has an EXACT external check for the first time; the membership band
was the weak supplement the audit correctly called out.
This is stronger tier-4 evidence than the R fixture for a reason unrelated to authority: it
runs unconditionally in the suite, whereas the R comparison reads a frozen file that only
regenerates on a machine with R. 5 tests added, 305 passing.
OPEN, a genuine scope question for the human: `chatterjee_xi` could now simply CALL scipy's,
with zero numerical effect (proven identical on four cases). Kept as ours for now - it is
validated, it is ten lines, and swapping an implementation during a review pass is not a
review-pass change. Recorded so the option is not forgotten.

---

## Increment C addendum - the cross-dataset independence survey

NEW | analyzers/independence_survey.py, plots/independence_survey_plot.py,
jobs/composite/independence_survey.py, tests/test_independence_survey.py (20 tests) |
Seven figures - one per `battery.ROW_KEYS` entry - each a grid of all 34 6D2S datasets by the
10-rung T2* ladder, both clocks stacked. Answers a question the per-dataset ledger cannot:
whether iid fails on THIS record or on the device.
Coloured by VERDICT, not p-value. A p-value heatmap is the obvious design and is wrong: on a
short window "did not reject" carries no information, and a pale cell would read as evidence
of independence when it is evidence of nothing. The ledger's three-way verdict already
encodes that against the bench, so the survey reuses it and prints p inside the cell.

REMOVED | jobs/active/check_ledger_q1_070423.py | Redundant, and redundant in the way its own
neighbour warns against. `jobs/composite/check_ledger_q1.py` already ledgers both q1 datasets
and gets its windows by `include` + `.ref` from the T2* jobs, so there is ONE carve wiring;
the job I added in this pass re-derived the carve, which that file's docstring explicitly
calls out as "fifty lines that must stay byte-identical to those jobs forever". Deleted. Its
output directory stays - output/ is append-only - so the figure already produced is intact.

NEW | analyzers/check_ledger.py `include_c3` | The survey needs C3 off: it is out of process,
costs 130 s at n=355 and grows steeply, and has no bench cell so it cannot be scored on the
grid anyway. Rows are OMITTED rather than written `not computed`, because a blank row would
claim the check was attempted. Threaded through `make_inputs_from_windows` so it is reachable
from a job and lands on the provenance label - the defect class the reviewer caught with
`xi_seed` and which I nearly repeated: the first version added the field to the dataclass only.
Also caught before it shipped: the first guard was `if not inputs.include_c3: return results`,
which returns the WRONG TYPE - the function goes on to build row dicts from `results` and
returns those. Fixed to a positive guard around the block.

FOUND, and PRE-EXISTING | the threshold ladder was not bit-identical between the survey and
the panels | The composites build the ladder as `k * 1e-6`; the T2* jobs write the literals
`1e-6 ... 10e-6`. These differ: 1e-6 is not exactly representable, so `5 * 1e-6` is
4.9999999999999996e-06 and `10 * 1e-6` is 9.999999999999999e-06 - two of the ten rungs land on
a different double from the panel's. `k / 1e6` reproduces all ten literals exactly.
`jobs/composite/check_ledger_q1.py` has had this since it was written, so its ledger has been
scoring a marginally different ladder than the panel it describes.
IMPACT IS NIL IN PRACTICE and this says so rather than dressing it up: the gap is about 4e-22
and no measured T2* value will fall inside it, so no window has ever changed. It is fixed
because the survey CLAIMS to score the same ladder, and a claim true only to 15 significant
figures is not the claim being made. Both composites now use `k / 1e6`, and
`test_the_survey_scores_the_same_threshold_ladder` asserts exact float equality against the
T2* job's literals - which is how it was found.

THE RESULT | 34 datasets x 10 thresholds x 2 clocks, C3 excluded. Rejection share of DECIDED
cells (pass + fail as the denominator, not all cells - an underpowered grid must not look
reassuring):
    instrument                        in_spec   calendar
    C1 LR asymptotic                  8/8 = 1.000    31/67 = 0.463
    C1 LR permutation                 8/22 = 0.364   25/118 = 0.212
    C2 AD asymptotic                  9/18 = 0.500   33/51 = 0.647
    C2 AD permutation                 9/22 = 0.409   32/118 = 0.271
    C5 rank studentized               37/168 = 0.220 16/167 = 0.096
    C5 rank raw                       36/168 = 0.214 15/167 = 0.090
    C6 exchangeability                49/168 = 0.292 18/167 = 0.108
Two readings worth carrying, both provisional: the IN-SPEC clock rejects at roughly twice the
calendar rate on every rank check, which is the same clock disagreement the single-dataset
ledger showed, now across 34 records; and C1/C2 have almost no usable cells on the in-spec
clock at all (179 of 208 `not computed`), because `tau == T_N` there makes eq (7) singular -
exactly what `include_tau_checks` documents. The asymptotic C1 in-spec cell reads 8/8 = 1.000
and MUST NOT be read as "C1 always rejects": it is 8 decided cells out of 208.

---

## Increment C addendum - REVIEWER PASS (independent, fresh context, 2026-08-14)

SCOPE: the cross-dataset independence survey + the unreviewed Increment C audit fixes.
COMMIT: 2f1bd95 (working tree; nearly all files in scope are untracked)

### Manifest
- [ ] analyzers/independence_survey.py (278)
- [ ] plots/independence_survey_plot.py (289)
- [ ] jobs/composite/independence_survey.py (261)
- [ ] tests/test_independence_survey.py (254)
- [ ] analyzers/check_ledger.py include_c3 (503)
- [ ] analyzers/instrument_validation.py tier verdicts (635)
- [ ] tests/test_tier3_calibration.py SEEDS (226)
- [ ] tests/test_r_cross_implementation.py scipy (288)
- [ ] ladder k/1e6 sweep
- [ ] claim sweep for a fifth failed edit

### Reviewed
- [x] src/quebra/core/job.py (_load_dataset, _load_dataframe_raw)
- [x] src/quebra/core/dataset.py
- [x] src/quebra/schemas/ramsey_series.py
- [x] src/quebra/schemas/null.py
- [x] tests/test_load_dataset_contract.py
- [x] docs/WRITING_A_SCHEMA.md
- [x] G1 .gitignore / Makefile / pytest.ini  (R0.1, R0.2 acceptance: all 5 criteria PASS)
- [x] G2 jobs/{reference,rscripts,bench} move
- [x] G3 spec/ledger/ (R0.3)
- [x] G6 .claude/ tracking + repo weight
- [x] G4 rename (R0.4)
- [x] G7 tests/ correctness spot-checks
- [x] G5 project identity + em dash sweep (R0.5)

### Findings (this pass)

IMPORTANT | analyzers/independence_survey.py:117-131 (`counts`) and the ledger's THE RESULT
table | The census silently drops every ledger row the survey cannot key. Measured on the
shipped artifact
(output/independence_survey_340722_20260814_030436/independence_survey.pkl): every grid is
34x10 = 340 cells but `sum(counts.values())` is 208 (in_spec) / 207 (calendar) - 132/133
cells are NaN and land in NO verdict bucket. Cause: `check_ledger._blank_row` writes
`check_id="(all)"`, which matches no `ROW_KEYS` entry, so a rung the carve declined
entirely vanishes instead of counting as `not computed`. It is not scattered: the 1 µs and
10 µs columns are 34/34 NaN in EVERY grid, i.e. two of the ten ladder rungs are wholly
absent from all seven figures and from every denominator. The ledger's prose "(179 of 208
`not computed`)" and "8 decided cells out of 208" therefore quote a denominator over the
subset that produced keyed rows; the honest statement is 311 of 340. `rejection_share_of_
decided` itself is unaffected (its denominator is pass+fail), but the "out of 208" framing
is not. Minimal fix: map `(all)` blank rows onto every ROW_KEYS entry for that
(dataset, threshold, clock), or record the count in the (currently always-empty) `dropped`
field and quote 340 as the grid size.

IMPORTANT | analyzers/independence_survey.py:225-226 | A second, differently-caused silent
blank: `run_battery` drops the C2 asymptotic row when m > 1 segments (Kvaloy-Lindqvist
4.2). Measured: the C2-asymptotic/calendar grid has 35 cells that are NaN while the C2-
permutation/calendar grid has a verdict at the same (dataset, threshold) - 168 NaN vs 133.
A reader comparing the two C2 figures sees cells disappear with no marking and no note. Same
class as above: absence is drawn and counted as nothing rather than as a stated reason.

MINOR | analyzers/independence_survey.py:142,244 | `dropped: dict[str, int]` is constructed
as `{}` unconditionally and never written - the one field that would have carried the two
findings above is dead. Either populate it or remove it.

CRITICAL | tests/test_independence_survey.py:98-113 | The ladder guard does not guard the
survey. `test_the_survey_scores_the_same_threshold_ladder` RE-TYPES the survey's ladder as a
literal in the test (`survey_ladder = [(f"{k} µs", k / 1e6, True) for k in range(1, 11)]`)
instead of reading `THRESHOLDS` out of `jobs/composite/independence_survey.py`; it cannot
read it, because `_module_constants` uses `ast.literal_eval`, which raises on the ListComp
and silently skips it (`"THRESHOLDS" not in SURVEY`). So the assertion compares the T2* job
to a copy of the reference. MUTATION TESTED in a scratch mirror (no source touched), 20 tests
each:
  - survey ladder `k / 1e6` -> `k * 1e-6` (the exact drift this pass was written to catch):
    20 passed, NOTHING FAILED.
  - survey ladder `range(1, 11)` -> `range(1, 9)` (drops two rungs): 20 passed, NOTHING FAILED.
  - T2* job `5e-6` -> `5.5e-6`: 1 failed (so the guard is real in ONE direction only).
The ledger's claim that this test "asserts exact float equality against the T2* job's
literals - which is how it was found" is therefore half true: it pins the panel side, not the
survey side, and would not have caught the defect had it been introduced in the survey.
This is CLAUDE.md's recorded "positive control must fail when the thing it guards is broken /
assert against the REAL output, not a copy of the reference". Minimal fix: import THRESHOLDS
from the job module (or `ast.literal_eval` the elements after unrolling the comprehension)
and compare THAT to the T2* literals.

VERIFIED | the ladder fix itself | `k / 1e6` reproduces all ten T2* literals to exact float
equality (checked element-wise against the AST-read literals); `k * 1e-6` differs at k=5
(4.9999999999999996e-06) and k=10 (9.999999999999999e-06), as claimed. `grep -rn "\* 1e-6"`
over the repo finds only the two explanatory comments - no construction site remains. Both
composites use `k / 1e6`.

VERIFIED by mutation | the other survey guards DO bite: GAP_MULT (survey side and T2* side),
MIN_EVENTS_PASS, INCLUDE_C3 each fail exactly one test when mutated. `_module_constants`
handles AnnAssign, so `DATASET_FILES` is really read.

IMPORTANT | jobs/composite/independence_survey.py:73-80 | False comment, and the one constant
in the block that is NOT pinned by a test. The block is headed "Verdict parameters ... Equal
to jobs/composite/check_ledger_q1.py so the two ledgers are comparable" and `SEED = 20260814`
sits inside it, but check_ledger_q1 uses `SEED = 20260812`. The seed selects the permutation
draws, hence the permutation p-values, hence the verdicts - so the two ledgers are NOT
comparable cell-for-cell on the one dataset they share (q1 070423), which is what the comment
asserts. The parametrized equality test lists ALPHA/MIN_EVENTS_PASS/TIE_CUTOFF_DISTINCT/
LAG_MAX/N_PERMUTATIONS and omits SEED, so nothing catches it. Minimal fix: either use the
same seed and add SEED to the parametrize list, or say plainly that the seed differs and the
permutation p-values are therefore not cell-comparable.

VERIFIED | plots/independence_survey_plot.py | The overview's labelling is lookup-based, not
positional: `grid = pd_.grid(key, clock)` and the row label comes from `grid.label`, the
column title from the same `clock` used in the lookup, so the off-by-one class asked about
cannot occur. `_draw_transposed` uses `grid.verdicts.to_numpy().T` with y ticks from
`.columns` (thresholds) and x ticks from `.index` (datasets) - genuinely the same data,
correctly relabelled. Dropping p-values in the overview is stated in the docstring AND in the
caption ("Read a p-value off the per-instrument figures, not this one"). `theme.verdict_color`
raises on an unknown verdict, so no silent colour fallback.

IMPORTANT | plots/independence_survey_plot.py:112-118,272-276 vs docs/FIGURE_STANDARD.md:84-91
| "Every panel that drops data says how much it dropped, in the panel." These panels drop 132
of 340 cells per grid (the unkeyable `(all)` rows - the whole 1 µs and 10 µs columns) and, for
C2 asymptotic, 35 more for m > 1 segments, and paint all of them with the `not computed`
colour without saying anywhere in the panel that they were dropped rather than computed-and-
declined. The caption's counts ("N of M decided cells reject") do not include them either. An
annotation of the form `no rows: 132 of 340 cells` is what the standard asks for.

MINOR | plots/independence_survey_plot.py:93-94, 251-252 | Caption points at the wrong colour.
Both captions say "Grey cells carry no information - a non-rejection on too few events is not
evidence of independence", but the class that sentence describes is `underpowered`, which
`theme.VERDICT_COLORS` renders #B9A87A (tan). The greys are #9A9A9A (`not interpretable
(ties)`) and #D8D8D8 (`not computed`). A reader following the caption will read the wrong
cells. Name the verdict, not the colour.

MINOR | plots/independence_survey_plot.py:140-143, 284-285 vs FIGURE_STANDARD.md:59-62 |
"Units go in the axis label, in parentheses. Never in the tick labels." The threshold ticks
are "1 µs" ... "10 µs" under an axis labelled "Threshold" (and in the overview, under no axis
label at all). Should be `Threshold (µs)` with bare numeric ticks.

MINOR | plots/independence_survey_plot.py:26,133 | Two small couplings: the plot imports the
private `_fmt_p` from the analyzer, and the text-colour test compares against the string
literals `("pass", "fail")` rather than `VERDICT_PASS`/`VERDICT_FAIL`, so renaming a verdict
constant silently changes label contrast instead of failing.

IMPORTANT | bench/results/instrument_report.md:35 | THE FIFTH STALE ARTEFACT, and it is the
one the pass was written to fix. The shipped report's xi tier-4 cell still reads
"pass - matches XICOR exactly tie-free; inside its tie-break distribution when tied" - the
weak tie-break-membership claim the audit called out. `analyzers/instrument_validation.py:501`
now says "matches `scipy.stats.chatterjeexi` to 1e-10 ... INCLUDING the tied one", but
`bench/instrument_report.py` was never re-run, so the generated file on disk does not carry
the fix. Verified by regenerating into a scratch copy and diffing: exactly that one line
changes; the file was restored byte-identical (mtime included) and NOT left modified. Minimal
fix: re-run `python bench/instrument_report.py`.
VERIFIED alongside it: the three tier-3 numbers in the shipped .md (0.0342 / 0.0292 / 0.0383)
are exactly what `measure_all_asymptotic_sizes()` returns on this checkout, so that edit did
land and the numbers are this run's.

IMPORTANT | analyzers/instrument_validation.py:413-416,427-429 | The bench cells are named by
shape only, and shape is not the factor that moves them. "the bench measures 0.0640 and
0.0740 on Weibull shapes 0.75 and 1.50" holds `arm`, `quantised`, `censoring_target` and
`clock` silently fixed at (A_iid_weibull, False, 0.00, in_spec). Measured: at n=20, C1
asymptotic, shape 0.75 there are TEN table rows spanning 0.0630 to 0.1762 - the
censoring_target=0.25 sibling is 0.1415, more than double the quoted cell. Same for C2
(0.0000 to 0.4000 at shape 0.75). The quoted values are correct as cells; the sentence
identifies them by the one factor that varies least. This is the recorded recurrence
("never call a mean the cell" / name the factor): fixed once for the shape pooling, now
under-specified on the remaining three factors. Minimal fix: say "primary null cell
(A_iid_weibull, uncensored, unquantised, in-spec)".

IMPORTANT | analyzers/instrument_validation.py:117-119,397 | A seed that decides a REPORTED
number is a module default, not a declared parameter. `ASYMPTOTIC_SIZE_SEED = 777` and
`ASYMPTOTIC_SIZE_REPLICATES = 1200` are defaults; `build_instrument_validation` calls
`measure_all_asymptotic_sizes()` with no arguments, and neither the seed nor the replicate
count reaches `meta` or the rendered table. The tier-3 verdicts (0.0342 / 0.0292 / 0.0383)
are functions of that seed. This is CLAUDE.md's rule verbatim and the same defect class the
previous reviewer caught as `xi_seed`. Minimal fix: make them arguments of
`build_instrument_validation` and put them in `meta`.

MINOR | analyzers/instrument_validation.py:413 etc. | The tier-3 rows quote four decimals with
no Monte-Carlo error. At 1200 replicates the MC SE is ~0.005, so "0.0342" against nominal
0.05 is ~3 SE - real, but the row does not let a reader see that. Carry the SE.

MINOR | analyzers/instrument_validation.py:151 | `rejected += int(p is not None and p <= alpha)`
folds a `p_value=None` into the denominator as a non-rejection, so a check that silently
declined would report a LOW size rather than a problem. Measured 0 of 300 Nones for each of
C1/C2/CvM at n=20, so the defect is latent, not active. Minimal fix: count Nones and raise or
report them.

IMPORTANT | tests/test_instrument_validation.py:246-261 | The pin that is supposed to prove
the tier-3 rows quote a number THIS RUN produced only checks that the sentence contains a
four-decimal number. `re.search(r"size measured at n=\d+: 0\.\d{4}", row.detail)` matches
`"size measured at n=20: 0.9999"` just as happily (demonstrated). If the f-string were
replaced by a hardcoded literal - which is exactly the defect the docstring names - the test
passes. CLAUDE.md: "assert against the REAL output". Minimal fix: call
`measure_all_asymptotic_sizes()` in the test and assert the formatted value appears in the row.

MINOR | tests/test_instrument_validation.py:241-243 | The direction guard is a literal-string
guard: it asserts only that the exact old sentence "anti-conservative at small n as the source
reports" is absent. Any newly written direction claim ("conservative at small n, as the theory
predicts") passes untouched. Pin the property, not the sentence that failed.

VERIFIED | tests/test_tier3_calibration.py:15-28 | Here the four bench cells ARE fully
specified - n, alpha, clock, arm A_iid_weibull, both shapes, and mc_se - and every number
checks out against bench/results/size_table.csv (0.0640/0.005473, 0.0740/0.005853,
0.0650/0.005512, 0.0865/0.006286), as does "anti-conservative in all four". The exponential
figures 0.0342 / 0.0292 / 0.0383 reproduce exactly. So the under-specification finding above
applies to the REPORT's rows, not to this docstring.

VERIFIED | tests/test_tier3_calibration.py:83-88 | The `hash()` flakiness is genuinely gone:
`SEEDS` is deterministic across processes and the file passes 12/12 on three runs with
PYTHONHASHSEED=random. (Nit: the comment says "written out rather than derived" while the
table is in fact derived from `enumerate`, so inserting a layout at the front silently
renumbers every seed. Deterministic, but not what the comment says.)

VERIFIED by mutation | tests/test_r_cross_implementation.py:248-288 | The scipy cross-check
really would fail on a wrong estimator. scipy 1.17.1 is installed, `chatterjeexi` is present,
and neither of the 5 tests is among the suite's 2 skips. Ours vs scipy agrees to 0.0e+00 on
all four cases; two plausible wrong estimators (the pre-P5 tie-free reduction, and
`rankdata(method="min")`) differ from scipy by 1.7e-2 to 4.3e-2 on `xi_tied` and
`xi_tied_y_only`, i.e. 8 orders of magnitude above TOL=1e-10. The tied path does now have an
exact external check, as claimed. (The tie-free cases cannot discriminate - all estimators
agree there - so the discrimination rests entirely on the two tied cases; that is what the
docstring says.)

MINOR | analyzers/shape_stats.py:78-110 | `chatterjee_xi` gives the source (Chatterjee 2021,
JASA 116(536):2009-2022) and writes the formula out, but cites no EQUATION NUMBER, which
CLAUDE.md requires of a transcribed formula. The ledger's own prose refers to "eq (8)".

MINOR | analyzers/independence_survey.py:254-271 | `survey_summary` is never called by the
job - nothing materializes it - so the seven-row table quoted in the ledger is not a run
product; it has to be recomputed by hand from the pickle. It does reproduce exactly (all 14
numbers checked against the shipped artifact), but a `job.materialize` of the summary would
make the quoted table an artifact rather than a scratch computation.

MINOR | analyzers/independence_survey.py:245-250 | `meta` is asserted, not derived:
`"c3_excluded": True` is hardcoded regardless of what the ledgers actually ran, and `alpha`
is a step kwarg that is only recorded - the verdicts came from each ledger's own alpha. Run
the survey over ledgers built with `include_c3=True`, or with a different alpha, and the
artifact still claims otherwise. (`CheckLedger` carries no `include_c3` field to derive it
from, which is the underlying gap.) Every figure caption repeats the C3 sentence
unconditionally for the same reason.

MINOR | analyzers/instrument_validation.py (635 lines) | Growing unwieldy, and mixed-purpose:
a 20-row hand-maintained tier table, a Monte-Carlo size measurement, three artifact builders
and a markdown renderer in one module. The tier table in particular is prose-in-code that a
reviewer must read line by line. Candidate split: the measurement and the tier table.

VERIFIED | provenance and removals | `include_c3=False` and `seed=20260814` both appear on
the step labels in
output/independence_survey_340722_20260814_030436/provenance/*.prov.{json,md}, so the C3 flag
really is threaded to the label and not just onto the dataclass. The two
output/check_ledger_q1_070423_* directories from the deleted job are intact (output/ not
touched). Every source file in scope predates the shipped run (latest source 03:04:29, run
03:06:13), so the shipped figures do reflect this code - with the single exception of
bench/results/instrument_report.md noted above.

### Reviewed (this pass)
- [x] analyzers/independence_survey.py
- [x] plots/independence_survey_plot.py
- [x] jobs/composite/independence_survey.py
- [x] tests/test_independence_survey.py (mutation tested, 7 mutants)
- [x] analyzers/check_ledger.py include_c3
- [x] analyzers/instrument_validation.py tier verdicts
- [x] tests/test_tier3_calibration.py SEEDS
- [x] tests/test_r_cross_implementation.py scipy (mutation tested, 2 mutants)
- [x] ladder k/1e6 sweep
- [x] claim sweep - one stale artefact found (instrument_report.md)

VERDICT: DO NOT SHIP until the ladder guard actually reads the survey's THRESHOLDS and
bench/results/instrument_report.md is regenerated. Baseline reproduced: 324 passed, 2 skipped;
ruff at the one known pre-existing F841.

---

## Survey review - fixes

FIXED, CRITICAL | tests/test_independence_survey.py | The ladder guard compared the T2* job to
a copy of the reference RETYPED INSIDE THE TEST, so it could not see the survey at all.
`_module_constants` uses `literal_eval`, which raises on the ListComp, so `THRESHOLDS` was
never in the dict and the test quietly built its own. Reviewer's mutation table: survey ladder
`k / 1e6` -> `k * 1e-6` and `range(1, 11)` -> `range(1, 9)` BOTH left 20/20 green. It would not
have caught the very defect the ledger says it found. This is CLAUDE.md's positive-control rule
verbatim - "assert against the REAL output, not a copy of the reference".
Fixed with `_evaluated_constant`, which unparses the assignment and evaluates it. Re-mutated
against the real code line (the first attempt patched the COMMENT, which is why one mutant
still looked to survive): `k*1e-6` CAUGHT, 8 rungs CAUGHT, `k/2e6` CAUGHT.

FIXED | analyzers/independence_survey.py | 132 of 340 cells per grid were silently blank.
`check_ledger._blank_row` writes `check_id="(all)"` for a rung that produced no windows or was
declined at segmentation - one row standing for all seven instruments - so matching on
`check_id` dropped them and the cell rendered as "no row". Measured: the 1 us and 10 us columns
were 34/34 absent in EVERY one of the seven figures. Blanket rows are now folded into each
instrument's grid, and any remaining hole is filled with `not computed` (35 of them are
C2-asymptotic, which `run_battery` legitimately drops for m > 1 segments). Every grid now sums
to 340 of 340 with zero blanks.
The consequence for the numbers already reported: the ledger's "179 of 208 not computed" quoted
a denominator over the keyed subset. The honest figure is 311 of 340 for C1 in-spec. The
rejection SHARES are unchanged, because they always divided by pass+fail.

FIXED | bench/results/instrument_report.md | The fifth stale artifact of this pass. The shipped
xi tier-4 cell still read "inside its tie-break distribution when tied" - the weak claim the
previous audit rejected - because the code changed and the generator was never re-run.
Regenerated; it now reads the scipy claim.

FIXED | analyzers/instrument_validation.py | Bench cells named by SHAPE only. Verified: ten
rows share n=20 / C1 asymptotic / shape=0.75 and span 0.0630 to 0.1762, so "the bench measures
0.0640" silently fixed arm, clock, quantised and censoring as well. Both tier-3 rows now name
the full cell (arm=A_iid_weibull, clock=in_spec, quantised=False, censoring=0.00) and state the
range the other rows cover. Same rule, third instance this pass.

FIXED | analyzers/instrument_validation.py, jobs/active/instrument_validation.py |
`ASYMPTOTIC_SIZE_SEED` decided the three reported tier-3 numbers from a bare default, reaching
neither `meta` nor the provenance label. Now a builder parameter, declared in the job, and both
it and `n` are recorded in `meta`. Same class as the `xi_seed` defect.

FIXED | tests/test_instrument_validation.py | The pin only regex-matched four decimals, so
"size measured at n=20: 0.9999" would have passed. It now re-measures and compares.

FIXED | jobs/composite/independence_survey.py | `SEED` sat under a comment claiming equality
with `check_ledger_q1.py` while differing from it. The difference is CORRECT and now says so:
a shared permutation seed would make two jobs' Monte Carlo error identical rather than
independent, so a rung near alpha would land the same way in both and read as corroboration.
New `test_the_two_ledgers_use_different_seeds` pins the difference; the five genuine
comparability parameters keep their equality test.

FIXED | plots/independence_survey_plot.py | FIGURE_STANDARD requires a panel that drops data to
say how much IN the panel. Each per-instrument caption now carries `N of 340 not computed`
alongside the decided count.

NEW | plots/independence_survey_plot.IndependenceSurveyOverviewPlot | One image, all seven
instruments: 7 rows of instruments x 2 columns of clocks, transposed so thresholds are on y and
datasets on x. p-values are dropped at this density and the caption says to read them off the
per-instrument figures. Reviewer verified the labelling is lookup-based (no off-by-one) and the
transpose is the same data.

STATE | 325 passed, 2 skipped. ruff at the one known pre-existing F841. Ratchet 17.

---

# Phase 0 review (2026-08-23)

SCOPE: complete Phase 0 change set (staged + untracked), whole repo
COMMIT: aff894e (HEAD, nothing committed)
AUTHORITY: spec/spec01hygiene.md, AGENTS.md

## Manifest

## Reviewed
- [x] G1 .gitignore / Makefile / pytest.ini  (R0.1, R0.2 acceptance: all 5 criteria PASS)
- [x] G2 jobs/{reference,rscripts,bench} move
- [x] G3 spec/ledger/ (R0.3)
- [x] G6 .claude/ tracking + repo weight
- [x] G4 rename (R0.4)
- [x] G7 tests/ correctness spot-checks
- [x] G5 project identity + em dash sweep (R0.5)

## Findings

### G1 - ignore rules, Makefile, pytest.ini

VERIFIED PASS | .gitignore | all five R0.1 acceptance checks pass: `git check-ignore -v` returns nothing for `tests/`, `README.md`, `.github/workflows/ci.yml`; `git ls-files | grep -c '\.pyc$'` = 0; R0.2 pattern count = 0 in both `git ls-files` and the staged set; no deletions in `git status --porcelain`. | none

VERIFIED PASS | pytest.ini:14 | `pythonpath = .` is confirmed to be the fix. Re-running the same suite with that one line stripped gives `Interrupted: 24 errors during collection`; with it, 330 passed / 2 skipped. Claim in the file's own comment is accurate. | none

IMPORTANT | Makefile:12 (`check: lint types arch test`) | `make check` exits 2 today. `ruff format --check .` reports 21 files, 10 of which are `tests/*` that THIS commit is what makes tracked - they have never been formatted because `tests/` was gitignored. Committing a phase named "Hygiene" with its own checkpoint gate red, and newly tracking 10 unformatted files, is the defect the phase exists to remove. (The other 11 are pre-existing: verified DIRTY at HEAD for the 4 sampled, so the change set introduces no NEW non-test violations.) | `ruff format tests/` before commit; the 11 pre-existing paths can stay for a follow-up.

MINOR | Makefile:16-27 | `types` and `arch` `exit 0` when mypy / lint-imports are absent, and both ARE absent here, so `make check` reports two gates green that never ran. It echoes SKIPPED, so it is not silent, but the exit code a CI job reads cannot tell "passed" from "not installed". Once SPEC 0002 adds CI this becomes a false-green. | Have the skip branch exit a distinct non-zero, or gate on a `SKIP_OPTIONAL=1` env var so CI cannot skip by accident.

MINOR | Makefile:8 | `.PHONY` lists `check lint types arch deps test test-all test-r docs clean` but omits `test-real`, which is a real target at line 45. | Add `test-real` to `.PHONY`.

MINOR | pytest.ini:14 | `pythonpath = .` puts the repo root on `sys.path` unconditionally. After SPEC 0002 adds `pyproject.toml` and `pip install -e .`, the suite will still import the source tree rather than the installed distribution, so the install itself is never exercised by any test. Documented as a stopgap but nothing records that it must be removed. | Add a one-line `# remove in SPEC 0002` to the comment, or note it as a SPEC 0002 requirement.

MINOR | .gitignore:36-41 | R0.1.3 says use Appendix A verbatim and, on a conflict with something already tracked, "stop and report the conflict rather than resolving it". `output/`, `output_backup/`, `output_backup2/` were added (correctly - Appendix A only has `outputs/`, and using it verbatim would have un-ignored 361 MB) but that is a resolution, not a report. Outcome is right, procedure was not followed. | Mention the deviation in the checkpoint banner so it is on record.

### G2 - the bench/reference/rscripts move under jobs/

VERIFIED PASS | jobs/bench/* | all six modules import cleanly as `jobs.bench.*`; `parents[2]` in `instrument_report.py:26,33` and `xi_ties.py:46` is correct for `jobs/bench/<name>.py`; `RESULTS_DIR = Path(__file__).parent / "results"` in `runner.py:59`, `report.py:43`, `xi_ties.py:55` is move-invariant; all five `Dataset(path=...)` sites updated (`jobs/active/check_calibration.py:43,44`, `jobs/composite/check_ledger_q1.py:52`, `jobs/composite/independence_survey.py:146`, `jobs/active/instrument_validation.py:48,50,52,53,54`). `main.py:119,126` globs only `jobs/active` and `jobs/archived`, so `run --all` does not sweep `jobs/bench/`. | none

IMPORTANT | jobs/bench/probe_unresolved.py:44 | `OUT_CSV = repo_root() / "bench" / "results" / "probe_unresolved_out.csv"` was MISSED by the move sweep. It is the one hardcoded absolute-from-root path in `jobs/bench/`; every sibling uses `Path(__file__).parent / "results"`. Running the script writes to a `bench/` that no longer exists, so it either raises or resurrects a stray top-level `bench/results/` next to the real one. | `Path(__file__).resolve().parent / "results" / "probe_unresolved_out.csv"`, matching runner.py:59.

MINOR | jobs/bench/probe_unresolved.py:20-21 | Docstring `Usage: PYTHONPATH=. python bench/probe_unresolved.py` / `Writes: bench/results/...` - both stale after the move, and `PYTHONPATH=.` is now unnecessary given pytest.ini/`jobs` being a package. | Update to `jobs/bench/`.

MINOR | analyzers/windows.py:22 | Points a reader at `scripts/probe_unresolved.py`. That file is now `jobs/bench/probe_unresolved.py`, and `scripts/` no longer exists on disk nor in `.gitignore`. A dangling pointer in a load-bearing analyzer docstring. | Repoint to `jobs/bench/probe_unresolved.py`.

MINOR | jobs/active/instrument_validation.py:17 | `bench/results/promotion_report.md` missing the `jobs/` prefix; the same file gets it right at lines 44-54. | Add `jobs/`.

MINOR | jobs/bench/instrument_report.py:3-4 | Claims `jobs/bench/results/*.md` is tracked because of "the `!jobs/bench/results/*.md` rule in .gitignore". The new .gitignore has exactly three `!` rules (lines 28, 32, 33) and none of them is that one - `*.md` is no longer ignored at all. The stated reason for the file's location is now false. | Reword: markdown is tracked by default now.

MINOR | jobs/bench/results/promotion_report.md:2 | The generated report says "Generated by `bench/report.py`", from `jobs/bench/report.py:529`. Stale path baked into a TRACKED artifact, so it ships wrong until the bench is re-run. | Fix report.py:529 and regenerate, or edit the one line.

MINOR (scope) | jobs/bench/probe_unresolved.py + jobs/bench/results/probe_unresolved_out.csv | Neither existed at HEAD (`git ls-tree -r HEAD` has 14 `bench/` files, not these two). 225 lines of new study script plus a 1021-row generated CSV enter git under a commit whose spec is R0.1-R0.5 hygiene. Its own docstring says "Wired into nothing, imported by nothing." | Either commit it separately with an honest message, or leave it untracked until it has a consumer.

IMPORTANT | CLAUDE.md:98-104, 113-114, 133-135 | CLAUDE.md is now TRACKED for the first time (the old `*.md` rule hid it) and its Layout section is stale on five counts, against its own rule "Every claim in a doc must match the code": (a) `reference/` and `rscripts/` are documented at top level, they are under `jobs/`; (b) ":100 It exists because tests/ is gitignored" - tests/ is no longer gitignored, which is this very phase; (c) ":113 UNTRACKED (gitignored) tests/" - false; (d) ":114 scripts/, monoliths/" - absent from disk AND from the new .gitignore; (e) ":133-135 "The tracked .md files are README.md, docs/**, and jobs/bench/results/promotion_report.md ... because `*.md` is ignored by default" - `*.md` is not ignored any more, and AGENTS.md, CONTRIBUTING.md, CODE_OF_CONDUCT.md, spec/** and every tests/*.md are now tracked too. | Rewrite the Layout and Docs sections, or the em-dash sweep will have shipped a doc that contradicts the commit it rides in.

MINOR | CLAUDE.md:6-11 | The Commands block still carries three `<FILL: confirm...>` placeholders and does not mention `make check` / `make test`, which AGENTS.md 6 makes the canonical entry points. Now tracked, so the placeholder ships. | Fill them from the Makefile.

MINOR | CLAUDE.md vs AGENTS.md | Two agent guides are now tracked side by side and disagree: project name "QRE tool" vs "QUEBRA"; layout `core/ analyzers/ panels/` vs a table whose package row is `src/quebra/`; `output/` vs `outputs/`. AGENTS.md 2 describes a tree that does not exist yet. A new contributor cannot tell which is normative. | State in one line at the top of each which is authoritative and for what.

### G3 - spec/ledger (R0.3)

VERIFIED PASS | spec/ledger/INDEX.md | R0.3 acceptance holds: 1 `.jsonl`, 1 index row. All 503 records parse as JSON after redaction (checked, not assumed - the README's own claim). All six redaction counts in the README table are exact: `<repo>` 671, `<workspace>` 7, `<home>` 20, `<dataset>.pickle` 70, `<dataset>_freq_log` 16, and `<dataset>` totals 646 = 70+560+16. The `Commit` column is right: `git log --diff-filter=A` puts all five named files' first add at `aff894e`. No API key, bearer token, private-key block, `sk-ant`, `ghp_`, or email address in the file. | none

IMPORTANT | spec/ledger/7c6ea3af-...jsonl (5 occurrences) + spec/ledger/README.md:31-46 | The redaction MISSED the dash-mangled scratchpad path. `/tmp/claude-1000/-home-sera-Desktop-polimi-thesis-code-912days-qre-tool/...` survives 5 times, and it encodes the username and the entire directory chain the README says was replaced 671 + 7 + 20 times. `912days` appears 6 times, `polimi` 5, `qre_tool` 2. The substitutions matched only the slash-separated spelling. The README's "Audited before copying, and clean" and its `<home>` row are therefore claims the artifact does not support - the exact defect class CLAUDE.md's "Claims discipline" and AGENTS.md 8 forbid. | Redact `-home-sera-Desktop-polimi-thesis-code-912days-qre-tool` too, then re-verify JSON parse; or drop the sentence.

IMPORTANT | spec/ledger/*.jsonl | 10.54 MB for ONE 1-hour session, and git history is append-only: this cannot be undone later without a rewrite of every downstream clone. 9.97 MB of it is `user` records (n=153, ~65 kB each) - tool results echoed back. 5.10 MB of THAT is the `toolUseResult` field, which duplicates the `tool_result` block already present in the same record's `message.content`. So roughly half the payload is literal intra-record duplication. gzip -9 gets it to 6.0 MB, which is about what the packfile will carry, permanently. At this rate the ledger grows ~10 MB per working session. | Either strip the redundant `toolUseResult` field before committing (halves it at zero information loss), or keep the ledger in a separate repo / release asset and reference it. Do not let the first commit set an unbounded-growth precedent.

MINOR | spec/ledger/README.md:23-29 | Explains that the second, still-open transcript is deferred because "copying a live file would capture a mid-sentence snapshot". Reasonable, but nothing schedules the copy, and the same retention sweep R0.3.1 warns about applies to it. | Add it to INDEX.md's "Not yet indexed" as a dated TODO, or note the deadline.

MINOR | spec/ledger/README.md | 3.2 kB of prose for a 3-file directory, with a "Reading a transcript" section that re-explains newline-delimited JSON. R0.3.2 asks for three things: what the files are, why retained, and that they are the JOSS/thesis evidence base. Sections "What survived, and what did not" and "Redaction" earn their place; "Reading a transcript" (last 5 lines) does not. | Cut the "Reading a transcript" section.

### G6 - .claude/ now tracked

CRITICAL | .claude/.claude/** (11 files staged) | A nested duplicate of the agent configuration, mtime 2026-06-29, is staged for commit. `.claude/.claude/settings.json` is an obsolete copy whose permission set is WEAKER than the live one: it allows `Bash(git add:*)` and `Bash(git commit:*)`, which AGENTS.md 1 makes a hard prohibition ("Never run `git add`, `git commit`..."), and it lacks the current `Edit(data/**)` / `Edit(outputs/**)` denials. `.claude/.claude/agents/code-reviewer.md` (2362 B) and `plan-critic.md` (1659 B) are likewise stale, shorter forks of the live 4009 B / 1394 B files. Committing this publishes a second, contradictory, permissive agent config that someone can copy back over the real one. | `git rm --cached -r .claude/.claude` and delete the directory. It is a copy accident, not configuration.

IMPORTANT | .claude/review-findings.md (180 kB, 1824 lines before this section) | AGENTS.md 2, added in this same change set: "Prose belongs in `spec/`, never in `.claude/`, which holds configuration only." This file is 1800 lines of review prose - the largest single non-transcript addition, bigger than any source file in the repo. It is a working ledger for a review that is already finished, and its content (P5 findings on shape_stats, check_ledger, permutation nulls) is not configuration by any reading. | Move to `spec/` if the findings still matter, or leave it untracked. Add `.claude/review-findings.md` to `.gitignore`.

IMPORTANT | .claude/qre_checks_reference.tex (70 kB) | Two AGENTS.md 2 rules at once: it is prose in `.claude/`, and "Write specs and docs in Markdown, not LaTeX. What makes a spec work is numbered requirements ... not the markup." AGENTS.md also says thesis and paper LaTeX live outside this repository. It is cited by `.claude/review-findings.md` as "Authority for what it should be", so it IS load-bearing reference material - which is an argument for `spec/`, not for `.claude/`. | Move to `spec/` (as .md, or as-is if converting is too costly), or keep it out and cite the PDF the review already names.

MINOR | .claude/agents/*, .claude/commands/*, .claude/hooks/* | These 11 files are genuine configuration and belong under version control by AGENTS.md 2's own reading. No objection - listed so the CRITICAL above is not read as covering them. | none

MINOR | .gitignore:24 | Only `.claude/settings.local.json` is excluded. That is correct today (verified: it is the one `.claude` file absent from the staged set), but nothing stops the next `review-findings.md`-shaped file from being auto-staged. | Consider `.claude/*.md` and `.claude/*.tex` exclusions, or an allowlist of `agents/ commands/ hooks/ settings.json`.

### G4 - the vocabulary rename (R0.4)

VERIFIED | R0.4.5 acceptance | `grep -rn "repairable" --include="*.py" . | wc -l` returns **11**, not 0 and not 14. Sites: `panels/across_calibration.py` 1,4,8,11; `panels/within_calibration.py` 47; `jobs/active/mtbc_q6.py` 3; `jobs/active/mtbf_q1.py` 3; `tests/test_artifact_guard.py` 240,242,243,269. R0.4.4 PASSES outright: zero `availability` in `analyzers/ panels/ plots/ core/ jobs/`.

VERIFIED PASS | rename mechanics | All 8 file moves register as `R`/`RM` in porcelain, history preserved. All panel modules import and expose the new symbols; `panels/comparison.py` -> `CompareSeriesData`/`CompareSeriesPanel` with its only consumer (`jobs/composite/compare_t2star_0704_vs_1004.py:25,29,49,68`) updated; `python main.py inspect` on that job resolves the DAG. `tests/test_within_calibration_builder.py` + `test_across_calibration_builder.py` + `test_artifact_guard.py`: 22 passed. Full suite 330 passed / 2 skipped, unchanged from the pre-rename baseline as R0.4 requires.

JUSTIFIED (7 of 11) | panels/across_calibration.py:1-14, panels/within_calibration.py:47, jobs/active/mtbc_q6.py:3, jobs/active/mtbf_q1.py:3 | Decision D2b holds for these. The mapping note is real, argues the distinction (calibration boundary vs restoration-after-failure), cites Ascher and Feingold / Rigdon and Basu, and the three other sites cross-reference it. Same for `README.md:5`, `AGENTS.md:4`, `CITATION.cff:7` ("repairable-systems statistics", the field's name) and `docs/iid_checks/README.md:116` (a paper TITLE - renaming that would falsify a citation). | none

IMPORTANT | AGENTS.md:87 vs panels/across_calibration.py:11-14 | AGENTS.md 4, "Locked vocabulary", says flatly: "Do **not** use `repairable` / `non_repairable`. They are being removed." D2b says the term survives in prose that cites the field. Both files are tracked in THIS commit and they contradict each other on the one rule AGENTS.md marks as a scientific error to violate. A contributor reading AGENTS.md would delete the mapping note. | Add the D2b carve-out to AGENTS.md 4 in one line: kept only in prose naming the literature, nowhere in identifiers.

IMPORTANT | docs/PANEL_CONTRACT.md:274 | "the standard **non-repairable** degradation analysis surface" - this is describing OUR panel's coverage, three lines below `WithinCalibrationPanel` and inside the normative panel contract, with no mapping note. Unlike the 11 code sites it does not cite the field; it just uses the retired name for our own tier. The rename sweep touched this file (three `+/-` pairs) and missed this line because it is hyphenated prose, not the `non_repairable` identifier. This is the laziness case, not a D2b case. | Change to "within-calibration degradation analysis surface".

IMPORTANT | docs/PANEL_CONTRACT.md (same paragraph region) | The doc says "`panels/within_calibration.py` is ~650 lines, exceeding the project's 200-line ..." - the rename updated the PATH on that line and left the NUMBER. Measured: **1051 lines**, 5x the stated budget and 60% above the figure the doc quotes. CLAUDE.md's own rule: a measured number in a doc must come from the artifact in the state it ships. | Update to 1051, and treat the file as the flagged-unwieldy one (`panels/_within_calibration_compute.py` is a further 492).

MINOR | tests/test_artifact_guard.py:243 | Quotes the spec ("R0.4.5 of SPEC 0001 forbids the word `repairable` remaining in the source tree") and in doing so is itself one of the 11 hits that make R0.4.5 fail. Lines 240 and 269 are forced - the pickle stream literally contains `panels.non_repairable...` and the assertion must match it - but 243 is avoidable prose. | Say "the rename spec" without the word.

MINOR | panels/comparison.py | `CompareNonRepairableData/Panel` became `CompareSeriesData/Panel`, not `CompareWithinCalibration*`. R0.4.2 gives one mapping and says to flag any site that does not fit rather than renaming it. Choosing a third name is defensible (the panel really is generic) but it is an unrecorded deviation from the table. | Note it in the checkpoint banner.

MINOR | R0.4.5, non-.py scope | R0.4.5 says "no occurrence ... outside `spec/`, `docs/adr/`, and `spec/ledger/`" - not "no occurrence in .py". Only the acceptance LINE is .py-scoped. 6 non-.py occurrences remain (README.md:5, AGENTS.md:4,87, CITATION.cff:7, docs/iid_checks/README.md:116, docs/PANEL_CONTRACT.md:274). Five are justified; PANEL_CONTRACT.md:274 is not. Worth stating so the deviation is scoped honestly. | Report the non-.py set in the banner too.

### G7 - correctness spot-checks

IMPORTANT | tests/test_artifact_guard.py:231-274 | The widened guard is not merely looser, it is now **vacuous for its stated purpose**. Measured over the 72 real candidates the test globs (`output/*/subjobs_output/*/t2star_panel_data.pkl` + `output_backup2`): **72/72 raise ModuleNotFoundError, 0 raise ValueError.** The `ValueError` branch - the only one that exercises `StaleArtifactGuard` on a real artifact - is unreachable and will stay unreachable, because any pickle written from here on is post-split by construction. The docstring says "TWO failure modes now" as if both fire; only one does. The test would pass with `StaleArtifactGuard.__setstate__` deleted. (Not CRITICAL: lines 92-205 give the guard synthetic coverage including `_StalePickle` through the real `pickle.load` path, so the guard itself is still tested.) | Either assert that at least one ValueError was seen and regenerate one pre-split artifact under the new module name, or rename the test to what it now checks and stop calling it a guard test.

IMPORTANT | tests/test_bench_isolation.py:55-97 and 174 | Mutation-tested, three planted files, each removed after. (A) `import jobs.bench.runner` in `panels/` -> FAILS `test_pipeline_packages_do_not_import_the_bench[panels]`. Correct. (B) `importlib.import_module("jobs.bench.runner")` in `jobs/active/` -> FAILS `[jobs]`. Correct. (C) **`from ..bench import runner` in `jobs/active/` -> 19 passed. The guard MISSES it.** `_imported_roots` discards `node.level`, so it yields the bare string `"bench"`, and `_imports_bench` compares against `"jobs.bench"` - no match. That import resolves at runtime to `jobs.bench.runner`: a real, working bench import the guard admits. This is a REGRESSION introduced by this diff: under the old top-level-name matching, `"bench"` would have matched. Worse, line 174 of the positive control lists exactly this form, `("from ..bench import runner\n", "bench")`, and passes - because the control asserts on `_imported_roots`, not on `_imports_bench`, so it certifies a form the predicate does not catch. This is the CLAUDE.md instance "the AST walk missed importlib/relative imports its own docstring named" recurring in the same file. | Resolve `node.level` against the file's package before adding, or add `name == "bench" or name.startswith("bench.")` to `_imports_bench`; and re-point the positive control at `_imports_bench`.

VERIFIED PASS | pytest.ini:14 | Confirmed by measurement, not assertion: stripping `pythonpath = .` and re-running gives `Interrupted: 24 errors during collection`; restoring it gives 330 passed / 2 skipped. It does not mask a packaging defect that exists today, because no `pyproject.toml` exists and SPEC 0001 explicitly forbids adding one. See the G1 MINOR for the SPEC 0002 hazard it creates.

MINOR | tests/test_bench_isolation.py:113-115 | `_python_files` drops any path with `"bench"` as a path component. Correct for `jobs/bench/`, but it would also silently exclude a future `analyzers/bench/` or `plots/bench/` from every scan. | Anchor on the relative prefix `jobs/bench` rather than the component name.

### G5 - project identity files and the em dash sweep (R0.5)

VERIFIED PASS | em dash sweep | `grep -rlP '\x{2014}'` over the whole tree (excluding .git and output*) returns exactly ONE file: `spec/ledger/*.jsonl`, 79 hits, which must not be edited because it is a verbatim record. Every source file, doc and new prose file is clean. Also clean of every banned word in AGENTS.md 8 (surfacing / brings into view / data-driven / delve / leverage / firstly / excellent) - the only hits are AGENTS.md quoting its own rule.

VERIFIED PASS | the non-em-dash half of the unstaged diff | Isolated it mechanically: every remaining `+/-` pair outside the em-dash swap is ruff line-wrapping forced by the longer names (`WithinCalibration` > `NonRepairable`), e.g. `panels/across_calibration.py` `_draw_scatter/_draw_rolling/_draw_histogram` signatures and `jobs/active/mtbf_q1.py` import block. No semantic change.

VERIFIED PASS | CODE_OF_CONDUCT.md | Verbatim Contributor Covenant 2.1 with correct attribution and the `[INSERT CONTACT METHOD]` placeholder properly resolved to `CITATION.cff`. 106 lines is the standard's own length. NOT AI bloat - this is the right artifact to copy rather than write. Nothing to cut. | none

VERIFIED PASS | CITATION.cff | 21 lines, valid CFF 1.2, parses as YAML, has all four required keys. Correctly omits `version`/`date-released`/`doi` because R0.5.4 forbids tagging. Minimal, not padded. | none

IMPORTANT | README.md:56 | "Requires Python 3.11 or newer, plus numpy, scipy, pandas, matplotlib, plotly and pandera." **Four dependencies are missing.** The pipeline also imports `allantools` (`analyzers/allan.py:6`), `sklearn` (`analyzers/tlf.py:5`), `yaml` (`loaders/registry.py:7`) and `joblib` (`jobs/bench/runner.py:33`). Anyone following the install block on a clean environment hits ImportError. Installability is the first thing a JOSS reviewer checks, and it is the whole point of R0.5. | Add allantools, scikit-learn, PyYAML, joblib to the list.

IMPORTANT | README.md:16-17 and AGENTS.md:80 | "the loss measured here is roughly 25 to 45 per cent of real crossings at working thresholds". "Measured here" asserts the repo contains the measurement. `grep -rn "25 to 45|25-45|25 and 45"` over the whole tree returns exactly two lines: this one and AGENTS.md's. There is no artifact, no table, no test, no locator. It is the headline number of the Statement of Need - the paragraph a JOSS reviewer reads first - and it cites itself. AGENTS.md 8: "Never invent a section number, equation number, figure number, or citation. If you do not know the locator, write 'no source located'." | Either point at the artifact that measured it, or reword to what the code enforces ("QUEBRA never resamples, because resampling destroys and invents crossings") and drop the number.

IMPORTANT | CONTRIBUTING.md:47-48 | "Tests in `tests/statistical/` and `tests/jobs/` must name their oracle." Neither directory exists. `tests/` is flat: `tests/test_*.py` plus `tests/fixtures/`. AGENTS.md 5's six-tier table describes the same non-existent layout. A contributor following CONTRIBUTING.md would create directories nothing collects. | Say "statistical and end-to-end tests" until the tiers exist, or label the tier table as the target layout.

MINOR | CONTRIBUTING.md:16 | "the run directory name under `outputs/`" - this repo writes `output/` (singular). `.gitignore:36-39` says so explicitly: "`outputs/` is the name SPEC 0002 moves to; `output/` ... is what this repository writes today". A bug reporter is sent to a directory that does not exist. | `output/`.

MINOR | CONTRIBUTING.md:36 | Tells contributors "Run the gates before you push: `make lint`". `make lint` exits 1 at this commit (21 files fail `ruff format --check`). The first contributor who follows the instruction fails on code they did not write. Same root cause as the G1 IMPORTANT. | Format `tests/` at minimum.

MINOR | README.md:12 | "Doing that correctly is harder than it looks, for three reasons this toolkit addresses directly." Pure announcer sentence; the three paragraphs that follow each open with their own thesis and need no herald. It is also the one line in the file that reads generated. | Cut the sentence.

MINOR | AGENTS.md:4-6 | "applies repairable-systems statistics (Kaplan-Meier, Nelson-Aalen, log-rank, RMST, MCF) to them". `grep -rli "nelson|log_rank|logrank|rmst|mean_cumulative" --include="*.py" analyzers/` returns NOTHING. README.md:40 gets this right ("Nelson-Aalen, log-rank, RMST and MCF are **not** implemented yet"); AGENTS.md, tracked in the same commit, states the opposite. | Match the README's wording.

MINOR | README.md:1 vs repository-code | The project is named QUEBRA throughout but `git remote -v` and every URL in README/CITATION/CONTRIBUTING say `github.com/seraconti/qretool`. Precondition P2 records the repository as renamed; it is not, or not to this name. The README's clone command is at least CORRECT (it matches the remote), so this is cosmetic, but a reader meets two names for one project in the first 50 lines. | Rename the GitHub repo, or say once in the README that the historical slug is `qretool`.

MINOR | CITATION.cff:14 | A bare public email address in a machine-readable, indexed file. Sera's call, but it is the one irreversible personal-data item in this commit and there is no ORCID to substitute. | Consider ORCID plus the GitHub handle instead.

IMPORTANT | panels/within_calibration.py:6-9 | The mechanical rename rewrote a sentence whose entire subject was the OLD name, and inverted it into a falsehood: "That re-export is load-bearing and must not be removed: artifacts materialized before the split name `panels.within_calibration.WithinCalibrationPanelData` in their pickle stream". They do not. They name `panels.non_repairable.NonRepairablePanelData` - measured, 72/72 pickles under `output/` and `output_backup2/` raise ModuleNotFoundError on exactly that module. The docstring now claims a re-export protects artifacts that in fact cannot load at all, which is the opposite of what `tests/test_artifact_guard.py:239-247` records. Searched for siblings of this class (`pre-split`, `written before`, `used to be`, `formerly`); this is the only one. | Restore the old dotted name in that sentence and state that the re-export no longer covers pre-rename artifacts.

MINOR | R0.5 acceptance "No file contains an em dash" | Technically fails: `spec/ledger/*.jsonl` has 79. Editing a verbatim transcript to satisfy a style rule would be worse than the violation. | Scope the criterion to source and prose in the banner.

BLOCKING-BY-SPEC | "Done when" | `git ls-files` must show `.github/` tracked, and `.github/` does not exist. The same section says "Explicitly not in this phase: ... no CI." The spec contradicts itself: `.github/` has no non-CI content to hold. The R0.1 acceptance line only requires that `.github/workflows/ci.yml` NOT be ignored, which passes. | Sera's call to amend "Done when" to the acceptance line's weaker form; do not create `.github/` just to satisfy it.

---

# Phase 1 / SPEC 0002 review (installability)

SCOPE: whole repo, Phase 1 / SPEC 0002 (staged + untracked)
COMMIT: 9057803 (nothing committed; phase is staged/untracked)

## Manifest
- [x] src/quebra/core/paths.py
- [x] src/quebra/cli.py
- [x] pyproject.toml + quebra.toml
- [x] scripts/acceptance.sh + .github/workflows/ci.yml
- [x] docs/WRITING_A_JOB.md + README.md
- [x] conftest.py + tests/ deltas
- [x] src/quebra/recipes.py + __init__.py files
- [x] jobs/ import rewrites
- [x] baselines reproduction

## Reviewed
- jobs/ import rewrites
- baselines reproduction
- src/quebra/recipes.py + __init__.py files
- conftest.py + tests/ deltas
- docs/WRITING_A_JOB.md + README.md
- scripts/acceptance.sh + .github/workflows/ci.yml
- pyproject.toml + quebra.toml
- src/quebra/cli.py
- src/quebra/core/paths.py

## Findings

### src/quebra/core/paths.py

IMPORTANT | src/quebra/core/paths.py:3-7 | The module docstring, the first thing a reader meets, still says the repo root is "anchored off this file's location, the same pattern as provenance.get_git_commit". That is exactly what `repo_root` (:44-48) stopped doing in this phase - it walks up from `Path.cwd()`. The header also still calls the project `qre_tool/` and asserts the dataset root is "one level above" the repo, which is now only true because `quebra.toml` happens to say `data_root = ".."`. Header and function docstrings contradict each other on the single largest behaviour change in the phase. | Rewrite lines 3-7 to describe cwd-anchored marker discovery and the four-mechanism data root.

IMPORTANT | src/quebra/core/paths.py:150-156 | An explicit root that does not exist SILENTLY falls through to the environment, then `quebra.toml`, then platformdirs. In this checkout `quebra.toml` always resolves, so `--data-root /tpyo` yields `912days/` with no warning, and the run then hashes a different dataset tree than the operator asked for. AGENTS.md 3 makes fall-through the wrong direction ("errors are raised, not swallowed"); R1.3.4 only requires a raise when NOTHING resolves, it does not ask for this. `tests/test_data_root_resolution.py:96-107` pins the fall-through as intended behaviour, so the test locks the defect in. | Raise when `explicit is not None and not root.is_dir()`; keep the fall-through for the three non-explicit mechanisms.

IMPORTANT | src/quebra/core/runner.py:390 and src/quebra/cli.py:130-132 vs paths.py:139 | R1.3.5 says the `--data-root` flag maps to mechanism 1. It does not: the CLI resolves it itself and `run_job` does `Path(data_root).resolve() if data_root else default_dataset_root()`. `resolve_data_root(explicit=...)` therefore has NO production caller - mechanism 1 exists only for the test that exercises it. The two routes also disagree: the runner accepts a nonexistent root verbatim (and fails later at dataset resolution), the tested route falls through. | `run_job` should call `resolve_data_root(data_root)` so one chain governs.

MINOR | src/quebra/core/paths.py:22 | `PROJECT_MARKERS` includes `pyproject.toml` and `.git`, so ANY python project or git repository above the cwd is accepted as the QUEBRA project root - a git-tracked `$HOME`, a sibling project, a vendored subrepo. Measured: `cd /tmp && repo_root()` returns `/tmp`; `cd / ` returns `/`, after which `resolve_dataset_path`'s second candidate is `/<relative path>`. No live defect (the dataset fallback then simply misses and raises), but the marker set makes "some other project" indistinguishable from "this project". | Prefer `quebra.toml` as the sole positive marker, or document that the other two are heuristics.

MINOR | src/quebra/core/paths.py:65-92 | `resolve_dataset_path` is now cwd-dependent through its `repo_root()` fallback candidate, so which file gets content-hashed into the run identity depends on where the process was launched. Inside the checkout this is stable; outside it the candidate silently changes. The docstring argues the fallback (D4) but does not say it moved from `__file__` to the cwd. | One sentence in the docstring stating the fallback now follows the caller's project.

MINOR | R1.3.2 | `src/quebra/_fixtures/` was not created and nothing anywhere uses `importlib.resources` (`grep -rn "importlib.resources" src/` is empty). The requirement's first mechanism, packaged resources, is simply unimplemented. Defensible - there are no packaged resources yet - but it is an undeclared deviation, and the phase claims R1.3 complete. | Say so in the checkpoint banner, or add the empty package with a README.

### src/quebra/cli.py (and the death of main.py)

VERIFIED PASS | R1.1.4 | `grep -rn "sys.path" src/ tests/ | wc -l` = 0. The only hit anywhere in the tree is the word appearing in `conftest.py:16`'s prose. No `jobs.*` import in `jobs/active/` or `jobs/composite/` (only stdlib + `quebra.*`). `main.py -> src/quebra/cli.py` and `jobs/common.py -> src/quebra/recipes.py` both register as `RM` renames, so R1.1.2 holds.

VERIFIED PASS | src/quebra/cli.py:31-36 | The output root cannot land in site-packages: `_output_root` is `Path.cwd()/"output"` (or an explicit flag) and `_jobs_dir` is `Path.cwd()/"jobs"`. Nothing in the module touches `Path(__file__)`.

IMPORTANT | AGENTS.md:45-47 and AGENTS.md:200 | This change set DELETES `main.py` (renamed to `src/quebra/cli.py`) and, in the same uncommitted working tree, ADDS a Commands block telling every agent to run `python main.py run jobs/active/<job>.py`, `python main.py run --all` and `python main.py inspect`. All three now fail with "can't open file". Line 200's layout block still lists `main.py   CLI` while the block directly above it already gained the new `src/quebra/recipes.py` entry - a half-applied edit, the exact failure class this review was asked to sweep for. | Replace with `quebra run ...` / `quebra inspect ...`, and change line 200 to `src/quebra/cli.py`.

IMPORTANT | AGENTS.md:180-205 | The whole "Layout (real)" block still describes the pre-phase top-level layout (`core/`, `loaders/`, `transforms/`, `analyzers/`, `plots/`, `panels/`, `schemas/` at the repository root). After R1.1.1 every one of those lives under `src/quebra/`. The block is labelled "(real)", which is now false for eight of its ten entries. | Prefix the tree with `src/quebra/`.

IMPORTANT | jobs/composite/compare_t2star_0704_vs_1004.py:3,7,8 and jobs/composite/independence_survey.py:12,33 | Both files were edited in this phase (import rewrites) and both still document `PYTHONPATH=. python main.py run jobs/composite/<job>.py` as the way to run them. `main.py` no longer exists and `PYTHONPATH=.` is no longer needed. A composite is the ONE job family a reader cannot discover by `run --all`, so its docstring is the only instruction that exists. | `quebra run jobs/composite/<job>.py`.

IMPORTANT | src/quebra/__init__.py (empty) vs SPEC 0002 R1.3 acceptance | The spec's acceptance line is `QUEBRA_DATA_ROOT=/tmp/qd python -c "import quebra; print(quebra.__version__)"`. Measured from `/tmp`: `AttributeError: module 'quebra' has no attribute '__version__'`. `scripts/acceptance.sh:50-51` prints `quebra.__file__` instead, so the script passes while the stated criterion does not hold. | Add `__version__ = importlib.metadata.version("quebra")` to `src/quebra/__init__.py`, or amend the criterion.

MINOR | src/quebra/cli.py:119 | `--data-root` help still says "(default: the repo's parent directory)". After R1.3 the default is the four-mechanism chain; the repo's parent is merely what `quebra.toml` happens to declare. | "default: QUEBRA_DATA_ROOT, then [tool.quebra] data_root, then the user data dir".

MINOR | src/quebra/cli.py:134-139 | `quebra run --all` from a directory with no `jobs/active/` glob-matches nothing, prints nothing and exits 0. Silent success for a command that did no work, and this is now reachable because the CLI is anchored on an arbitrary cwd rather than on the repo. | Raise if `_jobs_dir("active")` does not exist, or if `job_files` is empty.

MINOR | src/quebra/cli.py:16-19, 190 | `h5py` is wrapped in `try/except ImportError: h5py = None`, but `pyproject.toml:37` makes `h5py>=3.11` a mandatory dependency and `loaders/registry.py` imports it unguarded. The guard can no longer fire, and if it somehow did, `schema-wizard` on an `.h5` file falls silently through to `load(file_path)` rather than saying h5py is missing. | Drop the try/except now that the dependency is declared.

### pyproject.toml + quebra.toml + Makefile

VERIFIED PASS | pyproject.toml:28-40 vs R1.2.2 | Derived independently with an AST walk over all 76 files under `src/quebra/`. The third-party import set is exactly {allantools, h5py, matplotlib, numpy, pandas, pandera, platformdirs, plotly, scipy, sklearn, yaml} - 11 names, and all 11 are declared (sklearn -> scikit-learn, yaml -> PyYAML). Nothing declared is unimported. `joblib` really is absent from `src/` (it appears only in `jobs/bench/runner.py:33`), so the comment at :25-27 is accurate. `make deps` (deptry over src/quebra) exits 0, 76 files scanned, no issues.

VERIFIED PASS | pyproject.toml:81-90 | `select = ["E4","E7","E9","F"]` is ruff's documented default, so this phase changes no lint outcome. Confirmed the claim at :84-85 is in the right ballpark and in the right direction: `ruff check .` under the shipped selection reports exactly 1 finding (the pre-existing F841 at `src/quebra/plots/interpolation_stage_plot.py:63`), and widening would land a reformat nobody reviewed. Correct call for an installability phase.

VERIFIED PASS | pyproject.toml:66-72, 92-103, 105-125 | R1.2.4 satisfied: `[tool.ruff]`, `[tool.pytest.ini_options]` with `addopts = "--strict-markers"` and all four markers, `[tool.hatch.build.targets.wheel] packages = ["src/quebra"]`. R1.2.6 satisfied: no `[tool.importlinter]`, and `make arch` correctly SKIPS with that reason (exit 0). R1.2.3 satisfied (`dev` has all seven named tools plus `build`; `r = []`).

VERIFIED PASS | LICENSE vs pyproject.toml:11 | `license = "GPL-3.0-or-later"` matches the GPLv3 text actually in `LICENSE`.

IMPORTANT | .github/workflows/ci.yml:33-34 | CI runs `ruff check .` as its second step. Measured at this commit: `ruff check .` exits **1** on the pre-existing `F841` in `src/quebra/plots/interpolation_stage_plot.py:63`. The workflow therefore CANNOT be green on the first push, which is R1.4.2's own acceptance criterion ("The CI workflow passes on a push"). The spec asks CI to run the fast selector; the lint step is an addition that fails. | Fix the F841 (it is one dead assignment) or drop the lint step to match the spec.

IMPORTANT | Makefile:18-20 | The skip-guard on `types` was DELETED, so `make types` now hard-fails (measured: 11 mypy errors, exit 2), and `make check` (`lint types arch test`) fails with it. AGENTS.md:50-52, edited in this same working tree, still tells agents `make types` and `make arch` "do not pass yet" and that a failure "tells you which phase you are in" - but a hard failure in `check` no longer distinguishes "wrong phase" from "you broke the build". Same for `lint`: `ruff check .` exits 1 and `ruff format --check .` reports 11 files. `make check`, the documented pre-checkpoint gate, cannot pass at this commit. | Restore the guard (or gate `types` on a spec marker the way `arch` now is), and say in AGENTS.md which of the four are expected red.

MINOR | Makefile:9 | `.PHONY` lists `check lint types arch deps test test-all test-r docs clean` but not `test-real`, which is a real target at :45. A file named `test-real` would shadow it. | Add it.

MINOR | Makefile:46 | `pytest -m "real or regression"` names a `regression` marker that `[tool.pytest.ini_options] markers` does not define. `--strict-markers` does not validate `-m` expressions, so this silently selects nothing rather than erroring. | Drop `or regression`, or register the marker.

MINOR | pyproject.toml:69-72 vs R1.2.5 | R1.2.5 says "Exclude from the wheel: spec/, docs/adr/, tests/, data/, outputs/". No wheel `exclude` is written; the exclusion is implicit in `packages = ["src/quebra"]`. That is stronger in practice, but the sdist `exclude` at :72 is explicit while the wheel one is not, so a reader checking R1.2.5 finds nothing where they look. | One comment saying the wheel excludes by construction.

MINOR | quebra.toml (untracked) | The file is the sole reason mechanism 3 resolves in this checkout and is currently UNTRACKED. If it is not staged with the rest of the phase, a fresh clone has no data root at all and every dataset job dies with `DataRootNotFound`. | Confirm it is added before the commit.

### scripts/acceptance.sh + .github/workflows/ci.yml

VERIFIED PASS | scripts/acceptance.sh | Ran it: exit 0, seven steps, all 0. The venv IS outside the repo (`$(mktemp -d)/venv`) and the wheel IS what is under test - verified independently by running `env -C <repo> <venv>/bin/python -c "import quebra; print(quebra.__file__)"`, which prints `/tmp/tmp.*/venv/lib64/python3.13/site-packages/quebra/__init__.py`, not `src/`. So the pytest step, despite running with the repo as cwd, exercises the installed distribution. The wheel itself is clean: 82 members, top level exactly `quebra/` and `quebra-0.1.0.dev0.dist-info/`, zero entries under `spec/ data/ outputs/ tests/ jobs/ output/`. R1.2 acceptance holds.

IMPORTANT | scripts/acceptance.sh:57-58 and .github/workflows/ci.yml:36-39 | **The fast selector deselects nothing.** Measured: `pytest -m "slow or heavy or real or r"` collects 0 tests and deselects all 339; `grep` finds no `@pytest.mark.{slow,heavy,real,r}` anywhere in `tests/`. The acceptance step and CI both report 337 passed / 2 skipped, identical to the full suite. So ci.yml:36-37's stated mechanism - "The `real` marker is what keeps those tests out; CI never sees the tree they need" - is not a mechanism at all. R1.4.3 passes today only because `data/real_private/` does not yet exist and no test happens to read a private tree; nothing enforces it. | Say the marker set is declared but unused, or mark the tests that need `tool/datasets/` (`tests/test_paths_dataset_fallback.py` is the closest) before SPEC 0003 creates the private tree.

IMPORTANT | pyproject.toml:69-72 (sdist exclude) | The sdist ships `.claude/` - all of it, including `review-findings.md` (2000+ lines of internal review notes on unreleased statistics) and `qre_checks_reference.tex`, plus `settings.json`, the hooks and the agent definitions. Verified against the built `dist/quebra-0.1.0.dev0.tar.gz`. `spec/ledger/*.jsonl` was correctly excluded; `.claude/` was not considered. | Add `.claude/` to `[tool.hatch.build.targets.sdist] exclude`.

IMPORTANT | .github/workflows/ci.yml | CI does not run `scripts/acceptance.sh`, does not run `python -m build`, and installs from the source tree with the repository present. The single claim this phase exists to make - a wheel installs and works with no checkout - is therefore verified only by a script a human has to remember to run. R1.4 acceptance lists the two as separate criteria so this is not a spec violation, but "CI is green" will not catch a packaging regression. | Add one step running `bash scripts/acceptance.sh`, or say in the banner that CI does not cover packaging.

MINOR | scripts/acceptance.sh:37 | `WHEEL="$(ls -t "${REPO}"/dist/*.whl | head -1)"` picks the newest wheel in a directory that ACCUMULATES. If the build step fails, the script keeps going and installs a stale wheel from a previous run; every later step then passes and only `STATUS=1` from the failed build makes the overall exit nonzero. The steps after a failed build prove nothing about the current tree while printing "exit 0". | `rm -rf "${REPO}/dist"` before building, or exit immediately when the build step fails.

MINOR | scripts/acceptance.sh:44-46 | The venv is created but never removed, and neither is the `mktemp -d` parent. Each run leaks a full venv (numpy, scipy, matplotlib, plotly) into `/tmp`. | `trap 'rm -rf "$(dirname "${VENV}")"' EXIT`.

MINOR | scripts/acceptance.sh:15 | `set -u` without `set -o pipefail`. The one pipeline (`ls -t ... | head -1`) is guarded by the emptiness check at :38, so nothing is masked today. | Note only.

MINOR | .github/workflows/ci.yml:6-8 | No `concurrency` group and no `permissions:` block. Every push to every branch starts a full 2-interpreter matrix, and the job runs with the default (write) token. | `permissions: contents: read` and a concurrency group.

### docs/WRITING_A_JOB.md + README.md

VERIFIED PASS | docs/WRITING_A_JOB.md:13-40 | The example job is real, not illustrative. Wrote it verbatim to a scratch directory outside the repo and ran `quebra inspect my_first_job.py`: exit 0, two nodes printed. Every API it names exists with the signature shown - `Dataset(path, schema, qubit, device)` (core/dataset.py:9-15), `Job("name")` positional (core/job.py:279), `track912Schema` (schemas/track912.py:9), `t2star.run` / `t2star.make_inputs_from_norm` (analyzers/t2star.py:46,72), `job.load` / `job.step` / `job.materialize` / `job.figure` / `job.include(..., alias=)` / `.ref`. The run-directory format at :87 matches `runner.py:428` exactly (`<job>_<identity6>_<timestamp>`). The data-root order at :63-69 matches `resolve_data_root`. 124 lines, single-purpose: not doc bloat.

IMPORTANT | docs/WRITING_A_JOB.md:8 | "`pip install quebra`, put a `.py` file anywhere, and run it by path." Measured: `pip index versions quebra` -> "No matching distribution found"; `https://pypi.org/pypi/quebra/json` -> 404. The distribution is unpublished (`version = "0.1.0.dev0"`, and SPEC 0001 R0.5.4 forbids tagging). The first instruction in the new doc is a command that cannot succeed, and it is the one a JOSS reviewer would try first. | "`pip install .` from a checkout, or `pip install <wheel>`" until the name is claimed on PyPI.

IMPORTANT | README.md:105-106 vs scripts/acceptance.sh:35-58 | R1.4.4's acceptance is that the README install commands and the script "agree line for line", and README:107-109 asserts flatly "Nothing is documented here that the script does not do". They differ on four points: (a) the script builds a wheel with `python -m build` and installs THAT, the README says `pip install .`; (b) the script runs `pip install pytest` as its own step, the README omits it; (c) the script runs `pytest -m "not slow and not heavy and not real and not r" -q tests/`, the README says `pytest tests/`; (d) the script's venv is outside the repo, the README's `.venv` is inside it. (b) is the one that bites: after a plain `pip install .` pytest is NOT installed (it lives only in the `dev` extra), so a reviewer following the README literally gets "command not found" on the next line. | Add `pip install pytest` (or `pip install ".[dev]"`) to the README block, and soften the "line for line" claim to what is actually shared.

IMPORTANT | README.md:97-98 | The clone command changed from `https://github.com/seraconti/quebra.git` to `git@github.com:seraconti/qretool.git`. The SSH form requires a key registered with GitHub; a reviewer on a machine that has never seen this repository - the governing constraint of the whole spec - cannot run it. The repository-name mismatch (quebra vs qretool) was already logged in the Phase 0 section; this change fixes the name and breaks the protocol. | `git clone https://github.com/seraconti/qretool.git`.

MINOR | README.md:87-91 | The "Feedback is welcome, particularly from anyone who runs long characterization campaigns..." sentence appears twice: once inside the Status paragraph at :87 and again as its own paragraph at :89-91. Pre-existing (untouched by this diff), but the README is in R1.4.4's scope and this is in the section a reviewer reads. | Delete one.

MINOR | README.md:141-142 | The Documentation section lists "time and clock semantics, the panel contract, and the figure standard" and does not mention `docs/WRITING_A_JOB.md`, the doc this phase added and the only entry point for someone writing their own job. | Add it.

MINOR | docs/WRITING_A_JOB.md:31 | The showcase example wires its step as a `lambda`, nine lines after the doc explains that "anything that changes the answer belongs in a kwarg rather than a closure" (:50-52). It works (verified), but the first job a reader writes is now a closure, and `quebra inspect` renders the node as `<lambda>` with no kwargs, which is exactly the opaque graph the doc sells inspect against. | Use a named `def`.

MINOR | docs/WRITING_A_JOB.md:24 | `path="tool/datasets/6D2S/070423_6D2S_qubit1.pickle"` hard-codes this researcher's private tree into the generic getting-started doc, and SPEC 0003 replaces `tool/datasets/` with the three-way `data/` split. A reader without those files gets a FileNotFoundError from step one. | Use an obviously-placeholder path, e.g. `my_device/run_001.pickle`.

### conftest.py + tests/ deltas

VERIFIED PASS | baselines | Reproduced all of them. `pytest -q` -> **337 passed, 2 skipped**. `pytest --collect-only -q | tail -1` -> **339 collected**; 339 - 7 (the new `tests/test_data_root_resolution.py`) = **332**, the P2 baseline, so D6(b)'s collect-count condition holds. `ruff check .` -> exactly 1 finding, the pre-existing F841 at `src/quebra/plots/interpolation_stage_plot.py:63`. `bash scripts/acceptance.sh` -> exit 0, seven steps, all 0. `mypy src/quebra/core` -> 11 errors in 2 files (out of phase). `make arch` -> SKIPPED, exit 0, correct under R1.2.6. `make deps` -> exit 0.

VERIFIED PASS | tests/test_bench_isolation.py:80-96, 105-116 | The relative-import hole logged in the Phase 0 section is genuinely FIXED, confirmed by mutation against `_imports_bench` on planted files: `from ..bench import runner`, `from .bench import runner` and `from ..bench.carve import carve_windows` all now return True, while `from quebra.recipes import ...` and `import jobs` return False. `node.level` is recorded and matched.

IMPORTANT | tests/test_bench_isolation.py:76-79 | The SAME class of hole remains one line away: for a non-relative `ImportFrom` only `node.module` is recorded, never `module + "." + alias`. Measured: `from jobs import bench` yields `{'jobs'}` and `_imports_bench` returns **False**. That is an ordinary, working absolute import of the bench that the guard admits. AGENTS.md 4: "Fix every site of a class, not the one that failed" - the relative half was fixed and the absolute half was not. | `roots.update(f"{node.module}.{a.name}" for a in node.names)` for the `node.level == 0` branch, and add `from jobs import bench` to the control table.

IMPORTANT | tests/test_bench_isolation.py:41-50 | `PIPELINE_PACKAGES` names the seven `src/quebra/` subpackages but not the package's TOP-LEVEL modules. `cli.py` gets its own assertion at :148-149; `src/quebra/recipes.py` and `src/quebra/provenance.py` get none. `recipes.py` is the file this phase moved INTO the distributed package, so it is precisely the new module that must not depend on the study, and it is now the one that is unscanned. | Replace the seven entries with `"src/quebra"`; rglob covers all of them and the top-level modules too.

IMPORTANT | tests/test_artifact_guard.py:277-290 | The ModuleNotFoundError assertion was WIDENED from `"non_repairable" in exc or "repairable" in exc` to a five-token any() including the bare strings `"panels"`, `"analyzers"` and `"core"`. D6(b) says no assertion may change; this one did, and in the direction that weakens it. Concretely: if the wheel ever fails to ship `quebra/panels/`, unpickling raises `No module named 'quebra.panels'`, the token `"panels"` matches, and the test passes - the packaging break that this whole phase exists to prevent is now indistinguishable from the accepted rename loss. The comment two lines above claims the opposite ("this clause cannot swallow an unrelated packaging break"). | Match on the REMOVED module names only: `panels.non_repairable`, `panels.repairable`, or a bare `panels.`/`analyzers.`/`core.` prefix that cannot appear under `quebra.`.

MINOR | tests/test_bench_isolation.py:180-208 | The positive control still asserts `expected in _imported_roots(...)` rather than `_imports_bench(...)`, so it certifies the collector, not the predicate. This was already logged in the Phase 0 section and carried forward unchanged. It is why the `from jobs import bench` miss above is invisible: the control would happily accept `"jobs"`. | Assert `_imports_bench(planted)` as well.

MINOR | tests/test_path_resolution.py:62-66 | `test_resolve_repo_path_is_cwd_independent` is now a TAUTOLOGY. `resolve_repo_path(p)` is `(repo_root()/p).resolve()` and the assertion's right-hand side is `repo_root()/"jobs"/"active"`; both sides read the same cwd-derived value, so it holds for every cwd including the `tmp_path` it chdirs to. Under the old `__file__` implementation it was a real assertion. The name and the change it was written to guard are now opposites of the shipped behaviour. | Rename to what it checks, or pin against an absolute path captured before the chdir.

MINOR | tests/test_path_resolution.py:31-39, tests/test_data_root_resolution.py:110-116 | Three tests now depend on the process cwd being inside this checkout: `repo_root()` must find a marker, and `default_dataset_root()` must resolve through the repo's `quebra.toml`. Run the suite from anywhere else and they fail with `DataRootNotFound`. `scripts/acceptance.sh:57` avoids this only because it passes `env -C "${REPO}"`, which is also what stops the "clean environment" step from being as clean as its comment claims. | Note the cwd requirement in the module docstring so a future runner does not chase it.

MINOR | tests/test_data_root_resolution.py:63,78 | Raw `os.chdir(deep)` inside two tests instead of `monkeypatch.chdir`. Safe today only because the `isolated` fixture already called `monkeypatch.chdir`, whose undo restores the pre-test cwd; drop the fixture from either test and the suite leaks a cwd into every test that follows. | Use `monkeypatch.chdir`.

MINOR | conftest.py | Docstring-only, no code, which is the whole trick (pytest's prepend import mode puts the rootdir on `sys.path`). It works - verified, `jobs.*` and `tests.fixtures` both import - but the mechanism is entirely implicit and one `--import-mode=importlib` away from silently breaking. | Consider `sys.path` being explicit here instead; R1.1.4 scopes its ban to `src/` and `tests/`, and this file is neither.

### src/quebra/recipes.py + __init__.py files + jobs/ import rewrites

VERIFIED PASS | R1.1.1 / R1.1.2 / R1.1.3 | `git status --porcelain` = 136 entries: 71 `R`/`RM` renames, 0 `A`, 1 unstaged `D` (pytest.ini). No file appears as add-plus-delete. All 9 directories under `src/quebra/` carry an `__init__.py` (9 of 9, matching the acceptance criterion), plus `jobs/`, `jobs/active/`, `jobs/composite/`, `jobs/bench/`. Every module under `quebra.` imports cleanly (`pkgutil.walk_packages` over all 76: zero failures), and every file in `jobs/active/`, `jobs/composite/` and `jobs/bench/` imports cleanly.

VERIFIED PASS | jobs/bench/instrument_report.py, jobs/bench/xi_ties.py | The two `sys.path.insert(0, ...parents[2])` lines were REMOVED, not preserved, as R1.1.4 demands. Both documented invocations still work after the removal: `python jobs/bench/instrument_report.py` regenerates `jobs/bench/results/instrument_report.md` byte-identically (checked - `git status` on that path stays empty), and `python jobs/bench/xi_ties.py` starts its Monte Carlo. Zero un-rewritten imports anywhere: `grep` for a top-level `from|import (analyzers|core|panels|plots|loaders|schemas|transforms|provenance)` over every `.py` in the tree returns nothing.

IMPORTANT | src/quebra/provenance.py:23, :46 | `get_git_commit` and `is_tree_clean` still anchor on `Path(__file__).resolve().parent`. Before this phase that was the repository root; now it is `<...>/site-packages/quebra/`. So the two roots in one run disagree by construction: `paths.repo_root()` follows the caller's cwd while provenance asks git about wherever the LIBRARY is installed. Installed from a wheel this yields `"nogit"` and `is_tree_clean() == False` (reuse permanently disabled, silently); installed editable inside some other checkout it records THAT repository's commit into `.prov.json`. Provenance is a hard rule in AGENTS.md 1 and this is the one field a reader uses to reproduce a figure. | Anchor both on `paths.repo_root()`, and record explicitly when no repository was found rather than the string "nogit".

MINOR | src/quebra/recipes.py:55-77 | `run_start_unix_s_from_hdf5` opens an HDF5 file. It is called at job-build time rather than inside a step, so it does not violate "a step never reads disk" - but this phase moved the function INTO the distributed package on the grounds that it is reusable library code, and AGENTS.md 3 puts disk access in `loaders/`. The package now has a second I/O site outside the loader registry. | Either move it behind the loader registry or say in the docstring why it is exempt.

MINOR | src/quebra/analyzers/reliability_band.py:111, distinguish_band.py:143, check_ledger.py:133 | Three user-facing error messages instruct the reader to "construct via `analyzers.reliability_band.run()`" (and siblings). Those dotted paths no longer import after R1.1.1; the working form is `quebra.analyzers.*`. An actionable error message that gives an unimportable path is worse than none. About a dozen more docstring references carry the old prefix, but these three are the ones a user is told to type. | Prefix the three with `quebra.`.

MINOR | tests/test_checks_cvm.py:163, jobs/bench/probe_unresolved.py:20 | Two more run instructions that do not run: `python jobs/bench/runner.py` fails with `ModuleNotFoundError: No module named 'jobs'` (its `from jobs.bench.arms import ...` needs the repo root on the path; `python -m jobs.bench.runner` works), and probe_unresolved's usage line still says `PYTHONPATH=. python bench/probe_unresolved.py` from before the bench moved under `jobs/`. Both pre-date this phase, but both files are in this change set. | `python -m jobs.bench.runner`, `python jobs/bench/probe_unresolved.py`.

MINOR | the index vs the working tree | 48 files are `RM`: staged as renames, with the import rewrite left UNSTAGED. The staged tree alone does not import - `git show :src/quebra/recipes.py` still says `import analyzers.fidelity`. `pytest.ini`'s deletion is also unstaged (` D`), and `quebra.toml`, `conftest.py`, `pyproject.toml`, `.github/`, `scripts/`, `docs/WRITING_A_JOB.md`, `tests/test_data_root_resolution.py` and six `__init__.py` files are untracked. A `git commit` without `git add -A` produces a broken commit. | `git add -A` before committing; the phase is one commit, not two.

### baselines reproduction and the Done-when clause

CRITICAL | SPEC 0002 "Done when" vs scripts/acceptance.sh:57 | **The fast test suite does NOT pass from an unrelated working directory.** The Done-when clause is: "`pip install .` succeeds in a fresh virtualenv outside the repository, `import quebra` works from an unrelated working directory, the fast test suite passes THERE, and CI is green." Measured, using the very venv the acceptance script built:

```
cd /tmp && <venv>/bin/python -m pytest -q -m "not slow and not heavy and not real and not r" <repo>/tests
6 failed, 331 passed, 2 skipped
```

The six: `test_data_root_resolution.py::test_the_repository_resolves_through_its_own_quebra_toml`, `test_path_resolution.py::test_repo_root_is_the_tool_repo`, `::test_default_dataset_root_is_repo_parent`, `test_paths_dataset_fallback.py::test_a_tracked_in_repo_table_resolves`, `::test_without_the_fallback_that_path_would_not_exist`, `::test_an_absolute_path_is_still_existence_checked`. All six die on `DataRootNotFound` or a missing marker, because `repo_root()` and `resolve_data_root()` now both read `Path.cwd()`.

`scripts/acceptance.sh:57` runs that step with `env -C "${REPO}"`, so the one cwd under which these six pass is the one the script chooses. The script's own header claims it "is the only thing that can falsify" the installability claim; on this point it does the opposite. This is not a packaging bug - the wheel is correct - it is the phase's headline claim being checked under the single condition that makes it true. | Either fix the six tests to be cwd-independent (capture an anchor from `Path(__file__)` in the TEST, which is legitimate there), or amend Done-when to say the suite runs from the checkout because `tests/` is not shipped. Do not leave the current wording and the current script together.

VERIFIED PASS | baselines | `pytest -q` -> 337 passed, 2 skipped. `pytest --collect-only -q | tail -1` -> 339; minus the 7 new tests in `tests/test_data_root_resolution.py` = 332, the P2 baseline, so D6(b) holds on the count. `ruff check .` -> 1 pre-existing F841. `bash scripts/acceptance.sh` -> exit 0, seven steps. `mypy src/quebra/core` -> 11 errors. `lint-imports` -> skipped, no contract. `make deps` -> exit 0. Wheel: 82 members, `quebra/` + dist-info only.

VERIFIED PASS | src/quebra/core/paths.py:31 | "the four jobs that declare `jobs/bench/results/*.csv` as a Dataset" - exactly four: `jobs/active/check_calibration.py`, `jobs/active/instrument_validation.py`, `jobs/composite/check_ledger_q1.py`, `jobs/composite/independence_survey.py`. Accurate.

VERIFIED PASS | tests/test_artifact_guard.py:246-250 | "77 pickles under `output/` and `output_backup*/` are affected, plus 5 naming `panels.repairable`" - measured over all 1012 pickles in the three output trees: exactly 77 contain `panels.non_repairable` and exactly 5 contain `panels.repairable` without it. Accurate. (0 contain `quebra.panels`, which confirms the `ValueError` branch of that test is still unreachable - already logged in the Phase 0 section.)

VERIFIED PASS | .github/workflows/ci.yml:16-18 | The 3.11 claim checks out: `ast.parse(..., feature_version=(3,11))` over every `.py` in `src/`, `jobs/` and `tests/` reports zero syntax failures, so the nested f-string it cites is genuinely gone and the 3.11 matrix leg will not die on import.

MINOR | pyproject.toml:84-85 | "Measured on this tree, `select = ["E","F","I","UP","B"]` reports 1193 findings, 1116 of them line-too-long." Measured now: **1212 findings, 1138 E501**. The numbers were taken before the phase finished adding files. AGENTS.md 4 names this exact defect class ("a measured number in a docstring must come from the artifact it cites, in the state it ships"). The argument is unaffected; the numbers are wrong. | 1212 / 1138, or drop the counts and keep the reasoning.

MINOR | AGENTS.md:213 | "`monoliths/`, `scripts/` and `jobs_old/` were deleted and no longer exist." This phase CREATED `scripts/`, and `scripts/acceptance.sh` is R1.4.1's deliverable. Same half-applied-edit class as the `main.py` lines. The layout section also never gained `pyproject.toml`, `quebra.toml`, `conftest.py`, `scripts/` or `.github/`. | Drop `scripts/` from that sentence and add the five new root files to the layout.

---

## db61799 slice review - generic loading path / schema dispatch

SCOPE: src/quebra/core/job.py (_load_dataset + _load_dataframe_raw only), src/quebra/core/dataset.py, src/quebra/schemas/ramsey_series.py, src/quebra/schemas/null.py, tests/test_load_dataset_contract.py, docs/WRITING_A_SCHEMA.md
COMMIT: db61799 (reviewed at HEAD b14bbc7)

NOTE: a second reviewer was appending to this ledger concurrently, so this slice's
entries sit at the file tail interleaved with theirs rather than under the heading
below. Every entry of mine names one of the six scoped paths above; grep those.
STATUS: complete, all six files reviewed.

### Manifest
(empty - all reviewed)

### Reviewed

### Findings

---

## b14bbc7 slice review - DataUnavailable / packaged fixtures / data manifest

SCOPE: src/quebra/core/paths.py (new parts only: DataUnavailable, _is_manifested,
_MANIFEST_RELATIVE, hoisted tomllib, changed raise in resolve_dataset_path),
src/quebra/_fixtures/__init__.py, scripts/make_fixtures.py, scripts/make_data_manifest.py,
tests/test_data_unavailable.py, tests/test_packaged_fixtures.py, tests/test_data_manifest.py
COMMIT: b14bbc7

### Manifest
(reviewed) src/quebra/core/paths.py
(reviewed) src/quebra/_fixtures/__init__.py
(reviewed) scripts/make_fixtures.py
(reviewed) scripts/make_data_manifest.py
(reviewed) tests/test_data_unavailable.py
(reviewed) tests/test_packaged_fixtures.py
(reviewed) tests/test_data_manifest.py

### Reviewed

### Findings

IMPORTANT | src/quebra/core/paths.py:86 | `text.endswith(str(record.get("path","\0")))` is a raw STRING suffix match, so it matches across a path component boundary. Measured: with record `path = "6D2S/x.pickle"`, `_is_manifested` returns True for `data/other/evil_6D2S/x.pickle` AND for `/tmp/zzz6D2S/x.pickle`. Both are wrong paths that the new message will call "an EMBARGOED record rather than a wrong path" - the exact inversion the class exists to prevent. | Anchor on a separator: `text == p or text.endswith("/" + p)`.

IMPORTANT | src/quebra/core/paths.py:79-86 | The docstring promises the manifest read "must degrade to 'not manifested' rather than replace a missing-dataset message with a manifest-parsing one", but `except (OSError, tomllib.TOMLDecodeError)` does not cover the shapes that actually occur. Measured, all RAISING out of the error path: `[record]` written as a table instead of `[[record]]` (a one-bracket typo) -> AttributeError 'str' has no attribute 'get'; `record = "hello"` or `record = [1,2]` -> AttributeError; a manifest containing a non-ASCII filename read under `LC_ALL=C` -> UnicodeDecodeError (a ValueError, not OSError). Any of these replaces the useful FileNotFoundError with a traceback from the diagnostic helper. | `records = ...` then `if not isinstance(records, list): return False`, keep only `isinstance(r, dict) and isinstance(r.get("path"), str)`, and read with `manifest.read_bytes()` + `tomllib.loads(... .decode("utf-8"))` or `tomllib.load(open(...,"rb"))`.

MINOR | src/quebra/core/paths.py:86 | The `"\0"` sentinel does work (no POSIX path can contain NUL, so a record with no `path` never matches) but it is load-bearing obscurity, and `str()` around it silently accepts a non-string `path`: measured, `path = 5` makes `_is_manifested(Path("foo5"))` return True. | Drop the sentinel in favour of the isinstance filter above.

MINOR | src/quebra/core/paths.py:127-171 | `DataUnavailable.__init__` takes four required positional args but calls `super().__init__(one_string)`, so `self.args == (msg,)`. OSError's `__reduce__` reconstructs as `type(self)(*self.args)`, which raises TypeError - the exception cannot be pickled or copied, and multiprocessing/joblib propagation would surface a TypeError instead of the real error. | `self.args = (raw, candidates, dataset_root, manifested)` is wrong for str(); simplest is to accept defaults for the last three params.

MINOR | src/quebra/core/paths.py:69 | `_is_manifested` reads the manifest from `repo_root()`, not from `dataset_root`, so a reviewer running with `--data-root` pointed at a private tree that carries its own MANIFEST.toml gets "not in data/real_private/MANIFEST.toml" for a record that is in one. Defensible (the manifest is committed repo metadata) but undocumented. | One sentence in the docstring saying the manifest is always the checkout's.

IMPORTANT | tests/test_packaged_fixtures.py:29-36 | `test_the_fixture_path_is_not_computed_from_the_source_tree` compares `fixture_path(name).parent` with `module.fixture_path(name).parent` - the SAME function object called twice, second call served from `_RESOLVED`. It asserts `x == x` and cannot fail for any implementation, including one written as `Path(__file__).parent / name`. The claim in its own docstring ("it must come from the imported package") is unguarded, and this is the bug class `repo_root()` was already burned by. | Assert the path is under `resources.files("quebra._fixtures")`, or that it survives with `quebra` imported from a copied install tree.

MINOR | src/quebra/_fixtures/__init__.py:3-5 | "point a `Dataset` at `fixture_path(...)` and the whole pipeline runs" is not true as shipped: `ramsey_synthetic.csv` carries no `DDMMYY_` prefix and its timestamps start at 0, so `RamseySeriesSchema.to_norm` reaches level 3 of the run-start resolution and RAISES `ValueError: Cannot determine run start time`. All four tests in tests/test_packaged_fixtures.py pass `extra={"run_start_unix_s": 1.7e9}`; the docstring does not mention it. | Add the required `extra` to the docstring recipe (renaming the fixture to a DDMMYY_ prefix would instead buy the midnight UserWarning).

MINOR | src/quebra/_fixtures/__init__.py:41-44 | The `_RESOLVED` cache returns the path unconditionally on a hit, and the temp file's lifetime is `_STACK`'s. The leak itself is fine and documented (a wheel install is a directory, so `as_file` is a no-op), but the cache has no guard if anything ever closes the stack, and two threads racing the first call each `enter_context`, leaking one extraction. | `if name not in _RESOLVED or not _RESOLVED[name].is_file():` costs one stat on the error-free path.

MINOR | scripts/make_fixtures.py:26-27, 72 | "A rerun must reproduce the committed bytes" is TRUE here - regenerated into a clean tree, the CSV is byte-identical to the committed one (numpy 2.4.4, pandas 2.3.3, `cmp` clean) - but nothing enforces it and one dependency default breaks it off-platform: `DataFrame.to_csv` defaults `lineterminator` to `os.linesep`, so the same script on Windows writes CRLF and every byte after the first line differs. | Pass `lineterminator="\n"`, and add a test that regenerates into tmp_path and compares.

MINOR | scripts/make_fixtures.py:70-71 | "`%.9g` keeps ... the round-trip exact enough for a fixture" is correctly hedged - 9 significant digits is NOT a lossless float64 round-trip (17 are needed) - but the module docstring one screen above says the script exists "so the committed files are reproducible", which reads as value-exact. Measured, the loss is immaterial: 1 Hz quantisation against a 5 kHz noise sigma, and ~1e-13 s against a 4 us sigma. | Say "byte-reproducible; values are stored to 9 significant digits" once, at the top.

GREEN | scripts/make_fixtures.py:27-31 | Seed is a written-out constant, not derived from a clock or a hash, and every physical constant carries its unit suffix (`CADENCE_S`, `QUBIT_FREQUENCY_HZ`, `T2STAR_MEAN_S`, `t_s`, `drift_hz`, `t2star_error_s`). The unsuffixed CSV column names (`timestamp`, `frequency`, `T2star`) are imposed by `RamseySeriesSchema`, not chosen here.
GREEN | pyproject.toml:72-73 + .gitignore | The fixture really does ship: hatchling's wheel target takes all of `src/quebra`, and `git check-ignore` does not match `src/quebra/_fixtures/ramsey_synthetic.csv`, so the "inside the wheel" claim holds.

IMPORTANT | src/quebra/schemas/null.py:6-8 | Docstring: "`Dataset.schema=None` runs `RamseySeriesSchema`, which ten in-repo jobs depend on; `NullSchema` has to be asked for by name." The ten is real but it counts the wrong thing. There are exactly ten `schema=None` Datasets in `jobs/` (check_calibration 2, instrument_validation 5, ramsey_q1_100423 1, check_ledger_q1 1, independence_survey 1) and exactly ten `job.load_df()` calls, 1:1 - every one of those ten is a raw-DataFrame load, and `load_df` never consults `dataset.schema`. NOT ONE in-repo `job.load()` call site passes `schema=None`: all 8 pass `track912Schema` or `CalibrationLogSchema`. So changing the default tomorrow would alter the behaviour of zero jobs, not ten. What 7 job files actually depend on is `RamseySeriesSchema` running as the FALLTHROUGH after a `validate` schema (the track912 path), which is a different guarantee and the one worth stating. | Replace with the true statement: "seven in-repo jobs reach it as the fallthrough after `track912Schema.validate`; no in-repo caller relies on `schema=None`."

IMPORTANT | src/quebra/schemas/null.py:45-47 | Duplicate column labels silently produce a 2-D array. `frame[column]` on a duplicated label returns a DataFrame, so `.to_numpy()` is 2-D, and the loop then assigns the same key twice. Measured: `pd.DataFrame([[1,2,3],[4,5,6]], columns=["a","a","b"])` -> `norm["a"]` has shape `(2, 2)`, no error, no warning. Every step downstream expects 1-D under a Norm key. The `meta` collision two lines up proves the author was thinking about key collisions; this is the same class, unguarded. | Raise on `frame.columns.duplicated().any()` alongside the `meta` check.

MINOR | src/quebra/schemas/null.py:46 | `str(column)` can map two distinct labels onto one key. Measured: `pd.DataFrame({0: [1,2], "0": [3,4]})` -> a single key `"0"` holding `[3 4]`; the integer-labelled column is gone with no error. | Fold into the duplicate-label guard: check for collisions after `str()`, not before.

MINOR | src/quebra/schemas/null.py:37-44 | The `meta` it builds is `dataset_id` plus `**dataset.extra`, in that order, so a caller's `extra={"dataset_id": ...}` silently overrides the computed value. `RamseySeriesSchema` has the opposite precedence (`dict(extra)` then `.update({...})`, computed wins). Two schemas in one package disagreeing on who wins is a trap. Also absent here: `run_name`, `qubit`, `device`, `duration_h`, `n_points` - so a NullSchema Norm cannot say which qubit it came from even though the Dataset declared one. | Match RamseySeriesSchema's precedence, and carry qubit/device/n_points into meta.

MINOR | src/quebra/schemas/null.py:1 | "Every column of the file" is not what it returns. `_load_dataset` injects `qubit_id` and `device` into the frame before dispatch, so a NullSchema Norm carries up to two keys that are not columns of the file. | "Every column of the frame handed to it, including the qubit_id/device columns job.load injects."

MINOR | src/quebra/schemas/ramsey_series.py:103 | `stacklevel=2` was correct one frame out; the move added a frame and it was not re-counted. Measured: the DDMMYY warning is now attributed to `src/quebra/core/job.py:183` - the `return RamseySeriesSchema.to_norm(...)` line inside the library - where before it pointed at `_load_dataset`'s caller. Every such warning now blames the same library line regardless of which job triggered it. | `stacklevel=3`.

VERIFIED PASS | src/quebra/schemas/ramsey_series.py:50-145 vs 9057803:core/job.py:161-258 | The move is behaviour-preserving. Textual diff of the two bodies (indentation normalised) shows ONLY black re-wrapping of two `pd.to_numeric(...)` calls plus the relocated `timestamp`/`frequency` KeyError. The `[valid][order]` index chain is intact and identical on all three passthrough columns (`normalised chi-square`, `T2star`, `T2star error`), the meta dict construction order is unchanged, and the three-level run-start resolution including the no-`t_raw[0]`-fallback comment is verbatim. Live check on a synthetic `010423_demo.csv` gives `run_start_resolution='date_only_midnight'`, `run_start_unix_s=1680307200.0`, keys `t_rel_s/delta_hz/raw_frequency_hz/meta`.

VERIFIED PASS | src/quebra/core/job.py:140-190 | Dispatch is correct and total: `to_norm` before `validate` (same precedence as before), and the third branch is a raised `TypeError` naming both required methods and the doc - no silent passthrough. `import warnings`, `import numpy as np` and `check_unix_s` were removed from job.py and nothing there still references them. No import cycle: `transforms/lookup_prior.py` imports only numpy/pandas/typing.

VERIFIED PASS | src/quebra/core/job.py:150-155 | The conditional `qubit_id`/`device` injection is a real behaviour change but not a silent one. Old: `int(dataset.qubit)` on `qubit=None` raised `int() argument must be...`, naming neither field. New: the column is simply absent, and `BaseQubitSchema` declares `qubit_id` required, so the track912 path still fails - now with a pandera SchemaError that names the column. All 8 in-repo `job.load` datasets declare both fields, so nothing in-repo changes.

VERIFIED PASS | src/quebra/core/dataset.py:12-30 + loaders/registry.py | The `loader_kwargs`/`extra` split is complete. All four registered loaders audited: `_load_csv` does `pd.read_csv(path, **meta)`, `_load_hdf` reads `meta.get("key")`, `_load_yaml` and `_load_pickle` both `del meta`. Both `_load_dataset` and `_load_dataframe_raw` now pass `dict(dataset.loader_kwargs)`; `git grep` finds no remaining site that forwards `extra` to a reader. Every in-repo `extra=` is `{"run_name": ...}` on a `.pickle`, whose loader discarded it, so the split is a no-op for existing jobs - the docstring's account of why it was a defect is accurate.

IMPORTANT | tests/test_load_dataset_contract.py:128 | Same false claim as null.py, restated in a test docstring that presents it as the reason the test exists: "it is the behaviour ten in-repo call sites rely on". Zero in-repo `job.load` call sites pass `schema=None`. The ten `schema=None` Datasets are all `load_df` inputs and never touch a schema. | Same correction as the null.py entry.

IMPORTANT | tests/test_load_dataset_contract.py (whole file) | The contract file never tests the `validate` branch, which is the one every real job takes. All 8 in-repo `job.load` sites pass `track912Schema` (validate-shaped, 6 sites) or `CalibrationLogSchema` (to_norm-shaped, 2). The file tests `to_norm`, `schema=None`, and the neither-method TypeError - but the docstring's own middle claim, "a schema with `validate` cleans a frame and then the default normaliser runs on the result", has no test here and none anywhere else (`git grep` for `track912|_load_dataset|CalibrationLogSchema` over `tests/` returns only this file and `test_reference_model.py`). Also untested: a schema carrying BOTH methods, where the code silently prefers `to_norm`. | Add a two-line fake `validate` schema that renames a column, and assert the resulting Norm has `t_rel_s` - that pins the fallthrough the 6 track912 jobs live on.

MINOR | tests/test_load_dataset_contract.py:86-94 | `test_loader_kwargs_do_reach_the_file_reader` asserts only the NEGATIVE: that omitting `sep` raises KeyError. It never asserts anything about a run that supplies `loader_kwargs`, so the name overstates it, and a `KeyError` raised for any unrelated reason would satisfy it. The positive half exists but is in a different test. | `pytest.raises(KeyError, match="day")`, or assert the frame width with `sep` supplied.

MINOR | tests/test_load_dataset_contract.py:17-19 | "checked separately on four real 6D2S records (three `Norm` shapes, 7 keys each, arrays compared elementwise)". Three distinct shapes that all have exactly 7 keys is at best odd phrasing, and the harness is not in the repo, so no reader can re-derive the claim or re-run it after the next edit to `RamseySeriesSchema`. A number nobody can check is the class of claim AGENTS.md 4 exists to stop. | Either commit the comparison as a `@pytest.mark.real` test, or drop the parenthetical and say only that it was checked out of band.

VERIFIED PASS | tests/test_load_dataset_contract.py | 11 passed locally. No `filterwarnings = error` in `[tool.pytest.ini_options]`, so the two unasserted DDMMYY UserWarnings do not fail the gate - and their pytest output independently confirms the stacklevel finding above, printing `src/quebra/core/job.py:183` as the warning site.

IMPORTANT | tests/test_data_manifest.py:89-105 | `test_a_corrupted_record_would_be_caught` is declared "a positive control on the check above", but it never touches the check above, or any project code: it hashes a literal with `hashlib` and asserts the digest changes when a byte is appended. It is a test of `hashlib.sha256`, and it would still pass if `test_every_manifest_entry_matches_the_file_on_disk` were deleted, or if its assertion were inverted. The repo's own standard is that a positive control must FAIL when the thing it guards is broken. | Build a synthetic private root + manifest under `tmp_path`, run the SAME comparison (extract it into a helper), and assert it raises when the digest is wrong.

MINOR | scripts/make_data_manifest.py:67-68 | `_quote` escapes only `\` and `"`, but TOML basic strings forbid raw control characters. Measured through `tomllib`: a filename containing a newline or `\x07` produces a manifest that FAILS to parse (`Illegal character`), and since `_is_manifested` then degrades to False, the affected record silently stops being recognised as embargoed. Non-ASCII and tab are fine. | Escape the control range (`\n` -> `\\n`, else `\\u%04X`), or use `tomli_w`.

MINOR | scripts/make_data_manifest.py:99, 111 | `str(rel)` writes the OS separator, so a manifest generated on Windows records `2x2\030723_....pickle` while `paths._is_manifested` compares against `raw.as_posix()` - no record would ever match. `MANIFEST.write_text(...)` also uses the locale encoding, and TOML is defined as UTF-8. | `rel.as_posix()` and `write_text(..., encoding="utf-8")`.

MINOR | scripts/make_data_manifest.py:90 | `generated = date.today()` makes the output non-reproducible: rerunning on any later day produces a diff even when not one record changed, which is exactly what erodes trust in a "do not hand-edit" generated file. There is also no `--check` mode, so nothing detects a manifest that has drifted from `build()`. | Drop the date (git already dates the commit), or add `--check` comparing `build()` to the committed text.

MINOR | tests/test_data_manifest.py:7-9 + :26-39 | "These tests skip rather than fail when `data/real_private/` is absent" - the skips are correct where they exist (`_records`, `_present_files`, the `checked == 0` guard, and the conditional-expression `pytest.skip` at :80-85, which does work because `skip()` raises before the assignment). But the premise is wrong for a checkout: `MANIFEST.toml` is COMMITTED inside `data/real_private/`, so the directory always exists and the guards fire only for a wheel install. Verified in a synthetic checkout carrying only MANIFEST.toml: 4 passed, 1 skipped (`no manifested file is present to verify`) - nothing fails, which is the requirement. | Restate as "skip when the RECORDS are absent"; nothing to change in the code.

MINOR | tests/test_data_manifest.py | Nothing checks `manifest.total_bytes`, nothing checks the sha256 fields are 64 hex chars (`:46` only tests truthiness), and nothing checks the committed manifest still equals `scripts/make_data_manifest.build()`. The last is the one that would catch a hand-edit of a file whose header forbids hand-editing. | Add a skip-guarded regeneration test.

MINOR | src/quebra/core/paths.py:106 + quebra.toml:16 | SPEC 0003 set `data_root = "."`, so `dataset_root == repo_root()` for every in-checkout run and `candidates` holds the same resolved path twice. Every message now reads `Tried '/…/x.pickle' or '/…/x.pickle'` (measured, verbatim). `test_the_error_names_every_location_it_tried` passes trivially on the duplicate. | Dedupe candidates while preserving order before constructing the message.

MINOR | tests/test_data_unavailable.py | The suite covers absent, unparseable, manifested and unmanifested, but not the two shapes that actually break `_is_manifested`: a `[record]` table instead of `[[record]]` (raises AttributeError) and a path that matches only as a mid-component suffix (`/elsewhere/mine2x2/030723_2x2_qubit1.pickle` is reported True against the REAL manifest). Both are one-line additions to `fake_repo`. | Add the two cases alongside the fixes to paths.py.

GREEN | src/quebra/core/paths.py:110 + tests/test_data_unavailable.py | The intended distinction works end to end on the real manifest: a job's declared `data/real_private/6D2S/070423_6D2S_qubit1.pickle` is recognised as manifested, a sibling typo is not, `DataUnavailable` is still a `FileNotFoundError`, and the embargoed branch names the fixtures as the way forward. Nothing is swallowed - the raise is unconditional and only the WORDING depends on the manifest.
GREEN | tests/test_data_manifest.py:49-77 | The two checks that matter both run and are correctly directional: an unmanifested file present on disk fails, and a manifested file present on disk must match digest AND byte count, skipping per-entry only for records that are genuinely absent.

---

# SPEC 0003 Phase 1B review

SCOPE: scripts/promote_run.py, tests/test_promote_run.py, Makefile (promote target),
tests/test_data_root_resolution.py, tests/test_path_resolution.py,
tests/test_paths_dataset_fallback.py, tests/test_windows_not_interpolated.py,
quebra.toml, .gitignore, the nine jobs/ dataset repoints
COMMIT: b14bbc7

## Reviewed
- quebra.toml, .gitignore, Makefile, jobs/ (9 files), scripts/promote_run.py,
  tests/test_promote_run.py, and the four modified tests. Manifest empty.

## Findings

IMPORTANT | scripts/promote_run.py:29 | `_quote` escapes only `\` and `"`. A --note containing a newline or any control character writes a PROMOTED.toml that is not valid TOML; measured: note="line one\nline two" produces `note = "line one\nline two"` with a raw newline and `tomllib.loads` fails with "Illegal character '\n'". The script exits 0, so the committed audit record is silently unparseable. | Escape the TOML basic-string set (\b \t \n \f \r plus \uXXXX for other control chars), or refuse a note with control characters.

IMPORTANT | scripts/promote_run.py:69-71 | Destination is `published/{job_stem}_{identity[:6]}` with `mkdir(exist_ok=True)` and `shutil.copy2` on top. Measured: two DIFFERENT runs (identities abcdef1111 and abcdef2222) of the same job write into one directory - run A's `node_0.prov.json` is OVERWRITTEN, run A's other records survive alongside run B's, and PROMOTED.toml then says identity=abcdef2222, record_count=1 while six files sit in provenance/. A committed provenance record is destroyed and the manifest misdescribes what is there, silently. | Refuse to write into an existing destination unless the manifest's identity matches exactly, or use the full identity in the directory name.

IMPORTANT | scripts/promote_run.py:67 | `records[0][1].get("identity", "")[:6] or "unknown"` raises `TypeError: 'NoneType' object is not subscriptable` when the record has `"identity": null`, which `build_prov_record` (src/quebra/provenance.py:78) produces by default and which 142 of 400 sampled records under output/ actually carry. The `or "unknown"` branch is dead for the representation that occurs. Not reachable on a clean run in this checkout today (all clean runs have an identity), but reachable via --allow-dirty. | `(records[0][1].get("identity") or "")[:6] or "unknown"`.

MINOR | .gitignore:34-36 | The new comment says the two backups "are what this repository already holds (305M, 33M and 23M)". Measured today: `output_backup/` and `output_backup2/` do not exist, and `output/` is 309M, not 305M. A measured claim in a comment that is false as shipped. | State that the backup names are kept as a guard against re-creation, and drop or refresh the sizes.

MINOR | .gitignore:22-25 | The private-data comment still reads "These two trees do not exist yet - the dataset root is currently one level ABOVE the repository" and cites SPEC 0002. This commit is precisely what made both sentences false: data/real_private/ exists and data_root is ".". | Update the comment to SPEC 0003 and the present tense.

MINOR | tests/test_data_root_resolution.py:70 | Docstring of `test_a_relative_data_root_anchors_on_the_file_that_declares_it` still says "This repository's own `quebra.toml` says `data_root = \"..\"`". It now says ".". The test body is synthetic and unaffected, but the stated motivation is false. | Update the sentence.

MINOR | tests/test_data_root_resolution.py:111-126 and tests/test_path_resolution.py:38-50 | Both rewrites are honest re-expressions, not weakenings - but with declared `data_root = "."` the assertion collapses to `resolve_data_root() == repo_root()`, so a resolver that ignored quebra.toml entirely and returned the repo root would still pass. The docstring claim "through mechanism 3 rather than `__file__`" is not what the assertion demonstrates. Mechanism 3 itself is still properly covered by the synthetic-toml tests at lines 55-77. | Also assert the resolver reacts to a CHANGED declared value (write a temp quebra.toml with a different root), or soften the docstring claim to what is asserted.

MINOR | Makefile:54 | `python scripts/promote_run.py` uses bare `python`; the `arch` target three rules above uses `python3`. On a distro without the unversioned alias (this is Fedora) `make promote` fails with "python: command not found". | Use `python3`, or a `PYTHON ?= python3` variable.

MINOR | Makefile:54 | `"$(NOTE)"` is inside double quotes, so a note containing a backtick, `$(...)` or a `"` is interpreted by the shell rather than passed through; `$(PROMOTE_FLAGS)` is unquoted and reaches the shell whole. Self-inflicted only, but the guard rails above suggest care was intended. | Single-quote the expansions, or document PROMOTE_FLAGS as trusted input.

MINOR | scripts/promote_run.py:50 | `root: Path = Path(".")` makes the published tree cwd-relative, so running the script from a subdirectory silently creates a second `published/` there. The Makefile always runs from the repo root, so this only bites direct invocation. | Default to the repo root, or refuse when `root/PUBLISHED` is not under the checkout.

GREEN | scripts/promote_run.py:73-77 | `path.with_suffix("").with_suffix(".prov.md")` traced for `node_0.prov.json`, and for node names containing dots (`filter_0.5`, `a.1.5`, `a.prov`), a trailing dot, and multi-dot stems: it yields the matching `.prov.md` in every case and can never name a file outside `provenance/*.prov.md`. Copying is restricted to the `*.prov.json` glob plus that sibling, so no artifact can travel. | none

GREEN | scripts/promote_run.py:57 | `r.get("tree_clean", False)` (missing key = dirty) is the right default and matches src/quebra/core/runner.py:147, which uses the identical expression for reuse eligibility. 198 of 400 sampled legacy records lack the key and are correctly refused. | none

GREEN | tests/test_paths_dataset_fallback.py:33-46 | `test_without_the_fallback_that_path_would_not_exist` is STRONGER than what it replaced. The old form asserted only that a file was absent from this checkout's dataset root and never called the resolver; the new form drives an empty tmp_path root and asserts the resolver returns the repo-root path. Not tautological - it fails if the fallback is removed. | none

GREEN | jobs/ (9 files, 10 path lines) | All ten repoints resolve: every path exists under data/real_private/, and all 34 filenames in independence_survey.DATASET_FILES are present under data/real_private/6D2S/. `grep -rn "tool/datasets\|FOR ZENODO" jobs/ src/` returns 0. The survey's carve parameters were untouched, so tests/test_independence_survey.py's equality with the active T2* jobs is unaffected. Stale references remain only in AGENTS.md, docs/TIME_SEMANTICS.md and docs/WRITING_A_JOB.md - outside this slice. | none

GREEN | tests/ (all five) | `pytest` on the five files: 45 passed. No test was deleted and no assertion count dropped. | none

IMPORTANT | docs/WRITING_A_SCHEMA.md:77-78 | "It is also the schema `job.load` uses when `Dataset.schema` is `None`, which is why the jobs in `jobs/active/` that load Ramsey records do not name a schema at all." The second half is flatly false and a reader will act on it. EVERY Ramsey-record job in `jobs/active/` names a schema explicitly: `ramsey_2x2_q1_030723.py:24`, `ramsey_q1_100423.py:25`, `t2star_q1_070423.py:32`, `t2star_q1_100423.py:32` and `km_poster_6d2s.py` all pass `schema=track912Schema`; `mtbf_q1.py:32` and `mtbc_q6.py:35` pass `CalibrationLogSchema`. They reach `RamseySeriesSchema` as the fallthrough AFTER `track912Schema.validate`, which is the interesting mechanism and the one the doc should show. This is the third restatement of one wrong belief (see null.py and the test docstring). | Delete the clause, or replace with: "the six track912 jobs reach it as the fallthrough after their `validate` schema returns."

MINOR | docs/WRITING_A_SCHEMA.md:14-18 | The first code block does not run: it imports `NullSchema` and then uses `Dataset`, which is never imported. Measured - `exec` of the block verbatim raises `NameError: name 'Dataset' is not defined`. It is the first thing a new user copies. | Add `from quebra.core.dataset import Dataset`.

MINOR | docs/WRITING_A_SCHEMA.md:89 | "`duration_h=None,  # optional; run length, copied into meta`" is presented as a property of `Dataset`, but only `RamseySeriesSchema` and `CalibrationLogSchema` copy it; `NullSchema` - the schema the same document tells you to start with - drops it, along with `qubit` and `device`. A reader who follows the page top to bottom loses all three without notice. | Say "copied into meta by the shipped schemas; your own `to_norm` must do it itself."

VERIFIED PASS | docs/WRITING_A_SCHEMA.md | Every other runnable or measurable claim checks out. The `ReadingsSchema` block executes as written through `_load_dataset` and yields `{'t_rel_s': ..., 'value': ..., 'meta': {'dataset_id': 'r'}}` - `dataset.path.stem` works because `__post_init__` coerces `path` to `Path`. The `Dataset(...)` field list matches the dataclass in order and defaults, all seven fields. "Longest schema in the repo": ramsey_series 144 lines vs track912 83, null 48, calibration_log 46, base 11 - true. `t2star.make_inputs_from_norm` exists at `analyzers/t2star.py:46` (and in fidelity/mtbf/windows, so "and friends" holds). `calibration_log.py` does build `CalibrationEvent` objects. `pd.read_csv(path, **loader_kwargs)` and "an HDF5 key is loader_kwargs" both match `loaders/registry.py`. "Runs in about a second": measured 1.19 s for the 11 tests. The doc is 129 lines for a new public extension point - proportionate, not a dumped essay.

---

SCOPE: SPEC 0004 Phase 2 "Boundaries" - pyproject.toml, Makefile, .github/workflows/ci.yml,
src/quebra/core/job.py, src/quebra/core/runner.py, src/quebra/core/_artifact_guard.py (untracked),
src/quebra/panels/_artifact_guard.py (deleted), 7 repointed import files, spec/specboundaries04.md (untracked)
COMMIT: 1f58d37

## Manifest
- [ ] pyproject.toml (+71/-?)
- [ ] Makefile (+/-10)
- [ ] .github/workflows/ci.yml (+/-18)
- [ ] src/quebra/core/job.py (+17)
- [ ] src/quebra/core/runner.py (+47)
- [ ] src/quebra/core/_artifact_guard.py (untracked)
- [ ] src/quebra/panels/_artifact_guard.py (deleted)
- [ ] 7 repointed import files
- [ ] spec/specboundaries04.md (untracked)

## Reviewed
- [x] pyproject.toml
- [x] Makefile
- [x] .github/workflows/ci.yml
- [x] src/quebra/core/job.py
- [x] src/quebra/core/runner.py
- [x] src/quebra/core/_artifact_guard.py
- [x] src/quebra/panels/_artifact_guard.py (deleted)
- [x] 7 repointed import files

## Findings

IMPORTANT | pyproject.toml:178-181 | The comment justifying `follow_imports = "silent"` states a measurement that does not reproduce: "At the default (normal), mypy reports errors in followed modules too - plots/base.py, loaders/registry.py, schemas/ramsey_series.py, provenance.py". Measured today with a fresh cache and `--no-incremental` on three configs (shipped, `normal`, and a minimal `[tool.mypy]` with only python_version+normal+ignore_missing_imports): `mypy src/quebra/core` reports **0 errors / 9 source files in every case**, and `-v` shows **zero `followed=True` sources**. The real mechanism scoping errors to core is mypy's DEFAULT silencing of site-packages: `quebra` is an editable install (`_editable_impl_quebra.pth`) with no `py.typed`, so followed `quebra.*` modules are silenced regardless of `follow_imports`. Proof: adding `--no-silence-site-packages` turns the same run into "Found 628 errors in 36 files". So `follow_imports = "silent"` is a no-op for the shipped gate, and it is a global setting that will suppress real errors the moment `py.typed` is added or the scope widens. | Correct the comment to the measured behaviour, or drop `follow_imports` and state that site-packages silencing is what scopes the gate.

IMPORTANT | src/quebra/core/job.py:252-256 | The new comment says "a load node carrying a non-Dataset here would be a real defect worth failing on" - but the code `if not isinstance(ds, Dataset): continue` SKIPS it, silently omitting that dataset's content hash from `build_identity`'s `data` tuple. Before the change a non-Dataset reached `ds.path` and raised AttributeError. runner.py's `_dataset_of` (same condition, same phase) correctly RAISES. Comment contradicts code, and the direction is loud->silent in the identity computation. | Raise (reuse `_dataset_of`-style TypeError), or change the comment to say it is skipped.

MINOR | src/quebra/core/runner.py:44 | `_dataset_of(node: object)` uses `getattr(node, "kwargs", {})` - a silent `{}` fallback for a node with no `kwargs` - when `_DAGNode` is already imported in this module and the two call sites both pass a `_DAGNode`. | Type the parameter `_DAGNode` and use `node.kwargs.get("dataset")`.

MINOR | src/quebra/core/runner.py:216-221 | The comment says "they are always populated (run_job sets all three)" but FOUR fields are read here; `job_out_dir` is deliberately left un-`_require`d. If it were None the `job_out_dir not in d.parents` self-read guard silently disables (None is a legal `in` operand) - the one case where a None would produce a wrong reuse rather than a crash, and it is the one not guarded. | `_require` `job_out_dir` too, or say in the comment why it alone is exempt.

MINOR | pyproject.toml:189-190 | `[tool.deptry]`'s comment still opens "Scope the scan to the distributed package." With `deptry .` in the Makefile the scan is now repo-wide (107 files, was `src/quebra`); `known_first_party` declares first-party-ness, it does not scope. | Reword to "declare `quebra` first-party; the scan is repo-wide".

MINOR | pyproject.toml:129-133 | The layer comment says `cli`, `recipes` and `provenance` "are listed explicitly. Omitted, they sit outside every layer and are unconstrained in both directions" - but `quebra._fixtures` is a fourth top-level module and IS omitted, so it is unconstrained in exactly the way the comment warns against. (Harmless today: it imports nothing from quebra.) | Add `quebra._fixtures` at the bottom of `layers`, or say why it is exempt.

MINOR | .github/workflows/ci.yml:38-39 | "The `real` marker is what keeps those tests out; CI never sees the tree they need." Measured: `grep @pytest.mark` over tests/ yields only `parametrize` - ZERO tests carry `real`, `slow`, `heavy` or `r`, so `FAST` excludes nothing and the marker keeps nothing out. What actually keeps CI green without the private tree is per-test `pytest.skip` in `tests/test_data_manifest.py` (4 sites) and `tests/test_artifact_guard.py:305`. The claim is inherited from the deleted step, but it is now the sole stated safety argument for a step that runs the WHOLE suite. | Say that the guard is the per-test `pytest.skip` on an absent `data/real_private/`, and that the markers are currently unused.

VERIFIED | pyproject.toml [tool.importlinter] | `lint-imports` exits 0 (79 files, 159 deps, 1 contract kept). Exactly 10 `ignore_imports`. Ratchet is load-bearing both ways, measured on copies via `--config` (pyproject.toml left untouched, `git diff --stat` still 70 insertions): deleting `quebra.core.job -> quebra.loaders.registry` gives "Layered architecture BROKEN ... quebra.core.job -> quebra.loaders.registry (l.20)"; adding a fabricated `quebra.core.job -> quebra.loaders.fabricated` exits 1 with "No matches for ignored import". `unmatched_ignore_imports_alerting = "error"` is stated explicitly as claimed. All 10 entries match real imports (implied by the unmatched=error pass).

VERIFIED | src/quebra/core/job.py TYPE_CHECKING | `import quebra.core.job` loads no matplotlib/plotly/pylab module (measured over `sys.modules`). `_FigureSink` is still a dataclass with fields `plot_class, input, targets, name`; nothing calls `typing.get_type_hints` on it (`grep get_type_hints src/` is empty), so the string annotation `type[BasePlot]` is never resolved at runtime. `@dataclass(slots=True)` is unaffected - slots are built from `__annotations__` keys, not resolved types.

VERIFIED | src/quebra/core/runner.py `_require`/`_dataset_of` | Not a converted working path. `ResolutionContext(` is constructed at exactly ONE site repo-wide (runner.py:480, `grep` over src/ jobs/ tests/ scripts/), and it sets `subjobs_dir`, `job_out_dir`, `pool_root` and `dataset_root` unconditionally; `dataset_root = Path(data_root).resolve() if data_root else default_dataset_root()` can never be None. Both `load()` and `load_df()` always put a `Dataset` in `kwargs["dataset"]`. The `kwargs["dataset"]` -> `_dataset_of(node)` substitution at line 503 is equivalent: `kwargs` is still `node.kwargs` at that point.

VERIFIED | src/quebra/core/_artifact_guard.py | Move is byte-identical apart from the TYPE_CHECKING import and the cast. `from _typeshed import DataclassInstance` is under `if TYPE_CHECKING` and the cast target is a STRING literal, so nothing resolves at runtime. The cast is honest: all 10 concrete subclasses (SizeVsN, PowerVsDependence, ValidationCurve, ReadDependence, CheckLedger, SignalBand, ReliabilityBand, DistinguishBand, WithinCalibrationPanelData, AcrossCalibrationPanelData) are `dataclasses.is_dataclass() == True` and non-slots. No stale reference to `quebra.panels._artifact_guard` anywhere in src/, jobs/, tests/, docs/, AGENTS.md or spec/. 666 `.pkl` files under `output/`: 0 contain the string `_artifact_guard`, so the guard is never named in stored state (pickle records the concrete class's module).

VERIFIED | Makefile | `PKGROOT` has no remaining reference in the repo (only two prose mentions in spec/specboundaries04.md). `deptry .` exits 0, "Scanning 107 files... Success!".

VERIFIED | test suite | `pytest -m "not slow and not heavy and not real and not r"` = 378 passed, 2 skipped in 25s. The 4 tests that would flip without `data/real_private/` are all `pytest.skip` sites in `tests/test_data_manifest.py`, none ERROR: 378/2 -> 374/6 is arithmetically the reported figure. `testpaths = ["tests"]` so the pathless `pytest $(FAST)` in `make test` collects the same set the old CI line did.


--- SECOND PASS (fix layer + spec/specboundaries04.md), appended 2026-08-30 ---

## Reviewed (second pass)
- [x] spec/specboundaries04.md (untracked, 324 lines) - NEVER REVIEWED before this pass
- [x] pyproject.toml (fix layer: follow_imports dropped, deptry comment, _fixtures layer)
- [x] src/quebra/core/job.py (fix layer: isinstance -> raise)
- [x] src/quebra/core/runner.py (fix layer: _dataset_of typed, job_out_dir _require)
- [x] .github/workflows/ci.yml (fix layer: real-marker comment rewritten)

## Findings (second pass)

IMPORTANT | Makefile:22 (`arch:` -> `lint-imports`) | `.import_linter_cache/` can serve a STALE dependency graph across exactly the kind of file move R4.2.1 performed, and `make arch` reads it. Reproduced: `lint-imports --config <copy>` reports "Analyzed 80 files, 160 dependencies" and five violations naming `quebra.analyzers.* -> quebra.panels._artifact_guard` - a module deleted in this diff, at pre-move line numbers - while `--no-cache` on the identical config reports "79 files, 161 dependencies" and none of them. The stale entry (`.import_linter_cache/fa73d0*.data.json`, mtime 2026-08-30 00:32) is NEWER than the sources it misdescribes (`analyzers/reliability_band.py`, mtime 2026-08-29 15:31), so it survived both the edits and a rewrite; repeat cached runs keep returning it. The shipped contract is currently correct (79/159/KEPT with AND without the cache) and CI is unaffected (fresh checkout, and grimp self-ignores the dir via its own `.gitignore`), so this is a local-gate reliability defect, not a wrong result today. But the phase's whole claim is that direction is a build failure, and a build failure computed from a graph containing a deleted module is not one. | `arch: lint-imports --no-cache`, or add `rm -rf .import_linter_cache` to `make clean`.

MINOR | pyproject.toml:167 | Measured number in a shipped comment is wrong: "a mixin appearing in zero of 400 archived .pkl files". There are 666 `.pkl` files under `output/`, which is the figure `spec/specboundaries04.md` R4.2.1 uses. Zero is right (`find output -name '*.pkl' -print0 | xargs -0 grep -l _artifact_guard` = 0 of 666); 400 is not. Note `grep -rl --include='*.pkl' output` silently returns nothing here - it also reports 0 for `panels`, which 151 files do contain - so anyone re-measuring must use find|xargs. | s/400/666/.

MINOR | .github/workflows/ci.yml:44-47 | The rewritten comment fixes the marker claim but misattributes one of its two guards. `tests/test_artifact_guard.py:305` does NOT skip on an absent `data/real_private/`; it globs `output/*/subjobs_output/*/t2star_panel_data.pkl` and skips "no pre-split artifact present on this machine" - an absent `output/`, which is gitignored and therefore always absent in CI. The four `data/real_private/` skips are all in `tests/test_data_manifest.py`. Relatedly, "378 passed / 2 skipped becomes 374 passed / 6 skipped" is the local rename experiment, not what CI will print: a fresh runner also has no `output/` (one more skip) and no Rscript, which makes `tests/test_checks_c3_bridge.py:63` RUN rather than skip (it skips locally *because* Rscript is present). | Say `test_data_manifest.py` guards the private tree and `test_artifact_guard.py` guards an absent `output/`; label 374/6 as the local measurement, not the expected CI line.

MINOR | spec/specboundaries04.md:125 | "The move is a file move plus import-line edits in ten modules." Measured: seven modules carry the import (`analyzers/{calibration_summary,check_ledger,distinguish_band,reliability_band,signal_band}.py`, `panels/{_within_calibration_data,across_calibration}.py`), which is also what the first-pass reviewer scoped as "7 repointed import files". Nine files change in total if you count the added and deleted module. | s/ten modules/seven modules/.

MINOR | spec/specboundaries04.md:216-219 (R4.3.4) | The `.prov.json` re-run comparison is called "the only thing protecting the phase's behaviour-neutrality claim", yet R4.3 is marked **DONE 2026-08-29** and neither the requirement, the Acceptance block nor the Review outcome records that it was performed. Every other load-bearing claim in this document carries "MEASURED", "Verified before writing this requirement" or a re-measurable command. Behaviour neutrality does hold by inspection (`kwargs["dataset"]` -> `_dataset_of(node)` is value-identical, and `_require` only narrows a field that is unconditionally set), but the doc's own acceptance criterion is unrecorded. | Record the comparison, or downgrade R4.3.4 to "argued by inspection" and say so.

MINOR | spec/specboundaries04.md:275 | "The reviewer covered 8 of 9 scoped files before its session ended; the ninth was this spec, checked separately." Written by the document's own author before any independent check of it existed. As of this pass the statement is true; when committed it was an assertion about a review that had not happened. | Attribute the check ("second review pass, 2026-08-30") rather than asserting it in the passive.

MINOR | spec/specboundaries04.md:106, 185, 217 | Line citations are anchored to three different revisions. `job.py:35`/`job.py:385` and `runner.py:65`/`runner.py:365` are correct at `1f58d37`; `job.py:249` is correct in the CURRENT file; `job.py:254 narrows by isinstance` is now `job.py:261` (the fix layer inserted a 7-line comment) and `runner.py:362-366` is now `_emit_prov`'s signature, not the `fn.__name__` test. | State the revision each citation is against, or cite symbols instead of lines.

VERIFIED | spec/specboundaries04.md R4.0 "16 illegal imports in 6 rule violations" + the per-pair table | Reconstructed from `1f58d37` source, not from the tool: analyzers->panels 10 (5 `_artifact_guard` in calibration_summary/check_ledger/distinguish_band/reliability_band/signal_band, 4 `_within_calibration_compute` in fidelity/reliability_band/signal_band/t2star, 1 `within_calibration` in t2star), analyzers->plots 1 (fidelity->plots.theme), core->plots 2 (job->plots.base, runner->plots.targets), core->loaders 1, core->schemas 1, schemas->transforms 1. Sum 16 across 6 ordered pairs, and 16 - 1 (R4.1.4) - 5 (R4.2.1) = the 10 frozen entries. Table is exact.

VERIFIED | spec/specboundaries04.md:96 "813 modules against 648" | Exact, both figures. `python -c "import sys, quebra.core.job; print(len(sys.modules))"` = 648 with `matplotlib`/`plotly` both absent; adding `quebra.plots.base` = 813 with both present.

VERIFIED | spec/specboundaries04.md:123 "all 666 archived .pkl files, zero contain the string _artifact_guard" | 666 `.pkl` under `output/`; 0 match `_artifact_guard` (find|xargs, since `grep -r --include` is unreliable on this tree). Corroborated by what the archive DOES contain: 151 files name a `panels.*` module, and the recorded module paths are pre-src-layout (`analyzers.fidelity`, `panels.non_repairable`), so pickle really does record the concrete class's module and never the mixin.

VERIFIED | spec/specboundaries04.md:49-53 "242 errors in 51 files" at 1f58d37 vs "208 errors in 41 files" now | BOTH figures re-measured without checking out the old commit: `git archive HEAD src/quebra | tar -x -C <scratch>` then `mypy --config-file <bare> <scratch>/src/quebra` gives exactly "Found 242 errors in 51 files (checked 79 source files)"; `mypy src/quebra` on the working tree gives exactly "208 errors in 41 files". Presenting both is honest AND the stated reason for the drift is exact, checked file by file: under a bare config the working tree is 232/49, i.e. HEAD minus 7 errors in `core/runner.py`, 2 in `core/job.py` and 1 in `panels/_artifact_guard.py` (10 errors, 2 files - `job.py` stays on the list because its pandas `import-untyped` survives a bare config); the `[tool.mypy]` pandas override then removes exactly 24 `import-untyped` errors, one per module that imports pandas, of which 8 modules had no other error. 232-24 = 208, 49-8 = 41. R4.0's category breakdown for the 10 core errors is also exact: 4 arg-type, 2 attr-defined, 2 union-attr, 1 type-var, 1 import-untyped.

VERIFIED | spec/specboundaries04.md:128 "MEASURED: 10 entries" | Exactly 10 `ignore_imports`; `lint-imports` KEPT with and without `--no-cache` (79 files, 159 dependencies).

VERIFIED | spec/specboundaries04.md:55-65, the three tool behaviours | All three reproduced. (1) Moving `exclude_type_checking_imports` into the contract block makes the TYPE_CHECKING `job.py -> plots.base` import a violation again (161 deps vs 159), so it is top-level-only; a deliberately bogus `totally_bogus_option = "banana"` in the same block drew no complaint - unknown contract options ARE silently accepted. (2) `"quebra.core : quebra.provenance"` KEPT, `"quebra.core | quebra.provenance"` BROKEN naming `quebra.core.identity -> quebra.provenance (l.32)` - exactly the stated example. (3) With `unmatched_ignore_imports_alerting` deleted, a fabricated entry still fails: "No matches for ignored import ...", exit 1. R4.2.3's explicitness is belt-and-braces, as the text claims.

VERIFIED | spec/specboundaries04.md:141-146 (R4.2.4) | `panels/_within_calibration_data.py:36` does define `WithinCalibrationPanelData`, and 2 of the 666 archived pickles contain `quebra.panels._within_calibration_data`, `within_calibration` and the class name - so the move really would break them, unlike `_artifact_guard` (0 of 666). The four frozen `_within_calibration_compute` imports point at a module that appears in 0 pickles, but that module imports `_within_calibration_data` at its top level (and imports three analyzers back), so the argument transfers rather than failing. The rejected alternative reason is correctly labelled false: `Identity.code` is `hash_string(job_file.read_text())`, job file only.

VERIFIED | spec/specboundaries04.md "Not in this phase" | Honest on all four. No `jobs/validation|survey|comparison` exists (`jobs/` is active, bench, composite, reference, rscripts), so `acyclic_siblings` really would assert over absent packages. `grep @pytest.mark tests/` yields only `parametrize` and there is no `pytestmark` anywhere - zero tests carry slow/heavy/real/r. `tests/test_bench_isolation.py` is genuinely not subsumed: its `PIPELINE_PACKAGES` includes `"jobs"`, which `root_package = "quebra"` cannot see. The quote of quebraplan 2.1 is verbatim (quebraplan.md:140) and 2.3 does condition `acyclic_siblings` on the jobs/ split.

VERIFIED | spec/specboundaries04.md "Review outcome" vs the fixes actually present | Accurate. It reports 2 IMPORTANT + 5 MINOR, which is exactly what the first pass recorded, and all seven fixes are in the tree: `follow_imports` gone with a comment stating the measured no-op; `build_identity` raises `TypeError`; `_dataset_of(node: _DAGNode)` with no `getattr` fallback; `job_out_dir` inside `_require`; the deptry comment now says "the scan itself is repo-wide"; `quebra._fixtures` in `layers`; the CI marker comment rewritten (with the residual misattribution above). Its "verified and unchanged" list also still holds: one `ResolutionContext(` construction site (runner.py:483) setting all four fields.

VERIFIED | src/quebra/core/job.py:249-266 | The skip->raise change adds no failure mode any legitimate DAG can reach. A load node exists only via `Job.load`/`Job.load_df`, both annotated `dataset: Dataset` and both the only writers of `kwargs["dataset"]`; `Job.step` refuses any fn named `_load_dataset`/`_load_dataframe_raw`, so a user cannot forge one. `if ds is None: continue` is correctly restored AHEAD of the isinstance, and it agrees with `runner._dataset_load_nodes`, which filters `kwargs.get("dataset") is not None` - so the two modules treat a None-dataset load node identically (ignored) and a wrong-typed one identically (raise). Pre-change the same input raised `AttributeError` on `ds.path`, so this is loud->loud with a better message.

VERIFIED | src/quebra/core/runner.py:213-232, 443, 507 | `_require(context.job_out_dir, ...)` breaks no path that previously tolerated None. `ResolutionContext(` is constructed at exactly one site repo-wide (runner.py:483) and `job_out_dir` is already dereferenced two lines earlier (`job_out_dir.relative_to(out_dir)` at line 470), so it cannot be None there; `_locate_artifact` has exactly one caller, the `locate=` injection at that same site. `_dataset_of(node: _DAGNode)` is only ever called on members of `_dataset_load_nodes(job)` (line 443) or on a node whose id is in `resolved_datasets` (line 507), both of which already exclude a missing/None dataset, so its `TypeError` is unreachable from a valid DAG - it is a type narrowing, not a new gate.

VERIFIED | pyproject.toml layers - `quebra._fixtures` at the bottom | Correct for the stated purpose. `_fixtures` is a data package (`FIXTURES`, `fixture_path`, one CSV) importing nothing from `quebra`; bottom placement is what forbids that changing, which is exactly what the comment claims and what R4.1.3 exists to close. Worth knowing that it constrains one direction only: nothing now stops a production module importing packaged fixtures, since everything sits above it. Today only `tests/test_packaged_fixtures.py` imports it, and tests are outside `root_package` either way.

VERIFIED | gates | `make lint`, `make types`, `make arch`, `make deps` and `make check` each exit 0; `pytest -q` = 378 passed, 2 skipped (380 collected), matching P2 and "Done when". `mypy src/quebra/core` = "Success: no issues found in 9 source files". The `dev` extra carries ruff, mypy, import-linter>=2.13 and deptry, so CI's `make check` + `make deps` are installable; `requires-python = ">=3.11"` matches the 3.11/3.12 matrix and `[tool.mypy] python_version = "3.11"`.

---

# SPEC 0005 / Phase 3 review (identity closure)

SCOPE: src/quebra/core/closure.py, src/quebra/core/job.py (Job.code_hash only),
src/quebra/recipes.py (T2STAR_THRESHOLDS + configure_t2star_job only),
jobs/active/t2star_q1_070423.py, jobs/active/t2star_q1_100423.py,
src/quebra/analyzers/within_calibration_data.py,
src/quebra/analyzers/within_calibration_compute.py,
tests/test_identity_closure.py, tests/test_independence_survey.py
COMMIT: 1f58d37 (working tree, uncommitted)

## Manifest

- [ ] src/quebra/core/job.py (+75/-?)
- [ ] src/quebra/recipes.py (+113)
- [ ] jobs/active/t2star_q1_070423.py (-132)
- [ ] jobs/active/t2star_q1_100423.py (-132)
- [ ] src/quebra/analyzers/within_calibration_data.py (moved, 70)
- [ ] src/quebra/analyzers/within_calibration_compute.py (moved, 492)
- [ ] tests/test_identity_closure.py (new, 207)
- [ ] tests/test_independence_survey.py (+46/-?)

## Reviewed
- [x] src/quebra/core/closure.py

## Findings

---

## SPEC 0005 / Phase 3 review (identity + discovery slice)

SCOPE: src/quebra/core/discovery.py, src/quebra/core/job.py (Job.include only),
src/quebra/cli.py, JOB_ID/JOB_FAMILY in 12 jobs/active+jobs/composite files,
scripts/make_job_manifest.py, docs/JOBS.md, tests/test_job_discovery.py,
tests/test_job_manifest.py, spec/specidentity05.md
COMMIT: 1f58d37 (working tree, uncommitted)

### Manifest
- [ ] src/quebra/core/discovery.py (new, 125)
- [ ] src/quebra/core/job.py (Job.include only)
- [ ] src/quebra/cli.py (+33/-?)
- [ ] jobs/active/*.py x9 + jobs/composite/*.py x3 (JOB_ID/JOB_FAMILY)
- [ ] scripts/make_job_manifest.py (new, 134)
- [ ] docs/JOBS.md (new, 52)
- [ ] tests/test_job_discovery.py (new, 115)
- [ ] tests/test_job_manifest.py (new, 71)
- [ ] spec/specidentity05.md (new, 486)

### Reviewed

### Findings

### src/quebra/core/closure.py

CRITICAL | src/quebra/core/closure.py:1-25 (docstring) + effect on every job | The closure MISSES the pandera schema class that actually validates and cleans the loaded frame. `_load_dataset` (core/job.py:148) dispatches on `dataset.schema`, a class object carried on the `Dataset` at runtime; `Dataset` kwargs are skipped by `parameter_row` and grimp never sees the edge. Measured in an isolated copy of the tree: appending a line to `src/quebra/schemas/track912.py` left `t2star_q1_070423`'s code_hash byte-identical at `320f903b...` (whereas the same mutation to `analyzers/t2star.py` moved it). `track912Schema.validate` filters and coerces the frame, so this is code that changes the numbers while the identity asserts nothing changed - the exact defect the module was written to fix, one layer down. | Seed the traversal with `type(value).__module__` for every `Dataset` kwarg's `schema` (and, generally, for any kwarg whose value is a class or instance defined in `quebra.*`).

IMPORTANT | src/quebra/core/closure.py:6-11 ("keeps the radius narrow") and pyproject/spec claim R5.0.4 | The claim that the closure excludes render modules is FALSE as shipped. Measured closure of `t2star_q1_070423` (35 modules) contains `quebra.panels._within_calibration_render`, `quebra.plots.base`, `quebra.plots.theme`, `quebra.plots.allan_plot`, `quebra.plots.tlf_plot`, `quebra.plots.fidelity_helpers`. Mutation-verified: appending a line to `plots/theme.py` moves the digest (`b4192242...`), and to `panels/_within_calibration_render.py` moves it (`fe3e3996...`). Cause: the step functions now live in `quebra.recipes`, which function-locally imports `quebra.panels.within_calibration`, which imports the render module. So the recipe collapse (R5.3) re-created precisely the blast radius the docstring says was avoided. | Either drop the "narrow radius" claim from the docstring and the spec, or exclude `quebra.plots.*` / `*_render` from `_reachable`.

IMPORTANT | src/quebra/core/closure.py:120-127 `_reachable` seeded per-module | Seeding at MODULE granularity, not function granularity, means the whole of `quebra.recipes` seeds every t2star job. The measured closure therefore also contains `analyzers.allan`, `analyzers.tlf`, `transforms.interpolate`, `transforms.lookup_prior` - modules the T2* graph provably never executes (the docstring itself says it "adds NO interpolate node"). Editing the Ramsey/Allan path moves the T2* identity. Over-claiming is safer than under-claiming, but the docstring's "keeps the radius narrow" is not what the code does. | Say what it is (module-level, coarse), or seed from `fn.__code__.co_names`.

IMPORTANT | src/quebra/core/closure.py:178-184 `parameter_row` repr guard | A `set` kwarg breaks the cross-process determinism claim outright. Measured, three separate processes: `a(labels={'qq','y','ab','x','zz'})`, `a(labels={'y','zz','ab','x','qq'})`, `a(labels={'y','ab','zz','qq','x'})` - three different rows, three different digests, same job. String hashing is randomised per process and no ` at 0x` appears, so the guard never fires. | Render containers canonically: `sorted(repr(v) for v in value)` for `set`/`frozenset`, and reject anything else that is not a scalar/str/tuple/list/dict of scalars.

IMPORTANT | src/quebra/core/closure.py:178-184 | A `pathlib.Path` kwarg renders as `PosixPath('/home/sera/x/y.csv')` - an ABSOLUTE path with no ` at 0x` - so it enters the identity and breaks install-independence exactly as this module's own docstring (hazard 4) forbids. Realistic today: any future `rscript_path=` / sidecar-file kwarg. | Reject `Path` kwargs, or store them repo-relative.

IMPORTANT | src/quebra/core/closure.py:178-184 | A numpy array kwarg longer than the print threshold renders truncated: measured `array([   0,    1,    2, ..., 1997, 1998, 1999], shape=(2000,))`. Two arrays differing only in the elided middle produce an identical row and an identical identity - a silent collision, which is the failure `parameter_row` was added to close. | Hash array bytes rather than repr, or reject `ndarray`.

MINOR | src/quebra/core/closure.py:180 | `if "0x" in rendered and " at 0x" in rendered` - the first conjunct is subsumed by the second and can never change the outcome. | Drop it.

MINOR | src/quebra/core/closure.py:180 | The guard catches only the CPython default `<... at 0x...>` form. `functools.partial(len)` renders `functools.partial(<built-in function len>)` and passes; any custom `__repr__` embedding `id()`, a timestamp or a hostname without the literal " at 0x" also passes. The guard is a heuristic presented as a contract. | Allowlist accepted types instead of denylisting one repr shape.

GREEN | src/quebra/core/closure.py | Headline claim VERIFIED by mutation in an isolated tree copy: `analyzers/t2star.py` +1 line moves the digest 320f903b -> 2b8d8ab5 and restores exactly; `analyzers/kaplan_meier.py` (unreached) +1 line leaves it byte-identical. Repo left untouched (`git status` unchanged).
GREEN | src/quebra/loaders/registry.py | The runtime loader dispatch is NOT a closure hole: `_LOADER_REGISTRY` is populated by decorators inside `registry.py` itself, no `importlib`/entry points, and `registry` is in the measured closure.

### src/quebra/core/closure.py (second pass: `_graph()` and grimp)

CRITICAL | src/quebra/core/closure.py:51-59 `_graph()` | `grimp.build_graph(PACKAGE)` uses grimp's DEFAULT `cache_dir`, i.e. `.grimp_cache/` in the process CWD, and that cache is keyed on MTIME, not content. Measured in an isolated tree copy: append `import quebra.analyzers.kaplan_meier` to `analyzers/t2star.py` and restore its mtime (size changed 6894 -> 6940) and the closure stays at 35 modules with kaplan_meier ABSENT; the identical run with `cache_dir=None` gives 36 modules with kaplan_meier PRESENT. So any workflow that preserves mtime while changing content (`cp -p`, `rsync --times`, `tar -x`, a container layer restore, some `git` operations) yields a wrong import graph and therefore an identity that silently stops covering a module. That is the exact failure this module exists to prevent, reintroduced by a default argument. | `grimp.build_graph(PACKAGE, cache_dir=None)`.

IMPORTANT | src/quebra/core/closure.py:51-59 | Computing an identity now has a filesystem SIDE EFFECT and a hard failure mode. Measured: running `Job.code_hash()` from an empty scratch directory creates `.grimp_cache/` there; running it from a read-only CWD dies with `PermissionError: Permission denied (os error 13)` out of `grimp/adaptors/caching.py:229`, with no quebra-level message. Identity computation writing to CWD is I/O in the graph runtime, and read-only CWDs are normal in CI containers. | Same fix: `cache_dir=None`.

IMPORTANT | src/quebra/core/closure.py:61-75 `_quebra_imports_of_file` | `except (OSError, SyntaxError): return set()` SWALLOWS the error and returns an EMPTY seed set, which makes the closure empty and silently restores the pre-SPEC-0005 behaviour (job file text only). Project rule is errors raised, not swallowed - and this module's own `_seed_modules` raises for the analogous "cannot resolve a module" case two functions below, so the file is internally inconsistent about it. | Re-raise with the job path named.

IMPORTANT | src/quebra/core/closure.py:122-124 `_reachable` | `stack = [s for s in seeds if s in known]` silently discards any seed grimp does not know. It is needed because `_quebra_imports_of_file` deliberately emits attribute names (`quebra.recipes.RAMSEY_CONFIG`) alongside real modules, but the filter cannot tell "that was an attribute" from "that module is missing from the graph". A stale cache or a package-layout change therefore drops coverage with no signal. | Resolve `from X import y` candidates against `graph.modules` explicitly and raise if NEITHER the module nor the attribute parent is known.

MEASURED | cost | `import grimp` 27 ms, `build_graph("quebra")` 3 ms warm / 81 modules; per-job `code_hash()` 1.4-4 ms for simple jobs. Not a hot path. `instrument_validation` 545 ms and `independence_survey` 141 ms are dominated by job import, not by the closure. Cost is not a finding.

VERIFIED | determinism | Claim 3 parts 1 and 3 hold. Twelve jobs hashed all-in-one-process twice: byte-identical. The same twelve hashed one-per-subprocess: byte-identical to the all-in-one-process run, so `run --all` does not contaminate. Running from a copied tree at a different absolute path reproduced `320f903b...` exactly, which is the install-independence claim for the tree location.

IMPORTANT | src/quebra/core/discovery.py:72-75 | `read_declaration` swallows `OSError` and `SyntaxError` and returns None, so a job file that does not parse becomes INVISIBLE to discovery. Measured: a dir with `good.py` and a syntactically broken `broken.py` both declaring `JOB_ID = "dup"` returns `{'dup': good.py}` with NO `DuplicateJobId`. Consequence for `run --all`: the old code globbed `jobs/active/*.py` and let the import raise into the `failures` list (exit 1); the new code drops the file before anything can raise, so a broken job is skipped with exit 0. Loudness regression. | Re-raise (or collect) parse/read errors instead of returning None; at minimum let `discover` fail on a `*.py` under `jobs/active|composite` that will not parse.
IMPORTANT | src/quebra/core/discovery.py:87-88 + tests/test_job_discovery.py:100 | A job file that simply forgets `JOB_ID` is silently excluded from `run --all` forever, and no test catches it: `test_every_in_repo_job_declares_an_id_and_a_family` iterates over what discovery already FOUND, so a file with no `JOB_ID` cannot fail it. | Add a test asserting every non-`__init__` `.py` under `jobs/active/` and `jobs/composite/` yields a declaration.
MINOR | src/quebra/core/discovery.py:50-67 | `_string_constant` silently ignores a non-literal `JOB_ID` (`JOB_ID = "pre"+"fix"` -> job vanishes). Documented and tested as intended, but the failure mode is invisibility, not a message. | Raise when the name is assigned but not a plain string literal.
GREEN | src/quebra/core/discovery.py | No-import rule holds: `ast.parse` only, no `import`/`importlib`/`exec` anywhere in the module, and a file whose body is `raise RuntimeError(...)` is discovered fine (measured). | none
GREEN | src/quebra/core/discovery.py:97-104 | Duplicate `JOB_ID` raises `DuplicateJobId` and the message names BOTH paths (measured on two well-formed files). | none

IMPORTANT | src/quebra/cli.py:105-106 | `--family` help says "see `quebra jobs`". There is no `jobs` subcommand: `quebra --help` lists `{run,inspect,schema-wizard}`. Shipped help text pointing at a command that does not exist. | Point at `docs/JOBS.md`.
IMPORTANT | src/quebra/cli.py:136-154 | `--all` silently changed meaning. Old: `jobs/active/*.py` = 9 jobs. New (measured with `run_job` mocked): 12 jobs, newly including the three composites `check_ledger_q1`, `compare_t2star_0704_vs_1004`, `independence_survey`. Composites re-run their own sub-jobs, so `run --all` now executes the two t2star jobs three times over. Neither the spec's R5.4 section nor the CLI comment mentions the widening. | Either keep composites out of a bare `--all` (e.g. require `--family composite`) or state the widening in the spec and the help.
IMPORTANT | src/quebra/cli.py:103-107,136 | `--family` is accepted and silently ignored when `--all` is absent. Measured: `quebra run jobs/active/mtbf_q1.py --family totally-bogus` runs the job and exits 0. A mistyped selector silently does the wrong thing. | `parser.error` when `--family` is given without `--all`.
MINOR | src/quebra/cli.py:142-143 | `declared = discovery.discover(jobs_root) if jobs_root.is_dir() else {}` -> from a cwd with no `jobs/`, `quebra run --all` prints nothing, runs nothing and exits 0 (measured, RC=0). Same shape as the old dead glob, but `--all` is now the only selection route. | Error when `--all` selects zero jobs.
MINOR | src/quebra/cli.py:142 vs src/quebra/core/job.py:324 | Two different roots for the same tree: the CLI uses `Path.cwd()/"jobs"`, `Job.include` uses `repo_root()/"jobs"`. From a subdirectory the CLI finds nothing while `include` still resolves. | Use `repo_root()` in both.
GREEN | src/quebra/cli.py | `--include-archived` removal is clean: `jobs/archived/` does not exist, and no reference to the flag survives in tests, Makefile, CI, README, AGENTS.md or docs (grepped). `--all --family t2star` -> exactly the 2 t2star jobs; `--family composite` -> exactly the 3 composites; an unknown family calls `parser.error` and lists the declared families. | none

GREEN | src/quebra/core/job.py:312-333 | ID-first/path-second dispatch verified end to end: with `t2star_q1_070423.py` COPIED to `jobs/some/deep/place/renamed_file.py` in a throwaway project root, `quebra inspect jobs/cmp.py` still wires `t2star_q1_070423:t2star_panel_data`. (Done in a sandbox; the tracked tree was never moved.) Edge cases all loud: `''` -> ValueError listing declared IDs; `bare.py` -> path branch, resolves; `jobs\a\b.py` -> FileNotFoundError naming the resolved path; `/abs/x.py` -> FileNotFoundError; an ID containing a dot (`sub.dotted`) resolves as an ID. | none
MINOR | src/quebra/core/job.py:322-326 | When `jobs/` is absent, an argument that was meant as an ID falls through to `resolve_repo_path` and reports a path-not-found error naming a path the caller never wrote. | Mention "no jobs/ directory found" in that branch.

IMPORTANT | scripts/make_job_manifest.py:26-37 | `INTERESTING` omits parameters that decide what a run means, including two seeds: `VALIDATION_SEED=20260811`/`VALIDATION_REPLICATES=4000`/`VALIDATION_N` (check_calibration), `ASYMPTOTIC_SIZE_SEED=777`/`N_PERM_IN_TIE_STUDY=999`/`DIVERGENCE_THRESHOLD=0.02` (instrument_validation), `TIE_CUTOFF_DISTINCT=5` (check_ledger_q1, independence_survey), `C3_N_NULL_SIM=200`. Measured: editing `VALIDATION_SEED` leaves `build()` byte-identical, so the manifest and its staleness gate are both blind to it. The docstring claims the list is "the ones that decide what a run MEANS". | Surface every UPPER_CASE module-level literal, or subtract an explicit `BORING` set so a new constant is included by default.
IMPORTANT | tests/test_job_manifest.py:60-70 | `test_regenerating_is_idempotent` SHELLS OUT to the generator, whose `main()` writes `docs/JOBS.md`. Measured: `pytest tests/test_job_manifest.py` changes the file's mtime. The staleness gate therefore self-heals - the run that reports "JOBS.md is stale" also rewrites it, so the next run is green with nothing fixed. | Have the test call `build()` twice, or run the generator against a tmp copy; never let the suite write a committed artifact.
IMPORTANT | tests/test_job_manifest.py:40-52 | `test_the_generator_does_not_import_the_jobs` does not test the generator. It writes an explosive job into `tmp_path`, then asserts `discover(tmp_path)` (the wrong function) and `assert gen.build` (truthiness of a function object). `gen.build()` is never called, and cannot be pointed at `tmp_path` because `JOBS` is a module constant. A test asserting coverage it does not have. | Make `build(jobs_root=JOBS)` take the root and call it on `tmp_path`.
MINOR | scripts/make_job_manifest.py:53-56 | `except (ValueError, TypeError): pass` drops a non-literal INTERESTING constant from the row with no marker, so the reader sees an absent parameter rather than "computed". | Render `<computed>`.
MINOR | tests/test_job_manifest.py:20 | Docstring "Fails when a job changed without regeneration" overclaims: measured, appending a module-level constant or changing `VALIDATION_SEED` leaves `build()` identical. It fails only when JOB_ID/JOB_FAMILY/Dataset literals/include literals/an INTERESTING constant change. | Say what it actually covers.
MINOR | tests/test_job_discovery.py:24-26 | `pytest.importorskip("pathlib")` to obtain `Path` - stdlib, cannot be missing; if it ever were, the whole module would skip silently. | `from pathlib import Path`.

IMPORTANT | jobs/composite/independence_survey.py:60 | `JOB_FAMILY = "composite"` on a job that contains NO `job.include` at all - its own docstring (line 14) says "this one loads its datasets directly rather than through `job.include`". The family was copied from the directory, which is precisely what R5.4 claims to have stopped. | Give it a subject family (e.g. `validation` or `independence`), or drop "composite" as a family.
MINOR | jobs/ (all 12) | "composite" is a mechanism, not a subject; the other five families (t2star, ramsey, survival, interval, validation) are subjects. A reader cannot predict where a new job goes - `check_ledger_q1` is as much "validation" as `check_calibration` is. | Either make family purely subject-matter and expose "has includes" from the manifest, or document the rule in docs/JOBS.md.
MINOR | jobs/active/ramsey_q1_100423.py, jobs/active/t2star_q1_100423.py | Both declare `PREFIX = 'q1_13h_1004_dataset'` (visible in docs/JOBS.md rows). Distinct JOB_IDs now, but the shared PREFIX is what names node/artifact labels. | Confirm the two cannot collide in `output/`.

### Reviewed
- [x] src/quebra/core/discovery.py
- [x] src/quebra/core/job.py (Job.include)
- [x] src/quebra/cli.py
- [x] jobs/active + jobs/composite JOB_ID/JOB_FAMILY (12 files)
- [x] scripts/make_job_manifest.py
- [x] docs/JOBS.md
- [x] tests/test_job_discovery.py
- [x] tests/test_job_manifest.py

CRITICAL | src/quebra/cli.py:136-154 vs docs/WRITING_A_JOB.md:229 and AGENTS.md:200 | The committed docs state the invariant this change breaks. WRITING_A_JOB.md:229: "`quebra run --all` sweeps `jobs/active` and would otherwise re-run every sub-job" - given as the REASON composites live in `jobs/composite/`. AGENTS.md:200: "`jobs/composite/` ... NOT swept by `run --all`". Measured, `run --all` now selects all 12 including the composites, so it re-runs every sub-job and additionally launches `independence_survey`, which this spec's own budget section calls "~6.5 h wall-clock at INCLUDE_C3=True over 415 cells". Neither doc was updated and the spec never mentions the widening. | Exclude family `composite` from a bare `--all`, or update both docs and say so in R5.4.
IMPORTANT | docs/WRITING_A_JOB.md:210-232 | The job-authoring guide never mentions `JOB_ID` or `JOB_FAMILY`, and still teaches `job.include("other_job.py")`. A job written by following the guide declares no `JOB_ID` and is therefore silently invisible to `run --all` forever. | Add the two constants to the guide and switch the include example to an ID.
IMPORTANT | docs/WRITING_A_JOB.md:173,185 | R5.1.6 required a `docs/` sweep and it did not happen. Line 173 still says the identity is "a content hash of the job file's source, the datasets it loaded, and the identities of any sub-jobs" - R5.1 added the step-source closure and the parameter row. Line 185's example run directory is `output/t2star_q1_070423_43e8d4_...`, the exact stale prefix R5.1.6 named; the spec's own R5.3 outcome says that digest is now `ffba1a`. | Sweep both.
MINOR | src/quebra/cli.py:145 | `quebra run --all --family ""` - empty string is falsy, so the guard is skipped and all 12 jobs run instead of erroring. | Test `args.family is not None`.

IMPORTANT | spec/specidentity05.md:346 | "the two **133-line** files were byte-identical". Measured at HEAD: both `jobs/active/t2star_q1_0{7,10}0423.py` are **128** lines, which is also what this same spec says at R5.0.2 line ~80 ("128 lines each"). The document contradicts itself and the artifact. | Say 128.
IMPORTANT | spec/specidentity05.md:347 | "They are 34 lines each now — **198 lines removed**". Measured: 34 lines each is correct; the removal count is not. `git diff --numstat` gives 113 removed + 19 added per file = 226 removed / 38 added, net 188 fewer lines in the two job files, against 113 lines ADDED to `recipes.py` (net -75 repo-wide). 198 corresponds to no measurement. | Quote 188 (net, two files) or 226/38 and the +113 in recipes.
IMPORTANT | spec/specidentity05.md:150-152 (R5.0.4 acceptance) vs pyproject.toml:163,166-177 | Acceptance says "no entry names a phase that will not remove it". All five surviving `ignore_imports` entries name SPEC 0005, and the spec's own "Not in this phase" section defers them. The comment above the list also still reads "The frozen baseline: **10** imports" over a list of 5, and an orphaned comment block about `panels/_within_calibration_data.py` survives with a dangling "# SPEC 0005." and no entry. Measured numbers in a shipped comment that are false. | Fix the comment and retarget the five entries at the spec that will remove them.
MINOR | spec/specidentity05.md:296-299 (R5.2 acceptance) / :338-348 | R5.2 has no `**Outcome.**` section; the paragraph that opens R5.3's Outcome ("`core/discovery.py` reads `JOB_ID` and `JOB_FAMILY` statically ... the acceptance holds") is R5.2's outcome pasted under R5.3. Same shift puts R5.4's outcome inside R5.5. A reader checking R5.2 or R5.4 finds no record. | Move each Outcome under its own requirement.
MINOR | spec/specidentity05.md:290-294 (R5.2.4) vs src/quebra/core/job.py:322-333 | R5.2.4 says "No deprecation window ... If an external caller ever appears, add the compatibility path then, against a real case." The shipped `include` keeps the path branch permanently and the `FileNotFoundError` text advertises it ("a repo-root-relative path still works"). Reasonable engineering, but the spec asserts the opposite of what was built. | Reword R5.2.4 to state that the path form is retained as a fallback, and why.
MINOR | spec/specidentity05.md:411 (R5.5 acceptance) | "its test fails when a job changes without regeneration" is false in general. Measured against a copy of `jobs/`: changing `ALPHA` moves `build()`; changing `VALIDATION_SEED` or adding a module-level constant does not. The gate covers JOB_ID, JOB_FAMILY, `Dataset(path=...)` literals, `include` literals and the ten `INTERESTING` names only. | State the covered surface.
MINOR | spec/specidentity05.md:430-435 ("Done when") | "The collect count is recorded with its delta explained." R5.1's outcome records 380 -> 393; the tree as shipped collects **408** (`406 passed, 2 skipped`) and no line of the spec records that number or its delta. | Record 408 and the +15 from R5.4/R5.5's tests.
MINOR | spec/specidentity05.md:382-384 (R5.4 acceptance) | "A job with no `JOB_FAMILY` is reported, not silently skipped" is satisfied by a test over the repo, not by the tool: at runtime a family-less job is silently omitted from every `--family` selection and nothing prints. | Have `--family` (or a `quebra jobs` listing) name the uncategorised jobs.
GREEN | spec/specidentity05.md:75-79 (R5.0.2) | "20 non-`__init__` files — 9 active, 3 composite, 8 bench" - counted, exact. The "byte-identical once the date and the run duration were normalised" claim also holds: normalising `070423`/`100423`, `27h`/`13h`, `duration_h`, `0704`/`1004` makes the two HEAD files `diff`-clean. | none
GREEN | spec/specidentity05.md:140-144 (R5.0.4 outcome) | "`analyzers.t2star`'s transitive closure went from 14 quebra modules with 4 render to **9 with 0**" - re-measured with `grimp.build_graph("quebra").find_upstream_modules("quebra.analyzers.t2star")`: 9 modules, zero `plots.*` or `*_render`. | none
GREEN | spec/specidentity05.md:296 | "All 4 in-repo `include` sites use `JOB_ID`" - grepped: exactly 4, in `check_ledger_q1.py` and `compare_t2star_0704_vs_1004.py`, all IDs. The remaining `.include(` calls are tests passing tmp_path absolutes. | none
NOT VERIFIED | spec/specidentity05.md:326-335, 349-353 | The digests `5c35194470316077`, `43e8d4 -> ffba1a`, the "byte-identical `q1_27h_0704_dataset_windows.pkl`" and the "2 bytes" `t2star_panel_data.pkl` delta all require running the T2* jobs over `data/real_private/`, which this review is forbidden to touch. Unchecked, and they are the spec's load-bearing numbers. | Another reviewer with data access must confirm.

- [x] spec/specidentity05.md

---

SCOPE: src/quebra/recipes.py (T2* family section only), jobs/active/t2star_q1_070423.py,
jobs/active/t2star_q1_100423.py, tests/test_independence_survey.py (diff vs HEAD).
COMMIT: 1f58d37

## Manifest
- [ ] src/quebra/recipes.py (T2* section)
- [ ] jobs/active/t2star_q1_070423.py (+19/-113)
- [ ] jobs/active/t2star_q1_100423.py (+19/-113)
- [ ] tests/test_independence_survey.py (+34/-12)

## Reviewed

## Findings

### Verified by measurement (behaviour preservation, claim 1)

GREEN | recipes.py + both jobs | DAG is byte-for-byte the same shape as HEAD. Built HEAD's two job files (from `git show`) and the new two in-process and dumped `dag`: identical node ids in identical order (`load`, `t2star_filter`, `t2star_final_filter_stage`, `t2star`, `windows`, `t2star_panel_data`), identical `inputs` edges (`t2star_panel_data <- [t2star, windows]`), identical sinks in identical order (`_MaterializeSink q1_..._windows` then `_FigureSink q1_..._t2star` with `targets=['static','academic']`, `WithinCalibrationPanel`). `quebra inspect` agrees. | none

GREEN | claim 2, provenance visibility | All five land in `node.kwargs`, measured not read: `windows` -> `{'gap_mult': 10.0, 'k': 1.0, 'use_uncertainty': True}`; `t2star_panel_data` -> `{'shape_min_reads': 5, 'use_uncertainty': True, 'xi_seed': 20260813}`. Identical dicts at HEAD. `runner._label` builds from `node.kwargs`, so the label is unchanged. | none

GREEN | claim 3, threshold ladder | `[(f"{k} µs", k / 1e6, True) for k in range(1, 11)]` reproduces the ten old literals EXACTLY, compared by IEEE-754 bit pattern (`struct.pack('>d')`) on all ten rungs, labels included. `k * 1e-6` would have differed at k=5 (`4.9999999999999996e-06`) and k=10 (`9.999999999999999e-06`); the shipped form does not. | none

GREEN | claim 4, the two jobs still differ | They differ in exactly path / PREFIX / duration_h and nothing else; the two DAGs differ only in those three places, and `Job(name)` differs, so identities differ. | none

### Findings

IMPORTANT | tests/test_independence_survey.py:88-105 | The rewritten `_effective_t2star_carve` is NOT strictly stronger; it is weaker on one axis and the weakening is silent. It pre-seeds `effective` with the recipe defaults and then overwrites from the job file inside `try: ... except (ValueError, SyntaxError): pass`. A job that overrides with anything non-literal keeps the swallowed exception AND keeps the recipe default, so the test compares the survey against a value the job does not use. MEASURED: injecting `gap_mult=MY_GAP` into the job file yields `{'GAP_MULT': 10.0, ...}` and the test PASSES while the job carves with `MY_GAP`; `k=2.0 * 1.5` likewise yields `K: 1.0` and passes. At HEAD the same mutation left `found` without the key and `assert name in found` fired loudly. That is precisely the "survey and panel describe different windows" defect the module docstring says this file exists to catch. | Re-raise instead of `pass`: on `ast.literal_eval` failure, fail the test naming the kwarg, e.g. `pytest.fail(f"{kw.arg} is overridden with a non-literal; the control cannot read it")`.

IMPORTANT | src/quebra/recipes.py:297-302 vs src/quebra/analyzers/t2star.py:117-128 | The new comment claims "One ladder, one place". False as shipped: `analyzers/t2star.py` still holds `T2STAR_DEFAULT_LADDER`, the same ten rungs as explicit literals. The count went from 3 copies to 2, not to 1, and nothing pins the two equal - no test compares `recipes.T2STAR_THRESHOLDS` to `T2STAR_DEFAULT_LADDER` (grepped). Project rule: a claim in a shipped comment must be true as shipped. | Either import the analyzer ladder, or say "two ladders" and add an equality assertion.

IMPORTANT | src/quebra/recipes.py:319, jobs/active/t2star_q1_070423.py:5, jobs/active/t2star_q1_100423.py:5 | "two 133-line job files" / "these two files were 133 lines each". MEASURED: `git show HEAD:jobs/active/t2star_q1_070423.py | wc -l` = 128, same for 100423. The same wrong number is already logged against spec/specidentity05.md:346 and has now been copied into three shipped source files. | Say 128.

IMPORTANT | src/quebra/recipes.py:311, 338 | `thresholds` is the one parameter of `configure_t2star_job` that does NOT reach provenance: `ladder` is a closure capture of `_windows_run`/`_t2star_panel_data`, never a step kwarg. The function's own docstring states the rule it breaks - "an argument left to its default would be invisible to provenance" - and applies it to the other five. The ladder is the most meaning-bearing parameter in the graph (it defines every rung of the panel); a job that passes `thresholds=` gets no record of it on the node label. | Forward `thresholds=ladder` as a step kwarg on `windows` and `t2star_panel_data`, as the other five are.

MINOR | src/quebra/recipes.py:338 | `thresholds=[]` is accepted silently: `[] is not None`, so the empty ladder is used, the panel scores zero rungs and nothing raises. Repo rule is errors raised, not silent degradation. | `if thresholds is not None and not thresholds: raise ValueError(...)`.

MINOR | src/quebra/recipes.py:311 vs :271 | `gap_mult: float = 10.0` hardcodes a number that is already imported into this module as `DEFAULT_GAP_MULT` (line 16) and used by `configure_ramsey_job` at line 271. Equal today (`windows.DEFAULT_GAP_MULT = 10.0`, checked); if the constant moves, ramsey follows and t2star silently does not. | `gap_mult: float = DEFAULT_GAP_MULT`.

MINOR | src/quebra/recipes.py:334-336 | The three function-body imports are cargo-cult, not cycle-breaking. MEASURED: `import quebra.analyzers.t2star` does not pull in `quebra.recipes` (no cycle); `WithinCalibrationPanel` is ALREADY imported at module scope (recipes.py:19-22), so that line is dead; and `recipes` sits above `panels`/`analyzers` in the import-linter layer list (pyproject.toml:147-152), so a module-level import is permitted. Function-level imports are also visible to grimp, so they buy nothing from the contract either. | Move all three to module scope alongside the existing `windows`/`WithinCalibrationPanel` imports.

MINOR | src/quebra/recipes.py:305-311 vs :204-215 | Inconsistent with `configure_ramsey_job` directly above: ramsey takes `dataset` POSITIONALLY and annotates `job: Job` bare; t2star makes `dataset` keyword-only and quotes `"Job"` / `"Dataset"` even though both are module-level imports and `from __future__ import annotations` is in force. Two sibling recipes with two conventions. | Match the neighbour.

MINOR | src/quebra/recipes.py:300 | `T2STAR_THRESHOLDS` holds SI seconds and carries no unit suffix, against the project's unit-suffix rule; the comment carries the unit instead. (Same defect in `T2STAR_DEFAULT_LADDER`, pre-existing.) | `T2STAR_THRESHOLDS_S`.

MINOR | tests/test_independence_survey.py:71 | `T2STAR_JOB = _module_constants(...)` is now dead - its only consumer was the `_T2STAR_THRESHOLDS` lookup this diff removed (grepped: one occurrence, the assignment). | Delete it.

MINOR | tests/test_independence_survey.py:166-176 | The ladder control lost its anchor. At HEAD it compared the survey's `k / 1e6` comprehension against the job's ten explicit LITERALS - a genuine cross-form check. It now compares the survey comprehension against `recipes.T2STAR_THRESHOLDS`, which is the character-for-character identical expression, so the surviving comment "`k * 1e-6` fails this at k = 5 and k = 10 ... `k / 1e6` reproduces the literals exactly" describes a comparison the test no longer performs. The only remaining explicit-literal ladder, `T2STAR_DEFAULT_LADDER`, is now pinned by nothing. | Compare the survey ladder against `analyzers.t2star.T2STAR_DEFAULT_LADDER` (the literals) as well, or drop the claim.

MINOR | tests/test_independence_survey.py:96 | The override scan reads only `jobs/active/t2star_q1_070423.py`. There are now two declared members of the `t2star` family and the recipe makes per-job overrides a first-class feature, so an override in `t2star_q1_100423.py` alone is invisible to every carve control. (Pre-existing narrowness, newly load-bearing.) | Loop over both `jobs/active/t2star_*.py`.

MINOR | jobs/composite/independence_survey.py:70-76 | Stale as shipped after this diff: "the literal `5e-6` the T2* jobs write". The T2* jobs write no literal any more; the ladder lives in `recipes.T2STAR_THRESHOLDS` as the same comprehension. | Retarget the comment at `recipes.T2STAR_THRESHOLDS` / `T2STAR_DEFAULT_LADDER`.

MINOR | jobs/bench/probe_unresolved.py:17 | "Every job file declares a ladder byte-identical to T2STAR_DEFAULT_LADDER" is false after the collapse - no job file declares a ladder. | Reword.

## Reviewed
- [x] src/quebra/recipes.py (T2* section)
- [x] jobs/active/t2star_q1_070423.py
- [x] jobs/active/t2star_q1_100423.py
- [x] tests/test_independence_survey.py

--- PATCH COMMIT 1 of 4 (build/parse breakage + Makefile incoherence), appended 2026-09-02 ---

SCOPE: pyproject.toml, .gitignore, Makefile, .github/workflows/ci.yml, src/quebra/cli.py
COMMIT: 1241725

## Manifest
- [ ] pyproject.toml (+2/-1)
- [ ] .gitignore (+4/-1)
- [ ] Makefile (+8/-1)
- [ ] .github/workflows/ci.yml (+5/-28)
- [ ] src/quebra/cli.py (+1/-5)

## Reviewed

## Findings

### pyproject.toml

VERIFIED | pyproject.toml:80 | The dot fix is real and load-bearing, and worse than reported. MEASURED: `tomllib` on `git show HEAD:pyproject.toml` -> "Expected newline or end of document after a statement (at line 79, column 17)"; the working-tree file parses and `tool.ruff.lint` comes back as a table with key `select`. Additionally ruff itself HARD-FAILS on the HEAD form (reproduced on a two-line copy: "ruff failed / Cause: Failed to parse ... TOML parse error at line 1, column 17"), so at HEAD `make lint` and therefore all of `make check` were unrunnable, not merely `pip install`. Because nothing ran before, the fix cannot have changed any lint outcome.

VERIFIED | pyproject.toml:19 | 3.13 classifier is consistent with every other version claim in the tree, checked: `requires-python = ">=3.11"`, `[tool.ruff] target-version = "py311"` (floor, correct), `[tool.mypy] python_version = "3.11"` (floor, correct), README.md:95 "Python 3.11 or newer", AGENTS.md:5 "Python 3.11+". Nothing else needed updating. All five gates pass on this machine's 3.13 (`ruff check`, `ruff format --check` 147 files, `mypy` 11 files, `lint-imports --no-cache` 1 kept, `deptry .` 110 files), so the new matrix leg is very likely green.

## Reviewed
- [x] pyproject.toml

### .gitignore

VERIFIED | .gitignore:51 | Pattern syntax and position are correct (trailing-slash directory match, filed under the existing "# Tool caches" heading beside `.grimp_cache/`), and it shadows nothing. Load-bearing for the sdist, MEASURED both ways: built an sdist from a `git archive HEAD` tree with the real `.import_linter_cache/` copied in and the dot fix applied -> ships exactly 5 entries (`.gitignore`, `CACHEDIR.TAG`, `quebra.meta.json`, and two mtime-keyed `*.data.json` import-graph caches); built from the working tree -> `grep -Ei "import_linter|grimp_cache|__pycache__|\.pyc"` over 194 entries returns nothing. Also confirms the `.grimp_cache/` claim by construction: hatchling reads only the ROOT `.gitignore`.

MINOR | .gitignore:50-51 | Comment states the redundant reason first. "Ignored so it cannot be committed and so hatchling keeps it out of the sdist" - the commit half was already true and is not what this rule buys: `.import_linter_cache/.gitignore` contains `*` ("# Automatically created by Grimp."), and `git check-ignore -v` confirms every path inside was already blocked, which is why `git status` looked clean. The ONLY live leak was the sdist, because hatchling honours the root `.gitignore` and not nested ones. A reader who believes the stated reason may delete the nested Grimp file as superfluous, or may not learn the nested/root asymmetry that caused the leak. | Say: hatchling honours only the root `.gitignore`, so the nested one Grimp writes keeps the cache out of git but not out of the sdist.

## Reviewed
- [x] .gitignore

### Makefile

VERIFIED | Makefile:47 | Both halves of the premise reproduce, and the fix is correct GNU make. MEASURED on a `git show HEAD:Makefile` copy: `make -f <head> -n test-r` prints TWO lines, `pytest -m "r"` then `pytest -m "real or regression"` (a non-tab `#`/`##` comment line does not terminate a recipe); `make -f <head> -n test-real` prints "make: Nothing to be done for 'test-real'." exit 0, i.e. the silent no-op, because `test-real` is in `.PHONY` with no rule. Working tree: one line each. No other target has the defect - `make -n` on all of check/lint/types/arch/deps/test/test-all/test-r/test-real/docs/clean prints exactly its own recipe, and `promote` carries its own header. `.PHONY` lists all 12 targets.

VERIFIED | Makefile:44 | The "exit 5" claim is accurate. MEASURED: `pytest -m "r"` -> "419 deselected", exit 5; `pytest -m "real or regression"` -> "419 deselected", exit 5.

IMPORTANT | Makefile:48 | `regression` is not a declared marker. `[tool.pytest.ini_options].markers` declares only `slow`, `heavy`, `real`, `r`, and `addopts = "--strict-markers"`. So half this selector is a dead name in two directions: today `-m "real or regression"` silently evaluates `regression` false (`-m` expressions are not strict-marker-checked), and the moment anyone writes `@pytest.mark.regression` collection ERRORS under `--strict-markers`. This diff is promoting the recipe to a live target for the first time, so it is shipping a selector the project's own config forbids completing. | Declare `regression` in `markers`, or drop it from the selector.

IMPORTANT | Makefile:46 | "## Markers are assigned in the test-architecture phase." is project chronology in a comment - the exact category the author's hard rule bans (no spec numbers, no phase references). It also dates: the line becomes false and unnoticeable the moment markers land. | Delete it, or restate as state rather than plan: "no test carries these markers yet, so both selectors are empty."

MINOR | Makefile:41-45 | The two paragraphs undercut each other. Para 1 makes the missing header sound like a live leak - "that is how the private-data selector reaches a public runner". Para 2 says `pytest -m "r"` exits 5. Since make aborts on the first failing recipe line, at HEAD the attached second line could NEVER have executed: the leak only becomes live once an `r`-marked test exists and line 1 starts passing. As written a reader concludes an active exposure was just closed. | State that the exposure is latent, live from the first `r`-marked test.

MINOR | Makefile:47-48 | The commit claims "no behaviour change"; this target changes one. `make test-real` moves from exit 0 (nothing to be done) to exit 5. That is the intended fix and the comment owns it, but the commit message as stated is wrong. `make test-r` is fine - its exit code was already 5 from line 1, so dropping line 2 is observationally inert. | Say in the message that `make test-real` now fails loudly instead of passing vacuously.

MINOR | Makefile:40 | "Local only. Touches data/real_private/, so it must never run on GitHub Actions" justifies only the `real` half of `real or regression`. Nothing states what `regression` needs or why it may not run in CI. | Cover both, or split the targets.

## Reviewed
- [x] Makefile

### .github/workflows/ci.yml

VERIFIED | ci.yml:34-38 | CI cannot reach the private-data selector. The workflow has exactly two run steps after install: `make check` and `make deps`. `make -n check` expands to `ruff check .`, `ruff format --check .`, `mypy src/quebra/core`, `lint-imports --no-cache`, `pytest -m "not slow and not heavy and not real and not r"`. Neither `test-r` nor `test-real` is a prerequisite of anything CI invokes. `actions/setup-python@v5` supports 3.13.

IMPORTANT | ci.yml:17-18 | "it is the interpreter that produced every shipped number" is a provenance claim that no shipped artifact can support. MEASURED: a `prov.json` under `output/` carries keys `dataset_hash, dataset_hashes, dataset_path, dataset_paths, figure_node_label, git_commit, includes, job_file, job_file_hash, node_name, pipeline_steps, targets_rendered` - no interpreter field; `grep -rl "3\.11\|3\.12\|3\.13" output/*/provenance/` matches nothing; `grep` over `src/quebra/` finds no `sys.version` / `platform.python_version` capture anywhere. The rule is that a measured claim in a comment comes from the artifact it cites in the state it ships; this one cites nothing and rests on memory. | Drop the clause, or record `python_version` in provenance in a later commit and then make the claim.

IMPORTANT | ci.yml:19 | "and it went untested while the classifiers omitted it" is past-tense project chronology - "it used to" - in a diff whose other half is the author removing exactly that from this file. | Delete the clause; "3.13 is the interpreter this is developed on" already carries the reason.

IMPORTANT | ci.yml:35-56 (deleted block) | The deletion takes two things with it that are not chronology and are now recorded nowhere. (a) The prohibition "no step fetches data/real_private/, and none may be added" - the only statement of that invariant in the workflow. (b) The corrected mechanism for why CI is green without the private tree: per-test `pytest.skip` in `tests/test_data_manifest.py` (4 sites) and `tests/test_artifact_guard.py`, NOT the marker selector. That correction was the output of a prior review gate (logged in this ledger against ci.yml:38-39), and it remains true today - re-measured: `pytest -m "r"` reports 419 deselected, zero tests carry any of the four markers, so `FAST` excludes nothing. Deleting it means the next reader re-derives the false version. | Keep two lines: CI must never fetch `data/real_private/`, and the guard is the per-test `pytest.skip` on its absence, the markers being currently unused.

MINOR | ci.yml:3 vs :35 (deleted) | Inconsistent trim. The diff removes "SPEC 0004 R4.4.1" from the step comment but leaves "SPEC 0002 R1.4.2" in the file header nine lines above, so the file now applies the no-spec-numbers rule to one comment and not to the other. | Trim both or neither.

## Reviewed
- [x] .github/workflows/ci.yml

### src/quebra/cli.py

VERIFIED | cli.py:9,207 | The removed branch was unreachable AND was not a fallback. `h5py>=3.11` is a hard `dependencies` entry, and the stale local editable metadata DOES list it - `importlib.metadata` gives `Requires-Dist: [... 'h5py>=3.11' ...]`; what the stale metadata is missing is `grimp`, plus it says `import-linter>=2.0` where pyproject now says `>=2.13`. So the staleness does not touch h5py and no reachable environment loses it. Separately, the old `and h5py` guard did not provide a fallback: with h5py absent an `.h5` path fell through to `_schema_stub(load(file_path))`, and `loaders/registry._load_hdf` does its own bare `import h5py`, so the identical ImportError was raised a few frames later with a worse message. Deleting it removes a silent-fallback shape, per the project rule.

VERIFIED | cli.py:9 | No import invariant is broken and `code_hash()` is structurally unable to move. The "must not load a rendering stack" rule exists only as a comment at `src/quebra/core/job.py:15`, scoped to that module; grep over `tests/` finds no `sys.modules`-based import-weight assertion for job.py or cli.py, so nothing covers cli.py. And cli.py is already the heaviest module in the package: it imports `quebra.plots.targets`, so `import quebra.cli` already pulls matplotlib AND plotly - MEASURED 860 modules, 0.62 s. h5py adds ~0.16 s on top (`-X importtime`: 170 ms cumulative). Nothing in `src/` or `jobs/` imports `quebra.cli` (grep clean), and the import-linter layers contract puts `quebra.cli` at the top and passes, so cli.py can never enter a step's static import closure; `core/closure.py` hashes only `quebra.*` modules reachable from step functions. The byte-identical `code_hash()` claim holds for a structural reason, not luck.

IMPORTANT | src/quebra/recipes.py:64-70 (omitted from the diff) | Self-consistency gap. `run_start_unix_s_from_hdf5` still carries the identical dead guard on the identical mandatory dependency: `try: import h5py / except ImportError as exc: raise ImportError("h5py is required to read measurement_time from HDF5") from exc`. That is a catch-and-rewrap of an error that cannot occur. The stated reason for the cli.py edit - "the fallback branch could never execute" - applies to it verbatim. The repo now holds three conventions for one hard dep: module-level (cli.py:9), bare function-level (loaders/registry.py:73), function-level-with-rewrap (recipes.py:65). | Delete the try/except in recipes.py in this same commit, so the rule is applied once rather than half-applied.

MINOR | cli.py:9 | Module-scope import widens the blast radius in a broken environment: an h5py-less install now fails at `import quebra.cli`, taking `quebra run` and `quebra inspect` down with it rather than only `schema-wizard` on an `.h5` path. Loud is the right direction and this is unreachable via a supported install, but note `make docs` (`sphinx-build -W --nitpicky`; sphinx is in no extra) would autodoc-import `quebra.cli` in whatever minimal env a docs build has. | No change needed; do not add sphinx without adding h5py to that env.

MINOR | cli.py:207 | The only code line the commit touches has no test. Grep over `tests/` and `scripts/` for `schema-wizard`, `schema_wizard`, `_schema_stub`, `cli.main` returns nothing; the three tests that touch this module (`test_nesting.py:18`, `test_artifact_guard.py:174`, `test_path_resolution.py:198`) only import `_module_from_path`. They do exercise the new top-level import on every FAST run, which is the part most likely to break, so the exposure is small - but the no-behaviour-change guarantee is asserted, not tested. | Note it; a `schema-wizard` test is its own commit.

MINOR | cli.py:207-227 (pre-existing, newly exposed) | The `if` now reads as a total handler for `.h5`/`.hdf5` and is not: the `len(keys) > 1` arm returns, the `<= 1` arm falls out of the `with` and re-opens the same file through `load(file_path)`. With the `and h5py` short-circuit gone there is no longer any visual cue that the branch is partial. | Add an explicit `if len(keys) <= 1` comment, or restructure. Not introduced here.

## Reviewed
- [x] src/quebra/cli.py

--- PATCH COMMIT 2 of 4 (the provenance record describes the run that actually happened) — SLICE A, appended 2026-09-02 ---

SCOPE: src/quebra/provenance.py, src/quebra/core/runner.py, src/quebra/core/job.py, src/quebra/core/closure.py, scripts/promote_run.py, tests/test_reuse_completeness.py (new/untracked), tests/test_identity_closure.py, tests/test_shape_stats.py
COMMIT: a00894d (working tree, uncommitted)

## Manifest
- [ ] src/quebra/provenance.py (+39/-33)
- [ ] src/quebra/core/runner.py (+66/-21)
- [ ] src/quebra/core/job.py (+20/-19)
- [ ] src/quebra/core/closure.py (+19/-4)
- [ ] scripts/promote_run.py (+18/-0)
- [ ] tests/test_reuse_completeness.py (NEW, 181 lines)
- [ ] tests/test_identity_closure.py (+8/-9)
- [ ] tests/test_shape_stats.py (+4/-4)

## Reviewed (slice A)

## Findings (slice A)

--- PATCH COMMIT 2 of 4 — SLICE B (paths/cli/recipes + spec amendment), appended 2026-09-02 ---

SCOPE: src/quebra/core/paths.py, src/quebra/cli.py, src/quebra/recipes.py, tests/test_data_root_resolution.py, spec/specinstallabity02.md, spec/specdatalayout03.md
COMMIT: a00894d (working tree, uncommitted)

## Manifest (slice B)
- [ ] src/quebra/core/paths.py (+~70/-~30)
- [ ] src/quebra/cli.py (+3/-5)
- [ ] src/quebra/recipes.py (+2/-5)
- [ ] tests/test_data_root_resolution.py (+43/-12)
- [ ] spec/specinstallabity02.md (+20/-1)
- [ ] spec/specdatalayout03.md (+6/-0)

## Reviewed (slice B)

## Findings (slice B)

### src/quebra/core/paths.py

CRITICAL | src/quebra/core/paths.py:234-238 (and the whole `explicit` demand branch, :246-255) | The motivating failure is NOT fixed, and the docstring now claims it is impossible. `resolve_data_root(explicit=...)` has NO production caller: `cli.py:130-131` does `Path(args.data_root).expanduser().resolve()` itself and hands the result to `run_job`, which at `runner.py:465` does `dataset_root = Path(data_root).resolve() if data_root else default_dataset_root()` - never `resolve_data_root(data_root)`. So the new raise fires only for direct-API/test callers. MEASURED on the working tree: `resolve_dataset_path("data/real_private/6D2S/070423_6D2S_qubit1.pickle", Path("/mnt/typo"))` -> `/…/qre_tool/data/real_private/6D2S/070423_6D2S_qubit1.pickle`, and `jobs/bench/results/xi_tie_experiment.csv` likewise. Every job in `jobs/active/` declares a repo-relative `Dataset.path` (`data/real_private/…` or `jobs/…`), and `quebra.toml` declares `data_root = "."`, so `resolve_dataset_path`'s repo-root second candidate silently reproduces EXACTLY the described failure: `quebra run jobs/active/t2star_q1_070423.py --data-root /mnt/typo` completes against the repository tree and records that tree's dataset hashes. The docstring sentence "which is the failure `DataRootNotFound` exists to make impossible" is therefore an overclaim about the shipped CLI. (Already logged against this repo at ledger:2025 in the SPEC 0002 pass and still open.) | Make `run_job` call `resolve_data_root(data_root)` and drop the ad-hoc `Path(data_root).resolve()`, in this commit - otherwise the commit buys nothing a user can reach. If that is deferred, delete the "makes impossible" clause and say the CLI does not route through here yet.

IMPORTANT | src/quebra/core/paths.py:192-198 | `DataRootNotFound`'s own docstring - "No data root could be resolved, **with every location that was tried**" - is now false at two of its three raise sites. The two new demand raises name only the failed demand; the explicit one does not even mention that `QUEBRA_DATA_ROOT`, `quebra.toml` and platformdirs exist, so a user who hits it learns nothing about the mechanism set. The class is the contract carrier for R1.3.4 and it was not touched. | Narrow it: "either a root somebody named does not exist - in which case the message names that one - or nothing resolved, in which case it lists every location tried."

IMPORTANT | src/quebra/core/paths.py:246-255 vs docs/WRITING_A_JOB.md:145 (omitted from the slice) | The user-facing documentation of this exact contract still says "resolved in this order, **first hit wins**" over all four mechanisms, and "If none resolves, QUEBRA raises `DataRootNotFound`" - i.e. it documents the fall-through this commit deleted. The spec files were amended; the doc a user reads was not. | Add one sentence to WRITING_A_JOB.md:145-155: mechanisms 1 and 2 raise if they name a non-directory.

MINOR | src/quebra/core/paths.py:1 vs :8-10 | The new title line - "Path anchors and resolvers, **all of them relative to the caller's project**" - is contradicted seven lines below: "It may be inside the project or outside it." The dataset root is precisely the anchor that is NOT relative to the project. | "Path anchors and resolvers. Neither is computed from this file's location."

MINOR | src/quebra/core/paths.py:4-7 | Ambiguous rewrite lands on a false reading. "…never computed from this file's location, which is the same anchor `provenance.get_git_commit` uses." Nearest antecedent of "which" is "this file's location", so the sentence asserts `get_git_commit` anchors on `__file__`. It does not: `provenance._git` runs `subprocess.run([...], cwd=Path.cwd())` (provenance.py:43, unchanged by the slice-A diff), so it anchors on the cwd - the same anchor `repo_root` uses, which is presumably the intended point. The old text asserted the same falsehood in the other direction, so the rewrite did not fix it. | "…never computed from this file's location. The cwd is also what `provenance.get_git_commit` anchors on."

MINOR | src/quebra/core/paths.py:260-271 | `QUEBRA_DATA_ROOT=""` is treated as unset (`if env:`) and the failure message then prints the literally false line `QUEBRA_DATA_ROOT: not set`. Setting it empty is the plausible way a shell script "clears" it (`export QUEBRA_DATA_ROOT=` / a CI variable defined-but-blank), and a docstring that calls this mechanism a DEMAND makes silently ignoring a set value the odd branch. No test or doc sets it empty (grepped `tests/`, `docs/`, `scripts/`, `.github/`). Note `Path("").resolve()` is the cwd, so honouring "" literally would be worse - the issue is only the diagnostic. | `tried.append(f"{QUEBRA_DATA_ROOT_ENV}: set but empty, treated as unset")` when `env == ""`.

MINOR | src/quebra/core/paths.py:251 | The message hard-codes a CLI flag for a general parameter: a library caller doing `resolve_data_root(Path("/data"))` from Python is told "--data-root was given as ...". The env-var message at :266 gets this right by naming the actual mechanism. | "an explicit data root was given as '{explicit}' (the `--data-root` flag arrives here)".

MINOR | src/quebra/core/paths.py:250 | `DataRootNotFound` subclasses `RuntimeError`, so once the CRITICAL above is fixed a typo'd `--data-root` reaches the user as a bare RuntimeError traceback rather than argparse's `parser.error` two-liner. The sibling failure `DataUnavailable` deliberately subclasses `FileNotFoundError` for handler compatibility; "you named a path that is not a directory" is the same species. | Catch it in `cli.main` and route to `parser.error`, or subclass `NotADirectoryError` for this case.

VERIFIED | src/quebra/core/paths.py:246-270 | No legitimate caller relied on the fall-through. Every caller of `resolve_data_root` with an argument is in `tests/`; the only production entry is `default_dataset_root()` -> `resolve_data_root()` with `explicit=None` (grepped src/, jobs/, scripts/, conftest.py). Nothing passes a path that is created later: `run_job(data_root=tmp_path)` in the test suite bypasses this function entirely. A path to an existing FILE raises with the correct wording ("that is not a directory"), and `Path("")` -> cwd, which is a directory, so `resolve_data_root("")` returns the cwd rather than raising - the CLI never gets there because `if args.data_root` is falsy for `""`.

VERIFIED | src/quebra/core/paths.py:25,280-284 | The `platformdirs` guard removal is correct and the comment is accurate: `platformdirs>=4.0` is a hard `dependencies` entry in pyproject.toml. Module-scope import is safe for the identity closure (`core/closure.py` hashes only `quebra.*` modules) and platformdirs is import-light.

VERIFIED | src/quebra/core/paths.py:63-64 | The measured claim is true in the state it ships: working-tree `quebra.toml` says `data_root = "."` and `default_dataset_root()` returns `/…/qre_tool` (executed). Spec-number and `qre_tool/` removals at :115,:137,:154 are complete for this file - `grep -n "SPEC \|R1\.3\|qre_tool" src/quebra/core/paths.py` is empty.

## Reviewed (slice B)
- [x] src/quebra/core/paths.py

MINOR | src/quebra/core/paths.py:148,154-155 (addendum) | The chronology sweep is half-applied inside one docstring. The diff removes "SPEC 0003 R3.4." from :154 but leaves ":148 …the two situations that **used to** look identical" three lines above - the "it used to" form the rule bans - and leaves the deleted prefix's hole unreflowed, so :154 is a 63-column line in an 88-column file. | Say "the two situations that otherwise look identical" and reflow the paragraph.

### src/quebra/cli.py

IMPORTANT | src/quebra/cli.py:119 | The one user-visible statement of the data-root default is still the pre-`data_root = "."` one: `--data-root` help says "(default: **the repo's parent directory**)". `quebra run --help` therefore tells every user the opposite of what `default_dataset_root()` returns (measured: the repository itself), and it names a mechanism the resolver no longer has. This commit exists to make exactly this claim accurate in `paths.py`; it left the copy a user actually reads untouched, in a file the commit already edits. (Logged at ledger:2047 in the SPEC 0002 pass, still open.) | help="Dataset root for relative dataset paths; defaults to QUEBRA_DATA_ROOT, then [tool.quebra] data_root in a quebra.toml at or above the cwd, then the platformdirs user data dir. Used by both loading and provenance hashing."

VERIFIED | src/quebra/cli.py:9,207 | The guard removal is correct and cannot regress anything reachable: `h5py>=3.11` is a hard dependency, nothing in `src/` or `jobs/` imports `quebra.cli`, and the old `and h5py` short-circuit was not a fallback (it fell through to `load(file_path)` -> `loaders/registry._load_hdf`, which does its own bare `import h5py`). The `make docs` worry raised against this line in PATCH COMMIT 1 is moot: `docs/` contains no `conf.py` and no `automodule`/`autodoc` anywhere, so `sphinx-build -W --nitpicky -b html docs docs/_build/html` fails before importing any project module.

## Reviewed (slice B)
- [x] src/quebra/cli.py

### src/quebra/recipes.py

MINOR | src/quebra/recipes.py:65 | The surviving inline `import h5py` states no reason, so the repo now carries two unexplained conventions for one hard dependency: module scope in `cli.py:9`, function scope here and in `loaders/registry.py:73`. There IS a good reason here - `quebra.recipes` is in the step-import closure that `core/closure.py` hashes and every job imports it, so keeping h5py out of that import keeps job-build light - but a reader who applies the cli.py precedent will "tidy" it to the top of the file. | One line: "Imported here, not at module scope: every job imports this module and only this function needs h5py."

MINOR | src/quebra/recipes.py:58 | The function the commit edits is dead. `run_start_unix_s_from_hdf5` has zero call sites and zero tests repo-wide (grepped src/, jobs/, tests/, scripts/, docs/ - the only hit is its own `def`). Editing it is harmless, but it is the only behaviour-bearing line in this file's diff, and it moved the code hash of all six measurable jobs for a function nothing executes. | Either delete the function or give it a caller/test; note in the commit message that the recipes.py hash movement comes from dead code.

VERIFIED | src/quebra/recipes.py:61-67 (removed) | Nothing depended on the removed message: `grep "h5py is required"` over the tree is empty. The rewrap was loud, not a silent fallback, so behaviour changes only from a custom ImportError text to `ModuleNotFoundError: No module named 'h5py'` - unreachable under a supported install. File is 413 lines, single-purpose, not unwieldy.

## Reviewed (slice B)
- [x] src/quebra/recipes.py

### tests/test_data_root_resolution.py

VERIFIED | tests/test_data_root_resolution.py:100-136 | Both inverted tests are load-bearing, MEASURED by mutation without touching any source file: a pytest plugin re-bound `quebra.core.paths.resolve_data_root` to a verbatim copy of the pre-change fall-through implementation, and both new tests FAILED ("DID NOT RAISE DataRootNotFound") while the other 6 passed. On the working tree, `pytest tests/test_data_root_resolution.py` -> 8 passed. No env leak: ambient `QUEBRA_DATA_ROOT` is unset here, and `finally: del os.environ[...]` runs before `isolated`'s `monkeypatch.delenv` undo, so the two orders compose; running the env test followed by `test_when_nothing_resolves...` in one process passes.

IMPORTANT | tests/test_data_root_resolution.py (whole file) | The suite pins the FUNCTION and leaves the shipped BEHAVIOUR untested, which is why `make check` is green while the motivating hazard is live (see the CRITICAL against paths.py). There is no test that `run_job(job, out, data_root=<nonexistent>)`, or `quebra run --data-root <typo>`, refuses to run - and it does not refuse; it silently analyses the repo tree. A test at that level is what would have caught that the demand branch has no production caller. | Add one test: `run_job` with a nonexistent `data_root` must raise, not resolve datasets through the repo-root fallback.

MINOR | tests/test_data_root_resolution.py:131-136 | Hand-rolled `os.environ[...] = ...` + `try/finally` where the file's three sibling tests (:45, :56, :117) use `monkeypatch.setenv`, and where `monkeypatch` was deliberately dropped from the signature to make it possible. It is not leaking today, but it re-introduces exactly the anti-pattern already logged against this file at ledger:2129 (raw `os.chdir` at :67,:82), and a future edit that adds a second `os.environ` write or an early `return` inside the block loses the undo. Nothing in the test needs the env var to survive `monkeypatch`'s undo, so there is no reason for the deviation. | `monkeypatch.setenv(QUEBRA_DATA_ROOT_ENV, str(tmp_path / "gone"))` and delete the try/finally.

MINOR | tests/test_data_root_resolution.py:118,133 | The two `match=` strings pin prose at different strengths and neither pins the mechanism. `match="demand, not a candidate"` binds a copy-editable clause: any rewording of that sentence turns a passing test red for no behavioural reason. `match="not a directory"` is the opposite problem - BOTH demand raises contain that phrase, so it does not discriminate which mechanism fired (it only separates the demand raises from the catch-all, which has no such phrase). | Match on the stable identifiers instead: `match=r"--data-root"` for the explicit test, `match=QUEBRA_DATA_ROOT_ENV` for the env test.

MINOR | tests/test_data_root_resolution.py:3-4 vs :87,:142 | The chronology/oracle sweep is half-applied here too. The module docstring drops "SPEC 0002 R1.3.3 ... R1.3.4" in favour of "the documented resolution order" - which now cites nothing checkable, so a test file loses its oracle - while :87 still says `"""R1.3.4. ..."""` and :142 still says "The declared value changed in SPEC 0003 - `data_root` moved from `".."` to `"."`", which is the "it used to" form. Either tests may cite specs or they may not; right now this one file does both. | Pick one. If specs stay out of tests, point the oracle at `resolve_data_root`'s docstring and restate :142 as "asserts the checkout resolves to whatever its own quebra.toml declares".

MINOR | tests/test_data_root_resolution.py:3-4 | "the requirement that a failure name every location tried" is now only true of the catch-all raise; the two new raises name one location. The docstring states the old, unnarrowed contract as the file's oracle two lines above the tests that break it. | Say "a failure names either the demand that was not met or every candidate tried."

MINOR | tests/test_data_root_resolution.py:100-136 | Untested edge case introduced by this change: an explicit root that is an existing FILE rather than a missing path. It takes the same branch and the message ("that is not a directory") is correct for it, but nothing pins that - and it is the likelier operator mistake (`--data-root ./data/file.csv`, tab-completion). | Add `(tmp_path/"f").write_text(""); pytest.raises(DataRootNotFound, match="not a directory")` on `resolve_data_root(tmp_path/"f")`.

## Reviewed (slice B)
- [x] tests/test_data_root_resolution.py

### spec/specinstallabity02.md

IMPORTANT | spec/specinstallabity02.md:113-117 (R1.3.3a rationale) vs :128 (R1.3.5) | The amendment overclaims at the system level. Its rationale is written about `--data-root` ("With `quebra.toml` here declaring `data_root = "."`, the tree it silently reaches is the repository itself"), which only holds if R1.3.5 - "`--data-root` maps to mechanism 1" - is actually implemented. It is not: `cli.py:130` resolves the flag itself and `runner.py:465` does `Path(data_root).resolve()`, so `resolve_data_root`'s explicit branch is never reached from the CLI and the described silent-wrong-tree run still happens (MEASURED, see the CRITICAL against paths.py). So the amendment is accurate about the FUNCTION and false about the tool. R1.3.5 is stated as a requirement in the same section and is unmet, unmarked. | Either wire `run_job` through `resolve_data_root` in this commit, or add one sentence to R1.3.3a: "R1.3.5 is not yet implemented - the CLI bypasses `resolve_data_root`, so the demand rule does not yet protect `--data-root`."

MINOR | spec/specinstallabity02.md:107-110 vs the code at paths.py:261 | R1.3.3a says "If either names a path that is not a directory, `resolve_data_root` raises". The code has an unstated carve-out: `if env:` means `QUEBRA_DATA_ROOT=""` is treated as UNSET and falls through to mechanisms 3 and 4, so an env var that is set-but-empty is the one demand that silently does not raise. | State the carve-out: "an empty `QUEBRA_DATA_ROOT` counts as unset, because `Path("")` is the cwd and honouring it literally would be the very guess this rule forbids."

MINOR | spec/specinstallabity02.md:127-128 (R1.3.4, unamended) | R1.3.3 got an "AMENDED — see R1.3.3a" marker; R1.3.4 did not, although its text ("If no root resolves, raise a named exception stating ... every location tried") is exactly what the two demand raises no longer satisfy. The narrowing lives only inside R1.3.3a, and R1.3.4 is what other artifacts cite by number - `tests/test_data_root_resolution.py:87` is a docstring reading `"""R1.3.4. ..."""`. A reader arriving at R1.3.4 gets the unnarrowed contract. | Add "**NARROWED — see R1.3.3a**: applies to mechanisms 3 and 4" to R1.3.4.

VERIFIED | spec/specinstallabity02.md:118-119 | The claim "R1.3.4's requirement is unchanged for mechanisms 3 and 4: when nobody named a root and none is found, the exception still lists every location tried" is TRUE of the code: the catch-all at paths.py:286-291 still appends all four entries, including "explicit argument: none passed" and "QUEBRA_DATA_ROOT: not set", and `test_when_nothing_resolves_it_raises_and_names_every_location` passes (executed). The R1.3 acceptance bullet "its message lists all four locations tried" therefore still holds unmodified.

## Reviewed (slice B)
- [x] spec/specinstallabity02.md

### spec/specdatalayout03.md

MINOR | spec/specdatalayout03.md:36-37 vs spec/specinstallabity02.md:107-108 | The two specs state the same new rule with two different predicates: here "a path they name **that does not exist** raises", there "names a path **that is not a directory**". The code tests `root.is_dir()`, so an existing FILE also raises - which this file's wording says it should not. Whichever spec a future reader consults first, one of them is wrong. | Use "is not a directory" here too.

VERIFIED | spec/specdatalayout03.md:35-37 | The cross-reference is correct and the closed decision is re-opened explicitly rather than quietly: R1.3.3a exists at the cited location, the order and the four mechanisms are genuinely unchanged (only the miss-handling of 1 and 2 changed), and the "Do not re-open it" sentence at :32 is now scoped rather than contradicted. Both spec additions are short amendments to existing planning docs, not new unasked documents.

## Reviewed (slice B)
- [x] spec/specdatalayout03.md

## Slice B verdict: DO NOT SHIP (1 CRITICAL: the shipped CLI still does the thing the commit
## claims to have made impossible; the fix is one line in runner.py:465).

### src/quebra/provenance.py

VERIFIED | provenance.py:31-49 | The `_git` collapse is mechanically sound. `TimeoutExpired.__mro__` and `CalledProcessError.__mro__` both go through `SubprocessError`; `FileNotFoundError` is an `OSError`; `SubprocessError` is NOT an `OSError`, so both arms are needed and neither is redundant. MEASURED: `subprocess.run(check=True, capture_output=True, timeout=0.2)` on `sh -c 'sleep 5'` raises `TimeoutExpired` at 0.2 s and is caught by `(OSError, subprocess.SubprocessError)`. `cwd=Path.cwd()` is evaluated INSIDE the try, so the docstring's "a deleted cwd" claim is real.

VERIFIED | provenance.py:63 | The empty-string case is read correctly, not as failure. `_git` returns `completed.stdout` unstripped, so a clean tree yields `""`, and `out is not None and out.strip() == ""` -> True. MEASURED end-to-end by the new test file (`test_an_untracked_file_does_not_make_the_tree_dirty`, passes). `get_git_commit`'s `(out.strip() or "nogit") if out is not None else "nogit"` handles both None and `""`.

IMPORTANT | provenance.py:28-29 | The timeout's stated mechanism cannot occur for either command that routes through `_git`. "`git` can block indefinitely waiting on a credential prompt" applies to commands that contact a remote; `rev-parse --short HEAD` and `status --porcelain --untracked-files=no` are purely local and reach no credential helper. A reader who believes this will not think to check the reachable hangs (a stalled network filesystem holding the worktree, an fsmonitor daemon, a very large worktree). The rule is that a comment says why the code does what it does; this one gives a false why. | Name a reachable hang, or keep only "a provenance helper that hangs stops the run it exists to describe".

IMPORTANT | provenance.py:52-53, 62-64 with core/runner.py:227,231 | "nogit" is a sentinel that MATCHES ITSELF, and this collapse newly makes `(git_commit="nogit", tree_clean=True)` a reachable, self-consistent state, which turns the commit half of the reuse gate into a tautology. MEASURED in a fresh `git init -q` directory: `git status --porcelain --untracked-files=no` exits 0 with EMPTY output (-> `is_tree_clean() is True`) while `git rev-parse --short HEAD` exits 128, "fatal: Needed a single revision" (-> `"nogit"`). Before the collapse `is_tree_clean` read `Path(__file__).parent`, so from a wheel this pair was unreachable (site-packages is not a repo -> False). Now: run twice in an unborn repo and `_reuse_eligible_dir` sees `rec_commit == git_commit == "nogit"` and `rec_tree_clean == True` -> SKIP, with zero code-version guarantee. `job_code_hash` folding the closure narrows but does not close it: any module outside the step closure (`plots/*`, every plot class) can be edited with the identity byte-identical, and the stale figure is reused. | `_reuse_eligible_dir` returns None when `git_commit == "nogit"` (or `_read_prov_reuse_fields` treats a recorded "nogit" as ineligible).

IMPORTANT | provenance.py:44 | The widened except writes a provenance FALSEHOOD in a case that was previously loud. Under `(FileNotFoundError, CalledProcessError)` a `PermissionError` on the cwd propagated and there was no timeout to expire; under `(OSError, SubprocessError)` both collapse to `None` -> `get_git_commit()` returns `"nogit"`, which is then recorded as this run's `git_commit`. A run that really did happen at a real commit ships a record asserting it happened outside version control, and nothing in the record distinguishes "not a repo" from "git timed out / git not executable". In the commit whose title is "the provenance record describes the run that actually happened", that is the wrong direction. `tree_clean=False` keeps it reuse-safe, so this is a record-accuracy defect, not a wrong-artifact one. | Either keep the narrow excepts for the commit read, or record a distinct value (e.g. "git-unavailable") so the two are distinguishable in the record.

MINOR | provenance.py:37 | `cwd=Path.cwd()` is the `subprocess.run` default. The only thing the explicit form buys is that `Path.cwd()` is evaluated inside the try, which is what makes the docstring's "deleted cwd" OSError path real - nothing says so, so a future reader deletes it as redundant and loses that path. | One clause on the line.

MINOR | provenance.py:24-28 | The comment justifies `Path.cwd()` entirely by describing what a `__file__`-anchored helper would do, but `__file__` no longer appears anywhere in the module. It is rationale rather than chronology, so it does not break the hard rule, but the reader has no code to attach half the paragraph to. | Shorten to the invariant: both halves of the gate must describe the repository the run is happening in, which is the cwd.

VERIFIED | provenance.py:88-94 (job_file_hash -> job_code_hash) | The rename breaks NO reader. MEASURED: 909 files under `output/**/provenance/` carry the old key, 0 carry the new one, and grep over `src/`, `scripts/`, `tests/`, `docs/`, `README.md`, `AGENTS.md`, `CONTRIBUTING.md` finds no reader of either key - the only consumers of a record are `_read_prov_reuse_fields` (identity, git_commit, tree_clean), `promote_run._load_records` (identity, job_file, tree_clean) and `_mermaid_graph` (dataset_*, includes, pipeline_steps, git_commit, figure_node_label). `published/` holds 0 records, so no COMMITTED audit trail carries the old name. All 909 old-key records are also already reuse-ineligible after this commit because `git_commit` moves. The rename is a clean break.

MINOR | provenance.py:88-94 | Nothing in the tree states that records under `output/` written before this commit carry `job_file_hash`. Harmless today (no reader, see above), but the record format has no version field, so the next reader of an old directory has no way to learn why the key differs. | One line in the field comment, or a `record_format` key.

## Reviewed (slice A)
- [x] src/quebra/provenance.py

### src/quebra/core/runner.py

VERIFIED | runner.py:173-193 | `_expected_sink_pkls` derives the RIGHT set, checked against the sink loop line by line. Figure sink, `render_figures=True`: runner.py:568 writes `job_out_dir / f"{sink.input.node_id}.pkl"`. Figure sink, `render_figures=False`: runner.py:575-577 sets `prov_name = sink.input.node_id` and writes `f"{prov_name}.pkl"` - the SAME name. Materialize: runner.py:601 writes `f"{sink.name}.pkl"`. So the docstring's "mode-independent" claim holds. MEASURED over all 11 real jobs: `_expected_sink_pkls` matches the on-disk pkl set for every one (e.g. `instrument_validation` -> {instrument_validation.pkl, instrument_validation_report.pkl} for 5 sinks, correctly deduping 4 figures off one node).

VERIFIED | runner.py:494-503 | The "render_figures=False run reused by a render_figures=True call, producing no PDFs" scenario is UNREACHABLE, not a regression and not pre-existing-live. `render_figures` has exactly one non-default call site, runner.py:326 (`render_figures=inc.figures`), on the nested sub-job path, and that call also passes `force=True`, so the gate is skipped there. `cli.py` never passes it (default True). And nested runs are written under `job_out_dir / "subjobs_output"`, which is not a DIRECT child of `out_dir`, so the non-recursive `out_dir.glob(...)` at runner.py:498 cannot see them. Separately, the new filter is a pure conjunction on the candidate generator, so it can only REJECT candidates the old glob accepted - it is structurally incapable of introducing a new acceptance.

VERIFIED | runner.py:494-503 | The filter is also the right side of safe when it disagrees with identity. A materialize sink renamed on the same node leaves the identity unchanged (`parameter_row`/closure do not fold sink names) but changes the expected set, so the old directory is now incomplete and the job re-runs. Old behaviour reused it and served the artifact under the wrong name.

IMPORTANT | runner.py:492-493 | Project chronology in a comment, the exact category the author's hard rule bans: "the standalone path used to glob on directory name alone, so a run that died mid-sink-loop stayed eligible". Same rule PATCH COMMIT 1 was pulled up on at ci.yml:19. | State it as an invariant: "a candidate must hold every pkl a finished run writes; a run that died mid-sink-loop leaves a matching provenance record and would otherwise stay eligible forever."

MINOR | runner.py:181-187 | Same category, softer form: "so a run that raised on its third sink still holds a perfectly readable provenance record from its first ... without this test the gate reads a half-written run as reusable". This is a counterfactual about absent code rather than history, which is defensible, but "`run --all` makes it likely, because it catches per job and carries on" is an unverified behavioural claim about a code path in another module. | Cite it (`cli.py` line) or drop the sentence.

CRITICAL | runner.py:141-158 with jobs/composite/check_ledger_q1.py:112-118 | The docstring names a live provenance-overwrite defect and then EXEMPTS the only instance of it in the repo, calling that instance "harmless". The docstring says "BOTH write `provenance/{name}.prov.json`. So a figure title and a materialize name that coincide overwrite one another's provenance record even though their artifacts do not collide" - and then "Duplicates resolving to the same source node are harmless: that is one artifact requested twice." `check_ledger_q1` declares, per dataset, a FIGURE and a MATERIALIZE with the SAME `name` off the SAME node (`q1_27h_0704_dataset_check_ledger`, and the 1004 twin), so `existing == source_id` and the check passes - yet they are NOT one artifact requested twice: they emit two DIFFERENT records to one path. The materialize is declared second, so it wins. MEASURED on a shipped run, `output/check_ledger_q1_7174d1_20260812_095729/`: the directory holds `q1_27h_0704_dataset_check_ledger_static.pdf`, `..._academic.pdf` and the 1004 pair - four rendered PDFs - while BOTH surviving records say `targets_rendered = []` and `figure_node_label = None`. The provenance record does not describe the run that actually happened, which is this commit's title. | Key on `sink.name` for COLLISION (as done) but reject same-name duplicates whose PROV CONTENT differs - i.e. treat a figure and a materialize sharing a name as a collision regardless of source node - or give the figure record a distinct prov name. Also correct the "harmless" sentence.

IMPORTANT | runner.py:141-171 | This change turns a check that was structurally dead for figures (`x != x`) into one that can now `raise ValueError` at run start for any figure sink, and there is NO test for it anywhere: grep over `tests/` for `_check_sink_artifact_names` and for the message "artifact name collision" returns nothing. A slice that ships a new test file for the other two changes leaves the one newly-live rejection path untested, and a false positive here aborts `run --all` at job import. I did the falsification by hand - MEASURED: all 11 jobs under `jobs/` pass the new check, and the only verdict that changed is the correct one (two figures off DIFFERENT nodes with the same name: previously keyed on `node_id`, so two distinct keys and NO error; now one key, two sources, rejected). Two figures off the SAME node with different names is unchanged (allowed) under both keyings. | Add two tests: different-nodes-same-name rejects, same-node-different-names passes.

MINOR | runner.py:167-169 | The error message lost the information that made it actionable. It was "'{basename}.pkl' would be written for two different nodes"; it is now "'{sink.name}' would be written for two different nodes" - a bare name with no indication of WHAT would be written, when the whole point of the new keying is that one name spans a `.pkl`, `{name}_{target}.pdf` files and a `.prov.json`. | "'{sink.name}' names the artifacts and provenance record of two different nodes (...)".

MINOR | runner.py:153-155 | The `_safe_name` claim is accurate but incomplete, and the omission makes the collision LESS visible than it is. `_safe_name` (job.py:201-204) also `.strip("_")` and falls back to `"node"`, so beyond `"a b"`/`"a_b"` it also collapses `"_a_b_"` and maps `""`, `"***"`, `"---"`... to a shared literal `"node"`. Confirmed applied at both constructors (job.py:440, 450), so `sink.name` really is pre-sanitised as claimed. | Mention the strip and the `"node"` fallback.

MINOR | runner.py:494-503 vs :263-292 | The commit's stated invariant ("a candidate must be COMPLETE") is applied only to the standalone gate. `_cached_runs` still filters on the single `{node_name}.pkl` it needs, so a composite under `--reuse-deps` will happily reuse the one artifact it wants out of a directory left by a run that died on a LATER sink. That is content-correct for the artifact requested and I am not asking for a change, but the comment at runner.py:492 asserts the composite path "already filters this way", which conflates "requires the artifact it needs" with "requires the run to have finished". | Say which of the two the composite path checks.

## Reviewed (slice A)
- [x] src/quebra/core/runner.py

### src/quebra/core/job.py

VERIFIED | job.py:218,222,240-251,299 | Rename is complete and internally consistent: `_code_hash` -> `_job_code_hash` (3 sites), `code_hash` -> `job_code_hash` (def + 1 internal caller at :299 + error message at :242). Grep over `src/`, `scripts/`, `tests/`, `jobs/` finds no surviving `\.code_hash\(` or `_code_hash`. Memoisation semantics unchanged; the rewritten docstring's three-item list matches the three `parts` contributions at :243, :249, :250 exactly.

VERIFIED | job.py:222-234 | The docstring edit is the correct handling of the chronology rule: "Before SPEC 0005 R5.1 this hashed only the file" became "The job file alone WOULD NOT BE the code that produced the result", and "R5.1.7:" was dropped from the inline comment at what is now :249-250. This is the pattern the other two docstrings in the slice should have followed.

MINOR | job.py:222 | "Named for the three things it folds" is a name-justifies-itself claim that the name does not actually carry: `job_code_hash` says "code hash of the job", which covers items 1 and 2 but says nothing about item 3 (the step arguments). The old name's problem was `file`; the new name is better but the docstring oversells it. | "The identity's `code` contribution. It folds three things:".

## Reviewed (slice A)
- [x] src/quebra/core/job.py

### src/quebra/core/closure.py

NOTE | scope | This file carries a SEVENTH substantive change not in the six-item brief: `_module_digest` switches from `__import__` to `importlib.util.find_spec`. Reviewed here.

VERIFIED | closure.py:167-178 | The switch achieves its stated goal for the heavy dependencies. MEASURED on a cold interpreter: `importlib.util.find_spec("quebra.analyzers.checks.c1_lewis_robinson")` returns `origin` = that file and leaves `scipy`, `sklearn` and `matplotlib` all unloaded. Error handling is not a silent fallback - `(ImportError, AttributeError, ValueError)` is re-raised as a named `ValueError` with the module name, and `spec is None` falls into the existing "no source file to hash" raise.

IMPORTANT | closure.py:161-172 | "Located rather than imported" is false for ancestor packages, and the docstring's own failure argument therefore still applies. MEASURED, same run as above: `find_spec` on that leaf newly imported `quebra`, `quebra.analyzers`, `quebra.analyzers.checks`, `quebra.analyzers.checks.result` AND `numpy` (159 modules total) - because `find_spec` must import each ancestor to read its `__path__`, and `quebra/analyzers/checks/__init__.py` re-exports from `.result`. So the claim "a module-scope failure in a reachable-but-unused module would surface as a failure of IDENTITY COMPUTATION" is still true for any ancestor package, and the property is not structural: it holds only while every `quebra` package `__init__.py` stays trivial (two of the four non-empty ones already are not: `analyzers/checks` imports `result`, `_fixtures` is 2 kB). | Say "the leaf module is located rather than imported; its ancestor packages are still imported, so package `__init__.py` files must stay import-light", and consider an import-weight test on the package inits.

MINOR | closure.py:161-172 | The docstring names sklearn/scipy/plotting as what importing would pull, which is a measured-shaped claim with no locator. It is directionally right (verified above that they are avoided) but the reader cannot check it. | Name the module that does it, e.g. "sklearn via `analyzers/tlf.py`".

MINOR | closure.py:68-69 | The replacement comment drops a specific measured number ("81 modules in ~80 ms cold") for a vague one ("tens of milliseconds cold"). Dropping a number that can go stale is the right instinct, but "tens of milliseconds" is still a measurement with no source and no way to check it. | Either drop the cost claim entirely or cite how it was measured.

MINOR | closure.py:174-176 | Newly reachable hard failure that the old code could not produce: `find_spec` raises `ValueError` when a name is already in `sys.modules` with `__spec__ is None`, and the handler converts that to "could not be located to hash". Not reachable for real `quebra.*` modules (grimp only yields importable ones), so this is a note, not a defect - but the message will be misleading if it ever fires, because the module WAS located, it just has no spec. | Distinguish the ValueError arm in the message.

## Reviewed (slice A)
- [x] src/quebra/core/closure.py

### scripts/promote_run.py

VERIFIED | promote_run.py:90-105 | The new check CANNOT reject a legitimate run, checked both ways. (a) Within one run, `_emit_prov` is called from exactly two sites (runner.py:583, 606) and both pass the same `identity` and `git_commit` locals computed once at runner.py:466-468, so every record of one run agrees by construction. (b) A composite's nested sub-job records land in `job_out_dir / "subjobs_output" / <subjob-dir> / "provenance"`, and `_load_records` globs `run_dir / "provenance" / "*.prov.json"` NON-recursively, so nested records are never in the set being compared. Confirmed on disk: `output/check_ledger_q1_7174d1_20260812_095729/provenance/` holds only the composite's own 2 records, with `subjobs_output/` beside it. The all-records-missing-`identity` case degrades correctly: `{None}` has len 1, passes, and the existing `or ""` at :107 handles it.

MINOR | promote_run.py:96-98 | "so it checks the stronger thing it can see" is backwards. Mutual consistency of the records present is strictly WEAKER than verifying every sink ran: a directory holding one record of a five-sink run passes this check trivially. The sentence tells a reader the promotion gate is stronger than it is. | "so it checks the only thing it can see from here: that the records present are mutually consistent."

MINOR | promote_run.py:90-91 | The comment says records of one run share "the same identity, commit and tree state" and the loop then checks two of the three. `tree_clean` disagreement is mostly caught by the `dirty` check above, but not under `--allow-dirty`, where a mixed directory promotes silently. | Add `tree_clean` to the tuple, or drop it from the sentence.

MINOR | promote_run.py:90-93 | The check does not catch the most likely way a directory acquires two runs' records: `job_out_dir` is created with `mkdir(exist_ok=True)` on a name keyed by `{job.name}_{identity_short}_{timestamp}` at one-second resolution, so two runs of the SAME identity in the same second merge into one directory - and they agree on identity and commit, so this passes. The comment's stated coverage ("records from two runs have been merged") is therefore partial. | Note that same-identity merges are invisible here; the timestamp collision is the real hole.

## Reviewed (slice A)
- [x] scripts/promote_run.py

### runner.py addenda (found while checking targets.py naming)

MINOR | runner.py:144-147 | The namespace enumeration that the whole collision argument rests on is incomplete and one entry is wrong. Checked against `plots/targets.py`: `render_static` -> `{name}_static.pdf`, `render_academic` -> `{name}_academic.pdf`, `render_poster` -> `{name}_poster.png` (NOT `.pdf`), `render_interactive` -> `{name}.html` with NO target suffix. So the docstring's "a figure renders `{name}_{target}.pdf`" misses `.png` and misses the one filename that collides between two same-named figures regardless of their target lists. | "a figure renders `{name}_{target}.pdf`/`.png` and `{name}.html`".

MINOR | runner.py:149-151 | "Keying a figure on its SOURCE NODE instead cannot detect this: ... the mismatch test is `x != x` and the check is inert for every figure." This is a description of the implementation being replaced. It is phrased as an alternative design rather than as history, which is the defensible form, but it is the third comment in this slice whose subject is the previous code. | Keep; if the rule is applied strictly, compress to "keying a figure on its source node makes the mismatch test `x != x`".

### tests/test_reuse_completeness.py (NEW)

CRITICAL | tests/test_reuse_completeness.py:159 | The test is FLAKY and fails with a message that asserts the opposite of the truth. It asserts `(run_dirs[0] / "two.pkl").is_file()` - i.e. that the re-run wrote back into the FIRST run's directory - which only holds if both `run_job` calls land in the same wall-clock SECOND, because `job_out_dir` is `f"{job.name}_{identity_short}_{timestamp}"` at one-second resolution (runner.py:513). Cross a second boundary and the re-run correctly creates a SECOND directory, correctly produces `two.pkl` there, and the assertion fails saying "the partial run was treated as reusable and the job was skipped" - blaming the gate for behaving right. REPRODUCED deterministically: sleeping to the next second between the two calls gives dirs `[twosink_e29ff0_20260902_141126, twosink_e29ff0_20260902_141127]`, `run_dirs[0]/two.pkl` -> False, `two.pkl` present at `twosink_e29ff0_.._141127/two.pkl`. MEASURED window: t0->t2 is 15.2 ms mean over 20 iterations, so ~1.5% failure rate per run on this machine, higher on a loaded CI runner. The comment at :156-158 identifies the same-second behaviour and then relies on it instead of defending against it. | Assert on the artifact anywhere in the pool: `assert list(out.glob("twosink_*/two.pkl"))`. That is still a positive control (with the filter removed nothing is produced at all) and is timestamp-independent.

IMPORTANT | tests/test_reuse_completeness.py:1-18, 139 | The module docstring is written almost entirely as project chronology, the category the author's hard rule bans, and it is the largest single block of it in the slice. "Both of these guard defects that survived because nothing exercised them"; "It read `Path(__file__).parent` while its partner `get_git_commit` read `Path.cwd()`"; "It happened to work in this checkout because..."; "The standalone reuse gate globbed on directory name alone"; ":139 Before the filter this was indistinguishable from a finished run." All past tense about the implementation being replaced. | Restate as the invariants under test: both git helpers must describe the repository the run is happening in; a candidate run must hold every pkl a finished run writes. The "Oracle for both" paragraph at :16-17 is already in the right form and is the model.

IMPORTANT | tests/test_reuse_completeness.py:129-131 | The only test of `_expected_sink_pkls` uses a job with TWO MATERIALIZE SINKS AND NO FIGURE, so the figure branch - the one with the documented subtlety, and the only one whose correctness is non-obvious - is untested. `_expected_sink_pkls`'s docstring makes a specific claim ("a figure sink persists its input under the SOURCE NODE's id whether or not a PDF is rendered") that no test in the repo exercises. Regression direction is safe-but-silent: if the figure branch drifted to `sink.name`, the expected pkl would never exist and the job would re-run on every invocation forever, with no test failing and no message. | Add a figure sink to `_two_sink_job` (or a third job) and assert `{node_id}.pkl` is expected, under both `render_figures` values.

MINOR | tests/test_reuse_completeness.py:43-54 | `tiny_repo` is exposed to the developer's GLOBAL git config in ways that will fail a clean CI runner, and the fixture's own `_git` helper (:33-40) has `capture_output=True` and NO timeout, so a config that prompts hangs the suite silently rather than failing. `user.email`/`user.name` are correctly set locally, so those are covered - the uncovered ones are `commit.gpgsign=true` (fails, or blocks on a passphrase agent), `core.hooksPath` / `init.templateDir` pointing at hooks that run on commit. Verified none are set on THIS machine (`git config --get` on all three returns nothing), so it passes here and will pass on a bare runner; the risk is a contributor with a signing setup. `git init -q` without `--initial-branch` is NOT a portability risk: git 2.51 prints the `init.defaultBranch` hint to stderr, which is captured, exits 0, and no test reads the branch name. | `_git(repo, "-c", "commit.gpgsign=false", "-c", "core.hooksPath=/dev/null", "commit", "-qm", "initial")`, and a `timeout=` on the fixture's `_git`.

MINOR | tests/test_reuse_completeness.py:179 | `assert "Skipping twosink" in capsys.readouterr().out` couples the negative control to a `print` string in `run_job` (runner.py:506). Rewording that message turns this into a failure that reads "a complete run was not reused; the completeness filter is too strict" when nothing about the filter changed. Acceptable given the same-second constraint the docstring explains, but a structural oracle exists. | Prefer asserting `_reuse_eligible_dir(...) is not None`, or that the run-dir mtime did not change; or at minimum import the message from one place.

VERIFIED | tests/test_reuse_completeness.py | 7 passed in 0.82 s. `test_outside_a_repository_reports_nogit_and_dirty` and `test_a_dirty_tracked_file_makes_the_tree_dirty` both fail if `_git` stops returning None on error or if `is_tree_clean` stops treating None as dirty, so those are real. `test_both_git_helpers_describe_the_cwd_repository` fails if either helper is re-anchored on `__file__` (a temp repo is not the package's repo), so it is a genuine pin on the change.

## Reviewed (slice A)
- [x] tests/test_reuse_completeness.py

### tests/test_identity_closure.py

VERIFIED | test_identity_closure.py:1-11, 148 | Correct on both counts: the rename reaches the only two call sites, and the docstring is the ONE place in this slice where the chronology rule was applied properly - "SPEC 0005 R5.1. Before this, `code_hash` was..." became the counterfactual "`job_code_hash` over the job file alone WOULD leave...", and "the blast radius the plan warned about" became "an unacceptable blast radius". No spec number, no phase reference, no past tense. Nothing else in the file references the old name.

## Reviewed (slice A)
- [x] tests/test_identity_closure.py

### tests/test_shape_stats.py

IMPORTANT | tests/test_shape_stats.py:1-8 | The docstring's justification for the file's existence is now FALSE in the shipped state, and this diff rewrote the sentence rather than fixing it. It reads "with `Job.job_code_hash` covering only the job file, an analyzer edit left the run identity byte-identical ... Without this file the shape statistics could be changed with no signal anywhere in the repo." MEASURED: `quebra.analyzers.shape_stats` IS in the code closure of the jobs that use it - `code_closure` for `jobs/composite/independence_survey.py` returns 52 modules including `quebra.analyzers.shape_stats`, and for `jobs/active/instrument_validation.py` 23 modules, also including it. So editing that analyzer now DOES move `job_code_hash`, the identity, and the run directory name: there IS a signal, and it is the very mechanism the commit next door is about. The docstring also keeps the past-tense chronology the hard rule bans ("Written BEFORE...", "at the time", "left", "kept loading"). | State what is actually load-bearing now: the identity says the numbers CHANGED, it cannot say they are RIGHT; these tests are the only thing pinning the eq (8) values themselves. And drop the tense.

## Reviewed (slice A)
- [x] tests/test_shape_stats.py

### Cross-cutting / omissions

IMPORTANT | omission | The `job_file_hash` -> `job_code_hash` rename is not complete across the tree: `spec/specidentity05.md:39` and `:282` still name `Job.code_hash`, and that file is NOT among the modified files, so the rename leaves two dangling references in a document the repo treats as the authority for this mechanism. (Called out only as rename completeness - `spec/*.md` is excluded from this slice's review.) | Rename in the spec in the same commit, or note the alias.

VERIFIED | omission check | `_expected_sink_pkls`'s guarantee for the composite/`--reuse-deps` path was checked and needs nothing: nested runs are written to `subjobs_dir` and always with `force=True`, and `_cached_runs` (runner.py:263-292) applies its own artifact-existence filter. No third call site of the gate exists.

VERIFIED | measured-claim audit | Every measured claim in the added comments was checked. The only false ones are provenance.py:28-29 (credential prompt), closure.py:161-172 ("located rather than imported"), closure.py:68-69 ("tens of milliseconds", unsourced), runner.py:145 (`{name}_{target}.pdf` only), runner.py:157-158 ("harmless"), promote_run.py:97 ("stronger"), and test_shape_stats.py:1-8. All filed above.

MINOR | evidence for the "nogit" finding above | The exposure the commit+clean-tree gate is the ONLY defence for is real and measured: `jobs/active/instrument_validation.py` declares 4 figures whose plot classes all live in `quebra.plots.instrument_validation_plot`, and that module is NOT in the job's 23-module code closure (checked directly: `in closure? False` for all four). `job.figure(PlotClass, ...)` passes a class, not a step function, so the plot module is never seeded. Editing it changes every shipped PDF with the identity byte-identical. That is why "nogit" == "nogit" defeating the commit half of the gate matters.

VERIFIED | gates | `pytest tests/test_promote_run.py tests/test_reuse_gate.py tests/test_identity_closure.py tests/test_nesting.py -q` -> 45 passed. In particular no existing promote fixture writes a run directory with disagreeing `identity`/`git_commit`, so the new promote check does not break the existing suite. `pytest tests/test_reuse_completeness.py -q` -> 7 passed.

SLICE A VERDICT: DO NOT SHIP (2 blockers: runner.py:141-158 prov-record overwrite exempted by the same-node carve-out, MEASURED on a shipped run; test_reuse_completeness.py:159 flaky assertion, reproduced).

--- PATCH COMMIT 3 of 4 — BEHAVIOUR, appended 2026-09-02 ---

SCOPE: src/quebra/transforms/interpolate.py, src/quebra/transforms/filter.py,
src/quebra/analyzers/within_calibration_compute.py, src/quebra/analyzers/checks/_multiprocess.py,
src/quebra/analyzers/tlf.py, src/quebra/analyzers/instrument_validation.py,
src/quebra/plots/tlf_plot.py, src/quebra/recipes.py, jobs/active/instrument_validation.py,
tests/test_tier3_calibration.py, tests/test_instrument_validation.py, docs/JOBS.md (generated check only)
COMMIT: c0056f4 (working tree, uncommitted)

## Manifest (patch 3)
- [ ] src/quebra/transforms/interpolate.py  (+133/-?)
- [ ] src/quebra/transforms/filter.py  (+70)
- [ ] src/quebra/analyzers/within_calibration_compute.py  (+26)
- [ ] src/quebra/analyzers/checks/_multiprocess.py  (+19)
- [ ] src/quebra/analyzers/tlf.py  (+44)
- [ ] src/quebra/recipes.py  (+12)
- [ ] src/quebra/plots/tlf_plot.py  (+51)
- [ ] src/quebra/analyzers/instrument_validation.py  (+49)
- [ ] jobs/active/instrument_validation.py  (+4)
- [ ] tests/test_tier3_calibration.py  (+56)
- [ ] tests/test_instrument_validation.py  (+4)
- [ ] docs/JOBS.md  (generated)

## Reviewed (patch 3)

## Findings (patch 3)

--- PATCH COMMIT 3 of 4 — PROSE (comments/docstrings only), appended 2026-09-02 ---

SCOPE: prose (comments, docstrings, markdown) in the added/changed lines of:
src/quebra/transforms/interpolate.py, src/quebra/transforms/filter.py,
src/quebra/analyzers/within_calibration_compute.py, src/quebra/analyzers/checks/_multiprocess.py,
src/quebra/analyzers/tlf.py, src/quebra/analyzers/instrument_validation.py,
src/quebra/plots/tlf_plot.py, src/quebra/recipes.py, jobs/active/instrument_validation.py,
tests/test_tier3_calibration.py, tests/test_instrument_validation.py, docs/JOBS.md
COMMIT: c0056f4 (working tree, uncommitted)
NOTE: working tree is being edited concurrently - interpolate.py gained a paragraph between
my `git diff` and my `sed` of the same file. All findings below are against the CURRENT file
contents (what would ship), not against the diff text I first read.

## Manifest (prose)
- [ ] src/quebra/transforms/interpolate.py
- [x] src/quebra/transforms/filter.py MOVED-TO-REVIEWED
- [ ] src/quebra/analyzers/within_calibration_compute.py
- [ ] src/quebra/analyzers/checks/_multiprocess.py
- [ ] src/quebra/analyzers/tlf.py
- [ ] src/quebra/analyzers/instrument_validation.py
- [ ] src/quebra/plots/tlf_plot.py
- [ ] src/quebra/recipes.py
- [ ] jobs/active/instrument_validation.py
- [ ] tests/test_tier3_calibration.py
- [ ] tests/test_instrument_validation.py
- [ ] docs/JOBS.md

## Reviewed (prose)
- [x] src/quebra/transforms/interpolate.py

## Findings (prose)

### src/quebra/transforms/interpolate.py

MUST FIX | interpolate.py:83-85 | `run` docstring invariant 2 - "no value on the grid was invented to stand in for one the instrument did not produce" - is FALSE as written, and the module's own meta contradicts it: `n_interpolated_points` / `pct_interpolated_points` count exactly the grid values that stand in for reads the instrument did not produce at that time. Interpolation invents values by design; the distinction meant is fabricated-zero vs interpolated. Also "enforced rather than assumed" overclaims - nothing asserts either invariant on the way out, they hold by construction. | Replace all three lines with: "Row-aligned columns come out at the grid length. Reads whose value did not fit are dropped before the pchip, never zero-filled."

SHOULD FIX | interpolate.py:54-55 | "Dropping the pair instead bridges the gap from the real neighbours either side, which is interpolation doing what it says." - contradicted by the paragraph directly under it (a dropped leading/trailing value has no neighbour either side, and the code NaNs those grid points), and "which is interpolation doing what it says" is rhetorical filler. Para 3 subsumes it. | Delete both lines.

SHOULD FIX | interpolate.py:45-60 | 15-line, four-paragraph docstring on a 15-line function; the mechanism is stated in the summary and again twice in the body. Caps-for-emphasis three times (FINITE, THROUGH, plus RAISE at :197). | Keep the summary, one sentence of para 1, and para 3 trimmed: `"""pchip through the finite (t, value) pairs only, plus the count dropped.` / `Zero-filling instead would put a real, wrong measurement on the grid - for a detuning, exactly 0 Hz - and pchip passes through its data, so it would drag the neighbours too.` / `Grid points outside the finite support are NaN, not extrapolated: dropping a leading or trailing value shortens the support while the grid still spans the original interval."""` Drops ~7 lines.

SHOULD FIX | interpolate.py:114-115 | "a non-finite timestamp sorts last, so it becomes `t_rel_s[-1]`" is true only for NaN; `-inf` sorts FIRST under `np.argsort`, becomes `t_rel_s[0]`, and breaks the grid through the `t_rel_s - t_rel_s[0]` shift instead. Claim FALSE for half the values the guard rejects. Sentence 1 also duplicates the raise message two lines below ("The time axis has to be real before anything can be resampled onto it"). | Replace both lines with: `# A NaN timestamp sorts last and becomes t_rel_s[-1]; an infinite one survives the shift. Either way the whole uniform grid is unusable, and the time axis has no interpolation to fall back on.`

SHOULD FIX | interpolate.py:197-200 | "would leave it at the input length ... nothing downstream can detect because the lengths still look plausible" - FALSE. `x_uniform = np.linspace(..., num=len(t_rel_s))`, so the grid length EQUALS the input length by construction; there is no length discrepancy to look plausible, only a sampling one. "failure must RAISE" also describes the removed blanket `except`, not the code present. | Replace with: `# Row-aligned from here. Passing the array through unresampled would keep the input sampling at the grid length - a column misaligned against its own time axis, at a length nothing downstream can use to detect it.`

OPTIONAL | interpolate.py:191-192 | "scalars, lookup tables, per-dataset annotations" is a rhythm triad of speculative examples, and "Passed through untouched" restates the two lines of code under it. | Keep only: `# Not row-aligned with the time axis, so there is nothing to resample.`

OPTIONAL | interpolate.py:186-187 | "There is no time axis to resample it against" is the wrong reason - the time axis exists; what is missing is an array. | Keep only: `# Not array-shaped at all - ragged, or an object numpy cannot box.`

OPTIONAL | interpolate.py:217-219 | "Recorded rather than printed: a step is pure compute, and this belongs in the artifact a reader audits, not in a log line nobody keeps" narrates the removal of the `print` and restates a project rule; only the absent-key semantics is load-bearing. | Keep: `# Which columns lost reads to failed fits, and how many. Absent, nothing was dropped.`

OPTIONAL | interpolate.py:131-133 | Three lines where one does; "matching `rabi_hz`" is visible in the branch immediately above. | Keep: `# Raises rather than dropping the column: a missing key reads downstream as "this dataset has no raw frequency", not "its length was wrong".`

MINOR | interpolate.py:48 vs :217 | Vocabulary drift inside new prose: "at that read" (:48) vs "lost points" (:217) for the same thing. `read` is on the fixed list. | Say "read" in both.

VERIFIED | interpolate.py:51-52 | "Downstream finite-mask guards cannot reject any of it" - the guards exist (windows.py:371, shape_stats.py:116, t2star.py:81, filter.py:193). Claim true.
VERIFIED | interpolate.py:57-60 | NaN-outside-support paragraph matches the `outside` mask in the code, and pchip_interpolate does extrapolate by default.

### src/quebra/transforms/filter.py  [reviewed]

MUST FIX | filter.py:11-13 | "Named explicitly because the length test alone cannot tell them apart: a lookup table can coincidentally match the mask length, and a row-aligned column can coincidentally not." - the first hazard is NOT addressed by the code. The set is consulted only in the `elif`; a lookup table whose length happens to match the mask is still filtered, in the set or out of it. The comment credits the frozenset with a guarantee it does not provide. | Cut to two lines: `# Columns that are one value per read by definition. Consulted only to turn a length that` / `# disagrees with the mask into a raise instead of a silent pass-through.`

SHOULD FIX | filter.py:191-192 | `raises "cannot run on an empty dataset"` is not a verbatim string anywhere in the repo. The real strings are `"Cannot run filter on empty dataset."` (filter.py:115) and `"Cannot run T2* analysis on empty dataset."` (t2star.py:77); the comment quotes a composite of the two. A quoted error message must match the artifact it cites. | Drop the quotation: `...the mask all-False, and the dataset empty by the time an analyzer complains about something else.`

SHOULD FIX | filter.py:206-209 | Four lines explaining IEEE NaN comparison rather than the code. The load-bearing content is that the term is redundant now and defends a future rewrite. | Keep two: `# Redundant today - each NaN comparison is already False - but it states the intent, and it` / `# survives the bounds being rewritten as a negation, which would invert the NaN case.`

SHOULD FIX | filter.py:44-47 | "nothing detects it: the lengths still look plausible and the values are all real" - the stated reason is wrong for this function. If a row-aligned column is left unfiltered while `t_rel_s` is masked, its length visibly DIFFERS from `t_rel_s`; the hazard is that no consumer raises on it (interpolate.py:191 silently treats a length mismatch as not-row-aligned), not that the lengths look plausible. | `A row-aligned column left at full length is index-misaligned against `t_rel_s` for every step downstream, and no consumer raises on it - interpolate treats a length mismatch as simply not row-aligned.`

OPTIONAL | filter.py:49-50 | "scalars, lookup tables, per-dataset annotations" is the same speculative rhythm triad used at interpolate.py:191; "There is nothing to filter in them" restates the clause before it. | `Columns that are not row-aligned are left untouched.`

OPTIONAL | filter.py:42-51 | 10-line docstring on a 25-line function that the summary line already describes accurately. With the two trims above it is summary + 2 lines. | Net delete ~5 lines.

MINOR | interpolate.py:217-218 vs filter.py:118,222-225 | The justification written in interpolate.py - "Recorded rather than printed: a step is pure compute" - is contradicted by its sibling in the same commit: filter.py, also a transform step, keeps two `print` calls. The asserted rule is not true of the shipped state. | Either drop the rule from the interpolate comment or the prints are a code finding for the correctness reviewer.

VERIFIED | filter.py:202-203 | "Every finite value identical: no outlier is definable, so keep them and drop only the non-finite reads" matches `mask = finite`. Short, correct, keep as is - the best comment in the file.
VERIFIED | filter.py:189-190 | np.mean/np.std do propagate NaN; filter's own empty guard is on input only (:113), so an all-False mask does return an empty Norm. Substance of the claim is right; only the quoted string is wrong.

### interpolate.py

VERIFIED | interpolate.py:39-75 | Change 1, the numerics. pchip's node derivative is Fritsch-Carlson, a function of the two adjacent intervals only, so dropping one node perturbs the interpolant on at most the 4 intervals around it and nothing further. Measured: 20 nodes, drop index 10, max |full-reduced| away from index 10 is exactly 0.0. Interpolating from reduced support is the right operation and does NOT move values far from the dropped read.

VERIFIED | interpolate.py:71-74 | Change 1, ends. First/last non-finite does NOT extrapolate: `outside` NaNs the grid beyond `t_support`. Measured: delta_hz with index 0 and -1 NaN returns `[nan, 1..8, nan]`.

IMPORTANT | interpolate.py:209-214 | `_real_vs_interpolated_counts(t_rel_s, out["t_rel_s"])` is still passed the FULL input time axis, not the per-column finite support, so `pct_real_points` counts a grid point as "real" when the read at that timestamp was dropped as non-finite and its grid value was in fact interpolated - or is NaN. Measured: 10 reads, delta_hz[4]=NaN -> `pct_real_points=100.0`, `n_nonfinite_dropped={'delta_hz':1}`; both ends NaN -> `pct_real_points=100.0` with two NaN grid points that are neither real nor interpolated. `plots/interpolation_stage_plot.py:132` annotates that number onto a published figure. | Compute the counts against the delta_hz support (`t_rel_s[np.isfinite(delta_hz_sorted)]`), or subtract the dropped count from `n_real`.

MINOR | interpolate.py:69-70 | Change 2 regression: a row-aligned 2-D numeric column now raises `IndexError: too many indices for array` (`t_rel_s[finite]` with a 2-D `finite`). The old code resampled it correctly - `pchip_interpolate` broadcasts y along axis 0. Measured with a (6,2) column. Not reachable today (NullSchema refuses duplicate labels precisely to keep everything 1-D, RamseySeriesSchema emits only 1-D), but the raise the diff advertises is a ValueError and this is an unhandled IndexError from a private helper. | Guard `values.ndim != 1` in `_interpolate_finite` with a message, or mask along axis 0.

MINOR | interpolate.py:64-68 | Change 2, reachability of the new raise. The only row-aligned columns any shipped schema emits are `chi_squared`, `T2star_s`, `T2star_error_s` (ramsey_series.py:136-146, built with `errors="coerce"`, so NaN-capable). A dataset where one of those columns is present but entirely unparseable used to yield an all-zero resampled column; it now raises `cannot interpolate 'T2star_error_s': 0 finite value(s)`. For `track912Schema` the `dropna(subset=...)` at track912.py:79 makes it unreachable; for a bare `RamseySeriesSchema` load it is reachable. The raise is the right call, but this is a load path that used to succeed and now does not. | None needed; note it in the commit message.

MINOR | interpolate.py:128-134 | Change 3, `raw_frequency_hz` length raise. Nothing relied on the drop: the only producer is `ramsey_series.py:130`, which sets it to `f_hz`, the same masked-and-sorted array as `t_rel_s`, so the lengths cannot disagree. The old `len(...) == len(order)` was dead in every shipped path. `filter.py:15` lists it in `_ROW_ALIGNED_KEYS` and `fidelity.py:102` re-checks the length itself, so no consumer wanted the silent drop. Guard is correct but currently unfirable.

MINOR | interpolate.py:112-119 | Change 3, non-finite `t_rel_s` raise. Also unreachable from any shipped path: `ramsey_series.py:52-55` masks on `np.isfinite(timestamp) & np.isfinite(frequency)` before building `t_rel_s`. Correct as a contract check; it is a claim, not a live check.

MINOR | interpolate.py (whole file) | None of change 1/2/3 has a test. `tests/` has no unit test for `interpolate.run` at all (only `test_windows_not_interpolated.py`, which asserts DAG shape, and `test_load_dataset_contract.py`). The NaN-drop, the NaN-edge and `meta["n_nonfinite_dropped"]` are all unpinned, so a future revert to `nan_to_num` passes `make check`. | Add three asserts: dropped interior read bridges from neighbours, leading/trailing NaN stays NaN, `n_nonfinite_dropped` counts.

VERIFIED | interpolate.py:215-221 | Removing the `print` is the everything-is-a-step rule applied correctly (no I/O in a step) and the information is not lost - it moves to `meta`. `filter.py:144,153` still print, which is a pre-existing violation this diff did not touch.

## Reviewed (patch 3)
- [x] src/quebra/transforms/interpolate.py

### filter.py

IMPORTANT | filter.py:14-16 | Change 4: `_ROW_ALIGNED_KEYS` is INCOMPLETE and carries one dead entry. Missing, all one-value-per-read and all emitted by shipped code: `chi_squared`, `T2star_s`, `T2star_error_s` (ramsey_series.py:136-146) and `qubit_frequency_hz` (lookup_prior.py:133-135, injected by jobs/active/ramsey_q1_100423.py:48 in the same call that injects `rabi_hz`, which IS in the set). Dead: `t_unix_s` is not a Norm key any producer in src/ or jobs/ writes - it exists only as a field on `types.Measurement`. So the guard omits exactly the columns that arrive through the NaN-capable `errors="coerce"` path and includes one that cannot occur. | Add the four; drop `t_unix_s` or point it at a producer.

IMPORTANT | filter.py:10-13,67-76 | Change 4: the stated reason for a named set over the length test - "a lookup table can coincidentally match the mask length" - is not what the code does. Line 67 still masks ANY column whose length equals the mask, named or not, so a coincidentally-matching lookup table is still silently sliced. The set only converts silent-passthrough into a raise for the five named keys; it buys nothing for the first hazard. Behaviourally that is fine, but the guard is half of what it claims, and the half it does not do is the one that corrupts values rather than merely misaligning them. | Either drop the first clause of the justification, or add the reciprocal check (a NON-row-aligned key whose length coincidentally matches must not be masked - which needs a positive declaration, i.e. the set inverted).

MINOR | filter.py:57-64 | New raise: `np.asarray(value)` failure is now fatal where it used to `continue`. `interpolate.py:181-187`, changed in the same diff, takes the OPPOSITE decision for the identical input class (ragged/unboxable -> pass through like a scalar). Since filter runs before interpolate in `configure_ramsey_job`, interpolate's branch is dead in every shipped DAG. Two transforms disagreeing on the same input is a contract split. | Pick one. Raising in both is the stricter and more defensible half.

MINOR | filter.py:194-199 | Change 5: the new all-non-finite raise also fires for an EMPTY `delta_hz` (`np.any(np.zeros(0,bool))` is False), which is reachable when the chi stage removes every read (threshold 0.8 on normalised chi-square). The message then reads "0 read(s), all non-finite", which mis-describes an empty stage. Behaviour is still an improvement over the old path (NaN mu/sigma -> all-False mask -> a downstream "cannot run on an empty dataset"), but the diagnosis is wrong. | `if delta_hz.size == 0: raise ValueError("the sigma stage received 0 reads; an earlier stage removed them all")` before the finite test.

VERIFIED | filter.py:200-203 | Change 5, `sigma == 0.0` -> `finite`. Not a silent row-count change on any real dataset: `sigma == 0.0` requires every finite `delta_hz` identical, and `delta_hz` is built at ramsey_series.py:126 as `f_hz - mean(f_hz)` from an array already masked to finite (ramsey_series.py:52-55), so `finite` is all-True and `mask = finite` is identical to the old `np.ones(...)`. The only other reachable case is a single surviving read, where `np.std` is 0.0 and both old and new keep it. The change is correct and inert.

MINOR | filter.py:186-213 | Change 5 as a whole is unfirable on shipped data for the same reason: `delta_hz` cannot be non-finite downstream of `RamseySeriesSchema`. The `finite &` in the mask, the all-non-finite raise, and the `sigma == 0.0` branch are all contract checks against a norm no in-repo schema produces. Correct, but claims rather than live checks.

MINOR | filter.py (whole file) | `tests/` contains no test that imports `quebra.transforms.filter`. Every change 4 and change 5 behaviour is unpinned. | Two tests: a row-aligned key at the wrong length raises; a non-finite `delta_hz` is dropped by the sigma stage rather than emptying the mask.

## Reviewed (patch 3)
- [x] src/quebra/transforms/filter.py

### src/quebra/analyzers/within_calibration_compute.py  [reviewed]

SHOULD FIX | within_calibration_compute.py:219-223 | 12-line docstring on a 3-line function, and the first paragraph explains why a bug was a bug (an `ax.plot` crash) instead of stating the invariant (results are returned at row length so they align with `t_h`). "which is a rendering crash rather than a wrong number, but only because matplotlib happens to check" is editorial and adds nothing a caller can act on. This is also a compute module justifying its return shape by matplotlib's internal error text - the layer boundary the project keeps in the code should hold in the prose. | Keep summary + one sentence: `"""Put a finite-subset result back on the full row index, NaN where rows were dropped.` / `Callers hold the unmasked row index; a subset-length return would be silently misaligned against it. NaN rather than 0.0 because these are cumulative quantities - a zero would drop the curve back to the origin at each missing read."""` Deletes 6 lines.

MINOR | within_calibration_compute.py:220-221 | Quoted matplotlib text `x and y must have same first dimension` is a truncation of the real message (which continues ", but have shapes ..."). If the quote stays anywhere, it should not look verbatim. | Drop the quote with the paragraph.

MINOR | within_calibration_compute.py:220 | Caps-for-emphasis (`UNMASKED`); same tic as interpolate.py (FINITE/THROUGH/RAISE) and _multiprocess.py (LAST). | Lower-case it; the sentence carries the emphasis.

VERIFIED | within_calibration_compute.py:225-227 | NaN-not-zero rationale is correct for cumulative series, and matplotlib does break a line at NaN.

### src/quebra/analyzers/checks/_multiprocess.py  [reviewed]

SHOULD FIX | _multiprocess.py:296-306 | 11 lines of comment on a 7-line guard, and lines 303-305 restate the raise message immediately below them ("without gap_spans_s a gap flanked by out-of-spec reads leaves a censored death mid-block, and every event time after it would be wrong"). | Cut to the invariant, the mechanism, the reachability, and the authority - about 7 lines: `# Only the LAST window of a block may be censored: blocks are split at read gaps, so a` / `# gap_start death lands at a block end. x drops a censored window's duration while tau` / `# keeps it, so an interior one shifts every later T_i earlier - eq (4) and eq (7) both run` / `# on the wrong ones and report the censored duration as a residual that is really zero.` / `# Reachable: check_ledger.make_inputs_from_windows takes gap_spans_s from` / `# diagnostics.get(...), so a caller that omits it splits on birth types alone. Segment` / `# documents n_censored_dropped as 0 or 1 and nothing enforces it.`

SHOULD FIX | _multiprocess.py:305-306 | "Raising is right rather than defensive" argues with an imagined reader who has not objected yet. The `Segment` citation that follows is the actual justification and stands without it. | Delete the clause; keep the citation.

MINOR | _multiprocess.py:298 | "it is what makes the event times below correct" - the event times are not below; they are formed in c1_lewis_robinson and c2_anderson_darling from `Segment.x`. | Say "the T_i that eq (4) and eq (7) form from `x`".

VERIFIED | _multiprocess.py:296-302 | `x = durations[complete]` vs `tau = float(durations.sum())` - x drops, tau keeps. Blocks are split via `_segment_starts(births, t_birth_all, gap_spans_s)`. eq (4) = C1 Lewis-Robinson, eq (7) = C2 Anderson-Darling, both consuming `Segment.x`. All as claimed.
VERIFIED | _multiprocess.py:303-304 | `check_ledger.make_inputs_from_windows` does set `gap_spans_s=diagnostics.get("gap_spans_s")` (check_ledger.py:196), so None is reachable. Module lives at analyzers/check_ledger.py, not analyzers/checks/ - the unqualified name in the comment is still unambiguous.
VERIFIED | _multiprocess.py:305-306 | `Segment.n_censored_dropped` is documented "0 or 1 per segment" (result.py:88-90) and no code anywhere asserts it. "nothing else enforces it" is true.

### within_calibration_compute.py

VERIFIED | within_calibration_compute.py:216-231 | Change 6 is a real crash fix and it is correctly wired. `build_within_calibration_panel_data` (line 464, 486) hands the SAME `t_arr` to `signal_band.run` and `reliability_band.run`, and `signal_band.run` stores it unfiltered (signal_band.py:114,171), so `pd_.signal.t_h` really is the full row index and `_to_full_length(..., mask)` really does match it. Before: a single non-finite `t_h` or `primary_series` value made `cumulative_time_per_threshold[label]` shorter than `signal.t_h` and `within_calibration.py:899-906` / `:949-956` raised `x and y must have same first dimension`.

VERIFIED | within_calibration_compute.py:248-250,290-292 | The `len(t_f) < 2` short-circuit is correct. `mask.sum() == len(t_f)` by construction, so `full[mask] = np.zeros(len(t_f))` is shape-consistent for len(t_f) in {0,1}; the result is length `len(mask)`, zero at the surviving row and NaN elsewhere, which is what the panel needs.

VERIFIED | within_calibration_compute.py | Change 6 is COMPLETE for the two arrays. The other finite-masked helpers in the file (`_ttf`:317, `_threshold_in_spec_frac`:345, `_threshold_summary`:377) all return scalars or None, never a row-indexed array, so none of them needed the same padding. Every reader of the two arrays was checked: `reliability_band.py:224-225` (assignment), `reliability_band.py:100-101` (`to_dict`, not written to any file under `jobs/bench/results/`), `panels/within_calibration.py:897,948` (`ax.plot` only, no sum/max/index). Nothing sums, maxes or indexes them, so the NaN cannot propagate into a number.

MINOR | tests/test_within_calibration_builder.py:78-85 | Change 6 is untested. `test_cumulative_time_is_monotonic_and_bounded` already asserted `len(arr) == len(d.signal.t_h)` but runs on a dense all-finite synthetic series, so it passed before the fix and passes after; it never enters the padding path. Worse, its own assertions are NaN-hostile (`np.all(np.diff(arr) >= -1e-12)` and `arr[-1] <= total_h` are both False once a NaN exists), so the obvious extension of it would fail. | Add a case with one interior NaN in `primary_series` asserting `len(arr) == len(t_h)` and `np.isnan(arr[i])`, using `np.nanmax`/`np.diff` on the finite subset.

## Reviewed (patch 3)
- [x] src/quebra/analyzers/within_calibration_compute.py

### src/quebra/analyzers/tlf.py  [reviewed]

MUST FIX | tlf.py:48-50 | FALSE claim about another module: "the seed arrives here the way `xi_seed` does, as a step argument, so it reaches both the identity and the provenance label." It does not. `xi_seed` is a real step kwarg (`recipes.py:282: xi_seed=xi_seed`), which `closure.py:259-262` folds into the identity and `provenance.py:130` renders into the step label. `tlf_seed` is captured in a factory closure - `recipes.py:293: job.step(_tlf_step(tlf_seed), final_filtered, name="tlf")` passes no kwargs - so it appears in NEITHER the kwarg rows of the identity nor the step's provenance label. | Either delete the sentence or state what is true: `Passed to the step factory, so its value lives in the recipe text the code hash covers.` The structural fix is a code finding for the correctness reviewer.

MUST FIX | tlf.py:34-36 | "Without this, a failure returned bic_delta=0.0 with `gmm2` aliased to `gmm1`, so a consumer read ..." is a description of the code that was here, in the past tense, which rule 1 forbids. It explains why a bug was a bug instead of stating the invariant. | `# True when the 2-component fit raised. Consumers must read it before is_bimodal: False there means no comparison was made, not that one lobe was established.`

SHOULD FIX | tlf.py:15 | `"the two-lobe fit did not converge"` names the wrong failure mode. `fit_failed` is set by `except Exception` around `gmm2.fit`, i.e. the fit RAISED; `GaussianMixture` does not raise on non-convergence, it warns and returns. | `"the two-lobe fit did not run"`.

SHOULD FIX | tlf.py:13-15 | Three lines where sentence 2 and sentence 3 both restate sentence 1, and the whole thing is repeated at :34-36. | Keep: `# None when the 2-component fit raised. A number here always means the comparison happened.`

SHOULD FIX | tlf.py:43-50 | Eight-line paragraph for one fact plus one false claim; the same mechanism is written a second time at recipes.py:183-185. | Keep three lines: `` `seed` is required, not defaulted. `GaussianMixture` initialises by k-means, so an unseeded fit makes `bic_delta`, `is_bimodal`, the state assignment and every dwell statistic a fresh random variable per call while the run identity stays byte-identical. `` Drops 5 lines and the false sentence with them.

OPTIONAL | tlf.py:76-77 | "`fit_failed` is the field that says which" is the third statement of the same thing in one file (:13, :34, :76). | Keep only the first sentence: `# False because nothing established two lobes - not because one lobe was established.`

VERIFIED | tlf.py:43-45 | sklearn `GaussianMixture` does default to `init_params="kmeans"` and `random_state=None` draws on the global numpy RNG. Mechanism claim is right.
VERIFIED | tlf.py:47 | "The same defect is refused outright in `checks/_permutation`" - _permutation.py:92-93 does raise on an unseeded generator with the identity argument. True, but it is authority-by-analogy that the reader does not need; goes with the trim above.

### src/quebra/recipes.py  [reviewed]

MUST FIX | recipes.py:184-185 | "Passed as a step argument so it folds into that identity." FALSE for the same reason as tlf.py:48 - `TLF_SEED` reaches the step through a closure over `_tlf_step(seed)`, not through `job.step(..., tlf_seed=...)`, so `closure.py`'s kwarg rows never see it. | `# Its value lives in this file's text, which the code hash covers.` (accurate today) - or make it a real kwarg and keep the sentence.

SHOULD FIX | recipes.py:183-185 | The k-means mechanism is now written in full in two places (here and tlf.py:43-46). One of them should be a pointer. | Keep the one at the parameter (`tlf.run`); reduce here to `# Seeded because GaussianMixture initialises by k-means. See tlf.run.`

OPTIONAL (adjacent, unchanged line) | recipes.py:180 | "it is a design change, not a rename, and it is Increment C work" - internal chronology of exactly the kind rule 1 bans, sitting three lines above the new comment. Not in this diff, but it will be read as part of it. | Drop "and it is Increment C work".

### src/quebra/plots/tlf_plot.py  [reviewed]

OPTIONAL | tlf_plot.py:86-88 | Three lines; the second and third argue the hypothetical rather than state the rule. | Keep two: `# No GMM(2) curve when the two-component fit failed: a single-lobe curve under a "GMM(2)" label reads as evidence against bimodality rather than as an absent fit.`

VERIFIED | tlf_plot.py:243-253 | The summary block does distinguish absent from computed - "is_bimodal: n/a (2-component fit failed)" and "bic_delta: n/a" - so the downstream consumer really does see the difference the tlf.py comments promise. Claim consistent with what ships.

### checks/_multiprocess.py

VERIFIED | _multiprocess.py:308-309 | Change 7 index arithmetic is right. For a single-window block `complete[:-1]` is empty and `np.any` of an empty bool array is False, so a lone censored window (the normal `scan_end` case) does not raise - measured. For length n, `complete[:-1]` is exactly indices 0..n-2, i.e. every window but the last. An interior `gap_start` death with `gap_spans_s=None` raises; supplying `gap_spans_s=[(5.0, 99.0)]` splits the block and yields `Segment(x=[3.], tau=3.0, n_censored_dropped=0)` - measured both.

VERIFIED | _multiprocess.py:308-316 | The guard cannot fire on any shipped path, and the claim holds for the right reason. `windows.run` unconditionally writes `gap_spans_s` into `diagnostics` (windows.py:494-509, empty list when there are no gaps), and both shipped consumers take it from there (recipes.py:196, recipes.py:380 -> check_ledger.py:200). The only censoring death types are `gap_start` and `scan_end` (windows.py:45-47); `_segment_starts` splits at the first birth at-or-after each `t_after_s`, so a `gap_start` death is always the last window of its block, and `scan_end` only ever belongs to the last window of the record. `n_censored = int((~complete).sum())` is therefore in {0,1}, which is what `Segment` documents. It is a contract check on `CheckLedgerInputs(gap_spans_s=None)` (check_ledger.py:154 defaults it to None), not a live check.

MINOR | _multiprocess.py:308-316 | The same defect corrupts the CALENDAR clock and is not guarded there. With an interior censored window and no gap list, `x = np.diff(t_birth)` includes an inter-birth interval that spans unobserved hours and `tau = t_death[-1] - t_birth[0]` spans them too (`block_end_s` stays None because `gap_starts_s` is empty), so eq (4)/(7) get at-risk time the instrument was not watching - and `n_censored = 1` under-reports the two censored windows. The guard's placement inside `if clock == CLOCK_IN_SPEC` makes it half a check. | Hoist the `(~complete[:-1]).any()` test above the `if clock ==` branch.

MINOR | _multiprocess.py:308 | The guard sits BEFORE the `len(x) < min_events` drop, so a block that the old code silently discarded into `n_dropped` now aborts the whole ledger. Concretely: a two-window block `[censored, censored]` produced `x = []`, was dropped, and `n_dropped` recorded it; it now raises. | Move the test after the `min_events` drop, or accept it and say so - either is defensible, but the current order turns a counted drop into a hard failure.

## Reviewed (patch 3)
- [x] src/quebra/analyzers/checks/_multiprocess.py

### src/quebra/analyzers/instrument_validation.py  [reviewed]

MUST FIX | instrument_validation.py:118 | "which is the whole point" - banned construction, and it is doing no work: the sentence before it already says the count is random. Same phrase is already in the repo at checks/result.py:117 ("it is the whole point of time censoring"), so it is becoming a tic. | `# A truncation time, not an event count. Gaps are unit-mean, so about this many events land, and how many is random.` (2 lines -> 2, but drops the tell; also drop the caps on "A TRUNCATION TIME".)

MUST FIX | instrument_validation.py:142-143 | "`checks/battery` refuses CvM-asymptotic when it is not met" is FALSE. battery.py gates CvM-asymptotic on `include_tau_checks` and `len(segments) == 1` (battery.py:170-186); it cannot detect where tau came from. The only tau refusal in the repo is `validate_segment`'s `tau` vs `T_N` clearance (result.py:43), and `tau = T_N(1 + 1/n)` clears it comfortably - the exact scheme being criticised would sail straight through. Nothing refuses it, which is the actual reason this function has to construct the null itself. | `...chosen WITHOUT reference to the events, which `checks/result` states as a requirement and nothing downstream can check for you.`

SHOULD FIX | instrument_validation.py:147 | "Eq (7)'s tail term depends on that leftover window" is applied to a function that measures C1, C2 and CvM. Eq (7) is C2's statistic; C1 uses eq (4) and CvM's integrand explicitly has no `1/(s(1-s))` weight (cvm_cramer_von_mises.py:177). The mechanism named accounts for one of the three instruments. | `Every one of these statistics reads the window past the last event; frozen, it contributes no variance...` - or name eq (7) as the clearest case rather than the cause.

SHOULD FIX | instrument_validation.py:141-152 | 12 lines for one requirement and one mechanism, sitting under a 6-line paragraph of repo history (:134-139, unchanged). | Keep 5: `The null simulated here has to be the one the theory assumes: the truncation time is chosen without reference to the events (`checks/result`). Deriving tau from the draw as `g.sum() + g.mean()` makes it `T_N(1 + 1/n)`, so `T_N/tau` is pinned at `n/(n+1)` where the null has it Beta(n, 1), the event count is fixed where it should be random, and the leftover window carries no variance. `_exponential_segment` fixes tau and lets the count fall where it falls.` Deletes ~6 lines.

MINOR | instrument_validation.py:449-452 | The rename left a dangling antecedent in a shipped artifact string: the C2 row now reads "size measured at tau=20: ... Bench 0.0650 and 0.0865 ... other rows at the same n and shape run to 0.176". After the change the row no longer states an n, so "the same n" refers to nothing, and a reader sees tau=20 and n=20 in one row as if they were the same quantity. | `other rows at n=20 and the same shape run to 0.176`.

MINOR | instrument_validation.py:117-119 vs calibration_summary.py:294-296 | The commit's insisted-on vocabulary ("A TRUNCATION TIME, not an event count") is contradicted by the function that implements it: `_exponential_segment(n: int, rng)` still calls the argument `n` and the new caller passes `tau` into it positionally. Rule 3 vocabulary should hold at the definition, not only at the call site. | Rename that parameter to `tau: float` (code finding, but the prose is what makes it visible).

OPTIONAL | jobs/active/instrument_validation.py:46 (unchanged line, now inaccurate) | "these two decide the three tier-3 numbers" - three inputs decide them. `ASYMPTOTIC_SIZE_REPLICATES = 1200` is equally decisive and is NOT declared in the job; it reaches meta from the module constant (instrument_validation.py:571). | `these decide the three tier-3 numbers, together with ASYMPTOTIC_SIZE_REPLICATES in the analyzer`.

VERIFIED | instrument_validation.py:117 | "The mean gap is 1" - `rng.exponential(size=...)` uses scale 1.0. True.
VERIFIED | instrument_validation.py:144-146 | `g.sum() + g.mean() == T_N(1 + 1/n)` and `T_N/tau == n/(n+1)` are correct algebra; `T_N/tau ~ Beta(n, 1)` is correct for fixed-tau truncation conditional on N = n.
VERIFIED | instrument_validation.py:151-152 | `_exponential_segment` (calibration_summary.py:294-299) does fix tau and take a random count. Claim true.
VERIFIED | instrument_validation.py:432,452,513 | "size measured at tau=20" is a literal, but it is effectively pinned: test_instrument_validation.py:268-277 parses the tau out of the row, re-measures, and compares to the quoted size, so a changed default breaks the test. (Note for the code reviewer: the regex is `tau=(\d+)`, which will not match if the label is ever interpolated as "tau=20.0".)

### tests/test_tier3_calibration.py  [reviewed]

MUST FIX | test_tier3_calibration.py:33 | "reports 0.0342 / 0.0292 / 0.0383" - MEASURED and reproduced (I ran it: n=20, seed=777, 1200 replicates, event-derived tau gives exactly 0.0342 / 0.0292 / 0.0383), but no code in the repo produces them any more: that generator was removed from `measure_asymptotic_size` in this very commit. Three numbers no shipped artifact regenerates and no test asserts, describing a code path that no longer exists - rule 2 and rule 1 at once. | Drop the three numbers and keep the direction: `...removes eq (7)'s tail variance and turns the measured rate conservative - the size of a sampling scheme the theory does not cover.`

MUST FIX | test_tier3_calibration.py:231 | Docstring summary "The property whose failure made the reported sizes the size of another null" narrates a past defect instead of naming the property under test. | `The size generator truncates on time: tau fixed in advance, event count random.`

MUST FIX | test_tier3_calibration.py:233-234 | Same FALSE battery claim as instrument_validation.py:143 - "`checks/battery` refuses CvM-asymptotic where that does not hold". It gates on `include_tau_checks` and segment count, and cannot see where tau came from. | `checks/result` requires tau to be chosen without reference to the events; nothing downstream detects a violation.`

SHOULD FIX | test_tier3_calibration.py:23-24 | "tau = 20 gives C1 0.0708, C2 0.0700, CvM 0.0642" - VERIFIED (I ran `measure_all_asymptotic_sizes(tau=20.0, seed=777)`: 0.0708 / 0.0700 / 0.0642 exactly). But NOTHING pins them: the only size test runs at `tau=30.0` and asserts the wide bound `0.005 <= size <= 0.20`. They will drift silently if `ASYMPTOTIC_SIZE_REPLICATES` (1200) changes or if `_exponential_segment`'s draw size formula changes the RNG stream, and they are not reproducible from the docstring alone because it names tau but not the seed or the replicate count. | Either say `at the shipped defaults (seed 777, 1200 replicates)`, or add an assertion at tau=20 with abs tolerance so the docstring has a keeper.

SHOULD FIX | test_tier3_calibration.py:236-238 | Half the docstring describes a generator that does not exist and that the test does not exercise; the assertions are `len(set(counts)) > 1` and `np.std(used) > 1e-6`, neither of which touches `n/(n+1)`. "- to every decimal place, on every replicate -" is amplification of "exactly". | `Both consequences are asserted: the event count varies, and the fraction of the window used up by the last event varies. A tau derived from the draw pins the second at n/(n+1) and the first at n.` (n/(n+1) is algebra, not a measurement, so it needs no pinning - and it is correct.)

MINOR | test_tier3_calibration.py:33 | "removes eq (7)'s tail variance" is asserted for all three instruments; only C2 uses eq (7). Same defect as instrument_validation.py:147. | Scope it to C2 or name the leftover window rather than the equation.

MINOR | test_tier3_calibration.py:30-35 vs :107 | The docstring condemns the event-derived tau, and `_iid_segments` at :107 in the same file still builds `Segment(x=g, tau=float(g.sum() + g.mean()))` for the C5/C6 permutation nulls. Defensible (tau does not enter a permutation null the same way) but a reader will hit the contradiction. | One clause: `...a scheme the ASYMPTOTIC theory does not cover; the permutation checks below are indifferent to it.`

OPTIONAL | test_tier3_calibration.py:19-21 (unchanged) | "An earlier draft of this docstring quoted '0.069' and '0.0757' ... the first instance CLAUDE.md's claims-discipline section records, reintroduced here." Internal chronology plus a pointer to a process document - rule 1, in the docstring the new prose sits inside. Not in this diff, but it is what a reader meets first. | Cut to the standing rule: `Pooled means over the two shapes ("0.069", "0.0757") are in no cell of the table and must not be quoted as one.`

### docs/JOBS.md  [reviewed]

VERIFIED | docs/JOBS.md:52 | Only the parameter row changed (`ASYMPTOTIC_SIZE_N=20` -> `ASYMPTOTIC_SIZE_TAU=20.0`). The file is regenerated by `scripts/make_job_manifest.py` and its staleness is asserted by `tests/test_job_manifest.py:27`. Regenerated and pinned - no prose finding.

### tests/test_instrument_validation.py  [reviewed]

VERIFIED | no new prose. The two changed lines are code (regex and call site); the surrounding comment at :272 is unchanged and still accurate.

### jobs/active/instrument_validation.py  [reviewed]

VERIFIED | no new prose. One constant renamed; the comment above it is unchanged (see the OPTIONAL note under instrument_validation.py about "these two").

PROSE VERDICT: DO NOT SHIP as prose. 6 MUST FIX, of which 3 are false or unsupportable
claims (tlf.py:48 and recipes.py:184 on where the seed lands; instrument_validation.py:143
and test_tier3_calibration.py:233 on battery refusing CvM-asymptotic) and 2 are rule-1
chronology (tlf.py:34, test_tier3_calibration.py:231). Numbers: 0.0708/0.0700/0.0642 and
0.0342/0.0292/0.0383 both reproduce exactly, but neither triple is asserted anywhere and the
second describes a deleted code path. ~45 lines deletable outright.
CORRECTION to the line count above: the trims listed total ~55 deletable prose lines, not 45
(interpolate ~16, filter ~11, within_calibration ~7, _multiprocess 4, tlf ~10, recipes 2,
tlf_plot 1, instrument_validation ~7, test_tier3 ~7).
Also measured across the added prose in scope: 22 dash-aside constructions (" - x - ") in
~58 added comment lines.

### tlf.py + recipes.py + tlf_plot.py

IMPORTANT | recipes.py:293 (with tlf.py:47-49, recipes.py:183-186) | Change 8: the seed does NOT reach `parameter_row`, and both the code comment and the `tlf.run` docstring claim it does ("Passed as a step argument so it folds into that identity" / "the seed arrives here the way `xi_seed` does, as a step argument, so it reaches both the identity and the provenance label"). It is captured in the `_tlf_step` CLOSURE, not passed as a step kwarg. Measured: `parameter_row` for a job with `job.step(_tlf_step(20260902), raw, name="tlf")` returns `['raw()', 'tlf()']` - no seed. Compare recipes.py:270, `job.step(_fidelity_panel_data, ..., xi_seed=xi_seed)`, which yields `fidelity_panel_data(xi_seed=20260813)`. Consequences: (a) `provenance.pipeline_steps` is built by `runner._format_step`, which renders `node.kwargs`, so the shipped provenance record for the tlf node names no seed at all - the record does not state which seed produced the reported `bic_delta`/`is_bimodal`/dwell numbers; (b) identity still changes when `TLF_SEED` changes, but only because `code_closure` hashes `quebra/recipes.py` bytes - i.e. any edit to recipes.py, not the seed specifically. This is exactly the failure `parameter_row` was written to close, and the project rule is that a seed affecting a reported number is a declared parameter. | Mirror `xi_seed`: `def _tlf_step()` returning `def step(norm, seed)`, called as `job.step(_tlf_step(), final_filtered, name="tlf", seed=tlf_seed)`.

VERIFIED | tlf.py:40 | `seed` is keyword-only and required, so no caller can silently take a default. The only caller in the repo is recipes.py:135. `jobs/` never overrides `tlf_seed`.

MINOR | tlf.py:43-46 | Change 8: the nondeterminism the docstring asserts is not exhibitable on this data. `GaussianMixture` is 1-D here, and sklearn 1.8.0's k-means++ init converges to the same optimum every time: `random_state=None` repeated 8x on a deliberately multi-optimum 1-D mixture (three clusters fitted with 2 components) gives byte-identical `bic_delta` and `means_`, and so do seeds 0..7. So `random_state=seed` changes no currently reported number - it is correct hardening, but the claim "become a fresh random variable on each call" overstates what could be measured. | Keep the change; soften the claim to a contract statement.

VERIFIED | tlf.py:16,37,74-92 | `bic_delta: float | None`, `gmm2=None`, `fit_failed=True` on the failure path is correct and complete. The success path (tlf.py:98) always assigns a `float`, so the widened annotation costs nothing. `tlf_plot.py` is the ONLY reader of `bic_delta` or `gmm2` anywhere in src/, jobs/, tests/ or scripts/, and it handles both (`result.bic_delta is not None` at :252, `result.gmm2 is None` at :90). The old fallback genuinely was a silent-fallback violation: `bic_delta=0.0` with `gmm2=gmm1` reported "no evidence of bimodality" from a comparison that never ran.

MINOR | tlf.py:9,37 | Change 8, the pickle question, answered definitively. `TLFResult` is `slots=True`, so a new defaulted field is NOT safe for unpickling: a slots dataclass pickles as `(None, slots_dict)` and `__setstate__` never applies field defaults, so an old pickle loads with the `fit_failed` slot UNSET and any direct access raises `AttributeError: 'TLFResult' object has no attribute 'fit_failed'` - reproduced. `StaleArtifactGuard` cannot help: it is explicitly for NON-slots dataclasses (`_artifact_guard.py:22`) and `TLFResult` does not inherit it, so the repo's stale-artifact policy does not cover this class. In practice the hazard is moot right now: all 37 `*_tlf.pkl` files under `output/` were written before the src-layout move and reference `analyzers.tlf`, so they already fail with `ModuleNotFoundError` and none is loadable. The only reader is defensive (`getattr(result, "fit_failed", False)` at tlf_plot.py:89), which is what keeps this MINOR rather than IMPORTANT - but the `getattr` is load-bearing and nothing says so. | Either give `TLFResult` the guard by dropping `slots=True`, or note at tlf_plot.py:89 that the `getattr` default is required because `TLFResult` is a slots dataclass.

VERIFIED | tlf_plot.py:86-116,243-253 | `fit_failed` is bound at :89 on the straight-line path through `render`, before its uses at :113 and :247 - no branch reaches :247 with it unbound. `means2_axis = np.asarray([])` in the failed branch correctly short-circuits :118, :146 and :155, so the lobe annotation, `lobe_khz` and the per-lobe lines all fall to their existing "n/a" paths, and `covs2_axis`/`weights2` (only referenced at :158-159, inside `len(means2_axis) >= 2`) are never touched unbound. Behaviour for a pre-change pickle with `gmm2 = gmm1` is unchanged: `fit_failed` is False, `means2` has one element, and the figure takes the same "n/a" branches it already took.

MINOR | tlf.py, tlf_plot.py, recipes.py | Change 8 is entirely untested. `grep -rln "tlf\|TLF" tests/` returns nothing: there is no test for `tlf.run` at all, so the required `seed`, the `fit_failed` field, the `bic_delta=None` path and the suppressed GMM(2) curve are all unpinned. | One test: `run(..., seed=1)` twice is identical; a values array that forces a singular 2-component fit gives `fit_failed is True and bic_delta is None and gmm2 is None`.

## Reviewed (patch 3)
- [x] src/quebra/analyzers/tlf.py
- [x] src/quebra/recipes.py
- [x] src/quebra/plots/tlf_plot.py

### instrument_validation.py + jobs/active/instrument_validation.py + tests

CRITICAL | jobs/bench/results/instrument_report.md:29,30,34 | Change 9 moves the three tier-3 numbers but the COMMITTED artifact that publishes them is not regenerated. The tracked file still reads "size measured at n=20: 0.0342" (C1), "0.0292" (C2), "0.0383" (CvM) - both the values AND the `n=`/`tau=` label are stale. Proven: `python3 jobs/bench/instrument_report.py` on this working tree rewrites exactly those three cells to "size measured at tau=20: 0.0708 / 0.0700 / 0.0642". The file's own header says "GENERATED by `jobs/bench/instrument_report.py`", and `measure_asymptotic_size`'s docstring justifies its whole existence as "the tier-3 rows quote a number this run produced rather than one a docstring remembers" - which is now false of the shipped artifact, in the opposite direction (it claims CONSERVATIVE where the code says ANTI-CONSERVATIVE). | Re-run `python3 jobs/bench/instrument_report.py` and commit the result in the same commit.

VERIFIED | instrument_validation.py:124-165 | Change 9 is statistically correct and the numbers reproduce. Measured on this tree: NEW = C1 0.0708, C2 0.0700, CvM 0.0642. Reimplementing the old generator (`g = rng.exponential(size=20); Segment(x=g, tau=g.sum()+g.mean())`, seed 777, 1200 reps) gives exactly 0.0342 / 0.0292 / 0.0383 - the before-values are confirmed. `_exponential_segment` IS a correct time-truncated HPP observation: `kept = gaps[:searchsorted(cumsum(gaps), tau, "left")]` keeps every gap whose cumulative endpoint is strictly below a tau fixed in advance, so N is random (measured mean 19.99, min 4, max 41 over 20000 draws at tau=20) and the residual `tau - sum(x)` is the one censored gap. The docstring's diagnosis of the old scheme is exactly right: conditional on N=n, the ordered event times of an HPP on [0,tau] are uniform order statistics, so `T_n/tau ~ Beta(n,1)`; `tau = T_n(1+1/n)` pins that ratio at `n/(n+1)`, i.e. at the mean of the correct distribution with zero variance, and eq (7)'s tail term depends on it. Using the same generator as `size_vs_n`/`validation_curve` also makes the tier-3 rows commensurate with the bench numbers quoted beside them.

VERIFIED | instrument_validation.py:154, calibration_summary.py:297 | The pool cannot exhaust. `size = int(tau + 10*sqrt(tau) + 50)` is 114 at tau=20; `P(Gamma(114,1) < 20) = 2.0e-47`. At the bench's largest, n=355, the pool is 593 and `P = 6.8e-31`. So the question "can it exhaust before reaching tau" answers NO for every value shipped. There is still no guard: if it ever did, `searchsorted` returns `len(gaps)` and the segment is returned silently truncated at the pool's end rather than at tau - a silent fallback with no marker. | Add `if cumsum[-1] < tau: raise` in `_exponential_segment` (that file is outside this diff's scope; note it for the next patch).

VERIFIED | instrument_validation.py:154 | `n_censored_dropped=1` is right for this use. A time-truncated observation has exactly one incomplete final gap (the residual `tau - sum(x)`), which is what the field documents at `checks/result.py:88-91,116-119`. It is a pure diagnostic: grep shows it is only summed into `CheckResult.n_censored_dropped` by c1/c3/c5/c6 and never enters a statistic, so it changes no number.

IMPORTANT | instrument_validation.py:432,452,513 | The report quotes 4 decimals on a number whose Monte Carlo standard error is 0.0074. At 1200 replicates and p ~ 0.07, `sqrt(p(1-p)/1200)` = 0.00741 / 0.00737 / 0.00707 for C1 / C2 / CvM - measured. So `:.4f` overstates the precision by roughly 75x; the third decimal is already noise and the fourth is pure formatting. Only the coarse verdict survives: (0.0708 - 0.05)/0.0074 = 2.8 sigma, so "anti-conservative" is real, but "0.0708" is not. The same file's bench rows beside these DO carry an SE column (`size_table.csv`). | Quote 2 decimals with the SE, e.g. "0.071 +- 0.007 (1200 replicates)", or raise `ASYMPTOTIC_SIZE_REPLICATES` to ~20000 (SE 0.0018) if the fourth decimal is wanted. Whichever, `test_instrument_validation.py:268`'s `(0\.\d{4})` regex and `abs=5e-5` tolerance move with it.

IMPORTANT | tests/test_tier3_calibration.py:226-252 | `test_the_size_generator_truncates_on_TIME_not_on_events` CANNOT FAIL on the change it claims to pin, and neither can its companion. It asserts properties of `_exponential_segment`, which lives in `calibration_summary.py` - a file this diff does not touch at all (`git diff --stat` on it is empty), so the test passes identically before and after and would have passed at HEAD. PROVEN: monkeypatching `measure_asymptotic_size` back to the old generator while KEEPING the `tau` parameter name and the "tau=20" label reproduces 0.0342 / 0.0292 / 0.0383, and both `test_the_size_generator_truncates_on_TIME_not_on_events` and `test_the_tier_3_rows_quote_a_number_this_run_produced` still pass (2 passed). The second test is self-consistent, not anchored: it recomputes `expected` from the same reverted code it is checking. `test_the_asymptotic_size_is_in_the_documented_range`'s 0.005..0.20 bound also admits 0.0342. So the only thing standing between a semantic revert and green gates is the literal string "n=" vs "tau=". | Assert on the function that changed: `measure_asymptotic_size` must produce a varying event count. E.g. capture the segments it builds (or expose the generator as an injectable parameter) and assert `len(set(counts)) > 1`; or pin the three sizes against the report's own numbers loaded from `jobs/bench/results/instrument_report.md` so the artifact and the code cannot drift apart silently.

MINOR | instrument_validation.py:432,452,513 | The tier-3 detail strings hardcode "tau=20" while the value comes from the `asymptotic_size_tau` parameter. A job passing `asymptotic_size_tau=30.0` would publish a size measured at 30 under a "tau=20" label. `test_instrument_validation.py:268-276` happens to catch it (it recomputes at the label's value and compares), but the string itself is a lie waiting for a parameter change. Pre-existing with `n=20`, but the diff rewrote these exact three lines. | `f"size measured at tau={asymptotic_size_tau:g}: ..."`.

MINOR | instrument_validation.py:154 | `_exponential_segment` is annotated `(n: int, ...)` and its parameter is named `n`, and the new call passes `tau: float = 20.0`. It works (`float(20.0)`, `int(20.0 + ...)`), but the helper's name and type now contradict the meaning the caller was just renamed to state. `mypy` will not catch it: `make types` is scoped to `src/quebra/core` only. | Rename the helper's parameter to `tau: float` in the next patch (that file is out of scope here), or the rename this commit exists for is only half done.

MINOR | instrument_validation.py:157 | `rejected += int(p is not None and p <= alpha)` silently counts a `p_value=None` replicate as a non-rejection, which deflates the measured size - the exact "silence is not a pass" case. Pre-existing, but change 9 makes N random, so a too-small draw is now possible where n was formerly fixed at 20. Measured at tau=20/seed 777: 0 of 1200 replicates return None for any of the three modules and min N is 7, so it does not fire today. The adjacent risk is worse: at N < 2 the checks RAISE (`this check needs at least 2 complete gaps in a segment`), aborting the whole job; `P(N<2) = e^-20 * 21 = 4.3e-8` per replicate, ~1.5e-4 over the 3600 replicates the report runs. Deterministic under seed 777 and it does not happen, but a seed or replicate-count change could hit it and the old generator could not. | `if p is None: raise` (nothing legitimate produces None here), and either guard `_exponential_segment` for N>=2 or state the risk.

MINOR | instrument_validation.py:143 | An analyzer now imports a PRIVATE name from a sibling analyzer (`from quebra.analyzers.calibration_summary import _exponential_segment`). `lint-imports --no-cache` passes (1 contract kept, layering is package-level), so no rule is broken, and sharing the generator is the right call for commensurability with the bench. But `_exponential_segment` is now a two-caller shared contract wearing a private name. | Drop the underscore, or move it to a shared `checks/_nulls.py`.

VERIFIED | jobs/active/instrument_validation.py:47,80 | The rename reaches the declared job parameters and therefore the identity and the provenance label: `ASYMPTOTIC_SIZE_TAU = 20.0` is a module constant in the job file (hashed into `job_code_hash`) AND a step kwarg (`asymptotic_size_tau=ASYMPTOTIC_SIZE_TAU`), so `parameter_row` and `runner._format_step` both carry it. Contrast the `tlf_seed` finding above - this is how it should be done, in the same diff.

VERIFIED | docs/JOBS.md:52 | Matches the generator. `python3 scripts/make_job_manifest.py` on this tree reproduces the diff's line byte-for-byte (`ASYMPTOTIC_SIZE_SEED=777`<br>`ASYMPTOTIC_SIZE_TAU=20.0`), and the regenerated file leaves `git diff --stat` at exactly 1 insertion / 1 deletion. File restored to the working-tree state after the check.

VERIFIED | tests/test_instrument_validation.py:268,276 | The regex/`float()` change is consistent: the label renders "tau=20", `(\d+)` captures "20", `float("20")` is 20.0 and `measure_all_asymptotic_sizes(tau=20.0)` is what the artifact used. No stale `n=` reference remains in either test file.

## Reviewed (patch 3)
- [x] src/quebra/analyzers/instrument_validation.py
- [x] jobs/active/instrument_validation.py
- [x] tests/test_tier3_calibration.py
- [x] tests/test_instrument_validation.py
- [x] docs/JOBS.md

## Manifest (patch 3) - empty, all files reviewed

### Empirical confirmation on the real shipped dataset (100423_6D2S_qubit1)

VERIFIED | end-to-end | `_load_dataset` -> `lookup_prior` -> `filter.run` -> `interpolate.run` -> `tlf.run` on the real 2024-read dataset: every row-aligned column has ZERO non-finite values (`T2star_error_s`, `T2star_s`, `chi_squared`, `delta_hz`, `qubit_frequency_hz`, `rabi_hz`, `raw_frequency_hz`, `t_rel_s`), `meta["n_nonfinite_dropped"]` is absent, and no new raise fires. Confirms that changes 1, 2, 3, 4, 5 are all inert on shipped data: with `finite` all-True, `_interpolate_finite` reduces to the old `pchip_interpolate(t_rel_s, nan_to_num(v), x)` bit-for-bit, so change 1 moves no committed number.

VERIFIED | end-to-end | The norm's actual key set is `['T2star_error_s','T2star_s','chi_squared','delta_hz','meta','raw_frequency_hz','t_rel_s']` plus `qubit_frequency_hz`/`rabi_hz` from `lookup_prior`. Four of the seven row-aligned columns are NOT in `filter._ROW_ALIGNED_KEYS`, and the set's `t_unix_s` is not among them - the incompleteness finding above, confirmed against real data rather than inferred.

VERIFIED | tlf | On the real values, `random_state=None` repeated 5x and seeds 0/1/20260902 all give `bic_delta = 1999.35945` and identical `means_`. So change 8 also moves no committed number; the seed is hardening, not a correction.

VERIFIED | gates | Could not falsify. `ruff check .` clean, `ruff format --check .` 148 files formatted, `mypy src/quebra/core` clean, `lint-imports --no-cache` 1 contract kept, `deptry .` clean, `pytest -m "not slow and not heavy and not real and not r"` 433 passed / 2 skipped in 30s. Note what the gates do NOT cover: nothing under `tests/`, `Makefile` or `.github/workflows/` references `instrument_report.md`, so the CRITICAL artifact drift above passes `make check` silently.

## Verdicts (patch 3)
1 interpolate NaN masking .......... CORRECT-BUT-INCOMPLETE (pct_real_points not recomputed against the reduced support)
2 generic-loop except removal ...... CORRECT-BUT-INCOMPLETE (2-D row-aligned column now IndexErrors; disagrees with filter.py on ragged values)
3 t_rel_s / raw_frequency_hz raises  CORRECT (unfirable on shipped data; nothing relied on the drop)
4 _ROW_ALIGNED_KEYS ................ CORRECT-BUT-INCOMPLETE (4 real row-aligned keys missing, 1 dead entry, first stated hazard not actually closed)
5 sigma stage finite handling ...... CORRECT (and inert: delta_hz cannot be non-finite downstream of RamseySeriesSchema)
6 cumulative NaN padding ........... CORRECT (real crash fix, correctly wired, complete; untested)
7 interior-censored guard .......... CORRECT-BUT-INCOMPLETE (calendar clock unguarded; fires before the min_events drop)
8 TLF seed + fit_failed ............ CORRECT-BUT-INCOMPLETE (seed absent from parameter_row and the provenance label; slots-unpickle hazard; untested)
9 time-truncated size null ......... CORRECT (the null is right and the numbers reproduce) but the COMMITTED ARTIFACT WAS NOT REGENERATED and the new test cannot fail

## Patch 3 verdict: DO NOT SHIP
## (1 CRITICAL: jobs/bench/results/instrument_report.md still publishes 0.0342/0.0292/0.0383
## under an "n=20" label, contradicting the code in the opposite direction. Plus one test
## proven unable to fail on the change it names, and the tlf seed missing from parameter_row.)

# PATCH COMMIT 4 of 4 - BEHAVIOUR
SCOPE: tests/test_transform_guards.py (new), tests/test_windows_not_interpolated.py,
analyzers/checks/_multiprocess.py, transforms/interpolate.py, plots/tlf_plot.py,
panels/_check_ledger_render.py (+check_ledger.py context), tests/test_instrument_validation.py
COMMIT: ff4a8c5

## Manifest (patch 4)

## Reviewed (patch 4)
- [x] tests/test_transform_guards.py  (NEW, untracked)
- [x] tests/test_instrument_validation.py  (+5)
- [x] src/quebra/panels/_check_ledger_render.py  (+1)
- [x] src/quebra/plots/tlf_plot.py  (+4)
- [x] src/quebra/transforms/interpolate.py  (+9)
- [x] src/quebra/analyzers/checks/_multiprocess.py  (+23)
- [x] tests/test_windows_not_interpolated.py  (+92/-?)

## Findings (patch 4)

### tests/test_windows_not_interpolated.py  (change 1)

VERIFIED | tests/test_windows_not_interpolated.py:58-67 | `discover(jobs/)` yields exactly the
12 files carrying a module-level `JOB_ID`; the 8 `jobs/bench/*.py` study modules, 4
`__init__.py` and `jobs/reference|rscripts` are skipped. Collection is 16 tests = 12 params +
3 controls + the non-vacuity test. `independence_survey` IS reached: 206 nodes, 34 windows
nodes, 34 labels, `_unexpected` empty on all 34. No import-time I/O in the composite, so the
widening does not make the test depend on `data/real_private/`.
VERIFIED | tests/test_windows_not_interpolated.py:116-123 | Loop terminates on every degenerate
name and fails CLOSED on all of them: `interpolate`, `_final`, `final_`, `_`, `__`, `___`,
`a_`, `_a`, `filter_`, `windows_x_` all land in `unexpected`. `base` loses at least one
segment per iteration and the guard is `while suffix and base`, so no infinite loop is
reachable.
VERIFIED | tests/test_windows_not_interpolated.py:116-121 | The `final` / `final_filter_stage`
overlap does NOT mis-strip here: the loop tries every right-to-left split and accepts if ANY
(base in whitelist, suffix in labels) pair holds, so `final_filter_stage_040423_6D2S_qubit1`
is accepted via the long prefix and `final_smoothed_q1_040423` is still rejected. Confirmed by
running km_poster_6d2s and the survey.
VERIFIED | src/quebra/recipes.py:145-152 | `_final_stage` unwraps `FilterResult.final_norm`;
it resamples nothing, so adding `"final"` to the whitelist is factually justified.

IMPORTANT | tests/test_windows_not_interpolated.py:89-95 | `_carve_labels` is self-justifying,
and demonstrably defeatable. `_windows_nodes` matches ANY node whose `fn_name` contains
"windows", so a resampling step named `windows_pchip` sitting between `filter` and `windows`
(a) enters `_windows_nodes`, (b) contributes its own label `pchip`, and (c) is then vouched
for by that label when the real carve's ancestry is checked. Executed: `_unexpected(job,
"windows") == set()` and `_unexpected(job, "windows_pchip") == set()` - the guard passes on a
graph that carves from a resampled series. Same mechanism without any "windows" in the
resampler's name: the pair (`filter_smoothed`, `windows_smoothed`) also passes, i.e. the exact
name the new control at :225 asserts is rejected becomes accepted the moment the carve node
echoes it, which is the natural thing for a job author to write. | Two parts. (i) Pass the
root's own label into `_unexpected` instead of the job-wide union, so one chain's label cannot
vouch for another chain's node - that alone kills the `windows_pchip` case. (ii) A label set
derived from the graph's own names can never be a check on those names: intersect the derived
labels with an external source (dataset stems, as before, UNION an explicit
`SURVEY_LABELS = {"q1_040423", ...}` constant listing the 34), so admitting a new label stays
a deliberate edit to this file.

IMPORTANT | tests/test_windows_not_interpolated.py:142-156 | The widened guard is vacuous on 6
of the 12 jobs (check_calibration, instrument_validation, mtbc_q6, mtbf_q1,
ramsey_2x2_q1_030723, check_ledger_q1, compare_t2star_0704_vs_1004 have zero windows nodes),
and nothing asserts the matcher still matches anywhere. Rename the carve step from `windows`
to e.g. `segments` and all 12 parametrised cases pass with an empty loop while the three
controls keep passing on their hardcoded names - all 16 green, zero real coverage. This is the
precise failure mode the module docstring says was already hit once by `fidelity_windows`. |
Extend `test_there_are_jobs_to_check` to assert the carving jobs are still seen, e.g.
`{"km_poster_6d2s", "independence_survey", "t2star_q1_070423", "t2star_q1_100423",
"ramsey_q1_100423"} <= {f.stem for f in _job_files() if _windows_nodes(_load_job(f))}`.

MINOR | tests/test_windows_not_interpolated.py:108-111 | The blacklist runs BEFORE the
whitelist membership test and `continue`s, so the second layer overrides the first: a name
containing a token can never be whitelisted. Executed: re-adding `"fidelity_interp"` to
`OBSERVED_ANCESTORS` still yields `_unexpected(...) == {"fidelity_interp"}`. The failure
message at :152-156 tells the reader to "add it to OBSERVED_ANCESTORS deliberately", a remedy
that silently does not work for any name containing interp/resample/regrid/uniform/decimate/
rebin - and "uniform" is a plausible substring of a non-resampling step (a uniform prior). |
Move the `name in OBSERVED_ANCESTORS` test above the token scan, so an explicit whitelist
entry wins over the heuristic layer.

MINOR | tests/test_windows_not_interpolated.py:92-94 | `_carve_labels` adds a label for EVERY
whitelisted prefix that matches, not the longest. Adding `"final"` made two entries prefixes
of a third, so a carve node named `final_filter_stage_q1` now contributes both `q1` and the
bogus `filter_stage_q1`, and the bogus one then vouches for any `{whitelisted}_filter_stage_q1`
name. Inert today (no carve node uses a `final*` prefix) but activated by any future rename. |
Strip only the longest matching prefix: `known = max((k for k in OBSERVED_ANCESTORS if
name.startswith(f"{k}_")), key=len, default=None)`.

### src/quebra/analyzers/checks/_multiprocess.py  (change 2)

VERIFIED | _multiprocess.py:294-331 | `interior_censored` is bound on every path reaching the
guard: `clock` is validated against the two constants at :238-241, so the if/else at :294/:311
is total, and both arms assign it. No `UnboundLocalError` is reachable.
VERIFIED | _multiprocess.py:305,325 | Condition meaning unchanged. Old
`bool((~complete[:-1]).any())` vs new `int((~complete[:-1]).sum())` in a boolean context are
identical, including the single-window block where `complete[:-1]` is empty (sum 0, any False).
Executed: a one-row `gap_start` block still yields a Segment, not a raise.
VERIFIED | _multiprocess.py:320-331 | Ordering behaves as intended. Executed: a 2-row block
with an interior `gap_start` and `min_events=2` now returns `([], 1)` instead of raising; a
5-row block with the same defect still raises naming position 0; the calendar clock with the
same frame does not raise and yields tau=9.0.
VERIFIED | _multiprocess.py:325 | The guard still cannot fire from a real job. `windows.carve`
derives BOTH `DEATH_GAP_START` (windows.py:184) and `diagnostics["gap_spans_s"]`
(windows.py:496-498) from the SAME `is_gap` mask, so every gap that can censor a window is in
the list; `_segment_starts` then splits at the first birth at/after each gap end, which puts
every `gap_start` death last in its block. `DEATH_SCAN_END` is only reachable for the final
window. All three real call sites (recipes.py:196, recipes.py:380, check_ledger.py:200) pass
`diagnostics.get("gap_spans_s")`, and an empty list is only possible when there are no gaps
and therefore no `gap_start` deaths. Defensive only, as claimed.
VERIFIED | _multiprocess.py:325 | The added `clock == CLOCK_IN_SPEC` conjunct is redundant with
`interior_censored = 0` on the calendar arm, but harmless and it keeps the invariant local.

MINOR | _multiprocess.py:320-331 | The reorder makes the SAME data-integrity defect fatal or
invisible depending on block size: a caller that forgets `gap_spans_s` and produces an
interior-censored block of 5 windows aborts the run, while the identical defect in a 2-window
block is folded into `n_dropped` alongside ordinary too-short blocks and is never surfaced. No
wrong number reaches an artifact (the block contributes nothing either way), so it is minor,
but the drop reason is now conflated. | Count them separately, e.g. return the interior-
censored drop count in the diagnostics, or `print` the position once when a dropped block also
had `interior_censored`.

MINOR | _multiprocess.py:325-331 | Untested on both sides of the move. `grep` over `tests/`
finds no test that provokes the interior-censored raise and none that asserts a small
interior-censored block is dropped rather than aborting - i.e. the behaviour this change exists
to create has no regression guard, and a future revert to the pre-drop position would be green.
| Add the two cases exercised above to `tests/test_checks_statistics.py`.

### src/quebra/transforms/interpolate.py  (change 3)

VERIFIED | interpolate.py:174-183 | No real path now raises where it previously passed
through. On the recipe path (recipes.py:239-243) interpolate is fed `_final_stage(FilterResult)`
and `filter._subset_norm` (filter.py:47-54) has already raised on any un-arrayable value, so
the new arm is unreachable there. On the filter-less user path a loaded Norm comes through
`RamseySeriesSchema`/pandera as ndarray columns, so it is unreachable in practice too. The
revert is the right call: fail-closed and, on the only path where it can fire, it is the only
guard there is.
VERIFIED | interpolate.py:176 | Executed: with numpy 2.4.4 only a ragged nested sequence
actually raises (ValueError). dict / generator / set / str / bare object / an object with a
broken `__len__` all box to a 0-d object array and take the `arr.ndim == 0` passthrough at
:186. So the arm is neither too broad in effect nor mislabelled for the case it catches.

IMPORTANT | interpolate.py:185-188 | The parity with `filter._subset_norm` that this change
invokes is only half achieved, and the missing half is the silent one. `_subset_norm` RAISES on
a row-aligned column whose length disagrees (filter.py:57-62); interpolate still passes it
through untouched. Executed on a filter-less Norm: `t_rel_s` of length 4 and `chi_squared` of
length 3 gives an output whose `chi_squared` is still length 3 against a 4-point uniform grid,
no error, no diagnostic. That is exactly the "user job wires interpolate alone" case cited to
justify the raise nine lines above, and it is a misalignment rather than an unresampled value.
Pre-existing, but the reasoning that reinstated the neighbouring guard applies verbatim here. |
Raise for `arr.ndim > 0 and len(arr) != len(order)` with the same message shape
`_subset_norm` uses, since after interpolate the grid length is `len(order)` too, so no
legitimate column can have another length.

MINOR | interpolate.py:176 | `TypeError` in the tuple is dead code under numpy 2.x for the
reason measured above, and a hostile `__array__` raising anything else (measured: RuntimeError)
escapes the named message. Harmless - it still fails loudly - but the arm is narrower than it
reads. | Catch `Exception` and re-raise named, or drop `TypeError`.

### src/quebra/plots/tlf_plot.py  (change 4)

VERIFIED | tlf_plot.py:88 | No construction path lacks the attribute. `TLFResult` has only two
construction sites (tlf.py:69, tlf.py:199), both keyword dataclass calls, so the
`fit_failed: bool = False` slot is always set; and `TLFPlot` has exactly one call site
(recipes.py:294) fed by the `tlf` node, which returns a `TLFResult`. Semantics unchanged:
`bool(x) or y is None` parses as `bool(x) or (y is None)` exactly as before.
VERIFIED | tlf_plot.py:88 | The `slots=True` unpickle hazard is not reachable from any artifact
on disk. `fit_failed` and `slots=True` both predate this commit (present at HEAD), and all 140
`output/**/*tlf*.pkl` files fail to load at all with `ModuleNotFoundError: No module named
'analyzers'` (they predate the src-layout move), so `runner.py:362`'s reuse loader can never
hand the plot a slot-less instance from the existing tree. Any newly written pickle carries the
slot. The `getattr` was defending nothing.

### src/quebra/panels/_check_ledger_render.py  (change 5)

VERIFIED | _check_ledger_render.py:24 | The key is right: `battery.ROW_KEYS` uses exactly
`"cvm_cramer_von_mises"` (battery.py:64, :76) and `check_ledger.py:497` writes
`check_id = result.check`, so the entry is live rather than dead.
VERIFIED | check_ledger.py:59-71,149-171 | Nothing breaks at six. Rendered the panel with all
six check ids: the ladder draws 6 lines + the alpha axhline, the legend reads
`['C1 LR','C2 AD','C3 cop','C5 rank','C6 exch','CvM']`, and `ordered_color(i, n)` is a
continuous colormap sample (theme.py:253-265) with no palette length to exhaust and no wrap.
Legend row count is unchanged: `ncol=2` gives ceil(5/2) == ceil(6/2) == 3 rows.
VERIFIED | _check_ledger_render.py:30 | The TABLE does not gain columns. `ledger_grid` pivots
every row present, so the two CvM columns were already drawn under the `check_id[:7]` fallback
label `cvm_cra a` / `cvm_cra p`; the change only renames them to `CvM a` / `CvM p`. Confirmed
in the render: 12 columns before and after.
VERIFIED | tests, artifacts | Nothing pins the old count or the old label: no test imports
`CHECK_SHORT`, `column_label` or `ledger_grid` (grep over `tests/`), and no figure artifact is
tracked in git.

IMPORTANT | panels/check_ledger.py:152-171 (activated by _check_ledger_render.py:24) | The
sixth line makes the ladder's y-axis claim wrong for the check the project just promoted.
`statistic_series` (_check_ledger_render.py:83-88) takes the MIN p-value over calibrations at
each rung, and CvM has BOTH an asymptotic and a permutation row in `ROW_KEYS`. So a point
labelled "CvM" on an axis labelled "p-value", read against the alpha reference line, can be the
asymptotic value - the calibration the battery itself refuses when tau is event-determined
(tests/test_checks_statistics.py:310) and the one the whole size bench exists to distrust at
small n. Nothing in the legend or the marker says which calibration a point came from, and the
line also ignores `verdict`, so an `underpowered` or `not interpretable` cell is drawn as
ordinary evidence. Pre-existing for C1/C2; this change makes it true of CvM. | Restrict the
ladder to `calibration == CALIB_PERMUTATION` rows (add the filter in `statistic_series`), or
key the series on (check, calibration) so the two are separate lines.

MINOR | panels/check_ledger.py:162 | Six categorical series sampled from a SEQUENTIAL ramp:
the rendered colours are `#002b62, #3b496c, #646770, #8b8778, #b6a96f, #e5cf52`, whose middle
four are near-indistinguishable greys, and `ordered_color`'s own docstring says it is for
quantities with a natural order. Marginal at five, worse at six. | Use a categorical palette
for the ladder (check identity has no order), or add per-check markers.

### tests/test_instrument_validation.py  (change 6)

VERIFIED | test_instrument_validation.py:213-215 (deleted) | The deleted line was
`assert ("2 levels" in render_tier_table_markdown(data) or True)`, unconditionally true, so
nothing relied on it. The preceding `assert data.divergence_levels == {35: 2.0, 355: 2.0}`
carries the whole check, and it fails if `divergence_levels` regresses (verified by reading the
FINEST-crossing sibling test, which pins the same field the other way).
VERIFIED | test_instrument_validation.py:26,151,171 | `render_tier_table_markdown` is still
imported and still called twice, so no unused import and no dead helper. The deleted call was
the only render of a data with FINITE divergence in that test, but `_build()` (rendered at :151)
also yields finite values (measured: `divergence_tie_fraction == {35: 0.9, 355: 0.9}`), so the
finite-crossing render path keeps its smoke coverage.

### tests/test_transform_guards.py  (NEW, 16 tests)

All 16 pass as shipped (`pytest tests/test_transform_guards.py -q` -> 16 passed). Each was
attacked by mutating the guard in a COPY of the source and re-running the assertion.

VERIFIED (guard removal proven red), 13 of 16:
- :52 failed-fit-not-zero-filled - mutant that zero-fills before the finite mask: no
  `n_nonfinite_dropped` key (KeyError) AND `grid == 0.0` present. RED.
- :81 outside-support-is-NaN - mutant with the `outside` NaN-ing deleted: `grid[0]` finite and
  `finite.min() == 99.99999999999997`, so both assertions fire. RED.
- :110 non-finite timestamp, :118 too-few-finite, :131 raw_frequency length - each keys on a
  distinct `pytest.raises(match=...)` against a raise that is the only thing producing that
  message. RED by construction.
- :140 / :147 / :170 / :194 filter - the `filter cannot mask` length raise, the
  all-columns-same-length control (regresses to `{t_rel_s: 9, T2star_s: 10}` if a column is
  left unmasked), the propagating-mean regression (mu = NaN makes the mask all-False and
  `kept == 0`), and the all-non-finite raise. RED.
- :180 both-transforms-reject-a-ragged-column - ran HEAD's `interpolate.run` on the same input:
  no raise, `out["odd"] == [[1, 2], [3]]` passed straight through. RED at HEAD, so this test is
  the actual regression guard for change 3.
- :210 the-seed-reaches-both-estimators - monkeypatched `GaussianMixture` to drop
  `random_state`: `gmm1.random_state` becomes None and the assertion FAILS. The replacement
  does what the old agreement form could not. RED.
- :245 failed-fit-reports-fit_failed - the three assertions all read fields set only in the
  `except` block, so either regression (no catch -> ValueError escapes; or the historical
  `gmm2 = gmm1, bic_delta = 0.0`) is red.
- :252-263 `_FailsOnTwo.__getattr__` delegation verified by direct exercise: `bic(X)`,
  `predict(X)`, `means_`, `covariances_`, `weights_`, `random_state` all reach `_inner`, and
  the n=2 arm raises before touching it. The monkeypatch target is right - `tlf.py:60` uses the
  bare module-global name, and patching `tlf_module.GaussianMixture` demonstrably reaches it
  (the seed experiment above used the same route).

IMPORTANT | tests/test_transform_guards.py:125-128 | `test_a_row_aligned_companion_column_is_
resampled_to_the_grid` CANNOT FAIL, so the file's opening claim ("Every assertion in this file
was checked by removing the guard it covers") is false for this one. Proven: a mutant whose
generic column loop is reduced to `out[key] = value` - literally the passthrough the test says
is forbidden - still satisfies the assertion, because `x_uniform` is built with
`num=len(t_rel_s)` (interpolate.py:130), so the grid length ALWAYS equals the input length and
`len(out["T2star_s"]) == len(out["t_rel_s"])` is `10 == 10` either way. Tightening it to a
value comparison does not help with this fixture either: `_norm`'s clock is already uniform
(`np.arange(n) * 10.0`), so resampling is the identity - measured, `out["T2star_s"]` is
`array_equal` to the raw input. | Give the fixture a NON-uniform clock (e.g.
`t_rel_s=[0,10,20,25,40,...]`) and assert the companion column's values equal the pchip of the
companion column on the grid, not merely its length. The same fixture weakness makes
`test_a_clean_input_is_untouched_by_the_masking` (:70) a control over an identity map.

MINOR | tests/test_transform_guards.py:227-235 | `test_the_same_seed_gives_the_same_verdict` is
decoration with respect to the seed, as suspected. Measured under the seed-dropping mutant:
all three assertions still pass, and six unseeded runs return a bit-identical
`bic_delta = 870.932251`. It retains only residual value as a determinism smoke test of the
dwell/BIC arithmetic. The docstring calls it "weaker ... a sanity floor", which is honest, so
this is a note rather than a defect. | Keep it, or fold it into the seed test as a second
assertion so the file does not carry a test that cannot fail on the property it names.

MINOR | tests/test_transform_guards.py:241 | `pytest.raises(TypeError)` with no `match`: the
test passes if `run_tlf` raises TypeError for ANY reason, e.g. a future `np.asarray` type
error on `values`. The real message is "missing 1 required keyword-only argument: 'seed'". |
`pytest.raises(TypeError, match="seed")`.

MINOR | tests/test_transform_guards.py:262-263 | `__getattr__` recurses instead of raising when
`_inner` is unset - measured: `_FailsOnTwo.__new__(_FailsOnTwo).means_` gives RecursionError.
Unreachable in this test (both constructions complete), but it turns any future construction
failure inside the wrapper into an unreadable stack. | Guard with
`if name.startswith("_"): raise AttributeError(name)`.

MINOR | tests/test_transform_guards.py | The file is UNTRACKED. `git add` it before committing
or the whole change-3 regression guard ships as nothing; and note `provenance.is_tree_clean`
uses `--untracked-files=no`, so the tree reports clean while it is missing.

### Cross-check: behaviour hiding in the files labelled "prose"

Method: for every modified `.py` outside this scope, parsed HEAD and worktree, stripped every
docstring and (via `ast.unparse`) every comment, and compared. Six files differ, and all six
differ ONLY in human-readable message strings:
- `_permutation.py`, `battery.py`, `test_checks_c3_bridge.py`, `test_r_cross_implementation.py`
  - error text, a `pytest.skip` message and an assertion message. No test matches on the old
  wording (`grep` for the removed phrases across `tests/ src/ jobs/` is empty). Inert.
- `battery.py`'s new counts check out against the code: `ROW_KEYS` has 9 rows, 6 of them
  `CALIB_PERMUTATION` ("Six of its nine"), and `include_tau_checks=False` leaves 4
  (pinned by tests/test_checks_statistics.py:312), i.e. drops 5. Correct.

VERIFIED (this is the class that sank patch 3) | jobs/bench/report.py, analyzers/
instrument_validation.py | The two changed strings are the ones that get WRITTEN into committed
artifacts, and this time code and artifact AGREE: `report.py` now emits
"# Promotion report - the six checks against the calibration bench" and "**Five checks are
assessed, not six.**", which is verbatim `jobs/bench/results/promotion_report.md` lines 1 and 5;
`instrument_validation.py` now emits "`promotion_report.md`, which scores five checks - C3 has
no bench cell.", verbatim `jobs/bench/results/instrument_report.md:23`. Both artifacts are
tracked and unmodified, so these edits close a code/artifact divergence rather than opening one.

MINOR | tests/test_transform_guards.py:202-207,219-220,230-232,240 | Unit-suffix rule: the new
`timestamps` local (and the returned `values`, which are Hz) carry no unit suffix, and
`_bimodal` builds them as `np.arange(...) * 24.0` seconds. Mirrors `tlf.run`'s existing
parameter names, so it is consistent rather than novel, but it is new code under a hard rule. |
`t_rel_s` / `values_hz` in the fixture, keeping the call keyword-free as now.

## Patch 4 verdict: see the final report. Behaviour of all 6 listed changes is sound; the
## defects are in the GUARDS (one test that cannot fail, one guard defeatable by naming, one
## ladder line that plots min-over-calibrations under a bare "p-value" label).
<!-- PATCH4-BEHAVIOUR-FINDINGS-END -->

--- PATCH COMMIT 4 of 4 — PROSE (comments/docstrings/markdown only), appended 2026-09-02 ---

SCOPE: prose in the ADDED lines of: AGENTS.md, README.md, CONTRIBUTING.md, CLAUDE.md (new
symlink), docs/PANEL_CONTRACT.md, docs/TIME_SEMANTICS.md, docs/FIGURE_STANDARD.md,
docs/WRITING_A_JOB.md, docs/iid_checks/*, docs/iid_checks/CvM_cramer_von_mises.md (NEW),
the chronology purge across src/**, jobs/**, tests/**, scripts/**, conftest.py,
pyproject.toml, and src/quebra/analyzers/checks/battery.py's COUNTS.
COMMIT: ff4a8c5 (working tree, uncommitted)
NOTE: Makefile, .gitignore and .github/** are UNCHANGED in this diff despite being named in
scope; nothing to review there. CITATION.cff / quebra.toml / this ledger skipped per scope.

## Manifest (patch 4 prose)
- [ ] CLAUDE.md (new symlink)
- [ ] AGENTS.md
- [ ] README.md
- [ ] CONTRIBUTING.md
REVIEWED battery.py
- [ ] docs/PANEL_CONTRACT.md
- [ ] docs/FIGURE_STANDARD.md
- [ ] docs/TIME_SEMANTICS.md
- [ ] docs/WRITING_A_JOB.md
- [ ] docs/iid_checks/iid_checks_basics.md
- [ ] docs/iid_checks/C2_anderson_darling.md
- [ ] docs/iid_checks/C6_exchangeability.md
- [ ] docs/iid_checks/LIMITATIONS.md
- [ ] docs/iid_checks/CvM_cramer_von_mises.md (NEW)
- [ ] src/** chronology purge (grouped)
- [ ] jobs/** chronology purge (grouped)
- [ ] tests/** chronology purge (grouped)
- [ ] scripts/**, conftest.py, pyproject.toml (grouped)

## Reviewed (patch 4 prose)

## Findings (patch 4 prose)

### CLAUDE.md (new)
VERIFIED | CLAUDE.md | `readlink CLAUDE.md` -> `AGENTS.md`. It is a relative symlink to
AGENTS.md, exactly as AGENTS.md:14 now claims. No content of its own. Nothing to trim.

### AGENTS.md
VERIFIED | AGENTS.md:14 | "`CLAUDE.md` is a symlink to it so it loads every session" - true,
and AGENTS.md is the only agent-facing file (`find . -name 'AGENTS.md' -o -name 'CLAUDE.md'`
returns only these two).
VERIFIED | AGENTS.md:15 | `spec/quebraplan.md` exists. Old `spec/PLAN.md` does not.
VERIFIED | AGENTS.md:25,28 | `data/real_private/` exists in the checkout and holds
`6D2S/`, `2x2/`, `companion/`, `calibration_logs/`, `MANIFEST.toml`.
FALSE | AGENTS.md:29 | "`data/simulated/` holds regenerable payloads" - the directory is
EMPTY on disk and `git ls-files data/` returns exactly one path (`real_private/MANIFEST.toml`).
Present tense claim about content that does not exist. Also silently drops `data/real_public/`,
which DOES exist, so the sentence describes two thirds of a three-way split. | "The embargoed
records live under `data/real_private/`, inside the checkout. `data/real_public/` and
`data/simulated/` are the redistributable and regenerable halves."
VERIFIED | AGENTS.md:69 | `analyzers/within_calibration_compute.build_within_calibration_panel_data`
- module exists, `def build_within_calibration_panel_data` at line 395. Old
`panels/_within_calibration_compute` is gone.
VERIFIED | AGENTS.md:193-197 | Layout block: `analyzers/within_calibration_compute.py`,
`analyzers/within_calibration_data.py`, `panels/_within_calibration_render.py`,
`panels/within_calibration.py` all exist; `panels/_within_calibration_compute.py` and
`panels/_within_calibration_data.py` do not.
VERIFIED | AGENTS.md:231-232 | `scripts/` holds exactly `acceptance.sh`, `promote_run.py`,
`make_data_manifest.py`, `make_fixtures.py`, `make_job_manifest.py`. List is complete and
correct.
VERIFIED | AGENTS.md:234-235 | provenance note: `_git` runs with `cwd=Path.cwd()` and both
`get_git_commit` and `is_tree_clean` go through it, so the cwd-anchoring claim is true.
VERIFIED | AGENTS.md:252-253 | `data/real_private/6D2S/` and `data/real_private/companion/`
both exist.
FALSE | AGENTS.md:159 | "Fourteen lines remain in `*.py` for that reason." Fourteen is the
count at HEAD. The shipping tree has FIFTEEN (`grep -rn repairable --include='*.py'`), because
this same commit added the `# FROZEN VOCABULARY.` comment at
tests/test_artifact_guard.py:297. The number was counted against the wrong tree state. It is
also a liability by the author's own rule (2): nothing regenerates or asserts it
(`grep -rn 'Fourteen' tests/ scripts/ Makefile` is empty), and it will be wrong again on the
next comment edit. | Delete the sentence. `panels/across_calibration.py` carries the canonical
note; a running tally of grep hits adds nothing a reader can act on.
IMPORTANT | AGENTS.md:189-190 | The `checks/` layout list is now incomplete in a way that
makes the count on the next line unreadable: it names "C1 ..., C2 ..., C3 copula-via-R ...,
C5 ..., C6 ..." and omits CvM, then line 191 says "the five permutation checks". A reader
counts the five listed and concludes C3 is one of them. C3 is not in `ROW_KEYS` at all and is
never run by the battery; CvM is. | Add ", CvM Cramer-von Mises" to the list after C2.
VERIFIED | AGENTS.md:191 | "the five permutation checks" - `ROW_KEYS` has six
CALIB_PERMUTATION rows across five distinct check names (c1, c2, c5, c6, cvm). "Five checks"
is right; "five rows" would have been wrong.
IMPORTANT | AGENTS.md:234 | Over-correction. The deleted note - "`provenance.py` reports the
tree clean while it changes" - was a live hazard for the reader this file addresses: an agent
writes untracked files constantly, and `is_tree_clean` passes `--untracked-files=no`, so a
tree carrying brand-new untracked modules still reports clean and artifact reuse still fires.
The replacement note answers a different question (which repo the helpers describe) and the
hazard now survives only in the `is_tree_clean` docstring. | Keep both: append "Untracked
files do not make the tree dirty, so reuse can fire on a tree holding new untracked code."
MINOR | AGENTS.md:154 | Chronology the purge missed, and the exact banned form: "The rename
landed on 2026-08-23 (SPEC 0001 R0.4)." Spec number plus requirement id plus date. | "The
rename is ours."
MINOR | AGENTS.md:171 | Chronology the purge missed: "SPEC 0002 moved the packages under
`src/quebra/`." | "The packages live under `src/quebra/`."

### README.md
FALSE | README.md:125,132 | "`--all` ... runs the ones declaring `JOB_SWEEP = True`" and the
inline comment "# every job declaring JOB_SWEEP = True". `discovery.py:126` reads
`sweep=True if sweep is None else sweep`, so a job that declares NOTHING is swept. The
condition is opt-OUT, not opt-in, and `discovery.py:167` says so in its own docstring
("everything not opting out via `JOB_SWEEP = False`"). Two occurrences. | comment: "# every
job that does not opt out"; prose: "runs every job that does not declare `JOB_SWEEP = False`
- composites opt out."
VERIFIED | README.md:130-131 | "`--all` discovers every job under `./jobs`" - cli.py:150-175
walks `jobs/` (not `jobs/active`), errors if `jobs/` is absent, and filters on `d.sweep`.
The old "sweeps `./jobs/active`" was the false claim; the directory half of the fix is right.
VERIFIED | README.md:72 | The stray `**` on its own line before "Identity is content" is
fixed and the bold now renders. Real defect, correctly repaired.
VERIFIED | README.md:87-89 | The duplicated "Feedback is welcome..." sentences are removed
from the run-on paragraph and survive once, below. Correct de-duplication.
SHOULD | README.md:113-116 | Four sentences to say one thing, and the last one argues with an
imagined reader. "`scripts/acceptance.sh` checks the same claim from a harsher angle: it
builds a wheel, installs THAT into a throwaway virtualenv outside the repository, and runs
the suite from a directory that is not the checkout. If a packaging mistake makes the steps
above work only from a git clone, that script is what catches it." The second sentence is the
first sentence with a justification attached - the recognisable shape. | "`scripts/
acceptance.sh` runs these steps against a built wheel, in a throwaway virtualenv outside the
repository, from a directory that is not the checkout." (Verified against the script below.)
MINOR | README.md:113 | "from a harsher angle" and the shouting "THAT" are register drift in
the one file aimed at strangers. The mechanism carries the point without either.

### CONTRIBUTING.md
VERIFIED | CONTRIBUTING.md:36-37 | `make check` is `check: lint types arch test` where
`lint` = ruff check + ruff format --check, `types` = mypy, `arch` = `lint-imports --no-cache`
(the import contract), `test` = `pytest $(FAST)`. `.github/workflows/ci.yml:43-47` runs
`make check` then `make deps`. "what CI runs" is accurate. `make deps` is `deptry .`.
VERIFIED | CONTRIBUTING.md:8 | Repository renamed qretool -> quebra consistently; matches
the `quebra` console script and package name.
MINOR | CONTRIBUTING.md:26 | "Response is best effort, be kind <3" replaces "Response is
best effort. This is thesis work, not a funded project." The deleted half was the REASON the
response is best effort, which is the only part a stranger opening an issue needs; the
emoticon is register drift in a contributor-facing contract. | "Response is best effort.
This is thesis work, not a funded project."

### src/quebra/analyzers/checks/battery.py
FALSE | battery.py:107-108 | NEW false claim, and the file contradicts itself 33 lines later.
"`include_tau_checks=False` drops the three asymptotic rows - C1, C2 and CvM - leaving the
six permutation-calibrated ones." The `if include_tau_checks:` block at battery.py:153-197
contains FIVE rows: c1-asymptotic, c1-PERMUTATION, c2-asymptotic, cvm-asymptotic and
c2-PERMUTATION. So the flag drops 5 of 9 and leaves 4 - c5-studentized, c5-raw, c6, and
cvm-by-permutation - not 6. battery.py:143 and battery.py:217-221 both state correctly that
only CvM's permutation row survives the flag. | "`include_tau_checks=False` drops C1 and C2
entirely, both calibrations, along with CvM's asymptotic row. Four rows survive: the three
rank rows and CvM by permutation."
FALSE | battery.py:9-10 | The row count was updated and the measurement it labels was not.
"sharing it cuts a full nine-row replicate at N = 355 from 362 ms to 197 ms, and the whole
6-point n sweep from 613 ms to 341 ms." 362/197/613/341 ms were measured on a SEVEN-row
replicate; relabelling them "nine-row" makes four numbers describe a configuration nothing
ever timed. Rule (2) exactly. | Drop the milliseconds and keep the ratio-free claim:
"sharing it cuts a full replicate at N = 355 and the 6-point n sweep to roughly half."
Re-measure only if the figures are wanted.
VERIFIED | battery.py:100 | "Up to nine results in `ROW_KEYS` order." `ROW_KEYS` has 9
entries. Correct.
VERIFIED | battery.py:106 | "the other eight rows are still wanted" - dropping
c2-asymptotic leaves 8 of 9. Correct.
VERIFIED | battery.py:128 | "Six of its nine rows are permutation-calibrated" - CALIB_
PERMUTATION appears 6 times in `ROW_KEYS` (c1, c2, c5-stud, c5-raw, c6, cvm). Correct.
MINOR | battery.py:8-9 | Adjacent to an edited line and now false by omission: "C1, C2,
C5-studentized, C5-raw and C6 all want the identical one". CvM wants it too - battery.py:144
passes `permuted=permuted` to `cvm.run` - so six things share the matrix, not five. | Add
"and CvM" to the list.
MINOR | battery.py:108-109 | The rewrite left a ragged wrap with an orphan "It" ending line
108 and "is for records where" starting 109. | Re-wrap the paragraph.
VERIFIED | battery.py:143 | "so this is not gated on that flag alone" replaces "can no longer
be gated on that flag alone" - chronology removed, present tense true, no information lost.

### tests/test_windows_not_interpolated.py  [reviewed]  - CHANGE 1

IMPORTANT | tests/test_windows_not_interpolated.py:74 | The blacklist is NOT sufficient compensation for dropping the stem requirement, and the gap includes this repo's OWN interpolator. Measured old-vs-new over a 45-name corpus: the rule went from "suffix must be a dataset stem" to "any suffix", and the token list catches only 5 of the 20 plausible resampling names I could construct. ACCEPTED BY THE NEW RULE, REJECTED BY THE OLD: `filter_pchip`, `filter_spline` (transforms/interpolate.py:63 literally calls `scipy.interpolate.pchip_interpolate`, so `pchip` is the most likely name a second resampling recipe would carry), `filter_downsample`, `filter_upsampled`, `windows_downsampled` (the list has `resample` but not `sample`, so two thirds of the standard verb family escape), `filter_binned` / `filter_bin_10ms` (`rebin` only matches `rebin`), `filter_grid` / `filter_gridded` / `filter_onto_grid` (`regrid` only matches `regrid`), `filter_smooth`, `filter_nearest`, `filter_asfreq`, `filter_reindexed`, `filter_zoh`, `filter_evenly_spaced`, `filter_rolling_mean`, `filter_lowpass`, `filter_ffill`, `filter_imputed`. | Minimal: `RESAMPLING_TOKENS = ("interp","sample","grid","decimate","rebin","bin","pchip","spline","smooth","uniform")` - note `sample` subsumes `resample`/`upsample`/`downsample`/`subsample` and `grid` subsumes `regrid`. `bin` and `fill` are unsafe as bare substrings (`combined`, `binary` contain `bin`), so if they are wanted, split on `_` and match TOKENS rather than substrings.

IMPORTANT | tests/test_windows_not_interpolated.py:84-88 | A name-independent check is available and was not used, which is what makes the loosening avoidable rather than necessary. Measured: `job.dag[n].fn.__qualname__` is `_interpolate_step.<locals>.step` for the interpolate node and `_filter_step.<locals>.step` for the filter node - the closure factories ARE distinguishable, contradicting the premise at :117-119 (which is only true of `fn.__name__`). So the one resampling step that exists in this repo can be caught by identity regardless of what a job calls the node, and the name blacklist would only need to cover FUTURE resamplers. | Add, beside the token check: `assert not any(job.dag[a].fn.__qualname__.startswith("_interpolate_step") for a in _ancestors(...))`, or more generally reject any ancestor whose `fn.__module__`/qualname reaches `quebra.transforms.interpolate`. Keep the token list as the secondary net.

MINOR | tests/test_windows_not_interpolated.py:86-88 vs :6-10,144 | The token check runs BEFORE the whitelist and `continue`s, so the whitelist is no longer primary and the documented escape hatch does not work. A legitimate step whose name contains a token can NOT be fixed by "add it to OBSERVED_ANCESTORS deliberately" (the failure message at :144) - the entry would be dead, which is exactly why `fidelity_interp` had to be deleted from the whitelist in this diff. C1's null is stated as "events uniform on [0, tau]" (c1_lewis_robinson.py:19), so `uniform` is the plausible false positive: a step named `filter_uniform_null` fails unfixably. Fails loudly and in the safe direction, hence MINOR. | Either check the whitelist first and the tokens only on the non-whitelisted remainder, or amend the failure message to say a token-bearing name must be RENAMED, not whitelisted.

MINOR | tests/test_windows_not_interpolated.py:111-112 | Dead code. The `while` exits either by `break` (which happens only when `base in OBSERVED_ANCESTORS`) or by exhausting the condition (which runs the `else` and `continue`s), so line 111 is reached only with `base` already whitelisted and the `if` is never true. | Delete :111-112.

VERIFIED | tests/test_windows_not_interpolated.py:103-110 | The new prefix loop terminates and is correct on the degenerate names. Each iteration strictly shortens `base` (`rpartition` on a present `_`), so no infinite loop. Measured: no underscore (`foo`) -> `while` body never runs -> `else` -> rejected; `"_"` -> base `""` -> rejected; `"___"` -> three iterations -> rejected; leading underscore `_filter` -> base `""` -> rejected (same as old). One trivial loosening: trailing underscore `filter_` is now ACCEPTED (old rejected it, because `suffix` was empty and the `while` guard failed). Harmless.

VERIFIED | tests/test_windows_not_interpolated.py:58-67 | `discover` really skips every non-job file. Only 12 files under `jobs/` declare `JOB_ID` (9 in `active`, 3 in `composite`); all 8 `jobs/bench/*.py` and the 3 `__init__.py` are skipped, and `discover` skips `__init__.py` by name before parsing (discovery.py:137). Parametrised cases: 9 BEFORE -> 13 NOW (adds `check_ledger_q1`, `compare_t2star_0704_vs_1004`, `independence_survey`); file total 12 -> 16 tests, all passing in 1.7 s. New coupling, acceptable: a SyntaxError in any `jobs/**.py` now raises `UnreadableJobFile` during COLLECTION of this file, where before only `jobs/active` could break it.

VERIFIED | tests/test_windows_not_interpolated.py:136 | The 34 survey chains are genuinely reached, not nominally. `independence_survey` builds 206 nodes, `_windows_nodes` finds 34 (`windows_q1_040423` ...), and their union of ancestors is 170 nodes spanning `load_*`, `filter_*`, `final_*`, `t2star_*`, `windows_*`. The chain is `load -> filter_{label} -> final_{label} -> t2star_{label} -> windows_{label}` (independence_survey.py:242-248), and `_final_stage` (recipes.py) only unwraps `FilterResult.final_norm` - no resampling, so `final` is a correct whitelist addition. Under the OLD rule 3 x 34 = 102 of those names were unmatchable, so the widening genuinely could not have landed without loosening something.

VERIFIED | tests/test_windows_not_interpolated.py:42 | Deleting `fidelity_interp` from the whitelist is inert AND forced. `grep -rn fidelity_interp` over the whole repo returns nothing, so no job names such a node; and with the token check placed first the entry could never have matched anyway. It was also a rule violation while it lasted (a whitelisted "observed, never resampled" entry whose name says otherwise). `ramsey_q1_100423`'s `fidelity_windows` ancestors are `{fidelity_raw, filter, final_filter_stage, load, load_df, lookup_prior}` - `interpolate` is NOT among them, confirmed by construction.

VERIFIED | tests/test_windows_not_interpolated.py:190-209 | The new positive control does fail on the hole it names: `_unexpected` returns `{"filter_then_interpolate"}` with the token list present, and returns `set()` (test fails) with `RESAMPLING_TOKENS = ()`. It is not decoration.

## Reviewed (patch 4)
- [x] tests/test_windows_not_interpolated.py

### docs/FIGURE_STANDARD.md
VERIFIED | FIGURE_STANDARD.md:6-8 | Axis-label attribution is now correct.
`panels/comparison.py:34` (`x_label: str = "Elapsed time (h)"`) and
`plots/interpolation_stage_plot.py:107,129` (`set_xlabel("Elapsed time (h)")`) are the two
places with that label; `panels/across_calibration.py:191` has
`set_ylabel("Inter-event interval (h)")`. `panels/within_calibration.py` has neither, so
dropping it from the list was right. Both labels are ruled out by the table at
FIGURE_STANDARD.md:49 (`scan clock | wall time, elapsed time`) and :45
(`window | interval, period, span`).
SHOULD | FIGURE_STANDARD.md:8 | The inventory is presented as complete and is not:
`panels/across_calibration.py:190,223` also `set_xlabel("Elapsed time (days)")`. The new
precision about `(h)` is what creates the gap - it buys nothing and excludes a real
violation. | "`panels/comparison.py`, `plots/interpolation_stage_plot.py` and
`panels/across_calibration.py` all label an axis `Elapsed time`, and
`panels/across_calibration.py` labels one `Inter-event interval (h)`."
MINOR | FIGURE_STANDARD.md:9 | The rewrite broke the wrap: line 9 is 152 characters where
the rest of the file wraps near 90. Also `; and` splices two clauses that were separate
sentences. | Re-wrap and split at "rules out. No panel yet reports how much data it dropped."

### docs/PANEL_CONTRACT.md
VERIFIED | PANEL_CONTRACT.md:56-58 | The Kaplan-Meier rewrite is true on both halves.
`src/quebra/analyzers/kaplan_meier.py` exists; its only importers are
`plots/km_survival_plot.py` and `jobs/active/km_poster_6d2s.py`, and
`panels/within_calibration.py` does not import it - so "the estimator exists but this panel
does not consume it" holds. `ReliabilityBand` declares `cumulative_hazard`, `band_lower`,
`band_upper` as `None` defaults (reliability_band.py:70-72) and the only construction site
is `ReliabilityBand(estimator=ESTIMATOR_CRUDE)` at :187, so "are `None` here" holds.
VERIFIED | PANEL_CONTRACT.md:202-203 | "`analyzers/kaplan_meier.py` uses the censored windows
rather than discarding them" - kaplan_meier.py:10 "Right-censored windows are kept, not
dropped", and it carries `censor_time_min` / `censor_survival` / `n_censored`. True. The
crude path does discard them (reliability_band.py:218). Correct contrast.
VERIFIED | PANEL_CONTRACT.md:290 | `analyzers/fidelity.py::make_panel_data(result:
FidelityResult, windows, reads, ...)` exists at fidelity.py:222 and returns
`WithinCalibrationPanelData`. The old `plots/fidelity_plot.py` does not exist. Correct fix.
IMPORTANT | PANEL_CONTRACT.md:310-311 | Half-fixed sentence, so it still names a file that
does not exist. The compute path was corrected to `analyzers/within_calibration_compute.py`
but the same sentence still says "`_within_calibration_data.py` is the typed contract".
`src/quebra/panels/` holds only `within_calibration.py` and `_within_calibration_render.py`;
the typed contract is `analyzers/within_calibration_data.py`, as AGENTS.md:196 now says. |
"...and `analyzers/within_calibration_data.py` is the typed contract."
SHOULD | PANEL_CONTRACT.md:310 | "The split is partly done:" is status narration, the
author's rule (4), and it replaced status narration ("The split it proposed has since
happened in part") rather than removing it. The paragraph does not need a progress verb at
all. | "`analyzers/within_calibration_compute.py` builds the artifact, ..."
SHOULD | PANEL_CONTRACT.md:17-19 | The rewrite left a sentence that is half rule and half
obituary, and the middle clause stayed in the past: "`FidelityPlot` and `T2StarPlot` must
not build panel data inside `build_matplotlib`, where no DAG node could supply the window
tables; both are gone and jobs use `WithinCalibrationPanel` directly". If both classes are
gone, a prohibition addressed to them is unreadable. | "No plot builds panel data inside
`build_matplotlib`: no DAG node can supply the window tables there. Jobs use
`WithinCalibrationPanel` directly off a panel-data step."
OPTIONAL | PANEL_CONTRACT.md:57-58 | "wiring it in changes only that band" argues a
counterfactual nothing can check. The two verified sentences before it carry the point. |
Delete the clause.
MINOR | PANEL_CONTRACT.md:307 | Chronology the purge missed, in the paragraph it edited:
"The number in this paragraph read "~650" until 2026-08-23; the file had grown and the doc
had not." | Delete the sentence.
OPTIONAL | PANEL_CONTRACT.md:306 | Adjacent typed number, off by one: "is 1057 lines";
`wc -l` says 1056. Nothing regenerates or asserts it. "still the largest module in the tree"
IS true (next largest is `analyzers/instrument_validation.py` at 704). | "is over five times
the project's 200-line guideline" - drop the exact count.

### docs/TIME_SEMANTICS.md
VERIFIED | TIME_SEMANTICS.md:34,64,94 | All three dataset paths corrected correctly.
`data/real_private/companion/` and `data/real_private/6D2S/` exist and hold the named
patterns. No `.hdf5` exists anywhere in the checkout, so "outside the checkout" is true.
VERIFIED | TIME_SEMANTICS.md:38-39 | Smart quotes in the `python` block replaced with
straight quotes. That block was not valid Python before; real fix.
VERIFIED | TIME_SEMANTICS.md:61 | "There is no `t_s` key; downstream steps use `t_rel_s`."
`grep -rn '"t_s"' src/ jobs/ tests/` is empty. True, and the chronology ("no longer
emitted") is gone without losing the fact.
VERIFIED | TIME_SEMANTICS.md:128 | Deleting "(`_check_unix_s` is a deprecated alias and will
be removed.)" loses nothing: only `check_unix_s` exists
(`transforms/lookup_prior.py:16`); the underscore alias is gone.
SHOULD | TIME_SEMANTICS.md:77 | Over-correction that breaks the section's own contract. The
heading is under "## Verified Dataset Types" and lists "**Verified fields**", but
`frequency_cal_q1.hdf5` was replaced with "core-tools `.hdf5`, outside the checkout" - so
nobody can tell which file the verification was performed against, and the verification
becomes unreproducible. "core-tools" also now appears twice in one heading. | "### 2)
core-tools HDF5 calibration sweeps (`frequency_cal_q1.hdf5`, outside the checkout)"
OPTIONAL | TIME_SEMANTICS.md:34 and AGENTS.md:253 | `data/real_private/companion/` holds
exactly one file, `qubit1.pickle`. `qubit*.pickle` is an honest glob; AGENTS.md's
`qubit{N}.pickle` implies a family that is not there. | AGENTS.md: `qubit1.pickle`.

### docs/iid_checks/iid_checks_basics.md
VERIFIED | iid_checks_basics.md:13,103 | "all six" and "the six checks" -
`src/quebra/analyzers/checks/` holds c1_lewis_robinson, c2_anderson_darling,
c3_serial_copula, c5_rank_autocorr, c6_exchangeability, cvm_cramer_von_mises. Six.
VERIFIED | iid_checks_basics.md:104 | "runs the five permutation checks off ONE permutation
set" - `battery.ROW_KEYS` has six CALIB_PERMUTATION rows over five distinct check names
(c1, c2, c5, c6, cvm). Five checks is correct.
VERIFIED | iid_checks_basics.md:112-113 | "There is no C4. Lin-Wei-Ying is the fourth check
in the source's numbering and is not implemented here" -
`.claude/qre_checks_reference.tex:790` is titled "C4 -- Lin-Wei-Ying, and why it is not in
the battery", and :92 says it "is documented in Section~\ref{sec:lwy} and not implemented".
Both halves true.
VERIFIED | iid_checks_basics.md:110-111 | All six per-check pages linked exist on disk,
including the new `CvM_cramer_von_mises.md`.
SHOULD | iid_checks_basics.md:113 | "; nothing in the battery depends on it." A check that
is not implemented cannot be depended on, so the clause carries no information - it is the
closing-clause-that-argues shape. The reference explains WHY it is out (it tests a fitted
regression model, a different tier); that would be worth a clause, this is not. | "There is
no C4. Lin-Wei-Ying is the fourth check in the source's numbering; it tests a fitted
regression model, not the durations, and is not implemented here."

### docs/iid_checks/LIMITATIONS.md
VERIFIED | LIMITATIONS.md:1 | "all six checks" - correct, six modules on disk.

### docs/iid_checks/C6_exchangeability.md
VERIFIED | C6_exchangeability.md:13 | "the natural reference among the six" - correct.

### docs/iid_checks/C2_anderson_darling.md
SHOULD | C2_anderson_darling.md:30 | The chronology went, an empty assertion took its place:
"That distinction matters and is load-bearing - see the limitations below." "Matters" and
"is load-bearing" are the same claim twice, and neither says anything the bolded sentence
before it did not. | "**Both force `gamma = 1`, so neither pins the shipped path.** See the
limitations below."

### docs/iid_checks/CvM_cramer_von_mises.md (NEW, 73 lines)
Nothing on this page is invented. Every functional, function name, test, row and gate checks
out against the code. Three claims are wrong and one pointer is wrong.

VERIFIED | CvM:5 | "The third functional" is CORRECT and the CODE is wrong.
`.claude/qre_checks_reference.tex:278-281` tabulates the class as LR (4), Kolmogorov-Smirnov
(5), Cramer-von Mises (6), Anderson-Darling (7). CvM is third by both table order and
equation number. `cvm_cramer_von_mises.py:3` and `battery.py:69` both say "the FOURTH
functional"; those two are the errors. | Fix the two code comments to "third", not the doc.
SHOULD | CvM:5-13 | The ordinal and the code block contradict each other for the reader: the
page says "third" and then shows a three-line block in which CvM is listed FIRST and
Kolmogorov-Smirnov (the second of the four) is absent. | Add `sup_s |W0(s)|   <- Kolmogorov-
Smirnov, not implemented` to the block, or reorder it (4)(5)(6)(7).
FALSE | CvM:15-17 | The equation citation was dropped rather than corrected, so a transcribed
formula now cites no equation - the project's own statistical hard rule. "The reference
prints its middle term as `- i N (T^2_{i+1} - T^2_i)/tau`". The reference's CvM statistic is
eq (6) (`\label{eq:cvm}` at tex:589). `cvm_cramer_von_mises.py:16` cites "its eq (7) middle
term", and eq (7) is Anderson-Darling - so the module cites the wrong equation and the page
cites none. | "The reference prints eq (6)'s middle term as ...".
FALSE | CvM:67-68 | The pointer is wrong and LIMITATIONS.md disclaims it in its own second
line. "Shared with C1 and C2, and set out in `LIMITATIONS.md`: the asymptotic calibration
divides by an estimated `gamma_hat` ...". `grep -n gamma docs/iid_checks/LIMITATIONS.md` is
EMPTY; that file's eight sections are cross-cutting ones and it opens "Per-check limitations
live on the per-check pages." The gamma_hat caveat is on `C2_anderson_darling.md:34-38` (with
measured figures) and `C1_lewis_robinson.md:47`. | "Shared with C1 and C2, and stated on
their pages: the asymptotic calibration divides by an estimated `gamma_hat`, so ..."
SHOULD | CvM:49-51 | Mechanism misattributed, and the mechanism is the point of the
paragraph. "On the in-spec clock of a carved record `tau == T_N`, which makes eq (7) singular
and silences C1 and C2". Eq (7) is C2's alone; C1 is eq (4) and has no `1/(s(1-s))` weight.
`c1_lewis_robinson.py:120-155` has no `tau == T_N` guard and always returns a p-value - only
`c2_anderson_darling.py:89,147,165` declines. C1 goes dark because the CALLER turns off
`include_tau_checks`, not because its formula fails. | "On the in-spec clock of a carved
record `tau == T_N`. That makes eq (7) singular, so C2 declines outright, and it makes the
truncation time event-determined, so C1's asymptotic null is void too. CvM's integrand
carries no `1/(s(1-s))` weight, so its statistic stays finite there."
SHOULD | CvM:49 | "It also survives a case C2 cannot" is true only of the permutation row,
and the page's own next section says so. Stated unqualified four lines earlier it is the
kind of claim this commit exists to remove. | "Its permutation row also survives a case C2
cannot."
MINOR | CvM:55-56 | "The asymptotic row is gated with C1's and C2's - `include_tau_checks`
and a single segment". C1's asymptotic row is NOT gated on a single segment
(`battery.py:154-171` emits it for any m); only C2's and CvM's are. | "gated with C2's -
`include_tau_checks` and a single segment - and, like C1's, dropped when the tau gate is off".
VERIFIED | CvM:33-34 | `analyzers/checks/cvm_cramer_von_mises.py` defines `statistic` (:172),
`statistic_batch` (:191), `cvm_limiting_cdf` (:107) and `run` (:208). All four names correct.
VERIFIED | CvM:36-39 | All four pinned directions exist in `tests/test_checks_cvm.py`:
`test_at_gamma_one_it_is_the_classical_cramer_von_mises` (:34, against
`scipy.stats.cramervonmises`), `test_the_bracket_is_the_bridge_integral` (:43, via
`scipy.integrate.quad`), `test_the_limiting_cdf_reproduces_the_published_critical_values`
(:61), `test_statistic_batch_reproduces_the_statistic_under_the_identity` (:84).
VERIFIED | CvM:28-29 | The `gamma_hat = 1` identity claim matches the test exactly:
`cramervonmises(np.cumsum(seg.x) / seg.tau, "uniform")`, i.e. `u_i = T_i/tau`.
VERIFIED | CvM:43-45 | The `m > 1` reversal and the quotation are the source paper's, per
`.claude/qre_checks_reference.tex:598-600`: "for several processes we do not include the
Anderson-Darling test as the Cramer-von Mises test had better level properties in this case".
`C2_anderson_darling.md:53-55` says the same. Quote is accurate word for word.
MINOR | CvM:43 | Attribution slip on the sentence that carries the quote. "The reference
reverses its own single-process preference" - the reference document reports that the PAPER
reverses ("But for m > 1 the paper reverses the preference", tex:597). Everywhere else on
this page "the reference" means the project's reference document. | "The source paper
reverses its own single-process preference for `m > 1`".
VERIFIED | CvM:55 | "Two of the nine in `battery.ROW_KEYS`: one asymptotic, one
permutation-calibrated" - `ROW_KEYS` entries at battery.py:64 and :76. Correct.
VERIFIED | CvM:69-70 | "`test_tier3_calibration.py` measures that level and asserts only a
wide bound, deliberately" - :212-226 parametrizes `["c1", "c2", "cvm"]` and asserts
`0.005 <= size <= 0.20`; the module docstring at :11-13 says no direction is asserted.
VERIFIED | CvM:73 | "`run_battery` raises rather than defaulting one" - battery.py:127.
SHOULD | CvM:61-63 | A hazard is stated and never resolved, so the page reads as a warning
about itself. "`ROW_KEYS` is the schema of the bench tables, so registering a check without
re-running the bench makes `bench_acceptance_at_n` return None and every ledger row read
`underpowered / no bench cell`." The reader is not told this was avoided. The bench WAS
re-run: `jobs/bench/results/size_table.csv` holds 219 CvM rows (counted; the figure in
`cvm_cramer_von_mises.py:38` is right). Also `underpowered / no bench cell` is backticked as
a literal and is a paraphrase - the code has `VERDICT_UNDERPOWERED = "underpowered"`
(check_ledger.py:56) and `"no bench cell for this check at this event count"` (:252). |
"`ROW_KEYS` is the schema of the bench tables, so a check must be registered and the bench
re-run together or every ledger row reads underpowered. `size_table.csv` carries CvM cells."
STRUCTURE | CvM | Diverges from its five siblings. C1, C2, C3, C5 and C6 all run
`The equation`/`What it computes` -> `Where it lives` -> `Limitations` -> `What we do`. This
page has no `What we do` section - the one section that tells a reader how to USE the row -
and adds two the siblings do not have (`Why it is in the battery`, `Its rows`). At 73 lines
it is longer than every sibling (48-68). | Fold `Why it is in the battery` into two sentences
under `The equation`, keep `Its rows` (it is genuinely new information), and add `What we do`.
That lands near 60 lines.
NOTE | CvM | Dash-asides: 3 constructions (CvM:46, :56 pair, :72 pair) in ~40 prose lines,
against 22 in 58 last pass. The shape is much reduced but not gone. Does not read as
generated: the derivation, the quote and the gating are all specific and all check out.

### CORRECTION - the file changed UNDER the review at 17:30:57

The findings above were written against the version of `tests/test_windows_not_interpolated.py`
that was on disk when this pass started (prefix-only match, `RESAMPLING_TOKENS` as the gate).
That version was replaced mid-review by a `_carve_labels` version that RESTORES a suffix
requirement. The four IMPORTANT/MINOR items above about the blacklist being the only gate no
longer apply. Re-reviewed below. (Working tree is being edited while under review - re-run
the manifest before trusting any part of this section.)

### tests/test_windows_not_interpolated.py  [re-reviewed, _carve_labels version]  - CHANGE 1

IMPORTANT | tests/test_windows_not_interpolated.py:71-95 | The docstring's claim that taking
labels from the carve nodes "gives up nothing" is false, and the names it gives up are the
exact ones the new control cites. The anchor moved from EXTERNAL (a dataset path the job
loads) to INTERNAL (a name the same author gave a sibling node), so a self-vouching PAIR now
passes. MEASURED: a job with `filter_pchip` -> `windows_pchip` gives `_carve_labels == {"pchip"}`
and `_unexpected == set()`; same for `filter_smoothed`/`windows_smoothed` and
`filter_downsampled`/`windows_downsampled`. Under the OLD rule all three were REJECTED, because
`pchip`/`smoothed`/`downsampled` are not stems of any dataset the job loads. The tokens do not
save it: the list has `resample` but not `sample`, so `downsampled` passes, and
`transforms/interpolate.py:63` calls `scipy.interpolate.pchip_interpolate`, so `pchip` is the
likeliest name a second resampler would carry. The control at :219-231 only proves those three
names fail when the carve is named `windows_x`; rename the carve to match and they pass. | Cheapest real fix, and it removes the whole naming game: check callable IDENTITY, not names. `job.dag[n].fn.__qualname__` is `_interpolate_step.<locals>.step` for the resampler and `_filter_step.<locals>.step` for the filter - MEASURED, so the closure factories ARE distinguishable and the premise at :144-146 (`fn.__name__` is always `step`) does not extend to `__qualname__`. The whole `independence_survey` graph, 206 nodes, uses only 8 distinct callables. Add `assert not any(job.dag[a].fn.__qualname__.startswith("_interpolate_step") for a in _ancestors(job, node_id))` and the name rules become a secondary net rather than the guard.

MINOR | tests/test_windows_not_interpolated.py:70-73 | `RESAMPLING_TOKENS` is completely unpinned: setting it to `()` and running all 16 tests in this file gives 16 PASSES. MEASURED. So the second layer can be deleted with green gates, and the control at :207 that reads as if it pins the token list actually passes for the suffix reason (`then_interpolate` is not a carve label). | One assert on the tokens alone: a job with carve `windows_q1` and step `filter_q1_interp` (suffix IS built from a real label, so only the token list can reject it) must be flagged.

MINOR | tests/test_windows_not_interpolated.py:97-113 | The label set is a UNION over the whole job, so a label minted by dataset A's carve vouches for a step in dataset B's chain. In a 34-iteration loop the classic bug is exactly that - `windows_q5_070623` carved from `filter_q1_040423` - and neither the old nor the new rule can see it. MEASURED clean today: for every one of the survey's 34 carves, no ancestor name ends in a different chain's label (0 cross-wired ancestors). Worth pinning now that the survey is in scope. | Per-carve rather than per-job labels: derive the label from THIS carve node's own name and require every suffixed ancestor of it to carry that same label.

VERIFIED | tests/test_windows_not_interpolated.py:97-113 | Equal-strength on every pre-existing job, so the widening costs nothing on today's graphs. MEASURED label sets: `km_poster_6d2s` -> exactly its 5 dataset stems, i.e. identical to the old `_dataset_stems`; every other `jobs/active` job and both other composites -> the EMPTY set, so their ancestors must be exact whitelist entries as before; `independence_survey` -> its 34 readable labels. `fidelity_windows` and bare `windows` contribute no label (no `{known}_` prefix matches), so no job accidentally widens itself.

VERIFIED | tests/test_windows_not_interpolated.py:216-231 | The new control DOES kill the prefix-only mutant. Replacing `_unexpected` with the pure-prefix version and re-running the three controls fails `test_the_check_catches_a_resample_hidden_behind_a_whitelisted_prefix`. So the suffix half of the rule is genuinely pinned. `_carve_labels` calling `_windows_nodes`, which is defined later in the file, is fine (resolved at call time).

## Reviewed (patch 4)
- [x] tests/test_windows_not_interpolated.py (re-reviewed after mid-review rewrite)

--- RESUME (prose pass, 2026-09-02 late) ---
NOTE ON RESUME: the previous prose pass wrote per-file findings for README.md, CONTRIBUTING.md,
battery.py, PANEL_CONTRACT.md, FIGURE_STANDARD.md, TIME_SEMANTICS.md, iid_checks_basics.md,
LIMITATIONS.md, C6, C2 and the new CvM page, but died before ticking the manifest. All those
files are UNCHANGED on disk since that pass (mtimes 17:17-17:25, ledger 17:37), so those
findings stand and are NOT redone. Only AGENTS.md moved (22:37, the three fixes). Manifest
below updated accordingly; remaining real work is AGENTS.md past 191, WRITING_A_JOB.md and the
chronology purge groups.

## Manifest (patch 4 prose) - UPDATED
- [x] CLAUDE.md (new symlink)
- [x] AGENTS.md (three findings fixed + rest swept, this resume)
- [x] README.md
- [x] CONTRIBUTING.md
- [x] src/quebra/analyzers/checks/battery.py
- [x] docs/PANEL_CONTRACT.md
- [x] docs/FIGURE_STANDARD.md
- [x] docs/TIME_SEMANTICS.md
- [ ] docs/WRITING_A_JOB.md
- [x] docs/iid_checks/iid_checks_basics.md
- [x] docs/iid_checks/C2_anderson_darling.md
- [x] docs/iid_checks/C6_exchangeability.md
- [x] docs/iid_checks/LIMITATIONS.md
- [x] docs/iid_checks/CvM_cramer_von_mises.md (NEW)
- [ ] src/** chronology purge (grouped)
- [ ] jobs/** chronology purge (grouped)
- [ ] tests/** chronology purge (grouped)
- [ ] scripts/**, conftest.py, pyproject.toml (grouped)

### AGENTS.md - re-verification of the three fixes, and the sweep past 191

VERIFIED | AGENTS.md:25-30 | Fix 1 landed and is true. `data/real_private/`, `data/real_public/`
and `data/simulated/` all exist; the new sentence asserts no contents, so the empty-directory
problem is gone. `real_private/` holding the embargoed records is true (6D2S/, 2x2/,
companion/, calibration_logs/, MANIFEST.toml).
VERIFIED | AGENTS.md:159 | Fix 2 landed: the tally sentence is deleted outright. HEAD said
"Twelve", the mid-review tree said "Fourteen"; neither is now asserted.
VERIFIED | AGENTS.md:189-192 | Fix 3 landed and is true on both halves. CvM is named; "C3
copula-via-R (needs Rscript; not in ROW_KEYS, so the battery never runs it)" matches
`battery.ROW_KEYS` (no c3 entry, no c3 import); "the five permutation checks" matches 6
CALIB_PERMUTATION rows over 5 distinct names (c1, c2, c5, c6, cvm).

FALSE | AGENTS.md:53-55 | Stale claim, now false in both halves, on the page that opens by
telling the reader not to describe unshipped things as shipping. "`make types` and `make arch`
exist in the Makefile but do not pass yet: they depend on the `src/` layout (SPEC 0002) and the
import-linter contract (SPEC 0004). A failure from those two tells you which phase you are in."
MEASURED: `make arch` -> "Layered architecture KEPT. Contracts: 1 kept, 0 broken."; `make types`
-> "Success: no issues found in 11 source files". Both pass. Also two spec numbers and a phase
pointer, rule (1). | "`make lint`, `make test` and `make clean` wrap the first two. `make types`
runs mypy over `core/`; `make arch` runs the import-linter contract. `make check` is all four."
FALSE | AGENTS.md:279 | Inventory presented as complete and is 3 of 7. "`docs/` holds reference
docs: `TIME_SEMANTICS`, `PANEL_CONTRACT`, `FIGURE_STANDARD`." On disk: also `JOBS.md`,
`WRITING_A_JOB.md`, `WRITING_A_SCHEMA.md` and `iid_checks/` (8 files, one added by this commit).
A reader told to "refresh docs" cannot refresh what the list hides. | "`docs/` holds the
reference docs, including `TIME_SEMANTICS`, `PANEL_CONTRACT`, `FIGURE_STANDARD`,
`WRITING_A_JOB` and `iid_checks/`."
MINOR | AGENTS.md:155-156 | Chronology the purge missed, exactly the banned form (date + spec +
requirement id): "The rename landed on 2026-08-23 (SPEC 0001 R0.4)." | "The rename is ours."
MINOR | AGENTS.md:172 | Chronology missed: "SPEC 0002 moved the packages under `src/quebra/`." |
"The packages live under `src/quebra/`."
MINOR | AGENTS.md:206-208 | Chronology missed, and it is the useful kind badly framed: "Moved
out of jobs/ in SPEC 0002: it is reusable library code, and leaving it in jobs/ forced the CLI
to put the caller's directory on sys.path, which R1.1.4 forbids." The hazard is worth keeping;
the move and the requirement id are not. | "Library code, not a job: a job here would force the
CLI to put the caller's directory on sys.path."
MINOR | AGENTS.md:216 | Chronology missed: "tests/  TRACKED since 2026-08-23 (SPEC 0001 R0.1)" |
"tests/  TRACKED"
MINOR | AGENTS.md:281-282 | Chronology missed: "The blanket `*.md` and `.*` ignore rules were
removed in SPEC 0001 R0.1, so `.md` needs no `!` exception." | "Nothing in `.gitignore` blanks
`.md`, so it needs no `!` exception."
IMPORTANT | AGENTS.md:234-235 | Over-correction, unfixed from the earlier pass and repeated here
because it is the one that costs a reader something. The deleted note - "`provenance.py` reports
the tree clean while it changes" - warned about a live hazard for exactly this reader: an agent
writes untracked files constantly and `is_tree_clean` passes `--untracked-files=no`, so a tree
carrying brand-new untracked modules reports clean and artifact reuse still fires. The
replacement answers a different question (which repo the helpers describe). | Keep the new
sentence and append: "Untracked files do not make the tree dirty, so reuse can fire on a tree
holding new untracked code."
OPTIONAL | AGENTS.md:262 | "The six-tier layout lives in `spec/quebraplan.md` Phase 5" - a phase
number survives, on a line this commit edited. | "...lives in `spec/quebraplan.md` and is not
built".

### src/** chronology purge  [reviewed]

MUST | src/quebra/schemas/ramsey_series.py:3-4 | LOGICAL INVERSION introduced by the rewrite:
the sentence now says the opposite of the truth. "Keeping this out of
`quebra.core.job._load_dataset` is what keeps the generic loading path Ramsey-specific".
Keeping it out is what keeps the path GENERIC. The old line ("used to live inside
`_load_dataset`, which made the generic loading path Ramsey-specific") had the subject the
predicate belonged to; the purge swapped the subject and left the predicate. The three
following clauses are also stranded in the past ("required", "was involved", "could not
use"). | "Keeping this out of `quebra.core.job._load_dataset` is what keeps the generic
loading path generic: inside it, `job.load` would require `timestamp` and `frequency`
columns and a resolvable run-start time whether or not a schema was involved, so a job with
nothing to do with frequency tracking could not use `job.load` at all."
MUST | src/quebra/core/discovery.py:43-45 | Botched rewrite: the replacement line duplicates
the tail of the line it replaced, so the comment reads "and dropping the distinction entirely
/ dropped the distinction entirely". Also two fragments in a row ("Not the DIRECTORY (...):"
then "A declaration rather than a directory") saying the same thing twice, and the surviving
half no longer says WHOSE mistake widened `--all`. | "# Declared, not encoded in the DIRECTORY
(`jobs/active` swept, `jobs/composite` not). Dropping the distinction entirely silently
widens `--all` from 9 jobs to 12." (9 -> 12 VERIFIED: 9 `JOB_ID` files in `jobs/active`, 3 in
`jobs/composite`, and only the composites declare `JOB_SWEEP = False`.)
SHOULD | src/quebra/core/closure.py:3-6 | Present-tense rewrite of a measured fact now reads
as a claim about CURRENT behaviour, and the behaviour it claims is the one this module exists
to prevent: "appending a line to `analyzers/t2star.py` leaves `t2star_q1_070423`'s identity
byte-identical." It does not, as of this module. The mixed tense in the same sentence
("would let ... while the digest asserted nothing had") is the tell. | "Identity over
`hash(job file) + dataset hashes + child identities` alone would let the analyzer that
produced the numbers change completely while the digest asserted nothing had: appending a line
to `analyzers/t2star.py` would leave `t2star_q1_070423`'s identity byte-identical."
SHOULD | src/quebra/core/closure.py:248-249 | Same tense break, plus a ragged wrap the rewrite
left behind (line 248 ends "two jobs differing", 30 characters short of the paragraph). "Without
it, two jobs differing only by `x=1` versus `x=999` both produced `93f7286634eef03b`." Also
`93f7286634eef03b` is asserted by nothing: `grep` finds it in this docstring, in
`tests/test_identity_closure.py:155`'s DOCSTRING, and in `spec/`. Rule (2). | "Without it, two
jobs differing only by `x=1` versus `x=999` collide on one digest." Re-wrap.
SHOULD | src/quebra/core/dataset.py:17-21 | The rewrite made the first clause hypothetical and
left the whole body in the past, so a reader cannot tell whether the crash is a live defect or
an averted one: "Merging them is a defect rather than a simplification: `extra` was forwarded
to the loader AND read back as metadata, so ... reached `pd.read_csv(path, **extra)` and raised
`unexpected keyword argument`. It only ever worked because every in-repo caller used
`.pickle`". Also a ragged wrap at :17. | "Merging them would be a defect rather than a
simplification: `extra` would be forwarded to the loader AND read back as metadata, so
`extra={'run_start_unix_s': ...}` - the value `_load_dataset`'s own error message tells you to
set - would reach `pd.read_csv(path, **extra)` and raise `unexpected keyword argument`, and
only `.pickle` callers would escape it, because that loader discards the dict."
SHOULD | src/quebra/core/discovery.py:8-9 | Grammatically broken conditional: "**`JOB_FAMILY`**
- the category. A directory would mean recategorising costs moving a file, which under logical
IDs must stop being meaningful." | "A directory would make recategorising a file move, which
under logical IDs means nothing."
SHOULD | src/quebra/core/runner.py:49 | Over-correction that invents a design intent. "the
value model is deliberately untyped" replaces "stays untyped this phase". Nothing supports
"deliberately"; the original said the opposite (a temporary state). | "the value model is
untyped".
SHOULD | src/quebra/analyzers/checks/_multiprocess.py:305-307 | The new comment gives the WRONG
reason for the guard's scope, and the right reason is checkable. "the guard below is in-spec
only, because what the calendar clock's tau should be when a block ends at a gap is an open
question." The guard's stated hazard is corrupted EVENT TIMES, not tau: on the calendar clock
`x = np.diff(t_birth)` (:315), so events come from birth times and an interior censored window
cannot shift them. Tau ambiguity is a different issue and does not explain the scope. |
"# Not applicable on this clock: `x` comes from birth times, so an interior censored window
cannot shift the event times."
SHOULD | src/quebra/analyzers/independence_survey.py:256 | The purge deleted the numbers'
provenance and kept the numbers, which is the wrong half. "Measured on the shipped artifact
before this fix: the 1 us and 10 us columns were 34/34 absent" became "Without this the 1 us
and 10 us columns are 34/34 absent in EVERY one of the seven figures - 132 of 340 cells per
grid silently blank." Now four measured figures hang off a counterfactual nothing can rerun.
Keep the word that anchors them. | "Measured: without the blanket rows the 1 us and 10 us
columns come back 34/34 absent in every one of the seven figures - 132 of 340 cells per grid
silently blank."
SHOULD | src/quebra/panels/within_calibration.py:7-9 | Chronology the purge MISSED, in the
paragraph directly above the one it rewrote, and it is the banned form (a bare date): "The
re-export was load-bearing until 2026-08-23: artifacts materialized before the band split
name `panels.non_repairable.NonRepairablePanelData` in their pickle stream, and output/ is
append-only, so dropping it would have broken reloading them." With the following paragraph
rewritten to the present, the two now contradict each other: para 1 says the re-export is
load-bearing, para 2 says nothing can reach those artifacts. | Delete para 1 and open with:
"Artifacts materialized before the band split name
`panels.non_repairable.NonRepairablePanelData` in their pickle stream. That module does not
exist, so they raise ModuleNotFoundError and no re-export here can reach them: the module they
name is gone, not the symbol."
MINOR | src/quebra/panels/within_calibration.py:11 | Information dropped that a reader acting
on this docstring wanted: "82 pickles under output/ and the backups are affected." The
paragraph tells the reader to recover by re-running jobs but no longer says how much there is
to re-run. Defensible under rule (2) (nothing asserts 82), so MINOR - but if any number in
this file survives, this is the one worth keeping.
MINOR | src/quebra/analyzers/checks/_permutation.py:76-78 | Tense break under the new
counterfactual: "Defaulting it to `np.random.default_rng()` would make every permutation
p-value a fresh random variable: three consecutive calls on identical input returned p =
0.3860 / 0.4040 / 0.3790." Those three numbers can no longer be regenerated by anything (the
default is now impossible), so they are rule-(2) liabilities under a "would". | "...a fresh
random variable: three consecutive calls on identical input measured p = 0.3860 / 0.4040 /
0.3790."
MINOR | src/quebra/analyzers/checks/_permutation.py:89-92 | User-facing error message with a
three-way tense break and a ragged wrap the edit left: "None would draw OS entropy, which made
permutation p-values irreproducible across runs while the run identity stayed unchanged." |
"None draws OS entropy, which makes permutation p-values irreproducible across runs while the
run identity stays unchanged." Re-wrap.
MINOR | src/quebra/plots/targets.py:50-51 | Tense break: "static and academic would otherwise
write the same {name}.pdf, so academic clobbered static while the prov record listed both
targets." | "...so academic would clobber static while the prov record listed both targets."
MINOR | src/quebra/analyzers/checks/battery.py:68 | Chronology the purge missed inside a file
this commit edited, and the banned form: "# CvM, promoted 2026-08-14." | "# CvM."
MINOR | src/quebra/cli.py:152-153 | Chronology left on the two lines under the one that was
purged: "The directory is presentation only now. `jobs/archived/` does not exist and its glob
was dead, so `--include-archived` is gone with it." | "The directory is presentation only.
There is no `jobs/archived/`, so there is no `--include-archived`."
MINOR | src/quebra/recipes.py:180 | "it is a design change, not a rename. Placeholder for
future work." replaced "and it is Increment C work" with status narration (rule 4), and it can
be misread as saying `XI_SEED`'s VALUE is a placeholder - it is the live seed. | Delete the
sentence; "it is a design change, not a rename." is the whole point.
MINOR | src/quebra/recipes.py:327-328 | Past tense left where present tense is both true and
checkable: "two job files that were byte-identical once the date and the run duration were
normalised". MEASURED today: `jobs/active/t2star_q1_070423.py` and `t2star_q1_100423.py` are
34 lines each and differ at exactly two lines, `PREFIX` and `duration_h`. Present tense turns
an obituary into an assertion a reader can check in one command. Also a ragged wrap. | "Collapses
the family: the two T2* job files are byte-identical once the date and the run duration are
normalised."
VERIFIED | src/quebra/core/job.py:150-153, :167-172 | The two best rewrites in the commit.
Both convert history into a live hazard in the present tense, keep the error string and the
mechanism, and lose nothing. No change wanted.
VERIFIED | src/quebra/analyzers/checks/c2_anderson_darling.py:82 | "a nan statistic would slip
past the positivity guard" - correct counterfactual, correct tense, hazard intact.
VERIFIED | src/quebra/analyzers/checks/result.py:45-47 | "`validate_segment` would check
`tau > np.sum(x)`" is a true counterfactual: the shipped guard is
`segment.tau <= total * (1.0 + TAU_MARGIN)` (:175) against `last_event_time`, not `np.sum`.
VERIFIED | src/quebra/plots/km_survival_plot.py:25 and mtbc_hist_plot.py:17 | "no longer
visible" -> "not visible". Chronology gone, claim unchanged, both still true (the counts are
on the materialized artifacts and absent from the images).
VERIFIED | src/quebra/analyzers/fidelity.py:188-190 | "Kept out of the render layer: building
panel data inside build_matplotlib happens at draw time, where no DAG node can supply the
window and read tables." Matches PANEL_CONTRACT.md:17-19 and the live layout; the deleted
pointer to `plots/fidelity_plot.py` named a file that no longer exists, so nothing was lost.
VERIFIED | src/quebra/panels/_within_calibration_render.py:11 and within_calibration.py:15 |
`analyzers.within_calibration_compute` is the real module and `build_within_calibration_panel_data`
is really there. Both pointers now correct.
VERIFIED | src/quebra/panels/_check_ledger_render.py:24 | `"cvm_cramer_von_mises": "CvM"` -
the key matches `CHECK_NAME` as written into `CheckResult.check` (cvm_cramer_von_mises.py:255)
and matches the `ROW_KEYS` entries.
VERIFIED | src/quebra/transforms/interpolate.py:176-179 | The cross-reference is true:
`filter._subset_norm` raises `"filter cannot mask '{key}': it is neither a scalar nor an array
({type})..."` on exactly this input class (filter.py:50-54), and the new message mirrors it.
OPTIONAL | src/quebra/panels/_within_calibration_render.py:8-9 | Adjacent pre-existing breakage
in the docstring the commit edited: "They do arithmetic - `adaptive_ylim` takes
`observed_slices` searches for cut indices - but only to decide where to put ink". Two verb
phrases collided; the line is also ~110 characters against a file that wraps at 90. | "They do
arithmetic - `adaptive_ylim` picks limits, `observed_slices` searches for cut indices - but
only to decide where to put ink."
OPTIONAL | src/quebra/panels/_within_calibration_render.py:9 | "which CLAUDE.md places on the
renderer's side of the line" - `CLAUDE.md` is now a symlink to `AGENTS.md`, and AGENTS.md is
the file that carries the rule. Source should cite the real file. | "which AGENTS.md places...".

## Reviewed (patch 4 prose) - src/** chronology purge: DONE

### jobs/**, scripts/**, conftest.py, pyproject.toml, docs/WRITING_A_JOB.md  [reviewed]

MUST | jobs/bench/report.py:299 | NEW FALSE CLAIM, and the code contradicts it twice. "the
"0.25" arm realises 0.167-0.169 throughout, so the realised value is printed, never the label."
MEASURED: report.py:322-324 builds column headers as
`f"target {c:g} (realised {realised...:.3f})"` - the label IS printed, beside the realised
value - and report.py:210 prints `c={worst['censoring_target']}`, the label ALONE, in
`worst_cell`. | "the "0.25" arm realises 0.167-0.169 throughout, so every table prints the
realised value beside the target label." (Line is also 99 chars against a file wrapping at 88.)
MUST | jobs/active/t2star_q1_070423.py:5-6 AND jobs/active/t2star_q1_100423.py:5-6 | Broken
sentence in both files: deleting "128" left the noun behind. "Everything shared lives in the
recipe; these were lines each and byte-identical once the date and the run duration were
normalised." "these were lines each" is not a sentence. | "Everything shared lives in the
recipe: the two files are byte-identical once the date and the run duration are normalised."
MEASURED: both files are 34 lines and differ at exactly two lines (`PREFIX`, `duration_h`), so
the present tense is checkable in one command and the past tense is not.
MUST | jobs/rscripts/reference_values.R:76-77 | Botched rewrite: the replacement duplicates the
phrase from the line below it. Reads "Our estimator uses the tie-corrected eq (8), not the
tie-free reduction; / the tie-corrected eq (8), and no shipped window exercises it". Also
"the only external evidence that the change was made correctly" now refers to a change with no
antecedent. | "# THE CASE THAT MATTERS. Our estimator uses the tie-corrected eq (8), not the
tie-free / # reduction, and no shipped window exercises it - all 309 are tie-free. This is the
only / # external evidence that the tie-corrected path is right." (Line 76 is 100 chars.)
SHOULD | jobs/bench/probe_unresolved.py:45 | The rewrite inverted the warning and lost the
hazard. "Writing through repo_root() must point at `jobs/bench/results/`." `repo_root()` does
not point there and this code deliberately does not use it; the original said writing through
`repo_root()` landed in a `bench/results/` that does not exist. MEASURED: there is no
top-level `bench/` in the checkout, so the mistake is still available to make. | "# Beside
this file, like every other bench output: `repo_root()` would resolve to a top-level
`bench/results/`, which does not exist."
SHOULD | jobs/bench/report.py:557-558 | Verb deleted with the chronology, leaving a noun
phrase and a dangling conjunction: "# Computed, not typed: hand-written ranges ("44-90 null
cells", "z_crit 3.3-3.5") and both were wrong." | "# Computed, not typed: hand-writing these
two ranges got both of them wrong ("44-90 null cells", "z_crit 3.3-3.5")."
SHOULD | conftest.py:50 | "so it is not a small edit. Placeholder for future work." The
docstring already opens with "FUTURE:", so the added sentence is the third marker of the same
status in five lines, and it is rule (4) narration. Second occurrence of this exact phrase in
the commit (see recipes.py:180). | Delete "Placeholder for future work."; end at "so it is not
a small edit."
SHOULD | scripts/acceptance.sh:4 | Chronology/status the purge missed on the line under the one
it fixed: "The claim this phase makes is ..." | "The claim is "a reviewer runs the documented
install steps on a machine that has never seen this repository"."
MINOR | jobs/active/*.py, jobs/composite/*.py (8 files) | The `# SPEC 0005 R5.2/R5.4:` fix left
a ragged two-line wrap in all eight ("# The logical name and category. `include` resolves
JOB_ID, so this" / "# file can move without breaking any composite; recategorising costs one
string edit."), and the same 16 lines of boilerplate are duplicated across eight job files to
explain a convention `docs/WRITING_A_JOB.md` already documents under "What a job file
declares". | Delete the comment from all eight files (16 lines), or collapse each to one line:
"# Logical name and category; `include` resolves JOB_ID, so this file can move."
MINOR | scripts/make_data_manifest.py:3, make_fixtures.py:12, make_job_manifest.py:3,
promote_run.py:3 | Four module docstrings now open with a short orphan line where the spec id
was cut, against paragraphs that wrap at 88. Correct edits, unreflowed. | Re-wrap the four
paragraphs.
MINOR | jobs/bench/report.py:350 | Same unreflowed wrap: "Criterion 4 is scored here rather
than asserted in prose: for every" then "row, ...". Also "Criterion 4" is a bench requirement
id of the kind rule (1) targets, and it is now unattributed - nothing in the file says what
criterion 4 is. | Re-wrap, and name it: "This scores what the bench's criterion 4 asks for
rather than asserting it in prose: for every row, ...".
VERIFIED | pyproject.toml:57-59 | `qretool` -> `quebra` in all three project URLs, consistent
with `CITATION.cff:15-16` and `README.md:100`. No stale name left in the file.
VERIFIED | docs/WRITING_A_JOB.md:246 | "the directory decides nothing; the declaration does" -
true: `discovery` reads `JOB_ID`/`JOB_FAMILY`/`JOB_SWEEP` statically and `cli.py` walks all of
`jobs/`, so no code branches on `active` vs `composite`.
VERIFIED | scripts/make_job_manifest.py:6 | "That argument holds at 12 jobs" - 12 files under
`jobs/` declare `JOB_ID` (9 active, 3 composite). Correct.
VERIFIED | jobs/bench/report.py:449 | "a hand-typed driver string would print 0.999" - correct
counterfactual; the code takes the mean (`directed["rejection_rate"].mean()`), so 0.999 is not
printed today.
VERIFIED | jobs/bench/report.py:350 and scripts/*.py, acceptance.sh:2,72 | All seven spec-id
deletions removed only the id. No factual claim changed and no pointer went stale.
OPTIONAL | src/quebra/core/closure.py:240, core/discovery.py:16, scripts/make_job_manifest.py:3 |
The purge removed `SPEC 000N RN.N` everywhere but kept `quebraplan.md` section numbers
("3.1 asks for", "3.5 suggests", "3.4 makes the argument") in three places. Defensible - they
are citations to a live document, not chronology - but state the rule once in AGENTS.md so the
next pass does not delete them as a fourth defect class.

## Reviewed (patch 4 prose) - jobs/**, scripts/**, conftest.py, pyproject.toml,
## docs/WRITING_A_JOB.md: DONE

### THE CHECK COUNTS - the one the commit missed

MUST | src/quebra/analyzers/checks/c3_serial_copula.py:27-28 | The commit updated
"five checks"->six and "four permutation checks"->five in AGENTS.md, the docs and battery.py,
and left the OTHER count untouched, where it is now false in the direction that hides
evidence. "the other four checks must still run. The promotion report scores four checks, not
five." With CvM in `ROW_KEYS`, five other checks run (C1, C2, C5, C6, CvM) and the report
scores five, not four. | "the other five checks must still run. The promotion report scores
five checks, not six."
MUST | src/quebra/analyzers/instrument_validation.py:619 | Same count, in the string that
WRITES the shipped report, and it is measurably false: "`promotion_report.md`, which scores
four checks - C3 and CvM have no bench cell." MEASURED: `jobs/bench/results/size_table.csv`
carries 219 CvM rows and `power_table.csv` carries 348. CvM has bench cells. This is the exact
claim `tests/test_instrument_validation.py:136-140` was written to forbid, one file over. |
"`promotion_report.md`, which scores five checks - C3 has no bench cell."
IMPORTANT | jobs/bench/results/instrument_report.md:23 and
jobs/bench/results/promotion_report.md:1,5 | The committed artifacts carry the same stale
count: "which scores four checks - C3 and CvM have no bench cell", "# Promotion report - the
five checks", "**Four checks are assessed, not five.**". Six checks exist and five are
assessed. Not a prose edit - these regenerate - but nothing in this commit regenerated them,
so the promotion lands with three published sentences contradicting it. | Re-run the report
builders in the same commit as the count change.
MINOR | tests/test_checks_c3_bridge.py:64-65 | The skip message edited by this commit now
points at the stale count: "C3 is assessable, and the promotion report's four-check scope
needs revisiting on this machine." | "...and the promotion report's scope needs revisiting on
this machine." (Drop the number; the report owns it.)

### tests/** chronology purge  [reviewed]

MUST | tests/test_checks_statistics.py:132-133 | Subject deleted with the chronology, leaving
ungrammatical text: "not borrowed from a nearby run - a borrowed value would be / passed only
because the tolerance is 4 SE wide." | "not borrowed from a nearby run - a borrowed value
would pass anyway, because the tolerance is 4 SE wide." Re-wrap (line 133 currently ends on an
orphan "This").
MUST | tests/test_identity_closure.py:84-86 | Broken sentence with a capital mid-clause and a
doubled "so": "The T2* jobs were the original measurement for this, but / Both collapse onto
`recipes.configure_t2star_job`, so all their steps / live in a `quebra.*` module, so they do
not exercise this path." | "Uses `mtbf_q1`, not a T2* job: both T2* jobs collapse onto
`recipes.configure_t2star_job`, so all their steps live in a `quebra.*` module and they do not
exercise this path."
MUST | tests/test_identity_closure.py:154-155 | Broken: "so without the / parameter row was
folded these produced ONE digest — measured, both `93f7286634eef03b`." Also the digest is
asserted nowhere (this docstring, `closure.py:249` and `spec/` only). | "so without the
parameter row folded in, these two collide on one digest."
MUST | tests/test_independence_survey.py:136-137 | Tense collapse mid-sentence: "That is what
a test building its own copy of the ladder does: it asserts the / T2* job matched IT, so
mutating the SURVEY's ladder changed nothing and the test passed." | "That is what a test
building its own copy of the ladder does: it asserts the T2* job matches ITS copy, so mutating
the SURVEY's ladder changes nothing and the test still passes."
MUST | tests/test_instrument_validation.py:138-139 | Two sentences fused into a run-on by the
rewrite: "The tier-3 detail must not read "no bench cell exists for CvM" once the bench
carries CvM rows that sentence is false, and a report repeating it would understate the
evidence". | "The tier-3 detail must not read "no bench cell exists for CvM": once the bench
carries CvM rows that sentence is false, and a report repeating it would understate the
evidence."
MUST | tests/test_instrument_validation.py:136 | Chronology MISSED, and the identical sentence
one file over WAS purged - so the commit fixed one copy and left the other. Here:
"PROMOTED 2026-08-14; this test pinned the opposite and is inverted, not deleted."
`tests/test_checks_cvm.py:137` got "CvM is promoted; this test pins the promoted behaviour...".
| Delete the sentence: the test name already says `test_cvm_is_promoted_and_...`.
MUST | tests/test_job_discovery.py:143-145 | Verbless fragment: "Replacing the
`jobs/active/*.py` glob with "every discovered / job", which silently widened `run --all` from
9 jobs to 12 — pulling in the three / composites...". | "Replacing the `jobs/active/*.py` glob
with "every discovered job" silently widens `run --all` from 9 jobs to 12 — pulling in the
three composites that ...". (9 -> 12 VERIFIED against disk.)
MUST | tests/test_r_cross_implementation.py:20-21 | Botched rewrite that duplicates the phrase
below it and produces gibberish: "`chatterjee_xi` uses the tie-corrected eq (8), not the
tie-free / reduction to the tie-corrected eq (8)." Same defect as
`jobs/rscripts/reference_values.R:76`, from the same edit. "that change" at :21 and "the change
was made correctly" at :23 then refer to a change no longer named. | "**The headline case is
`xi_tied`.** `chatterjee_xi` uses the tie-corrected eq (8), not the tie-free reduction, and no
shipped window exercises it: all 309 window-rows on the T2* ladder are tie-free. This file is
the ONLY external evidence that eq (8) is right, and `XICOR::xicor(ties = TRUE)` is the
reference." (Line 20 is 100 chars.)
MUST | tests/test_tier3_calibration.py:205-207 | Broken beyond repair by clause-splicing: "The
measurement lives in `analyzers/instrument_validation.py` and is IMPORTED here, not / defined
here rather than in three docstrings and the report's / tier-3 verdicts claimed the report
imported it - it did not, and the pipeline cannot / import `tests/`, so those verdicts rested
on no number at all." | "# The measurement lives in `analyzers/instrument_validation.py` and
is IMPORTED here, not defined / # in `tests/`. The pipeline cannot import `tests/`, so a
measurement defined there would leave / # the report's tier-3 verdicts resting on no number at
all. Importing it in this direction is / # what makes "the figure and the test measure the same
thing" true rather than asserted."
SHOULD | tests/test_data_unavailable.py:3 | Uses a phrase on the author's own banned list, on
the line the purge edited: "Before this, both situations raised the same `FileNotFoundError`".
| "Without the distinction, both situations raise the same `FileNotFoundError`: a reviewer who
lacks an embargoed record and a user who mistyped a path get identical output, and only the
first is in a position to do anything about it."
SHOULD | tests/test_paths_dataset_fallback.py:36-39 | Banned phrase left in the paragraph the
purge edited: "The earlier / version asserted that THIS checkout's dataset root lacked the
table, which held only / while the root was the repo's parent; it is now the repo itself".
Three tenses and an obituary to say one thing. | "Driven from an explicit empty root rather
than `default_dataset_root()`: asserting that THIS checkout's dataset root lacks the table
holds only when the root is the repo's parent, and here it is the repo itself. `tmp_path`
tests the fallback itself, under every configuration."
SHOULD | tests/test_bench_isolation.py:3-4 | Chronology missed on the very line that was
rewritten, and the replacement clause is a non-sequitur. "The rule changed shape in P4b and it
is worth being precise about what it now protects, because a docstring is not a guard." "P4b"
is a phase id; "it is worth being precise" is filler; and if a docstring is not a guard, that
is an argument for LESS docstring, not more precision. | Delete both lines. The next paragraph
("**What moved.**") already carries the content.
SHOULD | tests/test_checks_statistics.py:219-220 | Present-tense rewrite now contradicts the
test body directly under it: "Clamping a materially negative eq (10) variance to zero reports
it as "every gap is identical", which is a different and false diagnosis." The shipped code
does not clamp, and the test asserts the two messages differ. | "Clamping a materially
negative eq (10) variance to zero would report it as "every gap is identical", a different and
false diagnosis."
SHOULD | tests/test_checks_cvm.py:137-138 | The purge removed the date and left the
chronology-about-itself: "this test pins the promoted behaviour and is inverted, not deleted,
so the change of decision is visible in the history rather than silent." A test does not need
to narrate its own edit history. | Delete the sentence; start the docstring at "Registration
alone is not the promotion."
SHOULD | tests/test_load_dataset_contract.py:127-129 | Over-correction: a false NUMBER was
replaced by a vague quantifier when the true number is one and is stated in the next clause.
"It is NOT relied on by many call sites: every `job.load` in `jobs/` names a schema, and
`tests/test_windows_not_interpolated.py` is the only site that reaches the default." VERIFIED
true: all six `job.load` sites in `jobs/` pass a Dataset carrying a schema, and
`test_windows_not_interpolated.py:190` is the only `Dataset(...)` with no schema. | "Exactly
one site relies on it: every `job.load` in `jobs/` names a schema, and
`tests/test_windows_not_interpolated.py` is the only site that reaches the default."
SHOULD | tests/test_checks_statistics.py:20-21 and :330 | The same rewritten sentence appears
twice in one file, and both copies are awkward: "The second is where a defect merging renewal
segments across unobserved read gaps is invisible to inspection." A location is not where a
defect "is invisible". | Keep one, at :20: "The second covers `segments_from_windows`, where a
gap-merge defect - renewal segments spliced across unobserved read gaps - cannot be caught by
reading the code." Delete the :330 copy.
MINOR | tests/test_bench_isolation.py:37-39 | Ragged wrap left by the id deletion, plus a
stranded past conditional: "`jobs/` deliberately sits at the root, because ... and moving it
would have put `output/` inside site-packages." | "...and moving it would put `output/` inside
site-packages." Re-wrap.
MINOR | tests/test_bench_isolation.py:92-93 | Half-purged: "which does not match "jobs.bench"
now that the package moved - so the guard admitted a real, working bench import." Present tense
then "now that ... moved" then "admitted". | "which does not match "jobs.bench", so the guard
admits a real, working bench import."
MINOR | tests/test_path_resolution.py:33-34 | The contrast now names a file that does not
exist anywhere in the tree: "the repo root is the / directory holding pyproject.toml, not the
one holding main.py." | "the repo root is the directory holding pyproject.toml." Re-wrap.
MINOR | tests/test_path_resolution.py:39, test_independence_survey.py:180 | "now" survives in
two purged sentences ("`quebra.toml` now declares", "which is where the T2* job now gets it").
| Drop both "now"s.
MINOR | tests/test_independence_survey.py:80-82 | Half-purged: "Both T2* jobs go through
`recipes.configure_t2star_job`, so these three moved from inline step kwargs into the recipe's
signature. Reading only the job file would now find nothing". | "...so these three live in the
recipe's signature rather than in inline step kwargs. Reading only the job file finds
nothing".
MINOR | tests/test_identity_closure.py:45 | Chronology left directly under the docstring line
the purge fixed: "Until the within-calibration compute was moved out of `panels/`,
`analyzers.t2star` reached `plots.base` ...". | "With the within-calibration compute in
`panels/`, `analyzers.t2star` reaches `plots.base` ...".
MINOR | tests/test_data_manifest.py:3, test_job_manifest.py:3, test_packaged_fixtures.py:3,
test_promote_run.py:3, test_paths_dataset_fallback.py:8, test_job_discovery.py:3,
test_path_resolution.py:33 | Seven more docstrings opening on a short orphan line where the
spec id was cut. Correct edits, unreflowed. | Re-wrap.
VERIFIED | tests/test_checks_statistics.py:14-16 | "The distinction matters: reading the gate
as proof of the whole implementation hides that the shipped asymptotic path is measurably
oversized at small n (0.0634 against 0.05 at n = 20)." Best rewrite in the tests group: the
chronology is gone, the claim is present tense and true, and 0.0634 is ASSERTED - :125
parametrizes `[(20, 0.0634), (50, 0.0514)]`. Both numbers regenerate.
VERIFIED | tests/test_bench_isolation.py:93 | "planting that line in jobs/active/ left all 19
tests passing" - `pytest --collect-only` reports exactly 19 tests in this file today. The
number is still current.
VERIFIED | tests/test_job_discovery.py:165, test_load_dataset_contract.py:60,
test_within_calibration_builder.py:30, test_r_cross_implementation.py:84,169,
test_identity.py:1, test_provenance_graph.py:1, test_reuse_gate.py:1, test_path_resolution.py:1 |
Clean rewrites: increment/spec labels removed, no claim changed, nothing stranded.
VERIFIED | tests/test_paths_dataset_fallback.py:8 | "This checkout's data lives under `data/`,
so its declared root is the repo itself and the two candidates coincide HERE." True, and the
"property of one configuration" caveat is preserved.

## Reviewed (patch 4 prose) - tests/** chronology purge: DONE
## ALL MANIFEST ITEMS COMPLETE

---

SCOPE: naming and structure only - check_outcome.py "dispatcher" terminology audit
COMMIT: ab885a0
(Appended, not overwritten: the P5 findings above are committed work.)

## Manifest
(all empty - review complete)

## Reviewed
- src/quebra/analyzers/check_outcome.py
- src/quebra/plots/check_outcome_plot.py
- src/quebra/analyzers/check_selection.py (terminology only)
- src/quebra/analyzers/assumptions.py (terminology only)
- jobs/active/km_with_checks_6d2s.py (terminology only)
- spec/specvalidity08.md (terminology only)
- src/quebra/plots/theme.py (diff)
- tests/test_check_outcome.py (docstring only)

## Findings

IMPORTANT | src/quebra/analyzers/check_outcome.py:1 | Module docstring calls the module "the dispatcher"; the repo's established name for exactly this entity is output-builder (CLAUDE.md:75). No entity in this codebase is called a dispatcher; "dispatch" is used only for schema selection in job.load. | Replace with "the output-builder".
IMPORTANT | spec/specvalidity08.md:322 | "**Data side, the dispatcher.** An output-builder resolves the display-set..." uses both names for one thing in one sentence. | Drop "the dispatcher"; keep "Data side."
MINOR | spec/specvalidity08.md:5 | "two selectors, a dispatcher and a render path". | "an output-builder".
MINOR | spec/specvalidity08.md:320 | Heading "R8.4 - The dispatcher and a separate, Phase-7-readable plot". | "The output-builder and ...".
MINOR | spec/specvalidity08.md:488 | "CHECKPOINT 8.4b - the dispatcher, the plot and the new job". | "the output-builder, ...".
MINOR | jobs/active/km_with_checks_6d2s.py:163 | Step wrapper docstring "The dispatcher: every drawn record's rows, reshaped ...". The step is not the builder; it is a wrapper. | Delete the word: "Every drawn record's rows, reshaped into one grid per shown check."
MINOR | tests/test_check_outcome.py:1 | "the dispatcher reshapes a ledger into grids". | "the output-builder reshapes...".
MINOR | src/quebra/analyzers/check_outcome.py:38 | `CELL_ABSENT = "no row"` is a sixth verdict string, but the five others are `VERDICT_*` in check_ledger.py:54-58 and theme.verdict_color raises on anything not in VERDICT_COLORS. Constant name and home both diverge from the family it joined. | Name it `VERDICT_ABSENT` and define it beside the other five in check_ledger.py.
MINOR | src/quebra/analyzers/check_outcome.py:64 | `OutcomeGrid.verdict_counts()` is a method; the same concept on `InstrumentGrid` is a property named `counts` (independence_survey.py:141). | Rename to `counts`, property or method, matching the twin.
MINOR | src/quebra/analyzers/check_outcome.py:31 | Imports CHECK_LABELS/CHECK_NULL from independence_survey; check_selection.py:34 imports C3_KEY/SURVEY_KEYS from it too. independence_survey is now the de-facto vocabulary module for check row keys and labels while being named as a figure builder. | Move the vocabulary (C3_KEY, SURVEY_KEYS, CHECK_LABELS, CHECK_NULL) to check_selection.py, which already owns RowKey and ALL_KEYS.

NOT findings, verified:
- Entity kind: same as the three cited builders. DataFrame/artifacts in, typed artifact out, in analyzers/, drawn by a BasePlot subclass.
- Module layout: types+builder in one module matches independence_survey.py and instrument_validation.py. The *_compute/*_data split is a two-instance panel-artifact pattern (within_calibration, across_calibration), described in CLAUDE.md 203-206 as layout for that panel, not as a rule.
- `build_check_outcome` matches `build_<module stem>` exactly, as all four other build_* output-builders do.
- `*Data` suffix is NOT the convention: 3 of 41 top-level classes in analyzers/ carry it, and only 2 of 11 StaleArtifactGuard artifacts repo-wide (both `*PanelData`). CheckOutcome sits with CheckLedger, SignalBand, ReliabilityBand, DistinguishBand, SizeVsN, ValidationCurve, ReadDependence. OutcomeGrid mirrors InstrumentGrid.
- Plot: `plots/check_outcome_plot.py` / `CheckOutcomePlot(BasePlot)` / `_draw_grid(ax, grid, ...)` all match independence_survey_plot.py.
- Pre-existing "dispatch" at core/job.py:156, schemas/ramsey_series.py:8, tests/test_load_dataset_contract.py:10,206,252 is schema selection. Correct usage, leave alone.

---

SCOPE: CHECKPOINT 7.1 - tracked working-tree diff vs HEAD (12 files). Excludes untracked
spec/specrboundary07.md and spec/specpresentation09.md.
COMMIT: 278a2fd
(Appended, not overwritten: prior findings above are committed work.)

## Manifest
- [ ] docs/iid_checks/C3_serial_copula.md (+/-52)
- [ ] docs/iid_checks/LIMITATIONS.md (+/-9)
- [ ] jobs/bench/report.py (+/-25)
- [ ] jobs/bench/results/promotion_report.md (+/-2)
- [ ] jobs/composite/independence_survey.py (+/-21)
- [ ] spec/quebraplan.md (+/-11)
REVIEWED battery.py
REVIEWED c3_serial_copula.R
- [ ] src/quebra/analyzers/checks/c3_serial_copula.py (+/-44)
- [ ] src/quebra/analyzers/independence_survey.py (+/-17)
- [ ] src/quebra/plots/independence_survey_plot.py (+/-24)
- [ ] tests/test_checks_c3_bridge.py (+/-17)
- [ ] CROSS-CUT: tenth-site hunt, number reproduction, em dashes

## Reviewed

## Findings

### src/quebra/analyzers/checks/battery.py
No finding. Verified: `len(ROW_KEYS) == 9`, no C3 entry; `grep -ic c3` on
jobs/bench/results/{size,power}_table.csv gives 0 and 0; `jobs/bench/arms.py` has no C3
mention; `check_ledger.py:175,199,223,463` carries `include_c3` and invokes the bridge
directly. The retired "(it does not here)" is gone. No em dash, no spec id.

### src/quebra/analyzers/checks/c3_serial_copula.R
No finding. "Historic: R was absent..." deleted; the surviving sentence states the
contract (non-zero exit rather than a plausible number) in the present tense. Site 5 of
nine, corrected.

### src/quebra/analyzers/independence_survey.py
No finding. Verified by import: `len(ROW_KEYS)==9`, `len(SURVEY_KEYS)==10`, `C3_KEY` last.
New prose is count-free ("One grid per entry of `SURVEY_KEYS`"), the CvM-appears-twice
addition is correct (cvm asymptotic + cvm permutation both in ROW_KEYS), and the
`build_independence_survey` comment is now count-free. "132 of 340 cells per grid" is
per-grid and unaffected by the instrument count (34 datasets x 10 rungs = 340). Retired
"the survey turns it off to stay affordable" is gone and matches
`jobs/composite/independence_survey.py` INCLUDE_C3=True.

### src/quebra/plots/independence_survey_plot.py
IMPORTANT | src/quebra/plots/independence_survey_plot.py:203 | TENTH SITE of the same
count-drift class the diff fixed three lines above: "each of the fourteen panels" is
7 instruments x 2 clocks. `SURVEY_KEYS` has 10 entries and `build_matplotlib` lays out
`len(keys) x len(clocks)` = 20 panels. The diff de-counted the class docstring's first line
and left the stale count in the same docstring. | "each of the twenty panels", or de-count
it as the first line was.
IMPORTANT | src/quebra/plots/independence_survey_plot.py:244 | Same stale count, same class:
"fourteen copies of a thirty-four-label axis". Twenty today. | De-count or say twenty.

### jobs/bench/report.py + jobs/bench/results/promotion_report.md
Verified: regenerating `build(runtime_note)` into a scratchpad file and diffing against the
tracked `promotion_report.md` gives IDENTICAL (no write to the repo). The R-absent branch
was rendered by monkeypatching `c3.rscript_path` and also reads true. No absolute path
survives. `ruff check jobs/bench/report.py` exits 0; `ruff format --check` on all seven
changed Python files reports "already formatted". No test regenerates and compares the
report, so a no-R CI machine will not fail a gate on this.

IMPORTANT | jobs/bench/results/promotion_report.md:5 | The abs-path correction removed the
path but the TRACKED generated artifact still encodes a machine fact: "`Rscript` was on
PATH when this report was generated". Regenerating on a machine without R rewrites this
tracked file, so `python jobs/bench/report.py` is no longer reproducible from the two
committed CSVs the report's own header line claims it is built from. | State it
unconditionally in both directions ("where R is present the bridge returns an uncalibrated
p-value; where absent it degrades to `not computed`") and drop the probe, or move the probe
result into an untracked runtime note.
MINOR | jobs/bench/report.py:544-547 | R-present branch reads "the bridge executes and
returns a real p-value" inside a paragraph about THIS bench, which never calls C3. The next
clause rescues it, but the first-read subject is ambiguous. | "...the bridge executes for
callers that request C3".

### tests/test_checks_c3_bridge.py
Behaviour: UNCHANGED. Read the full file: no assertion added, removed or altered; the skip
condition at :66 is still `rscript_path() is not None`; no marker or parametrize change.
Only the module docstring, one test docstring and the skip message moved.

IMPORTANT | tests/test_checks_c3_bridge.py:6-8 | The rewrite replaced one false environment
claim with a false claim about the file: "Every test here drives that behaviour directly,
monkeypatching the probe ... so the file asserts the same thing on a machine with R and on
one without". `test_r_is_genuinely_absent_here` (:60-72) does NOT monkeypatch and is exactly
the test whose outcome depends on the machine - it skips where R is present and asserts
where it is absent. | "All but one test monkeypatches the probe; `test_r_is_genuinely_absent_here`
reads the real machine and skips where R is present."
MINOR | tests/test_checks_c3_bridge.py:60 | The function NAME still asserts the retired
claim: `test_r_is_genuinely_absent_here`. Its own docstring was corrected to "Records
whether `Rscript` is absent"; the node id a reader sees in `-v` output was not. | Rename to
`test_records_whether_rscript_is_present`.

### jobs/composite/independence_survey.py
IMPORTANT | jobs/composite/independence_survey.py:104-109 | "Over this survey's 415
populated cells that BRACKETS the run between roughly a quarter of an hour ... and roughly
six hours". The arithmetic checks (415 x 2.1 s = 14.5 min, 415 x 54.6 s = 6.3 h) but it is
not a bracket: the endpoints are the smallest and largest MEASURED n (50 and 355), not the
survey's smallest and largest CELL, and nothing here establishes that no cell exceeds
n = 355. The one datum in the repo about a larger n was deleted from
`c3_serial_copula.py` in this same diff. | Say "brackets the run IF every cell falls between
n = 50 and n = 355, which has not been checked", or drop "brackets" for "suggests".
MINOR | jobs/composite/independence_survey.py:104-106 | Quotes 2.1 / 7.1 / 54.6 s as "the
measured cost" while `c3_serial_copula.py:212-216` says of those exact figures that "a third
significant figure here would be noise quoted as measurement" and gives a second run at
1.7 / 6.2 / 57.2. Two tracked sites disagree about how precisely the same measurement may be
quoted. | Quote the rounded form here too, or name it as one of two runs.

### docs/iid_checks/C3_serial_copula.md
Statistics and p-values REPRODUCED, all seven quoted values, exactly (see cross-cut below).
"Historic note" block removed; the `four-check` skip-message reference removed; the
"it reads the interpreter off the machine it is generated on" claim matches report.py.

### docs/iid_checks/LIMITATIONS.md
IMPORTANT | docs/iid_checks/LIMITATIONS.md:44-47 | The rewritten sentence ends "so large
windows can exhaust the timeout and report `not computed` for that reason alone", but the
timeout is 900 s (`c3_serial_copula.py:212`) and the largest measured cost is now 54.6-57.2 s
at n = 355. `c3_serial_copula.py:228` says the opposite in the same tree: "900 s covers every
n the T2* ladder has produced so far with margin". Two tracked sites now contradict each
other, and the one datum that supported "can exhaust" (n = 682 unfinished at 580 s) was
deleted in this same diff. | Either restore the n = 682 observation as the open question
CLAUDE.md sec 4 prescribes, or drop "can exhaust the timeout" for "costs minutes rather than
seconds".

### spec/quebraplan.md
No finding. Section 7 item 2 amended in place with the dated search, its result, and the
"a search is a negative result with a date on it" caveat. `statsmodels` claim not
independently checked (it is a negative search result, correctly labelled as one).

## CROSS-CUT

VERIFIED (CRITICAL check, passed) | All quoted C3 numbers reproduce. Ran the stated recipe
(one `numpy.random.default_rng(0)` drawing `exponential(1.0)` at n=50, 150, 355 in that
order) through `c3._invoke_rscript(x, 5, seed, ., 1000)` under R 4.5.3 / copula 1.1.7:
  seed=1 n=50  stat 0.00579077616284928 p 0.957542457542458  -> quoted 0.0057908 / 0.9575 OK
  seed=1 n=150 stat 0.00713259001510053 p 0.903596403596404  -> quoted 0.0071326 / 0.9036 OK
  seed=1 n=355 stat 0.00763312792149939 p 0.865634365634366  -> quoted 0.0076331 / 0.8656 OK
  seed=0 n=50  p 0.947552447552448                            -> quoted 0.9476 OK
  seed=0 n=150 p 0.911588411588412                            -> quoted 0.9116 OK
Every digit quoted at every site matches. No CRITICAL finding on the numbers.

IMPORTANT | src/quebra/analyzers/calibration_summary.py:252 | ANOTHER site of the corrected
class, outside the diff. "Drawing all seven put C1 and C2 on the page twice" - the power
table carries NINE distinct check/calibration/variant rows (measured from
`jobs/bench/results/power_table.csv`), and CvM appears twice as well. The diff corrected
this sentence's exact twin in `independence_survey.py` ("C1, C2 and CvM each appear twice")
and left this one. | "Drawing all nine put C1, C2 and CvM on the page twice".
MINOR | src/quebra/analyzers/checks/c3_serial_copula.py:210-213, docs/iid_checks/C3_serial_copula.md:14-16,
docs/iid_checks/LIMITATIONS.md:45, jobs/composite/independence_survey.py:104 | The module
docstring says of 2.1 / 7.1 / 54.6 that "a third significant figure here would be noise
quoted as measurement", then three other tracked sites quote exactly those three-significant-
figure values as "the measured cost", with no spread. A third run during this review measured
61.2 s at n = 355, outside the docstring's own "55 to 57 s" band. | Quote the rounded form at
all four sites, or state two runs at all four.
MINOR | src/quebra/analyzers/checks/c3_serial_copula.py:210-212 | "TWO RUNS, AND THE SPREAD IS
WHY THESE ARE ROUNDED" but the rounding NARROWS: observed 6.2-7.1 is quoted "6 to 7" and
observed 54.6-57.2 is quoted "55 to 57". Three of four endpoints are pulled inside the
observed hull, which is the wrong direction for a rounding justified by spread. | "6 to 8 s"
and "54 to 58 s", or drop the headline and keep the two raw triples.
MINOR | src/quebra/analyzers/checks/c3_serial_copula.py:228 | "900 s covers every n the T2*
ladder has produced so far with margin" now rests on no measurement above n = 355, because
the only larger-n datum was deleted from this docstring. | Bound the claim to the measured
range, or restore n = 682 as an open question.

CHECKED, NOT FINDINGS:
- Regression greps clean outside the spec and this ledger: "is absent on this machine",
  "four-check", "n^2.8" appear in no source, doc or test.
- `check_ledger.py:477` uses past tense ("while R was absent ... With R installed the call
  executes for real"). True, not a site.
- README.md:48,120 phrase R conditionally. Not a site.
- `jobs/bench/report.py:19,21` "all seven rows" / "factor of seven" are past-tense records of
  a measurement made before CvM was promoted, and CLAUDE.md sec 4 quotes the first verbatim.
- `plots/instrument_validation_plot.py:162` "four of seven rows" is the tier table, a
  different subject.
- No em dash in any added line. No SPEC/phase identifier in any added source line.
- `ruff check` exit 0 on jobs/bench/report.py; `ruff format --check` clean on all 7 changed
  Python files. `pytest tests/test_checks_c3_bridge.py` = 5 passed, 1 skipped.
- `test_there_is_exactly_one_figure_per_surveyed_instrument` exists at
  tests/test_independence_survey.py:253; the plot comment's citation is live.
- Regenerated promotion_report.md into the scratchpad: byte-identical to the tracked file.

## ALL MANIFEST ITEMS COMPLETE

---

SCOPE: FINAL SWEEP - whole working-tree diff vs 278a2fd (25 tracked files) PLUS untracked
spec/specrboundary07.md and spec/specpresentation09.md. Checkpoints 7.0-7.3 of SPEC 0007.
COMMIT: 278a2fd
(Appended, not overwritten.)

## Manifest
- [ ] CROSS-CUT P1: incomplete fixes (repo-wide grep of retired phrasings/numbers, count drift)
- [ ] CROSS-CUT P2: internal contradictions (spec vs docs vs docstrings vs promotion_report)
- [ ] CROSS-CUT P3: uncheckable claims (numbers without recipe/machine; absolutes)
- [ ] CROSS-CUT P4: writing rules (em dash, spec id in source, debugging history)
- [ ] Makefile
- [ ] conftest.py
- [ ] docs/iid_checks/C3_serial_copula.md
- [ ] docs/iid_checks/LIMITATIONS.md
- [ ] jobs/active/km_with_checks_6d2s.py
- [ ] jobs/bench/report.py
- [ ] jobs/bench/results/promotion_report.md
- [ ] jobs/composite/independence_survey.py
- [ ] pyproject.toml
- [ ] spec/quebraplan.md
- [ ] src/quebra/analyzers/calibration_summary.py
- [ ] src/quebra/analyzers/check_ledger.py
- [ ] src/quebra/analyzers/checks/__init__.py
- [ ] src/quebra/analyzers/checks/battery.py
- [ ] src/quebra/analyzers/checks/c3_serial_copula.R
- [ ] src/quebra/analyzers/checks/c3_serial_copula.py
- [ ] src/quebra/analyzers/independence_survey.py
- [ ] src/quebra/core/_artifact_guard.py
- [ ] src/quebra/core/discovery.py
- [ ] src/quebra/plots/independence_survey_plot.py
- [ ] tests/test_check_ledger.py
- [ ] tests/test_checks_c3_bridge.py
- [ ] tests/test_job_discovery.py
- [ ] tests/test_marker_discipline.py
- [ ] spec/specrboundary07.md (untracked, 524 lines)
- [ ] spec/specpresentation09.md (untracked, 490 lines)

## Reviewed

## Findings

---

SCOPE: COMMIT GROUP B - CHECKPOINTS 7.2 + 7.3 of SPEC 0007. Paths: conftest.py, Makefile,
pyproject.toml, .gitignore, src/quebra/analyzers/check_ledger.py,
src/quebra/analyzers/checks/c3_serial_copula.py, src/quebra/core/_artifact_guard.py,
tests/test_check_ledger.py, tests/test_marker_discipline.py, docs/iid_checks/LIMITATIONS.md
COMMIT: 278a2fd (working tree, uncommitted)
(Appended, not overwritten. Prior 7.1 pass above covered the PROSE-only half of
c3_serial_copula.py and LIMITATIONS.md; the 7.2/7.3 code is new and reviewed fresh here.)

## Manifest
- [ ] src/quebra/analyzers/checks/c3_serial_copula.py (+/-119)
- [ ] src/quebra/analyzers/check_ledger.py (+/-80)
- [ ] conftest.py (+/-56)
- [ ] tests/test_marker_discipline.py (+/-174)
- [ ] tests/test_check_ledger.py (+/-303)
- [ ] src/quebra/core/_artifact_guard.py (+/-2)
- [ ] Makefile (+/-3)
- [ ] pyproject.toml (+/-2)
- [ ] .gitignore (+/-6)
- [ ] docs/iid_checks/LIMITATIONS.md (+/-32)
- [ ] EXPERIMENTS: exception-escape matrix, state-collision matrix, mutation runs, cost

## Reviewed

## Findings

## Manifest ADDENDUM (scope given to this pass includes 3 paths the manifest above omitted)
- [ ] AGENTS.md (+21)
- [ ] tests/test_checks_c3_bridge.py (+42/-)
- [ ] spec/specrboundary07.md (untracked, 555 lines; R7.2, R7.3, Mutation ledger)

## Reviewed (Group B)
- [x] src/quebra/analyzers/checks/c3_serial_copula.py
- [x] src/quebra/analyzers/check_ledger.py
- [x] src/quebra/core/_artifact_guard.py
- [x] AGENTS.md
- [x] EXPERIMENTS: exception-escape matrix (DONE, all clean), sentinel leak (DONE, leaks)

## Findings (Group B)
IMPORTANT | src/quebra/analyzers/checks/c3_serial_copula.py:109-114 | the empty-probe sentinel
  added for the PROVENANCE consumer leaks into the FUNCTIONAL one: `_invoke_rscript:204`
  does `env["R_LIBS"] = os.pathsep.join(r_library_paths())`, so the branch hands R the
  literal string `probe returned nothing (exit 0)` as a library search path. Measured:
  R_LIBS handed to the bridge == 'probe returned nothing (exit 0)'. Fix: move the sentinel
  into `check_ledger._r_provenance` (which already owns the `probe failed (...)` sentinel)
  and leave `r_library_paths()` returning `()`.
MINOR | c3_serial_copula.py:113 and :158 | both empty-probe sentinels interpolate
  `{completed.returncode}` in a branch only reachable when returncode == 0 (non-zero exits
  earlier in both functions), so the code always emits "(exit 0)". Dead interpolation that
  invites a reader to think it varies. Fix: drop the parenthetical.
IMPORTANT | c3_serial_copula.py:262-266 | "No interval is quoted, and that is deliberate."
  is contradicted two sentences later by "have ranged from about 50 s to about 75 s across
  six runs". Also violates the AGENTS.md rule this same commit adds. Fix: delete the range
  sentence or delete the "no interval is quoted" claim.
IMPORTANT | AGENTS.md:419-423 | the new wall-clock rule's cited evidence is false as shipped:
  "the fourth attempt stopped quoting an interval at all" - the fourth attempt
  (c3_serial_copula.run, same commit) quotes 50-75 s, 130.2 s and 580 s. Fix: correct the
  docstring, or restate the evidence.
MINOR | c3_serial_copula.py:280-286 | cites `output/check_ledger_q1_070423_.../` for
  statistic 1.511 / p 0.0005 / n 682. Numbers VERIFY against the pickle's bytes (1.51123,
  0.0004995, 682 present), but output/ is gitignored AND the pickle names module
  `analyzers`, not `quebra.analyzers`, so the shipping code cannot unpickle it. A permanent
  docstring citing an artifact no reader can open. Fix: state the numbers with their job
  and date, not the ephemeral directory.
MINOR | c3_serial_copula.py:287 | "It still fails fast..." is glued to the end of the
  inserted paragraph with no blank line, and the line runs ~104 chars. Fix: reflow.
MINOR | check_ledger.py:143-145 | `c3_n_null_sim` is the one of the four new fields with no
  `not asked` state: it records `inputs.c3_n_null_sim` (module default 1000) even when
  `include_c3=False` and no C3 row exists. Its three siblings got a sentinel. Documented in
  the comment, so not a silent fallback, but the twin site is untreated.
MINOR | check_ledger.py:346 | the `r_version()` except clause omits `RuntimeError` while its
  twin at :338 includes it. Not reachable today (r_version returns a string on non-zero
  exit) but the asymmetry is unexplained at the call site.
VERIFIED CLEAN | `_r_provenance` exception escape: stubs for hang (120.1 s, TimeoutExpired),
  non-zero exit, vanished path (FileNotFoundError), non-executable and directory
  (PermissionError), silent, and undecodable bytes all returned normally. Nothing escaped.
VERIFIED CLEAN | `_artifact_guard.py` message change is accurate for the new break.

## Reviewed (Group B) - continued
- [x] conftest.py
- [x] tests/test_marker_discipline.py
- [x] tests/test_check_ledger.py
- [x] tests/test_checks_c3_bridge.py
- [x] Makefile
- [x] pyproject.toml
- [x] .gitignore
- [x] docs/iid_checks/LIMITATIONS.md
- [x] spec/specrboundary07.md (R7.2, R7.3, Mutation ledger)
## Manifest (Group B) - EMPTY, all items reviewed

## Findings (Group B) - continued
IMPORTANT | src/quebra/core/_artifact_guard.py:1,3,22,27 | the runtime message at :49 was
  generalised from "predates the builder/renderer split" to "predates a field added to this
  artifact", but the module and class docstrings still call the mixin's users "panel-data
  artifacts" / "panel-data dataclasses" and still frame every break as "pre-split".
  `CheckLedger` (check_ledger.py:120) is neither. One site of the class fixed, four left.
IMPORTANT | spec/specrboundary07.md:449-452 | the mutation ledger says "Every test these two
  checkpoints add or change" and "A test with no row here is not evidence", but
  `test_with_no_rscript_the_r_fields_read_absent_and_every_row_is_still_built` has no row -
  and it is the test R7.3's own Acceptance clause names. (Also
  `test_records_whether_rscript_is_absent_here`, renamed with no row.) Measured: an M12
  collapsing `absent` to `not asked` DOES turn it red, so the test is sound; the ledger's
  totality claim is not. Fix: add the row, or drop the totality sentence.
MINOR | src/quebra/analyzers/checks/c3_serial_copula.py:157-160 | the success value is
  `"R version 4.5.3 (2026-03-11) 1.1.7"` - the copula version is an unlabelled trailing
  token. Every failure string names itself; the success one does not. Fix: `cat(..., "
  copula", ...)` in the probe.
MINOR | tests/test_marker_discipline.py:154-157 | the test asserts
  `completed.returncode == 0` on a child that runs the WHOLE `-m policy` selection, so any
  unrelated policy failure turns this one red too. The file-existence assertion already
  covers the real condition. Measured: with `.git` absent, test_stale_references failing
  made this test fail with a child-run traceback.
MINOR | conftest.py:103-110 | the fixture hands back the LIVE stash list, not a copy; and
  the guard docstring (test_marker_discipline.py:87) says it "sees every item the run
  collected", slightly stronger than the fixture's own accurate "from under this conftest's
  tree". Measured: an item collected from outside the rootdir does NOT enter the record
  (n_record=1 while 2 tests ran).
MINOR | conftest.py:113-132 | `requires_rscript` ships with no consumer and no test; its
  first user arrives in R7.4. Nothing catches it breaking until then.
MINOR | pyproject.toml:110 / Makefile:62 | "needs Rscript and the copula package" is true
  today, but R7.4(b) adds a meta-row test checking XICOR, energy and randtests versions
  against the local R, which will make the description false again.
INFO | budget | 7.2+7.3 add 9 collected tests (2 marker + 6 ledger + 1 bridge) against
  stated envelopes of +1 (R7.2) and +2..+4 (R7.3). Under the 2x halt, but the banner should
  state it.

## MUTATIONS MEASURED BY THIS REVIEW (all 11 spec rows + 1 extra), in a scratch mirror
M1 RED test_a_ledger_that_never_asked... | M2 RED same | M3 RED (OSError escape kills the
ledger) test_a_version_probe_that_raises... | M4 RED test_an_unusable_r... | M5 RED
test_the_recorded_n... | M6 RED test_c3_n_null_sim... (+ M5's test) | M7 RED
test_the_record_is_pre_deselection_not_post (+ tryfirst test) | M8 RED
test_the_record_survives_a_competing_tryfirst_deselector | M9 RED
test_every_collected_test_carries_exactly_one_tier_marker (+2) | M10 RED / M11 RED
test_a_probe_that_answers_nothing... | M12 (extra) RED test_with_no_rscript...
SPEC MUTATION LEDGER IS HONEST: every row reproduces, caught by the test it names.
ALSO VERIFIED: the `-p` plugin competitor does NOT beat the rejected modifyitems form
(recorded=True) while the shipped subdirectory-conftest competitor DOES (recorded=False),
so AGENTS.md:321-329's second instance and the test comment are both accurate.

## GATES (real tree, read-only)
ruff check 0 | ruff format --check 0 (164 files) | mypy src/quebra/core 0 | lint-imports 0
mirror `pytest -m "not slow and not heavy and not r"`: 598 passed 3 skipped (the one extra
skip vs the stated 599/2 is test_artifact_guard.py:319, which needs output/, excluded from
the mirror). `pytest --collect-only -q`: 601. No planted files left in tests/ in either tree.

## Group B verdict: DO NOT SHIP (no CRITICAL; 5 IMPORTANT, all text or latent)

---

SCOPE: README.md, pyproject.toml, docs/iid_checks/C3_serial_copula.md, spec/quebraplan.md,
src/quebra/analyzers/instrument_validation.py, tests/test_instrument_validation.py,
jobs/bench/results/instrument_report.md  (CHECKPOINTS 7.5 + 7.7 slice)
COMMIT: cf55810
## Manifest
- [ ] pyproject.toml (+0/-2)
- [ ] README.md (+30/-0)
- [ ] docs/iid_checks/C3_serial_copula.md (+30/-0)
- [ ] src/quebra/analyzers/instrument_validation.py (+?/-?)
- [ ] tests/test_instrument_validation.py (+14/-?)
- [ ] jobs/bench/results/instrument_report.md (+1/-1)
- [ ] spec/quebraplan.md (+85/-0)
## Reviewed
## Findings

---

# CHECKPOINT 7.4 review - slice: the two r-marker test files

SCOPE: tests/test_r_cross_implementation.py, tests/test_checks_c3_bridge.py
COMMIT: cf55810
FREEZE: `git diff -- . ':(exclude).claude/review-findings.md' | md5sum` = 70349dcd1030b0194026ea0ad9c29567
        at start. (Bare `git diff | md5sum` moved only because a CONCURRENT reviewer
        appended 17 lines to this ledger. Source tree unmoved.)
BASELINE in scratch mirror: `pytest -m r -q` 3 passed, 604 deselected => collect 607.

## Manifest
- [ ] tests/test_r_cross_implementation.py  (+97/-0)  DONE
- [ ] tests/test_checks_c3_bridge.py        (+40/-0)  DONE

## Reviewed

## Findings

---

SCOPE: CHECKPOINT 7.6 slice - src/quebra/analyzers/check_ledger.py,
src/quebra/analyzers/independence_survey.py, tests/test_check_ledger.py
COMMIT: cf55810   TREE HASH AT START: 70349dcd1030b0194026ea0ad9c29567 (matches, frozen)

## Manifest
- [ ] src/quebra/analyzers/check_ledger.py  (+14/-0)
- [ ] src/quebra/analyzers/independence_survey.py  (+0/-1)
- [ ] tests/test_check_ledger.py  (+79/-0)

## Reviewed

## Findings

### src/quebra/analyzers/independence_survey.py  (+0/-1)
MINOR (verified clean) | independence_survey.py:305 | deleted `"c3_excluded": True`.
`git grep -n c3_excluded -- src/ jobs/ tests/` exits 1; `grep -rn c3_excluded .` finds it
only in spec/specrboundary07.md (the spec that ordered the deletion) and this ledger. No
test asserted it. Shipped artifact meta now reads
`{'alpha': 0.05, 'n_datasets': 34, 'n_instruments': 10}`. Correct and complete. | none

IMPORTANT | independence_survey.py:148-157 `rejection_share_of_decided` | For any instrument
whose `bench_accepted` is never True, `pass` is UNREACHABLE (`_verdict` sends every
non-rejection to `underpowered`), so `fail/(pass+fail)` is identically 1.0 - a constant, not
a measurement. MEASURED on the shipped artifact
output/independence_survey_d6fe38_20260927_095456/independence_survey_grids.pkl:
C3 in_spec pass=0 fail=46 share=1.000000; C3 calendar pass=0 fail=14 share=1.000000
(46+14 = the 60 fails this checkpoint is about). Same degeneracy at C1-asymptotic in_spec
(pass=0 fail=8 share=1.000000), where `bench_acceptance_at_n` returns bench_accepted=False
at ALL SIX grid n. This number reaches a shipped caption:
plots/independence_survey_plot.py:266-268 prints
`f"{g.counts['fail']} of {g.counts['pass']+g.counts['fail']} decided cells reject"`, so the
C3 figure states "46 of 46 decided cells reject". AGENTS.md section 4: that is a pooled/forced
100%, not a measured rate. | Return NaN (or a distinct sentinel) when `pass` is structurally
unreachable for the instrument, or drop the share for grids with no reachable `pass`.

MINOR | independence_survey.py:309-326 `survey_summary` | No caller anywhere in src/, jobs/
or plots/; the only reference is tests/test_independence_survey.py:382 asserting its row
count. A public function whose sole consumer is a test. | Either wire it to the sink that
needs it or note it as debug-only in the docstring.

### src/quebra/analyzers/check_ledger.py  (+14/-0)
VERDICT DECISION (weight 2): keeping `fail` is RIGHT. Demoting an uncalibrated rejection to
`underpowered` would bucket "the statistic landed in the tail" with "we had no power and saw
nothing", which is a strictly worse claim and hides a finding - against AGENTS.md section 5
("a reader is better served by a band they are told not to trust than by a missing one").
A sixth verdict is the principled answer but is not free: plots/theme.py:304-319
`VERDICT_COLORS` has six fixed keys and `verdict_color` RAISES on an unknown verdict, and
`independence_survey.VERDICT_ORDER` is a fixed 5-tuple keyed by `counts`, `_pivot_verdicts`
and two legends. Out of scope for 7.6.

IMPORTANT | check_ledger.py:301-307 vs panels/check_ledger.py:118 | FIX EVERY SITE OF THE
CLASS (AGENTS.md section 4). The defect class is "`bench_accepted is None` handled unlike
`is False`". The diff fixes it in `_verdict`'s note and leaves the live twin in the render
layer: `hatch=MISCALIBRATED_HATCH if bench_ok is False else None`. A C3 cell (None) draws
with NO hatch, pixel-identical to a calibrated cell, and the legend
(panels/check_ledger.py:136-142) has only "bench: miscalibrated at this n" - no entry for
"no bench cell". This is LIVE, not latent: `CheckLedgerInputs.include_c3` defaults to True
(check_ledger.py:190), jobs/composite/check_ledger_q1.py:100-116 does not override it and
draws `CheckLedgerPanel`, and panels/_check_ledger_render.py:21 carries an explicit
`"c3_serial_copula": "C3 cop"` column abbreviation. The panel's own comment at
check_ledger.py:113-117 makes the argument against this diff's chosen remedy: "The hatch says
so at a glance; burying it in a notes column would not." | Hatch (differently) when
`bench_ok is None` too, and add the legend entry. Note the docstring at
_check_ledger_render.py:53-55 already promises None "is a different thing from False".

MINOR | check_ledger.py:291-299 | The 9-line comment for a 5-line branch ends with change
history: "An earlier version named only the miscalibrated case, so an UNBENCHED check's
rejection was indistinguishable from a calibrated one." AGENTS.md section 10 bans exactly
this shape in prose you write ("debugging history: they date on the first refactor"). It also
restates the pre-existing 2-line comment above it. | Drop the last sentence and the
restatement; keep the constraint ("no bench cell means no measured level").

MINOR | check_ledger.py rejection branch reach | The new note reaches the materialized
per-dataset ledger (`notes` is in LEDGER_COLUMNS:115, materialized at
jobs/composite/check_ledger_q1.py:122) and `check_outcome.OutcomeGrid.notes`, but NOT the
survey: `InstrumentGrid` carries verdicts/p_values/n_events only, no notes
(independence_survey.py:135-141). The survey figures - the artifact that actually drew the
60 red cells - are told only by the caption at independence_survey_plot.py:99-103 and
:274-275. Those captions are accurate today, so this is reach, not error. | none required;
noted so nobody reads the note as covering the survey.

MINOR (latent, not live) | kaplan_meier.py:172-186 `check_summary()` | Tallies verdicts with
no calibration distinction, so an uncalibrated `fail` would count as a plain `fail` in a
band caption. Not reachable today: km_with_checks_6d2s.py:73 sets
`RUN_SET = PERMUTATION_KEYS`, and check_outcome.py:64-67 filters on
`key[1] == CALIB_PERMUTATION`, which excludes C3 (`r_copula`). | none now; will bite if a job
puts C3 in a band's run-set.

### tests/test_check_ledger.py  (+79/-0)

GUARD IS REAL (weight 1). The derivation is not a re-type: `benched` comes from
`pd.read_csv(size_table.csv)["check"].unique()`, `surveyed` from `SURVEY_KEYS`. MUTATIONS
RUN in a scratch mirror (PYTHONPATH=mirror/src, real tree untouched):
  M-A  append a c3_serial_copula row to size_table.csv -> RED
       `assert {'c3_serial_copula'} == set()`
  M-B  add ("c7_brand_new","permutation","") to SURVEY_KEYS -> RED
       `assert {'c3_serial_copula'} == {'c3_serial_c...c7_brand_new'}`
  M-C  delete the new `elif bench_accepted is None` branch -> RED, and ONLY the new test
       (`1 failed, 19 passed` over the whole file)
  M-D  add c3 to battery.ROW_KEYS unbenched -> RED (assert 3 fires)
All four reverted; mirror diff-clean against the real tree.

MINOR | tests/test_check_ledger.py:544 | `assert UNCALIBRATED_CHECKS <= {key[0] for key in
SURVEY_KEYS}` CANNOT FAIL. If assert 1 passed then `UNCALIBRATED_CHECKS == surveyed - benched`,
which is a subset of `surveyed` by construction. A dead assertion - the shape AGENTS.md
section 7 names ("a test can cite a source and still be unable to fail"). | Delete it.

MINOR | tests/test_check_ledger.py:545-547 | Assert 3 duplicates an existing test. M-D makes
BOTH this and `tests/test_independence_survey.py:235-241
::test_c3_is_not_in_the_row_schema_so_it_cannot_be_bench_scored` go red on the same mutation
(measured, `2 failed`). Same oracle, same subject, two files - AGENTS.md section 7 says that
is one file. The new test's unique contribution is assert 1 only, and its subject
(SURVEY_KEYS vs the bench tables) is test_independence_survey.py's, not this file's stated
subject ("the ledger's verdict rule", module docstring line 1). | Drop assert 3, or move the
whole test next to its twin.

MINOR | tests/test_check_ledger.py:533-542 | The guard's oracle is coarser than the property
it claims. The ledger's real uncalibrated condition is "no match in `bench_acceptance_at_n`
output on (check, calibration, variant, clock, n_target)" - that frame is FILTERED by
null_arms, `censoring_target <= ENVELOPE_MAX_CENSORING` and `n_used >= MIN_SUPPORT_FRACTION *
n_replicates` (calibration_summary.py:430-436). A check could have rows in size_table.csv that
are all filtered out, or be benched at only some grid n, and the test would call it "benched"
while every ledger cell still gets bench_accepted=None. I MEASURED the gap is empty today:
all 9 ROW_KEYS x 2 clocks x 6 grid n = 108 cells have an acceptance row, 0 missing. | Derive
`benched` from `bench_acceptance_at_n(size_table)` crossed with ROW_KEYS x clocks x grid
instead of from the raw `check` column; it is the same 8 lines and pins the real condition.

MINOR | tests/test_check_ledger.py:519 | `import pandas as pd` inside the test body shadows
the module-level `import pandas as pd` at line 15. | Delete the local import.

MINOR | tests/test_check_ledger.py:552 | Docstring says "Oracle: AGENTS.md section 5 - a
non-rejection is not evidence". Section 5 does not contain that sentence; its nearest claim
is "Empty means NOT ASSESSED ... a different claim from 'assessed and nothing rejected'",
which is about an empty outcome, not a non-rejection. AGENTS.md section 10 bans inventing a
locator. Also, section 7's oracle list is analytic value / reference implementation /
simulation truth - a policy doc is none of the three (arguably fine here, since the assertion
is on a note string, not a statistical result). | Cite the sentence that exists, or cite
section 4 (fix every site of a class), which is what this test actually guards.

## Reviewed
- [x] src/quebra/analyzers/check_ledger.py
- [x] src/quebra/analyzers/independence_survey.py
- [x] tests/test_check_ledger.py


### pyproject.toml - reviewed, no findings
Deleting `r = []` breaks nothing measured: `python -m deptry .` -> "No dependency issues
found"; `tomllib` parse leaves extras `['bench','dev']`; `.github/workflows/ci.yml:32`
installs `.[dev]`, `scripts/acceptance.sh:63-64` installs the wheel plus pytest; no
`quebra[r]` install route survives outside prose that denies it.

### README.md
IMPORTANT | README.md:142-143 | "`randtests`, `XICOR` and `energy` are needed only to
  regenerate the committed reference fixtures ... which the test suite never runs" is false
  for the `make test-r` route the next paragraph recommends.
  `tests/test_r_cross_implementation.py:333-368` (`r`-marked) runs
  `Rscript -e 'cat(... packageVersion("XICOR"), packageVersion("energy"),
  packageVersion("randtests") ...)'` with `check=True`. Measured:
  `Rscript -e 'packageVersion("definitelyNotInstalled")'` exits 1, so a missing package is a
  CalledProcessError -> test FAILURE, not a skip. A reader who follows the README (R +
  copula only) and runs `make test-r` gets a red gate. Makefile:62-63 carries the same wrong
  sentence (out of this slice). | Say the three extra packages are required by `make test-r`
  too, or make that one test skip on a missing package.
MINOR | README.md:150-151 | The pre-existing "The copula serial-independence check
  additionally needs R with the copula package. Without R that one check reports not
  computed" now sits at the END of the new section and repeats its opening paragraph
  verbatim in substance. | Delete the older paragraph.
Verified true: C3 is the only R check and goes out of process
  (`analyzers/checks/c3_serial_copula.py:43-46,208`); the rpy2 rationale matches that
  module's docstring lines 8-11; the package imports with no Rscript
  (`tests/test_checks_c3_bridge.py:113-147`); `jobs/rscripts/reference_values.R` is invoked
  by nothing in the suite (grep: only docstring mentions); `r` tests skip on a missing
  interpreter (`conftest.py:113-132`); `make test` excludes `r`, `make test-r` selects it.

### docs/iid_checks/C3_serial_copula.md
All four delegation claims verified against code, not prose:
  1. statistic not reimplemented: `src/quebra/analyzers/checks/c3_serial_copula.R:44,50`
     calls `copula::serialIndepTestSim` / `serialIndepTest`; the Python side only marshals.
  2. subprocess + CSV, not rpy2: `c3_serial_copula.py:43-46`, `:202-236`; the fresh-interpreter
     pin is `tests/test_checks_c3_bridge.py:113-147`, which hides `Rscript` in a subprocess.
  3. `r`-marked fixture pin: `tests/test_r_cross_implementation.py:303-329` asserts BOTH
     `serial_indep_global_statistic` and `serial_indep_global_p_value`;
     `jobs/reference/r_reference_values.csv:42-46` carries seed 707, N 1000, lag.max 5, and
     `jobs/rscripts/reference_values.R:163-164` matches (`set.seed(707)`, `N = 1000`,
     `lag.max = 5`).
  4. UNCALIBRATED: `battery.ROW_KEYS` has 9 entries, none C3 (`battery.py:63-79`);
     `grep -ci c3` is 0 in both `jobs/bench/results/size_table.csv` and `power_table.csv`;
     `check_ledger._verdict:300-306` writes the "UNCALIBRATED: no bench cell ..." note on
     every rejection; `tests/test_check_ledger.py:518` exists under the cited name and reads
     `size_table.csv`. No word in the section asserts calibration, power or level except as
     a denial. No em dashes.
MINOR | docs/iid_checks/C3_serial_copula.md:8-9 | "Porting it was considered and declined"
  reproduces the reasoning of `spec/quebraplan.md` row 6.2 (line 278) but cites no locator,
  in a repo whose rule is that a claim names its source. | Cite quebraplan 6.2.

### src/quebra/analyzers/instrument_validation.py + jobs/bench/results/instrument_report.md
IMPORTANT | src/quebra/analyzers/instrument_validation.py:566-575 | Tier 4's own definition
  is "does it agree with an independent implementation OF THE SAME STATISTIC"
  (instrument_validation.py:11, instrument_report.md:16). For C3 there is no second
  implementation: the fixture was produced by `copula::serialIndepTest`
  (`jobs/rscripts/reference_values.R:164-165`) and the bridge calls
  `copula::serialIndepTest` (`checks/c3_serial_copula.R:50`). What the new evidence pins is
  BRIDGE FIDELITY (argument order, CSV round trip), which is a real property but not tier
  4's question. The row directly above it, C2, is `partial` for a strictly weaker version of
  the same flaw ("tier-1 evidence wearing a tier-4 label") - C2 at least has a second
  implementation. `pass` is the strongest token in the table and the honesty now lives only
  in the detail string. | `TIER_PARTIAL` with the identical detail, or add one clause to the
  tier-4 legend saying a delegated instrument's tier-4 cell means bridge fidelity.
MINOR | src/quebra/analyzers/instrument_validation.py:571-575 | The detail names the CSV but
  not the test. `tests/test_stale_references.py` walks `src/` for test-module citations, so
  naming `tests/test_r_cross_implementation.py` in the detail would put this claim under an
  existing guard for free. Today nothing anywhere names that test (grep: one hit, its own
  `def`). | Add the module path to the detail string.
Verified: the shipped `jobs/bench/results/instrument_report.md` is byte-identical to what
  `python jobs/bench/instrument_report.py` produces from the current source (regenerated in
  a scratch mirror and diffed; no output).

### tests/test_instrument_validation.py - the guard that moved (mutations run)
Scratch mirror at scratchpad/m1 (src+tests+jobs+conftest+pyproject, never the real tree).
Baseline: `pytest tests/test_instrument_validation.py` 15 passed.
  M-B (verdict `pass`, detail "the bridge agrees with the R reference on the statistic and
  the p-value", no filename): RED at test_instrument_validation.py:253. The guard has teeth
  in the direction its docstring claims.
  M-A (verdict back to TIER_ABSENT with the exact removed text "circular as stated: 'it IS
  the R implementation' compares R to itself"): GREEN. Acceptable in principle - an `absent`
  cell claims no evidence - but the file's title property is
  "no_tier_verdict_claims_evidence_that_does_not_exist", and the C3 half of it is now
  conditional on one verdict token.
IMPORTANT | tests/test_instrument_validation.py:251-257 | The tier-4 `pass` rests on an
  `r`-marked test that the default gate never runs and that nothing names, so the claim can
  outlive its evidence with every gate green. MEASURED: renaming
  `test_the_bridge_reproduces_the_reference_serial_indep_values` and deleting its p-value
  assertion (the half the tier row calls load-bearing: "Both are needed") produced an
  IDENTICAL failure set to the unmutated mirror under
  `-m "not slow and not heavy and not r and not real"` (diff of the two FAILED lists: empty).
  `make test` excludes `r`; `make test-r` would still pass with the p-value half gone. | Name
  the test module in the tier detail so `test_stale_references` resolves it, and assert in
  this guard that the cited module exists.
MINOR | tests/test_instrument_validation.py:253-257 | The assertion is a text match on the
  detail string. It cannot distinguish "the artifact is named and the comparison was made"
  from "the artifact is named". That is inherent to a report test; worth one line saying so
  rather than leaving the docstring's "on real evidence" to carry it.

### spec/quebraplan.md
All measured numbers in 7B verified read-only against
`output/independence_survey_d6fe38_20260927_095456/independence_survey.pkl`:
  34 datasets (meta n_datasets=34); C3 grids 2 x (34 x 10) = 680 cells; 248 non-null
  p-values (124 in_spec + 124 calendar); verdict `fail` 60 = 46 in_spec + 14 calendar. Exact.
  `C3_N_NULL_SIM = 200` and `INCLUDE_C3 = True` at `jobs/composite/independence_survey.py:117-118`,
  and the run's own provenance records `include_c3=True, c3_n_null_sim=200`. 30 minutes matches
  the run dir (09:54:56) against its artifacts (10:24). `ROW_KEYS` nine entries, none C3: exact.
  Cited defects all reproduce: `jobs/bench/runner.py:209` counts `seen` unconditionally while
  `:210-211` guards `hits` on `p_value is not None`; `:169`, `:181`, `:200` catch only
  `(ValueError, KeyError)`. Section-6 "Adopt-lite, later" is quebraplan:602, section-4 optional
  image is :552, the ~30 min ACM/ICSE warning is :602, `r-lib/actions/setup-r` is
  specrboundary07.md:168,512. No em dashes.
IMPORTANT | spec/quebraplan.md:278 | Row 6.2 still PRESCRIBES "Keep `copula::serialIndepTest`
  behind `pip install quebra[r]`" and justifies it with "`[PRACTICE]`: core/extras splits ...
  are standard, and `rpy2` itself ships this pattern" - the two claims the new item 11 (:659)
  says can never work and that this same commit acted on by deleting the extra. The doc now
  contradicts itself, and item 11 lists its collisions with sections 6 and 4 while missing the
  one row that names the deleted install route. `spec/specrboundary07.md:126` already marks it
  "**Reshaped** by R7.4"; quebraplan does not. | One clause on row 6.2 pointing at 7.11.
MINOR | spec/quebraplan.md:662,678 | The new item is numbered 11 and inserted ABOVE item 10,
  so the source list reads 9, 11, 10. A renderer that renumbers will show the container item
  as 10 and the KM-gate item as 11, swapping the identifiers. | Put it after item 10 and
  renumber.
MINOR | spec/quebraplan.md:684 | "`report.py:658-673`" - the shared-permutation-set assertion
  is at `jobs/bench/report.py:666` (and again at :677); :658 is inside Criterion 3's prose
  about Kvaloy & Lindqvist. | Cite :665-673, or :666.
MINOR | spec/quebraplan.md:698-709 | 7B opens with an Implemented/Calibrated status table.
  CLAUDE.md section 8 bans status tables in agent-facing docs because they go stale, and this
  one duplicates CLAUDE.md's own "Implemented today / Not implemented" list, so there are now
  two places to update when Nelson-Aalen or MCF lands. | Keep the C3 prose, drop the table or
  make it a sentence.

## Reviewed
- [x] pyproject.toml
- [x] README.md
- [x] docs/iid_checks/C3_serial_copula.md
- [x] src/quebra/analyzers/instrument_validation.py
- [x] jobs/bench/results/instrument_report.md
- [x] tests/test_instrument_validation.py
- [x] spec/quebraplan.md

### line-number corrections to the block above
instrument_validation.py tier row: 567-575 (not 566-575).
quebraplan item 11 begins at :662 (not :659); the status table is :700-712; the
`report.py:658-673` citation is at quebraplan:748; the 9/11/10 ordering is item 11 at :662
ahead of item 10 at :681.
IMPORTANT | tests/test_r_cross_implementation.py:11 | module docstring still reads "**The
  suite never runs R.**" while this same file now holds three tests that do. Measured:
  `pytest -m r -q` runs 3 tests, all from this file, invoking /usr/bin/Rscript. The
  changeset promoted the instrument_validation.py TierRow from TIER_ABSENT to TIER_PASS for
  exactly this reason and missed the two in-file twins (AGENTS.md sec 4, fix every site of a
  class). | Rewrite: the non-`r` tests read committed fixtures and never run R; the `r`
  tier does. (The same paragraph carries a pre-existing garbled clause at :12-13, "which
  was the may be the state of a reviewer's machine" - fix while there.)
IMPORTANT | tests/test_r_cross_implementation.py:237 | test_the_c3_reference_is_pinned_with_
  its_simulation_seed's docstring says "C3's own agreement with this reference is a separate
  exercise: it needs R at test time, which this suite refuses to require." That exercise now
  lives 67 lines below at :304. | Point the sentence at the new test.
IMPORTANT | tests/test_r_cross_implementation.py:323 | test (a) drives the PRIVATE
  `c3._invoke_rscript`, not the public `c3.run` that check_ledger.py:540 calls, so the
  argument plumbing BETWEEN them is unpinned. MEASURED (mutation 11): change run()'s call to
  `_invoke_rscript(x, max_lag, n_null_sim, timeout_s, seed)` and
  test_r_cross_implementation + test_checks_c3_bridge + test_check_ledger are all GREEN
  (53 passed, 1 skipped). The shipped row would then carry notes="lag.max=5 N=1000 seed=707"
  beside a p-value R computed at seed 1000 / N 707. Same defect class as R7.4's own
  MUTATION 1, one frame up, on the production path. | Drive the test through
  `c3.run([Segment(x=x, tau=float(x.sum()))], max_lag=lag_max, seed=seed, n_null_sim=n_sim)`.
  VERIFIED free: that call returns statistic 0.00739025142047825 and p 0.791708291708292,
  matching the fixture exactly, and notes "lag.max=5 N=1000 seed=707".
IMPORTANT | tests/test_checks_c3_bridge.py:113 | named "the package imports", covers 19 of
  the 84 modules under src/quebra. MEASURED (mutation 10): an import-time
  `if shutil.which("Rscript") is None: raise ImportError` planted in
  src/quebra/analyzers/instrument_validation.py leaves the test GREEN. That module and
  independence_survey.py are the two most R-adjacent in the tree, both modified by this same
  checkpoint, and neither is reached. No other test imports with Rscript hidden. | Add both
  modules to the `-c` string, or walk `pkgutil.walk_packages`.
MINOR | tests/test_r_cross_implementation.py:378 | "the ledger passes no `timeout_s`, so its
  branch is unreachable from a real timeout" - no `timeout_s` means the 900.0 DEFAULT.
  c3_serial_copula.run's own docstring describes a row exceeding it becoming `not computed`,
  and check_ledger.py:555 catches TimeoutExpired. Reachable in production, just not cheaply
  in a test. | "cannot be driven from a test, because the ledger exposes no way to lower the
  900 s default".
MINOR | tests/test_r_cross_implementation.py:365 | non-strict `zip` over four package names
  and `out[2:]`: a short probe answer silently asserts fewer packages. Demonstrated with a
  two-version `out`: XICOR and energy asserted, randtests and copula silently skipped. Not
  live (check=True plus packageVersion's own error), and ruff's select = ["E4","E7","E9","F"]
  has no B905 to catch it. | `strict=True`.
MINOR | tests/test_r_cross_implementation.py:334 | compares major and minor only, because
  jobs/rscripts/reference_values.R:53-59 records only those. R7.4's three-way stop condition
  asks the reader to run this FIRST to rule out "R or copula having moved under it", and a
  `copula` PATCH bump is the likeliest cause of a p-value drift (the spec's own open question
  3 says so). This cannot see it. | Record full version strings in the fixture, or say in the
  docstring that patch drift is invisible here.
MINOR | tests/test_r_cross_implementation.py:389 | `pytest.raises(subprocess.TimeoutExpired)`
  does not pin WHICH subprocess timed out. `_invoke_rscript` first calls `r_library_paths()`,
  whose own `timeout=60` raises the same type, so a stalled libPaths probe passes this test
  for the wrong reason. | assert `excinfo.value.timeout == 0.5`.
MINOR | tests/test_checks_c3_bridge.py:141 | the child inherits no HOME. Harmless today
  (measured: platformdirs is not among the 19 modules imported) but quebra.core.paths uses
  platformdirs, so the first import reaching it turns this red for a non-R reason. | add
  `"HOME": os.environ.get("HOME", str(tmp_path))`.
MINOR | tests/test_checks_c3_bridge.py:143 | does not assert WHICH quebra the child imported.
  Measured: `'' in sys.path` is True in the pytest parent, so `":".join(sys.path)` carries an
  empty PYTHONPATH element, and `-c` puts CWD at sys.path[0] regardless. Resolves correctly
  today (child reported the tree under test). | have the child print quebra.__file__ and
  compare against the parent's.
MINOR | tests/test_r_cross_implementation.py:378 | "recorded as untested in the spec's Not
  done" is a phase-artifact pointer in a permanent file (AGENTS.md sec 10 bans spec/phase
  identifiers in docstrings). Softer than "SPEC 0008 R8.x" but it dates the same way. | state
  the constraint, drop the pointer.
INFO | tests/test_r_cross_implementation.py:304 | mutation 4 (write the input CSV at "%.6g"
  instead of repr) stays GREEN: serialIndepTest is rank-based, so serialisation precision is
  invisible to both assertions. Not a defect, but the bridge's write precision is unpinned
  by anything.

## MUTATIONS MEASURED BY THIS REVIEW (11), in a scratch mirror at scratchpad/rmir
M1 swap n_sim/seed argv in _invoke_rscript -> RED on the P-VALUE only (statistic assertion
   passed first). M1-counterfactual: a statistic-only copy of test (a) -> GREEN. R7.4's
   stated reason for asserting both quantities REPRODUCES exactly.
M2 delete env["R_LIBS"] -> RED ("the 'copula' package is required for C3"), and it fails
   rather than skipping, which is the `requires_rscript` fixture's stated contract.
M3 max_lag+1 on argv -> RED on the statistic.
M4 input CSV at "%.6g" -> GREEN (rank statistic; see INFO above).
M5 hardcode timeout=900 ignoring timeout_s -> RED (DID NOT RAISE).
M6 run() swallows TimeoutExpired into a declined row -> RED (DID NOT RAISE).
M7 fixture meta,copula_minor 1->9 -> RED for the new test (b); the pre-existing
   test_the_fixture_records_which_r_produced_it stayed GREEN. R7.4(b)'s premise that the
   `>= 0` rows are vacuous is confirmed, and keeping both is right.
M8 module-level Rscript probe in c3_serial_copula -> RED.
M9 module-level ImportError on missing R in check_ledger -> RED.
M10 same, in instrument_validation.py -> GREEN. Coverage gap, see IMPORTANT above.
M11 run() passes seed/n_null_sim swapped to _invoke_rscript -> GREEN. See IMPORTANT above.

## ACCEPTANCE CHECKS (all measured, scratch mirror)
`pytest -m r -q` with Rscript present: 3 passed, 604 deselected, exit 0.
With Rscript hidden behind a symlink farm (python present, no Rscript): 3 skipped, 604
  deselected, each naming the fixture's skip reason. The non-`r` import test still passes
  there (test_checks_c3_bridge.py: 8 passed).
With a MOCK Rscript on PATH (exit 0, writes nothing): all 3 FAILED. None passes with a mock.
Markers: the three `r` tests carry {'r','statistical'} - exactly one tier plus the cost
  marker. test_the_package_imports_with_no_rscript_on_path carries {'unit'} and NOT `r`, so
  it runs in CI. Correct.
Prose: no em or en dashes, no banned words in the added lines.
Budget: +4 collected from these two files against R7.4's stated +3..+5.
Collect restored to 607 in both the mirror and the real tree. Mirror diffed clean against
  the repo (src, tests, jobs/reference). No files planted in the repo.
FREEZE at end: `git diff -- . ':(exclude).claude/review-findings.md' | md5sum` =
  70349dcd1030b0194026ea0ad9c29567. Unmoved.

## Slice verdict: DO NOT SHIP (no CRITICAL; 4 IMPORTANT, 7 MINOR)

### docs/iid_checks/C3_serial_copula.md - two more
MINOR | docs/iid_checks/C3_serial_copula.md:5-6 | "the reference implementation of the
  Genest-Remillard empirical-copula test" is an unsourced superlative; the doc carries no
  year, DOI or section for Genest and Remillard anywhere (line 54 names them with no
  locator either, pre-existing). | "the only maintained implementation" with the citation,
  or drop "reference".
MINOR | docs/iid_checks/C3_serial_copula.md:3-33 | The new preamble restates two things the
  doc already said: claim 2 duplicates the "A bridge, not bindings" paragraph at :66-70
  near-verbatim, and claim 4 duplicates :46-50 ("no size or power evidence, because the
  bench never ran it and `battery.ROW_KEYS` has no C3 row"). The doc's own line 44 states
  the principle being broken: "`analyzers/checks/c3_serial_copula.py` is the one place the
  cost is quoted, so that a second copy cannot drift from it." | Make the lower sections
  point up at the preamble rather than repeat it.

## Slice verdict: DO NOT SHIP (no CRITICAL; 4 IMPORTANT, all text or drift-risk)
Freeze re-verified at the end: `git diff -- . ':(exclude).claude/review-findings.md' |
md5sum` = 70349dcd1030b0194026ea0ad9c29567, unchanged.

================================================================================
SCOPE: PHASE 7 COMMIT GROUP 1 (SPEC 0009 checkpoints 9.1/9.2/9.3) - paths:
  src/quebra/plots/targets.py, src/quebra/plots/theme.py,
  tests/test_style_baseline.py, spec/specpresentation09.md,
  tests/test_figure_determinism.py (untracked new), tests/test_poster_font.py (untracked new)
COMMIT: f134647
FREEZE at start: git diff -- . ':(exclude).claude/review-findings.md' | md5sum =
  8f12dc96c904ac918503953532f65632  (matches stated baseline)
## Manifest
- [ ] spec/specpresentation09.md  (+29/-?)
- [ ] src/quebra/plots/targets.py  (+21/-?)
- [ ] src/quebra/plots/theme.py  (+73/-?)
- [ ] tests/test_style_baseline.py  (+23/-?)
- [ ] tests/test_figure_determinism.py  (new, 113 lines)
- [ ] tests/test_poster_font.py  (new, 65 lines)
## Reviewed
## Findings

### MEASUREMENTS (all run by this review, mirrors under scratchpad/{mir,mirrepo})
Ratchet counts, repo state as it ships (glob src/quebra/{panels,plots}/*.py minus theme.py):
  bare-key `fontsize=`            -> 17     (claim verified)
  numeric-literal, both spellings -> 15     (claim verified)
  dict spelling with a literal    -> 0      ("adds none" verified)
  hex 6-digit                     -> 0
The TWO sites dropped by the narrowing are `fontsize="x-small"` at plots/allan_plot.py:123
  and :132. They are NOT `fontsize=theme.X`. The only two theme-reference sites are the DICT
  ones (check_outcome_plot.py:138, independence_survey_plot.py:144) and the bare-key pattern
  never matched them.
Planted-probe directions (mirrored repo root, real test file):
  `fontsize=9` -> CAUGHT (16>15).  `**{"fontsize": 9}` -> CAUGHT.
  `fontsize=theme.CAPTION["fontsize"]` -> not counted (correct).
  `fontsize="x-small"` -> NOT caught.  `**{'fontsize': 9}` single quotes -> NOT caught.
  `fontsize=_LOCAL_SIZE` -> not caught (pre-existing class).
Metadata mutations (mirror src, real test files):
  M1 drop metadata on static PDF   -> cross-process[static] RED, direct-property RED
  M2 drop metadata on academic PDF -> cross-process[academic] GREEN, direct-property RED
     repeated 6x under mutation: 5 RED / 1 GREEN. CreationDate is second-resolution, so the
     cross-process PDF assertion only fires when the two renders straddle a second boundary.
  M3 drop metadata on poster PNG   -> cross-process[poster] GREEN, direct-property RED
  M4 unpin Creator/Producer        -> cross-process GREEN, direct-property RED (b"Matplotlib")
Font-guard mutations:
  P1 check_fonts body disabled          -> test_poster_font RED
  P3 error message stops naming the face-> test_poster_font RED
  P2 DELETE `check_fonts(target)` from style_context -> FULL fast suite
     609 passed, 2 skipped, 3 deselected, exit 0. IDENTICAL to baseline. Nothing guards
     the wiring.
No-Roboto machine simulated (HOME redirected, fresh MPLCONFIGDIR, 132 families, Roboto absent):
  `pytest tests/test_figure_determinism.py tests/test_poster_font.py` -> 2 FAILED, 5 passed.
  Both failures are MissingFontError raised out of render_poster. They FAIL, they do not skip.

## Reviewed
- [x] spec/specpresentation09.md
- [x] src/quebra/plots/targets.py
- [x] src/quebra/plots/theme.py
- [x] tests/test_style_baseline.py
- [x] tests/test_figure_determinism.py
- [x] tests/test_poster_font.py

## Findings (Phase 7 commit group 1)
CRITICAL | tests/test_figure_determinism.py:75-113 | on a machine without Roboto the poster
  parametrisation and test_no_written_file_names_the_renderer_or_the_time FAIL, they do not
  skip (measured: 2 failed / 5 passed with the user font dir hidden); .github/workflows/ci.yml
  provisions no fonts and R9.2 itself says CI has no Roboto, so `make check` is red on CI, on
  every collaborator machine, and through scripts/acceptance.sh:103 | skip the poster cases
  when findfont(Roboto, fallback_to_default=False) raises, or land R9.2's vendor/declare/drop
  decision first.
CRITICAL | tests/test_style_baseline.py:13-14 and spec/specpresentation09.md:353-355 | the
  stated reason for 17->15 is false as measured: the two dropped sites are
  `fontsize="x-small"` at plots/allan_plot.py:123 and :132, NOT `fontsize=theme.X` (no such
  site exists in the counted tree; the two theme references are the DICT sites, which the
  bare-key pattern never matched). So the narrowing removed two genuine hardcoded sizes rather
  than two false positives, and blinded the ratchet to the whole named-size class (planted
  `fontsize="x-small"` measured NOT caught) without adding it to the :8-9 exclusion list |
  either count them (`fontsize=\s*[\d"']`) or state the real reason (named sizes scale with
  rcParams font.size) and list the exclusion.
IMPORTANT | src/quebra/plots/theme.py:413 | nothing guards the wiring: deleting
  `check_fonts(target)` from style_context leaves the full fast suite at 609 passed / 2
  skipped / 3 deselected, exit 0, identical to baseline | one test that enters
  style_context("poster") with a missing face patched in and asserts MissingFontError.
IMPORTANT | src/quebra/plots/theme.py:194-196 | the guard silently returns whenever
  font.family is not a list, so reverting the poster to its previous bare-string spelling
  disables it with no signal, and test_a_target_with_no_font_stack_is_not_checked pins that
  silence as correct | key on an explicit {target: required face} map, not on the value type.
IMPORTANT | tests/test_poster_font.py:32-37 | _clear_font_cache is a no-op on matplotlib
  3.10.8: font_manager has no module-level `_findfont_cached` (measured False; the lru_cache
  is on FontManager, i.e. font_manager.fontManager._findfont_cached), and getattr(...,None)
  swallows it while :16-17 claims the cache is cleared | clear fontManager._findfont_cached
  and let AttributeError surface, or drop the helper and the claim.
IMPORTANT | tests/test_poster_font.py | R9.2's acceptance "a test asserts the configured
  poster family resolves to itself" is not met; only structure is asserted | see CRITICAL 1:
  the positive test is what cannot ship until the font decision lands.
IMPORTANT | tests/test_figure_determinism.py:76-83 | the cross-process test is a coin flip for
  the PDF metadata (M2 measured 5 RED / 1 GREEN over six runs: CreationDate is second
  resolution) and fully inert for the PNG and for the Creator/Producer pin (M3, M4 GREEN).
  Every metadata mutation is caught only by the direct-property test. The docstring admits the
  PNG blindness, not the PDF's | state it, or keep the parametrised test only for what it
  uniquely guards.
IMPORTANT | src/quebra/plots/targets.py:17 | pinning Creator/Producer to "quebra" deletes the
  matplotlib version from the artifact, and provenance.py records no renderer version (grep:
  none), so a hash change can no longer be attributed to an upgrade - which is the exact
  payoff the test file's docstring claims | add the matplotlib version to build_prov_record.
MINOR | spec/specpresentation09.md:367-369 | the untouched "Why" still reads "the widened
  pattern gives 19 ... the true count", contradicted by the amendment above it, and cites
  independence_survey_plot.py:143 where the site is :144.
MINOR | tests/test_style_baseline.py:9 | "does not catch ... fontdict=" is now wrong
  (`fontdict={"fontsize": 9}` IS caught); single-quoted `**{'fontsize': 9}` is not caught
  (measured) and is not listed.
MINOR | tests/test_poster_font.py:61-65 | assertion-free test; stays green when the guard body
  is disabled (measured under P1).
MINOR | src/quebra/plots/theme.py:138 | the deleted block took the poster's design rationale
  with it ("everything scales up together ... same colours as the other targets"); the live
  block carries no equivalent.
MINOR | tests/test_poster_font.py:9-10 | "on the machine this was written on the font IS
  present" is drafting history in a permanent docstring (AGENTS.md s10).
MINOR | tests/test_figure_determinism.py | 9 subprocess interpreter starts, ~6s of a 47s fast
  suite, no `slow` cost marker.
MINOR | budget | the two files collect 7 (measured), against R9.1 +1..+2 and R9.2 +1..+2.
  R9.1 landed +4 = exactly 2x its ceiling, R9.2 +3. No halt, but the banner must state it.
MINOR | tests/test_poster_font.py:44 | the `is` identity pin forbids a future defensive copy.
  Harmless today: measured, matplotlib's validator copies the list into rcParams and savefig
  does not mutate the shared metadata dicts.

## Verified clean
Deleted 27 lines were genuinely the dead duplicate (its values differ from the live block:
  27.4/24.0 vs 28.0/30.0); no executable line went with it. POSTER_FONT_STACK:100 precedes
  RCPARAMS:105. No new module-level matplotlib import reaches a step (fidelity.py imports
  theme inside a function; provenance.py:118 deliberately does not import it). ruff check,
  ruff format --check, mypy (core), lint-imports all clean. Collect back to 614, nothing
  planted in the repo. No em dashes or banned words in added lines. No new .md doc.
`quebra run --all` isolates the failure (cli.py:188 catches per job); a single `quebra run`
  aborts with the actionable message; both poster jobs declare targets=["poster"] only, so no
  orphaned partial output today, but runner.py:589-591 is unguarded for a mixed-target sink.

## Slice verdict: DO NOT SHIP (2 CRITICAL, 6 IMPORTANT, 8 MINOR)

---

# Phase 7 commit group 2 (CHECKPOINTS 9.4-9.6), merged from four parallel reviewers

SCOPE: src/quebra/core/_artifact_guard.py src/quebra/core/runner.py src/quebra/panels/within_calibration.py src/quebra/panels/_across_calibration_compute.py tests/test_artifact_guard.py tests/test_within_calibration_builder.py tests/test_across_calibration_builder.py spec/specpresentation09.md .gitignore
COMMIT: 4b433e4 (uncommitted tree, frozen for the review)
Full ledgers: ~/.claude/projects/-home-sera-Desktop-polimi-thesis-code-912days-qre-tool/phase67/review_9_4_6/{A_9.4,B_9.5,C_9.6,D_cross}.md
Per-ledger verdicts: A DO NOT SHIP (0/6/6), B SHIP (0/1/12), C SHIP (0/0/7), D DO NOT SHIP (0/4/11).

## Findings (deduplicated)

IMPORTANT | src/quebra/core/_artifact_guard.py:17-18 | "moved AND gained a field still reaches StaleArtifactGuard" holds only for guard subclasses. Re-measured by main agent: 24 cached KaplanMeierComparison pickles (km_poster_6d2s) now load missing assumption_id/checks_asked/checks_unanswered/check_verdicts and read NOT ASSESSED silently; 40 cached TLF dict pickles now load holding a TLFResult without fit_failed (tlf_plot.py:88 raises mid-render). All failed at load before 9.4. (A1, D1) | refuse incomplete dataclasses anywhere in the loaded graph, or narrow the claim
IMPORTANT | src/quebra/core/_artifact_guard.py:53-57 | `except (ImportError, AttributeError): pass` hides an ImportError raised inside the moved module behind "No module named 'analyzers'", context lost. (A2, D minor) | catch only ModuleNotFoundError naming the moved module, re-raise the rest `from exc`
IMPORTANT | src/quebra/core/_artifact_guard.py:30-32 | comment says a single-module rename is a new entry; lookup matches the head component only, so such a key never matches. (A3) | longest-prefix match, or fix the comment
IMPORTANT | spec/specpresentation09.md:231-242, tests/test_artifact_guard.py:395-453 | "repairs 16" counts loads (12 complete, 4 incomplete of the 29-pickle set); "fixes the composite transport" is false because runner.py:243-254 reuses only same-commit clean-tree runs, which are written under quebra.* paths. The composite test builds an unreachable state and proves wiring only. (A4-A6, D2-D3) | restate the count and the transport claim; say "wiring only" in the test docstring
IMPORTANT | spec/specpresentation09.md:288-289 | says the gaps could be rebuilt from the window/read tables; a gap while out of spec leaves no trace there, spans exist only in WindowsResult.diagnostics. (B1) | "take them from the WindowsResult"
IMPORTANT | tests/test_across_calibration_builder.py:96-121 | R9.4b acceptance is "on the artifact"; the test calls private _binned_interval_stats, and a builder that strips NaN bins passes the whole suite (D M3). (D4) | assert through build_across_calibration_panel_data

## Mutations that leave the WHOLE suite green (test gaps)

- B M5: draw only the first stretch of each cumulative curve (drops every post-gap point).
- B M3 / D M8: label every stretch (duplicate legend entries).
- B M8: observed_slices honours only the first gap (fixture has one gap).
- C M6: np.digitize(right=True) drops the FIRST event of every record (twin of the fixed last-edge drop).
- C M7: delete the `len(edges) < 2` branch.
- D M3: builder strips NaN bins (see IMPORTANT above).
- B M9 (latent): comparing in seconds against hours makes the gap test vacuous; shipped test uses hours on both sides.

## MINOR

MINOR | _artifact_guard.py:48,57 | a renamed class raises "No module named 'panels'", reads as a packaging break. (A)
MINOR | _artifact_guard.py:52 | getattr skips find_class's dotted-qualname handling and audit hook; use super().find_class(moved, name). (A)
MINOR | tests/test_artifact_guard.py:265-269 | real-artifact test still uses pickle.load and its docstring contradicts the alias. (A, D)
MINOR | tests/test_artifact_guard.py | nothing pins "moved and stale" for a guard subclass. (A)
MINOR | tests/test_artifact_guard.py:399-402 | composite test docstring narrates history and calls the faked state "the real path". (A)
MINOR | spec:64,143,247,475 | still say 27 pickles; the amendment says 29 (the 2 extra are independence_survey). spec:145 cites _artifact_guard.py:43-51, now :80-102. spec:240 cites pickle.load, now load_artifact. (A, D)
MINOR | spec:325 | "2894 events" is the interval count; the log has 2895 events and the panel prints 2895. (C, D)
MINOR | spec:328-331 | a first or last bin needs only one empty neighbour to render as a lone point; a lone bin's IQR fill draws a faint vertical line, not nothing. mtbf_q1 has neither. (C)
MINOR | tests/test_within_calibration_builder.py:209 | cites AGENTS.md section 5 for the rule; it is docs/FIGURE_STANDARD.md:138. (B)
MINOR | tests/test_within_calibration_builder.py:211-212 | "if the reference fails the test is broken" is false: slice-helper mutations turn the reference red. (B)
MINOR | tests/test_within_calibration_builder.py:210-217, within_calibration.py:905-907,962-964 | history in docstring/comments ("already did this", "still drew"). (B)
MINOR | tests/test_within_calibration_builder.py:240 | x.min()/x.max() are NaN on NaN xdata, which would pass silently. (B)
MINOR | tests/test_within_calibration_builder.py:234-235 | for _draw_primary, "drew some line" is satisfied by the three threshold axhlines. (B)
MINOR | within_calibration.py:908-911,965-968 | slicing silently truncates an arr longer than t where the old call raised; latent (_to_full_length). (B)
MINOR | recipes.py:196,:380 (pre-existing, outside diff) | diagnostics.get("gap_spans_s") turns a missing key into no gaps; spec's "production is not affected" should mention it. (B)
MINOR | tests/test_across_calibration_builder.py:110 | comment "bin 1 (days 14-28)"; edges start at 1.0, so days 15-29. (C)
MINOR | tests/test_across_calibration_builder.py:1-5 | module docstring still claims array-equality with the pre-split values. (C, D)
MINOR | _across_calibration_compute.py:31-33 | docstring names only the exact-multiple case; the clamp also catches the last edge rounding below t1 (up to 9.1e-13 d). (C)
MINOR | budget | R9.3 and R9.4a each landed collect +3 against +1..+2: under the 2x stop, but a banner must say so. Group 613 -> 621. (B, D)
MINOR | .gitignore:63 | fixes the HEAD commit's probe, belongs to no requirement here; path matches test_marker_discipline.py:134 and the probe is removed in finally. (C, D)

## Verified clean (by at least one reviewer)
runner.py:363 is the only pipeline read site of a pickled artifact; no aliasing of core_utils/plots2-style names; every scan-clock artist split or justified; clocks consistent (t_h = t_rel_s/3600); both production callers pass gaps; NaN safe for every binned_interval_stats consumer (provenance JSON holds no artifact values, no plotly target, no nan-poisoned limits); clamp safe over adversarial and 300k random edge constructions; mtbf_q1 numbers reproduced (65 bins, one empty at day 77, 64 identical, none isolated); no live reuse risk (0 AcrossCalibrationPanelData pickles in output/, composites include only t2star jobs); no em dashes, no spec ids in source; every new test names an oracle and carries one tier marker.

## Group verdict: DO NOT SHIP (0 CRITICAL, 6 IMPORTANT after dedup, 20 MINOR, 7 whole-suite-green mutations)

## Fix pass applied (group 2), on Sera's quiz answers: refuse incomplete / keep runner wiring / gap list required / mark lone bins
All 6 IMPORTANT and the 7 whole-suite-green mutations addressed; MINORs applied except as noted in the checkpoint banner.
Mutations re-run after the fix (each RED, file restored byte-identical): 9.4 S1-S5, C1-C4, M3; 9.5 M1, M3, M5, M8, M9, M10a, M10b, Q3a-Q3d; 9.6 D-M3, C-M6, C-M7, M1, M2, L1-L4.
Not fixed, open: the same optional-gap-list pattern in analyzers/check_ledger.py:176,:231 and analyzers/checks/_multiprocess.py:160,:202 (bench calls it without gaps at jobs/bench/arms.py:314).
Re-review: not yet run on the fixed tree.

---

# Phase 7 commit group 2, RE-REVIEW after the fix pass (four reviewers; D2 stopped by Sera after leads 1-4)

Full ledgers: ~/.claude/projects/-home-sera-Desktop-polimi-thesis-code-912days-qre-tool/phase67/review_9_4_6_r2/{A2_9.4,B2_9.5,C2_9.6,D2_cross}.md
Verdicts: A2 DO NOT SHIP (0/1/9), B2 SHIP (0/0/11), C2 SHIP (0/0/5), D2 partial (0/1/7; leads 5 test sweep and 6 cross-file NOT reviewed).
Closure: every round-one IMPORTANT is FIXED; A6 (runner wiring unreachable) RECORDED; round-one MINORs FIXED or RECORDED, each verified by at least one reviewer.

## Findings (deduplicated)

IMPORTANT | tests/test_artifact_guard.py:355-395 | deleting the walk into dataclass FIELDS (_artifact_guard.py:149) leaves the whole suite green; no test nests a dataclass in a dataclass, and real artifacts do (IndependenceSurveyData.grids holds InstrumentGrid). (A2 MA1) | add an unguarded dataclass holding a stale one in a field
IMPORTANT | spec/specpresentation09.md:247-250 | the 6->8 / 27->29 explanation is false: the previous latest independence_survey run is a21fd1_20260830_110425 and is EMPTY; the latest-run rule cut at 278a2fd reproduces 27/6/21, so the one new run explains both moves. (D2, A2) | restate

## Mutations that leave the WHOLE suite green (round two)
- 9.4: MA1 fields walk removed; MA2 walk lists only; MA13 dict keys skipped; MA4 dotted-prefix clause dropped; MA8 prefix matched without "."; MA5 exc.name None read as a rename. The moved-package-itself-missing case is never exercised.
- 9.5: X1 y misaligned after the first stretch (test checks x only); X2 later thresholds draw only their first stretch (coverage is over the union of lines); X3 no labels at all (only duplicates checked); X4 t2star adapter forwards `gap_spans_s or []`; X5 length guard for the first threshold only.
- 9.6: N1-N4 lone marker at the wrong height / size 0 / labelled / under the fill (test reads x only); N5-N7 median->mean, p90->p80, q1->p20: NO test pins a statistic's value, so a wrong number reaches the artifact; N8 single-bin centre moved.

## MINOR
MINOR | tests/test_artifact_guard.py:1-8 | module docstring still says the real-artifact test uses pickle.load. (A2)
MINOR | _artifact_guard.py:100-119 | a non-slots field(default=X, init=False) is never in __dict__, so a FRESH artifact would be refused; none exists in src today. (A2)
MINOR | _artifact_guard.py:83 | nested qualname "Outer.Missing" fails loudly but without the "renamed" wording. (A2)
MINOR | spec:231-233 | 29 counts top-level pickles only; recursive gives 33 (compare job's subjobs_output adds 4). (A2)
MINOR | spec:144-146 | R9.0.5 says StaleArtifactGuard raises whenever a dataclass gains a field; only for its subclasses. (A2)
MINOR | within_calibration_compute.py:232,259,329,362 | four private helpers keep gap_spans_h=None with `or []` (same class as R9.4c; latent, reliability_band passes the list); SignalBand.gap_spans_h defaults []. (B2, D2)
MINOR | within_calibration.py:902,963 | the length guard sits after the `len(arr) == 0` skip. (B2)
MINOR | tests/test_within_calibration_builder.py:222 | "four ways" reads as exhaustive. (B2)
MINOR | spec:356-357 | "0 to 16 h, gaps=1" describes the pre-fix fixture; shipped fixture with the list omitted: one line 0-20 h, gaps=2, cumulative time 20.0 vs 14.0 h, 3 us occupancy 0.2915 vs 0.4165. (B2)
MINOR | spec R9.3/R9.4a/R9.4b Budget lines | measured +4 / +4 / +3 (files 3 / 2 / 4) against +1..+2; R9.3 and R9.4a are exactly at 2x. (A2, B2, C2, D2)
MINOR | across_calibration.py:210-211 | "a bin with data never vanishes": only the median is marked; a lone bin's p90 still draws nothing. (C2)
MINOR | spec:334 + lone-marker test docstring | "a record with no empty bin draws as before" is false for a one-bin record, which gains a marker. (D2)
MINOR | spec:355 | "marks a gap only where one falls inside an in-spec window": a window never spans a gap; "where an in-spec read borders it". (D2)
MINOR | quebraplan.md:780 | cite :225-231 (the getattr at :225). (D2)
MINOR | quebraplan.md:782-785 | "folded into tau ... nothing says so" overstates: on the in-spec clock an interior censored window raises; only the calendar clock folds it; the out-of-spec-flanked merge is silent on both. (D2)

## Verified clean in round two
No current artifact is refused (all 8 plain-loading pickles complete; 131 fresh instances of 47 dataclasses round-tripped complete); guard and walk raise the same type and prefix, nothing catches one and not the other; every caller of the three strict functions updated; windows.run has one return and always sets gap_spans_s; gap_spans_s=[] true for every gap-free fixture; float-exact coverage check; lone mask is render-side like observed_slices; legend, axis limits and PDF bytes unchanged for mtbf_q1 (0 lone bins); spec 29/8/24/20/12/4/5 reproduced; no requirement over 2x.

## Group verdict: DO NOT SHIP (0 CRITICAL, 2 IMPORTANT, ~20 MINOR, 20 whole-suite-green mutations); D2 leads 5-6 unreviewed

## Round-two fix pass applied, on Sera's answers: helpers fixed / p90 marked too / init=False exempt / D2 leads 5-6 only
Both IMPORTANTs fixed (nested-field, dict-key, set, tuple shapes in the stale test; 27->29 restated from a measurement cut at 278a2fd: 27/6/21, independence_survey a21fd1 empty).
All 20 whole-suite-green mutations now RED, plus 6 for new guards (26 total, real tree, each restored byte-identical): MA1 MA2 MA13 MA4 MA8 MA5 IF1 NQ1 | X1 X2 X3 X4a X4b X5 LG1 H1 SB1 | N1 N2 N3 N4 P1 N5 N6 N7 N8.
Scope additions: SignalBand.gap_spans_h made required (fifth site of the R9.4c class; tests/test_reuse_gate.py passes []); init=False exemption; nested-qualname walk; length guard moved before the empty skip; p90 lone marker.
Collect unchanged at 625. Gates: ruff 0, format 0, mypy 0, arch 0, deptry 0, make test 0 (620 passed, 2 skipped), make check-ci 0.
Harness note: scratchpad/mut.py was overwritten by reviewer C2's harness mid-session; mutation results above come from scratchpad/main_session_mutate.py, which targets the real tree.

## D2 round two (leads 5-6) verdict DO NOT SHIP (0/1/5), fixes applied
IMPORTANT fixed: new test_the_read_gaps_are_excluded_from_observed_time_and_timelines pins cumulative time at 4 us = observed 14 h (hand sum 5+4+5) and the two unobserved timeline holes; O21 (reliability band handed []), O22 (distinguish band handed []), O8 (_observed_dt_h stops zeroing) now RED.
MINORs fixed: runner completeness check now tested (composite wiring test runs a stale seed; O6 runner-bypass RED, M3 RED); docs/PANEL_CONTRACT.md names gap_spans_h; StaleArtifactGuard exempts init=False and the walk skips absent fields (GI1, GI2 RED); "nine" causes; "one point too long".
Collect 626 (+13 vs HEAD 613: R9.3 +4, R9.4a +4, R9.4b +3, R9.4c +2). Gates: ruff 0, format 0, mypy 0, arch 0, deptry 0, make test 0 (621 passed, 2 skipped), make check-ci 0.

---

# Phase 7 commit group 3 (CHECKPOINTS 9.7-9.10), three parallel reviewers
Ledgers: ~/.claude/projects/-home-sera-Desktop-polimi-thesis-code-912days-qre-tool/phase67/review_9_7_10/{E1_9.7,E2_9.8,E3_9.9}.md
Verdicts: E1 DO NOT SHIP (0/4/9), E2 DO NOT SHIP (0/5/14), E3 DO NOT SHIP (0/1/4).

## IMPORTANT
IMPORTANT | docs/FIGURE_STANDARD.md:8-13, spec R9.5 amendment | panels/check_ledger.py reports "no usable windows" in the panel and is omitted: five renderers, not four; the spec's :193 cite was that panel file. (E1; re-measured by main agent)
IMPORTANT | docs/FIGURE_STANDARD.md:12-13 | "C3 excluded" is false: C3 is drawn and flagged uncalibrated (independence_survey_plot.py:99-104, :274-276). (E1; re-measured)
IMPORTANT | tests/test_km_survival_plot.py:28-38 | fixture leaves n_windows_carved and n_unobserved_birth_dropped at 0, so n = (carved or n_windows) and n = n_windows + dropped both pass. (E1 K4, K5)
IMPORTANT | spec R9.5 Not done | omits transforms/filter.py's per-stage count prints; "the one that printed a dropped-data count the panel did not carry" is false for kaplan_meier.py:365-371. (E1)
IMPORTANT | analyzers/instrument_validation.py:633 | CvM "within 0.007 of nominal from n=35 up" is false: 0.0670 at n=35 shape 1.50 (0.017 off). (E2; re-measured)
IMPORTANT | analyzers/instrument_validation.py:569 | C2 "other rows ... run to 0.176": C2's n=20 shape 0.75 rows span 0.0-0.4 (0.037-0.077 without n_used<=5); 0.176 is C1's max. (E2; re-measured)
IMPORTANT | analyzers/instrument_validation.py:540 | tier-2 row hardcodes "4 of 7 ... 3 differ"; derivable from `published`. (E2)
IMPORTANT | analyzers/instrument_validation.py:546-551, 566-569, 631-638, 827, 830 | remaining bench and tie-study numbers hardcoded from tables the report never loads (correct today, except the two above). (E2)
IMPORTANT | core/runner.py:142-158 | _check_render_targets ignores includes: a composite's sub-job with a bad target fails partway, after the composite's run dir exists; at HEAD include(figures=False) completed. (E3)

## Whole-suite-green mutations
E1 W5 (only the [reliability_band] prefix is checked), K4, K5; E2-M1 (tau := gap sum in the builder: only the currency test catches), E2-M3 (the new residual<0 raise untested); E3 M4 (ghost target in the targets.py docstring), M9 (apply_common_style body unchanged, verified by reading only).

## MINOR (selected; full lists in the ledgers)
E1: mtbc_hist_plot records the same deviation as km_survival_plot; tlf_plot drops non-finite values silently; the "no curve" entry shows a line sample; survival legend runs past the figure edge (saved by bbox_inches="tight"); KM legend on the poster not re-rendered.
E2: residual guard uses x.sum() and a hard < 0 against result.py's cumsum()[-1] / TAU_MARGIN convention; {tau_h:g} goes scientific past 6 digits; "time censored" wording when 0 gaps are censored; spec's equivalent-mutant reason is wrong (c1.statistic never reads n_censored_dropped); an old pickle passes load_artifact and KeyErrors mid-render; module docstring says four checks have bench cells, there are five.
E3: spec "the one layer exception" should be "the one core-to-plots exception"; R9.7 carries two different 15s; apply_common_style now sits in theme, inside fidelity's compute closure (closure hash change); plotly has no real use left.

## Group verdict: DO NOT SHIP (0 CRITICAL, 9 IMPORTANT, ~27 MINOR, 8 whole-suite-green mutations)

## Group 3 fix pass applied, on Sera's decisions (composite check recurses; plotly recorded; report cites, not copies; hold the commit)
All 9 IMPORTANTs addressed. E1: FIGURE_STANDARD names five renderers, C3 drawn-uncalibrated, mtbc/tlf listed; km fixture separates n from carved/dropped; spec counts and Not done corrected. E2: every report number computed (published rows, R reference values, tie study, own MC) or replaced by a pointer (tests by name, bench by size_table cell key); bench claims limited to what the bench audit (phase67/bench_audit) supports; a citation test resolves every pointer and checks the bench claim in the committed table; the report regenerated. E3: the target check recurses through includes; documented targets == registered set.
MINORs applied: tau convention (last_event_time, TAU_MARGIN), non-scientific tau format, "observed to ... ending on a failure" wording, no line sample on "no curve" entries, spec wording (core-to-plots exception, two 15s, equivalent-mutant reason), module docstring (five benched checks).
Mutations (real tree, each restored byte-identical), all RED: K2 K4 K5, W5, NC1, T1 T4 T5 T6, E2-M1 E2-M3, CW1 CW2, RD1-RD6.
Not addressed (recorded): plotly (Sera: later); the bench's own defects (bench spec, after Sera's information); gold-standard values await Sera's web check.
Gates: ruff 0, format 0, mypy 0, arch 0, deptry 0, make test 0 (626 passed, 2 skipped), make check-ci 0, collect 631.

## Small review of the group-3 fix layer (F1 report: DO NOT SHIP 0/4/11; F2 rest: SHIP 0/1/6), fixes applied
F1: C5 tier 4 cited a test that never calls C5, whose lag autocorrelation differs from R's by construction (measured 0.062395 vs 0.062654): regraded partial -> absent, guarded in the audit test. C3's "iid smoke test" now cites the bridge test that runs c3.run. "max |difference|" broke the markdown table: now "max abs difference". The citation test parses the full cell key, checks the permutation half against the promotion report's envelope verdict and the envelope's definition, and unit-tests _span. Bench wording "point estimates above nominal". Tie paragraph lines split. Caption test docstring states the property.
F2: FIGURE_STANDARD names allan_plot's silent drop; nested-include errors carry an "included by" breadcrumb; the all-censored case asserts a legend exists; spec citations fixed (lines at 278a2fd, kaplan_meier :366-372), the numeric-only ADOPTED bullet marked superseded.
Mutations (real tree, restored byte-identical), all RED on a specific test: F1-M1, F1-M2, F1-M4, C5-1, BC1. Equivalent on today's data, not caught: the xi/dcor case counts typed in.
Gates: ruff 0, format 0, mypy 0, arch 0, deptry 0, make test 0 (626 passed, 2 skipped), make check-ci 0, collect 631.

---

SCOPE: bench grid extension to n=1000 (uncommitted): jobs/bench/grid.py, jobs/bench/report.py, jobs/bench/runner.py, src/quebra/analyzers/calibration_summary.py, jobs/bench/results/{size_table.csv,power_table.csv,promotion_report.md,runtime.txt}, docs/iid_checks/{BENCH,LIMITATIONS,C1_lewis_robinson,C2_anderson_darling,C5_rank_autocorr,C6_exchangeability,iid_checks_basics}.md
COMMIT: 2bd30dc
## Manifest
## Reviewed
- [x] stale-355 sweep (docs/, src/, jobs/bench/*.py, tests/)
- [x] docs/iid_checks/iid_checks_basics.md  (16)
- [x] docs/iid_checks/LIMITATIONS.md  (8)
- [x] docs/iid_checks/C6_exchangeability.md  (10)
- [x] docs/iid_checks/C5_rank_autocorr.md  (7)
- [x] docs/iid_checks/C2_anderson_darling.md  (2)
- [x] docs/iid_checks/C1_lewis_robinson.md  (17)
- [x] docs/iid_checks/BENCH.md  (69)
- [x] jobs/bench/results/promotion_report.md, runtime.txt
- [x] jobs/bench/results/size_table.csv, power_table.csv (pandas)
- [x] src/quebra/analyzers/calibration_summary.py  (+3/-3)
- [x] jobs/bench/runner.py  (+2/-2)
- [x] jobs/bench/report.py  (+4/-4)
- [x] jobs/bench/grid.py  (+6/-6)
## Findings
MINOR | jobs/bench/report.py:17 | edited line keeps "a flat tolerance rejected all seven rows" (first-draft history at the old grid) beside the new 86-144; the table has nine rows, and at this grid a flat 0.01 band rejects all nine | state the property: "a flat tolerance rejects every row, because the max of 86-144 deviations is ~3 SE by chance"
(grid.py seed claim VERIFIED: 1122 size rows at n<=355 identical to HEAD in key order and every value; n is the outer loop of size_cells. REAL_DATA_NOTES 935 / 28 VERIFIED from output/check_ledger_6d2s_q*/check_ledger_q*_*.pkl n_events. report.py:423 0.195-1.000 VERIFIED from 8 D-arm cells. calibration_summary.py:370-372 VERIFIED via bench_acceptance_at_n: C5-unstud n=20 in_spec z=-2.971 vs 2.955; C2-asym in_spec accepted at 50..1000. 86-144 matches the envelope n_cells.)
IMPORTANT | jobs/bench/report.py:681 (-> promotion_report.md "Two findings" item 1) and jobs/bench/arms.py:20 | "never positive" is now false: the regenerated report prints "every cell between -0.089 and +0.001 - never positive"; four Arm C cells have positive mean_induced_lag1 (max +0.00134 at calendar n=500 rho=0.3; also n=700 rho=0.05, n=1000 rho=0.1, size n=1000 rho=0). At HEAD max was -0.00074 | replace with a computed bound ("never above +0.002", or "within 0.09 of zero and centred below it") in both sites
(runtime.txt VERIFIED: 216 size x 2000 + 288 power x 1000 = 720000 replicates, 504 cells. Wall-clock line is a generated artifact, not prose.)
IMPORTANT | docs/iid_checks/BENCH.md:108-111 (old item 4 deleted) | half of the deleted limitation is still literally true: 100 and 355 are still adjacent, so 228-354 events are still judged at 355 (and 851-999 at 1000, 601-699 at 700); on the 28 6D2S ledgers 372 of 560 record/threshold/clock rows are judged at a LARGER n than they have, 231 of them below 20 (judged at 20), 68 in 228-354 | restore the item minus the "above 355" clause: nearest-n rounds up between grid points and below 20; keep the proposal
MINOR | docs/iid_checks/BENCH.md:66-68 | "On other data, check the event counts against N_GRID first: a record outside 20-1000 ..." implies 6D2S is inside the grid; 231 of its 560 rows have < 20 events | drop "On other data" or say the 6D2S records include counts below 20
MINOR | docs/iid_checks/BENCH.md:24, 131 | reflowed lines run to 105 and 131 chars (factors paragraph; item 7 citation) | rewrap
(VERIFIED: 216/288 cells; 504/720,000; 935; item 5 36%/34% = 0.357/0.341; item 7 162 groups, ~8 expected, three flagged permutation rows exactly as named; item 9 18 cells at c=0.25, nine inside at n=20/35, three size + six power, 1-3 replicates; groups 2/8/14/16 and +/-0.011/0.013/0.014/0.014 unchanged and correct. Numbering 1-14 contiguous; continuation indent 3 sp for 1-9, 4 sp for 10-14, consistent. No file in docs/, spec/, src/, tests/, jobs/ cites a BENCH.md item by number.)
MINOR | docs/iid_checks/C1_lewis_robinson.md:66-67 | unchanged "a quantised Arm A cell with 819 usable replicates of 2000 reaches 0.973" was the extreme of the poorly supported cells at the old grid; now poorly supported quantised Arm A cells reach 0.9983 (n=500, 587 usable) and 1.0 (n=700, 330 usable) | name the n (355) or quote the new extreme
(VERIFIED C1: 0.9625 Arm D shape 1.5 n=1000; perm 0.0405-0.0645 on supported c=0.25 cells; 128 cells, z -2.15 vs 3.55, power 0.603, range 0.195-1.000; 0.103 at 355 and 0.2045 at 1000 for A/1.5/q/c=0.03; size growth with n holds for every c=0.03 multi-segment cell (A, D; quantised or not), flat at c=0. Unchanged "within 0.009 from 35 on" still holds at 500-1000 (max dev 0.0090). C2: 2.67 vs 3.55, 0.621 verified.)
IMPORTANT | docs/iid_checks/C5_rank_autocorr.md:64-65 | not updated: "mean absolute difference of 0.016 across 192 shared cells"; regenerated check_agreement gives C5 studentized vs C6 0.018 over 288 cells (unstudentized 0.017), and C6_exchangeability.md:37 already says 0.017-0.018 / 288, so the two pages now disagree | 0.018 across 288 shared cells (or 0.017-0.018 for both variants)
IMPORTANT | docs/iid_checks/C6_exchangeability.md:47-48 | not updated and now false: "Its worst null cell in the bench came from Arm D at b = 1"; the regenerated envelope table puts C6's worst cell at A_iid_weibull/in_spec n=700 q=False c=0.0 (0.0365, z=-2.77). At HEAD it was D_trp_power_law n=100. The "cross-check that the TRP generator is sound" conclusion hangs on it | drop the paragraph or restate from the current worst cell
(VERIFIED: C5 0.106/0.236/0.357/0.626/0.804, z -3.18 vs 3.58, 0.357; C6 0.087/0.208/0.332/0.599/0.788, below studentized C5 at every n, z -2.77 vs 3.58, 0.332, 0.017-0.018 / 288 / max 0.17; LIMITATIONS table all 18 values; "below n=100 only rho=0.5 exceeds 0.5" holds; unchanged section 5 censoring ranges 0.147-0.174 / 0.171-0.195 still hold.)
(iid_checks_basics.md CLEAN: 504/720,000; 128 cells, 0.60/0.62; z 31.7 and 5.4; 144 cells z 2.56; 4.41 over 86; 0.33-0.36; 21-24% / 33-36% / 60-63% / 79-80% over C5 stud, C5 unstud, C6 all recomputed from power_table.csv.)
IMPORTANT | docs/iid_checks/iid_checks_basics.md:40-41 | not updated: "C5 ... 0.014-0.016 ... over 192 shared power cells, while C1/C2 differ from it by 0.18-0.19 over the 156 they share"; regenerated check_agreement: C5 0.017-0.018 over 288, C1/C2 0.256-0.267 over 234 | 0.017-0.018 over 288; 0.26-0.27 over 234
MINOR | src/quebra/analyzers/calibration_summary.py:470 | nearest_bracketing_n docstring "The bench measured six event counts"; N_GRID has nine | "nine event counts", or drop the count
MINOR | jobs/bench/grid.py:8-9, jobs/bench/report.py:44 | "levels are set from what the real data shows, measured before the bench was written" / "measured before the bench was built" now cover event counts that were taken from the 6D2S ledgers after it (BENCH.md:31-33 already separates the two) | scope the "measured before" clause to censoring and quantisation
MINOR | jobs/bench/grid.py:16-18 | "Every power cell is re-seeded, because power cells follow all the size cells" is history of this change in a module docstring (AGENTS.md section 10) | state the constraint: "extend N_GRID only at the end; doing so re-seeds every power cell"
(Sweep CLEAN otherwise: no other "20-355", "largest grid n", "judged at 355", grid-length or 144/192/336 count in docs/, src/, jobs/bench/*.py, tests/. Remaining 355 mentions are C3/permutation-gather cost, record examples (check_ledger.py:4,85; checks/__init__.py:6; C5 doc:48), validation_curve's own n_values, xi_ties' own grid, and instrument_validation tests at 35/355, all fine. No em dashes, spec ids or wall-clock claims in added lines; runner.py's ~16x was removed. Gates run: ruff check 0, ruff format --check 0 on the four .py files; pytest on the four table-reading test files + test_bench_isolation: 84 passed, exit 0.)

## Verdict: DO NOT SHIP (0 CRITICAL, 5 IMPORTANT, 7 MINOR). All five IMPORTANTs are one-line text fixes.

---

# C6 (SD2) patch, CP1 review

SCOPE: CP1 of C6 SD2 patch, uncommitted diff vs 3a9ab92: src/quebra/analyzers/checks/c6_exchangeability.py, src/quebra/analyzers/checks/_rank_serial.py, jobs/rscripts/reference_values.R, jobs/reference/r_reference_values.csv, tests/test_r_cross_implementation.py, tests/test_checks_statistics.py, tests/test_tier3_calibration.py, src/quebra/analyzers/instrument_validation.py, tests/test_instrument_validation.py, src/quebra/plots/instrument_validation_plot.py, jobs/bench/results/instrument_report.md
COMMIT: 3a9ab92
## Manifest
## Reviewed
- [x] src/quebra/analyzers/checks/c6_exchangeability.py
- [x] src/quebra/analyzers/checks/_rank_serial.py
- [x] jobs/rscripts/reference_values.R
- [x] jobs/reference/r_reference_values.csv
- [x] tests/test_r_cross_implementation.py
- [x] tests/test_checks_statistics.py
- [x] tests/test_tier3_calibration.py
- [x] src/quebra/analyzers/instrument_validation.py
- [x] tests/test_instrument_validation.py
- [x] src/quebra/plots/instrument_validation_plot.py
- [x] jobs/bench/results/instrument_report.md
## Findings (C6 SD2 CP1)
IMPORTANT | src/quebra/analyzers/checks/c6_exchangeability.py:59-63 | Assumption 3 says a zero-spread lag makes `studentize` raise and the ledger write `not computed`, 'there is no fallback'; VERIFIED false: on records whose segments are each internally constant ([1,1,1,1]+[2,2,2]; [1x5]+[2x5]+[3x4]; [5x6]+[1x7]) every lag has zero permutation spread, the augmented std is ~1e-17 from rounding, nothing raises, every z_h is +/-1 on every row, and C6 AND C5-studentized report p = 1.000 (3 seeds each); same on all-size-2 layouts ([2,2,2]). A degenerate cell reads as a computed non-rejection, which is the 'silence is not a pass' case, and these level-differing tied records are the docstring's own signature case. Behaviour pre-existing (old C6 also gave p = 1 there), the false claim is new | minimal for this slice: state in assumption 3 that rounding can leave a zero-spread lag unraised and what then happens; the real fix (Sera's call, touches C5): studentize raises when scale <= rtol * max(1, |centre|), e.g. rtol 1e-12, plus a test on [1,1,1,1]+[2,2,2] asserting the raise
MINOR | src/quebra/analyzers/checks/c6_exchangeability.py:23-24 | 'A sum of uncentred squares reads that offset as dependence' implies inflated size; the permutation calibration keeps any statistic exact, so the uncentred sum never over-rejected a pure level shift (sd2-update.md s.8 says so), it only biased power by sign | say the uncentred sum is not symmetric about the offset, so dependence of one sign inflates it and the other deflates it, while size stays exact
MINOR | src/quebra/analyzers/checks/c6_exchangeability.py:26-30 | unchanged lines now contradict the new text: 'The null is exchangeability of the duration sequence' (on a gapped record it is within-segment exchangeability, which assumption 2 now states) and 'the one check in the set whose validity does not rest on an approximation' (C5 reads the same z_h on the same permutations and is equally exact, as the rewritten 'Against C5' paragraph says) | 'exchangeability within each segment'; 'one of the permutation-calibrated checks' or drop the uniqueness claim
MINOR | src/quebra/analyzers/checks/c6_exchangeability.py:14 | one docstring line is 108 chars, the rest wrap at 93 (E501 not selected, so ruff stays green) | rewrap
MINOR | src/quebra/analyzers/checks/_rank_serial.py:128-129 | 'handled by studentization, which C5 and C6 share' - C5 also ships VARIANT_RAW, which does not studentize, so the sentence is false for one of the three rows that read these values | 'which C6 and C5's studentized variant share'
MINOR | src/quebra/analyzers/checks/_rank_serial.py:51-52 | 'against the one alternative they studied' asserts Ljung and Box studied only one alternative; the plan's quote ('for this particular alternative') supports only 'the alternative studied in that section' (inferred, PDF not opened) | 'against the alternative studied in that section; not shown in general'
(clean) jobs/rscripts/reference_values.R | VERIFIED: acf and Box.test consume no RNG and every generator block re-seeds, r_reference_inputs.csv is byte-identical to 3a9ab92 (git diff empty), and the LB comment's n(n+2)/(n-h) weight is correct | -
(clean) jobs/reference/r_reference_values.csv | VERIFIED: only 12 added rows; all 10 acf_rank_lag values and both ljung_box_rank_lag5 values re-derived from r_reference_inputs.csv in a scratch R run (R_LIBS=/home/sera/R/library Rscript --vanilla) to <= 4.7e-16, and the acf values match a hand Box-Jenkins sum to 1.4e-17. Neither ljung_box_rank_lag5 nor any bartels_* row is read by any test or src module (grep), as the plan intends ("context") | -
(clean) tests/test_r_cross_implementation.py | T1 reads the quebra functions (lag_layout, global_ranks, autocorrelations) against independent R acf rows at atol 1e-12, asserts the layout has five lags, and r_value raises on a missing row; the rewritten lag-1 docstring is accurate (that test calls only scipy/numpy). Module pytestmark = statistical covers T1 | -
(clean) tests/test_checks_statistics.py | T3-T6 VERIFIED to go red under scratch mutants (sd2_work/review/mutplug.py, loaded with -p): pairing across gaps -> T4 x2 + T5 red; z without centring -> T6 x2 red; moments from null only (observed row excluded) -> T6 x2 red; per-lag n/(n-h) rescale of r_h -> T1 x2 + T3 x4 + T4 x2 + T5 red; uncentred unscaled C6 sum -> T6 x2 red. Hand values in T5 re-derived by hand (18/35, 1/7, 46/105, 23/105, -5/30, -4/30 all correct); closed form in T4 re-derived (a_s^2 - W_s/(n_s(n_s-1)) per pair) | -
(clean) tests/test_tier3_calibration.py | T7 VERIFIED red under two scratch battery mutants: (a) battery resolve_perm given one block of the total size (plan's M12 shape) -> T7 x2 red, the only reds in the four touched files; (b) a cross-segment shuffle relabelled with the per-segment sizes so no consistency check fires -> T7 x2 red on the RATE (c5 0.0950, c6 0.1425 at nominal 0.01, band +/-0.0059). New seeds 90_100/90_101 do not collide with SEEDS (90000-90012). Bonferroni count 6 as planned
MINOR | tests/test_tier3_calibration.py:163-166 (note for CP3, not a defect of this file) | the plan's M12 as worded ("hands the shuffler one block of the total size", battery.py:135) turns T7 red through the ValueError at _permutation.py:213 ("blocked as (30,) but the segments are (15, 15)") raised by a sibling row, not through the exactness rate T7 asserts, so it proves the wiring check rather than T7's oracle | in CP3 also run a relabelled variant (one-block indices wrapped as PermutationSet(sizes=per-segment)); it goes red on the rate, as measured here
MINOR | src/quebra/analyzers/instrument_validation.py:627,638 | the C5/C6 tier-4 verdict is `max(abs_difference) <= AGREEMENT_TOL`, and Python's max skips a NaN that is not the first element; VERIFIED: a NaN acf_rank_lag5 for durations_trend in the values frame gives 10 rows, max 4.2e-16, verdict pass, while that row's abs_difference is nan (a NaN reference is not None, so the new raise does not fire) | `rank_acf_verdict = TIER_PASS if all(r.abs_difference <= AGREEMENT_TOL for r in rank_acf) else TIER_FAIL` (NaN compares False), or raise on a non-finite reference beside the None raise
MINOR | src/quebra/analyzers/instrument_validation.py:629-636 | the detail says "match R's acf" unconditionally while the verdict beside it is computed and can be `fail`; on a failing row the artifact would print fail next to "match" (inferred from reading) | phrase it as a comparison ("compared with R's acf ... max abs difference X against tolerance Y"), or branch the verb on the verdict
MINOR | tests/test_instrument_validation.py:339-355 | the computed C5/C6 tier-4 verdict has no test of its fail branch; VERIFIED in a scratch mirror: replacing it with a constant `TIER_PASS` leaves test_instrument_validation, test_instrument_report_currency and test_r_cross_implementation green (49 passed), and grep finds no other test touching TIER_FAIL or the rank acf rows. Mitigated: test_the_exact_reference_rows_actually_agree still goes red on any disagreeing row | one unit test that feeds build_instrument_validation (or the verdict line, factored out) a cross list with one rank acf row above AGREEMENT_TOL and asserts `fail` for both C5 and C6
MINOR | tests/test_instrument_validation.py:349-355 | the guard requires the literal "Not pinned:" but not what follows, so "Not pinned: nothing" passes, although the plan requires the clause to name the studentization, the max/sum and the pairing across gaps | also assert "studentization" and "pairing across gaps" in row.detail
(verified) tests/test_instrument_validation.py:272-288 | the missing-row test goes red when the raise is replaced by `continue` (scratch mirror, 1 failed); the guard docstring's C3 sentence now matches the code (C3 tier 4 is TIER_PARTIAL) | -
(clean) src/quebra/plots/instrument_validation_plot.py:338 | caption matches the plan's wording; the figure sizes from len(rows) and labels from statistic/case, so the 10 rows need no other change (read). Pre-existing, not this slice: every exact row is drawn in verdict_color("pass") whatever its difference | -
(clean) jobs/bench/results/instrument_report.md | VERIFIED current: tests/test_instrument_report_currency.py passes on the frozen tree (1 passed); the diff is only the C5/C6 tier rows and the 10 rank acf rows, and no other line mentions C5 or C6 | -
MINOR | src/quebra/analyzers/checks/c6_exchangeability.py:7 | 'the same values C5 takes the max of' - C5 takes the max of |z_h| (c5_rank_autocorr.py:109) | 'the values whose largest absolute value C5 reports'
## Commands run (frozen tree, PYTHONPATH=$PWD/src, quebra.__file__ = this worktree's src/quebra/__init__.py)
- ruff check . exit 0 | ruff format --check . exit 0 (196 files) | make types exit 0 | make arch exit 0 (1 contract kept)
- pytest -p no:cacheprovider on the four touched test files: 133 passed, exit 0
- full suite pytest -p no:cacheprovider -x: 1130 passed, 10 skipped, exit 0 (241 s)
- tests/test_instrument_report_currency.py: 1 passed
- collect delta read from the diff: +14 (T1 2, T3 4, T4 2, T5 1, T6 2, T7 2, missing-row 1), no test removed; plan band +13..+14
- prose scan of added lines: no em/en dash, no "gate", no spec/phase id, no Tsay
## Verdict: DO NOT SHIP (0 CRITICAL, 1 IMPORTANT, 11 MINOR, one of them a CP3 note). The IMPORTANT is a text fix for this slice (C6 assumption 3); the behaviour behind it is pre-existing and shared with C5, and fixing it is Sera's call.

# C6 (SD2) patch, CP1.1 + CP2 review

SCOPE: CP1.1 (sd2_work/review_cp11_cp2/cp11.diff, 5 files) + CP2 (sd2_work/review_cp11_cp2/cp2.diff, 7 files), uncommitted on sd2-c6
COMMIT: 3a9ab92
## Manifest
## Reviewed
- [x] src/quebra/analyzers/checks/_rank_serial.py (CP1.1)
- [x] src/quebra/analyzers/checks/c6_exchangeability.py (CP1.1)
- [x] src/quebra/analyzers/instrument_validation.py (CP1.1)
- [x] tests/test_checks_statistics.py (CP1.1)
- [x] tests/test_instrument_validation.py (CP1.1)
- [x] plan_snapshot.md R8 + CP3 notes (CP1.1)
- [x] .claude/qre_checks_reference.tex (CP2)
- [x] docs/iid_checks/BENCH.md (CP2)
- [x] docs/iid_checks/C5_rank_autocorr.md (CP2)
- [x] docs/iid_checks/C6_exchangeability.md (CP2)
- [x] docs/iid_checks/LIMITATIONS.md (CP2)
- [x] docs/iid_checks/iid_checks_basics.md (CP2)
- [x] jobs/bench/report.py (CP2)
- [x] missed C6 sites grep (CP2)
- [x] commands: ruff, format, types, arch, full suite
## Findings (C6 SD2 CP1.1 + CP2)
(verified) src/quebra/analyzers/checks/_rank_serial.py:58-62 | ZERO_SPREAD_TOL comment VERIFIED (probe.py): on 10,857 lags of 3,000 random tied gapped records no nonzero ptp fell below 3/(N(N^2-1)) (min ratio 4.0), 126 were exactly 0.0; 3/(N(N^2-1)) is 3.0000e-9 at N=1000, 1.0047e-12 at N=14,400, and first drops under 1e-12 at N=14,423, so "above 1e-12 up to N = 14,400" holds. The bound also leans on the denominator being identical in every row (fixed multiset of global ranks), which the comment leaves implicit | -
(verified) src/quebra/analyzers/checks/_rank_serial.py:170 | CP1 IMPORTANT resolved: [1x4]+[2x3], [1x5]+[2x5]+[3x4], [5x6]+[1x7] now raise at B=199 and B=999 (probe.py) | -
MINOR | src/quebra/analyzers/checks/_rank_serial.py:171-175 | the message names one cause, "as it does when every segment is constant inside", but a lag also has no spread when every segment pairing at lag 1 has exactly two events, whatever their values; VERIFIED [1,2]+[3,4]+[5,6] raises with this message at B=199 and 999. A ledger with min_events_per_segment=2 on a heavily gapped record gets a reason pointing at a cause the record does not have (the same sentence is in c6_exchangeability.py:61-62) | "as it does when every segment is constant inside, or when every segment that pairs at that lag has only two events"
(note) src/quebra/analyzers/checks/_rank_serial.py:51-52 | the locator reads "Sect. 4.3, p. 302"; the plan's Resolved opens say "Sect. 4.3, pp. 301-302", while Sera's CP1.1 instruction says "p. 302". The quoted phrase is inside the Resolved-opens quote, so the narrower page is consistent with both; flagged only because the brief asks for an exact match. CP2's docs must use the same page (checked under CP2) | Sera's call; no change needed if p. 302 is where the phrase sits
(verified) src/quebra/analyzers/checks/c6_exchangeability.py:59-64 | "the ledger writes the whole cell as `not computed` with the reason" and "the bench drops the replicate from every row" both RUN, not read: check_ledger._rows_for on [1x4]+[2x3] (segments_from_windows monkeypatched) returns ONE row, check_id "(all)", verdict "not computed", notes "battery failed: a lag has no spread ..." (the inner retry at check_ledger.py:497-511 raises again); bench run_cell with every fifth replicate replaced by that record gives n_used 16 of 20 and n_failed 4 on every row of both an include_tau_checks=True and a False cell (probe_bench.py) | -
MINOR | src/quebra/analyzers/checks/c6_exchangeability.py:62-64 | "C1, C2 and CvM included" lists the battery rows but not C3: _rows_for returns the "(all)" row at check_ledger.py:509-512 before the C3 call at :539, so C3 is never run either, though it shares nothing with the permutation set; VERIFIED include_c3=True and False both give the single "(all)" row. New fact for open decision (a): the loss reaches past run_battery. C5's unstudentized row is also lost and unnamed | "every row of the cell, C3 included, since the ledger returns before calling it"; and put the C3 point in front of Sera with (a)
MINOR | src/quebra/analyzers/checks/c6_exchangeability.py:30 | new line "C5 reads the same `z_h`" repeats the class CP1 minor 4 fixed in _rank_serial.py:135: C5's unstudentized variant reads raw r_h, not z_h (c5_rank_autocorr.py:105 studentizes only one variant). Inferred from reading. Same imprecision in the unchanged line 32 | "C5's studentized variant reads the same `z_h`; both C5 variants are exact in the same sense"
(clean) src/quebra/analyzers/checks/c6_exchangeability.py:7-8,12-19,21-26,28-29 | CP1 minors applied: "whose absolute maximum C5 takes" matches c5:109; rewrap leaves no docstring line over 93 chars (awk); the uncentred-sum paragraph now says size stays exact and power is distorted by sign; "within each segment"; the uniqueness claim and the "natural reference for criterion 4" sentence are gone. Dufour-Roy locators (2.1) p.260, (2.2) p.260, p.261, Sect. 3, p.268, p.271 and Ljung-Box (1.1) p.297 match Resolved opens. No em dash, no "gate", no spec id | -
(verified) src/quebra/analyzers/instrument_validation.py:628-641 | CP1 minors on the NaN verdict and the verb applied; mutants on scratch copies (mut_iv.sh, PYTHONPATH=copy:$PWD, quebra.__file__ printed as the copy): verdict back to max() -> [nan] red only; verdict constant TIER_PASS -> both params red; verb constant "match" -> both red | -
MINOR | src/quebra/analyzers/instrument_validation.py:627,636 | the printed "max abs difference" is still `max()`, the NaN-skipping call the new comment at :628 warns about; VERIFIED (probe_iv.py): NaN at durations_trend/acf_rank_lag5 gives verdict fail with detail "r_h do not all match R's acf ... (max abs difference 4.2e-16)", a fail whose stated worst difference is float noise | `rank_acf_max = float(np.max([r.abs_difference for r in rank_acf]))` (propagates NaN) or report the count of non-finite rows beside the max
(verified) tests/test_checks_statistics.py:510-532 | new no-spread test: MUT=std0 (CP1 body, std-based `scale <= 0.0`, patched on rs/c5/c6 via sd2_work/review_cp11_cp2/mutplug.py) -> [199] and [999] both red; MUT=c6_only_std0 (C5 keeps the raise) -> both red on "c6 returned a result", so C6 is reached and checked independently. unit marker, analytic oracle named | -
(note) tests/test_checks_statistics.py:510-532 | MUT=ptp0 (ZERO_SPREAD_TOL = 0.0) leaves both params green. Expected, not a defect: by the comment's own bound no midrank lag range lies in (0, 1e-12], so the tolerance value cannot be told from 0.0 on any rank input; the test guards the switch from std to ptp, which is what Sera's "mutation back to <= 0.0" means. A banner reading "<= 0.0" should say it is the std form | say "std-based <= 0.0" in the banner
(verified) tests/test_instrument_validation.py:291-310,376-377 | new fail-branch test and the two added guard strings: dropping "the studentization" from C6's not_pinned turns the guard red; replacing "and the pairing across gaps" turns it red; the three verdict/verb mutants above turn the new test red as listed. Module pytestmark = statistical covers the new test; the 1e-6 offset is 1000x AGREEMENT_TOL (1e-9). The guard does not require "max"/"sum" in the not-pinned clause the plan R5 lists; CP1 asked for only the two strings | -
(verified) src/quebra/analyzers/checks/_rank_serial.py:171-175 (addendum to the message MINOR above) | through run_battery the message's one example is the case the CP1 code already caught: [1,2]+[3,4]+[5,6] and constant-inside records raised before too (std exactly 0, or C1/CvM gamma_hat). The case CP1.1 newly catches is all-size-2 records with real-valued durations: 8 random exponential size-2 segments gave C5, C5-raw and C6 p = 1.0 under the CP1 std check and raise now (probe_size2.py). So the reason a ledger row will actually carry names the wrong cause | as in the message MINOR
MINOR | plan_snapshot.md, R8 step 3, "C5 and the no-spread raise" | "Count them per cell (the rise in `n_failed`, and `first_failure` naming "no spread")": the table column is `fail_reason` (runner.py:254), it holds only the cell's FIRST failure (runner.py:155-159), and a retry that raises is never recorded (runner.py:199-201). VERIFIED (probe_bench.py): a constant-inside replicate in a tau-checks cell leaves fail_reason = "gamma_hat is zero ..." although the replicate never reached a computed C5. Any cell whose first failure predates its first no-spread replicate will not name "no spread" | "Count them by the rise in n_failed; fail_reason names only the cell's first failure, so attribute the rise by re-running the cell's failed replicates and reading each message"
(clean) plan_snapshot.md, R8 additions and CP3 "M12 note" | the R8 note says what Sera asked (C5 rows may differ only in cells with degenerate replicates; count them; any other C5 difference is a stop); the "Added by Claude, pending Sera" note is labelled as such and is true (probe_bench.py: every row of the cell loses the replicate). The CP3 M12 note carries CP1 minor 11 faithfully, including "re-measure before citing". The "9 hits" for `grep -rn pending docs/iid_checks/` is VERIFIED (9), and the listed line numbers are the pre-CP2 ones, as the note says | -
(verified) .claude/qre_checks_reference.tex | compiles: a copy with line 6 (newtxtext,newtxmath) removed, two pdflatex passes in sd2_work/review_cp11_cp2/tex/, exit 0 both; the same 7 overfull boxes and the same "Label(s) may have changed" warning as 3a9ab92's copy compiled the same way (all shifted by the 19 added lines; none in the C6 section or the C6 table row). PDF text of the C6 section reads correctly. Dufour-Roy (2.1) p.260, (2.2) p.260 with E[r_k] = -(n-k)/(n(n-1)), p.261, Sect. 3, p.268, p.271 and Ljung-Box (1.1) p.297 match Resolved opens; the Dufour-Roy bibitem volume/pages match (29(3):257-273). Old m_h statistic, "why m_h weighting", "natural reference" and "burden" are all gone from the file (grep) | -
MINOR | .claude/qre_checks_reference.tex:800 | "C5 reads the same $z_h$ and is exact in the same sense" - same class as the c6:30 MINOR: C5's unstudentized row reads raw r_h. The rest of the sentence ("as is every permutation row of the battery") is the randomisation hypothesis the file already cites for C5 (piece 2), and the bench does not contradict it (permutation rows' size z max 2.26-2.98 over 198-216 cells per row, consistent with chance) | "C5's studentized variant reads the same $z_h$"
MINOR | .claude/qre_checks_reference.tex:818-821 | "Degenerate lags": the one example ("segments each constant inside") is the case C1/CvM's gamma_hat raise already caught; the case CP1.1 newly catches, all-size-2 segments, is unnamed (see the _rank_serial message MINOR). The rest (C5-studentized and C6 raise, the battery raises, the ledger writes not computed) VERIFIED | add "or every segment pairing at that lag has two events"
(note) .claude/qre_checks_reference.tex:700 | the C5 section's "Null. Exchangeability of the duration sequence" now differs from C6's "within each segment" for the same permutation scheme. Not a C6 site, so outside R6 as written | optional: "within each segment" in C5's Null too
(verified) docs/iid_checks/BENCH.md:137-140 | "Three exact permutation rows are flagged in the current tables" VERIFIED: calibration_summary.bench_acceptance_at_n(size_table.csv, alpha=0.05) flags exactly C5 studentized in_spec n=700 (0.0345), C5 unstudentized in_spec n=20 (0.0355), CvM permutation calendar n=100 (0.0625); no C6 row | -
MINOR | docs/iid_checks/BENCH.md:139-140 | "so this list may grow" names only C6 as the source of change, but the R8 note expects C5 rows to change in cells with degenerate replicates, and the implementer's own count found its two zero-range replicates in n = 20 in_spec cells (sd2_work/mut_cp11/bench_spread.txt), which is where the flagged C5-unstudentized n = 20 row sits. Inferred: whether that row moves is unmeasured | "so this list may change" or "The C6 rows, and C5 rows in cells with degenerate replicates, are pending the re-run"
MINOR | docs/iid_checks/C5_rank_autocorr.md:69 | "C5 takes the maximum of the same `z_h` that C6 sums": C5 takes the maximum of |z_h| and C6 sums z_h^2 - the imprecision CP1 minor 11 removed from c6_exchangeability.py:7 ("Fix every site of a class") | "C5 takes the largest |z_h|, and C6 sums the squares of the same z_h"
MINOR | docs/iid_checks/C5_rank_autocorr.md:46-49 | "as on a record whose segments are each constant inside": same wrong-example issue as the _rank_serial message MINOR; on such a record the battery already raised through C1/CvM gamma_hat before CP1.1, and the case newly not computed is the all-size-2 layout (VERIFIED, probe_size2.py) | name the two-event-segment case
(note) docs/iid_checks/C5_rank_autocorr.md:64-66 | C5's own bench numbers (worst z = -3.18, mean power 0.357 at n = 100) are not marked pending although R8 expects C5 rows to change in degenerate-replicate cells. Likely unmoved (2 of 28,800 replicates at n <= 50 by the implementer's count, none at n = 100 checked), so a note for R8 step 6, not a defect | add C5_rank_autocorr.md:64-66 to R8 step 6's refresh list, conditional on any C5 row changing
(clean) docs/iid_checks/C6_exchangeability.md:5-31,57-67 | states Q = sum_h z_h^2 with z_h centred and scaled by the augmented permutation moments (matches c6.run:127-129); every old-statistic number (0.087..0.788, 0.017-0.018, 0.17, -2.77, 0.332, 0.0365) and the "natural reference ... no citation" claim are gone; Dufour-Roy (2.2) p.260, p.271 and Ljung-Box (1.1) p.297 match Resolved opens; "agreement with C5 ... check_agreement" is right (check_agreement compares every non-reference row, C5 included, against C6, report.py:322-364); "every permutation row is exact under the same null" agrees with _permutation.py:14-17 and is not contradicted by the size table | -
MINOR | docs/iid_checks/C6_exchangeability.md:52-55 | the "no spread" paragraph repeats two issues already logged against the docstring: its example (segments constant inside) is the case the battery already refused through C1/CvM gamma_hat, while the all-size-2 case it newly refuses is unnamed; and "C1, C2 and CvM included" leaves out C3, which the ledger also never runs (check_ledger.py:509-512 returns before :539). VERIFIED by probe.py and probe_size2.py | carry the same fix as c6_exchangeability.py:59-64
(verified) docs/iid_checks/LIMITATIONS.md:10-19,111-118 | C5 studentized row recomputed from power_table.csv (arm E, mean over rho and shape): 0.106/0.174/0.236/0.307/0.357/0.626/0.691/0.750/0.804, matches; "Below n = 100 C5 does not exceed 0.5 except at rho = 0.5" holds per shape too (only rho = 0.5 cells at n = 35, 50, 75 exceed 0.5, both shapes). Section 9 meets D12: names every permutation row, cites _permutation.py, whose docstring (:14-17) gives the same reason (moving a gap changes both tau_j and destroys exactness); "rejected no more often than the nominal size" is the exactness bound under within-segment exchangeability | -
(verified) docs/iid_checks/iid_checks_basics.md:100-110 | "C5 HOLD ... mean power at n = 100 is only 0.34-0.36" VERIFIED: studentized 0.3567, unstudentized 0.3412 (power_table.csv, arm E); 24% / 36% / 63% / 80% at n = 50/100/355/1000 match 0.236/0.357/0.626/0.804; "exceeds 50% only at the strongest dependence" holds per shape. C6 row reads "pending" | -
MINOR | docs/iid_checks/iid_checks_basics.md:39 | "They are not redundant." is now a flat claim about C5 and C6 with no number behind it: the old measured agreement (0.017-0.018, "two of them are close") is gone, and the new C6's agreement with C5 is pending by the next sentence's own account. Claims discipline (AGENTS.md s.4): write an unmeasured claim as the open question it is | "C1/C2 and the rank checks are not redundant. Whether C5 and C6 are, with the current C6, is pending the bench re-run"
(note) docs/iid_checks/iid_checks_basics.md:39-40 | "C5 and C6 read the same studentized rank autocorrelations" - C5's unstudentized row does not; same class as the c6:30 MINOR | "C5's studentized variant and C6"
(verified) jobs/bench/report.py:659-671 | new criterion-4 text says nothing false: "C6 is the reference row" (check_agreement default reference, :322-323), "exact under within-segment exchangeability, as every permutation row's is" (_permutation.py:14-17; the size table's permutation rows show worst z 2.26-2.98 over 198-216 cells, consistent with chance), "a convention, not a ranking of validity". `build()` RUNS on the frozen tree and the rendered text no longer contains "natural reference" or "no citation". No test reads report.py's prose (test_instrument_validation.py:618 reads promotion_report.md, which R8 regenerates) | -
MINOR | jobs/bench/report.py:58 | missed C6 site: DIRECTED_AT still gives C6's alternative as "any departure from exchangeability", rendered into promotion_report.md's criterion-2 heading, while CP2 rewrote the same claim in iid_checks_basics.md:36 to "rank autocorrelations at lags 1-5 are jointly extreme". Q sees only departures that move r_1..r_5, so "any" overclaims for both the old and new statistic; the two descriptions now disagree | "serial dependence spread over lags 1-5" (or the basics wording)
(verified) missed-site grep (CP2) | grep over *.py, *.md, *.tex, *.R in src, docs, jobs, tests, spec, AGENTS.md, README.md for n_pairs(h), m_h, pair-count weighting, 0.017-0.018, "nearly redundant", "natural reference", 33-36%, 21-24%, 0.33-0.36, -2.77, 0.0365, 60-63%, 79-80%: no hit outside promotion_report.md (exempt), the review ledgers and sd2-update.md. The 9 "pending" placeholders are the 9 the plan lists. cp2.diff is byte-identical to `git diff` of the seven files on the frozen tree. Only missed C6 site found: jobs/bench/report.py:58 (above). Remaining "portmanteau" uses (checks/__init__.py:12, c5:9, c6:1, tex:706,764) stay true of a sum over lags | -
## Commands run (frozen tree, PYTHONPATH=$PWD/src[:$PWD], quebra.__file__ = this worktree's src/quebra/__init__.py; mutants on scratch copies or by -p plugin only)
- ruff check . exit 0 | ruff format --check . exit 0 (196 files already formatted) | make types exit 0 (11 files) | make arch exit 0 (1 kept, 0 broken)
- full suite pytest -p no:cacheprovider: 1134 passed, 10 skipped, exit 0 (log sd2_work/review_cp11_cp2/full_suite.log). Against CP1's 1130 passed: +4 = 2 no-spread params + 2 tier-4 fail params, none removed
- pdflatex on a copy without line 6: exit 0 twice; same 7 overfull boxes as 3a9ab92
- mutants: std0 / c6_only_std0 red at B = 199 and 999; ptp0 green (expected); verdict-max, const-pass, verb-match, no-studentization, no-pairing-across-gaps all red as listed
- probes: probe.py (raises, ledger, bound), probe_bench.py (+ _old), probe_size2.py, probe_iv.py, all in sd2_work/review_cp11_cp2/
- added-line scan of both diffs: no em/en dash, no "gate", no spec/phase id, no banned word
## Verdict: SHIP (0 CRITICAL, 0 IMPORTANT, 13 MINOR, 5 notes; the 13 come from 9 distinct defects, some repeated across sites). The CP1 IMPORTANT is resolved and mutation-verified; all ten CP1 minors are applied and minor 11 sits in the CP3 notes. Every MINOR is wording in a docstring, doc, message, plan note or report label. One of them is new information for open decision (a): the ledger also never runs C3 on a no-spread cell.

# C6 (SD2) patch, CP3 + CP4 review

SCOPE: CP4 (sd2_work/review_cp3_cp4/cp4.diff == git diff HEAD: jobs/bench/{grid,arms,runner,report}.py, docs/iid_checks/BENCH.md) + CP3 (sd2_work/mut/: sd2_mutants.py, run_mutant.sh, env.sh, summarize.py, logs/, logs_seq/; claims in sd2_work/review_cp3_cp4/cp3_claims.md)
COMMIT: 888de36
## Manifest
## Reviewed
- [x] jobs/bench/grid.py
- [x] CP4 re-measure: existing cells unchanged (index, seed, output)
- [x] jobs/bench/arms.py
- [x] jobs/bench/runner.py
- [x] jobs/bench/report.py
- [x] docs/iid_checks/BENCH.md
- [x] CP3 sd2_work/mut/sd2_mutants.py
- [x] CP3 sd2_work/mut/run_mutant.sh, env.sh
- [x] CP3 sd2_work/mut/summarize.py
- [x] CP3 claims 1-10 against logs
- [x] commands: ruff, format, types, arch, full suite
## Findings (C6 SD2 CP3 + CP4)
(verified) jobs/bench/grid.py:302-304 + runner.py main | existing cells keep index, seed and output: all_cells()[:504] field-for-field identical to HEAD's (cells_dump.py, HEAD overlay from `git show HEAD:jobs/bench/*` + current src, since the mirror's src is pre-commit-1 and carries the old C6); run_cell rows (TABLE_COLUMNS, 4 replicates) identical on 17 cells over arms A, C, D, E, both clocks, both kinds (139 rows); runner.main() on 10 size/power cells + 5 serial cells (RESULTS_DIR and all_cells monkeypatched, n_jobs=1) writes size_table.csv and power_table.csv BYTE-identical to HEAD's main() on the same 10 cells (main_dump.py). Counts: 216 + 288 + 87 = 591 cells, serial replicates 96,000, total 816,000 = the plan's figures. jobs/bench/results/ untouched (mtimes 15:08) | -
(verified) jobs/bench/grid.py:48-52 | the single-shape reason holds: CV^2 1.8302 (shape 0.75) and 0.4610 (1.5) match "about 1.8 against 0.46" and the module's own CV 1.35/0.68; at expected 5 events per segment (n 50, m 10, rho 0) the per-segment count variance is 7.07 against 2.36, and 52.5% of shape-0.75 replicates realise fewer than 10 segments against 5.6% at 1.5 (probe_arms.py, 2000 reps). The monotone-map argument is right: for a fixed seed the latent draw does not depend on shape | -
(note) jobs/bench/grid.py:46, arms.py:12 | "`phi = 0` is a null cell at every `delta`" is unqualified, while arm_el (arms.py:515-516) and the report say null for the PERMUTATION rows; EL runs and prints every row, C1 asymptotic included. Not a defect on measurement: probe_el_null.py, the real EL cells at phi 0 with 600 replicates, gives C1 asymptotic 0.0583 (se 0.0096) at delta 0.4 and 0.0400 (0.0080) at 0.8, every permutation row 0.040-0.063, C6 0.045 and 0.048 | for consistency, "a null cell for every permutation row at every `delta`"
MINOR | jobs/bench/grid.py:254-257 | `censoring_for_segments` promises the round trip through `n_segments_for_censoring`, which silently caps m at n // 5; a future SERIAL_SEGMENT_GRID or SERIAL_N_GRID entry past the cap would yield fewer segments than the grid names, and the report (which inverts through the same function) would print the capped target with no sign that the requested one was refused. Holds today (VERIFIED: 3, 5, 10 at n 50 and 60 all round-trip) | in serial_cells, raise when n_segments_for_censoring(n, censoring) != n_segments
IMPORTANT | jobs/bench/arms.py:459-461 | arm_e2's "It is the alternative where a check that reads lag 1 alone sees nothing and the lag-2 coordinate carries all of it" is false. The latent lag-1 correlation is zero, so r_1 has mean ~0, but Bartlett's formula gives Var(r_1) ~ (1 + 2 rho + 2 rho^2)/n for an MA whose only nonzero autocorrelation is rho at lag 2 (2.5x the null variance at rho = 0.5), and every lag h > 2 gets (1 + 2 rho^2)/n. VERIFIED (probe_arms.py, single segment n = 300, 600 replicates, B = 399): a two-sided exact lag-1-only rank permutation test rejects 0.047 (se 0.009) at rho 0, 0.097 (0.012) at rho 0.2, 0.180 (0.016) at rho 0.5. Scenario: E2 power for C5 (max over lags 1-5) is read in the thesis as detection "at a lag other than 1", while part of it is the lag-1 and lags 3-5 coordinates firing on variance inflation | "the latent autocorrelation is zero at every lag but 2, so lag-1 statistics are centred at their null mean; their variance is not (Bartlett), so a lag-1-only test still rejects above its size". The header's "power for C5, C6 at a lag other than 1" and BENCH.md's "away from lag 1" can stay
(verified) jobs/bench/arms.py:428-475 | lag2_theta solves theta/(1+theta^2) = rho on the invertible root and raises past |rho| = 0.5 (RHO_GRID's 0.5 gives theta = 1, correlation 0.5); _latent_lag2 has nonzero autocorrelation at lag 2 only (shared innovation e_{t+2}), stationary from its first draw; arm_e2 is a faithful port of ../quebra-c6/c6_work/sim/arms_lag2.py (same theta, latent and build), with arm_e's marginal, truncation and censoring-as-segments | -
(verified) jobs/bench/arms.py:478-500 | the level-shift tau uses quadrature with norm.logsf as the plan requires: level_shift_mean_gap agrees with an independent quad over u in (0, 1) to <= 3.9e-13 at shapes 0.75 and 1.5, delta 0, 0.4, 0.8 (1.0000000000000010, 1.5360083375366, 2.2746910355368 at shape 0.75), so "1 at delta = 0" holds to 1e-15. "both segments expect the same number of events": mean counts per segment 150.63/149.60, 150.27/150.28, 150.58/150.03 (se 0.37) at delta 0, 0.4, 0.8 over 2000 replicates (probe_arms.py) | -
(verified) jobs/bench/arms.py:513-517 | "phi = 0 ... null cells for every permutation row": the time-truncated renewal sample is exchangeable given N, so the claim is right; measured on the real EL cells at 600 replicates (probe_el_null.py): every permutation row 0.040-0.063 at delta 0.4 and 0.8, C6 0.045 (0.0085) and 0.048 (0.0088) | -
(verified) jobs/bench/runner.py:71-97,107-128,374-383 | TABLE_COLUMNS is untouched by the diff; size and power are reindexed to it after an `==` split, and main() writes them byte-identical to HEAD (main_dump.py, above). serial_table.csv's header (main_dump.py output and the implementer's smoke table) carries every plan column: delta beside b; n_segments_min, mean_n_segments, n_segments_max; segment_size_min/q25/median/q75/max; n_not_computed beside n_failed_tau_only. An unknown kind raises instead of being written nowhere. n_not_computed = n_replicates - n_used matches "replicates a row has no p-value for": the battery never returns p_value None (only C3 does, c3_serial_copula.py:168, and run_battery does not call it) | -
IMPORTANT | jobs/bench/runner.py:15-16,322 (serial cells) | serial_table.csv inherits `mean_induced_lag1`, documented as "the duration-level dependence it actually produced, which is the quantity a check sees", but it is C5's raw r_1 (c5_rank_autocorr.py:128), centred on the whole record: exactly the uncentred-r_h bias this patch removes from C6. On the new cells it reads the layout as dependence. VERIFIED in the implementer's smoke table (60 reps): EL phi 0 reads 0.010 / 0.042 / 0.126 at delta 0 / 0.4 / 0.8, and EL phi -0.12 at delta 0.8 reads +0.034, the WRONG SIGN; gapped E at rho 0 reads -0.016 to -0.056. A reader of the tracked serial table takes a level shift for positive dependence | drop mean_induced_lag1 from SERIAL_COLUMNS, or record the studentized z_1 mean instead, or at least state in the runner docstring that on the serial cells it is r_1 centred on the whole record and moves with the layout and the level shift
MINOR | jobs/bench/runner.py:178-185 | the ARM_EL branch ignores cell.quantised and cell.censoring; a future EL cell built with quantised=True or a censoring target would run continuous and two-segment with its label and table row saying otherwise (inferred; no such cell today, all 9 EL cells have quantised False, censoring 0) | raise in the EL branch when cell.quantised or cell.censoring is set
(note) jobs/bench/runner.py:286-288 (pre-existing class) | a row that no replicate computed (e.g. every replicate of a full-battery cell going through the tau-only retry) is absent from the table, not present with n_not_computed = n_replicates, so serial_layouts' n_not_computed_max cannot see it; only the all-rows-missing case gets the fallback row. Unlikely on the continuous serial arms | none needed for CP4; record if R8 shows a serial cell with fewer rows than its siblings
(verified) jobs/bench/report.py:127-136,551 | the report raises when serial_table.csv is missing: with RESULTS_DIR monkeypatched to a scratch copy of the committed size/power tables only, load_serial() and build() both raise FileNotFoundError (probe_report.py). With the implementer's smoke serial table, build() renders the section, no (cell, row) pair is duplicated (so aggfunc="first" hides nothing), and every SERIAL_LAYOUT value is constant within a cell (so groupby.first() hides no disagreement). Known consequence, plan-intended (D13): build() fails on the committed tree until R8 writes serial_table.csv; nothing in tests/ or src/ calls it | -
MINOR | jobs/bench/report.py:781-812 | the serial tables print every row's rate whatever its support, while the same report excludes size/power rows below MIN_SUPPORT_FRACTION "because a rate computed from the surviving replicates is conditional on whatever killed the others" (:819-823). The serial n = 50, m = 10 cells are where dropped segments and degenerate replicates concentrate (52.5% of shape-0.75 replicates lose a segment, probe_arms.py), and n_not_computed_max sits in a different table from the rate (inferred: no serial row under-supported in the 60-rep smoke) | mark or move rows under MIN_SUPPORT_FRACTION as the size/power section does, or print n_used beside each rate
MINOR | jobs/bench/report.py:786-787 | "a target of 3, 5 and 10 segments at n = 50 and 60" is hand-typed grid content in a report whose docstring says regenerating "cannot leave a stale figure in the prose"; a grid change leaves it stale | format it from grid.SERIAL_SEGMENT_GRID and SERIAL_N_GRID (the report already imports N_GRID)
MINOR | jobs/bench/report.py:806-807 | "Segments short of two events are dropped" holds for Arms E and E2 (_direct_segments) but Arm EL raises instead (arms.py:529-533), failing the replicate | "In Arms E and E2 segments short of two events are dropped ...; Arm EL fails the replicate"
MINOR | jobs/bench/report.py:4 | edited docstring line is 106 characters; the rest of the docstring wraps near 90 (E501 not selected, ruff green) | rewrap
(verified) docs/iid_checks/BENCH.md:13,22-23,33-38 | "Seven arms", "87 more" (42 E + 36 E2 + 9 EL, grid.serial_cells), "150 expected events per segment", delta and phi grids, shape 0.75, and the column list all match the code; the EL row's "part of every permutation row's null (phi = 0)" is the qualified form. The "504 cells and 720,000 replicates" at :51 still describes the last run; R8 step 6 owns it | -
MINOR | docs/iid_checks/BENCH.md:33-34 | "which only `promotion_report.md` reads": a markdown file reads nothing; the reader is jobs/bench/report.py (grep: no other reader in src, jobs, tests, scripts) | "which only `jobs/bench/report.py` reads, for `promotion_report.md`"
MINOR | jobs/bench/runner.py:1 | missed site: the summary line still says "write the size and power tables"; main() now writes three | "Execute the grid and write the size, power and serial tables."
(note) AGENTS.md section 6 "five generating arms" | left for Sera per the brief | -
(note) jobs/bench/{arms,grid}.py | no test touches lag2_theta, level_shift_mean_gap, arm_e2, arm_el or serial_cells; the plan did not ask for one and the bench's existing arms are untested too (arms.py module docstring). The probes in sd2_work/review_cp3_cp4/ are the only evidence for the quadrature and the EL counts | Sera's call
(verified) CP3 sd2_mutants.py, every mutant against plan "CP3: mutation" | M1-M8 bodies identical in substance to ../quebra-c6/c6_work/mut/autocorr_mutants.py (difflib: docstrings and the ORIGINAL/ORIG_AUTOCORR name only); M9 is 3a9ab92's c6.run body line for line (git show 3a9ab92); M10 is HEAD's studentize with the centre dropped and the ptp raise kept; M11 fully permutes the concatenated vector and keeps the per-segment labels; M12 hands resolve_perm one block of the total; M12r relabels a one-block shuffle with the segment sizes only when perm is None, as the plan's CP3 note asks | -
(verified) CP3 sd2_mutants.py:270-297,378-382, binding sites | install happens at plugin import, before any test module or late-imported pipeline module, so every later `from ... import name` binds the patched object; the sweep at import and after collection rebinds any same-name binding to the ORIGINAL. grep over src, tests, jobs, scripts finds no alias import, default-argument capture, partial or container capture of autocorrelations, block_permutations, studentize or resolve_perm, and every C6 call goes through the c6 module attribute (battery.py:212, assumptions.py via `_c6.`). Calls records confirm reach: M1 in 49 tests over 7 files incl. test_r_cross_implementation (own import, 3 tests) and test_instrument_validation (9); M9/M10 in T7 (2) and test_check_ledger (6); M11 in test_check_ledger (7), test_checks_cvm and T7; M12/M12r in T7 and test_check_ledger. Every run collected 1144 tests, so no mutant broke collection. Hole, not affecting any claim: subprocess-based tests (test_check_ledger.py:136, test_check_attachment.py:102) run without the plugin; none of them is C5/C6 evidence | -
(verified) CP3 summarize.py, can it hide a survivor | no: a mutant is SURVIVED when its only reds are the two race tests (exclusion is conservative), "ok" needs every EXPECTED tag red, and no mutant's red set was empty. NONE/M0's "ok; race red 2" is only cleared by the sequential logs_seq runs, which summarize.py does not read (claim 1 covers it) | -
MINOR | sd2_work/mut/summarize.py:41 | EXPECTED["M11"] is {"T7"} although the plan's table also names test_permutations_stay_inside_their_segment; had that test stayed green the summary would still print "ok". It is red in M11.xml, so nothing is hidden this time | add a pre-existing-test expectation (by base name) for M11
MINOR | sd2_work/review_cp3_cp4/cp3_claims.md claim 3 ("T7 under M5b" read as a real detection) | false: M5b replaces the observed and null rows with iid noise, which leaves the permutation p-value exactly uniform, so T7 (an exactness test) cannot detect it. The red was c6 at 0.0365 against nominal 0.05, 2.77 null SE out against a 2.64-SE Bonferroni band, a chance exceedance at that stream position. VERIFIED: T7 alone under M5b is green (2 passed), and green at 6 other noise seeds on a scratch copy of the plugin (sd2_work/review_cp3_cp4/mutcopy). The other extra reds are real: NS under M5b/M9/M11 (no raise without studentize or with cross-segment spread), T6 under M8 and M11 (its lag-1-centre assertion, test_checks_statistics.py:505, fires on the sign flip and on the cross-segment null) | banner: "T7 under M5b was a chance exceedance of its band; T7 cannot see M5b by construction"
(verified) CP3 claims 1, 2, 4, 5, 6, 7, 8, 9, 10 | 1: 16 runs in batches of three (exit timestamps 21:35 x3, 21:40 x3, 21:45 x3, 21:51 x3, 21:56 x3, 22:00 x1); logs_seq NONE and M0 exit 0, 1134 passed, 10 skipped. 2: every plan-named red is red, incl. M5a/M5b's "as before" sets (0 old reds lost) and M11's test_permutations_stay_inside_their_segment. 4: old harness M1-M4, M6-M8 each 1126 tests, 0 failed, exit 0. 5: the parallel failures read "the probe test did not run" and both tests pass in logs_seq. 6: M12's two T7 reds are "permutation set is blocked as (30,) but the segments are (15, 15)"; M12r's are rate failures c5 0.0950, c6 0.1425 at nominal 0.01, band +/-0.0059 (M11 gives the same two numbers, as expected: a one-block shuffle draws the same indices). 7: M11's reds contain no C2 or CvM test. 8: see binding-site entry. 9: git status --ignored identical before and after; the only files outside sd2_work/ written between .start_marker and git_status_after.txt are 361 in .hypothesis/constants and 1 in .pytest_cache. 10: proof files read M9 11.6987 (no z) vs 15.1247 (z1 2.2248); M10 21.0438; M11 6 of 6 rows cross; M12 battery RAISED; M12r battery 17.4879 -> 22.1912 | -
(note) CP3 claim 7 | the multi-segment C2 and CvM permutation tests cannot see a cross-segment shuffle; within-segment shuffling stays guarded by test_permutations_stay_inside_their_segment, T7 and NS, so this is a known blind spot of those two tests, not a survivor | -
## Commands run (frozen tree, PYTHONPATH=$PWD/src[:$PWD], quebra.__file__ = this worktree's src/quebra/__init__.py; nothing outside sd2_work/review_cp3_cp4/ and this ledger written)
- ruff check . exit 0 | ruff format --check . exit 0 (196 files already formatted) | make types exit 0 (11 files) | make arch exit 0 (1 kept, 0 broken)
- full suite pytest -p no:cacheprovider, run once and alone: 1134 passed, 10 skipped, exit 0 (sd2_work/review_cp3_cp4/full_suite.log); identical to CP2's count, as expected for a bench-only diff
- probes: cells_dump.py, main_dump.py, probe_arms.py, probe_el_null.py, probe_report.py; T7 under M5b alone and at 6 seeds on a plugin copy
- jobs/bench/results/ untouched (all mtimes 15:08-15:23); .claude/review-findings.md is TRACKED, so `git status` now shows it modified by this ledger: leave it out of the CP4 commit or commit it separately
- added-line scan of cp4.diff: no em/en dash, no "gate", no banned word, no spec/phase id
## Verdict CP4: DO NOT SHIP (0 CRITICAL, 2 IMPORTANT, 8 MINOR). Both IMPORTANTs are text-sized: the arm_e2 docstring's lag-1 claim, and mean_induced_lag1 in serial_table.csv (drop it or document it). The mechanics the brief asked about all hold: existing cells keep index, seed and byte-identical output, TABLE_COLUMNS unchanged, serial columns complete, logsf quadrature correct, report raises on a missing table.
## Verdict CP3: claims hold except claim 3's "T7 under M5b" (a chance exceedance, not a detection); no mutant survived; the harness reaches every binding site and summarize.py cannot hide a survivor.

# C6 (SD2) patch, CP4.1 + R10 review

SCOPE: CP4.1 + R10: sd2_work/review_cp41_r10/delta.diff (working tree vs cp4_tree = 888de36 + cp4.diff; 21 files, 881 diff lines); requirements.md; claims.md 1-18
COMMIT: 888de36
## Manifest
## Reviewed
- [x] R10 schema consumers (ROW_KEYS, SURVEY_KEYS, ALL_KEYS, PERMUTATION_KEYS, row_key, C5-raw) across src, jobs, tests, scripts
- [x] src/quebra/analyzers/checks/battery.py
- [x] src/quebra/analyzers/checks/c5_rank_autocorr.py, _rank_serial.py, c6_exchangeability.py, cvm_cramer_von_mises.py
- [x] src/quebra/analyzers/independence_survey.py, check_outcome.py, calibration_summary.py
- [x] jobs/bench/runner.py, grid.py, arms.py
- [x] jobs/bench/report.py
- [x] bench control (claims 7, 13): evidence power and re-measure
- [x] claims 8, 9, 10 (banner reports, CP3 record)
- [x] docs/iid_checks/{BENCH,C3_serial_copula,C5_rank_autocorr,CvM_cramer_von_mises,LIMITATIONS}.md
- [x] tests/test_bench_serial_cells.py, test_checks_statistics.py, test_check_outcome.py, test_checks_cvm.py (+ mutations)
- [x] commands: ruff, format, types, arch, full suite
## Findings (C6 SD2 CP4.1 + R10)
(verified) R10 consumers, src | grep over src, jobs, tests, scripts for ROW_KEYS, BENCH_ROW_KEYS, SURVEY_KEYS, ALL_KEYS, PERMUTATION_KEYS, row_key, C5_RAW_KEY, VARIANT_RAW, "unstudentized", include_c5_raw. Each consumer: check_ledger.py:487,498 calls run_battery with the default, so neither the first call nor the tau-only retry produces the raw row; its acceptance lookup (:578-584) is keyed per result, so the bench's raw acceptance rows are simply never matched. check_outcome.py:44 ALL_KEYS = SURVEY_KEYS (9 = 8 + C3), :65 PERMUTATION_KEYS 6 -> 5, TAU_DEPENDENT_KEYS unchanged, normalise() now raises on the raw key; no job names it (jobs/active only uses PERMUTATION_KEYS and two explicit C1/C6 keys). independence_survey.py:53 SURVEY_KEYS 10 -> 9, the raw entries leave CHECK_LABELS and CHECK_NULL together (test_every_surveyed_key_has_a_label_and_a_stated_null is relative and still binds). plots/independence_survey_plot.py:186 SURVEY_PLOTS is generated from SURVEY_KEYS, so the raw survey figure class disappears; nothing imports it by name; ledger_recipe.py:189 loops SURVEY_PLOTS. panels/_check_ledger_render.py builds one column per row present, so absence is handled. calibration_summary.py:254 and plots/calibration_plot.py:104 read the BENCH tables, where the raw rows stay. scripts/: no consumer. No min-P anywhere ever read C5 (min_p is imported only by tellability.py) | -
(verified) R10 consumers, positional | the only positional uses are tests/test_check_outcome.py:306-337 (ROW_KEYS[0], [1], [3], [7]). [0], [1], [3] are unchanged (c1 asymptotic, c1 permutation, cvm asymptotic). ROW_KEYS[7] moved from c6_exchangeability to cvm_cramer_von_mises/permutation; test_the_two_are_independent only needs a key other than [1], so it still tests what it says. No src or jobs code indexes ROW_KEYS by position | -
(verified) R10 consumers, ordering | the bench runner keys rows with row_key and writes them by sorted(seen.items()) (runner.py:283,307), so moving the raw row from position 7 to last changes no table order. The battery raises all-or-nothing, and raw cannot raise where the studentized row did not (same validate_segment, no studentize), so the order change also cannot change a fail_reason string | -
(verified) R10 consumers, pins re-pinned? | every count change in tests is a deliberate re-derivation, not a copy of output: test_check_outcome's 5 / 8 / 7 / gap 2 follow from ROW_KEYS minus raw and battery_flags(PERMUTATION_KEYS) (include_tau_checks on, C2-asymptotic off, C3 off); the tau test's len == 4 became the three named keys in the docstring's order. Pins written relative to ROW_KEYS (test_check_ledger.py:328,503; test_independence_survey.py:344-401) did not need to change. Requirement listed test_independence_survey, test_instrument_validation:210 and test_assumptions as places to update: none hard-codes the row set (test_assumptions counts six check MODULES, unchanged), so leaving them is right | -
(note) requirements R10 "jobs/bench/runner.py ... uses the bench list" | the runner does not import BENCH_ROW_KEYS; it never used ROW_KEYS either, it derives every row from row_key(result). Nothing to change; BENCH_ROW_KEYS's only consumers are the new test and docstrings | -
(note) Known risk to carry into the banner | km_with_checks_6d2s's RUN_SET = PERMUTATION_KEYS is a step kwarg, so its run identity changes and every KM band's checks_asked loses "c5_rank_autocorr permutation unstudentized"; the survey job writes one fewer figure and n_instruments reads 9. Old outputs keep the raw row; R9 re-runs them | -
IMPORTANT | src/quebra/analyzers/checks/c5_rank_autocorr.py:13-18 (edited by the delta), _rank_serial.py:30-31 (edited sentence's context), docs/iid_checks/C5_rank_autocorr.md:41 | the docstring now says both things: the retained paragraph claims that on an ungapped record "the two variants are near-identical by construction and any bench cell showing them differing would indicate a bug rather than a finding", and the paragraph this delta adds (:30-41) says the raw variant is uncentred by -(n-h)/(n(n-1)) and its power leans with the sign of the dependence. The shipped power table sides with the new paragraph. VERIFIED on jobs/bench/results/power_table.csv: every Arm E power cell is uncensored (one segment), and raw is LOWER than studentized in 84 of 88 non-tied cells (sign test p 1.6e-20), up to 7.6 points at n 20, rho 0.5 (0.240 vs 0.317, conservative unpaired z -3.81); Arm D 63 of 65 lower. A reader who follows the docstring reads that difference as a bug. The 1.04 SD ratio argument covers the SCALING only; the centring is not a no-op at small n (-0.05 against a null SD near 0.22 at n 20). This is the class of claim the requirement's docs bullet and claim 17 asked to reconcile | in all three sites, keep the SD-ratio sentence but scope it to the scale: "the per-lag scales nearly agree on an ungapped record (SD ratio 1.04), so studentizing's SCALING is nearly a no-op there; its CENTRING is not, and the raw variant's power is lower for positive dependence (power table, Arm E: raw below studentized in 84 of 88 non-tied cells)". Drop "would indicate a bug"
(verified) src/quebra/analyzers/checks/c5_rank_autocorr.py:30-46 | the new uncentred paragraph: E[r_h] = -(n-h)/(n(n-1)) on one segment is right (the pairwise covariance of exchangeable centred ranks is -S/(n(n-1)), times n-h pairs), it matches the existing Eq. (2.2) citation in c6 and C6.md, and the pinned test test_the_gapped_null_mean_is_the_closed_form exists (tests/test_checks_statistics.py:440). c5.run draws no randomness when `autocorr` is passed (:111-118), so "turning it on changes no other row" holds. extra["r_lag1"] is observed_r[0], the uncentred r_1, as the comments now say; C6 reports signed z_h in extra["z"] (c6_exchangeability.py:148) | -
(verified) src/quebra/analyzers/checks/battery.py | ROW_KEYS 8 rows, BENCH_ROW_KEYS = ROW_KEYS + C5_RAW_KEY, run_battery's result order is ROW_KEYS order then raw last, matching BENCH_ROW_KEYS; "leaving three, or four with include_c5_raw" is right (8 minus the five tau rows; C2 and CvM asymptotic already absent for m > 1). C3 has no BENCH_ROW_KEYS entry and no row in any table (grep size, power, serial) | -
(verified) cvm_cramer_von_mises.py:44, calibration_summary.py:369-371, check_outcome.py:43-48, independence_survey.py:16-21,44-54 | each reworded sentence matches the code: every ROW_KEYS row is a bench row; "five of eight"; "no ledger carries that row"; "C5 appears once, studentized" | -
(note) battery.py:96 | row_key's docstring says it returns "the BENCH_ROW_KEYS entry a result belongs to"; a C3 result (check_ledger builds one) maps to a key in neither tuple. Harmless | "the (check, calibration, variant) key of a result"
(verified) jobs/bench/runner.py:111-133 | SERIAL_COLUMNS leaves out mean_induced_lag1 and keeps every other CP4 column in place; TABLE_COLUMNS is untouched, so the size and power tables keep it (as the requirement asks) | -
(verified) jobs/bench/runner.py:15-18 | "on a single-segment cell ... a cell with no dependence reads about -1/n": on the committed size table's single-segment cells (A, B, C, D, E), mean_induced_lag1 sits within 2.7 SE of -(N-m)/(N(N-1)) everywhere (SE 1/sqrt(N n_used)). runner.py:287 reads r_lag1 from the studentized row, which is observed_r[0], the same uncentred value the raw row carries | -
IMPORTANT | src/quebra/analyzers/calibration_summary.py:76-81,210-230 (outside the delta's lines; sibling site of the CP4 finding this delta answers) | power_vs_dependence plots C5 and C6 power against mean_induced_lag1 and its docstring still says it is "the dependence a check sees". The delta documents the -1/n offset at the producer (runner.py:15-18) but not at the one consumer that draws it. VERIFIED on power_table.csv, Arm E: at n 20 the plotted x for rho 0.05 is -0.0096 (negative) while -1/N is -0.0497, so the offset is 124% of the signal; 68% at rho 0.10; 2% at n 1000. The curves for different n are therefore shifted against each other by up to 0.05 on an axis whose smallest grid step is 0.05. AGENTS.md section 4, "fix every site of a class" | in power_vs_dependence add the single-segment null mean back (x + 1/mean_n_events, the power cells are all one segment, which the function can assert), or relabel the axis and docstring "uncentred r_1, null mean -1/n". Tables unchanged either way. Sera's call whether that lands in CP4.1 or is recorded as owed
MINOR | jobs/bench/runner.py:15-18, docs/iid_checks/BENCH.md:39-40 | the reason given for dropping the column from the serial table ("across several segments the uncentred r_1 also moves with the layout") holds just as strongly for the 96 multi-segment size cells (Arms A and D, censored), which keep the column with only the single-segment caveat. VERIFIED: on those iid-null cells mean_induced_lag1 reaches -0.0908 against -1/N of -0.061, up to 48 SE from the exchangeable mean -(N-m)/(N(N-1)) (the truncated segments are exchangeable within, not across). A reader of the size table takes that for negative dependence | one more sentence: "on a censored size cell it also moves with the layout, as on the serial cells, and is not a measure of dependence there"
MINOR | jobs/bench/runner.py:186-191 | the EL refusal raises inside generate(), which run_cell catches per replicate (:241-245), so a misconfigured cell does not stop the bench. It becomes one all-failed row (VERIFIED, probe_el_e2.py: n_used 0, n_failed 3, check None, fail_reason naming the refusal) that the report prints as an "excluded (n 0)" column labelled "nan [nan]". It is visible, so this is not silent, but a configuration error is not a replicate failure: it is counted once per replicate and would be read as an arm defect. The implementer disclosed it as a known risk | raise it once, before the replicate loop (top of run_cell, outside the try) or in grid.serial_cells, so a bad cell aborts main()
MINOR | jobs/bench/arms.py:460-463 | "which inflates the variance of a lag-1 statistic: a test that reads lag 1 alone still rejects above its size" holds only for rho > 0. arm_e2 accepts signed rho (lag2_theta allows |rho| <= 0.5), and Bartlett's factor 1 + 2 rho + 2 rho^2 is below 1 for rho in (-1, 0). VERIFIED (probe_el_e2.py, n 300, 1500 reps): var(sqrt(n) r_1) is 0.570 at rho -0.4 (Bartlett 0.520), 1.051 at 0, 2.113 at +0.4 (2.120). The grid uses only positive rho, so no table is misdescribed | "for rho > 0, which is every cell on the grid, this inflates ...; a negative rho deflates it"
(verified) jobs/bench/arms.py:460-463, rest | the corrected E2 claim matches the CP4 finding: a lag-1 pair is independent (zero latent lag-1 correlation, jointly Gaussian), successive pairs share a duration and are linked through lag 2 (Cov = rho), so Var(r_1) grows by 2 rho + 2 rho^2 | -
(verified) jobs/bench/grid.py:263-275 | censoring_for_segments raises on any mismatch, not only the cap (also rounding); it runs at grid construction, so a bad SERIAL_SEGMENT_GRID aborts all_cells() before anything runs. The shipped grid round-trips (the new test's first loop) | -
(verified) jobs/bench/report.py:176-199,797-863 | on real-shaped data the marking is right. probe_report.py builds the report from the committed size and power tables plus the CP4 60-rep smoke serial table with four injected rows (an E C6 row at n_used 3, an EL C1-asymptotic row at n_used 1, an E row at exactly 50%, an all-failed fallback row as runner.py:337-358 writes it) and diffs it against the CP4 report.py on the same tables. Result: exactly the two sub-threshold rows print "excluded (n 3)" / "excluded (n 1)"; the row at n_used 30 of 60 prints a rate, consistent with with_support's >=; every other rate gains ", n 60" and is otherwise unchanged; the excluded list keeps all 64 size rows (n_target still prints as an int, NaN rho and delta print "-") and appends the three serial rows. The index-based isin is safe: every caller passes a filtered read_csv frame with a unique index. The grid sentence renders "3, 5 and 10 segments at n = 50 and 60"; the EL sentence and MIN_EVENTS_PER_SEGMENT (arms.py:297 drops, :531 raises) match the code | -
MINOR | jobs/bench/report.py:116-119 with :839 | the all-failed fallback row (check, calibration, variant all None) now reaches the excluded-rows table as "nan [nan]" with rejection_rate "-" (VERIFIED, probe_report.py line 362 of the generated report). It is the one row that says a whole cell produced nothing, and it reads like a formatting accident; it is also absent from the rank-rows rate table, which filters by check | in _row_label, print "(every replicate failed)" when check is NaN
MINOR | jobs/bench/report.py:848-858 | serial rows in the excluded table are identified by censoring_target (0.0566 for 3 segments at n 50), while every serial table above keys the same cells by segments_target, so a reader cannot match an excluded row to its cell without inverting n_segments_for_censoring. Pre-existing for size rows and not new: the table has no shape or clock column, so two A_iid_weibull 355 / quantised / 0.25 rows that differ only in clock print as apparent duplicates | add segments_target (blank for size and power) or run the serial part through _serial_cells; add shape and clock
(verified, re-measured) bench control, claims 7 and 13 | the implementer's evidence could not detect the change that matters. rows_before.json and rows_after.json are byte-identical, but they hold run_cell AGGREGATES at 40 replicates: a rejection count out of 40 moves only when a p-value crosses 0.05, so a perturbed p-value that does not cross is invisible. The continuous columns (mean_induced_lag1, mean_n_events) catch a generator or stream change, not a p-value change. Both files also carry the same quebra path, so which tree wrote "before" is unverifiable after the fact. The main() pair (cp41/main_new vs r10/main_after, byte-identical, 15 cells x 3 reps) is weaker still. RE-MEASURED at the replicate level (rev/pvals_dump.py): run_cell's loop (seed, generate, perm, battery, tau-only retry) recording every result's repr(statistic), repr(p_value) and repr(extra r_lag1), run once under cp4_tree (pre-R10 battery, raw inside ROW_KEYS) and once under the live tree with include_c5_raw=True. Sample: 304 cells (every cell with n <= 60 plus every sixth larger one), all seven arms, both clocks, size/power/serial, 18 quantised, 126 censored, n 20-1000, 6 replicates each. Result: 1,824 replicates, 14,290 rows (1,822 of them C5-raw) identical to the last bit; the 2 failing replicates (cell 53, gamma_hat zero) fail identically. AST comparison with docstrings stripped: c5, c6, cvm, _rank_serial, calibration_summary, check_outcome and arms.py are code-identical to CP4, so battery.py and the runner's flag are the only code paths R10 changed | -
(note) bench control, the one path it cannot reach | the tau-only retry (runner.py:258-275) never succeeds anywhere: n_failed_tau_only is 0 in every row of the committed size and power tables and of the CP4 smoke serial table, and in the sample above the retry was entered twice and raised twice (the same gamma_hat-zero error, because CvM's permutation row, which the retry keeps, needs gamma_hat too). So the retry call's include_c5_raw=True has no observable effect in any table today, and no test covers it (mutation B2, below). The retry comment's "C5 and C6 never touch gamma_hat" omits CvM, which is in the retry and does; pre-existing, out of this delta's scope | record it; if the retry is meant to work, it has to drop CvM too or catch per check
(verified) claim 8, consumers of C5-raw before R10 | matches the grep in the consumer entry above: the ledger (verdict, and a column in _check_ledger_render), the survey (grid, generated figure class, overview row), the KM band (checks_asked only; km_with_checks_6d2s's DISPLAY_SET is C1 and C6), the calibration size figure (calibration_plot.py:104 dashed line, bench tables), the bench report (verdicts, check_agreement). min_p is imported only by tellability.py, so no min-P ever had it | -
(verified, with a recipe gap) claim 9, Arm C offset | my independent recomputation from the committed tables gives the same numbers: 126 Arm C cells, mean_induced_lag1 + 1/N from -0.0198 to +0.0167, max |z| 3.16 with one cell past 3, mean z -0.06, Bonferroni 3.54 (null SD 1/sqrt(N) per replicate, as stated). Re-running sd2_work/cp41/arm_c_offset.py reproduces its file only up to line 63: the two z lines in arm_c_offset.txt are NOT printed by the cited script, so the banner's recipe cannot regenerate them. All Arm C cells are one segment (mean_n_segments 1.0), so the -1/N correction is the exact one there. "C5 studentized rates flat at nominal" holds (0.048-0.049 mean by clock, rho > 0) | add the z computation to arm_c_offset.py so the cited script prints what the banner quotes
(verified) claim 10, CP3 record | sd2_work/mut/summarize.py now adds the pre-existing M11 expectation (EXPECTED_OLD) and labels T7 under M5b as chance; re-running it reproduces summary_after_review.md byte for byte (M5b "ok (+NS) (+T7 by chance, not a detection)", M11 "ok (+NS,T6)"). Positive control: with a second, never-red name added to EXPECTED_OLD["M11"] on a copy, M11 reads "expected red missing", so the new expectation can fail. cp3_claims.md carries the correction | -
(verified) docs/iid_checks/{BENCH,C3_serial_copula,CvM_cramer_von_mises,LIMITATIONS}.md, C5_rank_autocorr.md:44-52 | every changed claim matches the code: serial_table.csv's only reader is jobs/bench/report.py (grep over src, jobs, tests, scripts, docs); C3 has no row in any table; "two of the eight", and every ROW_KEYS row is a bench row; LIMITATIONS gains section 10 after 9, which C5.md cites correctly; "C6 and C5's studentized row centre each lag at its own permutation mean" matches studentize (_rank_serial.py:177-179, centre over observed plus permutations); the Dufour-Roy locator matches the one already in C6.md and c6's docstring. Sweep of docs, README, AGENTS.md, src, jobs/active, jobs/composite for "nine rows", "ships both", "both variants", "schema of the bench": nothing stale left (BENCH.md:139,149 name C5 unstudentized bench results, which still exist; AGENTS.md:248 "not in ROW_KEYS" is still true for C3) | -
(verified) writing rules, all 339 added lines of the delta | no em or en dash, no "gate" in any form, no surfacing/data-driven/delve/leverage/firstly/excellent, no spec or phase identifier (R10, CP4.1, SPEC) in a source docstring. The new statistical prose (c5 docstring, LIMITATIONS 10) cites its formula to an equation and page; c5's module carries no Validity-assumptions section, but none was added or removed by this delta | -
(note) src/quebra/plots/independence_survey_plot.py:44 (not in the delta) | "subclasses set KEY to one entry of battery.ROW_KEYS"; the C3 subclass's key is not in ROW_KEYS. Pre-existing | "one entry of independence_survey.SURVEY_KEYS"
(note) no new long .md doc; .claude/sd2-update.md is untracked and outside src/docs (the implementer's notes), not reviewed here | -
(verified, re-run) mutations against the SHIPPED tree | rev/mut/revmut.py loads a textually mutated copy of one module under its real name before any import (proof_<M>.txt written per run; nothing tracked edited); each run is the nine test files that touch battery, runner, grid, the ledger, the survey, check_outcome, CvM, IV and assumptions (231 tests). NONE: 231 passed, exit 0. D (include_c5_raw default True): 6 red: battery_returns_every_row, c5_raw_runs_only_when_the_bench_asks, battery_drops_only_the_tau_checks, switches_restrict, two check_ledger R-absence tests. K (raw back in ROW_KEYS): 11 red (the two battery tests, two check_outcome, two check_ledger, five independence_survey). RAWMID (raw computed at its old position, not last): 6 red. B (runner's first call drops the flag): 1 red, the new test. R (EL raise removed): 2 red, both parametrisations. R2 (censoring half removed): 1 red, [censoring]. G (censoring_for_segments raise removed): 1 red. Every claimed kill holds | -
MINOR | sd2_work/review_cp41_r10/claims.md claim 15 (banner) | the recorded counts describe a test file that never shipped: D "5 red", B "2 red", K "12 red" against the shipped tree's 6, 1 and 11. tests/test_checks_statistics.py was last written at 01:41:53, 29 s after the last of mut_D/B/K (01:41:04-01:41:24). At that moment test_battery_drops_only_the_tau_checks_when_asked still asserted len == 4, which the R10 battery (3 rows) fails with NO mutation: the baseline was red. That is why it is "red" under B and K (not a detection) and "green" under D (the mutation restored the fourth row). Conclusions survive; the record does not. AGENTS.md section 7: a recorded mutation must describe the tree it guards | banner: replace the counts with D 6, B 1, K 11 (this ledger), and run a NONE baseline in the same batch as any mutation
MINOR | jobs/bench/runner.py:271 (tau-only retry) | mutation B2, the retry call dropping include_c5_raw=True, SURVIVES: 231 passed. "Fix every site of a class" covers the guard too: the new test reaches only the first call. The scenario is a replicate that takes the retry and so writes every row but C5-raw, leaving the raw row's n_used below its siblings with no failure counted. The exposure today is nil, because the retry never succeeds on any committed cell (bench-control note above) | in test_c5_raw_runs_only_when_the_bench_asks_for_it, also run run_cell on a cell built to take the retry, or monkeypatch run_battery to raise once on include_tau_checks=True and assert C5_RAW_KEY is still written; or collapse the two calls so the flag is passed once
(note) jobs/bench/runner.py:116-117 | mutation SCOL, mean_induced_lag1 put back into SERIAL_COLUMNS, survives (231 passed). The requirement asked tests only for minors 1 and 2, so this is not a defect, but nothing pins the CP4 IMPORTANT's fix | optional: assert "mean_induced_lag1" not in runner.SERIAL_COLUMNS and in TABLE_COLUMNS
(verified) tests/test_bench_serial_cells.py | marker unit (one tier marker per test); the first test's monkeypatch of grid.SERIAL_SEGMENT_GRID reaches serial_cells because it reads the module global at call time; the EL test asserts the valid cell is served (n_segments 2) before the two refusals, so it cannot pass on a generate() that refuses everything. The file is a distinct subject from test_bench_uses_real_carve.py (contract of the grid and runner, not the carve), so section 7's one-file rule does not apply | -
(verified) tests/test_checks_statistics.py:557-601 and :651-659 | the new test's bit-identity loop compares every shared row (statistic, p_value) and the raw row against c5.run standalone, both tau settings; the tau test's explicit three-key list is the docstring's contract, in ROW_KEYS order, and fails under D and RAWMID. Neither was re-pinned to output. test_check_outcome.py's 5/8/7 and gap 2 are derived (consumer entry above). test_checks_cvm.py changes a docstring only | -
(note) commit split and budget | CP4.1 (commit 2) and R10 (commit 3) sit in one uncommitted tree on 888de36. Splitting them needs hunk-level staging in files both touch: runner.py (docstring and SERIAL_COLUMNS are CP4.1, the two include_c5_raw lines are R10), c5_rank_autocorr.py (the r_lag1 paragraph and comment are CP4.1, the "runs only in the bench" sentences are R10), LIMITATIONS.md section 10 (one section, both). R10 touches about 15 files (8 src or jobs, 3 tests, 4 docs) against the requirement's "about 6-9": over the envelope, under the 2x STOP; the banner should say which ones the envelope did not foresee. Collect: 1148 against CP3's 1144, +3 for test_bench_serial_cells.py (CP4.1) and +1 for the new battery test (R10, inside +1..+2) | -
## Commands run (frozen tree; PYTHONPATH=<worktree>/src:<worktree>, PYTHONPYCACHEPREFIX=sd2_work/.pycache; quebra.__file__ = this worktree's src/quebra/__init__.py; jobs.bench.runner from this worktree)
- delta check: sd2_work/review_cp41_r10/delta.diff equals a fresh diff of the live tree against cp4_tree (hunk headers aside), and cp4_tree differs from HEAD only in the five CP4 files
- ruff check . exit 0 | ruff format --check . exit 0 (197 files already formatted) | make types exit 0 (11 source files) | make arch exit 0 (1 kept, 0 broken)
- full suite, run once and alone: pytest -p no:cacheprovider, 1138 passed, 10 skipped, exit 0 (rev/full_suite.log); --collect-only 1148
- mutations (rev/mut/run.sh, nine files, 231 tests, sequential, never beside a full suite): NONE exit 0; D 6 red; K 11 red; RAWMID 6 red; B 1 red; B2 0 red (survives); R 2 red; R2 1 red; SCOL 0 red (survives, untested by design); G 1 red
- probes: rev/pvals_dump.py under cp4_tree and live tree (exit 0 both, 304 cells x 6 reps, identical); rev/probe_report.py (report build on injected tables, old vs new); rev/probe_el_e2.py; power and size table analyses (C5 raw vs studentized, mean_induced_lag1 by layout, Arm C z); sd2_work/cp41/arm_c_offset.py re-run (exit 0); sd2_work/mut/summarize.py re-run (exit 0, identical) and a positive-control copy
- git status --porcelain --ignored identical before and after; written outside sd2_work/review_cp41_r10/: this ledger, .hypothesis/constants (20), .mypy_cache (1), .pytest_cache (1); jobs/bench/results/ untouched (mtimes 27 Oct 15:08-15:23)
## Verdict CP4.1 + R10: DO NOT SHIP (0 CRITICAL, 2 IMPORTANT, 7 MINOR). The mechanics hold: the R10 row schema is consistent at every consumer, with no positional reader broken; the bench rows, C5-raw included, are identical at the replicate level; every claimed mutation kills on the shipped tree. Both IMPORTANTs are text or a one-line x-offset. (1) The C5 docstring, _rank_serial and C5.md still say the two variants are near-identical on an ungapped record and that a difference "would indicate a bug", which the new paragraph and the shipped power table (84 of 88 cells) contradict. (2) The calibration power figure still plots the uncentred r_1 as "the dependence a check sees", up to 124% offset at n 20; it is outside the delta's lines, so fixing it now or recording it as owed is Sera's call. Fix (1) and it ships, whichever way (2) goes.

# C6 (SD2) patch, CP5.1 review

SCOPE: CP5.1: sd2_work/review_cp51/delta.diff (working tree vs reviewed_tree = 888de36 + CP4 + CP4.1 + R10; 11 files, 825 diff lines); requirements.md (fix round only; the R8 trend figure is out of scope); claims.md 1-15
COMMIT: 888de36
## Manifest
## Reviewed
- [x] delta check: sd2_work/review_cp51/delta.diff equals a fresh diff of the live tree against reviewed_tree (hunk headers aside); the 11 changed files are exactly the files that differ from reviewed_tree
- [x] src/quebra/analyzers/calibration_summary.py (power_vs_dependence, read_dependence, _pipeline_rows)
- [x] src/quebra/plots/calibration_plot.py
- [x] consumers of ReadDependence / PowerVsDependence / renamed fields (src, jobs, tests, scripts, docs)
- [x] artifacts built from the committed tables (check_calibration job path, rendered)
- [x] jobs/bench/runner.py (+ before/after re-measure through the real run_cell)
- [x] jobs/bench/report.py, jobs/bench/arms.py
- [x] src/quebra/analyzers/checks/c5_rank_autocorr.py, _rank_serial.py
- [x] docs/iid_checks/BENCH.md, C5_rank_autocorr.md (+ recount 84/88, 63/65)
- [x] tests/test_calibration_artifacts.py, tests/test_bench_serial_cells.py (+ mutations)
- [x] claims 1-15
- [x] commands: ruff, format, types, arch, full suite
## Findings (C6 SD2 CP5.1)
(verified) calibration_summary.py:71-78, power_vs_dependence x | arm_e (arms.py:385-425) draws a stationary latent AR(1) with lag-1 correlation rho (_stationary_ar1, burn-in), maps it through Phi and the Weibull quantile, and truncates at tau. Neighbouring durations are a monotone map of a bivariate normal pair with correlation rho, so their Spearman correlation is (6/pi) arcsin(rho/2) exactly for the untruncated sequence, independent of shape (copula) and n. Truncation conditions only on the sum, and the residual between the measured r_1 + 1/N and the identity is the estimator's own small-n bias, which shrinks with n (sd2_work/cp51/rank_corr_check.txt reproduced: -0.120 at n 20, rho 0.5; -0.002 at n 1000). x is shared by every n by construction; the wrong-x alternatives fail on the committed table at n 1000 (x = rho: max |z| 19.3; uncentred measured: 3.91) | -
(verified) calibration_summary.py:276-291, ordering and coverage | rho is sorted once from the Arm E cells; each (row, n) is groupby("rho").mean().sort_index(), so power lists, rho and neighbour_rank_corr share one order. The coverage raise compares exact floats taken from the same column, so it fires only when a (row, n) is missing a rho (a row absent from every replicate of a cell). Built from the committed tables: 6 rho, x = 0.0478..0.4826, 8 ROW_KEYS rows, 5 drawn (all permutation, C5 studentized), every n present | -
MINOR | src/quebra/analyzers/calibration_summary.py:282 (context line, pre-existing) with :117-133 (docstring rewritten by the delta) | every power value is a MEAN over the two Weibull shapes, and PowerVsDependence neither names the pooled factor nor carries the spread, which SizeVsN does (pooled_over, rate_lo/hi_by_row). AGENTS.md section 4, "never call a mean the cell". VERIFIED on power_table.csv, Arm E permutation rows: the 0.75-vs-1.5 spread reaches 0.087 for CvM and 0.084 for C2 (n 1000, rho 0.5: 0.291 vs 0.378, about 4 MC SE of a difference), 0.050 for C6 and 0.039 for C5 (within noise). The x change does not cause it, since x is shape-free, but the rewritten docstring is where a reader looks | add pooled_over and lo/hi per (row, n) as size_vs_n does, or record it as owed to R8 with the trend figure, which uses the same _representative_rows
MINOR | src/quebra/analyzers/calibration_summary.py:51-68 (_pipeline_rows) | the all-failed fallback row that runner.run_cell writes on purpose (check, calibration and variant None) is treated as an unknown check. VERIFIED (sd2_work/review_cp51/rev/probe_pipeline_rows.py): one such row anywhere in size_table.csv, even in an Arm D censored cell that no artifact draws, makes size_vs_n raise "rows outside battery.BENCH_ROW_KEYS: [(nan, nan, '')]", and power_vs_dependence does the same, so jobs/active/check_calibration.py dies on a legitimate table with a message naming the wrong cause. With a second unknown key present the sorted() of mixed float/str tuples raises TypeError instead. read_dependence does not call _pipeline_rows; on an all-failed Arm C row it would raise "more than one segment" (mean_n_segments NaN != 1.0), also the wrong cause. The committed tables carry no such row (check NaN: 0 in both), so nothing breaks today; R8 regenerates them | recognise the fallback row (check NaN, n_used 0) before the unknown-key test and raise a message naming the cell and its fail_reason only if it falls inside what the artifact draws; sort with key=str
(verified) _pipeline_rows on the committed tables and on spelling variants | variant NaN (read_csv) and "" (in memory) both map to "", a float all-NaN variant column works, a capitalised calibration raises naming the key, an empty frame raises KeyError 'check' (not reachable from the job). Built from the committed tables: size_vs_n (8 rows, pooled over shape, 2 cells), power_vs_dependence, read_dependence on both clocks all return; the three figures render through the real plot classes (rev/render.py, rev/{size,power,reads}.png) | -
(verified) calibration_summary.py:323-374, read_dependence | the centring is per cell and on mean_induced_lag1 (C5's uncentred r_1, read from the studentized row: runner.py:296-297), by 1/mean_n_events of the same cell; on one segment E[r_1 | N] = -1/N under exchangeability, and the truncated sample is exchangeable given N. The dedupe is sound on the committed table: in-spec has exactly 4 rows per (n_target, rho) and calendar 9, 54 cells each, and every identifying column (shape, quantised, censoring_target, b, mean_n_events, mean_n_segments) is single-valued within a cell. The SE is std(ddof=1)/sqrt(9) over the 9 cells at a rho, as the field comment says; the new/old ratio is 2.092 on in-spec and 3.162 on calendar for every rho, so claim 3's "about 2x" holds (and the calendar default was 3.2x). Committed in-spec result: centred -0.0016..+0.0026, SE 0.0004..0.0020 | -
MINOR | src/quebra/analyzers/calibration_summary.py:349-350, 366-367, 284-288, 269-272 | four behaviours the delta added have no test; each mutation SURVIVES the six relevant test files (136 passed, sd2_work/review_cp51/mut/log_*.txt): DD (dedupe removed, averaging over rows again), SE (SE divided by sqrt(4n), the old scale), NU (cell-disagreement raise removed), COV (rho-coverage raise removed), ARM (non-E raise removed). DD and SE are the substance of claim 3 (the across-cell SE, a 2.09x change on the drawn error bars) and nothing pins them. AGENTS.md section 7 on sentinels with no test | in test_read_dependence_adds_back..., give the three cells different centred values (e.g. -1/N + {0.01, 0.02, 0.06}) and assert centred_lag1_se == std(ddof=1)/sqrt(3), which kills DD and SE; add one pytest.raises each for a disagreeing cell, a rho missing from one (row, n), and arm != ARM_E
MINOR | src/quebra/analyzers/calibration_summary.py:256-267 (Validity assumptions, Reference) | "no source located for the identity itself": the Gaussian-copula Spearman identity is classical (usually attributed to Pearson 1907, restated in Kruskal 1958, JASA 53). AGENTS.md permits "no source located", so this is compliant; I have not checked an equation locator and do not supply one | optional: locate the equation in Kruskal (1958) and cite it, or leave as is
MINOR | src/quebra/analyzers/calibration_summary.py:457-461 (comment rewritten by the delta) | "The other direction is a row the per-n rule flags at n = 20, which the report, scoring only n >= 35, reads as calibrated" now has no live instance among ROW_KEYS rows: the three ROW_KEYS rows bench_acceptance_at_n flags at n 20 (C1, C2, CvM asymptotic) are all REJECT in the report. VERIFIED on size_table.csv. The live instances of that direction are different in kind: CvM permutation flagged on calendar at n 100 (z 2.56 against a per-n 2.24) and C5 studentized on in-spec at n 700 (z -3.18 against 2.955), both read as size-calibrated by the report (PROMOTE, HOLD), because the per-n threshold covers fewer cells, not because of n >= 35. AGENTS.md section 4: a claim about the artifact must hold in the state it ships | name one of the two live rows and the reason (the per-n Bonferroni covers fewer cells), or drop the sentence
(verified) src/quebra/plots/calibration_plot.py | the power plot draws every n against the one neighbour_rank_corr list (n_values from power_by_row_and_n, as induced_lag1_by_n did before); the read plot reads centred_lag1 and centred_lag1_se and labels both "+1/N"; the dashed raw line is gone with nothing else keyed on "unstudentized" in the file. Neither plotly backend existed before. No test renders either plot (pre-existing); rendered here from the committed tables without error | -
(verified) consumers of the renamed fields and artifacts | grep over src, jobs, tests, scripts, docs for induced_lag1_by_n, induced_lag1_se, induced_lag1 (field), PowerVsDependence, ReadDependence, read_dependence, power_vs_dependence, neighbour_rank_corr, centred_lag1: the only consumers are calibration_plot.py (updated), jobs/active/check_calibration.py (calls the steps by name, no field access) and the new test. No old field name survives anywhere. Known risk for the banner: the step outputs change type shape, so check_calibration's cached identities change and its figures need regenerating (R8 does) | -
MINOR | jobs/bench/grid.py:22-24 (not in the delta, invalidated by it) | "The duration-level lag-1 each rho induces is recorded per cell (mean_induced_lag1), and the power curve is read against that." The power curve is no longer read against it: x is now neighbour_rank_correlation(rho). VERIFIED by reading power_vs_dependence | "... recorded per cell (mean_induced_lag1, C5's uncentred r_1); the power figure is drawn against the generator's own rank correlation, calibration_summary.neighbour_rank_correlation"
IMPORTANT | jobs/bench/report.py:344-363 (dependence_grid_table) and :598-609 with :793 (finding 1), docs/iid_checks/BENCH.md:60 and :185 (all outside the delta) | the class the delta fixed in power_vs_dependence and read_dependence (the uncentred r_1 read as "the dependence induced") is left at its two other consumers, which print it into the committed promotion_report.md. VERIFIED on the committed report: the Arm E grid's induced_lag1 column reads "-0.011 to 0.047" at rho 0.05, whose generator rank correlation is +0.0478, so a reader is told the dependence at that rho may be negative; finding 1 says Arm C's "induced duration-level lag-1 averages -0.013 at rho = 0 and -0.024 at rho = 0.5, with every cell between -0.089 and +0.001", where the centred in-spec values are -0.0016..+0.0026 and -0.089 comes from censored multi-segment size cells (the CP4.1 MINOR's layout effect). The conclusion of finding 1 survives; the numbers that support it do not. AGENTS.md section 4, "fix every site of a class"; the CP4.1 IMPORTANT this delta answers was itself a sibling-site finding. R8 regenerates the report, so the fix is cheapest before R8 runs | in dependence_grid_table print neighbour_rank_correlation(rho) beside (or instead of) the range, or the range of mean_induced_lag1 + 1/mean_n_events; in finding 1 use the centred value on single-segment cells only (read_dependence's numbers). Update BENCH.md:60/185 to say which. Sera's call whether this lands now or as the first step of R8, but it must land before R8 regenerates promotion_report.md
(verified, re-measured through the REAL run_cell) jobs/bench/runner.py, battery_kwargs and _refuse_unservable | sd2_work/review_cp51/rev/runcell_dump.py calls runner.run_cell itself (not a copy of its loop) with run_battery wrapped in the runner's namespace to record every call's kwargs (perm hashed from sizes + indices bytes) and every result's repr(statistic), repr(p_value), repr(r_lag1), plus every returned row. Run from a neutral cwd under reviewed_tree (runner and quebra both from reviewed_tree, checked) and under the live tree: 296 cells (every cell with n <= 60, every 7th larger), all 7 arms incl. 2 EL, 4 replicates each. natural mode: 1,187 battery calls, 9,095 results, 2,279 rows, identical in every byte, 0 differing cells; the 6 natural raises (cells 53 and 149, gamma_hat zero) raise identically. force mode (the first call of every second replicate raises, so the retry runs): 530 forced retries, 528 succeed, every one carries C5-raw in both trees, 0 differing cells. So neither the refusal's move nor the shared kwargs changes any bench row, and the retry now asks for C5-raw exactly as the first call does | -
(verified) jobs/bench/runner.py:258-265, retry comment | true as written: include_tau_checks=False drops C1, C2 and CvM asymptotic and keeps C5, C6 and CvM by permutation (battery.py:126-134, 155-157, 193-198); CvM's statistic divides by gamma_hat (cvm_cramer_von_mises.py:206, 222), and in the natural run both gamma_hat-zero replicates raise again on the retry call with the same message | -
(verified) jobs/bench/runner.py:151-162, 222-226, _refuse_unservable | it draws no randomness and runs before the seed is built, so served cells are unchanged (above); in run_cell it is outside the per-replicate try, so the raise propagates out of run_cell, and out of main() through joblib Parallel, which re-raises a worker's exception. No grid cell trips it (the shipped grid round-trips; the full natural sample includes both EL cells) | -
MINOR | jobs/bench/runner.py:162 (generate's call) | mutation E2, the refusal removed from generate while kept in run_cell, SURVIVES (136 passed): the edited EL test now reaches only run_cell. generate's only other caller is that test's first half (grep over tests, jobs, scripts), so the exposure is a future direct caller being served an EL cell without what it asks for. Claim 6 says the refusal is in both places; only one is guarded by a test | keep one generate assertion in the EL test beside the run_cell one (the pre-edit line), so each site has a test
(note) jobs/bench/runner.py:253-281 | mutation K9R, the retry overriding the shared flag (**{**battery_kwargs, "include_c5_raw": False}), survives; K9 (the flag dropped from the dict) is red. The shared dict makes the CP4.1 B2 drift structurally unlikely rather than tested, which is the alternative the CP4.1 finding offered, so this is not a defect | -
(verified) jobs/bench/report.py:122-126, 170-191, 857-879 | sd2_work/review_cp51/rev/probe_report.py (the CP4.1 probe, retargeted: committed size and power tables, the CP4 smoke serial table with the same four injected rows) builds the report with reviewed_tree's report.py and with the live one. The ONLY difference is the excluded-rows table: it gains clock, shape and segments_target; the all-failed row reads "(every replicate failed)" with rejection_rate "-"; the two apparent A_iid_weibull 355 / quantised / 0.25 duplicates now differ visibly by shape (0.75, 1.5), not clock; size and power rows print "-" for segments_target. _serial_cells' output is unchanged (same decode, factored into _segment_targets). Verdict and rate tables exclude the fallback row by support (n_used 0), so the new label appears nowhere else | -
MINOR | jobs/bench/report.py:857 | segments_target prints as "3.0" and "2.0" in the excluded table (VERIFIED, rev/rep/report_new.md:360-362) because the None of the size and power rows turns the column float, while every serial table above prints the same cells as 3 and 2. Cosmetic, but it is the column added so a reader can match an excluded row to its serial cell | format it as rho and delta are in _serial_cells ("-" or f"{int(v)}") before the table is built
(verified) jobs/bench/arms.py:460-464 | the qualified sentence matches the CP4.1 measurement (var(sqrt(n) r_1) 0.570 at rho -0.4, 1.051 at 0, 2.113 at +0.4) and the grid: every E2 cell has rho in {0.05, 0.1, 0.15, 0.2, 0.3, 0.5}, so "the only sign on the grid" holds | -
(verified, recounted) c5_rank_autocorr.py:13-21, C5_rank_autocorr.md:41-48 | from jobs/bench/results/power_table.csv, pivoting C5's two variants per cell: Arm E 108 cells, 88 differ, raw lower in 84, higher in 4; Arm D 72 cells, 65 differ, raw lower in 63, higher in 2 (every b: 15/15, 18/18, 14/16, 16/16); mean_n_segments is 1.0 in every Arm E and Arm D power cell, so "every one of them a single segment" holds. Both counts exact | -
(verified) _rank_serial.py:30-36, c5_rank_autocorr.py:13-21, C5_rank_autocorr.md:41-48 | IMPORTANT 1's two deletions landed: "near-identical", "would indicate a bug", "very nearly a no-op" and "should agree" are gone from src, docs and jobs (grep). The new wording scopes the 1.04 SD ratio to the scale and names the centre as what separates the variants on one segment, consistently with the uncentred paragraph (c5:33-46) and C5.md:50-56 | -
MINOR | src/quebra/analyzers/checks/c5_rank_autocorr.py:23-28 (context lines next to the edit) | the residue of IMPORTANT 1: the paragraph after the new one still says "Where they GENUINELY separate is short, numerous segments" and "the studentized form is a safeguard for the multi-process regime, not a general-purpose improvement", one sentence after the docstring says centring alone separates them on one segment in 84 of 88 cells. Read in sequence, the single-segment difference is again presented as not genuine. Not as strong as the deleted "would indicate a bug", and "not a general-purpose improvement" is defensible for the sign reason (the raw variant gains power on negative dependence), but the docstring does not give that reason | scope both to the scale: "Where their SCALES separate is short, numerous segments ..." and "So studentizing's scaling is a safeguard for the multi-process regime; its centring matters on any record"
(verified) docs/iid_checks/BENCH.md:40-41 and jobs/bench/runner.py:19-20 (MINOR 3) | true on the committed size table: its censored cells are several segments and keep mean_induced_lag1 (TABLE_COLUMNS unchanged), and the CP4.1 measurement (up to -0.0908 on iid-null cells against -1/N of -0.061) is the layout effect the new sentence names | -
(verified) writing rules, all 387 added lines of the delta | no em or en dash; no "gate" in any form; no surfacing/leverage/delve/data-driven/firstly/excellent/"trust"; no spec or phase identifier (CPn, Rn, SPEC) in a source docstring. power_vs_dependence and read_dependence each carry a Validity assumptions section with the four fields per assumption; the one Dufour-Roy locator matches the one already in c5, C5.md and C6.md | -
(note) no new long .md doc; .claude/sd2-update.md is untracked implementer notes, outside src and docs | -
(verified, re-run) mutations | sd2_work/review_cp51/mut/revmut51.py loads a textually mutated copy of calibration_summary or the runner under its real name before any import (proof_<M>.txt per run; nothing tracked edited); each run is six files (test_calibration_artifacts, test_bench_serial_cells, test_checks_statistics, test_bench_uses_real_carve, test_check_ledger, test_bench_isolation; 136 tests). NONE: 136 passed, exit 0. Claimed kills all hold: X (x = measured lag-1 per rho) 1 red, test_the_power_axis...; S (raw kept) 1 red and S2 (unknown-row raise removed) 1 red, test_no_calibration_artifact...; C (+1/N removed) 1 red, test_read_dependence...; M (multi-segment raise removed) 1 red, same test; E (run_cell's refusal removed) 2 red, both EL parametrisations; K9 (flag dropped from battery_kwargs) 1 red, test_c5_raw_runs_only_when_the_bench_asks_for_it. Extra kills: XRHO (x = rho itself) 1 red; CP (centring by the pooled event count) 1 red, so the per-cell centring is pinned. Survivors: DD, SE, NU, COV, ARM, E2, K9R (findings above) | -
(verified) tests/test_calibration_artifacts.py, markers and oracles | four tests, one tier marker each (unit, statistical, unit, unit; the full suite's marker guard passes, below). Each statistical assertion names its oracle in the first docstring line (the copula identity; the copula identity; Dufour-Roy Eq. (2.2) for -1/N); the raw-row test is a contract test, not a statistical one. A new file is right: no existing file shares its subject (test_tier3_calibration tests null_se and bonferroni_z_crit) | -
(note) tests/test_calibration_artifacts.py:53-68 | the unit test re-types the identity it checks (6/pi arcsin(0.15)), so it pins the WIRING (x is a function of rho, 0 at rho 0, one list for every n) and cannot catch a wrong identity; XRHO dies only because rho != the identity. The independent evidence for the identity is the n = 1000 test against the generator, which is the right split | -
MINOR | tests/test_calibration_artifacts.py:71-91 (n = 1000 agreement) | what it asserts is real but narrower than it reads. (a) It calls no src function: it recomputes the centring inline from the committed table, so it is the Diagnostic for power_vs_dependence's Validity assumption, not a test of read_dependence (mutation C leaves it green; test 4 catches C). Dropping its OWN +1/N turns it red (uncentred max |z| 3.91), and x = rho gives 19.3, so it does discriminate. (b) Its docstring states neither the multiplicity (12 cells, each at |z| <= 3) nor the known second-order bias left after the -1/N centring: the lag-1 autocorrelation's small-sample bias is about -(1 + 3 rho_s)/N to first order (standard for AR(1); Marriott and Pope 1954, locator not checked), i.e. about -1.45 SE at rho 0.5 at n 1000. VERIFIED on the table: adding (1 + 3x)/N instead of 1/N gives max |z| 2.77 at n 1000, 2.39 at 355 and 2.77 at 100, against 2.91, 3.98 and 7.29 with 1/N. The pass at 2.91 of 3 is on these seeds; PLAUSIBLE (back of envelope from the per-cell bias and Bartlett's SD): a re-seeded bench, which any grid reorder causes (BENCH.md limitation 14), goes red about 10% of the time with the identity exactly right | say "12 cells, |z| <= 3 each; the residual -3 rho_s/N bias is about 1.5 SE at rho 0.5" in the docstring, or compare against r_1 + (1 + 3x)/N and cite the bias; either way state the margin
(verified) tests/test_bench_serial_cells.py:40-58 | the EL test now reaches run_cell with n_replicates=1 and still asserts the valid cell is served (n_segments 2) first, so it cannot pass on a generate() that refuses everything; E turns both parametrisations red (with the refusal only inside generate, run_cell catches it per replicate and returns an all-failed row). The pre-edit generate() assertion was replaced rather than kept (E2 survives, finding above) | -
### Claims (sd2_work/review_cp51/claims.md)
- 1 VERIFIED, with a residue: the three sites no longer say near-identical / would indicate a bug; 84/88 and 63/65 recounted exact. Residue: c5:23-28 still says the variants "genuinely separate" only on short segments and that studentizing is "not a general-purpose improvement" (MINOR above).
- 2 VERIFIED: x = neighbour_rank_correlation(rho), one list, mean_induced_lag1 unused, arm != ARM_E raises (the raise is untested: ARM survives).
- 3 VERIFIED: +1/mean_n_events per cell, multi-segment raise, disagreement raise, mean and SE over the 9 cells; SE ratio new/old 2.092 in-spec (3.162 calendar). Only the centring and the multi-segment raise are tested (DD, SE, NU survive).
- 4 VERIFIED: _pipeline_rows keeps ROW_KEYS and raises outside BENCH_ROW_KEYS; dashed raw line gone; _representative_rows simplified. ReadDependence does not pass through _pipeline_rows but carries no row, so nothing reaches it. Caveat: the all-failed fallback row also raises (MINOR above).
- 5 VERIFIED (runner.py:19-20, BENCH.md:40-41).
- 6 VERIFIED that run_cell raises before its first replicate and that the raise leaves main(); both call sites exist. Only the run_cell site is tested (E2 survives).
- 7 VERIFIED (arms.py:460-464; grid rho all positive).
- 8 VERIFIED (probe: "(every replicate failed)").
- 9 VERIFIED (clock, shape, segments_target present; segments_target prints as float, MINOR).
- 10 VERIFIED: one dict; real run_cell before/after identical at the replicate level, natural and forced-retry; retry comment true.
- 11 VERIFIED: X 1, S 1, C 1, M 1, E 2, K9 1 red, re-run on the shipped tree with a NONE baseline in the same batch; n = 1000 max |z| 2.91 recomputed. Not claimed and surviving: DD, SE, NU, COV, ARM, E2, K9R.
- 12 PARTLY FALSE (number), conclusion VERIFIED: 54 cells, max |z| 3.16 (n 35, rho 0.05), mean z +0.03, 26 of 54 positive, recomputed. The Bonferroni bound is 3.31 (calibration_summary.bonferroni_z_crit(54) = 3.312), not 3.34; 3.16 is under either. UNVERIFIABLE recipe: no script under sd2_work prints sd2_work/cp51/arm_c_in_spec.txt (the CP4.1 claim-9 gap again).
- 13 VERIFIED: sd2_work/cp51/rank_corr_check.txt reproduced independently; gap -0.120 at n 20, rho 0.5; the residue is the estimator's -(3 rho_s)/N-type bias (bias-corrected max |z| 2.77 at n 1000).
- 14 VERIFIED: grep over src, jobs (results excluded), scripts, docs, README, AGENTS: every remaining hit is C5's definition, battery's bench path, a statement of exclusion, the studentization comparison, LIMITATIONS 10, BENCH.md's report and multiplicity notes.
- 15 VERIFIED: ruff 0, format 0 (198 files), types 0, arch 0, full suite 0 with 1142 passed and 10 skipped, collect 1152 (commands below).
## Commands run (frozen tree; PYTHONPATH=<worktree>/src:<worktree>, PYTHONPYCACHEPREFIX=sd2_work/.pycache; quebra.__file__ = this worktree's src/quebra/__init__.py; "before" runs used reviewed_tree for both quebra and jobs.bench, run from a neutral cwd and checked by printing both __file__s)
- delta check: delta.diff equals a fresh diff of the live tree against reviewed_tree (hunk headers aside); exactly the 11 listed files differ
- ruff check . exit 0 | ruff format --check . exit 0 (198 files already formatted) | make types exit 0 (11 source files) | make arch exit 0 (1 kept, 0 broken)
- full suite, run once and alone: pytest -p no:cacheprovider, 1142 passed, 10 skipped, exit 0 (rev/full_suite.log); --collect-only 1152
- mutations (mut/run.sh, six files, 136 tests, sequential, never beside the full suite): NONE 0 red; X 1; XRHO 1; S 1; S2 1; C 1; CP 1; M 1; DD 0 (survives); NU 0 (survives); SE 0 (survives); COV 0 (survives); ARM 0 (survives); E 2; E2 0 (survives); K9 1; K9R 0 (survives)
- probes, all exit 0: rev/runcell_dump.py x4 (before/after x natural/force, 296 cells x 4 reps, identical); rev/build_artifacts.py and rev/render.py (all three artifacts and figures from the committed tables); rev/probe_pipeline_rows.py; rev/probe_report.py (reviewed_tree vs live report.py); table analyses for the 84/88 and 63/65 recount, shape spread, Arm C z, the n = 1000 z under variants, bench_acceptance_at_n flags
- git status --porcelain --ignored identical before and after; written outside sd2_work/review_cp51/: this ledger only (sd2_work/ is excluded via .git/info/exclude); jobs/bench/results/ untouched (mtimes 27 Oct 15:08-15:23)
## Verdict CP5.1: SHIP (0 CRITICAL, 1 IMPORTANT, 10 MINOR). The delta does what Sera's answers asked, and every mechanical claim holds: the power x is the generator's own rank correlation, right for Arm E as built; read_dependence centres per cell on the right quantity, with an across-cell SE that is 2.09x the old one, as stated; C5-raw reaches no artifact; the runner refactor changes no bench row at the replicate level through the real run_cell, including the retry path; every claimed mutation kills. The one IMPORTANT is outside the delta's files: report.py and BENCH.md still print the uncentred r_1 as "the dependence induced". It must land before R8 regenerates promotion_report.md, which makes it the first step of R8 if not this commit. Of the MINORs, the cheapest to take now are the c5 docstring residue, the untested SE/dedupe and the 3.34 -> 3.31 correction in the banner.

# C6 (SD2) patch, R8 step 0 review

SCOPE: R8 step 0: sd2_work/review_r8step0/delta.diff (= git diff HEAD, HEAD 1308cbd; 8 files, +157/-41); requirements.md; claims.md 1-12
COMMIT: 1308cbd
## Manifest
## Reviewed
- [x] jobs/bench/report.py (rebuilt from the committed size/power tables + the probe serial table with HEAD and live report.py; rev/report_headbuild.md vs rev/report_new.md differ only at the grid table and finding 1)
- [x] delta check: delta.diff equals git diff HEAD (index lines aside)
- [x] src/quebra/analyzers/calibration_summary.py (_pipeline_rows, PowerVsDependence, power_vs_dependence, the :487-492 comment; probes rev/probe_failed_cell.py)
- [x] src/quebra/plots/calibration_plot.py (rendered all three figures from the committed tables: rev/size.png, rev/power.png, rev/reads.png)
- [x] tests/test_calibration_artifacts.py, tests/test_bench_serial_cells.py (mutations: sd2_work/review_r8step0/mut/revmut8.py + run.sh, six files, 137 tests, sequential; proof_<M>.txt and log_<M>.txt per run)
- [x] src/quebra/analyzers/checks/c5_rank_autocorr.py, jobs/bench/grid.py, docs/iid_checks/BENCH.md, and their unedited siblings (jobs/bench/arms.py, runner.py); writing-rule sweep of all 157 added lines
- [x] claims 1-12, banner correction, commands
## Findings (C6 SD2 R8 step 0)
IMPORTANT | jobs/bench/report.py:739-743 (context lines, not edited) | the prose that introduces the dependence grid still says "with the range of duration-level lag-1 each `rho` induces across its cells", but the column under it is now ONE closed-form generator value per rho (`neighbour_rank_corr`, the same at every n and both shapes). VERIFIED in the rebuild (rev/report_new.md:65-67). R8 regenerates promotion_report.md with a sentence that calls 0.048 a measured range across cells: the third site of the CP5.1 IMPORTANT's class, missed by the commit whose purpose is to fix that class before R8 | "... as the mean over the two Weibull shapes, with the rank correlation of neighbouring durations each `rho` sets in the generator, (6/pi) arcsin(rho/2), the same at every n:"
MINOR | jobs/bench/report.py:428-446 (centred_arm_c_cells) | groupby(clock, n_target, rho).first() takes the first non-null value of EACH column independently and checks nothing about the cell, where its twin read_dependence raises when a cell's rows disagree (AGENTS 4, fix every site of a class). The key also assumes Arm C has one shape, quantised, censoring and b per (clock, n, rho). VERIFIED that the committed tables satisfy both: every identifying and centring column has nunique 1 per cell; 63 cells per clock; the size table holds only Arm C rho 0.0 (never NaN, so the fillna is a no-op) and the power table only rho > 0, so the concat and fill collide on no key. A future Arm C shape or a rho = 0 power cell would be merged silently into one row | before .first(), raise if mean_induced_lag1, mean_n_events, mean_n_segments, shape, quantised or censoring_target has nunique(dropna=False) > 1 within a cell, as read_dependence does
MINOR | jobs/bench/report.py:816-817 | "centred by adding back its exchangeable mean -1/N" reads as adding -1/N, the wrong sign (the code adds +1/mean_n_events); and "averages +0.000 at rho = 0 and -0.003 at rho = 0.5" pools both clocks and every n without naming them (the every-cell range does carry the spread). VERIFIED: rho 0.5 is -0.0046 on calendar, -0.0012 on in-spec | "centred by subtracting its exchangeable mean -1/N (adding 1/N per cell) ... averages, over both clocks and every n, ..."
(verified) claim 1 numbers | rebuilt: grid rho 0.05 reads 0.048; finding 1 reads +0.000 at rho 0, -0.003 at rho 0.5, every cell -0.020..+0.017. 126 cells (63 per clock, 7 rho x 9 n), all mean_n_segments 1.0 and censoring 0, so the raise never fires | -
(verified) claim 1's correction of the CP5.1 ledger | the old minimum -0.089 is a POWER cell, calendar, n 20, rho 0.5, mean_n_segments 1.0, censoring_target 0.0, mean_n_events 14.428 (-1/N = -0.069). The CP5.1 statement that it came from censored multi-segment size cells was wrong; the implementer is right | -
(verified) the centred values are consistent with zero, so finding 1's conclusion is supported | z = centred / (1/sqrt(mean_n_events) / sqrt(n_used of the C5 studentized row)): mean -0.06, SD 1.06, max |z| 3.16 (in-spec n 35 rho 0.05) against bonferroni_z_crit(126) = 3.54; the -0.020 cell is z -2.38. Note: the centring uses -1/mean(N), not mean(-1/N) (Jensen gap about Var(N)/E[N]^3, not measurable from the tables, which carry no Var(N)); nothing in the z pattern needs it | -
(verified) MINOR 10 | probe rev/tables_thin (probe serial table + one injected under-supported serial row): the excluded table prints segments_target "3" with the live report.py and "3.0" with HEAD's; size/power rows print "-". The same row's rho still prints "0.0" where the serial tables print "0" (cosmetic, not asked) | -
IMPORTANT | src/quebra/analyzers/calibration_summary.py:59 (_pipeline_rows, `frame = frame[frame["check"].notna()]`) | the fix drops run_cell's all-failed row UNCONDITIONALLY, where the CP5.1 finding asked to raise when it falls inside what an artifact draws. runner.py:352-353 says that row exists so a failed cell is not "silently omitted"; the drop restores exactly that hole in the calibration figures, and pooled_over is computed over the whole arm frame, so nothing notices a group with fewer cells. VERIFIED (rev/probe_failed_cell.py, fallback row written through CSV as the job reads it): one shape of a DRAWN cell failing in every replicate -> power_vs_dependence returns, the drawn CvM permutation point at n 1000, rho 0.5 moves from the pooled 0.3345 to 0.378 (shape 1.5 alone), lo = hi = 0.378, pooled_over still ['shape']; size_vs_n returns C5 at n 20 as 0.0475 (one shape; pooled 0.0513) with the description still "mean over 2 cells varying shape" and n_cells_pooled 2. Both shapes failing -> power_vs_dependence raises "covers rho [...]" (wrong cause). Before the delta all of these raised, wrongly worded but loudly; the hard rule is errors raised, not swallowed. A fallback row in an undrawn cell (Arm D censored) now passes, as intended. Not reachable on the committed tables (no check-NaN row in either); R8 regenerates them | keep the drop for the unknown-key test, but in size_vs_n and power_vs_dependence raise, after each filter to the drawn region, if any fallback row (check NaN) lies inside it, naming the cell and its fail_reason; and assert every (row, n[, rho]) group holds the full pooled cell count, which also catches a lost shape. read_dependence's "more than one segment" on an all-failed Arm C cell is the same class (unchanged by the delta)
(verified) PowerVsDependence pooled_over and lo/hi on the committed power_table.csv | through the real power_vs_dependence: pooled_over ['shape'] (clock, quantised, censoring_target single-valued on Arm E); all 432 (row, n, rho) groups hold exactly shapes (0.75, 1.5); mean, lo and hi match an independent groupby recompute with 0 mismatches over 8 rows x 9 n; lo <= mean <= hi everywhere. Max hi-lo: CvM perm 0.087, C2 perm 0.084, C6 0.050, C5 0.039 (CP5.1's numbers reproduced); CvM perm n 1000 rho 0.5 lo/hi 0.291/0.378 | -
(verified) _pipeline_rows key=str | an unknown key now raises naming it ("rows outside battery.BENCH_ROW_KEYS: [('zz_new', 'asymptotic', '')]"); the NaN-tuple case that caused the TypeError no longer reaches the sort | -
(verified) calibration_summary.py:487-492 comment | logically right: bench_acceptance_at_n's Bonferroni covers the cells at one n, the report's covers every in-envelope cell of a row, so the per-n z_crit is lower and a cell can exceed it inside the report's band; CP5.1's two live rows (CvM perm calendar n 100, C5 studentized in-spec n 700) are instances. No example, as claim 7 says | -
(note) PowerVsDependence carries no n_cells_pooled or cell_description as SizeVsN does; the count is implicit (2 per point on the committed table). Not asked; the plot footnote is the only place a reader is told (see calibration_plot.py) | -
MINOR | src/quebra/plots/calibration_plot.py:171-179 | the power figure names the pooled factor but does not draw the spread the artifact now carries (power_lo/hi_by_row_and_n is read by nothing; grep), where the size figure fills between rate_lo/hi (:89-94). AGENTS 4 is "name the factor and carry the spread"; on the drawn CvM panel the two shapes differ by up to 0.087 (n 1000, rho 0.5), about 4 MC SE, invisible in the figure. The footnote's "mean over shape" also uses the bare column name where the size figure says "mean over 2 cells varying shape" | draw lo/hi per curve (a thin fill or end-cap bars at each point), or put the largest hi-lo of the drawn rows in the footnote; say "the two Weibull shapes"
(verified) the footnote renders | rendered from the committed tables: "shown: permutation calibration - ...; each point a mean over shape", one line, inside the figure (rev/power.png); no plotly backend before or after | -
(verified, re-run) mutations | plugin loads a textually mutated copy under the real module name before any import (proof file per run; nothing tracked edited); six files (test_calibration_artifacts, test_bench_serial_cells, test_checks_statistics, test_bench_uses_real_carve, test_check_ledger, test_bench_isolation), 137 tests. NONE 137 passed, exit 0. Claimed kills all hold: PF 1, POOL 1, COV 1, ARM 1 (test_power_vs_dependence_names_its_pooling...); DD 1, SE 1, NU 1, C 1, M 1 (test_read_dependence_adds_back...); X 1 (axis test); S 1 (unstudentized-row test); K9 1 (test_c5_raw_runs_only...); E 2 and E2 2 (both EL parametrisations). Extra kills: LO 1, HI 1 (lo/hi set to the mean), SEDDOF 1 (ddof 0). Survivors: KEYSTR, RSIGN, RGRID, RRAISE (findings below) | -
MINOR | jobs/bench/report.py:345-370, 428-446 (no test) | mutations RSIGN (centring sign flipped: finding 1 would print about -2/N), RGRID (the grid column prints rho itself) and RRAISE (the multi-segment raise removed) each SURVIVE the 137 tests; grep finds no test importing jobs.bench.report except the isolation walk. Claim 12 discloses this, but the raise in centred_arm_c_cells is a new sentinel added in answer to a review finding with no test, which AGENTS 7 names as a recurring defect, and the grid column is the substance of the CP5.1 IMPORTANT | a toy-frame unit test beside the read_dependence one: Arm C cells at -1/N + d centre to d, a two-segment cell raises; dependence_grid_table's column equals (6/pi) arcsin(rho/2) and reads 0 at rho 0 (test_calibration_artifacts.py shares the subject, the dependence axes)
MINOR | src/quebra/analyzers/calibration_summary.py:68 (key=str) | mutation KEYSTR (key=str removed) SURVIVES: no test has two unknown keys of mixed type, so claim 4's "no TypeError" is unguarded. Reachable only with a malformed table now that check-NaN rows are dropped first | in the contract test, add two unknown rows, one with calibration NaN, and assert ValueError naming BENCH_ROW_KEYS
MINOR | tests/test_calibration_artifacts.py:79-82 (docstring of test_power_vs_dependence_names_its_pooling_and_refuses_what_it_cannot_draw) | says "another arm, a curve missing a `rho` and an unknown check raise", but the body tests only the first two; the unknown-check raise is asserted in test_no_calibration_artifact_carries_the_unstudentized_row. The test's own docstring overstates what it is evidence for | drop "and an unknown check" or move the stranger assertion here
(note) tests/test_calibration_artifacts.py:97-99 | the all-failed assertion builds a state run_cell never writes (a fallback row beside its own cell's real rows: runner.py:352 writes it only when the cell has no other row) inside the drawn region, and pins the silent drop the IMPORTANT on _pipeline_rows describes. It changes with that fix: a fallback row inside a drawn cell should raise, one outside should pass | -
(verified) tests/test_calibration_artifacts.py:111-117 (n = 1000 docstring, claim 6) | 12 cells (2 shapes x 6 rho at n 1000, n_used 1000 each), recomputed; calls no quebra function (true, it reads the CSV only); the residual bias stated qualitatively and in the right direction (mean z -0.85, low). The margin is not stated: max |z| 2.913 against 3. R8 regenerates power_table.csv, which this test reads, so the CP5.1 estimate (a re-seeded bench goes red about 10% of the time with the identity exact, PLAUSIBLE) belongs in the R8 banner's Known risk | -
(verified) tests/test_calibration_artifacts.py read_dependence test (claim 3) | offsets +0.03, 0, -0.03 at N 20, 50, 400, two rows per cell: centred mean 0 and SE std(ddof=1)/sqrt(3) = 0.03/sqrt(3) exactly; DD (6 row-values: SE 0.0110 not 0.0173), SE, SEDDOF and NU all red. Each statistical assertion names its oracle (Dufour-Roy Eq. (2.2)); the new contract test is marked unit, the oracle rule does not apply | -
(verified) tests/test_bench_serial_cells.py (claim 9) | the EL test asserts the refusal through run_cell AND generate on the same refused cell, after confirming generate serves the valid cell (n_segments 2); E and E2 each turn both parametrisations red | -
MINOR | jobs/bench/arms.py:23-25, :375-377, :404-405 (not in the delta, invalidated by it) | three more sites of the class the delta fixed, in the file grid.py's MINOR 8 fix sits beside. (a) :23 cites promotion_report.md for "the induced duration-level lag-1 never exceeds +0.002 over the whole rho grid"; that is the UNCENTRED maximum (+0.0013, calendar n 500 rho 0.3), and the report R8 regenerates with this delta prints cells up to +0.017 (20 of 126 centred cells exceed +0.002; VERIFIED in rev/report_new.md:204). AGENTS 4: a number in a docstring must come from the artifact it cites in the state it ships. (b) :375-377, Arm C: "the runner reports the INDUCED duration-level lag-1 for every cell, and the power curve is read against that, not against rho", the twin of the grid.py sentence the delta rewrote; Arm C's power is read against read-level rho (report.arm_c_flatness) and no figure draws Arm C power. (c) :404-405 "induces NO positive duration-level dependence at any rho": the evidence was the uncentred r_1 sitting below 0, which is the -1/N offset; centred, 49 of 108 rho > 0 cells are positive, all consistent with zero (max |z| 3.16 against a Bonferroni 3.54 over 126) | (a) "the centred duration-level lag-1 is consistent with zero at every rho on the grid"; cite the report section without a number, or with the regenerated one; (b) "the runner records C5's uncentred r_1 per cell (mean_induced_lag1); the report centres it by 1/N"; (c) "no duration-level dependence distinguishable from zero"
MINOR | jobs/bench/arms.py:399-401 (not in the delta; the identity it describes is now the power axis) | "(6/pi)*arcsin(rho/2), which is within 4% of rho itself over the grid used here": VERIFIED false at 5 of the 6 grid points; the identity is 4.50%, 4.47%, 4.42%, 4.35% and 4.15% below rho at 0.05-0.3, 3.48% at 0.5 (small-rho limit 3/pi = 0.955, 4.5% below) | "within 5%", or "about 4.5% below rho at small rho"
MINOR | docs/iid_checks/BENCH.md:60-61, jobs/bench/grid.py:24-25, jobs/bench/report.py:349-350 | "reads about -1/n with no dependence": the offset is -1/N with N the cell's realised event count, and only on one segment. On Arm C's calendar clock at n 20, N is 14.4, so -1/n = -0.050 against -1/N = -0.069; on the size table's censored cells the column reads down to about -0.09 (BENCH.md:39-41 says so two paragraphs up) | "about -1/N, N the cell's event count, on a single segment"
(verified) c5_rank_autocorr.py:23-28 (claim 2) | "Where their SCALES genuinely separate" and "the scaling half of studentizing is a safeguard for the multi-process regime, and the centring half matters on every record" scope the two statements as the CP5.1 fix proposed, and agree with :13-21 (84 of 88 single-segment Arm E cells differ, the centre being what separates them on one segment). The numbers (1.64, 43%, 21-29%) are unchanged context | -
(verified) grid.py:22-25 (claim 8) and BENCH.md:185-187 | the power curve (power_vs_dependence) and the report's grid column are both neighbour_rank_correlation(rho) = (6/pi) arcsin(rho/2), as stated (calibration_summary.neighbour_rank_correlation; report.py:365) | -
(verified) writing rules, all 157 added lines | no em or en dash; no "gate" in any form; no surfacing/leverage/delve/data-driven/firstly/excellent; no spec, phase or checkpoint identifier in a source docstring. One "taking the claim on trust" (c5:28) is carried over from the pre-existing sentence and is not about a check outcome. No new .md file; .claude/sd2-update.md is untracked implementer notes | -
(verified) banner correction | sd2_work/cp51/arm_c_in_spec.py now exists; re-run pinned from the repo root it reproduces sd2_work/cp51/arm_c_in_spec.txt byte for byte (rev/arm_c_in_spec.rerun.txt): 54 cells, max |z| 3.16, Bonferroni 3.31 | -
### Claims (sd2_work/review_r8step0/claims.md)
- 1 VERIFIED as stated (mechanics and numbers: 0.048; +0.000, -0.003, -0.020..+0.017; BENCH.md :60 and :185 changed; the -0.089 correction of the CP5.1 ledger is right, it is a single-segment calendar n 20 rho 0.5 power cell). INCOMPLETE as a fix of the class: report.py:739-743 still introduces the column as "the range of duration-level lag-1 ... across its cells" (IMPORTANT), and arms.py:23, :375-377, :404 (MINOR).
- 2 VERIFIED.
- 3 VERIFIED: DD, SE, NU, COV, ARM each 1 red; SEDDOF also red.
- 4 VERIFIED as stated; but the drop is unconditional, which silences an all-failed cell inside a drawn region (IMPORTANT), and key=str is untested (KEYSTR survives).
- 5 VERIFIED on the committed table (pooled_over ['shape'], lo/hi exact); the figure names the pooling but does not draw the spread (MINOR).
- 6 VERIFIED (12 cells, residual bias direction, no source function called); margin 2.913 of 3 not stated (note for R8's Known risk).
- 7 VERIFIED.
- 8 VERIFIED.
- 9 VERIFIED: E 2 red, E2 2 red.
- 10 VERIFIED by probe ("3" live, "3.0" at HEAD).
- 11 VERIFIED: every named mutation kills at the claimed count on a 137-test six-file set (the 81-test baseline is a different subset, UNVERIFIABLE as a count, immaterial). "Restored byte for byte" VERIFIED: sd2_work/r8step0/cs.orig and runner.orig equal the live files, and the tree diff equals delta.diff.
- 12 VERIFIED: report.py has no test; RSIGN, RGRID, RRAISE survive (MINOR).
## Commands run (frozen tree; PYTHONPATH=<worktree>/src:<worktree>, PYTHONPYCACHEPREFIX=sd2_work/.pycache; quebra.__file__ = this worktree's src/quebra/__init__.py; jobs.bench.report from this worktree)
- delta check: delta.diff equals git diff HEAD (index lines aside), 8 files
- ruff check --no-cache . exit 0 | ruff format --check --no-cache . exit 0 (198 files) | make types exit 0 (11 source files) | make arch exit 0 (1 kept, 0 broken)
- full suite, once and alone: pytest -p no:cacheprovider -q, 1143 passed, 10 skipped, exit 0 (rev/full_suite.log); --collect-only 1153, exit 0
- mutations (mut/run.sh, six files, 137 tests, sequential, never beside the full suite): NONE 0 red; PF 1; POOL 1; DD 1; SE 1; NU 1; COV 1; ARM 1; X 1; S 1; C 1; M 1; E 2; E2 2; K9 1; LO 1; HI 1; SEDDOF 1; KEYSTR 0 (survives); RSIGN 0 (survives); RGRID 0 (survives); RRAISE 0 (survives)
- probes, all exit 0: rev/build_report.py (live) and rev/build_report_head.py (HEAD's report.py) on sd2_work/r8step0/tables, and on rev/tables_thin; rev/probe_failed_cell.py; rev/render.py (three figures from the committed tables); table analyses (Arm C cells, z, pooled spread, n 1000 margin, identity vs rho); sd2_work/cp51/arm_c_in_spec.py re-run
- git status --porcelain --ignored identical before and after except this ledger ( M .claude/review-findings.md); written outside sd2_work/review_r8step0/: this ledger only; jobs/bench/results/ untouched (mtimes 27 Oct 15:08-15:23); no stray file at the repo root
## Verdict R8 step 0: DO NOT SHIP (0 CRITICAL, 2 IMPORTANT, 9 MINOR). The mechanics are right and every claim holds as written: the grid column is the generator's rank correlation, Arm C is centred per cell with no key collision, the old -0.089 correction is right, pooled_over and lo/hi are exact, and every claimed mutation kills. Two IMPORTANTs, both small: (1) report.py:739-743 still introduces the new column as "the range of duration-level lag-1 ... across its cells", so the report R8 regenerates misdescribes it. This is the class this commit exists to close, and the fix is one sentence. (2) _pipeline_rows now drops the all-failed row unconditionally, so a drawn cell that fails in every replicate silently narrows a pooled mean (VERIFIED by probe), where before the delta it raised. That breaks the hard rule on silent fallbacks, just before R8 regenerates the tables the figures read. Fix both and this ships; of the MINORs, arms.py:23 (a number the regenerated report will contradict) and a toy test for centred_arm_c_cells are the cheapest to take with them.

# R8 trend figure and doc numbers, 4525662

SCOPE: R8 trend figure and doc numbers, 4525662: sd2_work/review_r8/delta.diff (b9839cf..4525662, 10 files: calibration_summary.py, calibration_plot.py, jobs/active/check_calibration.py, tests/test_calibration_artifacts.py, docs/iid_checks/*.md, AGENTS.md; 389 changed lines)
COMMIT: 4525662
## Manifest
## Reviewed
- [x] delta check: sd2_work/review_r8/delta.diff equals git diff b9839cf..4525662 on the 10 files (index lines aside); quebra.__file__ = this worktree's src/quebra/__init__.py
- [x] src/quebra/analyzers/calibration_summary.py (refactor compared against b9839cf on real and 6 perturbed tables: rev/refactor_eq.py; Q2: rev/asym_vs_perm.py, rev/sign_pattern.py; Q3: rev/arm_d_direction.py)
- [x] src/quebra/plots/calibration_plot.py (rendered from the committed tables, both styles: rev/trend_default.png, rev/trend_paper.png, rev/power_*.png; spread rev/trend_spread.py, rev/trend_by_shape.py)
- [x] jobs/active/check_calibration.py (quebra inspect, rev/inspect.txt)
- [x] tests/test_calibration_artifacts.py (13 mutants, rev/mut/summary.txt, logs and proofs per mutant; margin rev/stat_margin.py)
- [x] docs/iid_checks/BENCH.md (runtime.txt; bench_acceptance_at_n recomputed, rev/accept.py)
- [x] docs/iid_checks/C5_rank_autocorr.md (check_agreement re-run + independent pairing, rev/docnums.py)
- [x] docs/iid_checks/C6_exchangeability.md (power by n, agreement, verdict, per-n flag recomputed: rev/docnums.py, rev/accept.py)
- [x] docs/iid_checks/LIMITATIONS.md (C6 row and below-100 claim recomputed, rev/below100.py, rev/docnums.py)
- [x] docs/iid_checks/iid_checks_basics.md (every changed number and its unchanged neighbours checked against the re-run report)
- [x] AGENTS.md (five artifacts counted; sweep for stale four-figure lists)
## Findings (R8 trend figure and doc numbers, 4525662)
MINOR | src/quebra/analyzers/calibration_summary.py:49-52 (REPRESENTATIVE_REASON, now on both power figures) | two defects in the reason a figure prints. (a) "the bench rejected both" is stale and misreads: three asymptotic rows exist (C1, C2, CvM; all REJECT in promotion_report.md:28-35) and the permutation rows PROMOTE, so "both" reads as "both calibrations", the opposite of why permutation is shown (origin db61799, "studentized where a variant exists", when it meant C1 and C2). (b) "within Monte Carlo error" is measured UNPAIRED only: max |z| 2.50 on Arm D's 270 trend cells (Bonferroni 3.74) and 2.21 on Arm E's 324 (3.78), so it holds cell by cell; but the rows share replicates, so the comparison is PAIRED, and the paired SE is not computable from the tables (no discordance count; the paired-z upper bound reaches 6.7). The sign pattern shows a real gap: Arm D b 1.5 shape 0.75, asymptotic below permutation in 15 of 15 nonzero cells (n 20-100, by 0.010-0.043; sign test p 6e-5); Arm E shape 1.5 rho 0.3, above in 26 of 27. Largest gap on a drawn (shape-pooled) point 0.028 (Arm D, C2 n 35 b 1.5) and 0.0225 (Arm E) | "permutation calibration - on the drawn points the asymptotic curves differ from it by at most 0.03, and the bench rejected every asymptotic row on size"
MINOR | src/quebra/analyzers/calibration_summary.py:190-191 (PowerVsTrend docstring) | "event intensity is `b t^(b-1)`" is the TRP's trend function, not its event intensity: for Weibull renewal shapes 0.75 and 1.5 the conditional intensity is z(Lambda(t) - Lambda(T_prev)) * b t^(b-1); arms.py:270 says "trend intensity". The direction wording at :193 is right (see verified) | "whose trend function is `b t^(b-1)`" (as arms.py)
MINOR | src/quebra/analyzers/calibration_summary.py:193-194 | "the curves start from the nominal level there": b = 1 is the middle of the axis (0.7, 0.85, 1, 1.2, 1.5), not a start, and the value there is the row's empirical size, not the nominal level (Arm D b 1 C1 permutation n 20 shape 0.75 reads 0.058). power_vs_trend's own docstring (:428-429) says it right | "so each curve passes through the check's size there"
MINOR | src/quebra/analyzers/calibration_summary.py:420-434 (power_vs_trend) | no Validity assumptions section, where its sibling power_vs_dependence carries one (AGENTS 10). The assumption the figure rests on: the b = 1 point is the same generator as the b != 1 cells (Arm D size cells, 2000 replicates and separate seeds, beside 1000-replicate power cells), and nothing checks it beyond the arm name | add the four-field block: assumption (b = 1 cells are Arm D at b = 1, uncensored, same clock), diagnostic (the arm/censoring filter and _require_full_pooling), consequence (the curve's centre point is drawn from another configuration), reference (jobs/bench/arms.py arm_d)
(verified, Q4) power_vs_dependence after the _power_grid extraction | rev/refactor_eq.py loads b9839cf's module beside the live one: every dataclass field IDENTICAL on the 4525662 power_table and on b9839cf's power_table; and identical exception type and message on 6 perturbed tables (a rho missing for both shapes, one shape cell missing, a failed cell, an unknown check, a clock varying in 5 rows, no Arm E), plus identical output with clock confounded with shape (pooled_over grows). Error texts render identically ("covers rho [...], not [...]"; "power_vs_dependence <row> at n=<n>") | -
(verified, Q3) direction "b < 1 falling event rate, b > 1 rising" (:193) | right for arms.py:270 lambda(t) = b t^(b-1), and measured through the real arm_d (rev/arm_d_direction.py, 400 replicates, n 100): mean duration first vs last quarter 4.07 -> 9.55 at b 0.7 (rate falls), 0.347 -> 0.150 at b 1.5 (rate rises), flat at b 1, both shapes. Agrees with arms.py:276-277 "b < 1 stretches later durations (improving); b > 1 compresses them (degrading)" | -
(verified) power_vs_trend on the committed tables | Arm D: power cells b 0.7/0.85/1.2/1.5 and size cells b 1, all in_spec, quantised False; censoring 0 size cells carry all 9 rows (censored ones lack C2/CvM asymptotic, unread here); pooled_over ['shape'], 2 cells per point. Unfiltered power-side censoring/quantised would raise in _require_full_pooling rather than pool silently | -
IMPORTANT | src/quebra/plots/calibration_plot.py:224-232 (CalibrationTrendPlot draws only power_by_row_and_n) | the new figure draws a mean over the two Weibull shapes and never the spread the artifact carries (power_lo/hi_by_row_and_n read by nothing), and on Arm D the shape decides the answer: max spread per drawn row C1 0.522, C2 0.554, CvM 0.525, C5 0.730, C6 0.742; 129 of 225 drawn points spread > 0.05. C6 at n 355, b 1.5 is 0.216 at shape 0.75 and 0.958 at 1.5, drawn as 0.587, a power neither shape has; C1 at n 100, b 1.2 is 0.195 vs 0.657. This is the AGENTS 4 "never call a mean the cell" class at 10x the power figure's spread (0.087, MINOR in the R8 step 0 review); the footnote "each point a mean over shape" names the factor but a reader cannot see that the factor is larger than the effect. MEASURED: rev/trend_spread.py, rev/trend_by_shape.py on the committed tables | draw the spread as CalibrationSizePlot does (fill_between lo/hi, :74-94), or draw one line style per shape (the pooled factor has two levels); apply the same to CalibrationPowerPlot to close the class
MINOR | src/quebra/plots/calibration_plot.py:243 (footnote) and calibration_summary.py PowerVsTrend | the figure never says what configuration it draws (in-spec clock, uncensored) or that its b = 1 points are size cells from another table (2000 replicates vs 1000), and it drops Arm D's censored b = 1 cells (censoring 0.03, 0.25) without saying so (FIGURE_STANDARD "Dropped data" binds figures you add). SizeVsN carries a cell_description for exactly this | add cell_description to PowerVsTrend ("D_trp_power_law, in_spec, censoring 0; b = 1 from the size table") and print it in the footnote
(verified) rendering | rev/render.py renders both power figures from the committed tables in "default" and "paper" styles, saved as the targets save them (bbox_inches tight): no warning; one-line footnote inside the figure; "no trend" label clear of the curves; legend outside the last panel, not clipped; title one line and names the object (FIGURE_STANDARD); no plotly backend before or after. rows_drawn: the five permutation rows (C5 studentized), b [0.7, 0.85, 1.0, 1.2, 1.5], pooled_over ['shape'] | -
(verified) _power_footnote extraction | the power figure's footnote text is byte-identical to before (selection_reason unchanged, see Q4); only its suptitle changed ("induced duration dependence" -> "serial dependence") | -
(note) calibration_plot.py:204 | `rows = list(result.rows_drawn) or sorted(...)` copies the power figure's fallback: an empty rows_drawn would draw every row under a footnote that still says "shown: permutation calibration". Unreachable today (rows_drawn is empty only when power is, which raises first) | -
(verified) jobs/active/check_calibration.py:60 wiring | `quebra inspect` (rev/inspect.txt, exit 0, nothing written under output/): "power_vs_trend <- [load_df_2, load_df]", i.e. (power_table, size_table) in signature order; five job.figure calls, as the docstring now says; both tables already Datasets, so the new step adds no input outside provenance | -
(verified, Q5) mutations | harness rev/mut/revmut9.py + run.sh: a -p plugin loads a textually mutated calibration_summary under its real name before any import (each edit asserted to match exactly once; proof_<M>.txt per run), or for TABLE patches the test module's RESULTS at setup to rev/mut/tables_bad (one value changed, verified). Test set: test_calibration_artifacts, test_tier3_calibration, test_bench_isolation (42 tests), sequential, never beside the full suite. NONE 42 passed exit 0. Red (test that went red): T_CENS b = 1 censoring filter removed 2 (unit, statistical); T_NOCONCAT b = 1 cells never joined 2 (unit, statistical); T_NORAISE no-null raise removed 1 (unit); T_TRENDARM power side not filtered to Arm D 1 (statistical); T_ALLROWS rows_drawn = every row 1 (statistical: an asymptotic row breaks the band); G_MEANMIN point = min not mean 2 (power_vs_dependence contract, statistical); G_POOL pooled_over emptied 5; G_COV coverage check removed 1 and G_FULL _require_full_pooling removed 1 (both by the pre-existing power_vs_dependence contract test, which now guards the shared _power_grid); TABLE (Arm D b 1 C6 permutation n 100 shape 0.75 0.049 -> 0.085) 1 (statistical). Survivors: T_B1 and T_NULLARM, see MINOR | -
MINOR | tests/test_calibration_artifacts.py:270-305 (unit test) | the `b == 1.0` and `arm == ARM_D` filters on the size side are unguarded: mutants T_B1 (b filter removed) and T_NULLARM (arm filter removed) SURVIVE all 42 tests, because the synthetic size table holds only Arm D b = 1 cells and the committed one has Arm D only at b = 1 and NaN b elsewhere. Equivalent on every table today; they stop being so the day the bench adds an Arm D size cell at b != 1 or another arm carries b | add one Arm A b = 1 row and one Arm D b = 0.7 size row to the synthetic size table, and assert the b = 1 point is unchanged
(verified) the statistical test's numbers | independent recompute from the CSV alone (rev/stat_margin.py): 45 points (5 drawn rows x 9 n), z_crit 3.261, max |z| 2.974 (C6 n 500), mean z -0.29, SD 1.13; margin 0.29 SE, deterministic on the committed table, so no flake but a narrow pass worth a Known risk line when the bench is re-run. Sensitivity: at n_used 2000 the pooled point resolves a mean deviation of 3.261 x 0.00345 = 0.011, so a single shape cell must be about 0.022 off alpha to fail it alone. Oracle named on the docstring's first line; marker statistical; the SE formula for the mean of two independent cells is right | -
(verified) the unit test's "a censored size cell (rate 0.9) must not enter it" | true. Under T_CENS what turns it red is _require_full_pooling raising; with that guard also removed (double mutant T_CENS_G_FULL) the unit test is still red, on `pooled_over == ["shape"]` (it reads [shape, censoring_target]), so the test holds without leaning on the pooling guard | -
(verified, Q1) docs/iid_checks/BENCH.md:54 "591 cells and 816,000 replicates" | jobs/bench/results/runtime.txt at 4525662: cells=591, replicates=816000 | -
(verified, Q1) docs/iid_checks/BENCH.md:150-156 finding 8 | bench_acceptance_at_n on the committed size_table.csv (rev/accept.py, rev/acceptance.csv): 162 groups (9 rows x 18), 24 flagged, of them exactly four permutation rows: C5 studentized in_spec n 700 at 0.0345 and C6 in_spec n 700 at 0.0355, whose worst cell is the same (Arm A, shape 0.75, continuous, uncensored, n_used 2000), both below alpha; C5 unstudentized in_spec n 20 (0.0355, Arm A q True c 0.03); CvM permutation calendar n 100 (0.0625, Arm B). "Same rank autocorrelations off the same permutations": C6 sums squared studentized rank autocorrelations at lags 1-5 and the battery runs the permutation checks on one shared set (AGENTS 6). "Demotes a ledger non-rejection to underpowered": check_ledger._verdict:313-319, bench_accepted False appends to unmet and returns VERDICT_UNDERPOWERED | -
MINOR | docs/iid_checks/BENCH.md:148-149 (context line beside the new count) | "over 162 groups, so about eight false flags are expected even for exact tests" is now read against "Four exact permutation rows are flagged", but only the 90 permutation groups (5 rows x 18) are exact, so the like-for-like expectation is 4.5 (90 x 0.05), not eight; 162 is still the right group count (recomputed). The reader is invited to compare 4 with 8 and conclude the exact rows run below chance | "about eight false flags are expected over all 162, about four and a half over the 90 permutation groups"
(note) docs/iid_checks/BENCH.md:156 | "A flagged group demotes a ledger non-rejection" is true for ledger rows; C5 unstudentized, one of the four, is in BENCH_ROW_KEYS only and has no ledger row, so its flag demotes nothing | -
(verified, Q1) docs/iid_checks/C5_rank_autocorr.md:81-83 "288 power cells ... 0.019 on average ... at most 0.172 (promotion_report.md, criterion 4)" | promotion_report.md criterion 4 row c5 studentized: 288 / 0.019 / 0.172; report.check_agreement re-run on the committed tables gives the same; independent pairing on the full cell key: 288 shared, mean 0.0186, max 0.1720 (rev/docnums.py) | -
MINOR | docs/iid_checks/C5_rank_autocorr.md:81-83 and docs/iid_checks/C6_exchangeability.md:62-63 | "over the 288 power cells they share" does not say which: check_agreement pairs over EVERY power cell, Arm C 108 + Arm D 72 + Arm E 108, and the 0.172 maximum is an Arm D trend cell (in_spec n 1000 shape 1.5 b 0.85), not a dependence cell. In C6.md the sentence closes a paragraph about Arm E power, so it reads as the dependence grid; on Arm E alone the figures are 0.017 mean and 0.077 max. AGENTS 4: state the aggregation. (promotion_report.md:177 says "across the shared power cells of the arm both are directed at", which the code does not do; outside this slice) | "over the 288 power cells of Arms C, D and E they share (largest on an Arm D trend cell; on Arm E alone 0.017 and 0.077)"
(verified, Q1) docs/iid_checks/C6_exchangeability.md:57-61 power by n | Arm E, mean over rho and shape, from power_table.csv: C6 0.220 at n 50, 0.343 at 100, 0.605 at 355, 0.791 at 1000 (22%, 34%, 60%, 79% as written); C5 studentized 0.236, 0.357, 0.626, 0.804, above C6 at all nine n by 0.008-0.022 ("slightly below at every n" holds) | -
(verified, Q1) docs/iid_checks/C6_exchangeability.md:68-74 verdict and size | report.verdicts on the committed tables: C6 HOLD, worst 0.0355 z -2.98 of 3.58, mean power 0.343 at n 100; promotion_report.md size-per-row: 144 cells, worst A_iid_weibull/in_spec n=700 q=False c=0.0. Per event count: bench_acceptance_at_n flags C6 in_spec n 700, z -2.975 against 2.955 over 16 cells (prints as -2.98 and 2.96). "BENCH.md, section 5, item 8" exists and is the finding-8 paragraph | -
(note) docs/iid_checks/C6_exchangeability.md:69 "the 144 null cells inside the envelope" | report.scored_cells also drops n = 20 (min_n = 35), so "inside the envelope" here means censoring <= 0.03 AND n >= 35, while bench_acceptance_at_n (the per-event-count reading in the next sentence) includes n = 20. Harmless for C6 (no flag at n 20); it is why BENCH.md's finding 8 can list C5 unstudentized at n = 20 while promotion_report.md gives that row's worst |z| as 2.56 | -
(verified, Q1) docs/iid_checks/LIMITATIONS.md:10-19 | the new C6 row (0.097 0.166 0.220 0.295 0.343 0.605 0.674 0.732 0.791) equals the Arm E mean over rho and shape from power_table.csv at all nine n, as does the unchanged C5 row; "below n = 100 neither exceeds 0.5 except at rho = 0.5" holds for both checks on the shape mean AND on every single-shape cell (rev/below100.py: the only cells above 0.5 below n 100 are rho 0.5 at n 35, 50, 75); "0.48" is neighbour_rank_correlation(0.5) = 0.483 | -
(verified, Q1) docs/iid_checks/iid_checks_basics.md:41-43, :93, :102, :108-109 | 288 / 0.019 / 0.172 (criterion 4 and re-run); 591 cells, 816,000 replicates (runtime.txt); C6 HOLD, worst z -2.98 over 144 cells, power 0.34 at n 100 (0.343); C6 22/34/60/79% (0.220, 0.343, 0.605, 0.791). The unchanged neighbours still match the regenerated report: C1/C2 permutation 128 cells, 0.60/0.62; C1/C2 asymptotic z 31.7 and 5.4; CvM 144 cells z 2.56, asymptotic 4.41 over 86; C5 0.34-0.36; C5 24/36/63/80% | -
(note) docs/iid_checks/iid_checks_basics.md:41-43 | the 288-cell agreement sentence repeats the unnamed aggregation of the C5/C6 MINOR above (Arms C, D, E; max on an Arm D trend cell); same fix | -
(note) docs/iid_checks/iid_checks_basics.md:109 "under the same condition" | ambiguous between "the same averaging" and "exceeds 50% below n = 100 only at the strongest dependence"; both readings are true on the committed table (rev/below100.py), so it is wording only | "with the same averaging, and the same pattern below n = 100"
(verified) writing rules, all 317 added lines of the slice | no em or en dash; no "gate" in any form; none of surfacing / brings into view / data-driven / delve / leverage / firstly / excellent / trust; no spec, phase or checkpoint identifier; no new .md file; no unsuffixed time or frequency name added | -
MINOR | docs/iid_checks/BENCH.md:91-93 (unchanged lines in a file this commit edits) | "draws calibration_summary.size_vs_n, power_vs_dependence, read_dependence and validation_curve" now omits power_vs_trend, the fifth step the job wires; the commit updated "four" to "five" in calibration_summary, calibration_plot, check_calibration and AGENTS.md but not this list (AGENTS 4, every site of a class; AGENTS 8, a doc claim must match the code). grep finds no other stale "four" for the calibration figures | add `power_vs_trend` to the list
(verified) AGENTS.md:244 "five typed artifacts" | calibration_summary.py defines five public artifact dataclasses (SizeVsN, PowerVsDependence, PowerVsTrend, ValidationCurve, ReadDependence) plus the private _PowerGrid; calibration_plot.py five plot classes; the job five figures | -
(verified) promotion_report.md in its shipped state | report.build() re-run on the committed tables with runtime.txt as main() passes it, written to rev/promotion_report_rebuilt.md: byte-identical to jobs/bench/results/promotion_report.md (diff exit 0), so every doc number cited to the report is also a number the committed tables produce | -
## Commands run (frozen tree at 4525662; PYTHONPATH=<worktree>/src:<worktree>, PYTHONPYCACHEPREFIX=sd2_work/.pycache; quebra.__file__ = this worktree's src/quebra/__init__.py; jobs.bench.* from this worktree)
- delta check: delta.diff equals git diff b9839cf..4525662 on the 10 files, exit 0
- ruff check --no-cache . exit 0 | ruff format --check --no-cache . exit 0 (199 files) | make types exit 0 (11 source files) | make arch exit 0 (1 kept, 0 broken)
- full suite once, alone: pytest -p no:cacheprovider -q exit 0, 1148 passed, 10 skipped (rev/full_suite.log)
- mutations, sequential (rev/mut/run.sh, 42 tests in 3 files): NONE 0 red; T_CENS 2; T_NOCONCAT 2; T_NORAISE 1; T_TRENDARM 1; T_ALLROWS 1; G_MEANMIN 2; G_POOL 5; G_COV 1; G_FULL 1; TABLE 1; T_CENS_G_FULL 2; T_B1 0 (survives); T_NULLARM 0 (survives)
- probes, all exit 0: rev/refactor_eq.py, rev/asym_vs_perm.py, rev/sign_pattern.py, rev/arm_d_direction.py, rev/render.py, rev/trend_spread.py, rev/trend_by_shape.py, rev/stat_margin.py, rev/accept.py, rev/docnums.py, rev/below100.py, report.build() rebuild; quebra inspect jobs/active/check_calibration.py exit 0
- git status --porcelain --ignored: only this ledger modified; nothing written outside sd2_work/review_r8/rev/ and this ledger; jobs/bench/results/ untouched; output/ unchanged
## Verdict R8 trend figure: DO NOT SHIP (0 CRITICAL, 1 IMPORTANT, 9 MINOR). The mechanics hold: power_vs_dependence is identical to b9839cf on real and perturbed tables, error messages included; the trend wiring, the b direction and every doc number recompute from the committed artifacts, and the committed promotion_report.md is a byte-exact rebuild of the committed tables; both new tests go red under the mutations they guard (two equivalent-on-data filters excepted). The IMPORTANT is the figure itself: it draws the mean over two Weibull shapes whose powers differ by up to 0.74 (C6 n 355 b 1.5: 0.216 vs 0.958), the AGENTS 4 class at ten times the spread the power figure was flagged for. Drawing the lo/hi the artifact already carries, or one line style per shape, fixes it. Of the MINORs, the REPRESENTATIVE_REASON wording ("rejected both"), the BENCH.md:91-93 figure list and naming the arms behind "288 power cells" are the cheapest to take with it.

# Patch A: Q1.1 band maximum (uncommitted on 24a46d8)

SCOPE: src/quebra/analyzers/tellability.py src/quebra/analyzers/spectrum.py docs/TELLABILITY_SEMANTICS.md AGENTS.md tests/test_tellability_metric.py tests/test_tellability_result.py tests/test_spectrum.py
COMMIT: 24a46d8 (change uncommitted; delta.diff == git diff, cmp exit 0)
## Manifest
## Reviewed
- [x] src/quebra/analyzers/spectrum.py  (+16/-9)
- [x] src/quebra/analyzers/tellability.py  (+68/-11)
- [x] docs/TELLABILITY_SEMANTICS.md  (+13/-3)
- [x] AGENTS.md  (+3/-3)
- [x] tests/test_spectrum.py  (+15/-1)
- [x] tests/test_tellability_metric.py  (+79/-6)
- [x] tests/test_tellability_result.py  (+2/-2)
## Findings (Patch A)
(verified) src/quebra/analyzers/spectrum.py:61-70 | chunking is bit-identical to the old one-block formula `abs(exp(-2j*pi*outer(f,t)).sum(1))**2/N**2`: 930 layouts (N in 1..101, 1..97 frequencies, every block size 1..F+2, short last block included), rows=1 when _CHUNK_ELEMENTS < N, and the default 3-block path on 3000x3000, all np.array_equal (rev/p1_grid_chunk.py) | -
(verified) src/quebra/analyzers/spectrum.py:63-64 "about 1.5 GB" | largest record 6D2S/220423_6D2S_qubit1: N 4827, 20130 grid points, one complex128 block 1.55 GB (rev/p3_real.py) | -
MINOR | src/quebra/analyzers/spectrum.py:46-47 | "Whether anything above the band edge can fold onto the band is read over the whole first harmonic's cell" overclaims: a W peak above 3 f_s/2 (any higher harmonic, or a schedule rate outside the cell, which window_band_maximum's own second assumption names) also folds onto the band, and the cell maximum does not see it | "Whether the first harmonic can fold onto the band is read over its whole cell"
IMPORTANT | src/quebra/analyzers/tellability.py:240-241 | f_s is a grid point only when n_steps = ceil(4 T f_s) is even, and on a lattice the parity is set by float rounding: 4*T*f_s = 1196.0000000000018 for 300 reads 0.1 s apart, so ceil gives 1197. 15 of 112 perfect lattices probed (spacings 0.1 s to 3600 s, N 2 to 2001) read 0.942-0.950 at f*/f_s != 1, where the old row read exactly 1 (rev/p1_grid_chunk.py, rev/p4_tests_geom.py). The new row can read BELOW W(f_s), the very ordinate the docstring calls "not enough"; the half-step test's own record shows it: W(f_s) 0.9960, band max 0.9472. test_the_window_band_maximum_is_one_at_f_s_on_a_lattice passes only because 300 x 24 s happens to round even. No real record reads below W(f_s) today (53 of 53, rev/p3_real.py), but 050423_6D2S_qubit1 peaks exactly at f_s (0.4298) on an even grid | n_steps = 2 * math.ceil(2.0 * T * first_alias_hz): still a step <= 1/(4T), and f_s is always the middle point, so band max >= W(f_s) by construction and is 1.000000 on all three failing lattices (probed). Add a lattice whose spacing rounds odd (0.1 s x 300) to the lattice test, and rebuild the half-step test around a peak placed off f_s
IMPORTANT | src/quebra/analyzers/tellability.py:226-230 | the second assumption's Diagnostic, `window_band_peak_over_fs`, cannot detect its violation. Probe: reads on a 10 s lattice with two thirds of spacings 20 s, so median 20 s and the schedule's rate is 2 f_s: W(1/d) = 1.0000, the band max reads 0.0136 (near the (ln M)/N noise scale) at f*/f_s = 0.9973, an ordinary-looking location (rev/p4_tests_geom.py). The row then reads "nothing folds" while a perfect alias sits at 2 f_s = 4 f_e | "Diagnostic: none shipped; W at 2 f_s, or the maximum over the next cell, would show it. f*/f_s does not: the cell's lesser maximum can sit anywhere in it"
MINOR | src/quebra/analyzers/tellability.py:211-213 | "Every frequency lies within half a step, 1/(8T), of a grid point" is false within half a step of either cell end: the ends are excluded, so a frequency just above f_s/2 or just below 3 f_s/2 is up to a full step, 1/(4T), from the nearest point (step*T = 0.2500 probed) | "Every frequency in the cell, except within half a step of its ends, lies within ..."
MINOR | src/quebra/analyzers/tellability.py:212-213 and :232-235 | sinc^2(1/8) = 0.9496 is the large-N limit for reads spread evenly over the span; a lattice of N reads with T = (N-1) Delta keeps 0.9387 at N 10, 0.9011 at N 3, 0.8536 at N 2 (probed). The third assumption's cos^2(pi/8) = 0.8536 (verified) is more than an example: it is the floor for ANY fully coherent peak, since at 1/(8T) off every read's phase lies in an arc of pi/4, so the mean phasor has modulus >= cos(pi/8). Bursts of width 1% and 5% of T read 0.856 and 0.867 (probed) | state 0.95 as the many-evenly-spread-reads value and cos^2(pi/8) = 0.85 as the floor for a fully coherent peak, reached by two bursts at the ends
(verified) src/quebra/analyzers/tellability.py:215-218 | the (ln M)/N scale: uniform read times, 20 seeds each, mean band max / ((ln M)/N) = 1.26, 1.23, 1.18, 1.19 at N 50, 200, 1000, 3000 (1.09-1.11 against (ln M + 0.577)/N); M = f_s T is the band width over the independent spacing 1/T. A fair reading scale, and it is framed as one (rev/p2_lnM.py) | -
(verified) src/quebra/analyzers/tellability.py:238-241 | grid: step f_s/ceil(4 T f_s) <= 1/(4T); ends excluded; T >= the median positive spacing (every positive spacing <= span), so 4 T f_s >= 4 and the grid has >= 3 points; 2 reads give n_steps 4, points 0.75/1.00/1.25 f_s, band max 1 at f_s for spacings 0.1, 1, 7.3, 24 s; t far from 0 (t0 = 3.3e7 s) reads the same as t0 = 0 (rev/p1_grid_chunk.py) | -
(verified) src/quebra/analyzers/tellability.py:1292 | moving the computation into `_attempt` changes no surfaced exception: `_attempt` catches only NotComputed; band_edges raises it before window_band_maximum is evaluated (tuple elements left to right), as before when the window call sat behind `if edges is not None`; window_band_maximum raises only ValueError/numpy errors, which propagate as before. The status is still set only after both succeed. Old pickles lack window_band_max, so StaleArtifactGuard.__setstate__ refuses them | -
(verified) real records, read-only (rev/p3_real.py) | 34 6D2S + 19 2x2 records: band max >= W(f_s) on every one; W(f_s) is 0.0000-0.43 while the band max is 0.007-0.88 at f*/f_s 0.94-1.05, so Q1.1 changes by up to three orders of magnitude on real data. Per-record cost is of order a second on the largest record (load-dependent), so "small for the record sizes here" holds | -
(note) src/quebra/analyzers/tellability.py:242-247 | argmax takes the first of tied maxima; on a perfect lattice W(f_s) and W near the excluded harmonics never tie inside the cell, so the location is unambiguous there | -
MINOR | docs/TELLABILITY_SEMANTICS.md:135-136 | repeats the docstring's "keeps about 95% of a 1/T-wide peak" with neither caveat: no floor for few or bursty reads (cos^2(pi/8) = 0.85, the docstring's third assumption) and, until the parity IMPORTANT is fixed, a perfect lattice can read 0.95 rather than 1. A reader of the doc alone gets the best case only; AGENTS 4, every site of a class | carry the 0.85 floor beside the 95%, and change this sentence in the same commit as the grid fix
MINOR | docs/TELLABILITY_SEMANTICS.md:137-139 | `(ln M)/N`, `N` and `M = f_s·T` sit outside backticks, the only formulas in the paragraph that do, and `·` is the only non-ASCII character the patch adds | `` `(ln M) / N` ``, `` `N` ``, `` `M = f_s T` ``
(verified) docs/TELLABILITY_SEMANTICS.md:101-104 | "about 1/N at any one frequency for scattered times" is right (E W = 1/N for independent uniform phases, and the module's own scattered-times test); the folding sentence is the converse-free reading. :35 "fewer than <n> reads ... window (2)" still matches band_edges' NotComputed; :79-80 "the spectral window ... count every read" still matches window_band_maximum (obs.t_rel_s, no collapse) | -
(verified) leftovers | grep over src, tests, docs, jobs, scripts, README, AGENTS: no `window_at_first_alias`, no `spectral_window_at_first_alias`, no sentence still saying W(f_s) near 0 means nothing folds; the only W(f_s) mentions left are the new "W(f_s) alone is not enough" ones | -
(verified) AGENTS.md:292-294 | seven arms A, B, C, D, E, E2, EL match jobs/bench/arms.py:92 `ARMS`; jobs/bench/results/serial_table.csv exists; no other "five arms" left in AGENTS.md or docs | -
MINOR | AGENTS.md:294 | 129 characters, the only line in the file over 100; the paragraph wraps near 95 | rewrap the paragraph
IMPORTANT | tests/test_spectrum.py:63-73 with src/quebra/analyzers/spectrum.py:66 | the block test cannot fail for a chunk loop that stops early. Mutant BLOCK_FIRST_ONLY (`range(0, min(f.size, rows), rows)`, only the first block written) SURVIVES the three touched files: `power = np.empty(f.size)` is handed back the buffer the test's own one-block call freed a moment earlier, which already holds the right unnormalised values, so the output is bit-equal. Verified directly: equal right after the one-block call, wrong from index 7 once the small-buffer cache is churned (rev/p7_empty_reuse.py). The two dropped-block mutants that do go red go red only because they also break the one-block path. The test is the stated evidence for "computed in blocks, it equals the same window computed in one block" (AGENTS 4, a positive control must fail when the thing it guards is broken) | `power = np.full(f.size, np.nan)` in spectral_window, so an unwritten block is NaN; and in the test compare against an inline one-block formula computed BEFORE the chunked call, not after
MINOR | tests/test_spectrum.py:63-73 | the `rows = max(1, ...)` floor is untested: mutant ROWS_NO_MAX1 survives (rows 0 makes `range` raise only when N > _CHUNK_ELEMENTS) | add a case with _CHUNK_ELEMENTS below t.size, so rows is 1
(verified) tests/test_spectrum.py:65-66 | the docstring's layout is right: (7*25+3)//25 = 7 rows, 40 = 5 x 7 + 5 | -
IMPORTANT | tests/test_tellability_metric.py:70-138 | the cell's lower edge, the one the docstring calls what "keeps the main lobe, W(0) = 1, out", is unguarded: mutant LOW_EDGE_TO_0+ (grid from f_s/n_steps instead of f_s/2) SURVIVES all three files. Its consequence is large: 210423_6D2S_qubit1 would ship 0.809 at f*/f_s = 0.0003 instead of 0.478 at 0.976, and 1000 uniform read times 0.815 instead of 0.0096 (rev/p6_lowedge.py). Every new test uses a lattice-like record whose alias peak (0.88-1.0) outranks the main-lobe sample at 1/(4T) (about 0.81), so none can see it | add a structureless-times test: uniform reads, band max below a small multiple of (ln M)/N (measured 1.18-1.26 x) and f*/f_s inside (1/2, 3/2). It also turns the (ln M)/N sentence into a tested claim
MINOR | tests/test_tellability_metric.py:70-138 | the other edges are unpinned too: LOW_EDGE_TO_fs/4, HIGH_EDGE_TO_5fs/2 and ENDS_INCLUDED survive (near-equivalent on these records, so lower stakes). The rate-outside-cell schedule (10 s lattice, median 20 s: W(2 f_s) = 1, band max 0.0136) would pin the upper edge and document the second assumption's limit in one test | add it, asserting the band max stays near the noise scale
MINOR | tests/test_tellability_metric.py:71-77 | test_the_window_band_maximum_is_one_at_f_s_on_a_lattice passes only because 300 x 24 s rounds to an even n_steps; 300 x 0.1 s reads 0.9494 at 0.999582 (see the tellability.py:240-241 IMPORTANT). The analytic oracle holds for this one spacing, not for "a perfect lattice" | parametrize over a spacing that rounds odd; it goes red today and green after the even-grid fix
MINOR | tests/test_tellability_metric.py:126-138 | the half-step test re-derives n_steps with its own copy of the grid formula, so `assert n_steps % 2 == 1` describes the test's grid, not the code's: change the code's grid (the even fix, or any other) and the precondition still passes while f_s is no longer half a step off, and the docstring's "the worst offset" silently stops being true. The 0.94 tolerance itself is justified: measured 0.9472 (sinc^2 at N 502 less the extra read's share), and it pins density, 1/(2T) and 1/(3T) grids read 0.8095 and 0.9099 (both mutants red here) | read the offset from the code's output (|f*/f_s - 1| against the reported step), or build the off-grid peak independently of f_s
MINOR | tests/test_tellability_metric.py:106-123 | "within two grid steps" is twice to four times the needed slack, because the reference is median/mean over the truncated 499 steps (0.980739) while the docstring names 24.48 s, the full-period mean, where the continuous peak lies (0.980400 against 24/24.48 = 0.980392). Measured: grid peak 0.78 steps from the test's reference, 0.09 steps from 24/24.48. With 2 steps, an off-by-one location mutant passes this test (caught only by the lattice and repeated-times tests) | reference 24.0 / 24.48 and tolerance half a step plus a small slack
MINOR | tests/test_tellability_metric.py:80-103 | the jittered test's window is loose: band mean 0.8922 sits 18.5 SE above its lower bound and 9.9 SE below its upper, so 1/(2T) and 1/(3T) grids pass it (it catches the old row, a coarse 1/T grid and an upper-half-only cell). Both bounds are sound: W at 1/Delta is 1.34 SE from the analytic value; the continuous maximum exceeds W(1/Delta) by only 0.00017 (0.08 SE), so the upper bound holds; seeds are fixed, so no flake. Across seeds 1/Delta falls uniformly within half a step of the grid (mean |offset| 0.25 steps), so the expected value is expected x mean sinc^2 over the offset, 0.983, i.e. 0.891 against 0.892 measured (rev/p5_jitter.py) | if this test is to pin the step, assert that value within 3 SE; otherwise say in the docstring that the half-step test pins the step and this one pins the location of the cell
(verified) tests/test_tellability_metric.py:256-275 | the repeated-times test goes red under COLLAPSE_REPEATS and under PEAK_LOC_OFF_BY_1/2; rel=1e-9 is right (same frequency up to an ulp of f*/f_s x f_s) | -
(verified) reviewer mutations on a copy (rev/mut.py, rev/mut.log): 20 mutants. RED: FS_ONLY (old row), GRID 1/T, 1/(2T), 1/(3T), UPPER_HALF_ONLY, PEAK_LOC_OFF_BY_1, _BY_2, COLLAPSE_REPEATS, WIRE_SWAP, WIRE_PEAK_FROM_MAX, BLOCK_DROPS_LAST, BLOCK_DROPS_PARTIAL, BLOCK_WRITES_FIRST_SLOT, NORM_BY_N. SURVIVED: LOW_EDGE_TO_0+, LOW_EDGE_TO_fs/4, HIGH_EDGE_TO_5fs/2, ENDS_INCLUDED, BLOCK_FIRST_ONLY, ROWS_NO_MAX1 | -
(verified) tests/test_tellability_result.py:303-317 | run() is checked field by field against the row functions: WIRE_SWAP and WIRE_PEAK_FROM_MAX both go red here; the [two] case of test_a_legal_record_of_any_shape_ships_a_result exercises the 2-read grid and goes red under a 1/T step (n_steps 1, empty grid) | -
## Commands run (frozen tree at 24a46d8 + uncommitted delta; PYTHONPATH=<worktree>/src:<worktree>, PYTHONPYCACHEPREFIX under rev/; quebra.__file__ = this worktree's src/quebra/__init__.py; scratch in /home/sera/.claude/jobs/cbef08f8/tmp/patchA/rev/)
- delta check: git diff == delta.diff, cmp exit 0 (before and after the review)
- ruff check --no-cache . exit 0 | ruff format --check --no-cache . exit 0 (199 files) | mypy src/quebra/core (make types scope) exit 0 | lint-imports --no-cache exit 0 (1 kept, 0 broken) | mypy on the two touched analyzers exit 1, 18 errors, none on a changed line (tellability is outside make types)
- targeted: pytest the three touched test files exit 0, 111 passed
- full suite once, alone: pytest -p no:cacheprovider -q exit 0, 1153 passed, 10 skipped (rev/full_suite.log)
- mutations on a copy of the tree (rev/mut.py -> rev/mutroot, PYTHONPATH pinned to the copy and printed), 20 mutants, results in rev/mut.log, exit 0, copy restored
- probes, all exit 0: rev/p1_grid_chunk.py (chunk bit-identity, lattice parity, 2 reads, sinc/cos bounds, ends), p2_lnM.py, p3_real.py (34 6D2S + 19 2x2 records, read-only by absolute path), p4_tests_geom.py, p5_jitter.py, p6_lowedge.py, p7_empty_reuse.py
- git status: only this ledger modified; .claude/sd2-update.md was untracked before the review
## Verdict Patch A: DO NOT SHIP (0 CRITICAL, 4 IMPORTANT, 12 MINOR). The statistic is the right one and it matters: on 53 real records W(f_s) reads 0.0000-0.43 where the band maximum reads 0.007-0.88. The grid step, the >= 3-point minimum, the 2-read case, sinc^2(1/8) = 0.9496, cos^2(pi/8) = 0.8536, the (ln M)/N scale, the chunking (bit-identical for every layout) and the exception behaviour all hold. The four IMPORTANTs are small: (1) f_s is on the grid only when float rounding makes ceil(4 T f_s) even, so 13% of perfect lattices read 0.95 instead of 1, below the old W(f_s); n_steps = 2 ceil(2 T f_s) fixes it. (2) The second assumption's diagnostic, f*/f_s, reads 0.997 on a schedule that violates it. (3) The block test passes a loop that computes only the first block, because np.empty recycles the previous call's buffer; np.full(nan) fixes it. (4) Letting the main lobe into the cell survives the suite and would ship 0.81 at f*/f_s 0.0003 on a real record; one structureless-times test catches it.

# Patch B: underpowered -> inconclusive (uncommitted on f826388)

SCOPE: /home/sera/.claude/jobs/cbef08f8/tmp/patchB/delta.diff (== `git diff` in this worktree, cmp exit 0) plus every reader of the renamed values
COMMIT: f826388
## Manifest
## Reviewed
- [x] src/quebra/analyzers/check_ledger.py
- [x] src/quebra/analyzers/independence_survey.py
- [x] src/quebra/analyzers/check_outcome.py
- [x] src/quebra/analyzers/checks/cvm_cramer_von_mises.py
- [x] src/quebra/plots/theme.py
- [x] src/quebra/plots/instrument_validation_plot.py
- [x] tests/test_check_ledger.py
- [x] tests/test_check_outcome.py
- [x] tests/test_checks_cvm.py
- [x] tests/test_independence_survey.py
- [x] docs/iid_checks/*.md (6 files, +7/-7)
- [x] repository-wide residue search (Q1) and readers (Q3)
## Findings (Patch B)
IMPORTANT | src/quebra/analyzers/check_ledger.py:8 | the module docstring still says the `n_events` condition is "enough events for the check to have power". That is the claim the rename exists to withdraw, and it now contradicts the new `_verdict` docstring (:284, "Neither `pass` nor `inconclusive` measures power") and `ledger_recipe.py:51-53` ("A policy value: it is not derived from power"). AGENTS 4, fix every site of a class | `n_events  at least min_events_pass events (a policy floor, not a power threshold)`
MINOR | src/quebra/analyzers/check_ledger.py:284-285 | "how much a `pass` excludes is read from the bench's power at the record's n" points at a number the ledger does not carry (LEDGER_COLUMNS has bench size only, no power column), at an n the bench never measured (the ledger judges at the nearest grid point, `bench_n_target`), against alternatives the sentence does not name (the bench's generated arms, averaged over rho and shape). It overclaims the reach of the bench in the other direction | "... is read off `jobs/bench/results/power_table.csv` at `bench_n_target`, and holds only against the bench's generated alternatives"
(verified) src/quebra/analyzers/check_ledger.py:63-80,316-324 | constant, VERDICTS tuple and the return site renamed consistently; every reason string byte-identical to f826388; precedence unchanged (not computed > ties > fail > inconclusive > pass) | -
MINOR | src/quebra/analyzers/independence_survey.py:142-144 | the reworded example broke its own arithmetic. Before: "40 cells with no power and 4 rejections is not '10% rejection'" (readable as 40 cells, 4 of them rejections). Now: "40 inconclusive cells and 4 rejections", which is 44 cells, and 4/44 is 9.1%, not "10% rejection" | "36 inconclusive cells and 4 rejections" (40 cells, 4/40 = 10%)
IMPORTANT | src/quebra/analyzers/independence_survey.py:48-49 (and its rendered twin src/quebra/plots/independence_survey_plot.py:98-100) | says C3's `pass` "means 'did not reject' and cannot mean 'and it had the power to'", which (a) describes a cell that cannot occur: C3 has no bench cell, so `_verdict` sends every C3 non-rejection to `inconclusive` (check_ledger.py:594-604 sets accepted None; :319-320), as this file's own :147-150 says ("pass is unreachable", "46 fail, 0 pass"); and (b) implies a green cell on every OTHER instrument does mean "had the power to", the claim the new `_verdict` docstring (check_ledger.py:284) disowns. The plot version is a figure caption, so it reaches the survey figures | comment: "C3 has no bench cell, so `_verdict` can never return `pass` for it: a C3 non-rejection is `inconclusive`"; caption: drop the "green cell ... power to" clause and say C3 non-rejections read inconclusive because no level stands behind them
MINOR | src/quebra/analyzers/independence_survey.py:35-41 (also tests/test_check_ledger.py:20-27, tests/test_check_outcome.py:31-37, tests/test_independence_survey.py:22-26) | the import lists were alphabetical with UNDERPOWERED last; INCONCLUSIVE was dropped into the same slot, so they are no longer sorted. Cosmetic: ruff selects only E4/E7/E9/F, so nothing flags it | move VERDICT_INCONCLUSIVE after VERDICT_FAIL
MINOR | src/quebra/analyzers/independence_survey.py:316-333 | the renamed `survey_summary` count column is pinned by nothing: the only test (tests/test_independence_survey.py:464) asserts the row count. No reader exists (git grep: tests and this file only), so the rename breaks nothing, and a mis-keyed column would not be caught either | optional: assert the column set once
(verified) src/quebra/analyzers/independence_survey.py:112-118,136 | VERDICT_ORDER carries the new constant, so `counts` keys, both survey legends and the check-outcome legend pick it up; `counts[VERDICT_INCONCLUSIVE]` resolves | -
(verified) src/quebra/analyzers/check_outcome.py:17 | docstring only; the module reads verdicts as opaque strings (`str(cell["verdict"])`, :285) and tallies them (:178-183), so nothing there keys on the old value | -
MINOR | src/quebra/analyzers/checks/cvm_cramer_von_mises.py:45-46 (same sentence at docs/iid_checks/CvM_cramer_von_mises.md:70-71 and tests/test_checks_cvm.py:160-161) | "every CvM ledger row read `inconclusive / no bench cell`" was already too wide and the patch re-touched it: with no bench cell a REJECTION reads `fail` with the UNCALIBRATED note (check_ledger.py:303-312), and a None p-value reads `not computed`; only non-rejections read inconclusive | "every CvM non-rejection read `inconclusive / no bench cell`"
MINOR | src/quebra/plots/theme.py:313-317 | the reworded legend comment now reads inconsistently: (a) "The two 'no information' verdicts" lost its antecedent, since the `inconclusive` line no longer says "no information"; (b) the untouched line above it says `pass / fail` both mean "was calibrated at this event count", but an uncalibrated or miscalibrated rejection stays `fail` (check_ledger.py:303-312; the ledger panel hatches it for that reason) | "pass  calibrated at this event count and did not reject; fail  rejected (calibration in the note and the hatch)"; "The two desaturated verdicts (inconclusive, not interpretable) ..."
(verified) src/quebra/plots/theme.py:325-336 | Q2: `_INCONCLUSIVE_TONE` is the old hex, so no figure changes colour; "underpowered" now raises ValueError in `verdict_color`; every ledger verdict and VERDICT_ABSENT resolve (rev/render_probe.py) | -
(verified) src/quebra/plots/instrument_validation_plot.py:35-40,145 | Q2: both former "underpowered" uses (TIER_PARTIAL and a disagreeing published value) now resolve "partial"; TierMatrixPlot draws its partial cells and PublishedValuesPlot the forced-disagreeing row's bar in #B9A87A; all four instrument-validation plots render (rev/render_probe.py, render_probe2.py). TierMatrixPlot is rendered by no test in the suite, so this probe is the only evidence for that site | -
(verified) tests/test_check_ledger.py:66-82 | Q5, reviewer mutations on a copy (rev/mut.py, rev/mut.log; PYTHONPATH pinned to the copy and printed). RED in the renamed test: N_EVENTS_DROPPED, MISCAL_DROPPED, NOBENCH_DROPPED (each of the three conditions removed alone), UNMET_RETURNS_PASS, ALL_THREE_DROPPED. CONST_REVERT (VERDICT_INCONCLUSIVE = "underpowered") survives the renamed test, as it must since it compares to the constant, but goes red in the three verdict files through theme.verdict_color in test_check_outcome. Copy restored byte-identical | -
MINOR | tests/test_check_ledger.py:80-82 | the docstring promises "its reason names the condition", but the test checks only that the expected fragment is PRESENT. Mutant MISCAL_ABSORBS_NONE (`bench_accepted is False` -> `is not True`, so every unbenched row's note also claims "bench found this check miscalibrated") survives this test and all three verdict files (72 passed). A wrong reason in the ledger's notes column goes undetected | assert the whole reason: each case alone yields exactly one clause, so `assert reason == expected_reason` with the three full strings
MINOR | tests/test_check_ledger.py:49-56,66-82 | the `n_events` boundary is unpinned: mutant N_EVENTS_LE (`<` -> `<=`, so a record with exactly `min_events_pass` events reads inconclusive, against BENCH.md:87 "at least `min_events_pass` events") survives the renamed test and all three verdict files; git grep finds no test using `min_events_pass` at the boundary. Pre-existing, not introduced by the rename | add `({"n_events": 34}, "n_events 34 < 35")` to the parametrize and an `n_events=35` case to the pass test
(verified) tests/test_check_ledger.py:524-527,555-558 | docstring renames only; no assertion changed. The old test name is referenced nowhere (git grep) | -
(verified) tests/test_check_outcome.py, tests/test_checks_cvm.py, tests/test_independence_survey.py | constant and docstring renames only; no assertion's logic changed. test_check_outcome's parametrized render test is now the one test that goes red if VERDICT_INCONCLUSIVE's value is reverted (rev/mut.log, CONST_REVERT) | -


IMPORTANT | docs/iid_checks/C5_rank_autocorr.md:76-79 | the reworded sentence keeps a causal "so": C5 is underpowered (0.357 at n = 100), "so `check_ledger.py` marks a non-rejection below `min_events_pass` as `inconclusive`". The floor is a policy value "not derived from power" (ledger_recipe.py:51-53), and it does not answer C5's power: the bench ACCEPTS C5 studentized at every grid n from 35 to 500 on both clocks (rev/accept_probe.py; only in_spec 700 is flagged), so a C5 non-rejection at n = 35-100 reads `pass` at Arm E mean power 0.174-0.357 (the C5.md:69-71 figures; 0.174 at n = 35, recomputed as the mean of all 12 Arm E studentized rows at that n, which reproduces the doc's 0.106/0.236/0.357 at 20/50/100). The paragraph's own lead, "require more than its p-value before calling anything a pass", is therefore not what the ledger does for C5. Same class as check_ledger.py:8; the power prose itself is right and is kept as the brief intends | replace "so ... never `pass`" with: "The ledger does not correct for this. A non-rejection below `min_events_pass` reads `inconclusive`, but at or above it, with the bench accepting C5 at that n, it reads `pass`, and what that `pass` excludes is the power above."
MINOR | docs/iid_checks/LIMITATIONS.md:25-27, docs/iid_checks/iid_checks_basics.md:115-117 | both now say the 070423 q1 ledger's non-rejections "read `inconclusive`", but the artifact they cite still says `underpowered` in its shipped state: the newest q1 run in the main checkout (output/check_ledger_6d2s_q1_687bef_20261028_204812, q1_070423_check_ledger_data.pkl, read-only) holds 56 non-rejections = 40 underpowered + 4 ties + 12 pass, so "44 of 56" is right and the label is not, until the ledger jobs are rerun (AGENTS 4, "in the state it ships"). Out of scope by the brief; recorded so the regeneration is not forgotten | rerun check_ledger_6d2s_q* (and the survey) in the commit that follows, or land the docs with it
(verified) docs/iid_checks/BENCH.md:82-89,154-157 | the rule and its order match `_verdict` (check_ledger.py:286-324); 35 is the `CheckLedgerInputs` default (:180) | -
(verified) docs/iid_checks/C6_exchangeability.md:74-76 | "flagged at n = 700 (z = -2.98 against 2.96 over 16 cells)" recomputes from size_table.csv (z -2.975, crit 2.955, 16 cells, bench_accepted False), so a non-rejection there does read `inconclusive` | -
(verified) Q1 residue | `git grep -i underpower` over the tracked tree (minus .claude) leaves only docs/iid_checks/C5_rank_autocorr.md:68 and :78 (power prose, kept as briefed), spec/quebraplan.md:668 and a session transcript under spec/ledger/. No residue in src/, tests/, jobs/ (bench report included), scripts/, README/AGENTS, or the tracked .claude/qre_checks_reference.tex. No schema, ordering tuple or legend names a verdict by literal: VERDICTS (ledger legend), VERDICT_ORDER (survey legends, counts, check-outcome legend) and VERDICT_COLORS are the only lists, and all three carry the new key. The survey-summary column is listed nowhere else | -
MINOR | spec/quebraplan.md:668 | "a stated mapping for `underpowered` and `not computed`" names the old verdict as a label in a live open question (section 7 item 10), which a reader of the ADR it asks for will search the code for. AGENTS 10 lets spec/ keep older wording as history, so this is optional | "`inconclusive` (formerly `underpowered`)"
(verified) Q3 readers | kaplan_meier.KaplanMeierComparison.check_verdicts and check_attachment.summarise treat verdicts as opaque strings (tally by value, check_attachment.py:43-46), filled in jobs/active/km_with_checks_6d2s.py:215-220 from CheckOutcome grids; panels/check_ledger.py:132, independence_survey_plot.py:142 and check_outcome_plot.py:136 key only on pass/fail. No job, script or test pins the literal "underpowered"; nothing outside the test reads the summary column. Cached artifacts cannot reach the new code through the runner: reuse needs a matching commit and a clean tree at both ends (core/runner.py:256-284) | -
MINOR | src/quebra/plots/independence_survey_plot.py:120-127 (twin at :291-296) | the colour lookup sends ANY string outside VERDICT_ORDER to the `not computed` tone, and `InstrumentGrid.counts` (independence_survey.py:136) counts only VERDICT_ORDER, so an old-vocabulary verdict is silently redrawn as `not computed` and counted nowhere. Measured on the pre-rename q1 070423 ledger (main checkout output, read-only): 40 of 180 survey cells counted under no verdict and drawn grey, no error, while CheckLedgerPanel raises `Unknown verdict: 'underpowered'` on the same artifact (rev/artifact_probe.py). Pre-existing; the rename is what makes the branch live. The runner's reuse rule keeps it off the normal path | let only non-strings (NaN, the no-row case) take the `not computed` tone and pass every string to `theme.verdict_color`, so a stale label raises as it does in the ledger panel
MINOR | src/quebra/plots/independence_survey_plot.py:94-96 | the rendered caption says "Grey cells carry no information - a non-rejection on too few events is not evidence of independence", which describes `inconclusive` by one of its three conditions, the narrowing the rename exists to remove (the bench-miscalibrated and no-bench-cell routes also produce it, and on C3 they are the only route). The `inconclusive` tone, #B9A87A, is also khaki rather than grey | "Khaki cells did not reject but cannot be read as a pass (too few events, or no accepted bench cell at this n); grey cells carry no answer."
## Commands run (frozen tree at f826388 + uncommitted delta; PYTHONPATH=<worktree>/src:<worktree>, PYTHONPYCACHEPREFIX under rev/; quebra.__file__ = this worktree's src/quebra/__init__.py; scratch in /home/sera/.claude/jobs/cbef08f8/tmp/patchB/rev/)
- delta check: git diff == delta.diff, cmp exit 0 before; after, git diff excluding this ledger == delta.diff, cmp exit 0
- ruff check --no-cache . exit 0 | ruff format --check --no-cache . exit 0 (199 files) | lint-imports --no-cache exit 0 (1 kept, 0 broken) | mypy src/quebra/core (make types scope) exit 0
- targeted: pytest test_check_ledger, test_check_outcome, test_independence_survey, test_checks_cvm, test_instrument_validation exit 0, 116 passed
- full suite once, alone (no other pytest running): pytest -p no:cacheprovider -q exit 0, 1156 passed, 10 skipped (rev/full_suite.log)
- render probes exit 0: rev/render_probe.py (CheckLedgerPanel, all 9 survey plots, overview, CheckOutcomePlot, all 4 instrument-validation plots incl. a forced disagreement; "underpowered" raises in verdict_color), rev/render_probe2.py (inconclusive and partial cells drawn in #B9A87A)
- mutations on a copy (rev/mut.py -> rev/mutroot, PYTHONPATH pinned to the copy and printed), 8 mutants, rev/mut.log, exit 0, copy restored byte-identical; CONST_REVERT re-run to name the catching test
- probes exit 0: rev/accept_probe.py (bench acceptance for C5/C6 by n; C6 n 700 numbers), rev/artifact_probe.py (pre-rename q1 070423 ledger from the main checkout's output/, read-only), power_table recompute for C5 Arm E
- git status: only this ledger modified; .claude/sd2-update.md was untracked before the review
## Verdict Patch B: DO NOT SHIP (0 CRITICAL, 3 IMPORTANT, 12 MINOR). The rename is mechanically complete and safe: no verdict-label residue in src/, tests/, jobs/ or docs/ beyond the briefed power prose; every draw site resolves, including TierMatrixPlot, which no test renders; nothing reads the old string or the survey column; reuse cannot feed a stale ledger to the new code; and the renamed test goes red when any one of the three conditions is removed. All three IMPORTANTs are one-sentence prose fixes in the class this patch exists to close, a verdict or the event floor read as a power statement: (1) check_ledger.py:8 still says the n_events condition means "enough events for the check to have power", in the same file as the new docstring that says it does not; (2) the C3 comment and the survey caption it feeds (independence_survey.py:48-49, independence_survey_plot.py:98-100) describe a green C3 cell that cannot occur and imply a green cell elsewhere means "it had the power to"; (3) C5.md:76-79 presents the policy floor as the answer to C5's low power, while the ledger prints `pass` for C5 at n 35-100 at mean power 0.17-0.36. Fix those three and it ships. Of the MINORs, the cheapest to take with them are the 40+4 arithmetic (independence_survey.py:143), the "every CvM ledger row" wording, and exact-reason assertions in the renamed test (MISCAL_ABSORBS_NONE survives it today).
