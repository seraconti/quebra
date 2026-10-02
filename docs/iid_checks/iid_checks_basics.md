# The iid checks - what they test, and what their answers are worth

Working reference for `analyzers/checks/` and the calibration study in `jobs/bench/`.

SPEC 0008 decided the licence wiring, and "how a check verdict reaches a panel" has folded into
`docs/PANEL_CONTRACT.md` as that sentence promised. The rest of this page is DELIBERATELY KEPT
rather than deleted, which is a stated deviation from the original plan: these pages are the
citation home for the assumption records in `analyzers/assumptions.py`, and they carry the
equation numbers and the per-check limitations that a reader checking a transcription needs.
Deleting them would move that material nowhere.

The decision itself, so this page does not have to be read for it: there is no licence and no
gate. A verdict annotates a figure and never suppresses one.

This is the first directory in `docs/` to carry citations. That is deliberate: these are
implementations of published statistics and the equation numbers are load-bearing - a
reader checking the transcription needs to know which equation, in which paper.

---

## The question all six are asking

The panel carves a metric into in-spec **windows**. Reliability arithmetic on those windows
- occupancy, survival, mean time between failures - assumes the durations behave like a
renewal process: independent, identically distributed, no trend. If they do not, the
arithmetic still produces numbers, and the numbers are wrong in ways nothing else catches.

Six checks ask whether that assumption survives contact with the data.

| check | null it tests | rejects when |
|---|---|---|
| C1 Lewis-Robinson | renewal process, no trend | events crowd early or late |
| C2 Anderson-Darling | renewal process | event times are not uniform on `[0, tau]` |
| C3 copula serial independence | full serial independence at all lags | any lag shows dependence |
| C5 rank autocorrelation | zero rank autocorrelation | consecutive durations covary |
| C6 exchangeability | order carries no information | any departure from exchangeability |
| CvM Cramer-von Mises | renewal process | event times are not uniform on `[0, tau]` |

They are not redundant, but two of them are close. Measured against C6 as the reference,
C5 agrees with it to a mean absolute difference of 0.017-0.018 in rejection rate over 288
shared power cells, while C1/C2 differ from it by 0.256-0.267 over the 234 they share. C5 and
C6 are largely measuring the same thing; C1/C2 measure a different one.

## Two clocks, because the mapping is not unique

A window table can become a renewal process two ways, and they ask different questions.

- **In-spec clock** (`CLOCK_IN_SPEC`): an event is a window dying by `down_crossing`, `x`
  is its duration, `tau` is the segment's total in-spec time. Asks whether successive
  in-spec lifetimes look renewal.
- **Calendar clock** (`CLOCK_CALENDAR`): an event is a window BIRTH, `x` is wall-clock time
  between births, `tau` is the observed length. Asks whether failures arrive as a renewal
  process in real time - the question a maintenance schedule poses.

Both are reported. Reporting one alone would hide that the answer depends on the choice,
and on the real record they DO disagree: at 4 µs the calendar clock passes every check
while the in-spec clock cannot compute C1 or C2 at all.

That last point is structural, not a bug. On the in-spec clock, time stops accruing the
moment the record ends out of spec, so `tau == T_N` and eq (7) (Kvaloy and Lindqvist 2020, p. 104) is `+inf`: its
summand `ln((tau - T_{N-1})/(tau - T_N))` has a zero denominator. Measured on synthetic iid reads: ~73% of replicates. Those rows read `not computed`
with the reason attached.

## Segments, because a read gap is not an interval

A gap in the reads means unobserved time. Treating the stretch across it as one inter-event
interval invents evidence. So a record is split at its gaps into independent time-censored
processes and combined by the multi-process forms (Kvaloy and Lindqvist 2020, eqs (13)-(16), p. 105).

The split uses `WindowsResult.diagnostics["gap_spans_s"]`, not the birth taxonomy alone:
`analyzers/windows.carve` only emits `gap_resume` next to an in-spec read, so a gap flanked
by out-of-spec reads leaves NO trace in the window table. Reading births alone merged two
processes separated by 488 s of unobserved time and handed that 493 s span to eq (4) as a
renewal interval.

## Permutation, not asymptotics

Every check ships a permutation calibration and it is the one to use. Permutation is
**exactly** valid here rather than asymptotically: `tau`, `N` and `gamma_hat` are all
invariant to reordering gaps within a segment, so the conditional null is exact at finite
`B`. The p-value is tie-corrected, `p = (1 + #{null >= observed}) / (1 + B)`, because this
project's duration vectors are heavily tied on quantised metrics and `#/B` would count a
tie as a refutation.

Permutation requires an explicit seed. `block_permutations` raises on `rng=None`: it used
to default to OS entropy, which made three consecutive runs on identical input return
p = 0.3860 / 0.4040 / 0.3790 while the provenance record stayed identical.

## What the bench establishes, and what it does not

`jobs/bench/` measures empirical size and power at the event counts this project actually has,
over 504 cells and 720,000 replicates. Its verdicts:

| | verdict | why |
|---|---|---|
| C1, C2 permutation | PROMOTE | size holds across 128 null cells; mean power 0.60/0.62 at n = 100 |
| C1, C2 asymptotic | REJECT | oversized inside the envelope (worst z = 31.7 and 5.4) |
| CvM permutation | PROMOTE | size holds across 144 null cells (worst z = 2.56) |
| CvM asymptotic | REJECT | oversized inside the envelope (worst z = 4.41 over 86 null cells) |
| C5, C6 | HOLD | correctly calibrated, but mean power at n = 100 is only 0.33-0.36 |
| C3 | RUNS, UNCALIBRATED | exercised under Rscript 4.5.3 with `copula`; smoke-tested on iid input only - no bench cell, no size or power evidence |

**The number to read before trusting a non-rejection**: averaged over the bench's dependence
grid, C5 and C6 have 21-24% power at n = 50, 33-36% at n = 100, 60-63% at n = 355 and 79-80%
at n = 1000, and below n = 100 they exceed 50% only at the strongest dependence the grid
tests (`promotion_report.md`, "Power across the dependence grid"). A non-rejection at a
threshold with 50 windows is close to uninformative unless the dependence is strong. A
rejection still means something; the silence does not.

That asymmetry is why `analyzers/check_ledger.py` requires three conditions for a `pass`
and not one. On the 070423 record of qubit 1, 49 of its 63 non-rejections read
`underpowered` or `not interpretable (ties)`, and would have printed `pass` under a p-value-only
rule (`jobs/active/check_ledger_6d2s_q1.py`, its `q1_070423_check_ledger_data`).

## Files

- `analyzers/checks/` - the six checks, the permutation harness, the segment mapping
- `analyzers/checks/battery.py` - runs the five permutation checks off ONE permutation set
- `analyzers/check_ledger.py` - scores each answer against event count, calibration, ties
- `jobs/bench/` - the calibration study; `jobs/bench/results/promotion_report.md` is its output
- `analyzers/calibration_summary.py` - the shared definition of "calibrated at this n"

## Per-check pages

[C1](C1_lewis_robinson.md) - [C2](C2_anderson_darling.md) - [C3](C3_serial_copula.md) -
[C5](C5_rank_autocorr.md) - [C6](C6_exchangeability.md) -
[CvM](CvM_cramer_von_mises.md) - [limitations](LIMITATIONS.md)

There is no C4. Lin-Wei-Ying is C4 in the numbering of `.claude/qre_checks_reference.tex`, the repository's own design note, and is not
implemented here; nothing in the battery depends on it.

## Sources

- J. T. Kvaloy and B. H. Lindqvist, *A class of tests for trend in time censored recurrent
  event data*, Technometrics 62(1):101-115, 2020, doi:10.1080/00401706.2019.1605936.
  Equations (4), (6), (7), (10), (11), (13)-(16); Sections 3.6, 4.2, 6.1, 8.1, 8.2; Appendix
  A.2. The Technometrics paper's preprint, arXiv:1802.08339v1 (2018), has the same title but different section numbers,
  a different second example, and DIFFERENT weights in eqs (14)-(16). `docs/GOLD_STANDARD.md`
  section 2 maps the locators this repository uses between the two versions.
- G. Marsaglia and J. Marsaglia, *Evaluating the Anderson-Darling Distribution*,
  Journal of Statistical Software 9(2):1-5, 2004, doi:10.18637/jss.v009.i02. The `adinf`
  limiting function (p. 3), and the correction of the published 1% point from 3.857 to
  3.878125 (p. 2).
- T. W. Anderson and D. A. Darling, *Asymptotic theory of certain "goodness of fit" criteria
  based on stochastic processes*, Annals of Mathematical Statistics 23:193-212, 1952. The
  limiting Cramer-von Mises distribution, eq (4.35) p. 202 and Table 1 p. 203; the limiting
  Anderson-Darling distribution, eq (4.38) p. 204.
- M. A. Stephens, *EDF statistics for goodness of fit and some comparisons*, JASA
  69(347):730-737, 1974. Table 1A part 1.0, p. 732: the tabled critical values 1.933, 2.492,
  3.070, 3.857 for A^2 and 0.347, 0.461, 0.581, 0.743 for the modified W^2.
- C. Genest and B. Remillard, *Tests of independence and randomness based on the empirical
  copula process*, Test 13:335-369, 2004, as the references of `copula::serialIndepTest`
  (copula 1.1-7) give it. UNVERIFIED: the paper was not opened. Implemented in R as
  `copula::serialIndepTest`.
- S. Chatterjee, *A new coefficient of correlation*, JASA 116(536):2009-2022, 2021,
  doi:10.1080/01621459.2020.1758115. Used by `analyzers/shape_stats.py`, not by these checks.
- B. H. Lindqvist, G. Elvebakk and K. Heggland, *The trend-renewal process for statistical
  analysis of repairable systems*, Technometrics 45:31-44, 2003. The model behind the bench's
  Arm D.
- D. R. Cox and P. A. W. Lewis, *The Statistical Analysis of Series of Events*, 1966, for
  the renewal-versus-trend framing.
