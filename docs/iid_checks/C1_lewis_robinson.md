# C1 - Lewis-Robinson test for trend

## The equation

Kvaloy & Lindqvist, *A class of tests for trend in time censored recurrent event data*,
Technometrics 62(1):101-115 (2020), eq (4), p. 103, single process:

```
LR = (1/gamma_hat) * sqrt(12)/(tau*sqrt(N)) * [ sum_i T_i - (N/2)*tau ]
```

Multi-process, eq (13) combines per-process statistics with weights `A_j`. The journal and
its preprint arXiv:1802.08339v1 use DIFFERENT weights, so for `m > 1` they are different
statistics. Both ship, selected by `weights=`:

```
WEIGHTS_TECHNOMETRICS (default)  journal eqs (14)-(16), p. 105, A_j ~ sqrt(N_j)/gamma_hat_j
LR_m = sqrt(12) / sqrt( sum_k N_k / gamma_hat_k^2 )
       * sum_j (1/gamma_hat_j^2) [ sum_i T_ij/tau_j - N_j/2 ]

WEIGHTS_ARXIV_V1                 preprint eqs (14)-(16), p. 8, A_j ~ gamma_hat_j tau_j sqrt(N_j)
LR_m = sqrt(12) / sqrt( sum_k gamma_hat_k^2 * tau_k^2 * N_k )
       * sum_j [ sum_i T_ij - (N_j/2) * tau_j ]
```

The journal weights are the ones behind the paper's `m = 6` example (Section 8.2, Table 4,
LR p = 0.019), which `tests/test_checks_published_values.py` reproduces; the preprint
weights give 0.0073 on the same data.

`sum_i T_i` is the total of the event times; `(N/2)*tau` is its expectation when events are
uniform on `[0, tau]`, which is the null of no trend. So LR is a centred, scaled "are the
events early or late" contrast - negative for a decreasing intensity, positive for an
increasing one, two-sided here.

The `1/gamma_hat` is what makes it a test of TREND rather than of the Poisson assumption:
dividing by the estimated coefficient of variation removes the renewal distribution's own
dispersion. That is the whole difference between Lewis-Robinson and the plain Laplace test
(Section 3.1, p. 103).

## Where it lives

`analyzers/checks/c1_lewis_robinson.py`. One formula covers both weightings and every `m`;
setting `m = 1` gives eq (4) exactly for either weighting, so a separate single-process
path would be a second copy of the same formula, free to drift.
`tests/test_checks_statistics.py::test_eq16_reduces_to_eq4_for_a_single_segment` asserts the
identity numerically rather than trusting the algebra (agreement to 1e-12).

## Limitations

**The asymptotic calibration is oversized at the event counts this project has.** It is
N(0,1) in the limit; on the bench's primary null configuration (arm A, in-spec clock,
unquantised, uncensored) it rejects at n = 20 at 0.058 for Weibull shape 0.75 and 0.0635 for
shape 1.5 against a nominal 0.05 (their pooled mean is 0.061), and stays within 0.009 of
nominal at every n from 35 on. Kvaloy & Lindqvist say the same: "most
of the tests are a bit nonconservative" for small samples (Technometrics Section 6.1,
p. 107; arXiv v1 Section 5.1, p. 10 words it "all tests are a bit non-conservative for small
samples in the underdispersed case"), and their Figure 1 (p. 107) shows LR near 0.08 at 10
expected events. The bench figures are from `jobs/bench/results/size_table.csv`
(check=c1_lewis_robinson, calibration=asymptotic, arm=A_iid_weibull, clock=in_spec,
quantised=False, censoring_target=0.0).

**The multi-process asymptotic collapses under segmentation.** With many short segments the
normal approximation fails outright: at a realised censoring of 0.17 the bench measured a
rejection rate of **0.8625** against a nominal 0.05 (Arm A, shape 1.5, n = 355,
unquantised). Among the adequately supported c = 0.25 cells the worst is 0.8695 (Arm D,
shape 1.5, n = 355); a quantised Arm A cell with 819 usable replicates of 2000 reaches
0.973. The
permutation form on the Arm A cell stayed at 0.0435, and within 0.0425-0.0645 on every
adequately supported c = 0.25 cell (`size_table.csv`, check=c1_lewis_robinson,
censoring_target=0.25). The
paper's asymptotics let every `tau_j` grow; many short segments is the other limit. Whether a
renewal start-up bias, carried by each short segment that starts at an event, accumulates
over segments is an open question: the bench has not measured it.

**`gamma_hat` is undefined on a fully tied duration vector.** A quantised metric where every
window is one read long gives zero dispersion, and the statistic it scales has no meaning.
This raises rather than returning a number.

## What we do

Use the permutation calibration. It is exactly valid rather than asymptotically - `tau`,
`N` and `gamma_hat`, hence the weights, are all invariant to reordering gaps within a
segment, so the only thing a permutation moves is `sum_i T_ij`. The bench PROMOTEd it:
calibrated across 80 null cells inside the censoring envelope (worst z = -2.15 against a
Bonferroni threshold of 3.42), with mean power 0.596 at n = 100 across the trend grid
(range 0.170-0.999 - the spread is over trend strength and Weibull shape, and mild trends
are genuinely hard).

The asymptotic form is REJECTed by the bench and the ledger hatches any cell that uses it
at an event count where it was miscalibrated.
