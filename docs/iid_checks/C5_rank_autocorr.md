# C5 - studentized rank autocorrelation

## What it computes

The rank autocorrelation of the durations at lags 1..H, studentized by the permutation
ensemble's own spread, with the MAXIMUM over lags calibrated by permutation:

```
r_h = sum_{pairs at lag h} (R_i - Rbar)(R_{i+h} - Rbar) / sum_i (R_i - Rbar)^2
statistic = max_h |r_h - mean_perm(r_h)| / sd_perm(r_h)
```

Pairs never cross a segment boundary: a lag-1 pair spanning a read gap would relate two
windows separated by unobserved hours.

Studentizing by the ensemble's own moments is this project's construction. It is not the
studentization of Romano and Tirlea (Journal of Time Series Analysis, 2022, as
`.claude/qre_checks_reference.tex` cites it; UNVERIFIED, the paper was not opened), which
divides by an estimate of the asymptotic standard error. The moments are computed from the observed
value TOGETHER WITH the permutations, which makes the scaling a symmetric function of the
augmented sample and keeps the test exactly valid at finite `B` rather than valid to
O(1/B).

## Where it lives

`analyzers/checks/c5_rank_autocorr.py`; the shared rank machinery is in
`analyzers/checks/_rank_serial.py`.

## Limitations

**Studentization does much less than expected, and the reason is measurable.** Every lag
shares one denominator (the full centred sum of squares) while its numerator sums only the
within-segment pairs at that lag, so the null spread of `r_h` GROWS with the pair count -
the opposite of the usual intuition. Measured:

| configuration | lags, pairs | null SD | ratio |
|---|---|---|---|
| one segment of 60 | 1-5, 59..55 | 0.1269..0.1224 | 1.04 |
| 10 segments of 5 | 1-4, 40..10 | 0.0985..0.0600 | 1.64 |

So on an ungapped record studentizing barely changes each lag's scale, and only on short,
numerous segments does the unstudentized max concentrate on one lag (43% of replicates on
lag 1, against 21-29% once each lag is on its own scale). Studentizing also moves each
lag's centre, and that separates the two variants even on one segment: in
`jobs/bench/results/power_table.csv` the unstudentized variant rejects less often than the
studentized one in 84 of the 88 Arm E cells where they differ, and in 63 of 65 on Arm D,
all single segments. The unstudentized variant runs in the bench so that these are numbers
rather than claims.

**The unstudentized variant runs only in the bench.** It is uncentred: under
exchangeability its `r_h` averages `-(n - h)/(n(n - 1))` on one segment, `-1/n` at lag 1
(Dufour and Roy 1985, Eq. (2.2), p. 260), and moves up with any difference in level between
segments. Its max then leans as an uncentred sum does, so its power depends on the sign of
the dependence. The pipeline's battery leaves it out (`battery.ROW_KEYS`), so no ledger,
survey or band carries it; the bench asks for it with `include_c5_raw=True`
(`battery.BENCH_ROW_KEYS`). `LIMITATIONS.md`, section 10, covers it and `extra["r_lag1"]`.

**A lag with no spread is not computed.** When every permutation gives the same `r_h` at
some lag, as when every segment holds two events, or every segment that pairs at that lag
is constant inside, the studentized variant raises rather than scale by rounding noise. The
battery raises with it, so the ledger writes the whole cell as one `not computed` row.

**It is a rank statistic, so ties degrade it.** On a quantised metric where most windows are
one read long the ranks are mostly ties and the statistic loses meaning. The ledger's tie
rule applies. On the T2* ladder it never binds: durations are wall-clock seconds and 355
windows gave 355 distinct values.

**Underpowered unless the dependence is strong.** This is the binding limitation. Averaged
over the bench's dependence grid (`jobs/bench/results/power_table.csv`, arm=E_copula_ar1_durations, mean over `rho` and shape), the studentized variant's power is
0.106 at n = 20, 0.236 at n = 50, 0.357 at n = 100, 0.626 at n = 355 and 0.804 at
n = 1000. Below n = 100 it
exceeds 0.5 only at the strongest grid point, `rho` = 0.5.

## What we do

Report it, and require more than its p-value before calling anything a pass. The bench put
C5 on HOLD - correctly calibrated (studentized variant: worst z = -3.18 against 3.58) but
underpowered (mean power 0.357 at n = 100) - so `analyzers/check_ledger.py` marks a
non-rejection below `min_events_pass` as `underpowered`, never `pass`.

The studentized variant takes the largest `|z_h|`, and C6 sums the squares of the same
`z_h`. Over the 288 power cells they share, their rejection rates differ by 0.019 on
average and by at most 0.172 (`promotion_report.md`, criterion 4).
