# C6 - permutation test of exchangeability

## What it computes

A Ljung-Box-type sum over the same rank autocorrelations C5 uses, each lag centred at and
scaled by its own permutation moments, calibrated by permutation:

```
z_h = (r_h - mean_perm(r_h)) / sd_perm(r_h)
statistic = sum_h z_h^2
```

The moments are taken over the observed value together with its permutations, so `z_h` is
exactly the value whose absolute maximum C5's studentized variant takes. Where C5 asks "is
any single lag extreme", C6 asks "is the whole autocorrelation profile further from its
permutation centre than reordering would produce". Its p-value is exact under
within-segment exchangeability at any number of permutations, as C5's is.

**Why each lag is centred.** Under exchangeability an autocorrelation does not average zero:
Dufour and Roy (1985, Eq. (2.2), p. 260) give its exact mean on one unbroken series, and
recommend the exact means and variances for tests of randomness (p. 271). Ljung and Box
(1978, Eq. (1.1), p. 297) sum squares that are not centred. On a gapped record whose
segments sit at different levels, every within-segment product leans positive whatever the
order inside each segment, so the null mean of `r_h` moves well away from zero. A sum of
uncentred squares is not symmetric about that offset: dependence of one sign inflates it
and the other sign deflates it. The size stays exact either way; what the offset distorts
is power. Dufour and Roy's closed forms do not cover a gapped record, so the permutation
moments take their place.

`extra["z"]` carries each lag's `z_h` with its sign, since the sum itself is symmetric in
the sign of every `z_h`.

## Where it lives

`analyzers/checks/c6_exchangeability.py`, sharing C5's harness in
`analyzers/checks/_rank_serial.py` and the same permutation set via
`analyzers/checks/battery.py`.

## Limitations

**Its null is weaker than it sounds.** Exchangeability is implied by iid but does not imply
it: a sequence can be exchangeable and still not independent. A non-rejection is evidence
that ORDER carries no information, not that the durations are independent.

**A level difference between segments is part of the null.** The permutations never move a
duration across a gap, so no permutation-based check (C1, C2 and CvM by permutation, C5,
C6) can see a record whose segments differ only in level. Here the centring absorbs such a
difference by construction.

**A rank statistic, so ties degrade it** - identically to C5. See that page.

**A lag with no spread is not computed.** When every permutation gives the same `r_h` at
some lag, as when every segment holds two events, or every segment that pairs at that lag
is constant inside, `studentize` raises. The battery raises with it, so the ledger writes
the whole cell as one `not computed` row and reports no other check for it, C3 included.

**Power is low below n = 100.** Averaged over the bench's dependence grid
(`jobs/bench/results/power_table.csv`, arm=E_copula_ar1_durations, mean over `rho` and
shape), C6 has 22% power at n = 50, 34% at n = 100, 60% at n = 355 and 79% at n = 1000,
slightly below C5 studentized at every n of the grid. On that grid its rejection rate
differs from C5 studentized's by 0.017 on average and by at most 0.077 (`power_table.csv`,
Arm E rows). Over all 288 power cells the two share (108 on Arm C, 72 on Arm D, 108 on
Arm E) the figures are 0.019 and 0.172, that maximum on an Arm D trend cell
(`jobs/bench/report.py::check_agreement`).

## What we do

Report it, and require the same three conditions before a `pass`. The bench compares the
other checks against it (`jobs/bench/report.py::check_agreement`); that choice of reference
is a convention, since every permutation row is exact under the same null. The bench puts it
on HOLD: correctly calibrated across the 144 null cells inside the envelope, worst 0.0355 on
Arm A, in-spec, n = 700, continuous and uncensored (z = -2.98 against 3.58), but with mean
power 0.343 at n = 100 (`promotion_report.md`, verdicts and size per row). Judged per event
count, as the ledger judges it, that cell is flagged at n = 700 (z = -2.98 against 2.96 over
16 cells): conservative, not oversized, and a non-rejection there reads `inconclusive`
(`BENCH.md`, section 5, item 8).
