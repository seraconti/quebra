# C6 - permutation test of exchangeability

## What it computes

A portmanteau over the same rank autocorrelations C5 uses, weighted by the pair count at
each lag, calibrated by permutation:

```
statistic = sum_h n_pairs(h) * r_h^2
```

Where C5 asks "is any single lag extreme", C6 asks "is the whole autocorrelation profile
larger than reordering would produce". It is the natural reference among the six: a pure
permutation test of exchangeability, resting on no asymptotic approximation and no citation
beyond the permutation principle itself.

## Where it lives

`analyzers/checks/c6_exchangeability.py`, sharing C5's harness in
`analyzers/checks/_rank_serial.py` and the same permutation set via
`analyzers/checks/battery.py`.

## Limitations

**Its null is weaker than it sounds.** Exchangeability is implied by iid but does not imply
it: a sequence can be exchangeable and still not independent. A non-rejection is evidence
that ORDER carries no information, not that the durations are independent.

**A rank statistic, so ties degrade it** - identically to C5. See that page.

**Underpowered unless the dependence is strong.** Averaged over the bench's dependence grid
(`jobs/bench/results/power_table.csv`, arm=E_copula_ar1_durations, mean over `rho` and shape), its power is 0.086 at n = 20, 0.207 at n = 50, 0.340 at n = 100 and
0.599 at n = 355, below the studentized C5 at every n, which is the price of testing a broad
null with a portmanteau.

**It is nearly redundant with C5.** Mean absolute difference in rejection rate against C5
is 0.014-0.016 across 192 shared power cells; the maximum is 0.14. They are measuring
nearly the same thing.

## What we do

Report it, treat it as the reference the others are compared against (that is what
`jobs/bench/report.py::check_agreement` does), and require the same three conditions before a
`pass`. The bench put it on HOLD: calibrated (worst z = -2.46 against 3.45), with mean
power 0.340 at n = 100 over its dependence grid.

Its worst null cell in the bench came from Arm D at b = 1 - the trend-renewal generator
with no trend - which is the cross-check that the TRP generator itself is sound.
