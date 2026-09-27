# C3 - copula-based serial independence

**Runs, but uncalibrated.** R is installed and C3 has executed. Smoke test on iid
exponential input at `lag.max=5`, `n.sim=1000`, C3 `seed=1`: n=50 gives statistic 0.0057908
/ p 0.9575, n=150 gives 0.0071326 / p 0.9036, n=355 gives 0.0076331 / p 0.8656.

**The data recipe is part of the number**, as much as the seed is. Those three series come
from ONE `numpy.random.default_rng(0)` drawing `exponential(1.0)` sequentially in the order
n=50, then n=150, then n=355 - not re-seeded per n, and not prefixes of one long draw. A
statistic quoted without its recipe cannot be reproduced. The C3 `seed` is a separate knob:
the R side simulates its own null, so the statistic is seed-invariant while the p-value is
not, and at `seed=0` the first two series give 0.9476 and 0.9116.

Wall clock at `n.sim=1000` runs from seconds at n=50 to about a minute at n=355.
`analyzers/checks/c3_serial_copula.py` is the one place the cost is quoted, so that a second
copy cannot drift from it. It is load-dependent: measure it for the n you actually have.

Failing to reject data that satisfies the null is the expected outcome and is a SMOKE TEST,
NOT calibration. C3 still has no size or power evidence, because the bench never ran it and
`battery.ROW_KEYS` has no C3 row. **Its silence is not a pass.**

## What it tests

Genest and Remillard's empirical-copula test of serial independence, based on the Mobius
decomposition of the independence hypothesis. Implemented in R as
`copula::serialIndepTest`, which needs a simulated null from `serialIndepTestSim(n, lag.max)`.

It is here because it tests a strictly STRONGER null than C5 or C6: full serial
independence at all lags jointly, not merely zero rank autocorrelation. Two duration
sequences can have zero rank autocorrelation at every lag and still be dependent; this
would see that and C5 would not.

## Where it lives

`analyzers/checks/c3_serial_copula.py` plus its sibling `c3_serial_copula.R`.

**A bridge, not bindings.** CSV out, `Rscript`, CSV back - no `rpy2`. rpy2 pins an ABI
against a specific R build and turns "R is missing" into an import-time failure of the whole
package; a subprocess turns it into a per-call `None`. The pipeline must import cleanly on
a machine without R.

## Limitations

**Uncalibrated in the only sense that matters.** `tests/test_checks_c3_bridge.py` pins the
decline paths: a missing interpreter returns `p_value=None` with `notes="R unavailable"`,
nothing raises, and the ordinary data guards still fire. It does not pin the numeric path.
The CSV round-trip, the exit-code handling and the result parsing have all been exercised by
the smoke test at the top of this file, but no test in the suite re-runs them, so a change to
either side of the bridge would go uncaught. Beyond that, what is missing is everything the
bench would have measured: size, power, and behaviour on data that is NOT iid.

**No multi-process form.** `serialIndepTest` takes one series. Concatenating segments would
manufacture lag pairs across read gaps - exactly the relation the carve refuses to assert -
and combining per-segment p-values needs a dependence-free combination rule that is not
established for this statistic. So `m > 1` returns `None` with that reason recorded rather
than a number nobody can defend.

**The continuity assumption.** The distribution-free property is derived for continuous
observations. On a quantised metric the durations are integers times the read interval and
the empirical copula is not what the test assumes. The ledger's tie rule applies to C3 for
this reason, though it has never fired on the T2* ladder.

## What we do

Ship the bridge, call it, and report `not computed` where R is absent. The promotion report
scores five checks, not six, and says so at the top. It asserts nothing about the interpreter
either way: the report is TRACKED and generated from two committed CSVs, so probing the machine
would make a committed artifact differ per machine. C3's absence from it rests on the bench, not
on whether R is installed.

`jobs/composite/independence_survey.py` sets `INCLUDE_C3 = True`, so where R is present the
survey draws a real C3 p-value in every eligible cell. Read a green C3 cell as "did not
reject" and nothing stronger: there is no bench cell behind it to say whether that cell had
the power to reject.
