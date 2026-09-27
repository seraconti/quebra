# Limitations that cut across all six checks

Per-check limitations live on the per-check pages. These are the ones that would mislead a
reader who took any single ledger row at face value.

## 1. A non-rejection is usually not evidence

This is the binding one. At the dependence this instrument actually shows (duration lag-1
0.12-0.15), C5 and C6 have:

| n | 20 | 35 | 50 | 75 | 100 | 355 |
|---|---|---|---|---|---|---|
| C5 studentized | 0.070 | 0.101 | 0.135 | 0.192 | 0.292 | 0.849 |
| C6 | 0.060 | 0.092 | 0.117 | 0.168 | 0.249 | 0.782 |

Below n = 75 a non-rejection is the expected outcome whether the durations are dependent or
not. A REJECTION is meaningful, because the checks are correctly calibrated. It is the
silence that carries no information.

`analyzers/check_ledger.py` enforces this: a `pass` requires the p-value AND a sufficient
event count AND that the bench found the check calibrated at that count. On the 0704
record, 31 of 68 non-rejections would have printed `pass` under a p-value-only rule.

## 2. Multiplicity across the ladder is not corrected

A ledger runs 5 checks x 2 calibrations x 2 clocks over up to 10 thresholds. Nothing
corrects for that. The thresholds are nested (a window in spec at 3 µs is in spec at 4 µs),
so the tests are strongly positively dependent and a Bonferroni correction would be far too
conservative - but "far too conservative" is not "unnecessary". **Read the ladder as a
pattern, not as ~200 independent decisions.** A single isolated rejection at one threshold
is much weaker evidence than a monotone trend down the ladder.

## 3. The two clocks can disagree, and neither is wrong

On the 0704 record at 4 µs, the calendar clock passes every check and the in-spec clock
cannot compute C1 or C2 at all. That is not a contradiction: they ask different questions
(see the README), and the in-spec clock's `tau == T_N` collapse is structural for a record
that ends out of spec.

## 4. C3 is unassessed

Where R is present C3 executes; where R is absent the bridge still degrades to
`p_value=None` and the row reads `not computed`. What remains missing is CALIBRATION, not the
interpreter: C3 has no bench cell, so its size and power are unmeasured. Two operational
limits: `--vanilla` implies `--no-environ` and so needs `R_LIBS` passed explicitly or
`copula` is invisible, and the run cost climbs steeply with n - seconds at n = 50 against
about a minute at n = 355, at `n.sim = 1000`; `analyzers/checks/c3_serial_copula.py` carries
the cost and is the one place it is quoted. Whether a much larger window can reach the 900 s
timeout is OPEN, and the two records that bear on it disagree: one unverified observation has
n = 682 still running at 580 s, while the shipped ledger in
`output/check_ledger_q1_070423_4fc898_20260813_170015/` has C3 completing at n_events = 682
on the in-spec clock (statistic 1.511, p 0.0005) with no row timing out.
No growth exponent is claimed: three points do not determine a power law, and measuring the
cost for the n at hand is the only reliable answer.
The promotion report scores five checks. No conclusion anywhere rests on C3.

## 5. The bench's censoring arm never reached its label

The grid asked for censoring 0.25. The generator caps the segment count so each segment
expects at least a few events, and that cap binds at every n, so the realised value is
0.1675-0.1685 throughout. The out-of-envelope findings (notably C1 asymptotic at 0.699) are
real, but they happened at c ≈ 0.17. `jobs/bench/report.py` now prints the realised value.

## 6. The tie cutoff is weakly determined

`jobs/bench/grid.py` sweeps quantisation as a BOOLEAN, not as a distinct-duration count, so the
evidence brackets the cutoff between roughly 5 and 20 distinct values and no more finely.
The ledger declares 5 and says so. On the T2* ladder it never binds - durations are
wall-clock seconds and effectively continuous. It will bind on a quantised metric, and
there the number is a placeholder rather than a calibrated threshold.

## 7. What the provenance record cannot hold

The record has a closed schema. `alpha`, the ladder, the seed and every step kwarg reach
the Mermaid label, and the bench table's sha256 enters `dataset_hashes` and the run
identity. **What R produced a C3 p-value cannot**: it is discovered at runtime, so it lives
on the materialized `CheckLedger` and nowhere else - `r_version` (the R build and the
`copula` version), `r_executable` and `r_library_paths`. They separate three states: `not
asked` when the run-set excluded C3, `absent` when it was asked for and no `Rscript` was
found, and a real version when it was probed. Any probe failure against a present but
unusable R - a non-zero exit, a hang, a vanished interpreter, a silent one, or one
emitting bytes that are not valid UTF-8 (which the probes decode with `errors="replace"`,
so the catch there is defence in depth) - is recorded
on the artifact rather than stopping the ledger: a version-probe failure lands in
`r_version`, a library-probe failure in `r_library_paths`.

`c3_n_null_sim` sits beside them for a different reason. It does reach the label wherever a
job passes it as a step kwarg, but the C3 row's `notes` holds the verdict reason rather than
the simulation count, so the artifact is the only place a reader can recover the `N` a C3
p-value was drawn against.

## 8. Bench numbers describe synthetic records

Every size and power figure comes from generated data. Arms B and C carve their reads with
the real `analyzers/windows.py` primitives, so the carve policy is not idealised - but the
underlying processes are Weibull and Gaussian, not this instrument. The bench answers "what
would this check do on a record of this shape and size", not "what is this record".

The one place that gap is visible: after segments began splitting at read gaps, the real
record moved into an m > 1 regime that the bench measured only OUTSIDE its censoring
envelope. That regime has no in-envelope bench support.
