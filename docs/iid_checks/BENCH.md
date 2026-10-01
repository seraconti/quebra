# The calibration bench

`jobs/bench/` measures how often each check rejects when its null is true (size) and when it
is false (power), on data generated with known properties. It is a study, not a pipeline
layer: no pipeline package and not `cli.py` imports it (`tests/test_bench_isolation.py`;
tests may), and its results reach the pipeline only through
`jobs/bench/results/size_table.csv` and `power_table.csv`, which jobs declare as Datasets.
The numbers quoted below are from those tables as committed with this page.

## 1. What is simulated

A cell is one generating configuration at one event count on one clock
(`jobs/bench/grid.py:Cell`). Five arms generate the data (`jobs/bench/arms.py`):

| arm | generator | stands for |
|---|---|---|
| A, `arm_a` | iid Weibull durations, unit mean, truncated at a pre-chosen `tau` | the null of every check |
| B, `arm_b` | iid Gaussian reads, thresholded at the occupancy-matched quantile, carved | the null, with intrinsically quantised durations |
| C, `arm_c` | AR(1) Gaussian reads at correlation `rho`, carved identically | read-level dependence; its power cells are null cells, because read correlation does not reach the durations |
| D, `arm_d` | Weibull trend-renewal process, cumulative intensity `t^b` | the trend C1, C2 and CvM are directed at; `b = 1` is a second null |
| E, `arm_e` | Weibull durations with a Gaussian-copula AR(1) dependence of strength `rho` | the serial dependence C5 and C6 are directed at; `rho = 0` is Arm A |

The factors (`jobs/bench/grid.py`): event count `N_GRID` = 20, 35, 50, 75, 100, 355; Weibull
shape 0.75 and 1.5; quantised or continuous durations (Arm A); target censoring 0, 0.03 and
0.25 (Arms A and D, through `arms.n_segments_for_censoring`); `rho` 0.05-0.5 (Arms C and E);
`b` 0.7, 0.85, 1.2, 1.5 (Arm D). Arms B and C are carved on two clocks, and their in-spec
cells run only the rows that stay defined when `tau == T_N`: C5, C6 and CvM by
permutation (`grid._carved_clocks`, `checks/battery.py:run_battery`). `grid.size_cells` and
`grid.power_cells` give 144 and 192 cells.

The levels were set from one record family, measured before the bench was written
(`grid.py` module docstring, `REAL_DATA_NOTES` in `jobs/bench/report.py`): 20-355 usable
windows, censoring 0.000-0.026 wherever n >= 20, and its T2* quantisation. Nothing in the
code compares another record's event counts, censoring or quantisation with the grid.

## 2. What is measured

Each replicate is generated, the whole battery runs on it, and a row counts a rejection
when `p < 0.05` (`jobs/bench/runner.py:run_cell`). Size cells have 2000 replicates, power
cells 1000, and permutation rows use B = 999 (`grid.py`: `N_REPLICATES_SIZE`,
`N_REPLICATES_POWER`, `N_PERM`). The last run covered 336 cells and 480,000 replicates
(`jobs/bench/results/runtime.txt`).

The Monte Carlo standard error of a size at a true 0.05 is `sqrt(0.05 * 0.95 / 2000)` =
0.0049 for a size cell and 0.0069 at 1000 replicates. Every row carries its own `mc_se`, its
`n_used` replicates and the first failure reason, the realised censoring, and the
duration-level lag-1 the cell actually induced (`runner.py` module docstring).

## 3. How the numbers are consumed

**The ledger's gate, `analyzers/calibration_summary.py:bench_acceptance_at_n`.**
- It groups size rows by (check, calibration, variant, clock, n).
- It keeps only null-arm cells inside the censoring envelope, `ENVELOPE_MAX_CENSORING` =
  0.03, that rest on at least half their replicates (`MIN_SUPPORT_FRACTION`).
- It takes the WORST cell of each group and accepts it when `|z| <= bonferroni_z_crit(k)`,
  with `z` from the null standard error and `k` the group's cell count.

The groups hold 2, 8, 14 or 16 cells, so at 2000 replicates the size a group can show and
still be accepted is 0.05 +/- 0.011, 0.013, 0.014 and 0.014 respectively
(`bonferroni_z_crit(k) * null_se(0.05, 2000)`).

**`analyzers/check_ledger.py`.**
- `run` computes the acceptance table once.
- `_rows_for` judges each real record at `calibration_summary.nearest_bracketing_n`: the grid
  n nearest its event count, with ties going to the smaller.
- `_verdict` applies the outcome, in this order:
  - no p-value is `not computed`, and a tie-sensitive check on fewer than
    `tie_cutoff_distinct` (5) distinct durations is `not interpretable (ties)`;
  - a rejection is `fail`, with a note when the bench found the check miscalibrated at that
    n or has no cell for it;
  - a non-rejection is `pass` only if the record has at least `min_events_pass` events
    (35 by default, `CheckLedgerInputs`) AND the bench accepted the check at that n;
  - otherwise it is `underpowered`.

**The figures.** `jobs/active/check_calibration.py` declares both tables as Datasets and draws
`calibration_summary.size_vs_n`, `power_vs_dependence`, `read_dependence` and
`validation_curve`.

**The reports.**
- `jobs/bench/report.py:verdicts` sets the verdicts in `promotion_report.md`. It REJECTs a
  row whose size fails the same z-test over every in-envelope cell at n >= 35. It PROMOTEs a
  row whose mean power at n = 100 over its alternative's grid is at least 0.5, and HOLDs the
  rest.
- `analyzers/instrument_validation.py` cites size-table cells in the instrument report's
  tier-3 column.

## 4. Why measuring once and reusing is legitimate

Size and power are properties of a procedure at a given event count and null model, not of
a dataset. A check's rejection rate under iid Weibull durations at n = 50 is the same number
whichever record is later tested. So it is measured once, on generated data where the truth
is known, and read off for a real record by its event count and clock.

That reading is valid only as far as the real record's generating conditions sit inside
what was simulated. The limits below are mostly about that condition.

## 5. What it cannot do, and the proposed fix for each

None of these fixes is implemented. Each is a proposal.

1. **Staleness is undetected.** The tables record no commit, code hash, weighting or
   divisor. A change to a check leaves `bench_acceptance_at_n` reading numbers from code
   that no longer exists, and nothing fails. *Proposal:* write those four columns, and make
   `bench_acceptance_at_n` raise when they differ from the running code.
2. **Acceptance means "no detected miscalibration".** Its tolerance is a property of the
   replicate count (section 3), not of what calibration is required. *Proposal:* declare an
   equivalence margin in advance, for example `|size - 0.05| <= 0.01`, and test against it.
3. **The envelope was fitted to one record family.** The ledger does not check a new
   record's censoring, segment count or coefficient of variation against it. *Proposal:* an
   envelope check in the ledger that marks a record outside the simulated range.
4. **Nearest-n extrapolates.** A record below 20 events is judged at 20 and one above 355 at
   355. Between grid points the nearest is taken by distance, so 228-354 events are judged at
   355, a larger n than the record has. That is optimistic for asymptotic rows, whose size
   improves with n. *Proposal:* judge at the largest grid n not above the record's count, and
   report records below 20 as outside the grid.
5. **Permutation rows are exact by construction.** Within-segment permutation of
   exchangeable gaps has the nominal size up to the `p < 0.05` grid (49/1000 at B = 999), so
   the bench tests only their implementation. No fix proposed: that is what such a row can
   show.
6. **Power never enters a ledger verdict.** `_verdict` reads the size acceptance only, so
   `pass` is a non-rejection by a calibrated check, not evidence of no trend or no
   dependence. Averaged over its dependence grid, C5 has 36% (studentized) and 35%
   (unstudentized) power at n = 100
   (`promotion_report.md`, verdict table, `power_n100_mean`). *Proposal:* carry the
   power at the record's n into the ledger row.
7. **A transcription error shared by bench and production is invisible.** The bench runs the
   same check code it calibrates, so a wrong equation is wrong in both and its size can
   still look right. Only published values catch it (`tests/test_checks_published_values.py`,
   `docs/GOLD_STANDARD.md`). No fix proposed.
8. **Multiplicity across groups is uncorrected.** The Bonferroni is per group, over 108
   groups, so about five false "miscalibrated" flags are expected even for exact tests. Two
   exact permutation rows are flagged today: C5 unstudentized, in-spec, n = 20, and CvM
   permutation, calendar, n = 100 (`bench_acceptance_at_n` computed from `size_table.csv`).
   *Proposal:* correct across groups, or exempt permutation rows from the size gate.
9. **Ties are under-counted on quantised data.** `_permutation.permutation_p_value` counts
   ties with a float-exact `>=`. On a read grid, shuffles that tie the observed statistic
   exactly can differ in the last bits and are missed, so p is slightly low. *Proposal:* a
   relative tolerance, itself to be calibrated, and a flag on the row when the share of
   near-ties is large.
10. **The tau-free retry never succeeds.** `run_cell` retries a failed replicate without the
    tau checks, but CvM's permutation row stays in that reduced battery and raises on the
    same degenerate segment. `n_failed_tau_only` is 0 in every row of both tables, so C5 and
    C6 lose those replicates too. All 12 cells that lose replicates this way are at
    censoring 0.25, outside the envelope; inside it, five n = 20 cells (three size, two Arm E
    power) lost one replicate each to a segment too short for any lag, which no retry can
    save. *Proposal:* drop CvM
    from the retry.
11. **Some of the null is not exercised.**
    - The carved arms have no gaps, so the carve's gap handling is never run under a null.
    - At n <= 35 a censoring target of 0.03 yields one segment, the same as 0
      (`n_segments_for_censoring`).
    - On the calendar clock the bench passes an observation end to the carve
      (`arms._carved_segments`), and the ledger passes none (`check_ledger._rows_for`).

    *Proposal:* a gapped carved arm, and the same observation end on both sides.
12. **The rule differs by one boundary.** The bench counts `p < 0.05` and the ledger rejects
    on `p <= alpha`, about +1/1000 on permutation rows. *Proposal:* use one rule.
13. **Nothing checks the table's alpha.** `bench_acceptance_at_n` z-scores the measured
    `P(p < 0.05)` against whatever `alpha` it is given. *Proposal:* raise when the table's
    `alpha` column differs.
14. **Seeds depend on grid order.** Each replicate is seeded from its cell's position in
    `grid.all_cells()`, which no table records, so reordering the grid silently re-seeds
    every later cell. *Proposal:* seed from the cell's key.
15. **The dependence grid is not calibrated to a record.** `rho` sets the strength of the
    generated dependence; no grid point is claimed to be a real record's. Power is reported
    over the whole grid, with the lag-1 each `rho` induced (`report.py:dependence_grid_table`).
    No fix proposed.
