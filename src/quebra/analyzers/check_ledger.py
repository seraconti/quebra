"""The check ledger: every check, on every threshold, with what its answer is worth.

A p-value alone is not a result here. The same 0.4 means "no evidence of dependence" from
a calibrated check on 355 windows and "this test cannot see anything" from the same check
on 20. So every row carries three things beside the p-value, and a `pass` requires all
three to hold:

    n_events            at least min_events_pass events (a policy floor, not a power
                        threshold)
    bench_size_at_n     the check was calibrated at that event count, per the bench
    tie_fraction        the durations are not so tied that a rank test is meaningless

If any fails, the verdict is NOT `pass`. That is the entire point of the ledger: a
non-rejection must never print as a pass on a p-value alone. A group row (below) is the one
exception to the bench condition.

**Group rows.** A window ledger also carries two rows per (threshold, side, clock) that are
no single check: `trend_group` over C1, C2 and CvM and `dependence_group` over C5 and C6, each
the single-step min-P (`min_p`) of its members that answered, on their shared permutation
draws. A member left out is named with its reason, and only a refusal of the record
(`checks.result.CheckRefused`) may leave one out: any other failure leaves the group `not
computed`, since a member set that could change with the shuffle would void the size. The
group's size rests on the permutation argument, exact under exchangeability, so it reads no
bench cell; it keeps the event floor, and the dependence group the tie rule of its members.
The trend group counts the events its members read, `n - 1` on a segment that ends at one.
Calibration cycles score through `score_segments` without them.

**Which bag, which clock.** Both bags, the in-spec and the out-of-spec windows, on both
clocks: a row is one (threshold, side, clock, check or group). The bags are scored by one rule
(`_multiprocess.segments_from_windows`), draw their permutations from one stream per
threshold and clock (`stream_for`), and read one bench cell: the bench calibrated the side
clock on arms that either carry no side or threshold a symmetric series at its median, and
the calendar clock on the median arms alone (B and C), so its cells do not depend on the
bag. That holds at the median only. At a rung where a fraction p of reads is in spec, the
out-of-spec bag has the law of the in-spec bag at 1 - p, which no cell measured either: the
cells transfer to both bags as far as they transfer to one, an open question for both. The
shared stream means the bags' p-values share their draws wherever their segment sizes
match, so they are not independent.

**Segments that end at an event.** On the side clock of a carved record, time on the side
stops accruing the moment the record ends on the other side, so `tau == T_N` (measured at
~73% of synthetic replicates on the in-spec bag). The same holds where the final window is
one read long, and on every calibration-cycles segment that ends with a complete cycle.
That is failure truncation: given `T_n`, under a constant rate the first `n - 1` event times
are uniform order statistics on `(0, T_n)`, so C1, C2 and CvM read them with `tau = T_n`
(Rigdon and Basu 2000, Sect. 4.3; `Segment.ends_at_event`). On all `N` events Kvaloy and
Lindqvist's eq (7) (Technometrics 2020, p. 104) would also be singular there. Their rows
count, judge and look up the bench by the events read, and say how many segments ended at an
event. No bench cell measured the `n - 1` form: C1's and C2's cells are all time truncated,
and CvM's side-clock cells include the carved arms' (B and C), which read all `n` events
where a segment ended at an event. The rank checks never touch `tau` and read every event.
The ledger carries no asymptotic row.

**What the provenance record cannot hold.** Its schema is closed, so what R produced the
C3 p-values is discovered at runtime and lives on the materialized artifact, not in the
`.prov.json`: `r_version` (the R build and the `copula` version), `r_executable` and
`r_library_paths`. With no `Rscript` on PATH the two strings read `absent` and the tuple is
empty. `alpha`, the ladder, the seed and the thresholds DO reach the label, because they
are step kwargs.

`c3_n_null_sim` is on the artifact for a different reason. It does reach the label wherever
a job passes it as a step kwarg, but a C3 p-value simulated against 200 draws is not the
same object as one simulated against 1000, and a ledger has to be readable without the job
that produced it.
"""

from __future__ import annotations

import subprocess
import zlib
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

import quebra.analyzers.checks.c1_lewis_robinson as c1
import quebra.analyzers.checks.c2_anderson_darling as c2
import quebra.analyzers.checks.c3_serial_copula as c3
import quebra.analyzers.checks.cvm_cramer_von_mises as cvm
from quebra.analyzers.calibration_summary import (
    bench_acceptance_at_n,
    nearest_bracketing_n,
)
from quebra.analyzers.checks._multiprocess import segments_from_windows
from quebra.analyzers.checks._permutation import block_permutations
from quebra.analyzers.checks._rank_serial import MAX_LAG_CAP
from quebra.analyzers.checks.battery import (
    REFUSED,
    ROW_KEYS,
    UNDEFINED,
    drop_unusable_segments,
    row_key,
    run_battery,
)
from quebra.analyzers.checks.result import (
    BENCH_CLOCK,
    CALIB_R_COPULA,
    CLOCKS,
    CheckResult,
    Segment,
)
from quebra.analyzers.event_table import SIDES
from quebra.analyzers.min_p import min_p
from quebra.core._artifact_guard import StaleArtifactGuard

VERDICT_PASS = "pass"
VERDICT_FAIL = "fail"
VERDICT_INCONCLUSIVE = "inconclusive"
VERDICT_TIES = "not interpretable (ties)"
VERDICT_NOT_COMPUTED = "not computed"
# Not produced by `_verdict` and deliberately NOT in `VERDICTS` below: the ledger never emits
# it. It is what a downstream builder writes when it has no row for a cell at all, and it
# lives here so the render contract has one vocabulary rather than two - `theme.verdict_color`
# raises on anything outside its palette, so "no row" is a verdict as far as drawing is
# concerned. Keeping it out of `VERDICTS` keeps `panels/check_ledger.py`'s legend honest,
# since that legend lists what the LEDGER can say.
VERDICT_ABSENT = "no row"
VERDICTS = (
    VERDICT_PASS,
    VERDICT_FAIL,
    VERDICT_INCONCLUSIVE,
    VERDICT_TIES,
    VERDICT_NOT_COMPUTED,
)

# HOW OFTEN THIS RULE BINDS, measured rather than assumed on the 1 to 10 µs T2* ladder
# (commit aff894e; not re-measured on the 2 to 5 µs ladder): it never did except where the
# event count already disqualified the row. T2* durations are wall-clock
# seconds and effectively continuous - 355 windows gave 355 distinct values - so ties are a
# property of QUANTISED metrics (a fidelity threshold where most windows are one read long),
# not of this one. The rule earns its place on those ladders, not here.
#
# Checks whose validity rests on the durations being effectively continuous: C5 and C6 are
# rank statistics, C3's distribution-free property is derived for continuous observations.
# C1 and C2 use the durations directly and degrade differently (a fully tied vector makes
# gamma_hat zero, which raises rather than misleads), so the tie rule does not apply to
# them.
TIE_SENSITIVE_CHECKS = ("c3_serial_copula", "c5_rank_autocorr", "c6_exchangeability")

# The checks that read event times against a segment's `tau`, so read only the first n - 1
# events of a segment that ends at one (`Segment.ends_at_event`). Their rows count, judge
# and look up the bench by the events they read.
EVENT_TIME_CHECKS = (c1.CHECK_NAME, c2.CHECK_NAME, cvm.CHECK_NAME)

# The two groups a window ledger calibrates by min-P, each the permutation rows of the checks
# that test one null (module docstring). A group row is keyed like a check row.
GROUP_TREND = "trend_group"
GROUP_DEPENDENCE = "dependence_group"
CALIB_MIN_P = "min_p"
_PERMUTATION_ROW = {key[0]: key for key in ROW_KEYS}
GROUPS: dict[str, tuple[tuple[str, str, str], ...]] = {
    GROUP_TREND: tuple(
        _PERMUTATION_ROW[check]
        for check in (
            "c1_lewis_robinson",
            "c2_anderson_darling",
            "cvm_cramer_von_mises",
        )
    ),
    GROUP_DEPENDENCE: tuple(
        _PERMUTATION_ROW[check] for check in ("c5_rank_autocorr", "c6_exchangeability")
    ),
}
# The checks `include_tau_checks=False` leaves out, which a group then reads `not asked`.
_TAU_CHECKS = ("c1_lewis_robinson", "c2_anderson_darling")
# How a member reads in a group row's notes.
MEMBER_SHORT = {
    "c1_lewis_robinson": "C1",
    "c2_anderson_darling": "C2",
    "cvm_cramer_von_mises": "CvM",
    "c5_rank_autocorr": "C5",
    "c6_exchangeability": "C6",
}
GROUP_KEYS = tuple((group, CALIB_MIN_P, "") for group in GROUPS)
# Members whose p reads |statistic| (`_permutation.two_sided_p_value`): min-P ranks them so.
TWO_SIDED_MEMBERS = frozenset({"c1_lewis_robinson"})
# What every figure that draws a group row says of its calibration.
GROUP_CALIBRATION_NOTE = (
    "group rows: size from the permutation argument, exact under exchangeability; no "
    "bench cell is used."
)

# The side of a row whose durations are no bag's windows (calibration cycles).
NO_SIDE = ""

LEDGER_COLUMNS = [
    "dataset_id",
    "threshold_label",
    "side",
    "clock",
    "check_id",
    "calibration",
    "variant",
    "statistic",
    "p_value",
    "verdict",
    "n_events",
    "n_windows_used",
    "n_censored_dropped",
    "n_segments",
    "n_distinct_durations",
    "tie_fraction",
    "bench_n_target",
    "bench_size_at_n",
    "bench_size_ci",
    "bench_size_cell",
    "bench_accepted",
    "notes",
]


@dataclass
class CheckLedger(StaleArtifactGuard):
    """One row per (threshold, side, clock, check or group, calibration), plus the settings
    that made it.

    The settings are ON the artifact rather than only in the job file because a ledger read
    six months from now has to be interpretable without opening the job that produced it -
    `alpha` and `min_events_pass` are what turn a p-value into a verdict.
    """

    rows: pd.DataFrame
    dataset_id: str = ""
    alpha: float = 0.05
    min_events_pass: int = 0
    tie_cutoff_distinct: int = 0
    lag_max: int = MAX_LAG_CAP
    n_permutations: int = 0
    seed: int = 0
    # Runtime-discovered, and the provenance schema has no free-form field to hold them.
    # `not asked` when the run-set excluded C3, `absent` when it was asked for and no
    # interpreter was found. `r_library_paths` is R's own report of its search path, so it
    # is only as trustworthy as the profile that produced it.
    r_version: str = "not asked"
    r_executable: str = "not asked"
    r_library_paths: tuple[str, ...] = ()
    # The value the inputs carried, recorded whether or not C3 ran. Not inferable from the
    # rows: the verdict reason replaces the note that carried it.
    c3_n_null_sim: int = c3.N_NULL_SIM
    thresholds: list[tuple[str, float, bool]] = field(default_factory=list)

    def check_thresholds(self, labels: list[str]) -> None:
        """Every threshold on the ladder appears in the ledger, or the artifact is partial."""
        present = set(self.rows["threshold_label"]) if len(self.rows) else set()
        missing = sorted(set(labels) - present)
        if missing:
            raise ValueError(
                f"incomplete CheckLedger: no rows for threshold(s) {missing} - construct "
                "via analyzers.check_ledger.run()"
            )

    def passing(self) -> pd.DataFrame:
        return self.rows[self.rows["verdict"] == VERDICT_PASS]

    def naive_passes(self) -> pd.DataFrame:
        """Rows a p-value-only rule would have called `pass`, and this rule does not.

        This is the ledger's own demonstration that the extra two conditions do work.
        """
        naive = self.rows["p_value"].notna() & (self.rows["p_value"] > self.alpha)
        return self.rows[naive & (self.rows["verdict"] != VERDICT_PASS)]


@dataclass(slots=True)
class CheckLedgerInputs:
    # Both bags' windows, one table, each row with its `side` (`make_inputs_from_windows`).
    # A bag with no rows at a threshold reads `no windows at this threshold` there.
    windows: pd.DataFrame
    bench_size_table: pd.DataFrame
    thresholds: list[tuple[str, float, bool]]
    dataset_id: str = ""
    gap_spans_s: list[tuple[float, float]] | None = None
    alpha: float = 0.05
    min_events_pass: int = 35
    tie_cutoff_distinct: int = 5
    lag_max: int = MAX_LAG_CAP
    n_permutations: int = 999
    seed: int = 0
    min_events_per_segment: int = 2
    # C3 is the only out-of-process check and by far the most expensive: about a minute at
    # n = 355 (see c3_serial_copula for the measured range and its machine). A ledger does
    # not run it unless asked: C3 has no bench cell, so its rows are `not computed` or an
    # uncalibrated p-value either way. When False the C3 rows are OMITTED rather than
    # written as `not computed`: a blank row would claim the check was attempted and
    # failed, when in fact it was never asked.
    include_c3: bool = False
    # Derived from a run-set by `check_outcome.battery_flags`, never set by hand once one is
    # in play. COARSE on purpose: `include_tau_checks` is all-or-nothing over C1's and C2's
    # permutation rows, so a run-set naming one still computes the other. The declared
    # run-set therefore names AT MOST what ran, not exactly what ran, and a reader comparing
    # the two will find the ledger carrying rows the run-set does not list.
    include_tau_checks: bool = True
    include_c2_asymptotic: bool = True
    # See c3_serial_copula.N_NULL_SIM: the dominant cost, and part of what a C3
    # p-value means, so it is carried rather than left to a module default.
    c3_n_null_sim: int = c3.N_NULL_SIM
    # Set when the record failed to load or carve (`core.types.RecordFailure.note`): every
    # threshold and clock then reads `not computed` with this note, so the record is in
    # the ledger and the survey rather than missing from them.
    record_failure: str = ""


def make_inputs_from_windows(
    window_result: object,
    bench_size_table: pd.DataFrame,
    *,
    thresholds: list[tuple[str, float, bool]],
    alpha: float,
    min_events_pass: int,
    tie_cutoff_distinct: int,
    lag_max: int,
    n_permutations: int,
    seed: int,
    include_c3: bool = False,
    include_tau_checks: bool = True,
    include_c2_asymptotic: bool = True,
    c3_n_null_sim: int = c3.N_NULL_SIM,
) -> CheckLedgerInputs:
    """Build inputs from a `WindowsResult`: both bags' windows, and the gap spans.

    `gap_spans_s` is not optional in spirit: without it `segments_from_windows` falls back
    to reading birth types alone, which cannot see a gap flanked by out-of-spec reads and
    silently merges two renewal processes separated by unobserved hours.
    """
    diagnostics = getattr(window_result, "diagnostics", {}) or {}
    return CheckLedgerInputs(
        windows=pd.concat(
            [window_result.windows_in_spec, window_result.windows_out_of_spec],
            ignore_index=True,
        ),
        bench_size_table=bench_size_table,
        thresholds=list(thresholds),
        dataset_id=str(getattr(window_result, "meta", {}).get("dataset_id", "")),
        gap_spans_s=diagnostics.get("gap_spans_s"),
        alpha=alpha,
        min_events_pass=min_events_pass,
        tie_cutoff_distinct=tie_cutoff_distinct,
        lag_max=lag_max,
        n_permutations=n_permutations,
        seed=seed,
        include_c3=include_c3,
        include_tau_checks=include_tau_checks,
        include_c2_asymptotic=include_c2_asymptotic,
        c3_n_null_sim=c3_n_null_sim,
    )


def stream_for(label: str, clock: str) -> int:
    """The permutation rng stream for one (threshold label, clock) cell.

    crc32, not `hash()`: Python randomises string hashing per process (PYTHONHASHSEED), so
    `hash(label)` would make every permutation p-value differ between runs while the
    provenance record stayed identical - the exact defect the explicit-rng rule exists to
    prevent. crc32 is stable across processes and versions.

    Public rather than private because the cross-process stability claim above is only
    testable if a test can call the derivation the shipped code uses. Inlined, it was
    guarded by a test that re-typed the formula as a string literal and therefore could not
    fail when this line changed.
    """
    return zlib.crc32(f"{label}|{clock}".encode()) & 0x7FFFFFFF


def _tie_stats(durations: np.ndarray) -> tuple[int, float]:
    """Distinct value count and the fraction of observations sharing a value."""
    if not len(durations):
        return 0, float("nan")
    _values, counts = np.unique(durations, return_counts=True)
    tied = float(counts[counts > 1].sum()) / float(len(durations))
    return int(len(_values)), tied


def _verdict(
    result: CheckResult,
    *,
    alpha: float,
    n_events: int,
    min_events_pass: int,
    n_distinct: int,
    tie_cutoff: int,
    bench_accepted: bool | None,
) -> tuple[str, str]:
    """The verdict and the reason it is not `pass`, in precedence order.

    Neither `pass` nor `inconclusive` measures power. How much a `pass` excludes is read
    from the bench's power table (`jobs/bench/results/power_table.csv`) at the row's
    `bench_n_target`, the nearest grid n, and only against the bench's own alternatives.
    """
    if result.p_value is None:
        return VERDICT_NOT_COMPUTED, result.notes
    if result.check in TIE_SENSITIVE_CHECKS and n_distinct < tie_cutoff:
        return (
            VERDICT_TIES,
            f"{n_distinct} distinct durations, below the cutoff of {tie_cutoff}",
        )
    if result.p_value <= alpha:
        # A rejection. Whether it is trustworthy still depends on the calibration, so the
        # note carries that even when the verdict does not.
        #
        # BOTH calibration states are named here, not just one. The non-rejection branch
        # below demotes `bench_accepted is None` to `inconclusive` because a non-rejection
        # from a check of unknown size is not evidence for the null. The same is true of a
        # rejection: with no bench cell there is no measured level, so `p <= alpha` does not
        # mean the test rejects at alpha. The verdict stays `fail` because the statistic did
        # land in the tail and suppressing that would hide a finding, but the note must say
        # the level behind it is unmeasured.
        note = "rejected"
        if bench_accepted is False:
            note += "; but the bench found this check miscalibrated at this event count"
        elif bench_accepted is None:
            note += (
                "; UNCALIBRATED: no bench cell for this check at this event count, so no "
                "level stands behind this rejection"
            )
        return VERDICT_FAIL, note
    unmet = []
    if n_events < min_events_pass:
        unmet.append(f"n_events {n_events} < {min_events_pass}")
    if bench_accepted is False:
        unmet.append("bench found this check miscalibrated at this event count")
    if bench_accepted is None:
        unmet.append("no bench cell for this check at this event count")
    if unmet:
        # A non-rejection that cannot be read as evidence for the null.
        return VERDICT_INCONCLUSIVE, "; ".join(unmet)
    return VERDICT_PASS, ""


def _r_provenance(asked: bool) -> tuple[str, str, tuple[str, ...]]:
    """What R produced the C3 rows: `(version, executable, library paths)`.

    THREE STATES, not two. `not asked` when the run-set excluded C3, so no R was consulted
    and none was needed; `absent` when C3 was asked for and no interpreter was found; and a
    real version string when it was asked for and probed. Collapsing the first two would put
    a complete R provenance block on a ledger that never ran a C3 row, which reads as
    evidence about work that did not happen.

    Probed once per ledger, above the threshold loop: inside it, the two subprocesses would
    be paid again on every row of every threshold and clock.

    NOTHING HERE RAISES, and the except clause is wide because the failure modes are the
    environment's, not the data's. A hung interpreter reaches the caller as
    `subprocess.TimeoutExpired`, a deleted one as `OSError`, and one emitting bytes that are
    not valid UTF-8 as a `ValueError`. That last one is defence in depth rather than a live
    path: both probes pass `errors="replace"`, so the decode cannot raise today, and the
    catch is what keeps that true if the decoding changes. None of the three is a
    `RuntimeError`, and any of them
    escaping would kill the ledger above the first row, which is the one thing a check outcome
    may never do. `check_ledger`'s per-row C3 clause already lists `ValueError` for the same
    reason. The reason is recorded in `r_version` instead.
    """
    if not asked:
        return "not asked", "not asked", ()
    executable = c3.rscript_path()
    if executable is None:
        return "absent", "absent", ()
    try:
        libraries = c3.r_library_paths() or ("probe returned nothing",)
    except (RuntimeError, ValueError, subprocess.TimeoutExpired, OSError) as exc:
        # A visible sentinel, not `()`. An empty tuple is also what `not asked` and `absent`
        # carry, so a silent fallback here would make a failed probe indistinguishable from
        # a run that never made one.
        libraries = (f"probe failed ({type(exc).__name__})",)
    try:
        version = c3.r_version()
    except (RuntimeError, ValueError, subprocess.TimeoutExpired, OSError) as exc:
        return f"probe failed ({type(exc).__name__})", str(executable), libraries
    return "absent" if version is None else version, str(executable), libraries


def run(inputs: CheckLedgerInputs) -> CheckLedger:
    """Run every check on every threshold and score each answer."""
    acceptance = bench_acceptance_at_n(inputs.bench_size_table, alpha=inputs.alpha)
    grid = sorted({int(n) for n in acceptance["n_target"]})
    r_version, r_executable, r_libraries = _r_provenance(inputs.include_c3)

    rows: list[dict[str, object]] = []
    if not inputs.record_failure and len(inputs.windows):
        if "side" not in inputs.windows.columns:
            raise KeyError("the ledger's window table has no side column")
        unknown = sorted({str(v) for v in inputs.windows["side"]} - set(SIDES))
        if unknown:
            raise ValueError(f"window rows of no bag: side {unknown}; bags: {SIDES}")
    for label, _value, _big_good in inputs.thresholds:
        for side in SIDES:
            for clock in CLOCKS:
                if inputs.record_failure:
                    rows.append(
                        _blank_row(
                            label,
                            side,
                            clock,
                            inputs.dataset_id,
                            inputs.record_failure,
                        )
                    )
                    continue
                table = inputs.windows[
                    (inputs.windows["threshold_label"] == label)
                    & (inputs.windows["side"] == side)
                ]
                rows.extend(
                    _rows_for(
                        table=table,
                        label=label,
                        side=side,
                        clock=clock,
                        inputs=inputs,
                        acceptance=acceptance,
                        grid=grid,
                    )
                )
    frame = pd.DataFrame(rows).reindex(columns=LEDGER_COLUMNS)
    ledger = CheckLedger(
        rows=frame,
        dataset_id=inputs.dataset_id,
        alpha=inputs.alpha,
        min_events_pass=inputs.min_events_pass,
        tie_cutoff_distinct=inputs.tie_cutoff_distinct,
        lag_max=inputs.lag_max,
        n_permutations=inputs.n_permutations,
        seed=inputs.seed,
        r_version=r_version,
        r_executable=r_executable,
        r_library_paths=r_libraries,
        c3_n_null_sim=inputs.c3_n_null_sim,
        thresholds=list(inputs.thresholds),
    )
    ledger.check_thresholds([label for label, _v, _b in inputs.thresholds])
    return ledger


def _blank_row(
    label: str,
    side: str,
    clock: str,
    dataset_id: str,
    note: str,
    n_windows: int = 0,
) -> dict:
    return {
        "dataset_id": dataset_id,
        "threshold_label": label,
        "side": side,
        "clock": clock,
        "check_id": "(all)",
        "calibration": "",
        "variant": "",
        "statistic": float("nan"),
        "p_value": None,
        "verdict": VERDICT_NOT_COMPUTED,
        "n_events": 0,
        "n_windows_used": int(n_windows),
        "n_censored_dropped": 0,
        "n_segments": 0,
        "n_distinct_durations": 0,
        "tie_fraction": float("nan"),
        "bench_n_target": None,
        "bench_size_at_n": float("nan"),
        "bench_size_ci": "",
        "bench_size_cell": "",
        "bench_accepted": None,
        "notes": note,
    }


def _rows_for(
    *,
    table: pd.DataFrame,
    label: str,
    side: str,
    clock: str,
    inputs: CheckLedgerInputs,
    acceptance: pd.DataFrame,
    grid: list[int],
) -> list[dict[str, object]]:
    if not len(table):
        return [
            _blank_row(
                label, side, clock, inputs.dataset_id, "no windows at this threshold"
            )
        ]
    try:
        segments, n_dropped = segments_from_windows(
            table,
            clock=clock,
            side=side,
            min_events=inputs.min_events_per_segment,
            gap_spans_s=inputs.gap_spans_s,
        )
    except (ValueError, KeyError) as exc:
        return [
            _blank_row(
                label,
                side,
                clock,
                inputs.dataset_id,
                f"segmentation failed: {exc}"[:200],
                len(table),
            )
        ]
    # `segments_from_windows` already dropped the short ones, by the same rule.
    segments, n_short, n_degenerate = drop_unusable_segments(
        segments, inputs.min_events_per_segment
    )
    return score_segments(
        segments,
        n_dropped=n_dropped + n_short,
        n_degenerate=n_degenerate,
        n_used=len(table),
        label=label,
        side=side,
        clock=clock,
        dataset_id=inputs.dataset_id,
        acceptance=acceptance,
        grid=grid,
        alpha=inputs.alpha,
        min_events_pass=inputs.min_events_pass,
        tie_cutoff_distinct=inputs.tie_cutoff_distinct,
        lag_max=inputs.lag_max,
        n_permutations=inputs.n_permutations,
        seed=inputs.seed,
        min_events_per_segment=inputs.min_events_per_segment,
        include_tau_checks=inputs.include_tau_checks,
        include_c2_asymptotic=inputs.include_c2_asymptotic,
        include_c3=inputs.include_c3,
        c3_n_null_sim=inputs.c3_n_null_sim,
        groups=True,
    )


def score_segments(
    segments: list[Segment],
    *,
    n_dropped: int,
    n_degenerate: int,
    n_used: int,
    label: str,
    clock: str,
    dataset_id: str,
    acceptance: pd.DataFrame,
    grid: list[int],
    alpha: float,
    min_events_pass: int,
    tie_cutoff_distinct: int,
    lag_max: int,
    n_permutations: int,
    seed: int,
    min_events_per_segment: int,
    include_tau_checks: bool,
    include_c2_asymptotic: bool,
    include_c3: bool,
    c3_n_null_sim: int,
    side: str = NO_SIDE,
    groups: bool = False,
) -> list[dict[str, object]]:
    """Ledger rows for durations already cut into segments, on one label, side and clock.

    `side` is the bag the segments came from, `NO_SIDE` for durations that are no bag's
    windows. The permutation stream and the bench cell are looked up under the clock's
    bench name (`checks.result.BENCH_CLOCK`) and do not depend on the side.

    The half of the ledger that does not care where the segments came from: the battery
    off one permutation set, C3, the bench match and the verdict. `_rows_for` reaches it
    after carving windows into segments; any other duration (calibration cycles) reaches it
    with its own segments and the same rule. `n_dropped` and `n_degenerate` count the
    segments the caller dropped through `battery.drop_unusable_segments`, as short and as
    degenerate, and `n_used` is what the `n_windows_used` column records. A nonzero
    `n_degenerate` is named in every row's notes. A check undefined on the segments that
    remain reads `not computed` with its own reason, and the others still answer, C1 and C2
    included. With no segment left it returns one `not computed` row with its reason and
    never raises.

    C1, C2 and CvM read `n - 1` events of a segment that ends at an event
    (`EVENT_TIME_CHECKS`): their rows' `n_events`, event floor and bench n count the events
    read, and their notes say how many segments ended at an event.

    `groups=True` appends the two group rows (`group_rows`): the window ledgers ask for
    them, and calibration cycles do not.
    """
    if not segments:
        note = (
            f"no segment left: {n_dropped} below {min_events_per_segment} events, "
            f"{n_degenerate} with every duration identical"
            if n_degenerate
            else f"no segment reached {min_events_per_segment} events "
            f"({n_dropped} dropped)"
        )
        return [_blank_row(label, side, clock, dataset_id, note, n_used)]
    degenerate_note = (
        f"{n_degenerate} segment(s) with every duration identical dropped"
        if n_degenerate
        else ""
    )

    durations = np.concatenate([s.x for s in segments])
    n_distinct, tie_fraction = _tie_stats(durations)
    n_events = int(len(durations))
    n_read = int(sum(s.n_events_read for s in segments))
    n_at_event = sum(s.ends_at_event for s in segments)

    stream = stream_for(label, BENCH_CLOCK[clock])
    rng = np.random.default_rng([seed, stream])
    perm = block_permutations([s.n_events for s in segments], n_permutations, rng)

    try:
        results = run_battery(
            segments,
            clock=clock,
            perm=perm,
            max_lag=lag_max,
            include_tau_checks=include_tau_checks,
            include_c2_asymptotic=include_c2_asymptotic,
            undefined_as_rows=True,
        )
    except (ValueError, KeyError) as exc:
        return [
            _blank_row(
                label, side, clock, dataset_id, f"battery failed: {exc}"[:200], n_used
            )
        ]

    # C3 is out-of-process and shares nothing with the permutation set, so it is called
    # separately. With R absent it returns p_value=None and the row reads `not computed`.
    if include_c3:
        try:
            results.append(
                c3.run(
                    segments,
                    clock=clock,
                    max_lag=lag_max,
                    seed=seed,
                    n_null_sim=c3_n_null_sim,
                )
            )
        # RuntimeError and TimeoutExpired are the R bridge's own failure modes - a
        # non-zero Rscript exit, a missing or unparseable result file, or a simulation
        # that outran its timeout. Neither is caught by the row count, and it does not show
        # while R was absent because `rscript_path()` returned None and the bridge never
        # ran. With R installed the call executes for real, and an uncaught RuntimeError
        # kills the whole ledger job instead of writing the `not computed` row this
        # except clause exists to write.
        except (ValueError, KeyError, RuntimeError, subprocess.TimeoutExpired) as exc:
            results.append(
                CheckResult(
                    check=c3.CHECK_NAME,
                    statistic=float("nan"),
                    p_value=None,
                    calibration=CALIB_R_COPULA,
                    clock=clock,
                    n_events=n_events,
                    n_segments=len(segments),
                    notes=f"declined: {exc}"[:160],
                )
            )

    # What a row that reads `n - 1` events says of it (`EVENT_TIME_CHECKS`).
    read_note = (
        f"{n_at_event} of {len(segments)} segment(s) end at an event: first n - 1 events "
        "read"
        if n_at_event
        else ""
    )
    out: list[dict[str, object]] = []
    for result in results:
        variant = ""
        for token in result.notes.split():
            if token.startswith("variant="):
                variant = token.split("=", 1)[1]
        reads_n_minus_1 = result.check in EVENT_TIME_CHECKS
        n_row = n_read if reads_n_minus_1 else n_events
        bench_n = nearest_bracketing_n(grid, n_row)
        match = acceptance[
            (acceptance["check"] == result.check)
            & (acceptance["calibration"] == result.calibration)
            & (acceptance["variant"].fillna("") == variant)
            & (acceptance["clock"] == BENCH_CLOCK[clock])
            & (acceptance["n_target"] == bench_n)
        ]
        if len(match):
            entry = match.iloc[0]
            accepted = bool(entry["bench_accepted"])
            size = float(entry["bench_size_at_n"])
            half = 1.96 * float(entry["bench_size_se"])
            ci = f"[{size - half:.4f}, {size + half:.4f}]"
            cell = str(entry["bench_size_cell"])
        else:
            accepted, size, ci, cell = None, float("nan"), "", ""

        verdict, reason = _verdict(
            result,
            alpha=alpha,
            n_events=n_row,
            min_events_pass=min_events_pass,
            n_distinct=n_distinct,
            tie_cutoff=tie_cutoff_distinct,
            bench_accepted=accepted,
        )
        if UNDEFINED in result.extra:
            reason = str(result.extra[UNDEFINED])
        # Only these rows read n - 1, so only they say so; the bench clause only where a
        # cell was read.
        at_event_note = read_note if reads_n_minus_1 else ""
        if at_event_note and cell:
            at_event_note += "; the bench cell was measured without the n - 1 form"
        notes = "; ".join(
            part for part in (reason, at_event_note, degenerate_note) if part
        )
        out.append(
            {
                "dataset_id": dataset_id,
                "threshold_label": label,
                "side": side,
                "clock": clock,
                "check_id": result.check,
                "calibration": result.calibration,
                "variant": variant,
                "statistic": float(result.statistic),
                "p_value": result.p_value,
                "verdict": verdict,
                "n_events": n_row,
                "n_windows_used": int(n_used),
                "n_censored_dropped": int(result.n_censored_dropped),
                "n_segments": int(result.n_segments),
                "n_distinct_durations": n_distinct,
                "tie_fraction": tie_fraction,
                "bench_n_target": bench_n,
                "bench_size_at_n": size,
                "bench_size_ci": ci,
                "bench_size_cell": cell,
                "bench_accepted": accepted,
                "notes": notes,
            }
        )
    if groups:
        out.extend(
            group_rows(
                results,
                dataset_id=dataset_id,
                label=label,
                side=side,
                clock=clock,
                alpha=alpha,
                n_events=n_events,
                n_read=n_read,
                read_note=read_note,
                n_used=n_used,
                min_events_pass=min_events_pass,
                n_distinct=n_distinct,
                tie_fraction=tie_fraction,
                tie_cutoff=tie_cutoff_distinct,
                include_tau_checks=include_tau_checks,
                n_segments=len(segments),
                n_censored_dropped=int(sum(s.n_censored_dropped for s in segments)),
                degenerate_note=degenerate_note,
            )
        )
    return out


def _first_clause(text: str, limit: int = 90) -> str:
    """A reason cut to its first clause, so the member names lead a group row's notes."""
    for end in (". ", "; "):
        if end in text:
            text = text.split(end, 1)[0]
    return text if len(text) <= limit else text[: limit - 3] + "..."


def _group_verdict(
    p_value: float | None,
    reason: str,
    *,
    alpha: float,
    n_events: int,
    min_events_pass: int,
    tie_sensitive: bool,
    n_distinct: int,
    tie_cutoff: int,
) -> tuple[str, str]:
    """A group row's verdict and its reason, in `_verdict`'s precedence without a bench."""
    if p_value is None:
        return VERDICT_NOT_COMPUTED, reason
    if tie_sensitive and n_distinct < tie_cutoff:
        return (
            VERDICT_TIES,
            f"{n_distinct} distinct durations, below the cutoff of {tie_cutoff}",
        )
    if p_value <= alpha:
        return VERDICT_FAIL, "rejected"
    if n_events < min_events_pass:
        return VERDICT_INCONCLUSIVE, f"n_events {n_events} < {min_events_pass}"
    return VERDICT_PASS, ""


@dataclass(frozen=True)
class GroupMembers:
    """A group's members on one (threshold, side, clock), ready for `min_p`.

    `columns` maps each answering member to `[observed, *null]`, a two-sided member through
    `|.|`. `absent` names each member left out and why. `blocking` is set when a member has
    no answer and no refusal of the record, which leaves the group `not computed`.
    """

    columns: dict[str, np.ndarray]
    absent: list[str]
    blocking: str


def group_members(
    results: list[CheckResult], group: str, *, include_tau_checks: bool
) -> GroupMembers:
    """Which of `group`'s members answer, which are left out and why (`group_rows`)."""
    by_key = {row_key(r): r for r in results}
    columns: dict[str, np.ndarray] = {}
    absent: list[str] = []
    for key in GROUPS[group]:
        check = key[0]
        short = MEMBER_SHORT[check]
        result = by_key.get(key)
        if result is None:
            if check in _TAU_CHECKS and not include_tau_checks:
                absent.append(f"{short} not asked")
                continue
            return GroupMembers(columns, absent, f"{short} has no row")
        if result.p_value is not None and result.null_statistics is not None:
            column = np.concatenate([[result.statistic], result.null_statistics])
            columns[short] = np.abs(column) if check in TWO_SIDED_MEMBERS else column
            continue
        reason = _first_clause(str(result.extra.get(UNDEFINED, result.notes)))
        if result.p_value is None and result.extra.get(REFUSED):
            absent.append(f"{short}: {reason}")
            continue
        return GroupMembers(
            columns,
            absent,
            f"{short} has no answer, and not by a refusal of the record: {reason}",
        )
    return GroupMembers(columns, absent, "")


def group_rows(
    results: list[CheckResult],
    *,
    dataset_id: str,
    label: str,
    side: str,
    clock: str,
    alpha: float,
    n_events: int,
    n_read: int,
    n_used: int,
    min_events_pass: int,
    n_distinct: int,
    tie_fraction: float,
    tie_cutoff: int,
    include_tau_checks: bool,
    n_segments: int,
    n_censored_dropped: int,
    read_note: str = "",
    degenerate_note: str = "",
) -> list[dict[str, object]]:
    """The two group rows of one (threshold, side, clock), from its check results.

    Each group is the single-step min-P (`min_p`) of its members that answered, on their
    kept draws (`CheckResult.null_statistics`), a two-sided member through `|.|` as its own
    p reads it. A member the run-set left out reads `not asked`; a member with no answer is
    named with its reason, and may be left out only when that reason is a refusal of the
    record (`extra[REFUSED]`). Any other failure, or no member answering, leaves the group
    `not computed`. The verdict follows `_group_verdict`; no bench cell is read. A group
    row's `statistic` is its smallest member p on the record, and its `p_value` the min-P
    adjusted p of that smallest p; `n_segments` and `n_censored_dropped` are the record's.
    A group whose members all read `n - 1` events of a segment that ends at one
    (`EVENT_TIME_CHECKS`: the trend group) counts and judges by `n_read`, as their rows do,
    and carries their `read_note`; the dependence group reads every event, `n_events`.

    Validity assumptions

    - Assumption: under the null the B + 1 rows of the shared permutation set are
      exchangeable, and the set of members that answer does not change with which row is
      taken as observed.
      Diagnostic: only refusals that read sizes, gap multisets and tau, or the augmented
      ensemble as a whole, may leave a member out (`checks.result.CheckRefused`).
      Consequence of violation: the group p loses its size.
      Reference: Westfall, P. H. and Young, S. S. (1993), Resampling-Based Multiple
      Testing, Wiley, the single-step min-P adjustment; no section locator located.
    """
    rows: list[dict[str, object]] = []
    for group, members in GROUPS.items():
        found = group_members(results, group, include_tau_checks=include_tau_checks)
        columns, absent, blocking = found.columns, found.absent, found.blocking
        outcome = None
        if not blocking and columns:
            outcome = min_p(columns)
        elif not blocking:
            blocking = "no member answered"
        tie_sensitive = any(key[0] in TIE_SENSITIVE_CHECKS for key in members)
        reads_n_minus_1 = all(key[0] in EVENT_TIME_CHECKS for key in members)
        n_group = n_read if reads_n_minus_1 else n_events
        p_value = None if outcome is None else outcome.p_value
        verdict, reason = _group_verdict(
            p_value,
            blocking,
            alpha=alpha,
            n_events=n_group,
            min_events_pass=min_events_pass,
            tie_sensitive=tie_sensitive,
            n_distinct=n_distinct,
            tie_cutoff=tie_cutoff,
        )
        parts = [] if outcome is None else [f"members {', '.join(columns)}"]
        if absent:
            parts.append("no answer: " + "; ".join(absent))
        if outcome is not None:
            parts.append(
                "member p "
                + ", ".join(f"{m} {v:.3g}" for m, v in outcome.member_p_values.items())
                + f"; smallest {', '.join(outcome.argmin_labels)}"
            )
        parts += [reason, read_note if reads_n_minus_1 else "", degenerate_note]
        rows.append(
            {
                "dataset_id": dataset_id,
                "threshold_label": label,
                "side": side,
                "clock": clock,
                "check_id": group,
                "calibration": CALIB_MIN_P,
                "variant": "",
                "statistic": (
                    float("nan")
                    if outcome is None
                    else float(min(outcome.member_p_values.values()))
                ),
                "p_value": p_value,
                "verdict": verdict,
                "n_events": n_group,
                "n_windows_used": int(n_used),
                "n_censored_dropped": int(n_censored_dropped),
                "n_segments": int(n_segments),
                "n_distinct_durations": n_distinct,
                "tie_fraction": tie_fraction,
                "bench_n_target": None,
                "bench_size_at_n": float("nan"),
                "bench_size_ci": "",
                "bench_size_cell": "",
                "bench_accepted": None,
                "notes": "; ".join(part for part in parts if part),
            }
        )
    return rows
