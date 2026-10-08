"""The check ledger: every check, on every threshold, with what its answer is worth.

A p-value alone is not a result here. The same 0.4 means "no evidence of dependence" from
a calibrated check on 355 windows and "this test cannot see anything" from the same check
on 20. So every row carries three things beside the p-value, and a `pass` requires all
three to hold:

    n_events            enough events for the check to have power
    bench_size_at_n     the check was calibrated at that event count, per the bench
    tie_fraction        the durations are not so tied that a rank test is meaningless

If any fails, the verdict is NOT `pass`. That is the entire point of the ledger: a
non-rejection must never print as a pass on a p-value alone.

**Which clock.** Both are run. On the in-spec clock of a carved record, time stops
accruing the moment the record ends out of spec, so `tau == T_N` - measured at ~73% of
synthetic replicates. Kvaloy and Lindqvist's eq (7) (Technometrics 2020, p. 104) is
singular there, and every asymptotic row's limiting null assumes a `tau` chosen
independently of the events. Those rows read `not computed` with
the reason, and the rank checks (which never touch `tau`) still run there.

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

import quebra.analyzers.checks.c3_serial_copula as c3
from quebra.analyzers.calibration_summary import (
    bench_acceptance_at_n,
    nearest_bracketing_n,
)
from quebra.analyzers.checks._multiprocess import segments_from_windows
from quebra.analyzers.checks._permutation import block_permutations
from quebra.analyzers.checks._rank_serial import MAX_LAG_CAP
from quebra.analyzers.checks.battery import run_battery
from quebra.analyzers.checks.result import (
    CALIB_ASYMPTOTIC,
    CALIB_PERMUTATION,
    CALIB_R_COPULA,
    CLOCK_CALENDAR,
    CLOCK_IN_SPEC,
    CheckResult,
)
from quebra.core._artifact_guard import StaleArtifactGuard

VERDICT_PASS = "pass"
VERDICT_FAIL = "fail"
VERDICT_UNDERPOWERED = "underpowered"
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
    VERDICT_UNDERPOWERED,
    VERDICT_TIES,
    VERDICT_NOT_COMPUTED,
)

# HOW OFTEN THIS RULE BINDS, measured rather than assumed: on the T2* ladder it never does
# except where the event count already disqualifies the row. T2* durations are wall-clock
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

LEDGER_COLUMNS = [
    "dataset_id",
    "threshold_label",
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
    """One row per (threshold, clock, check, calibration), plus the settings that made it.

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
    # n = 355 (see c3_serial_copula for the measured range and its machine). Skipping it
    # loses no calibrated evidence - C3 has no bench cell, so its
    # rows are `not computed` or an uncalibrated p-value either way. When False the C3 rows
    # are OMITTED rather than written as `not computed`: a blank row would claim the check
    # was attempted and failed, when in fact it was never asked.
    include_c3: bool = True
    # Derived from a run-set by `check_outcome.battery_flags`, never set by hand once one is
    # in play. COARSE on purpose: `include_tau_checks` is all-or-nothing over five rows, so a
    # run-set naming one asymptotic row still computes the other four. The declared run-set
    # therefore names AT MOST what ran, not exactly what ran, and a reader comparing the two
    # will find the ledger carrying rows the run-set does not list.
    include_tau_checks: bool = True
    include_c2_asymptotic: bool = True
    # See c3_serial_copula.N_NULL_SIM: the dominant cost, and part of what a C3
    # p-value means, so it is carried rather than left to a module default.
    c3_n_null_sim: int = c3.N_NULL_SIM


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
    include_c3: bool = True,
    include_tau_checks: bool = True,
    include_c2_asymptotic: bool = True,
    c3_n_null_sim: int = c3.N_NULL_SIM,
) -> CheckLedgerInputs:
    """Build inputs from a `WindowsResult`, taking the gap spans with it.

    `gap_spans_s` is not optional in spirit: without it `segments_from_windows` falls back
    to reading birth types alone, which cannot see a gap flanked by out-of-spec reads and
    silently merges two renewal processes separated by unobserved hours.
    """
    diagnostics = getattr(window_result, "diagnostics", {}) or {}
    return CheckLedgerInputs(
        windows=window_result.windows_in_spec,
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
    """The verdict and the reason it is not `pass`, in precedence order."""
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
        # below demotes `bench_accepted is None` to `underpowered` because a non-rejection
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
        return VERDICT_UNDERPOWERED, "; ".join(unmet)
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
    for label, _value, _big_good in inputs.thresholds:
        table = inputs.windows[inputs.windows["threshold_label"] == label]
        for clock in (CLOCK_IN_SPEC, CLOCK_CALENDAR):
            rows.extend(
                _rows_for(
                    table=table,
                    label=label,
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
    clock: str,
    inputs: CheckLedgerInputs,
    note: str,
    n_windows: int = 0,
) -> dict:
    return {
        "dataset_id": inputs.dataset_id,
        "threshold_label": label,
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
    clock: str,
    inputs: CheckLedgerInputs,
    acceptance: pd.DataFrame,
    grid: list[int],
) -> list[dict[str, object]]:
    if not len(table):
        return [_blank_row(label, clock, inputs, "no windows at this threshold")]
    try:
        segments, n_dropped = segments_from_windows(
            table,
            clock=clock,
            min_events=inputs.min_events_per_segment,
            gap_spans_s=inputs.gap_spans_s,
        )
    except (ValueError, KeyError) as exc:
        return [
            _blank_row(
                label, clock, inputs, f"segmentation failed: {exc}"[:200], len(table)
            )
        ]
    if not segments:
        return [
            _blank_row(
                label,
                clock,
                inputs,
                f"no segment reached {inputs.min_events_per_segment} events "
                f"({n_dropped} dropped)",
                len(table),
            )
        ]

    durations = np.concatenate([s.x for s in segments])
    n_distinct, tie_fraction = _tie_stats(durations)
    n_events = int(len(durations))
    bench_n = nearest_bracketing_n(grid, n_events)

    stream = stream_for(label, clock)
    rng = np.random.default_rng([inputs.seed, stream])
    perm = block_permutations(
        [s.n_events for s in segments], inputs.n_permutations, rng
    )

    results: list[CheckResult] = []
    battery_note = ""
    try:
        results = run_battery(
            segments,
            clock=clock,
            perm=perm,
            max_lag=inputs.lag_max,
            include_tau_checks=inputs.include_tau_checks,
            include_c2_asymptotic=inputs.include_c2_asymptotic,
        )
    except (ValueError, KeyError) as exc:
        battery_note = f"C1/C2 undefined here: {exc}"[:180]
        try:
            results = run_battery(
                segments,
                clock=clock,
                perm=perm,
                max_lag=inputs.lag_max,
                # Hard False, not the input: this branch is a mathematical refusal (tau is
                # singular here), and it must override a run-set that asked for those rows.
                include_tau_checks=False,
                include_c2_asymptotic=inputs.include_c2_asymptotic,
            )
        except (ValueError, KeyError) as inner:
            return [
                _blank_row(
                    label, clock, inputs, f"battery failed: {inner}"[:200], len(table)
                )
            ]
        # The retry dropped C1 and C2. Say so with their own rows: without them the reader
        # sees a ladder where the trend checks simply vanish at the loose thresholds and
        # cannot tell whether they were attempted, declined, or forgotten. The reason
        # (usually tau == T_N on the in-spec clock) belongs on the row it applies to.
        for dropped, calibrations in (
            ("c1_lewis_robinson", (CALIB_ASYMPTOTIC, CALIB_PERMUTATION)),
            ("c2_anderson_darling", (CALIB_ASYMPTOTIC, CALIB_PERMUTATION)),
        ):
            for calibration in calibrations:
                results.append(
                    CheckResult(
                        check=dropped,
                        statistic=float("nan"),
                        p_value=None,
                        calibration=calibration,
                        clock=clock,
                        n_events=int(sum(s.n_events for s in segments)),
                        n_segments=len(segments),
                        notes=battery_note,
                    )
                )

    # C3 is out-of-process and shares nothing with the permutation set, so it is called
    # separately. With R absent it returns p_value=None and the row reads `not computed`.
    if inputs.include_c3:
        try:
            results.append(
                c3.run(
                    segments,
                    clock=clock,
                    max_lag=inputs.lag_max,
                    seed=inputs.seed,
                    n_null_sim=inputs.c3_n_null_sim,
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

    out: list[dict[str, object]] = []
    for result in results:
        variant = ""
        for token in result.notes.split():
            if token.startswith("variant="):
                variant = token.split("=", 1)[1]
        match = acceptance[
            (acceptance["check"] == result.check)
            & (acceptance["calibration"] == result.calibration)
            & (acceptance["variant"].fillna("") == variant)
            & (acceptance["clock"] == clock)
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
            alpha=inputs.alpha,
            n_events=n_events,
            min_events_pass=inputs.min_events_pass,
            n_distinct=n_distinct,
            tie_cutoff=inputs.tie_cutoff_distinct,
            bench_accepted=accepted,
        )
        # The battery note explains why C1/C2 were dropped, so it belongs on THEIR rows.
        # Appending it to every row put "C1/C2 undefined here: tau ..." on a C3 row, which
        # is a different check with a different reason.
        relevant = battery_note if result.check.startswith(("c1_", "c2_")) else ""
        notes = "; ".join(part for part in (reason, relevant) if part)
        out.append(
            {
                "dataset_id": inputs.dataset_id,
                "threshold_label": label,
                "clock": clock,
                "check_id": result.check,
                "calibration": result.calibration,
                "variant": variant,
                "statistic": float(result.statistic),
                "p_value": result.p_value,
                "verdict": verdict,
                "n_events": n_events,
                "n_windows_used": int(len(table)),
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
    return out
