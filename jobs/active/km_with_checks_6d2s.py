"""Kaplan-Meier survival across five 6D2S records, with the independence checks beside it.

Two figures, deliberately: the survival band, and a separate grid of check verdicts. The band is ALWAYS computed and ALWAYS drawn - a rejected check annotates it and
never suppresses it, because control flow that depends on what the data happened to say is
unpredictable and a reader is better served by a band they are told not to trust than by no
band at all.

This re-derives `km_poster_6d2s`'s carve rather than including it, so the two carves are
held equal by `tests/test_windows_carve.py`, which reads both jobs' `gap_mult` and asserts
the two configurations agree on their window tables.

The run-set and the display-set are STEP KWARGS, so they enter the run identity and land on
the Mermaid label. `analyzers/check_outcome.py` explains why they are two independent tuples
over the full `(check_id, calibration, variant)` key rather than one list of checks.
"""

from __future__ import annotations

import pandas as pd

from quebra.analyzers import check_ledger, kaplan_meier, t2star, windows
from quebra.analyzers.assumptions import A1_RENEWAL_DURATIONS
from quebra.analyzers.check_outcome import (
    PERMUTATION_KEYS,
    CheckOutcome,
    battery_flags,
    build_check_outcome,
    resolve,
)
from quebra.analyzers.kaplan_meier import KaplanMeierComparison, KaplanMeierCurve
from quebra.analyzers.t2star import T2StarResult
from quebra.analyzers.windows import WindowsResult
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.plots.check_outcome_plot import CheckOutcomePlot
from quebra.plots.km_survival_plot import KMSurvivalPlot
from quebra.recipes import RAMSEY_CONFIG, _filter_step, _final_stage
from quebra.schemas.track912 import track912Schema

JOB_ID = "km_with_checks_6d2s"
JOB_FAMILY = "independence"
# Five ledgers at 999 permutations is long, and `discovery.swept` is opt-out. Same reason
# the `check_ledger` jobs opt out; selectable with `--family independence` or by path.
JOB_SWEEP = False

PREFIX = "km_with_checks_6d2s"

THRESHOLD_LABEL = "3.0 µs"
THRESHOLD: list[tuple[str, float, bool]] = [(THRESHOLD_LABEL, 3.0e-6, True)]

# Carve parameters. MUST equal `jobs/active/km_poster_6d2s.py` - the differential test is what
# holds them equal, and if they diverge the grids describe a different carve from the band.
GAP_MULT = 10.0

# Verdict parameters, equal to `ledger_recipe.LEDGER_KNOBS` so this job's verdicts and the
# `check_ledger` jobs' are comparable, except SEED, which is deliberately different: an
# identical permutation seed across two jobs would make their Monte Carlo error identical
# rather than independent.
ALPHA = 0.05
MIN_EVENTS_PASS = 35
TIE_CUTOFF_DISTINCT = 5
LAG_MAX = 5
N_PERMUTATIONS = 999
SEED = 20260903

# WHICH CHECKS RUN, and separately WHICH ARE SHOWN.
#
# Run-set: the permutation rows. At the event counts these records carry, the asymptotic
# calibrations are measured as oversized (`docs/iid_checks/C1_lewis_robinson.md`), and C3 is
# excluded because it has no bench cell, so its rows carry no calibrated evidence. It also
# shells out to R, which costs about a minute at the largest n this ladder produces.
# Those reasons are about entitlement to the data, not about convenience.
RUN_SET = PERMUTATION_KEYS
# Display-set: C1 and C6. One asks whether durations TREND, the other whether their ORDER
# carries information; between them they cover the two ways the renewal assumption fails that
# a reader of a survival curve most needs to know about. Everything in RUN_SET is still in the
# artifact and in provenance - this only chooses what the figure shows.
DISPLAY_SET = (
    ("c1_lewis_robinson", "permutation", ""),
    ("c6_exchangeability", "permutation", ""),
)

CLOCKS = ("in_spec", "calendar")

# The record in `analyzers/assumptions.py` this band rests on. Named on the artifact so a
# figure can say WHICH assumption the checks below are evidence about.
ASSUMPTION_ID = A1_RENEWAL_DURATIONS.id

DATASETS: list[tuple[str, int, str]] = [
    ("280623_6D2S_qubit2", 2, "280623 qubit 2"),
    ("040423_6D2S_qubit1", 1, "040423 qubit 1"),
    ("220423_6D2S_qubit1", 1, "220423 qubit 1"),
    ("090623_6D2S_qubit6", 6, "090623 qubit 6"),
    ("070723_6D2S_qubit4", 4, "070723 qubit 4"),
]

BENCH_SIZE_TABLE = Dataset(path="jobs/bench/results/size_table.csv", schema=None)

job = Job(PREFIX)

_bench = job.load_df(BENCH_SIZE_TABLE)


def _t2star_run(norm: object) -> T2StarResult:
    return t2star.run(t2star.make_inputs_from_norm(norm))  # type: ignore[arg-type]


def _windows_run(result: T2StarResult, gap_mult: float) -> WindowsResult:
    return windows.run(
        windows.make_inputs_from_frame(
            result.frame,
            time_col="t_rel_s",
            value_col="t2star_s",
            thresholds=THRESHOLD,
            dataset_id=str(result.meta.get("dataset_id", "")),
            gap_mult=gap_mult,
        )
    )


def _km_run(
    window_result: WindowsResult, threshold_label: str, label: str
) -> KaplanMeierCurve:
    return kaplan_meier.run(
        kaplan_meier.make_inputs_from_windows(
            window_result.windows_in_spec,
            threshold_label=threshold_label,
            label=label,
            dataset_id=str(window_result.meta.get("dataset_id", "")),
        )
    )


def _ledger_run(
    window_result: WindowsResult,
    bench: pd.DataFrame,
    *,
    run_set: tuple[tuple[str, str, str], ...],
    alpha: float,
    min_events_pass: int,
    tie_cutoff_distinct: int,
    lag_max: int,
    n_permutations: int,
    seed: int,
) -> check_ledger.CheckLedger:
    """Score one record's checks. Every battery switch is DERIVED from the run-set.

    The switches are coarser than the run-set, so this reduces what is computed without
    making the two match exactly: `include_tau_checks` is all-or-nothing over five rows.
    The ledger may therefore carry a row the run-set does not name.
    """
    flags = battery_flags(run_set)
    inputs = check_ledger.make_inputs_from_windows(
        window_result,
        bench,
        thresholds=THRESHOLD,
        alpha=alpha,
        min_events_pass=min_events_pass,
        tie_cutoff_distinct=tie_cutoff_distinct,
        lag_max=lag_max,
        n_permutations=n_permutations,
        seed=seed,
        include_c3=flags["include_c3"],
        include_tau_checks=flags["include_tau_checks"],
        include_c2_asymptotic=flags["include_c2_asymptotic"],
    )
    return check_ledger.run(inputs)


def _outcome(
    *ledgers: check_ledger.CheckLedger,
    run_set: tuple[tuple[str, str, str], ...],
    display_set: tuple[tuple[str, str, str], ...],
    threshold_label: str,
    alpha: float,
) -> CheckOutcome:
    """Every drawn record's rows, reshaped into one grid per shown check."""
    rows = pd.concat([led.rows for led in ledgers], ignore_index=True)
    return build_check_outcome(
        rows,
        threshold_label=threshold_label,
        run_set=run_set,
        display_set=display_set,
        # Ledger rows are keyed by `dataset_id`, the file STEM. The display label is a
        # different string, so passing it as the key finds nothing and renders as absent.
        datasets=tuple(stem for stem, _qubit, _display in DATASETS),
        dataset_labels=tuple(display for _stem, _qubit, display in DATASETS),
        clocks=CLOCKS,
        alpha=alpha,
    )


def _compare(
    *curves_and_outcome: object, threshold_label: str, assumption_id: str
) -> KaplanMeierComparison:
    """Rank the pairs, then ATTACH the check outcome to the band artifact.

    The outcome is the last input, not a kwarg: it is a node, and a typed object cannot be a
    step kwarg (`core/closure.py` renders only scalars and containers of scalars). Attaching
    here rather than drawing a sibling figure is what makes `AGENTS.md` section 5 structural -
    a reader who opens the band alone can still see what was checked.
    """
    *curves, outcome = curves_and_outcome
    comparison = kaplan_meier.compare(
        [c for c in curves if isinstance(c, KaplanMeierCurve)], threshold_label
    )
    assert isinstance(outcome, CheckOutcome)
    comparison.assumption_id = assumption_id
    comparison.checks_asked = tuple(
        " ".join(p for p in k if p) for k in outcome.run_set
    )
    comparison.checks_unanswered = tuple(
        " ".join(p for p in k if p) for k in outcome.asked_but_unanswered
    )
    comparison.check_verdicts = tuple(
        (grid.label, label, verdict)
        for grid in outcome.grids
        for label, row in zip(grid.dataset_labels, grid.verdicts, strict=True)
        for verdict in row
    )
    return comparison


# Resolved once, at build time, so a display-set outside the run-set fails when the graph is
# constructed rather than after five ledgers have been computed.
_RUN, _DISPLAY = resolve(RUN_SET, DISPLAY_SET)

_curves = []
_ledgers = []
for stem, qubit, display in DATASETS:
    dataset = Dataset(
        path=f"data/real_private/6D2S/{stem}.pickle",
        schema=track912Schema,
        qubit=qubit,
        device="6D2S",
        extra={"run_name": stem},
    )
    _raw = job.load(dataset)
    _filtered = job.step(_filter_step(RAMSEY_CONFIG), _raw, name=f"filter_{stem}")
    _final = job.step(_final_stage, _filtered, name=f"final_filter_stage_{stem}")
    _result = job.step(_t2star_run, _final, name=f"t2star_{stem}")
    _windows = job.step(
        _windows_run, _result, name=f"windows_{stem}", gap_mult=GAP_MULT
    )
    _curves.append(
        job.step(
            _km_run,
            _windows,
            name=f"kaplan_meier_{stem}",
            threshold_label=THRESHOLD_LABEL,
            label=display,
        )
    )
    _ledgers.append(
        job.step(
            _ledger_run,
            _windows,
            _bench,
            name=f"check_ledger_{stem}",
            run_set=_RUN,
            alpha=ALPHA,
            min_events_pass=MIN_EVENTS_PASS,
            tie_cutoff_distinct=TIE_CUTOFF_DISTINCT,
            lag_max=LAG_MAX,
            n_permutations=N_PERMUTATIONS,
            seed=SEED,
        )
    )

_check_outcome = job.step(
    _outcome,
    *_ledgers,
    name="check_outcome",
    run_set=_RUN,
    display_set=_DISPLAY,
    threshold_label=THRESHOLD_LABEL,
    alpha=ALPHA,
)

# Both materialized: the ranking behind the drawn pair, and the check outcome. The run-set is
# NAMED on the outcome, not resolved into it - rows that ran but are not displayed survive as
# key tuples, not as verdicts.
_comparison = job.step(
    _compare,
    *_curves,
    _check_outcome,
    name="km_compare",
    threshold_label=THRESHOLD_LABEL,
    assumption_id=ASSUMPTION_ID,
)

job.materialize(_comparison, name=f"{PREFIX}_comparison")
# Not `{PREFIX}_check_grid`: that is the figure's safe-named title below, and a materialize
# sharing it would write the same `provenance/<name>.prov.json`, so whichever sink ran second
# would overwrite the first's record. The runner refuses the collision outright.
job.materialize(_check_outcome, name=f"{PREFIX}_check_outcome")

job.figure(
    KMSurvivalPlot,
    _comparison,
    targets=["static", "academic"],
    title=f"{PREFIX}_km_survival",
)
job.figure(
    CheckOutcomePlot,
    _check_outcome,
    targets=["static", "academic"],
    title=f"{PREFIX}_check_grid",
)
