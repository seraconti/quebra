"""The record figures: the blocks a figure of one record draws, as typed artifacts.

One record, one PDF per block. A block is a typed artifact built here, with every number it
draws, and one plot class in `plots/record_figure_plot.py` that draws it; the jobs wire both
through `record_figure_recipe.BLOCKS`. Nothing here imports `plots/` or `panels/`.

A record's windows are the carve its check ledger scored, recomputed from the ledger job's
own carve kwargs and held to that ledger by `cross_check_carve`: at every rung and bag the
carve must hold as many windows as the ledger's `n_windows_used` says. A mismatch is a
wiring error and raises, so nothing is drawn on windows the ledger did not score. The
builders raise on wiring errors only; an empty bag or a ledger row that gave no answer is
drawn and labelled.

**Survival block** (`build_survival_block`). One Kaplan-Meier curve with its band per rung
and bag, both bags drawn, each band with the check outcome of exactly that band: the ledger
rows of that rung and that bag, on both clocks, over the run-set the ledger ran
(`ledger_run_set`). A ledger row `(all)`, written when a bag could not be scored at all,
stands for every check of the run-set; it is expanded into one row per check before the
outcome is built (`expand_all_rows`), so its verdict and note reach every cell instead of
reading `no row`. The bag rides on the outcome's dataset label (`{record}, {side} bag`), so
an attachment read alone says which band it belongs to. A record whose carve failed draws
the failure's note and no curve.

**Hazard block** (`build_hazard_block`). The Nelson-Aalen cumulative hazard per rung and
bag, with its band, from `survival.nelson_aalen_set` on the same carve. Each item's own
`checks` field is filled with the outcome of the same rung and bag as the survival band, so
the hazard band carries its check outcome as a field and a reader of the item alone sees
it. The ages are converted to minutes here, the unit the survival block draws in.

**Curve summaries** (`build_record_summaries`, `summaries_table`). Each record's quantiles
and restricted mean per rung and bag (`survival.curve_summaries_set`, tau each curve's own
longest age), each item's `checks` field filled as the hazard items' are, and the
restricted mean's Wald interval from `mtbf.rmst_interval`, the function the MTBC uses. The
job's table holds one row per drawn record, rung and bag, unrounded and keyed on
`threshold_label`; a ledger record outside the paper's tables gets one row with
`drawn = False` and the reason, and a record whose carve failed one row with its note. A
quantile or an interval end the curve never reaches is NaN, and so is every number of an
empty bag (`n_windows` 0). The interval is kept as computed, as the MTBC's is, and
`rmst_interval_note` says when it has no width (no death before tau, so SE 0) or its lower
end is at or below 0. `check_verdicts` lists each check's verdicts on the side clock, then
the calendar clock.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import numpy as np
import pandas as pd

from quebra.analyzers.check_attachment import A1_RENEWAL_DURATIONS_ID, CheckAttachment
from quebra.analyzers.check_ledger import VERDICT_ABSENT, VERDICTS, CheckLedger
from quebra.analyzers.check_outcome import (
    PERMUTATION_KEYS,
    CheckOutcome,
    RowKey,
    build_check_outcome,
    check_attachment_of,
    resolve,
)
from quebra.analyzers.checks.result import CLOCKS
from quebra.analyzers.event_table import SIDES, event_tables_from_carve
from quebra.analyzers.kaplan_meier import (
    KaplanMeierCurve,
    KaplanMeierSet,
    kaplan_meier_set,
)
from quebra.analyzers.mtbf import rmst_interval
from quebra.analyzers.survival import (
    DEFAULT_QUANTILES,
    ByThresholdSide,
    CurveSummaries,
    NelsonAalen,
    curve_summaries_set,
    nelson_aalen_set,
)
from quebra.analyzers.windows import WindowsResult
from quebra.core._artifact_guard import StaleArtifactGuard
from quebra.core.types import RecordFailure

# The check id of the one row the ledger writes for a rung, bag and clock it could not score.
ALL_ROWS = "(all)"
# The reason a ledger record outside the paper's dataset tables carries in the summaries.
EXCLUDED_REASON = "not in the paper's dataset tables, so not drawn"
# `rmst_interval_note` of a restricted mean whose SE is 0, and of one whose interval
# reaches 0.
NOTE_NO_WIDTH = "SE 0: no death before tau, so the interval has no width"
NOTE_LOWER_AT_ZERO = "lower end at or below 0, reported as computed"
# The verdicts a survival block's key names: what the ledger can say, then what a builder
# writes for a cell with no row.
VERDICT_KEY = (*VERDICTS, VERDICT_ABSENT)


def ledger_run_set(*, include_c3: bool) -> tuple[RowKey, ...]:
    """The run-set a ledger ran, from its node's kwargs: the battery's permutation rows.

    A ledger run with `include_c3` holds C3 rows, which no run-set can select
    (`check_outcome.ALL_KEYS`), so it raises. A window ledger's group rows
    (`check_ledger.GROUPS`) are no run-set's either: a band's check outcome reads single
    checks only.
    """
    if include_c3:
        raise ValueError(
            "the ledger ran C3, whose rows no run-set can select; a record figure reads "
            "a ledger of the permutation rows"
        )
    return PERMUTATION_KEYS


def expand_all_rows(rows: pd.DataFrame, run_set: tuple[RowKey, ...]) -> pd.DataFrame:
    """Every `(all)` row as one row per check of `run_set`, with its verdict, counts and
    note; every other row unchanged, in order."""
    expanded: list[dict[str, object]] = []
    for row in rows.to_dict("records"):
        if row["check_id"] != ALL_ROWS:
            expanded.append(row)
            continue
        expanded += [
            {**row, "check_id": check, "calibration": calibration, "variant": variant}
            for check, calibration, variant in run_set
        ]
    return pd.DataFrame(expanded, columns=rows.columns)


def cross_check_carve(
    carve: WindowsResult | RecordFailure, ledger: CheckLedger
) -> None:
    """Raise unless `carve` is the carve `ledger` scored.

    The same record; at every rung of the ledger and every bag, as many windows as each
    ledger row of that rung and bag used, `(all)` rows included; a failed carve against a
    ledger whose every row carries that failure's note.
    """
    record = (
        carve.dataset_id
        if isinstance(carve, RecordFailure)
        else str(carve.meta["dataset_id"])
    )
    if record != ledger.dataset_id:
        raise ValueError(
            f"the carve of {record!r} is checked against the ledger of "
            f"{ledger.dataset_id!r}"
        )
    rows = ledger.rows
    if isinstance(carve, RecordFailure):
        notes = set(rows["notes"].astype(str))
        if notes != {carve.note()}:
            raise ValueError(
                f"{record}: the carve failed ({carve.note()}) but the ledger rows say "
                f"{sorted(notes)[:3]}"
            )
        return
    windows = pd.concat([carve.windows_in_spec, carve.windows_out_of_spec])
    for label, _value, _big_good in ledger.thresholds:
        for side in SIDES:
            n = int(
                (
                    (windows["threshold_label"] == label) & (windows["side"] == side)
                ).sum()
            )
            on_cell = (rows["threshold_label"] == label) & (rows["side"] == side)
            used = {int(v) for v in rows.loc[on_cell, "n_windows_used"]}
            if used != {n}:
                raise ValueError(
                    f"{record}, {label}, {side}: the carve holds {n} windows and the "
                    f"ledger rows used {sorted(used)}; the figure would draw windows the "
                    "ledger did not score"
                )


def kaplan_meier_of(
    carve: WindowsResult | RecordFailure, *, conf_level: float
) -> KaplanMeierSet | RecordFailure:
    """The Kaplan-Meier set of a record's carve; a failed carve passes through."""
    if isinstance(carve, RecordFailure):
        return carve
    return kaplan_meier_set(event_tables_from_carve(carve), conf_level=conf_level)


def nelson_aalen_of(
    carve: WindowsResult | RecordFailure, *, conf_level: float
) -> ByThresholdSide[NelsonAalen] | RecordFailure:
    """The Nelson-Aalen set of a record's carve; a failed carve passes through."""
    if isinstance(carve, RecordFailure):
        return carve
    return nelson_aalen_set(event_tables_from_carve(carve), conf_level=conf_level)


def _band_outcomes(
    ledger: CheckLedger, ladder: tuple[str, ...], run: tuple[RowKey, ...]
) -> dict[tuple[str, str], CheckOutcome]:
    """The check outcome of every band of `ladder`, by (rung, side): that rung's and bag's
    ledger rows on both clocks, `(all)` rows expanded over `run` (module docstring)."""
    rows = expand_all_rows(ledger.rows, run)
    return {
        (label, side): build_check_outcome(
            rows,
            threshold_label=label,
            run_set=run,
            display_set=run,
            datasets=(ledger.dataset_id,),
            dataset_labels=(f"{ledger.dataset_id}, {side} bag",),
            clocks=CLOCKS,
            alpha=ledger.alpha,
            side=side,
        )
        for label in ladder
        for side in SIDES
    }


@dataclass
class SurvivalBlockData(StaleArtifactGuard):
    """Everything the survival block of one record draws.

    `curves`, `outcomes` and `attachments` are keyed by (rung, side), every rung of
    `ladder` and every side of `sides`. `outcomes[key]` is the check outcome of the band
    `curves[key]`: one grid per check of `run_set`, one row (the record) and one column
    per clock of `clocks`. `attachments[key]` is that outcome as the fields a band carries,
    naming a1. `failure` is the carve's failure note, empty when it carved; a failed record
    has no curves. `flag` is the text a record the paper's dataset tables shade carries,
    empty otherwise. `verdict_key` is the verdicts the key names, in order.
    """

    dataset_id: str
    flag: str
    failure: str
    ladder: tuple[str, ...]
    sides: tuple[str, ...]
    clocks: tuple[str, ...]
    run_set: tuple[RowKey, ...]
    conf_level: float
    curves: dict[tuple[str, str], KaplanMeierCurve]
    outcomes: dict[tuple[str, str], CheckOutcome]
    attachments: dict[tuple[str, str], CheckAttachment]
    verdict_key: tuple[str, ...]


def build_survival_block(
    km: KaplanMeierSet | RecordFailure,
    ledger: CheckLedger,
    *,
    run_set: tuple[RowKey, ...],
    conf_level: float,
    flag: str,
) -> SurvivalBlockData:
    """The survival block of one record: its curves and, per band, the check outcome of
    that rung and bag (module docstring)."""
    run, _display = resolve(run_set, run_set)
    record = ledger.dataset_id
    if isinstance(km, RecordFailure):
        failure = km.note()
        ladder = tuple(label for label, _v, _b in ledger.thresholds)
        curves: dict[tuple[str, str], KaplanMeierCurve] = {}
    else:
        failure = ""
        ladder = tuple(km.ladder)
        curves = {
            (label, side): km.curve(label, side) for label in ladder for side in SIDES
        }
        if km.conf_level != conf_level:
            raise ValueError(
                f"the curves were built at {km.conf_level}, the block declares {conf_level}"
            )
    outcomes = _band_outcomes(ledger, ladder, run)
    return SurvivalBlockData(
        dataset_id=record,
        flag=flag,
        failure=failure,
        ladder=ladder,
        sides=SIDES,
        clocks=CLOCKS,
        run_set=run,
        conf_level=conf_level,
        curves=curves,
        outcomes=outcomes,
        attachments={
            key: check_attachment_of(outcome, A1_RENEWAL_DURATIONS_ID)
            for key, outcome in outcomes.items()
        },
        verdict_key=VERDICT_KEY,
    )


@dataclass
class HazardBlockData(StaleArtifactGuard):
    """Everything the hazard block of one record draws.

    `hazards`, `time_min` and `outcomes` are keyed by (rung, side), every rung of `ladder`
    and every side of `sides`. `hazards[key]` is the Nelson-Aalen item with its `checks`
    field filled from `outcomes[key]`, naming a1; `time_min[key]` is its `time_s` in
    minutes. `failure`, `flag` and `verdict_key` read as on `SurvivalBlockData`.
    """

    dataset_id: str
    flag: str
    failure: str
    ladder: tuple[str, ...]
    sides: tuple[str, ...]
    clocks: tuple[str, ...]
    run_set: tuple[RowKey, ...]
    conf_level: float
    hazards: dict[tuple[str, str], NelsonAalen]
    time_min: dict[tuple[str, str], np.ndarray]
    outcomes: dict[tuple[str, str], CheckOutcome]
    verdict_key: tuple[str, ...]


def build_hazard_block(
    hazards: ByThresholdSide[NelsonAalen] | RecordFailure,
    ledger: CheckLedger,
    *,
    run_set: tuple[RowKey, ...],
    conf_level: float,
    flag: str,
) -> HazardBlockData:
    """The hazard block of one record: each Nelson-Aalen item with the check outcome of its
    rung and bag in its `checks` field (module docstring)."""
    run, _display = resolve(run_set, run_set)
    if isinstance(hazards, RecordFailure):
        failure = hazards.note()
        ladder = tuple(label for label, _v, _b in ledger.thresholds)
        items: dict[tuple[str, str], NelsonAalen] = {}
    else:
        failure = ""
        ladder = tuple(hazards.ladder)
        items = {
            (label, side): hazards.get(label, side)
            for label in ladder
            for side in SIDES
        }
        mismatched = {i.conf_level for i in items.values()} - {conf_level}
        if mismatched:
            raise ValueError(
                f"the hazards were built at {sorted(mismatched)}, the block declares "
                f"{conf_level}"
            )
    outcomes = _band_outcomes(ledger, ladder, run)
    filled = {
        key: dataclasses.replace(
            item, checks=check_attachment_of(outcomes[key], A1_RENEWAL_DURATIONS_ID)
        )
        for key, item in items.items()
    }
    return HazardBlockData(
        dataset_id=ledger.dataset_id,
        flag=flag,
        failure=failure,
        ladder=ladder,
        sides=SIDES,
        clocks=CLOCKS,
        run_set=run,
        conf_level=conf_level,
        hazards=filled,
        time_min={key: item.time_s / 60.0 for key, item in filled.items()},
        outcomes=outcomes,
        verdict_key=VERDICT_KEY,
    )


def curve_summaries_of(
    carve: WindowsResult | RecordFailure, *, conf_level: float
) -> ByThresholdSide[CurveSummaries] | RecordFailure:
    """The curve summaries of a record's carve, each restricted to its own longest age; a
    failed carve passes through."""
    if isinstance(carve, RecordFailure):
        return carve
    return curve_summaries_set(event_tables_from_carve(carve), conf_level=conf_level)


@dataclass
class RecordSummaries(StaleArtifactGuard):
    """One record's curve summaries, keyed by (rung, side), each with its `checks` field
    filled from its band's check outcome, naming a1, and `rmst_interval_s` its restricted
    mean's Wald interval at `conf_level`. `failure` is the carve's failure note, empty when
    it carved; a failed record has no item."""

    dataset_id: str
    flag: str
    failure: str
    ladder: tuple[str, ...]
    conf_level: float
    summaries: dict[tuple[str, str], CurveSummaries]
    rmst_interval_s: dict[tuple[str, str], tuple[float, float]]


def build_record_summaries(
    summaries: ByThresholdSide[CurveSummaries] | RecordFailure,
    ledger: CheckLedger,
    *,
    run_set: tuple[RowKey, ...],
    conf_level: float,
    flag: str,
) -> RecordSummaries:
    """One record's summaries with their check outcomes and intervals (module docstring)."""
    run, _display = resolve(run_set, run_set)
    if isinstance(summaries, RecordFailure):
        failure = summaries.note()
        ladder = tuple(label for label, _v, _b in ledger.thresholds)
        items: dict[tuple[str, str], CurveSummaries] = {}
    else:
        failure = ""
        ladder = tuple(summaries.ladder)
        items = {
            (label, side): summaries.get(label, side)
            for label in ladder
            for side in SIDES
        }
        mismatched = {i.conf_level for i in items.values()} - {conf_level}
        if mismatched:
            raise ValueError(
                f"the summaries were built at {sorted(mismatched)}, the block declares "
                f"{conf_level}"
            )
    outcomes = _band_outcomes(ledger, ladder, run)
    filled = {
        key: dataclasses.replace(
            item, checks=check_attachment_of(outcomes[key], A1_RENEWAL_DURATIONS_ID)
        )
        for key, item in items.items()
    }
    return RecordSummaries(
        dataset_id=ledger.dataset_id,
        flag=flag,
        failure=failure,
        ladder=ladder,
        conf_level=conf_level,
        summaries=filled,
        rmst_interval_s={
            key: rmst_interval(item, conf_level) for key, item in filled.items()
        },
    )


def _quantile_columns(p: float) -> tuple[str, str, str]:
    name = f"q{100 * p:g}"
    return f"{name}_s", f"{name}_lower_s", f"{name}_upper_s"


# The summaries table's columns, in order; the quantiles are `curve_summaries_of`'s.
SUMMARY_COLUMNS = (
    "record",
    "drawn",
    "reason",
    "flag",
    "threshold_label",
    "side",
    "conf_level",
    "n_windows",
    "n_deaths",
    "n_censored",
    "n_unobserved_birth_dropped",
    "max_observed_age_s",
    "rmst_tau_s",
    "rmst_s",
    "rmst_se_s",
    "rmst_lower_s",
    "rmst_upper_s",
    "rmst_interval_note",
    *(name for p in DEFAULT_QUANTILES for name in _quantile_columns(p)),
    "point_standing",
    "interval_standing",
    "assumption_id",
    "checks_asked",
    "checks_unanswered",
    "check_verdicts",
    "check_summary",
)
# Counts, written as integers with an empty cell where a row has none.
COUNT_COLUMNS = ("n_windows", "n_deaths", "n_censored", "n_unobserved_birth_dropped")


def _interval_note(item: CurveSummaries, lower: float) -> str:
    if item.rmst_se_s == 0.0:
        return NOTE_NO_WIDTH
    if lower <= 0.0:
        return NOTE_LOWER_AT_ZERO
    return ""


def _verdicts_text(checks: CheckAttachment) -> str:
    """`check: verdict, verdict; ...`, each check's verdicts in clock order."""
    by_check: dict[str, list[str]] = {}
    for check, _dataset, verdict in checks.check_verdicts:
        by_check.setdefault(check, []).append(verdict)
    return "; ".join(f"{check}: {', '.join(v)}" for check, v in by_check.items())


def summaries_table(
    *records: RecordSummaries, excluded: tuple[str, ...]
) -> pd.DataFrame:
    """One row per record, rung and bag of `records`, then one per record in `excluded`
    (module docstring). Floats are kept as computed; a quantile the curve never reaches is
    NaN."""
    rows: list[dict[str, object]] = []
    for record in records:
        head = {"record": record.dataset_id, "drawn": True, "flag": record.flag}
        if record.failure:
            rows.append({**head, "reason": record.failure})
            continue
        for (label, side), item in record.summaries.items():
            lower, upper = record.rmst_interval_s[(label, side)]
            row: dict[str, object] = {
                **head,
                "reason": "",
                "threshold_label": label,
                "side": side,
                "conf_level": item.conf_level,
                "n_windows": item.n_windows,
                "n_deaths": item.n_deaths,
                "n_censored": item.n_censored,
                "n_unobserved_birth_dropped": item.n_unobserved_birth_dropped,
                "max_observed_age_s": item.max_observed_age_s,
                "rmst_tau_s": item.rmst_tau_s,
                "rmst_s": item.rmst_s,
                "rmst_se_s": item.rmst_se_s,
                "rmst_lower_s": lower,
                "rmst_upper_s": upper,
                "rmst_interval_note": _interval_note(item, lower),
            }
            for p, point, low, high in zip(
                item.quantile_p,
                item.quantile_s,
                item.quantile_lower_s,
                item.quantile_upper_s,
                strict=True,
            ):
                names = _quantile_columns(p)
                for name, value in zip(names, (point, low, high), strict=True):
                    row[name] = np.nan if value is None else value
            row.update(
                {
                    "point_standing": item.point_standing,
                    "interval_standing": item.interval_standing,
                    "assumption_id": item.checks.assumption_id,
                    "checks_asked": "; ".join(item.checks.checks_asked),
                    "checks_unanswered": "; ".join(item.checks.checks_unanswered),
                    "check_verdicts": _verdicts_text(item.checks),
                    "check_summary": item.checks.summary("interval"),
                }
            )
            rows.append(row)
    rows += [
        {"record": record, "drawn": False, "flag": "", "reason": EXCLUDED_REASON}
        for record in excluded
    ]
    unknown = sorted({key for row in rows for key in row} - set(SUMMARY_COLUMNS))
    if unknown:
        raise ValueError(f"summary columns {unknown} are not in SUMMARY_COLUMNS")
    table = pd.DataFrame(rows, columns=list(SUMMARY_COLUMNS))
    return table.astype({column: "Int64" for column in COUNT_COLUMNS})
