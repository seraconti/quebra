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
"""

from __future__ import annotations

from dataclasses import dataclass

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
from quebra.analyzers.windows import WindowsResult
from quebra.core._artifact_guard import StaleArtifactGuard
from quebra.core.types import RecordFailure

# The check id of the one row the ledger writes for a rung, bag and clock it could not score.
ALL_ROWS = "(all)"
# The verdicts a survival block's key names: what the ledger can say, then what a builder
# writes for a cell with no row.
VERDICT_KEY = (*VERDICTS, VERDICT_ABSENT)


def ledger_run_set(*, include_c3: bool) -> tuple[RowKey, ...]:
    """The run-set a ledger ran, from its node's kwargs: the battery's permutation rows.

    A ledger run with `include_c3` holds C3 rows, which no run-set can select
    (`check_outcome.ALL_KEYS`), so it raises.
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
    rows = expand_all_rows(ledger.rows, run)
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
    outcomes = {
        (label, side): build_check_outcome(
            rows,
            threshold_label=label,
            run_set=run,
            display_set=run,
            datasets=(record,),
            dataset_labels=(f"{record}, {side} bag",),
            clocks=CLOCKS,
            alpha=ledger.alpha,
            side=side,
        )
        for label in ladder
        for side in SIDES
    }
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
