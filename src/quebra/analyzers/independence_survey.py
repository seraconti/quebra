"""Window independence across every dataset, one instrument at a time.

The check ledger answers "is this record's duration sequence iid?" for ONE dataset. This
reshapes many ledgers into one grid per instrument, bag and clock - datasets down,
thresholds across - so the question becomes "is the independence assumption safe on this
device at all, or only on the record we happened to look at first?" Both bags are surveyed,
the in-spec and the out-of-spec windows, each on both clocks.

**Coloured by VERDICT, not by p-value.** A p-value heatmap would be the obvious thing and
would be wrong: on a short window "did not reject" carries no information, and a pale cell
would read as evidence of independence when it is evidence of nothing. The ledger already
computes the three-way distinction - `pass`, `fail`, `inconclusive`/`not interpretable` -
against the calibration bench, so this reuses that verdict and prints the p inside the cell
as a secondary detail. That is the whole reason the survey is built on `check_ledger.run`
rather than on raw check results.

**One grid per entry of `SURVEY_KEYS`**, which is `battery.ROW_KEYS`: the five permutation
rows, each scored against a bench cell. The asymptotic rows of C1, C2 and CvM, and C5's
unstudentized variant, run only in the bench (`battery.BENCH_ROW_KEYS`). C3 is not
surveyed; the comment on `C3_KEY` below says why.

**Pure compute.** Ledger row tables in, typed artifact out. The job loads the datasets and
runs the ledgers; nothing here touches disk or matplotlib.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from quebra.analyzers.check_ledger import (
    VERDICT_FAIL,
    VERDICT_INCONCLUSIVE,
    VERDICT_NOT_COMPUTED,
    VERDICT_PASS,
    VERDICT_TIES,
)
from quebra.analyzers.checks.battery import ROW_KEYS
from quebra.analyzers.checks.result import CALIB_PERMUTATION
from quebra.analyzers.event_table import SIDES

# C3's key. It is NOT in `ROW_KEYS` and must not be: every row there is scored against a
# bench cell, and the bench does not run C3 - it is an out-of-process R call with its own
# simulated null, so `_verdict` could never return `pass` for it. The ledgers no longer run
# it and the survey draws no C3 grid; the key names C3's rows for a caller that asks a
# ledger for them.
C3_KEY = ("c3_serial_copula", "r_copula", "")

# What the survey draws: every row the pipeline's battery returns.
SURVEY_KEYS = tuple(ROW_KEYS)

# Display names, looked up by key rather than by position. Kept here rather than in the plot because they name
# the INSTRUMENT, which is a fact about the analysis and not about the axes.
CHECK_LABELS: dict[tuple[str, str, str], str] = {
    ("c1_lewis_robinson", "asymptotic", ""): "C1 Lewis-Robinson (asymptotic)",
    ("c1_lewis_robinson", "permutation", ""): "C1 Lewis-Robinson (permutation)",
    ("c2_anderson_darling", "asymptotic", ""): "C2 Anderson-Darling (asymptotic)",
    ("c2_anderson_darling", "permutation", ""): "C2 Anderson-Darling (permutation)",
    (
        "c5_rank_autocorr",
        "permutation",
        "studentized",
    ): "C5 rank autocorrelation (studentized)",
    ("c6_exchangeability", "permutation", ""): "C6 exchangeability",
    ("cvm_cramer_von_mises", "asymptotic", ""): "CvM Cramer-von Mises (asymptotic)",
    ("cvm_cramer_von_mises", "permutation", ""): "CvM Cramer-von Mises (permutation)",
    C3_KEY: "C3 serial copula (R, no bench cell)",
}

# What each rejection MEANS, so a reader does not have to remember that `fail` is a
# statement about the data and not about the test. FIGURE_STANDARD forbids a title that
# argues; this is a caption, and stating the null is description, not argument.
CHECK_NULL: dict[tuple[str, str, str], str] = {
    ("c1_lewis_robinson", "asymptotic", ""): "rejects when durations TREND over time",
    ("c1_lewis_robinson", "permutation", ""): "rejects when durations TREND over time",
    (
        "c2_anderson_darling",
        "asymptotic",
        "",
    ): "rejects when the process is NOT a renewal process",
    (
        "c2_anderson_darling",
        "permutation",
        "",
    ): "rejects when the process is NOT a renewal process",
    (
        "c5_rank_autocorr",
        "permutation",
        "studentized",
    ): "rejects when consecutive durations are RANK-CORRELATED",
    (
        "c6_exchangeability",
        "permutation",
        "",
    ): "rejects when the ORDER of the durations carries information",
    (
        "cvm_cramer_von_mises",
        "asymptotic",
        "",
    ): "rejects when the process is NOT a renewal process",
    (
        "cvm_cramer_von_mises",
        "permutation",
        "",
    ): "rejects when the process is NOT a renewal process",
    C3_KEY: "rejects when the durations are NOT serially independent at any lag",
}

VERDICT_ORDER = (
    VERDICT_PASS,
    VERDICT_FAIL,
    VERDICT_INCONCLUSIVE,
    VERDICT_TIES,
    VERDICT_NOT_COMPUTED,
)


@dataclass(frozen=True)
class InstrumentGrid:
    """One instrument's verdicts across every dataset and threshold, on one bag and clock."""

    key: tuple[str, str, str]
    label: str
    null_statement: str
    side: str
    clock: str
    verdicts: pd.DataFrame  # index = dataset, columns = threshold label
    p_values: pd.DataFrame
    n_events: pd.DataFrame

    @property
    def counts(self) -> dict[str, int]:
        flat = self.verdicts.to_numpy().ravel()
        return {v: int((flat == v).sum()) for v in VERDICT_ORDER}

    @property
    def rejection_share_of_decided(self) -> float:
        """`fail / (pass + fail)`. The denominator is DECIDED cells, not all cells.

        Dividing by every cell would let a mostly inconclusive grid look reassuring: 36
        inconclusive cells and 4 rejections is not "10% rejection", it is "4 of 4 decided cells
        rejected". NaN when nothing was decided, which is itself the finding.

        ALSO NaN when `pass` never occurred. For an instrument the bench has no cell for,
        `_verdict` sends every non-rejection to `inconclusive`, so `pass` is unreachable and
        this ratio is 1.0 by construction whenever anything rejected at all. Reporting that
        as a rejection SHARE would put a forced constant where a reader expects a measured
        rate: such a grid can read many fails and no pass, and 1.000 is not a finding about
        the data. The same degeneracy reaches any instrument the bench scored as
        miscalibrated at every event count.
        """
        c = self.counts
        decided = c[VERDICT_PASS] + c[VERDICT_FAIL]
        if decided == 0 or c[VERDICT_PASS] == 0:
            return float("nan")
        return c[VERDICT_FAIL] / decided


@dataclass(frozen=True)
class IndependenceSurveyData:
    """One grid per (instrument, bag, clock), plus what it took to build them."""

    grids: list[InstrumentGrid]
    datasets: list[str]
    thresholds: list[str]
    sides: list[str]
    clocks: list[str]
    dropped: dict[str, int] = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    def grid(
        self, key: tuple[str, str, str], side: str, clock: str
    ) -> InstrumentGrid | None:
        for g in self.grids:
            if g.key == key and g.side == side and g.clock == clock:
                return g
        return None


def _pivot(
    frame: pd.DataFrame, value: str, datasets: list[str], thresholds: list[str]
) -> pd.DataFrame:
    # One row per cell, or the pivot below would keep the first and drop the rest unseen.
    repeated = frame.duplicated(["dataset", "threshold_label"], keep=False)
    if repeated.any():
        cells = sorted(
            set(zip(frame["dataset"][repeated], frame["threshold_label"][repeated]))
        )
        raise ValueError(f"cells {cells} hold more than one row of one grid")
    out = frame.pivot_table(
        index="dataset", columns="threshold_label", values=value, aggfunc="first"
    )
    # reindex so every grid has identical axes even where a dataset produced no row
    return out.reindex(index=datasets, columns=thresholds)


def _pivot_verdicts(
    frame: pd.DataFrame, datasets: list[str], thresholds: list[str]
) -> pd.DataFrame:
    """Verdicts, with every remaining hole named rather than left as NaN.

    After the blanket rows are folded in, a hole means the BATTERY produced no row for this
    instrument here. Leaving those NaN made `counts` sum to
    less than the grid and let a reader take a blank cell for "nothing to say" when the
    truth is "this instrument does not apply to this record". `not computed` is that
    statement, and it is the tone the figure already uses.
    """
    return _pivot(frame, "verdict", datasets, thresholds).fillna(VERDICT_NOT_COMPUTED)


def build_independence_survey(
    *ledgers: object,
    dataset_labels: tuple[str, ...] = (),
    alpha: float = 0.05,
) -> IndependenceSurveyData:
    """Stack per-dataset `CheckLedger` artifacts into one grid per instrument, bag and
    clock. Raises ValueError on a row of no bag: every row of a window ledger has one.

    `dataset_labels` is passed as a step kwarg rather than read off the artifacts, so the
    row order of every grid is declared in the job and appears on the provenance label.
    Reading it off `dataset_id` would order the rows by whatever the loader happened to
    produce, and the figure's row order is part of how it is read.
    """
    if not ledgers:
        raise ValueError("the survey needs at least one ledger")
    if dataset_labels and len(dataset_labels) != len(ledgers):
        raise ValueError(
            f"{len(dataset_labels)} labels for {len(ledgers)} ledgers; the survey will "
            "not guess which is which"
        )
    # The grids pivot on the label, so a repeated label would draw its first record twice
    # and drop the other.
    repeated = sorted({x for x in dataset_labels if dataset_labels.count(x) > 1})
    if repeated:
        raise ValueError(
            f"dataset labels {repeated} repeat; each row must be one record"
        )

    frames = []
    for i, ledger in enumerate(ledgers):
        rows = getattr(ledger, "rows", None)
        if rows is None:
            raise TypeError(
                f"ledger {i} is a {type(ledger).__name__}, not a CheckLedger; the survey "
                "reshapes ledger rows and cannot take raw check results"
            )
        frame = rows.copy()
        frame["dataset"] = (
            dataset_labels[i]
            if dataset_labels
            else str(getattr(ledger, "dataset_id", "") or f"dataset {i}")
        )
        frames.append(frame)
    everything = pd.concat(frames, ignore_index=True)

    datasets = (
        list(dataset_labels)
        if dataset_labels
        else sorted(everything["dataset"].unique())
    )
    # Threshold order from the ladder as it appears, NOT sorted: "10 µs" sorts before
    # "2 µs" as a string and the ladder would be drawn out of order.
    thresholds = list(dict.fromkeys(everything["threshold_label"]))
    unknown = sorted({str(v) for v in everything["side"]} - set(SIDES))
    if unknown:
        raise ValueError(f"ledger rows of no bag: side {unknown}; bags: {SIDES}")
    present = set(everything["side"])
    sides = [side for side in SIDES if side in present]
    clocks = list(dict.fromkeys(everything["clock"]))

    # Rows that apply to EVERY check at a threshold. `check_ledger._blank_row` writes
    # `check_id="(all)"` when a rung produced no windows, or when segmentation declined it
    # - one row standing for every instrument. Matching on `check_id` alone therefore
    # drops them, and the cell renders as "no row" rather than as `not computed`: a rung
    # declined on every record would be blank in every per-instrument figure. A blank cell
    # reads as "nothing to say here"; `not computed` says the rung was reached and declined,
    # which is the truth and is a different statement.
    blanket = everything[everything["check_id"] == "(all)"]

    grids: list[InstrumentGrid] = []
    for key in SURVEY_KEYS:
        check, calibration, variant = key
        for side in sides:
            for clock in clocks:
                cell = (everything["side"] == side) & (everything["clock"] == clock)
                sub = everything[
                    cell
                    & (everything["check_id"] == check)
                    & (everything["calibration"] == calibration)
                    & (
                        everything.get("variant", pd.Series(dtype=str)).fillna("")
                        == variant
                    )
                ]
                # Fold in the blanket rows of this bag and clock, so a declined rung is
                # `not computed` in every instrument's grid rather than missing from all.
                sub = pd.concat([sub, blanket[cell[blanket.index]]], ignore_index=True)
                if sub.empty:
                    continue
                grids.append(
                    InstrumentGrid(
                        key=key,
                        label=CHECK_LABELS[key],
                        null_statement=CHECK_NULL[key],
                        side=side,
                        clock=clock,
                        verdicts=_pivot_verdicts(sub, datasets, thresholds),
                        p_values=_pivot(sub, "p_value", datasets, thresholds),
                        n_events=_pivot(sub, "n_events", datasets, thresholds),
                    )
                )

    return IndependenceSurveyData(
        grids=grids,
        datasets=datasets,
        thresholds=thresholds,
        sides=sides,
        clocks=clocks,
        dropped={},
        meta={
            "alpha": alpha,
            "n_datasets": len(datasets),
            "n_instruments": len(SURVEY_KEYS),
        },
    )


def survey_summary(data: IndependenceSurveyData) -> pd.DataFrame:
    """One row per (instrument, bag, clock): how many cells decided, and how many
    rejected."""
    rows = []
    for g in data.grids:
        c = g.counts
        rows.append(
            {
                "instrument": g.label,
                "side": g.side,
                "clock": g.clock,
                "pass": c[VERDICT_PASS],
                "fail": c[VERDICT_FAIL],
                "inconclusive": c[VERDICT_INCONCLUSIVE],
                "not_interpretable": c[VERDICT_TIES],
                "not_computed": c[VERDICT_NOT_COMPUTED],
                "rejection_share_of_decided": g.rejection_share_of_decided,
            }
        )
    return pd.DataFrame(rows)


def _fmt_p(value: object) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    v = float(value)
    return "<.001" if v < 0.001 else f"{v:.3f}".lstrip("0")


# --------------------------------------------------- the trend group, counted across ledgers

# The trend group: the three checks that read the same Brownian bridge of event times.
TREND_CHECKS = ("c1_lewis_robinson", "c2_anderson_darling", "cvm_cramer_von_mises")
_CELL = ["dataset_id", "threshold_label", "side", "clock"]
_TREND_COLUMNS = (*_CELL, "check_id", "calibration", "verdict", "notes")
# A number in a message, not a digit of a name such as `C1` or `T_N`.
_NUMBER = re.compile(r"(?<![A-Za-z_\d.])[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


@dataclass
class TrendGroupCounts:
    """Where the trend group answers, counted over many ledgers' cells.

    A cell is one (record, threshold, bag, clock). A check answers in a cell when its
    permutation row's verdict is not `not computed`; other calibrations are not read, so a
    ledger written while the asymptotic rows still ran counts the same way.

    `no_answer`: per bag, clock and message, the cells where neither C1 nor C2 answered,
    `n_cells`. The message is C1's note, or the cell's `(all)` row's when the cell has no
    C1 row (no windows, no segment), read as its first clause (up to the first `. ` or
    `, `), a note the ledger wrote twice as `X; X` read once, and numbers read `#`. The
    ledger truncates notes and prints numbers at varying widths, so the clause, not the
    whole note, is what groups the cells that failed for one reason.

    `cvm_alone`: per bag and clock, `n_cells`, `n_trend_answered` (any of the three answered),
    `n_cvm_alone` (CvM answered, C1 and C2 did not), and the share of CvM-alone cells among
    the answered cells and among all cells. A share is NaN when its denominator is 0.
    """

    no_answer: pd.DataFrame
    cvm_alone: pd.DataFrame


_CLAUSE_END = re.compile(r"(?<=\S)[.,] ")


def _message(note: object) -> str:
    parts = str(note).split("; ")
    half = len(parts) // 2
    if len(parts) % 2 == 0 and parts[:half] == parts[half:]:
        parts = parts[:half]
    clause = _CLAUSE_END.split("; ".join(parts), maxsplit=1)[0]
    text = " ".join(_NUMBER.sub("#", clause).split())
    return text or "(no message)"


def trend_group_counts(ledgers: Sequence[object]) -> TrendGroupCounts:
    """Count, over every cell of `ledgers`, where C1 and C2 gave no answer and where the
    trend group is CvM alone. Raises KeyError naming a missing row column, and ValueError
    on no ledger or no cell, on a record that appears in two ledgers, and on a cell where
    C1 and C2 were never asked: no row of theirs and no `(all)` row. That cell is not a
    silent C1 and C2, and counting it as one would read "not asked" as "no answer"."""
    if not ledgers:
        raise ValueError("no ledger to count")
    ids = [str(getattr(ledger, "dataset_id", "")) for ledger in ledgers]
    repeated = sorted({i for i in ids if ids.count(i) > 1})
    if repeated:
        raise ValueError(
            f"records {repeated} repeat across ledgers; each cell would be counted twice"
        )
    frame = pd.concat([ledger.rows for ledger in ledgers], ignore_index=True)
    missing = [c for c in _TREND_COLUMNS if c not in frame.columns]
    if missing:
        raise KeyError(f"the ledger rows have no column(s) {missing}")
    if not len(frame):
        raise ValueError("the ledgers hold no cell to count")

    cells = []
    for (_record, _label, side, clock), cell in frame.groupby(_CELL, sort=True):
        permutation = cell[cell["calibration"] == CALIB_PERMUTATION]

        def answered(check: str, permutation=permutation) -> bool:
            rows = permutation[permutation["check_id"] == check]
            return bool((rows["verdict"] != VERDICT_NOT_COMPUTED).any())

        asked = permutation["check_id"].isin(TREND_CHECKS[:2]).any()
        if not asked and not (cell["check_id"] == "(all)").any():
            raise ValueError(
                f"C1 and C2 were never asked in the cell {(_record, _label, side, clock)}: "
                "no row of theirs and no (all) row"
            )
        c1, c2, cvm = (answered(check) for check in TREND_CHECKS)
        source = permutation[permutation["check_id"] == TREND_CHECKS[0]]
        if not len(source):
            source = cell[cell["check_id"] == "(all)"]
        cells.append(
            {
                "side": side,
                "clock": clock,
                "no_c1_c2": not c1 and not c2,
                "message": _message(source["notes"].iloc[0]) if len(source) else "",
                "trend_answered": c1 or c2 or cvm,
                "cvm_alone": cvm and not c1 and not c2,
            }
        )
    table = pd.DataFrame(cells)
    no_answer = (
        table[table["no_c1_c2"]]
        .groupby(["side", "clock", "message"], sort=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    per_clock = table.groupby(["side", "clock"], sort=True).agg(
        n_cells=("clock", "size"),
        n_trend_answered=("trend_answered", "sum"),
        n_cvm_alone=("cvm_alone", "sum"),
    )
    answered_cells = per_clock["n_trend_answered"].where(
        per_clock["n_trend_answered"] > 0
    )
    per_clock["share_of_answered"] = per_clock["n_cvm_alone"] / answered_cells
    per_clock["share_of_all"] = per_clock["n_cvm_alone"] / per_clock["n_cells"]
    return TrendGroupCounts(no_answer=no_answer, cvm_alone=per_clock.reset_index())
