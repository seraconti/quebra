"""Window independence across every dataset, one instrument at a time.

The check ledger answers "is this record's duration sequence iid?" for ONE dataset. This
reshapes many ledgers into one grid per instrument - datasets down, thresholds across - so
the question becomes "is the independence assumption safe on this device at all, or only
on the record we happened to look at first?"

**Coloured by VERDICT, not by p-value.** A p-value heatmap would be the obvious thing and
would be wrong: on a short window "did not reject" carries no information, and a pale cell
would read as evidence of independence when it is evidence of nothing. The ledger already
computes the three-way distinction - `pass`, `fail`, `underpowered`/`not interpretable` -
against the calibration bench, so this reuses that verdict and prints the p inside the cell
as a secondary detail. That is the whole reason the survey is built on `check_ledger.run`
rather than on raw check results.

**One grid per entry of `SURVEY_KEYS`**, which is `battery.ROW_KEYS` plus C3. C1, C2 and
CvM each appear twice, once per calibration, because the asymptotic and permutation routes
are different instruments that happen to share a statistic - the tier work measured them
separately and they disagree. C5 appears twice for its two variants. C3 is drawn but is
the one instrument with no bench cell behind it, and the comment on `C3_KEY` below carries
what that costs a reader.

**Pure compute.** Ledger row tables in, typed artifact out. The job loads the datasets and
runs the ledgers; nothing here touches disk or matplotlib.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from quebra.analyzers.check_ledger import (
    VERDICT_FAIL,
    VERDICT_NOT_COMPUTED,
    VERDICT_PASS,
    VERDICT_TIES,
    VERDICT_UNDERPOWERED,
)
from quebra.analyzers.checks.battery import ROW_KEYS

# C3's key. It is NOT in `ROW_KEYS` and must not be: that tuple is the bench tables' row
# schema, and the bench does not run C3 - it is an out-of-process R call with its own
# simulated null. But the SURVEY draws it, because a survey that runs an instrument and
# then omits it from every figure answers a different question from the one it appears to.
# The consequence travels with it: C3 is the one grid here with no bench cell behind it, so
# its `pass` means "did not reject" and cannot mean "and it had the power to".
C3_KEY = ("c3_serial_copula", "r_copula", "")

# What the survey draws: everything the bench scores, plus C3.
SURVEY_KEYS = tuple(ROW_KEYS) + (C3_KEY,)

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
    (
        "c5_rank_autocorr",
        "permutation",
        "unstudentized",
    ): "C5 rank autocorrelation (raw)",
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
        "c5_rank_autocorr",
        "permutation",
        "unstudentized",
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
    VERDICT_UNDERPOWERED,
    VERDICT_TIES,
    VERDICT_NOT_COMPUTED,
)


@dataclass(frozen=True)
class InstrumentGrid:
    """One instrument's verdicts across every dataset and threshold, per clock."""

    key: tuple[str, str, str]
    label: str
    null_statement: str
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

        Dividing by every cell would let an underpowered grid look reassuring: 40 cells with
        no power and 4 rejections is not "10% rejection", it is "4 of 4 decided cells
        rejected". NaN when nothing was decided, which is itself the finding.

        ALSO NaN when `pass` never occurred. For an instrument the bench has no cell for,
        `_verdict` sends every non-rejection to `underpowered`, so `pass` is unreachable and
        this ratio is 1.0 by construction whenever anything rejected at all. Reporting that
        as a rejection SHARE would put a forced constant where a reader expects a measured
        rate: the shipped C3 grid reads 46 fail, 0 pass, and 1.000 is not a finding about
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
    """One grid per (instrument, clock), plus what it took to build them."""

    grids: list[InstrumentGrid]
    datasets: list[str]
    thresholds: list[str]
    clocks: list[str]
    dropped: dict[str, int] = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    def grids_for_clock(self, clock: str) -> list[InstrumentGrid]:
        return [g for g in self.grids if g.clock == clock]

    def grid(self, key: tuple[str, str, str], clock: str) -> InstrumentGrid | None:
        for g in self.grids:
            if g.key == key and g.clock == clock:
                return g
        return None


def _pivot(
    frame: pd.DataFrame, value: str, datasets: list[str], thresholds: list[str]
) -> pd.DataFrame:
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
    instrument here - `run_battery` drops C2-asymptotic when a record has more than one
    segment, which is 35 cells on the shipped data. Leaving those NaN made `counts` sum to
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
    """Stack per-dataset `CheckLedger` artifacts into one grid per instrument and clock.

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
        for clock in clocks:
            sub = everything[
                (everything["check_id"] == check)
                & (everything["calibration"] == calibration)
                & (everything["clock"] == clock)
                & (
                    everything.get("variant", pd.Series(dtype=str)).fillna("")
                    == variant
                )
            ]
            # Fold in the blanket rows for this clock, so a declined rung is `not
            # computed` in every instrument's grid rather than missing from all of them.
            sub = pd.concat(
                [sub, blanket[blanket["clock"] == clock]], ignore_index=True
            )
            if sub.empty:
                continue
            grids.append(
                InstrumentGrid(
                    key=key,
                    label=CHECK_LABELS[key],
                    null_statement=CHECK_NULL[key],
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
        clocks=clocks,
        dropped={},
        meta={
            "alpha": alpha,
            "n_datasets": len(datasets),
            "n_instruments": len(SURVEY_KEYS),
        },
    )


def survey_summary(data: IndependenceSurveyData) -> pd.DataFrame:
    """One row per (instrument, clock): how many cells decided, and how many rejected."""
    rows = []
    for g in data.grids:
        c = g.counts
        rows.append(
            {
                "instrument": g.label,
                "clock": g.clock,
                "pass": c[VERDICT_PASS],
                "fail": c[VERDICT_FAIL],
                "underpowered": c[VERDICT_UNDERPOWERED],
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
