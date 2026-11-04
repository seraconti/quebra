"""Render-local helpers for CheckLedgerPanel: functions of axes and theme only.

Nothing here decides a verdict or derives a statistic. Every cell's colour comes from a
verdict string that `analyzers/check_ledger.py` already resolved, and every number is read
from the ledger frame as-is. These functions only work out where to put ink: which rows to
show, how to shorten a label to fit a cell, and how to lay the table out.

Kept out of `panels/check_ledger.py` so that file is the drawing sequence and this one is
the arithmetic the drawing needs.
"""

from __future__ import annotations

import pandas as pd

# How a check id reads in a table cell. The ledger keeps the full name; a heat-table column
# header has ~8 characters before it starts overlapping its neighbour.
CHECK_SHORT = {
    "c1_lewis_robinson": "C1 LR",
    "c2_anderson_darling": "C2 AD",
    "c3_serial_copula": "C3 cop",
    "c5_rank_autocorr": "C5 rank",
    "c6_exchangeability": "C6 exch",
    # `check_ledger` iterates these keys to decide which lines the p-value ladder draws, and
    # `statistic_series` below takes the minimum p-value over calibrations per rung. A
    # ledger carries one calibration per check or group, the permutation rows and a window
    # ledger's `min_p` group rows, so that minimum is the row's own p-value.
    "cvm_cramer_von_mises": "CvM",
    # The two group rows of a window ledger (`check_ledger.GROUPS`), beside the checks.
    "trend_group": "trend group",
    "dependence_group": "dependence group",
}


def column_label(check_id: str, calibration: str, variant: str) -> str:
    """One table column per (check, calibration, variant), shortened to fit."""
    base = CHECK_SHORT.get(check_id, check_id[:7])
    if calibration == "asymptotic":
        base += " a"
    elif calibration == "permutation":
        base += " p"
    if isinstance(variant, str) and variant.startswith("unstud"):
        base += " (u)"
    return base


def ledger_grid(
    rows: pd.DataFrame, clock: str
) -> tuple[list[str], list[str], list[list]]:
    """Pivot the ledger into (row labels, column labels, cell records) for one clock.

    Cell records are `(verdict, p_value, bench_accepted, calibration)` so the renderer can
    colour by the first, annotate with the second and mark the third without going back to
    the frame. `bench_accepted` is None where the bench has no cell for that check at that
    event count, which is a different thing from False; a group row has none by
    construction, which its calibration (`min_p`) says.
    """
    frame = rows[rows["clock"] == clock].copy()
    if not len(frame):
        return [], [], []
    frame["column"] = [
        column_label(c, k, v)
        for c, k, v in zip(frame["check_id"], frame["calibration"], frame["variant"])
    ]
    frame = frame[frame["check_id"] != "(all)"]
    thresholds = list(dict.fromkeys(frame["threshold_label"]))
    columns = list(dict.fromkeys(frame["column"]))
    lookup = {
        (r["threshold_label"], r["column"]): (
            r["verdict"],
            r["p_value"],
            r["bench_accepted"],
            r["calibration"],
        )
        for _i, r in frame.iterrows()
    }
    grid = [[lookup.get((t, c)) for c in columns] for t in thresholds]
    return thresholds, columns, grid


def statistic_series(
    rows: pd.DataFrame, clock: str, check_id: str, ladder: list[str]
) -> list[float]:
    """One check's p-value at every rung of the ladder, NaN where it has none.

    NaN rather than omission is the point: matplotlib breaks a line at NaN, and
    `docs/FIGURE_STANDARD.md` requires that nothing be drawn across a gap. Filtering the
    missing rungs out instead drew C1 straight from 3 us to 7 us across three thresholds
    where it is undefined, implying evidence that does not exist.
    """
    frame = rows[(rows["clock"] == clock) & (rows["check_id"] == check_id)]
    best = (
        frame.dropna(subset=["p_value"])
        .groupby("threshold_label")["p_value"]
        .min()
        .to_dict()
    )
    return [float(best.get(label, float("nan"))) for label in ladder]


def format_p(value: object) -> str:
    """A p-value as it should appear in a table cell.

    Permutation p-values have a floor of `1/(B+1)`, so "0.001" at B = 999 means "the
    smallest value this test can report", not "0.001 exactly". Writing it as `<.001`
    says so without a footnote.
    """
    if value is None or pd.isna(value):
        return "-"
    number = float(value)
    if number <= 0.001:
        return "<.001"
    return f"{number:.3f}".lstrip("0")
