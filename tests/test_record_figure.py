"""Record figures: the carve held to its ledger, the survival and resolvability blocks, and
the figure jobs.

Oracles: hand tables for the cross-check; the real `ledger_step` on a synthetic record's
real carve for the run-set, the `(all)` expansion and each band's check outcome; the real
tellability of that record for the resolvability block, with G10 recomputed here from its
definition; the artifact's own fields read back off the drawn figure; the job graphs for the
one-ref rule, the record list and the one tellability node per record; and a run directory
for the nested-run count.
"""

from __future__ import annotations

import copy
import dataclasses
import io
import json
import pickle
import warnings
from pathlib import Path
from statistics import NormalDist
from types import SimpleNamespace
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import to_rgba
from matplotlib.text import Text

from quebra.analyzers.check_ledger import (
    CALIB_MIN_P,
    VERDICT_ABSENT,
    VERDICT_NOT_COMPUTED,
    CheckLedger,
)
from quebra.analyzers.check_outcome import PERMUTATION_KEYS
from quebra.analyzers.checks.result import CLOCKS
from quebra.analyzers.event_table import SIDES
from quebra.analyzers.kaplan_meier import CONF_LEVEL
from quebra.analyzers.event_table import EventTable, event_tables_from_carve
from quebra.analyzers.independence_survey import CHECK_LABELS
from quebra.analyzers.compliance_survey import record_compliance
from quebra.analyzers.record_figure import assemble_record_figure
from quebra.analyzers.record_figure import (
    ALL_ROWS,
    COUNT_COLUMNS,
    EXCLUDED_REASON,
    NOTE_LOWER_AT_ZERO,
    NOTE_NO_WIDTH,
    SUMMARY_COLUMNS,
    RecordSummaries,
    FractionBlockData,
    TimelineBlockData,
    build_fraction_block,
    build_hazard_block,
    build_record_summaries,
    LIMIT_INDEPENDENT,
    LIMIT_SMOOTH,
    build_resolvability_block,
    build_survival_block,
    cross_check_carve,
    curve_summaries_of,
    kaplan_meier_of,
    ledger_run_set,
    nelson_aalen_of,
    panel_data_of,
    summaries_table,
)
from quebra.analyzers.record_figure import build_timeline_block
from quebra.analyzers.survival import (
    curve_summaries,
    curve_summaries_set,
    nelson_aalen_set,
)
from quebra.core.job import Job, _import_job
from quebra.core.reference import ArtifactRef
from quebra.core.runner import run_job
from quebra.analyzers.tellability import METRIC_GROUPS, THRESHOLD_GROUPS
from quebra.core.types import RecordFailure
from quebra.ledger_recipe import LEDGER_KNOBS, ledger_step
from quebra.paper_tables import FLAG_TEXT, PAPER_RECORDS, SHADED
from quebra.plots.record_figure_plot import (
    BOUNDS_GID,
    CELL_GID,
    CROSSING_GID,
    FRACTION_GID,
    FractionBlockPlot,
    TimelineBlockPlot,
    HazardBlockPlot,
    LIMIT_DASHES,
    NOT_COMPUTED_GID,
    NO_BOUNDS_GID,
    RecordFigurePlot,
    ResolvabilityBlockPlot,
    SummariesTablePlot,
    SurvivalBlockPlot,
    ZERO_RATIO_GID,
)
from quebra.panels.within_calibration import WithinCalibrationPanel
from quebra.plots import record_figure_plot
from quebra.plots.targets import render_table
from quebra.plots.theme import state_color
from quebra.plots.theme import verdict_color
from quebra.record_figure_recipe import (
    BLOCKS,
    PAGE_ORDER,
    RecordNodes,
    check_tellability_kwargs,
    configure_record_figure_job,
    t2star_tellability_of,
)
from quebra.recipes import T2STAR_THRESHOLDS, XI_SEED
from quebra.tellability_recipe import TELLABILITY_PERMUTATIONS, TELLABILITY_SEED
from tests.test_compliance_survey import _carve, _track912_csv

REPO = Path(__file__).resolve().parents[1]
BENCH = REPO / "jobs" / "bench" / "results" / "size_table.csv"
RECORD = "q1_020423"
FIGURE_JOBS = sorted((REPO / "jobs" / "composite").glob("record_figures_6d2s_q*.py"))


@pytest.fixture(scope="module")
def scored(tmp_path_factory) -> SimpleNamespace:
    """A synthetic record's real T2* fit, its carve and the real ledger of it, at 99
    permutations."""
    root = tmp_path_factory.mktemp("record_figure")
    fit, carve = _carve(
        _track912_csv(root / "020423_r.csv", seed=2, gap=True, n=600), RECORD
    )
    ledger = ledger_step(
        carve,
        pd.read_csv(BENCH),
        thresholds=T2STAR_THRESHOLDS,
        **{**LEDGER_KNOBS, "n_permutations": 99},  # type: ignore[arg-type]
    )
    return SimpleNamespace(fit=fit, carve=carve, ledger=ledger)


def _block(km: object, ledger: CheckLedger, flag: str = "") -> object:
    return build_survival_block(
        km,  # type: ignore[arg-type]
        ledger,
        run_set=ledger_run_set(include_c3=False),
        conf_level=CONF_LEVEL,
        flag=flag,
    )


def _hand_ledger(rows: list[tuple[str, str, str, int, str]]) -> CheckLedger:
    """A ledger of (rung, side, check id, windows used, note) rows on the side clock."""
    return CheckLedger(
        rows=pd.DataFrame(
            [
                {
                    "dataset_id": "hand",
                    "threshold_label": rung,
                    "side": side,
                    "clock": "side",
                    "check_id": check,
                    "n_windows_used": n,
                    "notes": note,
                }
                for rung, side, check, n, note in rows
            ]
        ),
        dataset_id="hand",
        thresholds=[("a", 1.0, True), ("b", 2.0, True)],
    )


def _hand_carve(counts: dict[tuple[str, str], int]) -> SimpleNamespace:
    """A carve of `counts[(rung, side)]` windows per rung and side."""
    frames = {
        side: pd.DataFrame(
            {
                "threshold_label": [
                    r for (r, s), n in counts.items() if s == side for _ in range(n)
                ],
                "side": [
                    side for (_r, s), n in counts.items() if s == side for _ in range(n)
                ],
            }
        )
        for side in SIDES
    }
    return SimpleNamespace(
        meta={"dataset_id": "hand"},
        windows_in_spec=frames[SIDES[0]],
        windows_out_of_spec=frames[SIDES[1]],
    )


@pytest.mark.unit
def test_the_cross_check_holds_a_carve_to_the_ledger_that_scored_it() -> None:
    """Oracle: hand tables. The ledger used 3 and 2 windows at rung a and 0 and 1 at rung
    b, the out-of-spec bag of a in one `(all)` row; the carve that holds those counts
    passes, and one window more at any cell, a ledger of another record, or a failed carve
    against a scored ledger raises. A failed carve passes against a ledger of its failure."""
    in_spec, out_of_spec = SIDES
    ledger = _hand_ledger(
        [
            ("a", in_spec, "c1", 3, ""),
            ("a", in_spec, "c5", 3, ""),
            ("a", out_of_spec, ALL_ROWS, 2, "too few events"),
            ("b", in_spec, ALL_ROWS, 0, "no windows at this threshold"),
            ("b", out_of_spec, "c1", 1, ""),
        ]
    )
    counts = {
        ("a", in_spec): 3,
        ("a", out_of_spec): 2,
        ("b", in_spec): 0,
        ("b", out_of_spec): 1,
    }
    cross_check_carve(_hand_carve(counts), ledger)  # type: ignore[arg-type]
    for cell in counts:
        with pytest.raises(ValueError, match="ledger rows used"):
            cross_check_carve(_hand_carve({**counts, cell: counts[cell] + 1}), ledger)  # type: ignore[arg-type]
    other = _hand_carve(counts)
    other.meta = {"dataset_id": "other"}
    with pytest.raises(ValueError, match="is checked against the ledger of 'hand'"):
        cross_check_carve(other, ledger)  # type: ignore[arg-type]
    failure = RecordFailure(dataset_id="hand", stage="t2star", reason="ValueError: x")
    with pytest.raises(ValueError, match="the carve failed"):
        cross_check_carve(failure, ledger)
    failed = _hand_ledger([("a", in_spec, ALL_ROWS, 0, failure.note())])
    cross_check_carve(failure, failed)


@pytest.mark.unit
def test_the_run_set_holds_every_check_the_ledger_writes(scored) -> None:
    """Oracle: the real `ledger_step`. Every row of the synthetic record's ledger but an
    `(all)` row or a group row has a key of the derived run-set, and the run-set's every
    key has a row; a ledger run with C3 raises."""
    rows = scored.ledger.rows
    answered = rows[
        (rows["check_id"] != ALL_ROWS) & (rows["calibration"] != CALIB_MIN_P)
    ]
    keys = {
        (c, k, "" if pd.isna(v) else v)
        for c, k, v in answered[["check_id", "calibration", "variant"]].itertuples(
            index=False
        )
    }
    assert keys == set(ledger_run_set(include_c3=False)) == set(PERMUTATION_KEYS)
    with pytest.raises(ValueError, match="ran C3"):
        ledger_run_set(include_c3=True)


@pytest.mark.unit
def test_each_band_carries_the_check_outcome_of_its_own_rung_and_bag(scored) -> None:
    """Oracle: the real ledger's rows. Every cell of every band's outcome is the ledger row
    of that rung, bag, clock and check; the attachment lists those cells under a label
    naming the bag, the grid's verdicts in order; and the two bags differ somewhere, so a
    block giving both bands one bag's rows cannot pass."""
    block = _block(kaplan_meier_of(scored.carve, conf_level=CONF_LEVEL), scored.ledger)
    rows = scored.ledger.rows
    differ = False
    for (rung, side), outcome in block.outcomes.items():
        for grid in outcome.grids:
            check, calibration, variant = grid.key
            for clock, verdict in zip(CLOCKS, grid.verdicts[0], strict=True):
                cell = rows[
                    (rows["threshold_label"] == rung)
                    & (rows["side"] == side)
                    & (rows["clock"] == clock)
                    & (rows["check_id"].isin([check, ALL_ROWS]))
                ]
                assert len(cell) == 1, (rung, side, clock, grid.key)
                assert verdict == cell["verdict"].iloc[0]
        label = f"{RECORD}, {side} bag"
        assert block.attachments[(rung, side)].check_verdicts == tuple(
            (grid.label, label, v) for grid in outcome.grids for v in grid.verdicts[0]
        )
        other = block.outcomes[(rung, SIDES[1 - SIDES.index(side)])]
        differ |= [g.verdicts for g in outcome.grids] != [
            g.verdicts for g in other.grids
        ]
    assert differ, "premise: the bags' outcomes differ somewhere"
    assert set(block.ladder) == {label for label, _, _ in T2STAR_THRESHOLDS}


@pytest.mark.unit
def test_an_all_row_reaches_every_check_of_its_cell(scored) -> None:
    """Oracle: a hand edit of the real ledger. One cell's five rows replaced by one `(all)`
    row reads, in every grid of its band, that row's verdict and note, never `no row`; the
    other clock of that band keeps its own rows."""
    rows = scored.ledger.rows
    rung, side, clock = T2STAR_THRESHOLDS[2][0], SIDES[0], CLOCKS[0]
    cell = (
        (rows["threshold_label"] == rung)
        & (rows["side"] == side)
        & (rows["clock"] == clock)
    )
    blank = {**rows[cell].iloc[0].to_dict(), "check_id": ALL_ROWS, "calibration": "",
             "variant": "", "verdict": VERDICT_NOT_COMPUTED, "notes": "hand note"}  # fmt: skip
    edited = pd.concat([rows[~cell], pd.DataFrame([blank])], ignore_index=True)
    ledger = CheckLedger(
        rows=edited, dataset_id=RECORD, thresholds=list(T2STAR_THRESHOLDS)
    )
    block = _block(kaplan_meier_of(scored.carve, conf_level=CONF_LEVEL), ledger)
    outcome = block.outcomes[(rung, side)]
    assert len(outcome.grids) == len(PERMUTATION_KEYS)
    for grid in outcome.grids:
        assert grid.verdicts[0][0] == VERDICT_NOT_COMPUTED
        assert grid.notes[0][0] == "hand note"
        assert grid.verdicts[0][1] != VERDICT_ABSENT
    assert outcome.asked_but_unanswered == ()


@pytest.mark.unit
def test_the_survival_block_draws_its_artifact(scored) -> None:
    """Oracle: the artifact's fields read back off the figure. One curve axes per rung and
    bag, titled by rung, its counts under its age label; the check cells in each band's
    verdict colours, bag by bag and rung by rung; the flag in the caption."""
    block = _block(
        kaplan_meier_of(scored.carve, conf_level=CONF_LEVEL), scored.ledger, FLAG_TEXT
    )
    fig = SurvivalBlockPlot("survival").build_matplotlib(block)
    try:
        curves = [ax for ax in fig.axes if ax.get_xscale() == "log"]
        assert len(curves) == len(SIDES) * len(block.ladder)
        expected_cells = []
        for ax, (side, rung) in zip(
            curves, [(s, r) for s in SIDES for r in block.ladder], strict=True
        ):
            curve = block.curves[(rung, side)]
            assert ax.get_title() == rung
            assert (
                f"n {curve.n_windows}, d {curve.n_deaths}, c {curve.n_censored}"
                in ax.get_xlabel()
            )
            assert f"excluded {curve.n_unobserved_birth_dropped}" in ax.get_xlabel()
            (line,) = [
                line for line in ax.get_lines() if line.get_linestyle() != "None"
            ]
            assert np.max(line.get_xdata()) == pytest.approx(curve.max_observed_min)
            expected_cells += [
                to_rgba(verdict_color(v))
                for grid in block.outcomes[(rung, side)].grids
                for v in grid.verdicts[0]
            ]
        drawn_cells = [
            p.get_facecolor()
            for ax in fig.axes
            for p in ax.patches
            if p.get_gid() == CELL_GID
        ]
        assert drawn_cells == expected_cells
        drawn = "\n".join(t.get_text() for t in fig.findobj(Text))
        assert FLAG_TEXT in drawn
        assert f"Window survival by rung, {RECORD}" in drawn
    finally:
        plt.close(fig)


@pytest.mark.unit
def test_a_failed_record_draws_its_failure_and_no_curve() -> None:
    """Oracle: the failure's own note. A record whose carve failed has no curve; its bands'
    outcomes read `not computed` with the note, and the figure draws the note."""
    failure = RecordFailure(
        dataset_id=RECORD, stage="t2star", reason="ValueError: no reads"
    )
    ledger = ledger_step(
        failure,
        pd.read_csv(BENCH),
        thresholds=T2STAR_THRESHOLDS,
        **{**LEDGER_KNOBS, "n_permutations": 99},  # type: ignore[arg-type]
    )
    cross_check_carve(failure, ledger)
    block = _block(kaplan_meier_of(failure, conf_level=CONF_LEVEL), ledger)
    assert block.failure == failure.note() and block.curves == {}
    for outcome in block.outcomes.values():
        assert {v for g in outcome.grids for v in g.verdicts[0]} == {
            VERDICT_NOT_COMPUTED
        }
    fig = SurvivalBlockPlot("survival").build_matplotlib(block)
    try:
        drawn = "\n".join(t.get_text() for t in fig.findobj(Text))
        assert f"not computed: {failure.note()}" in drawn
    finally:
        plt.close(fig)


def _hazard_block(carve: object, ledger: CheckLedger, flag: str = "") -> object:
    return build_hazard_block(
        nelson_aalen_of(carve, conf_level=CONF_LEVEL),  # type: ignore[arg-type]
        ledger,
        run_set=ledger_run_set(include_c3=False),
        conf_level=CONF_LEVEL,
        flag=flag,
    )


@pytest.mark.unit
def test_each_hazard_item_carries_its_own_rung_and_bag_outcome(scored) -> None:
    """Oracles: `survival.nelson_aalen_set` called on its own, the raw ledger rows, and
    minutes = seconds / 60. Each item keeps the reference's numbers, which carry no check
    outcome; its `checks` field names a1 and lists, check by check and clock by clock, the
    ledger's verdict of its own rung and bag."""
    block = _hazard_block(scored.carve, scored.ledger)
    reference = nelson_aalen_set(
        event_tables_from_carve(scored.carve), conf_level=CONF_LEVEL
    )
    rows = scored.ledger.rows
    assert set(block.hazards) == set(reference.items)
    for (rung, side), item in block.hazards.items():
        plain = reference.get(rung, side)
        assert plain.checks.checks_asked == (), "premise: the step alone fills nothing"
        for field in (
            "time_s",
            "cumulative_hazard",
            "variance",
            "band_lower",
            "band_upper",
        ):
            np.testing.assert_array_equal(getattr(item, field), getattr(plain, field))
        assert (item.n_windows, item.n_deaths, item.n_unobserved_birth_dropped) == (
            plain.n_windows, plain.n_deaths, plain.n_unobserved_birth_dropped
        )  # fmt: skip
        np.testing.assert_allclose(block.time_min[(rung, side)], plain.time_s / 60.0)
        expected = []
        for key in PERMUTATION_KEYS:
            for clock in CLOCKS:
                cell = rows[
                    (rows["threshold_label"] == rung)
                    & (rows["side"] == side)
                    & (rows["clock"] == clock)
                    & (rows["check_id"].isin([key[0], ALL_ROWS]))
                ]
                expected.append(
                    (
                        CHECK_LABELS[key],
                        f"{RECORD}, {side} bag",
                        cell["verdict"].iloc[0],
                    )
                )
        assert item.checks.assumption_id == "a1_renewal_durations"
        assert item.checks.check_verdicts == tuple(expected), (rung, side)


@pytest.mark.unit
def test_the_hazard_block_draws_its_artifact(scored) -> None:
    """Oracle: the artifact's fields read back off the figure. One log-log axes per bag
    holding each rung's H after its first death, in ladder order; every rung's counts in
    that bag's legend; the check cells in each band's verdict colours, bag by bag, rung by
    rung, check by check; the flag in the caption. A failed record draws its note."""
    block = _hazard_block(scored.carve, scored.ledger, FLAG_TEXT)
    fig = HazardBlockPlot("hazard").build_matplotlib(block)
    try:
        axes = [ax for ax in fig.axes if ax.get_yscale() == "log"]
        assert len(axes) == len(SIDES)
        drawn = "\n".join(t.get_text() for t in fig.findobj(Text))
        expected_cells = []
        for ax, side in zip(axes, SIDES, strict=True):
            assert ax.get_title(loc="left").startswith(side.split("_")[0])
            steps = [line for line in ax.get_lines() if len(line.get_ydata())]
            items = [block.hazards[(r, side)] for r in block.ladder]
            with_deaths = [item for item in items if item.n_deaths]
            assert len(steps) == len(with_deaths)
            for line, item in zip(steps, with_deaths, strict=True):
                np.testing.assert_array_equal(
                    line.get_ydata(), item.cumulative_hazard[1:]
                )
            for rung, item in zip(block.ladder, items, strict=True):
                assert (
                    f"{rung}: n {item.n_windows}, d {item.n_deaths}, "
                    f"excluded {item.n_unobserved_birth_dropped}"
                ) in drawn
                expected_cells += [
                    to_rgba(verdict_color(v))
                    for grid in block.outcomes[(rung, side)].grids
                    for v in grid.verdicts[0]
                ]
        cells = [
            p.get_facecolor()
            for ax in fig.axes
            for p in ax.patches
            if p.get_gid() == CELL_GID
        ]
        assert cells == expected_cells
        assert FLAG_TEXT in drawn
    finally:
        plt.close(fig)
    failure = RecordFailure(
        dataset_id=RECORD, stage="t2star", reason="ValueError: no reads"
    )
    failed = _hazard_block(
        failure,
        ledger_step(
            failure,
            pd.read_csv(BENCH),
            thresholds=T2STAR_THRESHOLDS,
            **{**LEDGER_KNOBS, "n_permutations": 99},  # type: ignore[arg-type]
        ),
    )
    assert failed.hazards == {} and failed.failure == failure.note()
    fig = HazardBlockPlot("hazard").build_matplotlib(failed)
    try:
        assert f"not computed: {failure.note()}" in "\n".join(
            t.get_text() for t in fig.findobj(Text)
        )
    finally:
        plt.close(fig)


def _summaries(scored: SimpleNamespace, level: float, flag: str = "") -> object:
    return build_record_summaries(
        curve_summaries_of(scored.carve, conf_level=level),
        scored.ledger,
        run_set=ledger_run_set(include_c3=False),
        conf_level=level,
        flag=flag,
    )


@pytest.mark.unit
def test_each_summaries_row_is_its_items_fields(scored) -> None:
    """Oracles: `survival.curve_summaries_set` called on its own, the interval in closed
    form at 95% and at 90% (z from the standard library), and the raw ledger rows. One row
    per rung and bag holding every number of the item unrounded, the restricted mean's
    interval at the block's level with a note where it has no width or reaches 0, and the
    band's check outcome, each check's verdicts on both clocks; an excluded record one row,
    not drawn, with the reason; a failed record one row with its note."""
    rows = scored.ledger.rows
    failure = RecordFailure(
        dataset_id="q1_030423", stage="t2star", reason="ValueError: x"
    )
    failed = build_record_summaries(
        failure,
        ledger_step(
            failure,
            pd.read_csv(BENCH),
            thresholds=T2STAR_THRESHOLDS,
            **{**LEDGER_KNOBS, "n_permutations": 99},  # type: ignore[arg-type]
        ),
        run_set=ledger_run_set(include_c3=False),
        conf_level=CONF_LEVEL,
        flag="",
    )
    notes: set[str] = set()
    for level in (0.95, 0.90):
        z = NormalDist().inv_cdf(0.5 + level / 2)
        table = summaries_table(
            _summaries(scored, level), failed, excluded=("q1_040423",)
        )
        assert list(table.columns) == list(SUMMARY_COLUMNS)
        assert {str(table[c].dtype) for c in COUNT_COLUMNS} == {"Int64"}
        reference = curve_summaries_set(
            event_tables_from_carve(scored.carve), conf_level=level
        )
        drawn = table[table["record"] == RECORD]
        assert len(drawn) == len(reference.items)
        for row in drawn.to_dict("records"):
            plain = reference.get(row["threshold_label"], row["side"])
            assert row["drawn"] and row["reason"] == "" and row["conf_level"] == level
            for name in (*COUNT_COLUMNS, "max_observed_age_s", "rmst_tau_s", "rmst_s",
                         "rmst_se_s"):  # fmt: skip
                np.testing.assert_array_equal(
                    row[name], getattr(plain, name), err_msg=name
                )
            for i, p in enumerate(plain.quantile_p):
                q = f"q{100 * p:g}"
                for column, values in (
                    (f"{q}_s", plain.quantile_s),
                    (f"{q}_lower_s", plain.quantile_lower_s),
                    (f"{q}_upper_s", plain.quantile_upper_s),
                ):
                    expected = np.nan if values[i] is None else values[i]
                    np.testing.assert_array_equal(row[column], expected, err_msg=column)
            lower = plain.rmst_s - z * plain.rmst_se_s
            upper = plain.rmst_s + z * plain.rmst_se_s
            if plain.n_windows:
                assert row["rmst_lower_s"] == pytest.approx(lower, rel=1e-12, abs=1e-9)
                assert row["rmst_upper_s"] == pytest.approx(upper, rel=1e-12, abs=1e-9)
            note = row["rmst_interval_note"]
            notes.add(note)
            if plain.rmst_se_s == 0.0:
                assert (
                    note == NOTE_NO_WIDTH and row["rmst_lower_s"] == row["rmst_upper_s"]
                )
            elif plain.n_windows and lower <= 0.0:
                assert note == NOTE_LOWER_AT_ZERO
            else:
                assert note == ""
            cells = []
            for key in PERMUTATION_KEYS:
                verdicts = []
                for clock in CLOCKS:
                    cell = rows[
                        (rows["threshold_label"] == row["threshold_label"])
                        & (rows["side"] == row["side"])
                        & (rows["clock"] == clock)
                        & (rows["check_id"].isin([key[0], ALL_ROWS]))
                    ]
                    verdicts.append(str(cell["verdict"].iloc[0]))
                cells.append(f"{CHECK_LABELS[key]}: {', '.join(verdicts)}")
            assert row["check_verdicts"] == "; ".join(cells)
            assert row["assumption_id"] == "a1_renewal_durations"
            assert row["check_summary"].startswith("5 checks asked")
            assert row["interval_standing"] == "rests on exchangeability"
        (failed_row,) = table[table["record"] == "q1_030423"].to_dict("records")
        assert failed_row["drawn"] and failed_row["reason"] == failure.note()
        (excluded_row,) = table[table["record"] == "q1_040423"].to_dict("records")
        assert not excluded_row["drawn"] and excluded_row["reason"] == EXCLUDED_REASON
    assert notes == {""}, "premise: the synthetic record's intervals all have width"
    # Hand tables for the notes. One death at tau = 10: SE 0, no width. Ages 1 (death) and
    # 100 (censored): S = 1/2 from 1, the mean to 100 is 1 + 99 / 2 = 50.5, its SE
    # sqrt(49.5^2 / (2 x 1)) = 35.0, so the 95% lower end is about -18.1.
    hand = {
        ("a", SIDES[0]): ([10.0], [True]),
        ("a", SIDES[1]): ([1.0, 100.0], [True, False]),
    }
    items = {
        key: curve_summaries(
            EventTable(dataset_id="hand", threshold_label=key[0], side=key[1],
                       age_s=np.array(ages), event=np.array(event)),
            conf_level=0.95,
        )
        for key, (ages, event) in hand.items()
    }  # fmt: skip
    z95 = NormalDist().inv_cdf(0.975)
    table = summaries_table(
        RecordSummaries(
            dataset_id="hand",
            flag="",
            failure="",
            ladder=("a",),
            conf_level=0.95,
            summaries=items,
            rmst_interval_s={
                key: (
                    item.rmst_s - z95 * item.rmst_se_s,
                    item.rmst_s + z95 * item.rmst_se_s,
                )
                for key, item in items.items()
            },
        ),
        excluded=(),
    )
    no_width, below = table.to_dict("records")
    assert (no_width["rmst_s"], no_width["rmst_se_s"]) == (10.0, 0.0)
    assert no_width["rmst_interval_note"] == NOTE_NO_WIDTH
    assert below["rmst_s"] == pytest.approx(50.5)
    assert below["rmst_lower_s"] == pytest.approx(50.5 - z95 * 49.5 / np.sqrt(2.0))
    assert below["rmst_interval_note"] == NOTE_LOWER_AT_ZERO


@pytest.mark.unit
def test_the_table_target_writes_the_frame_unrounded(scored, tmp_path: Path) -> None:
    """Oracle: the frame itself. The CSV read back with round-trip float parsing holds the
    frame's every value bit for bit; a plot with no table raises before writing."""
    frame = summaries_table(
        _summaries(scored, CONF_LEVEL, FLAG_TEXT), excluded=("q1_040423",)
    )
    render_table(SummariesTablePlot("summaries"), frame, tmp_path)
    read = pd.read_csv(tmp_path / "summaries.csv", float_precision="round_trip")
    assert list(read.columns) == list(frame.columns) and len(read) == len(frame)
    for column in frame.columns:
        if pd.api.types.is_float_dtype(frame[column]) or pd.api.types.is_integer_dtype(
            frame[column]
        ):
            np.testing.assert_array_equal(
                read[column].to_numpy(dtype=float, na_value=np.nan),
                frame[column].to_numpy(dtype=float, na_value=np.nan),
                err_msg=column,
            )
        else:
            assert read[column].fillna("").astype(str).tolist() == (
                frame[column].fillna("").astype(str).tolist()
            ), column
    with pytest.raises(TypeError, match="no build_table"):
        render_table(SurvivalBlockPlot("survival"), frame, tmp_path)


@pytest.mark.unit
def test_the_fraction_block_is_the_compliance_surveys_numbers_for_its_record(
    scored,
) -> None:
    """Oracles: `record_compliance` called on its own, the bound's definition (it holds
    the fraction), and the carve's gap spans summed. Each rung's fraction, bounds and
    statuses are the survey's, in ladder order, with lower <= fraction <= upper; the gaps
    the observed time leaves out are counted; a failed record has every rung not computed
    and no gap count; the compliance of another record is refused."""
    survey = record_compliance(scored.carve, T2STAR_THRESHOLDS)
    block = build_fraction_block(survey, scored.carve, flag="")
    assert block.ladder == tuple(label for label, _, _ in T2STAR_THRESHOLDS)
    for i, rung in enumerate(survey.rungs):
        assert (block.fraction[i], block.lower[i], block.upper[i]) == (
            rung.in_spec_fraction, rung.in_spec_fraction_lower, rung.in_spec_fraction_upper
        )  # fmt: skip
        assert (block.status[i], block.bounds_status[i]) == (
            rung.status,
            rung.bounds_status,
        )
        assert block.lower[i] <= block.fraction[i] <= block.upper[i]
    assert len({*block.fraction}) == len(block.ladder), "premise: the rungs differ"
    spans = scored.carve.diagnostics["gap_spans_s"]
    assert block.n_gaps == len(spans) >= 1
    assert block.gap_h == pytest.approx(sum(hi - lo for lo, hi in spans) / 3600.0)
    failure = RecordFailure(dataset_id=RECORD, stage="t2star", reason="ValueError: x")
    failed = build_fraction_block(
        record_compliance(failure, T2STAR_THRESHOLDS), failure, flag=""
    )
    assert failed.failure == failure.note() and failed.n_gaps is None
    assert set(failed.fraction) == {None}
    with pytest.raises(ValueError, match="is drawn with the carve of"):
        build_fraction_block(
            survey,
            RecordFailure(dataset_id="q9_010101", stage="t2star", reason="x"),
            flag="",
        )


@pytest.mark.unit
def test_the_fraction_block_draws_its_artifact() -> None:
    """Oracle: a hand artifact read back off the figure. A computed rung with its bounds
    draws its bar and whisker; a rung whose bounds are not computed its bar and a dotted
    whisker from 0 to 1; a rung not computed a grey cell and no bar; the gaps are counted
    in the axes and the flag reaches the caption."""
    data = FractionBlockData(
        dataset_id=RECORD,
        flag=FLAG_TEXT,
        failure="",
        ladder=("a", "b", "c"),
        fraction=(0.7, 0.4, None),
        lower=(0.6, None, None),
        upper=(0.8, None, None),
        status=(
            "computed",
            "computed",
            "not computed: no segment with at least 2 reads",
        ),
        bounds_status=("computed", "not computed: no read", "not computed: no read"),
        unsure_k=1.0,
        gap_mult=10.0,
        n_gaps=2,
        gap_h=1.5,
    )
    fig = FractionBlockPlot("fraction").build_matplotlib(data)
    try:
        (ax,) = fig.axes
        bars = [p for p in ax.patches if p.get_gid() == FRACTION_GID]
        assert [(b.get_x() + b.get_width() / 2, b.get_height()) for b in bars] == [
            (pytest.approx(0.0), pytest.approx(0.7)),
            (pytest.approx(1.0), pytest.approx(0.4)),
        ]
        (cell,) = [p for p in ax.patches if p.get_gid() == NOT_COMPUTED_GID]
        assert cell.get_x() + cell.get_width() / 2 == pytest.approx(2.0)
        assert cell.get_height() == 1.0
        whiskers = {
            gid: [
                seg.tolist()
                for c in ax.collections
                if c.get_gid() == gid
                for seg in c.get_segments()
            ]
            for gid in (BOUNDS_GID, NO_BOUNDS_GID)
        }
        assert whiskers[BOUNDS_GID] == [[[0.0, 0.6], [0.0, 0.8]]]
        assert whiskers[NO_BOUNDS_GID] == [[[1.0, 0.0], [1.0, 1.0]]]
        drawn = "\n".join(t.get_text() for t in fig.findobj(Text))
        assert "gaps left out: 2 (1.5 h)" in drawn
        assert FLAG_TEXT in drawn
    finally:
        plt.close(fig)


def _panel(scored: SimpleNamespace) -> object:
    return panel_data_of(
        scored.fit,
        scored.carve,
        kaplan_meier_of(scored.carve, conf_level=CONF_LEVEL),
        shape_min_reads=5,
        xi_seed=XI_SEED,
        k=1.0,
        use_uncertainty=True,
        thresholds=tuple(T2STAR_THRESHOLDS),
    )


def _panel_drawn_rungs(panel: object) -> tuple[list[str], str]:
    """The rungs today's within-calibration panel draws on its 2-state timeline, read off
    its figure, and that axes' title."""
    fig = WithinCalibrationPanel("panel").build_matplotlib(panel)
    try:
        (ax,) = [
            a for a in fig.axes if a.get_title().startswith("Compliance by threshold")
        ]
        return [t.get_text() for t in ax.get_yticklabels()][::-1], ax.get_title()
    finally:
        plt.close(fig)


@pytest.mark.unit
def test_the_timeline_block_draws_the_rungs_todays_panel_draws(scored) -> None:
    """Oracle: today's within-calibration panel, drawn on the same panel data. The block
    draws exactly the rungs the panel's timelines draw, in ladder order, on the record as
    carved (one rung culled) and on a copy whose occupancy is edited to 0.07, 0.03 and
    unmeasured at three rungs; each drawn rung's states and first read out of spec are the
    panel data's own."""
    panel = _panel(scored)
    block = build_timeline_block(panel, record=RECORD, flag="")
    drawn, _title = _panel_drawn_rungs(panel)
    assert list(block.drawn) == drawn
    assert block.culled_below, "premise: the synthetic record culls a rung"
    for label in block.drawn:
        assert block.states_4[label] == tuple(
            panel.distinguish.state_series_per_threshold[label]
        )
        assert block.states_2[label] == tuple(
            panel.reliability.compliance_state_series[label]
        )
        assert block.first_out_h[label] == panel.reliability.ttf_per_threshold.get(
            label
        )
    edited = copy.deepcopy(panel)
    first, second, third = (label for label, _, _ in T2STAR_THRESHOLDS[:3])
    edited.reliability.occupancy.update(
        {first: 0.07, second: 0.03, third: float("nan")}
    )
    block = build_timeline_block(edited, record=RECORD, flag="")
    drawn, title = _panel_drawn_rungs(edited)
    assert list(block.drawn) == drawn and first in drawn
    assert second in block.culled_below and block.culled_unmeasured == (third,)
    assert "1 not measured" in title


@pytest.mark.unit
def test_the_timeline_block_draws_its_artifact(scored) -> None:
    """Oracle: the artifact read back off the figure. Both timelines hold the drawn rungs
    top down in ladder order, each rung's bars its own states in their colours, the four
    states above and the two below; each rung's first read out of spec sits beside it; the
    culled rungs are named on the title when they fit and listed under the lower axes when
    they do not."""
    block = build_timeline_block(_panel(scored), record=RECORD, flag=FLAG_TEXT)
    fig = TimelineBlockPlot("timelines").build_matplotlib(block)
    try:
        upper, lower = [
            a for a in fig.axes if a.get_yticklabels() and a.get_xlim()[1] > 1
        ]
        for ax, states in ((upper, block.states_4), (lower, block.states_2)):
            assert [t.get_text() for t in ax.get_yticklabels()][::-1] == list(
                block.drawn
            )
            for label in block.drawn:
                bars = [p for p in ax.patches if p.get_gid() == f"state:{label}"]
                assert [(b.get_x(), b.get_x() + b.get_width()) for b in bars] == [
                    (pytest.approx(t0), pytest.approx(t1))
                    for t0, t1, _ in states[label]
                ]
                assert [b.get_facecolor() for b in bars] == [
                    to_rgba(state_color(state)) for _, _, state in states[label]
                ]
        drawn = "\n".join(t.get_text() for t in fig.findobj(Text))
        for label in block.drawn:
            hour = block.first_out_h[label]
            assert (f"{hour:.1f} h" if hour is not None else "never") in drawn
        assert (
            f"(culled, below 5% in spec: {', '.join(block.culled_below)})"
            in upper.get_title(loc="left")
        )
        assert FLAG_TEXT in drawn
    finally:
        plt.close(fig)
    many = TimelineBlockData(
        dataset_id=RECORD, flag="", failure="", ladder=tuple("abcdefghijklmnop"),
        drawn=("a",), culled_below=tuple(f"rung {c}" for c in "bcdefghijklm"),
        culled_unmeasured=("rung n", "rung o", "rung p"), states_4={"a": ((0.0, 1.0, "in_spec"),)},
        states_2={"a": ((0.0, 1.0, "in_spec"),)}, first_out_h={"a": None}, t_min_h=0.0,
        t_max_h=1.0,
    )  # fmt: skip
    fig = TimelineBlockPlot("timelines").build_matplotlib(many)
    try:
        drawn = "\n".join(t.get_text() for t in fig.findobj(Text))
        assert "(15 culled, listed below)" in drawn
        assert (
            "Culled, below 5% in spec: rung b" in drawn
            and "not measured: rung n" in drawn
        )
    finally:
        plt.close(fig)


@pytest.mark.policy
def test_each_figure_job_takes_one_ref_loads_only_data_and_draws_the_listed_records() -> (
    None
):
    """Oracle: the job graphs and the paper's tables. Each 6D2S figure job includes one
    ledger job and takes exactly one ref into it; every load reads under `data/`; it draws
    a survival block for each of its qubit's records the tables list, flags the shaded
    ones, and draws none of the four ledger records outside them. Together the six draw
    the tables' 24 6D2S rows."""
    assert len(FIGURE_JOBS) == 6
    drawn_all: list[str] = []
    for path in FIGURE_JOBS:
        job = _import_job(path)
        refs = {
            id(i)
            for n in job.dag.values()
            for i in n.inputs
            if isinstance(i, ArtifactRef)
        }
        assert len(refs) == 1, (path.name, len(refs))
        assert len(job.includes) == 1
        for node in job.dag.values():
            if "dataset" in node.kwargs:
                assert str(node.kwargs["dataset"].path).startswith("data/"), (
                    node.node_id
                )
        for node in job.dag.values():
            if node.fn is build_survival_block:
                record = node.node_id.removeprefix("survival_block_")
                assert node.kwargs["flag"] == (FLAG_TEXT if record in SHADED else "")
                drawn_all.append(record)
    paper_6d2s = [record for device, record in PAPER_RECORDS if device == "6D2S"]
    assert sorted(drawn_all) == sorted(paper_6d2s)
    excluded = [
        record
        for path in FIGURE_JOBS
        for record in _import_job(path).dag["curve_summaries_table"].kwargs["excluded"]
    ]
    assert sorted(excluded) == ["q4_070723", "q4_100723", "q5_070723", "q6_070723"]
    assert not {"q4_070723", "q4_100723", "q5_070723", "q6_070723"} & set(drawn_all)


@pytest.mark.unit
def test_each_block_is_handed_its_own_records_fit_windows_and_carve_kwargs(monkeypatch):
    """Oracle: the job graphs. Every `RecordNodes` a block receives names the fit and
    windows nodes its own record's carve wired (`t2star_<record>`, `windows_<record>`, the
    latter the input the cross-check holds), and carries the carve kwargs of that record's
    windows node in the included ledger job. Its tellability is the one node of the record,
    on that fit and those windows, with those kwargs and the standalone jobs' permutation
    count and seed (`tellability_recipe`)."""
    import quebra.record_figure_recipe as recipe

    seen = []
    monkeypatch.setattr(recipe, "BLOCKS", (lambda job, nodes: seen.append(nodes),))
    job = _import_job(FIGURE_JOBS[0])
    ledger_job = job.includes[0].job
    assert seen, "premise: the job wires at least one record"
    for nodes in seen:
        assert nodes.fit.node_id == f"t2star_{nodes.record}"
        assert nodes.windows.node_id == f"windows_{nodes.record}"
        checked = job.dag[nodes.carve.node_id]
        assert checked.inputs[0].node_id == nodes.windows.node_id
        (ledger_windows,) = [
            n
            for n in ledger_job.dag.values()
            if n.kwargs.get("stage") == "windows"
            and n.kwargs.get("record") == nodes.record
        ]
        assert dict(nodes.carve_kwargs) == {
            key: ledger_windows.kwargs[key] for key in recipe.CARVE_KWARGS
        }
        assert nodes.metric == "t2star"
        tellability = job.dag[nodes.tellability.node_id]
        assert tellability.node_id == f"tellability_t2star_{nodes.record}"
        assert [i.node_id for i in tellability.inputs] == [
            nodes.fit.node_id,
            nodes.windows.node_id,
        ]
        assert tellability.kwargs == {
            "thresholds": ledger_windows.kwargs["thresholds"],
            "gap_mult": ledger_windows.kwargs["gap_mult"],
            "unsure_k": ledger_windows.kwargs["k"],
            "n_permutations": TELLABILITY_PERMUTATIONS,
            "seed": TELLABILITY_SEED,
        }
    tellability_nodes = [
        n for n in job.dag.values() if n.fn is recipe.t2star_tellability_of
    ]
    assert len(tellability_nodes) == len(seen), "one tellability node per record"


@pytest.mark.unit
def test_a_listed_record_its_ledger_job_does_not_hold_raises() -> None:
    """Oracle: the ledger job's own records. A table row of this qubit the job does not
    hold is a wiring error."""
    with pytest.raises(ValueError, match="which it does not hold"):
        probe = Job("record_figure_probe")
        configure_record_figure_job(
            probe,
            included=probe.include("check_ledger_6d2s_q1"),
            device="6D2S",
            paper_records=(*PAPER_RECORDS, ("6D2S", "q1_310123")),
            shaded=SHADED,
            flag_text=FLAG_TEXT,
        )


@pytest.fixture(scope="module")
def told(tmp_path_factory) -> SimpleNamespace:
    """A synthetic record's T2* fit, real carve and real tellability, at 19 permutations."""
    root = tmp_path_factory.mktemp("resolvability")
    fit, carve = _carve(
        _track912_csv(root / "020423_r.csv", seed=2, gap=True, n=600), RECORD
    )
    kwargs = {
        "thresholds": T2STAR_THRESHOLDS,
        "gap_mult": carve.meta["gap_mult"],
        "unsure_k": carve.meta["k"],
        "n_permutations": 19,
        "seed": 1,
    }
    tellability = t2star_tellability_of(fit, carve, **kwargs)
    assert not isinstance(tellability, RecordFailure), tellability
    return SimpleNamespace(fit=fit, carve=carve, kwargs=kwargs, tellability=tellability)


def _g10_by_definition(f: np.ndarray, s: np.ndarray, lower: float) -> float:
    """G10 = L(f_top) / L(lower), L(f) the integral from 0 to f of (2 pi f')^2 S(f') df'
    (docs/TELLABILITY_SEMANTICS.md, Q2.3): S held at its first value below the first grid
    frequency, trapezoids above it, and L at `lower` interpolated between grid points.
    Written out here, independent of `spectrum.derivative_variance`."""
    integrand = [(2.0 * np.pi * fi) ** 2 * si for fi, si in zip(f, s, strict=True)]
    L = [(2.0 * np.pi) ** 2 * s[0] * f[0] ** 3 / 3.0]
    for i in range(1, len(f)):
        L.append(L[-1] + 0.5 * (integrand[i] + integrand[i - 1]) * (f[i] - f[i - 1]))
    j = int(np.searchsorted(f, lower))
    if f[j] == lower:
        at_lower = L[j]
    else:
        w = (lower - f[j - 1]) / (f[j] - f[j - 1])
        at_lower = L[j - 1] + w * (L[j] - L[j - 1])
    return L[-1] / at_lower


@pytest.mark.unit
def test_g10_recomputed_from_the_stored_spectrum_and_top_decade_is_the_stored_g10(
    told,
) -> None:
    """Oracle: G10's definition, written out in the test. The block stores the spectrum and
    the decade G10 reads; G10 recomputed on them equals the G10 in the record's table, and
    the decade is the one below the last grid frequency."""
    block = build_resolvability_block(told.tellability, flag="")
    lower, upper = block.top_decade_hz
    assert upper == block.frequency_hz[-1]
    recomputed = _g10_by_definition(block.frequency_hz, block.spectrum_per_hz, lower)
    stored = float(block.per_record.loc[0, "g10"])
    assert recomputed == pytest.approx(stored, rel=1e-9)
    assert block.per_record.loc[0, "top_decade_lower_hz"] == lower


@pytest.mark.unit
def test_the_block_and_its_tables_carry_the_tellability_it_reads(told) -> None:
    """Oracle: the tellability's own fields, and the definitions in
    docs/TELLABILITY_SEMANTICS.md. The drawn arrays and marks, every rung's count ratio and
    crossing counts, and both tables equal the fields they come from; the edge level is the
    mean of S over the stored stretch, which ends at the band edge; the count ratio's limits
    are 1 (independent reads) and m (a smooth signal), Q2.2; G10 and its halves carry their
    statuses and DEC-12a's one sentence."""
    tellability = told.tellability
    metric = tellability.metric
    block = build_resolvability_block(tellability, flag="shaded")
    np.testing.assert_array_equal(block.frequency_hz, metric.spectrum_frequency_hz)
    np.testing.assert_array_equal(block.spectrum_per_hz, metric.spectrum_per_hz)
    assert block.band_edge_hz == metric.band_edge_hz
    assert block.edge_level_per_hz == metric.edge_level_per_hz
    assert block.noise_floor_per_hz == metric.noise_floor_per_hz
    lower, upper = block.edge_stretch_hz
    assert upper == metric.band_edge_hz
    on = (block.frequency_hz >= lower) & (block.frequency_hz <= upper)
    assert on.sum() >= 2, "premise: the stretch holds grid frequencies"
    # abs=0: S is of order 1e-12 here, the size of approx's default absolute tolerance
    assert block.edge_level_per_hz == pytest.approx(
        float(np.mean(block.spectrum_per_hz[on])), rel=1e-12, abs=0.0
    )
    m = block.decimation_ladder
    assert m == (1, 2, 3, 5, 8, 13, 20)
    assert block.count_limits == {
        LIMIT_INDEPENDENT: (1.0,) * len(m),
        LIMIT_SMOOTH: tuple(float(x) for x in m),
    }
    assert (block.value_key, block.threshold_value_key) == ("t2star_s", "t2star_s")
    record = block.per_record.iloc[0]
    assert record["g10"] == metric.g10
    assert record["g10_first_half"] == metric.g10_first_half
    assert record["g10_second_half"] == metric.g10_second_half
    assert record["g10_status"] == metric.status["g10"]
    assert record["g10_reading"] == "read one way: small means a red spectrum"
    assert record["window_band_max"] == metric.window_band_max
    assert block.ladder == tuple(row.threshold_label for row in tellability.thresholds)
    table = block.per_threshold.set_index("threshold_label")
    for row in tellability.thresholds:
        label = row.threshold_label
        assert block.count_ratio[label] == row.count_ratio
        assert block.crossing_counts[label] == row.crossing_counts
        for m, ratio in zip(row.decimation_ladder, row.count_ratio, strict=True):
            assert table.loc[label, f"count_ratio_m{m}"] == ratio
        for count in row.crossing_counts:
            for kind in ("clean", "merged", "unbracketed"):
                assert table.loc[label, f"{kind}_k{count.k:g}"] == getattr(count, kind)
        assert table.loc[label, "count_exponent"] == row.count_exponent
        assert table.loc[label, "clean_fraction_k1"] == row.clean_fraction_k1


def _not_computed(
    row: Any, groups: dict[str, tuple[str, ...]], reason: str, *names: str
) -> Any:
    """`row` with each named group not computed for `reason`: its fields None and its status
    the reason, as the tellability leaves a group it could not compute."""
    fields = {field: None for name in names for field in groups[name]}
    status = {**row.status, **dict.fromkeys(names, reason)}
    return dataclasses.replace(row, **fields, status=status)


@pytest.mark.unit
def test_a_group_the_tellability_did_not_compute_shows_its_status(told) -> None:
    """Oracle: the tellability's statuses. With the spectrum not computed (and so the noise
    floor and G10, as the tellability leaves them), one rung's count ratio and another's
    crossing counts not computed: the block holds no spectrum, decade or stretch, those
    rungs have no ratio or no counts and blank table cells, every other rung keeps its own,
    and the figure prints the reasons. With the spectrum computed but the band edge, noise
    floor and G10 not: no mark is drawn, no decade is stored, and the caption names no mark
    and gives each one's status."""
    tellability = told.tellability
    reason = "not computed: zero variance"
    metric = _not_computed(
        tellability.metric, METRIC_GROUPS, reason, "spectrum", "noise_floor", "g10"
    )
    thresholds = list(tellability.thresholds)
    placed, counted = thresholds[2], thresholds[3]
    thresholds[2] = _not_computed(
        placed,
        THRESHOLD_GROUPS,
        "not computed: no down-crossing at this threshold",
        "crossing_counts",
    )
    uncounted = "not computed: fewer than two rungs with a down-crossing"
    thresholds[3] = _not_computed(
        counted, THRESHOLD_GROUPS, uncounted, "count_exponent"
    )
    block = build_resolvability_block(
        dataclasses.replace(tellability, metric=metric, thresholds=tuple(thresholds)),
        flag="",
    )
    assert block.spectrum_status == reason and block.g10_status == reason
    assert block.frequency_hz.size == 0 and block.spectrum_per_hz.size == 0
    assert block.top_decade_hz is None and block.edge_stretch_hz is None
    table = block.per_threshold.set_index("threshold_label")
    assert block.crossing_counts[placed.threshold_label] == ()
    assert (
        table.loc[placed.threshold_label, ["clean_k1", "clean_k2", "clean_fraction_k1"]]
        .isna()
        .all()
    )
    assert block.count_ratio[counted.threshold_label] == ()
    assert block.count_status[counted.threshold_label] == uncounted
    ratio_cells = [f"count_ratio_m{x}" for x in block.decimation_ladder]
    assert (
        table.loc[counted.threshold_label, [*ratio_cells, "rungs_excluded"]]
        .isna()
        .all()
    )
    for other in tellability.thresholds:
        label = other.threshold_label
        if label != placed.threshold_label:
            assert block.crossing_counts[label] == other.crossing_counts
        if label != counted.threshold_label:
            assert block.count_ratio[label] == other.count_ratio
    fig = ResolvabilityBlockPlot("status").build_matplotlib(block)
    texts = [t.get_text() for t in fig.findobj(Text)]
    assert any("spectrum not computed: zero variance" in t for t in texts)
    assert any(t == "no down-crossing at this threshold" for t in texts)
    drawn = {line.get_label() for line in fig.axes[1].get_lines()}
    assert counted.threshold_label not in drawn and placed.threshold_label in drawn
    plt.close(fig)

    marks = _not_computed(
        tellability.metric,
        METRIC_GROUPS,
        "not computed: span shorter than 20 median spacings",
        "window",
        "g10",
    )
    marks = _not_computed(marks, METRIC_GROUPS, "not declared", "noise_floor")
    block = build_resolvability_block(
        dataclasses.replace(tellability, metric=marks), flag=""
    )
    assert block.frequency_hz.size and block.top_decade_hz is None
    assert (
        block.per_record.loc[0, ["top_decade_lower_hz", "top_decade_upper_hz"]]
        .isna()
        .all()
    )
    fig = ResolvabilityBlockPlot("marks").build_matplotlib(block)
    spectrum_ax = fig.axes[0]
    assert {line.get_label() for line in spectrum_ax.get_lines()} >= {"S of t2star_s"}
    assert not {"band edge", "edge level", "noise floor"} & {
        line.get_label() for line in spectrum_ax.get_lines()
    }
    assert not [p for p in spectrum_ax.patches if p.get_label() == "top decade (G10)"]
    caption = fig._supxlabel.get_text()
    assert "the spectrum of t2star_s; per rung" in caption, caption
    short = "not computed: span shorter than 20 median spacings"
    for note in (
        f"band edge ({short})",
        "edge level (no band edge to place it on)",
        "noise floor (not declared)",
        f"the top decade G10 reads ({short})",
    ):
        assert note in caption, (note, caption)
    plt.close(fig)


@pytest.mark.unit
def test_a_record_whose_carve_or_tellability_failed_reads_not_computed(told) -> None:
    """Oracle: decision 4 of #21, builders raise on wiring only. A failed carve passes
    through, the windows' failure first; a tellability that refuses the record's content
    (a negative sigma) becomes a failure of the stage `tellability`; the block of either
    draws the note and no panel."""
    failed = RecordFailure(dataset_id=RECORD, stage="windows", reason="ValueError: x")
    assert t2star_tellability_of(told.fit, failed, **told.kwargs) is failed
    frame = told.fit.frame.copy()
    frame.loc[frame.index[5], "t2star_error_s"] = -1.0e-7
    refused = t2star_tellability_of(
        dataclasses.replace(told.fit, frame=frame), told.carve, **told.kwargs
    )
    assert isinstance(refused, RecordFailure) and refused.stage == "tellability"
    assert "a reported sigma is negative" in refused.reason
    block = build_resolvability_block(refused, flag="")
    assert block.failure == refused.note() and block.ladder == ()
    assert block.per_record.loc[0, "failure"] == refused.note()
    fig = ResolvabilityBlockPlot("failed").build_matplotlib(block)
    texts = [t.get_text() for t in fig.findobj(Text)]
    assert f"not computed: {refused.note()}" in texts
    plt.close(fig)


@pytest.mark.unit
def test_the_resolvability_figure_draws_the_artifacts_values(told) -> None:
    """Oracle: the artifact's own fields, read back off the drawn figure: the spectrum, its
    band edge, edge level and noise floor, each named in the caption; the two limits, each
    in its own dash, and every rung's count ratio at its decimation; and every crossing
    count as one segment of its rung's stack, at k = 1 left of k = 2 as the caption says."""
    block = build_resolvability_block(told.tellability, flag="")
    fig = ResolvabilityBlockPlot("values").build_matplotlib(block)
    spectrum_ax, counts_ax, crossings_ax = fig.axes[:3]
    lines = {line.get_label(): line for line in spectrum_ax.get_lines()}
    np.testing.assert_array_equal(
        lines["S of t2star_s"].get_ydata(), block.spectrum_per_hz
    )
    assert lines["band edge"].get_xdata()[0] == block.band_edge_hz
    assert tuple(lines["edge level"].get_xdata()) == block.edge_stretch_hz
    assert lines["edge level"].get_ydata()[0] == block.edge_level_per_hz
    assert lines["noise floor"].get_ydata()[0] == block.noise_floor_per_hz
    caption = fig._supxlabel.get_text()
    assert (
        "with its band edge, edge level, noise floor and the top decade G10 reads"
        in caption
    )
    lines = {line.get_label(): line for line in counts_ax.get_lines()}
    m = list(block.decimation_ladder)
    for name in (LIMIT_INDEPENDENT, LIMIT_SMOOTH):
        guide = lines[f"limit: {name}"]
        assert list(guide.get_xdata()) == m
        assert tuple(guide.get_ydata()) == block.count_limits[name]
        assert guide.get_linestyle() == LIMIT_DASHES[name]
    assert (
        lines[f"limit: {LIMIT_INDEPENDENT}"].get_linestyle()
        != lines[f"limit: {LIMIT_SMOOTH}"].get_linestyle()
    )
    for label in block.ladder:
        ratio = block.count_ratio[label]
        kept = [
            (x, r) for x, r in zip(m, ratio, strict=True) if r is not None and r > 0
        ]
        assert list(lines[label].get_xdata()) == [x for x, _ in kept]
        assert list(lines[label].get_ydata()) == [r for _, r in kept]
    # One slot per (rung, k), left to right in ladder order and increasing k, each a stack
    # clean, merged, unbracketed from 0.
    slots: dict[float, list[tuple[float, float]]] = {}
    for bar in crossings_ax.patches:
        if bar.get_gid() == CROSSING_GID:
            centre = round(bar.get_x() + bar.get_width() / 2.0, 9)
            slots.setdefault(centre, []).append((bar.get_y(), bar.get_height()))
    expected = [
        (i, count)
        for i, label in enumerate(block.ladder)
        for count in sorted(block.crossing_counts[label], key=lambda c: c.k)
    ]
    assert len(slots) == len(expected)
    for centre, (i, count) in zip(sorted(slots), expected, strict=True):
        assert i - 0.5 < centre < i + 0.5
        clean, merged = count.clean, count.merged
        assert sorted(slots[centre]) == sorted(
            [(0, clean), (clean, merged), (clean + merged, count.unbracketed)]
        )
    plt.close(fig)


@pytest.mark.unit
def test_a_count_ratio_of_zero_sits_on_the_axis_floor_and_the_caption_names_it(
    told,
) -> None:
    """Oracle: the artifact. A count ratio of exactly 0, which a log axis cannot place, is
    drawn as one marker at its decimation on the axis floor in its rung's colour, left off
    its rung's line, and named in the caption with its rung and m."""
    tellability = told.tellability
    thresholds = list(tellability.thresholds)
    rung = thresholds[1]
    assert rung.status["count_exponent"] == "computed", "premise: the rung is counted"
    ratio = list(rung.count_ratio)
    ratio[3] = 0.0
    thresholds[1] = dataclasses.replace(rung, count_ratio=tuple(ratio))
    block = build_resolvability_block(
        dataclasses.replace(tellability, thresholds=tuple(thresholds)), flag=""
    )
    zero_m = block.decimation_ladder[3]
    fig = ResolvabilityBlockPlot("zero").build_matplotlib(block)
    counts_ax = fig.axes[1]
    (marker,) = [
        line for line in counts_ax.get_lines() if line.get_gid() == ZERO_RATIO_GID
    ]
    assert list(marker.get_xdata()) == [zero_m]
    assert list(marker.get_ydata()) == [counts_ax.get_ylim()[0]]
    lines = {line.get_label(): line for line in counts_ax.get_lines()}
    assert marker.get_color() == lines[rung.threshold_label].get_color()
    assert zero_m not in list(lines[rung.threshold_label].get_xdata())
    caption = fig._supxlabel.get_text()
    assert f"{rung.threshold_label} at m = {zero_m}" in caption, caption
    plt.close(fig)
    fig = ResolvabilityBlockPlot("none").build_matplotlib(
        build_resolvability_block(tellability, flag="")
    )
    assert not [ln for ln in fig.axes[1].get_lines() if ln.get_gid() == ZERO_RATIO_GID]
    assert "count ratio of 0" not in fig._supxlabel.get_text()
    plt.close(fig)


@pytest.mark.unit
def test_a_wiring_error_raises_when_the_job_is_built() -> None:
    """Oracle: the wiring rule. A metric outside the two a record figure draws, and carve
    kwargs or a permutation count every record would refuse, raise at build time instead of
    reading `not computed` record by record, and the recipe that builds a job checks them."""
    with pytest.raises(ValueError, match="is not one of"):
        RecordNodes(
            record="r",
            carve=None,
            ledger=None,
            run_set=(),
            flag="",
            fit=None,
            windows=None,
            carve_kwargs={},
            metric="fidelity",
            tellability=None,
        )
    good = {"k": 1.0, "gap_mult": 10.0, "thresholds": T2STAR_THRESHOLDS}
    check_tellability_kwargs("r", good, n_permutations=19, seed=0)
    for bad, match in (
        ({**good, "k": 0.0}, "k must be finite and positive"),
        ({**good, "gap_mult": float("nan")}, "gap_mult must be finite"),
        (
            {**good, "thresholds": [*T2STAR_THRESHOLDS[:2], T2STAR_THRESHOLDS[0]]},
            "repeated",
        ),
    ):
        with pytest.raises(ValueError, match=match):
            check_tellability_kwargs("r", bad, n_permutations=19, seed=0)
    with pytest.raises(ValueError, match="n_permutations >= 1"):
        check_tellability_kwargs("r", good, n_permutations=0, seed=0)
    probe = Job("record_figure_probe")
    with pytest.raises(ValueError, match="n_permutations >= 1"):
        configure_record_figure_job(
            probe,
            included=probe.include("check_ledger_6d2s_q1"),
            device="6D2S",
            paper_records=PAPER_RECORDS,
            shaded=SHADED,
            flag_text=FLAG_TEXT,
            tellability_permutations=0,
        )


_LEDGER_JOB = """
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.ledger_recipe import LEDGER_KNOBS, configure_check_ledger_job
from quebra.schemas.track912 import track912Schema

job = Job("ledger_probe")
configure_check_ledger_job(
    job,
    prefix="probe",
    records=[
        Dataset(path="{a}", schema=track912Schema, qubit=1, device="probe",
                extra={{"run_name": "q1_010423", "run_start_unix_s": 1.68e9}}),
        Dataset(path="{b}", schema=track912Schema, qubit=1, device="probe",
                extra={{"run_name": "q1_020423", "run_start_unix_s": 1.68e9}}),
    ],
    bench_size_table=Dataset(path="{bench}", schema=None),
    **{{**LEDGER_KNOBS, "n_permutations": 19}},
)
"""


@pytest.mark.integration
@pytest.mark.slow
def test_a_figure_job_runs_its_ledger_job_once_and_draws_each_record(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Oracle: the run directory and the blocks it persisted. A figure job over an invented
    two-record ledger job runs that job once, through its one ref, renders every block per
    record, and writes the curve summaries and the two resolvability tables, each table
    the drawn records' blocks stacked."""
    a = _track912_csv(tmp_path / "010423_r.csv", seed=1, gap=False)
    b = _track912_csv(tmp_path / "020423_r.csv", seed=3, gap=True)
    ledger_py = tmp_path / "ledger_probe.py"
    ledger_py.write_text(_LEDGER_JOB.format(a=a.name, b=b.name, bench=BENCH))
    (tmp_path / BENCH.name).write_bytes(BENCH.read_bytes())
    job = Job("figures_probe")
    job.job_file = (tmp_path / "figures_probe.py").resolve()
    job.job_file.write_text("# synthetic figure job\n")
    configure_record_figure_job(
        job,
        included=job.include(str(ledger_py)),
        device="probe",
        paper_records=(("probe", "q1_010423"), ("probe", "q1_020423")),
        shaded=("q1_020423",),
        flag_text=FLAG_TEXT,
        tellability_permutations=19,
    )
    monkeypatch.chdir(tmp_path)
    out = tmp_path / "out"
    run_job(job, out, force=True, data_root=tmp_path)
    (run_dir,) = out.glob("figures_probe_*")
    records = [
        json.loads(p.read_text()) for p in (run_dir / "provenance").glob("*.prov.json")
    ]
    nested = {
        json.dumps(include, sort_keys=True)
        for r in records
        for include in r["includes"]
    }
    assert len(nested) == 1, nested
    for record in ("q1_010423", "q1_020423"):
        for block in ("survival", "resolvability"):
            assert list(run_dir.glob(f"{record}_{block}_*.pdf")), sorted(
                p.name for p in run_dir.iterdir()
            )
    summaries = pd.read_csv(run_dir / "figures_probe_curve_summaries.csv")
    assert len(summaries) == 2 * len(T2STAR_THRESHOLDS) * len(SIDES)
    assert set(summaries["record"]) == {"q1_010423", "q1_020423"}
    # The two resolvability CSVs are every drawn record's tables, stacked, as the blocks
    # hold them (the persisted figure inputs); read as text, so no dtype is guessed.
    blocks = []
    for record in ("q1_010423", "q1_020423"):
        with (run_dir / f"resolvability_block_{record}.pkl").open("rb") as handle:
            blocks.append(pickle.load(handle))
    for part in ("per_threshold", "per_record"):
        written = pd.read_csv(
            run_dir / f"figures_probe_resolvability_{part}.csv",
            dtype=str,
            keep_default_na=False,
        )
        expected = pd.concat(
            [
                pd.read_csv(
                    io.StringIO(getattr(b, part).to_csv(index=False)),
                    dtype=str,
                    keep_default_na=False,
                )
                for b in blocks
            ],
            ignore_index=True,
        )
        assert len(expected) == (
            2 * len(T2STAR_THRESHOLDS) if part == "per_threshold" else 2
        )
        pd.testing.assert_frame_equal(written, expected)


def _five_blocks(scored: SimpleNamespace, told: SimpleNamespace) -> dict[str, object]:
    """The five blocks of the synthetic record, by name, in the page's reading order."""
    km = kaplan_meier_of(scored.carve, conf_level=CONF_LEVEL)
    return {
        "resolvability": build_resolvability_block(told.tellability, flag=""),
        "timeline": build_timeline_block(_panel(scored), record=RECORD, flag=""),
        "fraction": build_fraction_block(
            record_compliance(scored.carve, T2STAR_THRESHOLDS), scored.carve, flag=""
        ),
        "survival": _block(km, scored.ledger),
        "hazard": _hazard_block(scored.carve, scored.ledger),
    }


# Each block's figure size in inches on the synthetic record, measured on the block plots
# before their sizes moved into `<block>_size` (commit 65883e9), and by hand: survival
# 1.9 x 7 rungs + 1.7 wide, (2.4 + 0.3 + 0.19 x 5 checks) x 2 bags + 1.2 high; hazard
# max(4.5, 0.34 x 2 clocks x 7) x 2 + 1.7 wide, 3.0 + 0.17 x 4 + 0.55 + 0.19 x 5 + 0.9 high;
# fraction max(0.75 x 7 + 1.5, 6) wide, 3.2 + 1.2 high; resolvability 3 x 3.4 wide, 2.8 + 1.2
# high; timeline 12 wide, 0.45 + 2 x (0.28 x 6 drawn + 0.8) + 1.2 high.
_SIZES_BEFORE = {
    "survival": (15.0, 8.5),
    "hazard": (11.22, 6.08),
    "fraction": (6.75, 4.4),
    "resolvability": (10.2, 4.0),
    "timeline": (12.0, 6.61),
}
_SIZE_OF = {
    "survival": (record_figure_plot.survival_size, SurvivalBlockPlot),
    "hazard": (record_figure_plot.hazard_size, HazardBlockPlot),
    "fraction": (record_figure_plot.fraction_size, FractionBlockPlot),
    "resolvability": (record_figure_plot.resolvability_size, ResolvabilityBlockPlot),
    "timeline": (record_figure_plot.timeline_size, TimelineBlockPlot),
}


@pytest.mark.unit
def test_each_blocks_size_is_the_size_its_figure_had_before_the_page(
    scored, told
) -> None:
    """Oracle: the sizes measured on the block plots before the refactor, which the hand
    formulas above give. Each block's size function, and the figure its own plot now
    draws, are those sizes; the resolvability block's plot, which keeps its own formula,
    agrees with `resolvability_size`, the page's reading of it."""
    for name, block in _five_blocks(scored, told).items():
        size_of, plot = _SIZE_OF[name]
        assert size_of(block) == pytest.approx(_SIZES_BEFORE[name], rel=1e-12), name
        fig = plot(name).build_matplotlib(block)
        try:
            assert tuple(fig.get_size_inches()) == pytest.approx(
                _SIZES_BEFORE[name], rel=1e-12
            ), name
        finally:
            plt.close(fig)


@pytest.mark.unit
def test_the_page_stacks_each_block_once_in_reading_order(scored, told) -> None:
    """Oracle: each block's own figure. The page has one constrained row per block in the
    reading order it was handed, each row holding as many axes as that block's own figure,
    its width the widest block's and its height every block's plus its title; drawing it
    raises no layout warning. Every figure job hands its page each record's five blocks in
    `PAGE_ORDER`."""
    blocks = _five_blocks(scored, told)
    page = assemble_record_figure(*blocks.values(), record=RECORD)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fig = RecordFigurePlot("page").build_matplotlib(page)
        try:
            fig.canvas.draw()
            assert type(fig.get_layout_engine()).__name__ == "ConstrainedLayoutEngine"
            assert len(fig.subfigs) == len(blocks)
            for row, (name, block) in zip(fig.subfigs, blocks.items(), strict=True):
                own = _SIZE_OF[name][1](name).build_matplotlib(block)
                try:
                    assert len(row.axes) == len(own.axes), name
                finally:
                    plt.close(own)
            width = max(w for w, _ in _SIZES_BEFORE.values())
            height = sum(h for _, h in _SIZES_BEFORE.values()) + 0.5
            assert tuple(fig.get_size_inches()) == pytest.approx((width, height))
        finally:
            plt.close(fig)
    assert not [w for w in caught if "layout" in str(w.message).lower()], caught
    with pytest.raises(ValueError, match="handed blocks of"):
        assemble_record_figure(*blocks.values(), record="q9_010101")
    assert set(BLOCKS) == set(PAGE_ORDER) and len(PAGE_ORDER) == len(set(PAGE_ORDER))
    job = _import_job(FIGURE_JOBS[3])
    pages = [n for n in job.dag.values() if n.fn is assemble_record_figure]
    assert pages
    for node in pages:
        assert [job.dag[i.node_id].fn_name.rsplit("_", 2)[0] for i in node.inputs] == [
            "resolvability_block", "timeline_block", "fraction_block", "survival_block",
            "hazard_block",
        ]  # fmt: skip
