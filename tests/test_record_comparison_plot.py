"""The comparison figures: what each image states, and the data its builders hand it.

Synthetic carves only: noisy sines, and a record that rises into spec once and stays, so
it has no out-of-spec window with a seen birth and no in-spec death before tau.
"""

from __future__ import annotations

import zlib

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from quebra import comparison_recipe  # noqa: E402
from quebra.analyzers import event_table, kaplan_meier, record_comparison, windows  # noqa: E402
from quebra.analyzers.survival import ByThresholdSide  # noqa: E402
from quebra.analyzers.windows import SIDE_IN_SPEC, SIDE_OUT_OF_SPEC  # noqa: E402
from quebra.core.job import Job  # noqa: E402
from quebra.plots import record_comparison_plot as plot  # noqa: E402

pytestmark = pytest.mark.unit

LADDER = [("4", 4.0, True), ("6", 6.0, True)]
RULE = record_comparison.TAU_RULE_MIN_LARGEST_AGE
PER_GROUP = record_comparison.TAU_PER_GROUP


def _carve(record: str, phase: float) -> windows.WindowsResult:
    t = np.arange(120.0)
    rng = np.random.default_rng(zlib.crc32(record.encode()))
    v = 5.0 + 3.0 * np.sin(t / 3.0 + phase) + rng.normal(0.0, 0.8, t.size)
    return windows.run(
        windows.WindowsInputs(
            t_rel_s=t, values=v, thresholds=LADDER, dataset_id=record, gap_mult=10.0
        )
    )


def _rise(record: str) -> windows.WindowsResult:
    t = np.arange(120.0)
    return windows.run(
        windows.WindowsInputs(
            t_rel_s=t,
            values=np.where(t < 20.0, 2.0, 9.0),
            thresholds=LADDER,
            dataset_id=record,
            gap_mult=10.0,
        )
    )


def _texts(fig: plt.Figure) -> list[str]:
    return [t.get_text() for t in fig.texts] + [
        t.get_text() for ax in fig.axes for t in ax.texts
    ]


def _draw(plot_class: type, data: object) -> plt.Figure:
    return plot_class(name="unit").build_matplotlib(data)


def _flat(text: str) -> str:
    return " ".join(text.split())


@pytest.fixture(scope="module")
def pair_data() -> record_comparison.PairLadderData:
    a, rise = _carve("a", 0.0), _rise("rise")
    pairs = record_comparison.compare_pair_set(a, rise, tau_rule=RULE)
    return record_comparison.build_pair_ladder_data(
        pairs, a, rise, comparison="p", evidence="one synthetic run"
    )


@pytest.fixture(scope="module")
def forest_set() -> record_comparison.RandomEffectsSet:
    carves = [_carve(f"r{i}", 0.7 * i) for i in range(1, 6)] + [
        _rise("rise"),
        _rise("rise2"),
    ]
    return record_comparison.random_effects_set(
        *carves,
        groups={
            "big": ["r1", "r2", "r3", "rise"],
            "small": ["r4", "r5"],
            "lone": ["rise2"],
        },
        tau_scope=PER_GROUP,
        tau_rule=RULE,
    )


EVIDENCE = {"big": "line big", "small": "line small", "lone": "line lone"}


def test_a_pair_figure_has_a_panel_per_rung_and_side_and_says_where_nothing_was_computed(
    pair_data: record_comparison.PairLadderData,
) -> None:
    """Oracle: specification; each panel's numbers read off the pair it draws.

    The rise has no out-of-spec window, so both out-of-spec panels must say "not computed"
    with the reason, and every in-spec panel must print tau, the record that set it, the
    difference and the log-rank p.
    """
    fig = _draw(plot.PairLadderPlot, pair_data)
    try:
        panels = {ax.get_title(): ax for ax in fig.axes if ax.get_title()}
        assert len(panels) == 2 * len(LADDER)
        for (label, side), panel in pair_data.panels.items():
            words = "in spec" if side == SIDE_IN_SPEC else "out of spec"
            text = _flat(
                " ".join(t.get_text() for t in panels[f"{words}, {label}"].texts)
            )
            pair = panel.pair
            tau = f"τ = {pair.tau.rmst_tau_s / 60.0:.3g} min, set by {pair.tau.tau_set_by}"
            assert tau in text, text
            if side == SIDE_OUT_OF_SPEC:
                assert pair.not_computed_reason is not None
                assert _flat(f"not computed: {pair.not_computed_reason}") in text
            else:
                assert pair.not_computed_reason is None and pair.logrank_p is not None
                assert f"log-rank p = {pair.logrank_p:.2g}" in text, text
                assert f"a - b = {pair.rmst_difference_s / 60.0:.3g} min" in text
            assert (
                f"windows (deaths): a {pair.n_windows_a} ({pair.n_deaths_a}), "
                f"b {pair.n_windows_b} ({pair.n_deaths_b})"
            ) in text
            assert (
                f"birth unseen, excluded: a {panel.curve_a.n_unobserved_birth_dropped}, "
                f"b {panel.curve_b.n_unobserved_birth_dropped}"
            ) in text
    finally:
        plt.close(fig)


def test_a_pair_caption_says_each_p_is_unadjusted_and_carries_evidence_and_checks(
    pair_data: record_comparison.PairLadderData,
) -> None:
    """Oracle: specification. The p of one panel is one of many; the caption says so."""
    fig = _draw(plot.PairLadderPlot, pair_data)
    try:
        caption = _flat(" ".join(t.get_text() for t in fig.texts))
        legend = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    finally:
        plt.close(fig)
    assert legend == ["a: a", "b: rise", "comparison tau"]
    assert "unadjusted for the other thresholds and sides" in caption
    assert "Evidence: one synthetic run." in caption
    first = next(iter(pair_data.panels.values())).pair
    assert _flat(first.check_summary()) in caption


def test_the_pair_figure_curves_are_kaplan_meier_on_the_tables_the_pair_compared(
    pair_data: record_comparison.PairLadderData,
) -> None:
    """Oracle: `kaplan_meier.curve_from_event_table`, called here on each carve's tables.

    Carves handed in the other order than the pair compared them are refused.
    """
    a, rise = _carve("a", 0.0), _rise("rise")
    tables = {
        "a": event_table.event_tables_from_carve(a),
        "rise": event_table.event_tables_from_carve(rise),
    }
    for key, panel in pair_data.panels.items():
        for curve, record in ((panel.curve_a, "a"), (panel.curve_b, "rise")):
            expected = kaplan_meier.curve_from_event_table(
                tables[record].table(*key), conf_level=panel.pair.conf_level
            )
            np.testing.assert_array_equal(curve.band_lower, expected.band_lower)
            np.testing.assert_array_equal(curve.band_upper, expected.band_upper)
            assert curve.dataset_id == record
            np.testing.assert_array_equal(curve.time_min, expected.time_min)
            np.testing.assert_array_equal(curve.survival, expected.survival)
            assert curve.n_windows == expected.n_windows
    pairs = record_comparison.compare_pair_set(a, rise, tau_rule=RULE)
    with pytest.raises(ValueError, match="given for the pair"):
        record_comparison.build_pair_ladder_data(
            pairs, rise, a, comparison="p", evidence="x"
        )
    first_rung = ByThresholdSide(
        ladder=("4",), items={k: v for k, v in pairs.items.items() if k[0] == "4"}
    )
    with pytest.raises(ValueError, match="carves on the ladder"):
        record_comparison.build_pair_ladder_data(
            first_rung, a, rise, comparison="p", evidence="x"
        )
    with pytest.raises(ValueError, match="no threshold and side to draw"):
        record_comparison.build_pair_ladder_data(
            ByThresholdSide(ladder=pairs.ladder, items={}),
            a,
            rise,
            comparison="p",
            evidence="x",
        )


def test_a_group_too_small_to_pool_is_drawn_as_a_table_not_a_forest(
    forest_set: record_comparison.RandomEffectsSet,
) -> None:
    """Oracle: specification (the pooling floor). An unpooled group has no pooled numbers,
    so its panel is text: no line, no marker, each record named with its restricted mean."""
    data = record_comparison.build_forest_data(
        forest_set,
        comparison="f",
        evidence=EVIDENCE,
        threshold_label="4",
        side=SIDE_IN_SPEC,
    )
    fig = _draw(plot.ForestPlot, data)
    try:
        by_group = {
            ax.get_title(loc="left"): ax for ax in fig.axes if ax.get_title(loc="left")
        }
        small = by_group["small"]
        assert data.summaries["small"].pooled_s is None
        assert not small.lines and not small.collections and not small.patches
        text = _flat(" ".join(t.get_text() for t in small.texts))
        for row in data.summaries["small"].rows:
            assert f"{row.record}: {row.summary.rmst_s / 60.0:.3g} min" in text
        assert by_group["big"].lines, "the pooled group drew no forest"
    finally:
        plt.close(fig)


def test_the_forest_names_every_record_left_out_and_prints_tau_and_its_record(
    forest_set: record_comparison.RandomEffectsSet,
) -> None:
    """Oracle: the summary drawn: its tau, its excluded and zero-variance records.

    In spec the rise has no death before tau; out of spec it has no window at all.
    """
    for side, line in (
        (SIDE_IN_SPEC, "no death before τ, not pooled (1): rise"),
        (SIDE_OUT_OF_SPEC, "no windows (1): rise"),
    ):
        data = record_comparison.build_forest_data(
            forest_set,
            comparison="f",
            evidence=EVIDENCE,
            threshold_label="4",
            side=side,
        )
        fig = _draw(plot.ForestPlot, data)
        try:
            text = _flat(" ".join(_texts(fig)))
        finally:
            plt.close(fig)
        assert line in text, (side, text)
        for group, summary in data.summaries.items():
            tau = summary.tau
            if tau.rmst_tau_s is not None:
                assert (
                    f"τ = {tau.rmst_tau_s / 60.0:.3g} min, set by {tau.tau_set_by} "
                    f"(τ per group)"
                ) in text
            assert f"evidence: {EVIDENCE[group]}" in text
            assert _flat(f"check outcome: {summary.check_summary()}") in text


def test_figures_of_several_groups_say_they_are_described_not_tested(
    forest_set: record_comparison.RandomEffectsSet,
) -> None:
    """Oracle: specification. No step tests groups against each other, so a figure that
    sets them side by side says so; a figure of one group has nothing to say."""
    several = [
        (
            plot.ForestPlot,
            record_comparison.build_forest_data(
                forest_set,
                comparison="f",
                evidence=EVIDENCE,
                threshold_label="4",
                side=SIDE_IN_SPEC,
            ),
        ),
        (
            plot.PooledLadderPlot,
            record_comparison.build_pooled_ladder_data(
                forest_set, comparison="f", evidence=EVIDENCE
            ),
        ),
    ]
    for plot_class, data in several:
        fig = _draw(plot_class, data)
        try:
            assert "described, not tested" in _flat(" ".join(_texts(fig))), plot_class
        finally:
            plt.close(fig)
    one = record_comparison.random_effects_set(
        *[_carve(f"r{i}", 0.7 * i) for i in range(1, 4)],
        groups={"only": ["r1", "r2", "r3"]},
        tau_scope=PER_GROUP,
        tau_rule=RULE,
    )
    fig = _draw(
        plot.ForestPlot,
        record_comparison.build_forest_data(
            one,
            comparison="f",
            evidence={"only": "x"},
            threshold_label="4",
            side=SIDE_IN_SPEC,
        ),
    )
    try:
        assert "described, not tested" not in _flat(" ".join(_texts(fig)))
    finally:
        plt.close(fig)


def test_the_ladder_prints_tau_its_record_and_the_check_outcome_at_every_rung(
    forest_set: record_comparison.RandomEffectsSet,
) -> None:
    """Oracle: the summaries drawn; one line per rung names tau and the record that set it."""
    data = record_comparison.build_pooled_ladder_data(
        forest_set, comparison="f", evidence=EVIDENCE
    )
    fig = _draw(plot.PooledLadderPlot, data)
    try:
        text = _flat(" ".join(_texts(fig)))
    finally:
        plt.close(fig)
    for (group, label, side), summary in data.summaries.items():
        tau = summary.tau
        if tau.rmst_tau_s is None:
            assert f"{label}: no record with windows" in text
        else:
            assert (
                f"{label}: τ {tau.rmst_tau_s / 60.0:.3g} min, set by {tau.tau_set_by};"
            ) in text, (group, label, side)
        assert _flat(summary.check_summary()) in text
    small = data.summaries[("small", "4", SIDE_IN_SPEC)]
    assert f"not pooled, {small.k_pooled} with a variance" in text


def test_each_figure_draws_its_numbers_where_the_artifact_puts_them(
    pair_data: record_comparison.PairLadderData,
    forest_set: record_comparison.RandomEffectsSet,
) -> None:
    """Oracle: the artifact's own numbers, divided by 60 on every axis in minutes.

    One panel or group of each figure, at the second rung, so a figure drawn at the first
    rung fails: the pair's curves and tau line, the forest's squares, intervals, diamond,
    prediction interval and tau line, and the ladder's diamonds, bars and tau line.
    """
    panel = pair_data.panels[("6", SIDE_IN_SPEC)]
    fig = _draw(plot.PairLadderPlot, pair_data)
    try:
        ax = next(a for a in fig.axes if a.get_title() == "in spec, 6")
        steps = [line for line in ax.lines if line.get_drawstyle() == "steps-post"]
        assert len(steps) == 2
        for line, curve in zip(steps, (panel.curve_a, panel.curve_b), strict=True):
            # The steps, then at most one more point holding the last value out to the
            # curve's longest observed age (`KaplanMeierCurve.max_observed_min`).
            n = len(curve.time_min)
            x, y = line.get_xdata(), line.get_ydata()
            np.testing.assert_array_equal(x[:n], curve.time_min)
            np.testing.assert_array_equal(y[:n], curve.survival)
            assert len(x) <= n + 1 and y[-1] == curve.survival[-1]
            assert x[-1] == max(curve.max_observed_min, curve.time_min[-1])
        (tau,) = [line for line in ax.lines if list(line.get_ydata()) == [0, 1]]
        assert list(tau.get_xdata()) == [panel.pair.tau.rmst_tau_s / 60.0] * 2
        assert steps[1].get_linestyle() != tau.get_linestyle()
    finally:
        plt.close(fig)

    data = record_comparison.build_forest_data(
        forest_set,
        comparison="f",
        evidence=EVIDENCE,
        threshold_label="6",
        side=SIDE_IN_SPEC,
    )
    for group in forest_set.groups:
        assert data.summaries[group] is forest_set.get(group, "6", SIDE_IN_SPEC)
    big = data.summaries["big"]
    pooled = [row for row in big.rows if row.pooled]
    assert big.pooled_s is not None and pooled
    fig = _draw(plot.ForestPlot, data)
    try:
        ax = next(a for a in fig.axes if a.get_title(loc="left") == "big")
        squares = [line for line in ax.lines if line.get_marker() == "s"]
        assert [line.get_xdata()[0] for line in squares] == [
            row.summary.rmst_s / 60.0 for row in pooled
        ]
        bars = [
            list(line.get_xdata())
            for line in ax.lines
            if line.get_linestyle() == "-" and len(line.get_xdata()) == 2
        ]
        assert bars == [[row.lower_s / 60.0, row.upper_s / 60.0] for row in pooled]
        (dotted,) = [line for line in ax.lines if line.get_linestyle() == ":"]
        assert list(dotted.get_xdata()) == [
            big.prediction_lower_s / 60.0,
            big.prediction_upper_s / 60.0,
        ]
        (diamond,) = ax.patches
        np.testing.assert_array_equal(
            diamond.get_xy()[:4, 0],
            np.array(
                [big.pooled_lower_s, big.pooled_s, big.pooled_upper_s, big.pooled_s]
            )
            / 60.0,
        )
        (tau,) = [line for line in ax.lines if list(line.get_ydata()) == [0, 1]]
        assert list(tau.get_xdata()) == [big.tau.rmst_tau_s / 60.0] * 2
        text = _flat(" ".join(_texts(fig)))
        assert f"pooled, k = {big.k_pooled}: {big.pooled_s / 60.0:.3g} min" in text
        assert f"between-record variance {big.between_var_s2 / 3600.0:.3g} min²" in text
        assert "τ²" not in text
        dropped = sum(row.summary.n_unobserved_birth_dropped for row in big.rows)
        carved = dropped + sum(row.summary.n_windows for row in big.rows)
        assert dropped > 0
        assert f"birth unseen, excluded: {dropped} of {carved} windows" in text
    finally:
        plt.close(fig)

    ladder = record_comparison.build_pooled_ladder_data(
        forest_set, comparison="f", evidence=EVIDENCE
    )
    summaries = [
        ladder.summaries[("big", label, SIDE_IN_SPEC)] for label in ladder.ladder
    ]
    fig = _draw(plot.PooledLadderPlot, ladder)
    try:
        ax = next(a for a in fig.axes if a.get_title(loc="left") == "big, in spec")
        diamonds = [line for line in ax.lines if line.get_marker() == "D"]
        expected = [
            (float(i), s.pooled_s / 60.0)
            for i, s in enumerate(summaries)
            if s.pooled_s is not None
        ]
        assert expected
        assert [
            (line.get_xdata()[0], line.get_ydata()[0]) for line in diamonds
        ] == expected
        bars = [
            (list(line.get_xdata()), list(line.get_ydata()))
            for line in ax.lines
            if line.get_marker() != "D"
        ]
        assert bars == [
            pair
            for i, s in enumerate(summaries)
            if s.pooled_s is not None
            for pair in (
                ([i, i], [s.prediction_lower_s / 60.0, s.prediction_upper_s / 60.0]),
                ([i, i], [s.pooled_lower_s / 60.0, s.pooled_upper_s / 60.0]),
            )
        ]
        (twin,) = [
            a
            for a in fig.axes
            if a is not ax and a.get_position().bounds == ax.get_position().bounds
        ]
        (tau,) = twin.lines
        assert list(tau.get_ydata()) == [s.tau.rmst_tau_s / 60.0 for s in summaries]
        assert twin.get_ylabel() == "Comparison tau (min)"
        text = _flat(" ".join(_texts(fig)))
        for label, s in zip(ladder.ladder, summaries, strict=True):
            dropped = sum(row.summary.n_unobserved_birth_dropped for row in s.rows)
            carved = dropped + sum(row.summary.n_windows for row in s.rows)
            assert (
                f"{label}: τ {s.tau.rmst_tau_s / 60.0:.3g} min, set by "
                f"{s.tau.tau_set_by}; pooled {s.k_pooled} of {len(s.rows)}, "
                f"{len(s.tau.excluded)} without windows; birth unseen, excluded "
                f"{dropped} of {carved} windows"
            ) in text
    finally:
        plt.close(fig)


def test_the_drawn_rung_is_a_kwarg_of_the_figure_step_alone() -> None:
    """Oracle: specification. Drawing another rung must not change the comparison node,
    and a rung or side the ladder lacks is refused at graph build."""
    kwargs: dict[str, dict[str, object]] = {}
    for label in ("4", "6"):
        job = Job(f"rung_{label}")
        carves = {
            r: job.step(_carve, name=f"windows_{r}", record=r, phase=0.7 * i)
            for i, r in enumerate(("r1", "r2", "r3"))
        }
        refs = comparison_recipe.wire_record_forest(
            job,
            carves,
            name="f",
            groups={"g": ["r1", "r2", "r3"]},
            evidence={"g": "x"},
            tau_scope=PER_GROUP,
            tau_rule=RULE,
        )
        figures = comparison_recipe.wire_forest_figures(
            job,
            refs["forest"],
            name="f",
            evidence={"g": "x"},
            threshold_label=label,
            side=SIDE_IN_SPEC,
            thresholds=LADDER,
        )
        kwargs[label] = {
            "forest": job.dag[refs["forest"].node_id].kwargs,
            "figure": job.dag[figures["forest_data"].node_id].kwargs,
        }
        assert {sink.name for sink in job.sinks} >= {
            "f_forest_figure",
            "f_ladder_figure",
        }
        for bad in (
            {"threshold_label": "5", "side": SIDE_IN_SPEC},
            {"threshold_label": "4", "side": "inside"},
        ):
            with pytest.raises(ValueError, match="rung of the ladder|side must be"):
                comparison_recipe.wire_forest_figures(
                    job,
                    refs["forest"],
                    name="x",
                    evidence={"g": "x"},
                    thresholds=LADDER,
                    **bad,
                )
    assert kwargs["4"]["forest"] == kwargs["6"]["forest"]
    assert kwargs["4"]["figure"]["threshold_label"] == "4"
    assert kwargs["6"]["figure"]["threshold_label"] == "6"


def test_a_flagged_record_carries_a_dagger_and_its_flag_is_printed_beside_it(
    forest_set: record_comparison.RandomEffectsSet,
) -> None:
    """Oracle: specification. A flagged record's name carries the dagger wherever it is
    drawn (forest rows, a small group's table, the pair legend), an unflagged one does
    not, and each flag's line is printed with the records it marks."""
    flags = {"r1": "line one", "r4": "line one", "rise": "line rise"}
    data = record_comparison.build_forest_data(
        forest_set,
        comparison="f",
        evidence=EVIDENCE,
        threshold_label="4",
        side=SIDE_IN_SPEC,
        flags=flags,
    )
    fig = _draw(plot.ForestPlot, data)
    try:
        ticks = [t.get_text() for ax in fig.axes for t in ax.get_yticklabels()]
        text = _flat(" ".join(_texts(fig)))
    finally:
        plt.close(fig)
    assert "r1 †" in ticks and "r2" in ticks and "r2 †" not in ticks
    assert "rise † (no death before τ)" in ticks
    assert "r4 †: " in text and "r5 †" not in text
    for line in ("† line one: r1", "† line rise: rise", "† line one: r4"):
        assert line in text, line

    a, rise = _carve("a", 0.0), _rise("rise")
    pair = record_comparison.build_pair_ladder_data(
        record_comparison.compare_pair_set(a, rise, tau_rule=RULE),
        a,
        rise,
        comparison="p",
        evidence="one synthetic run",
        flags={"rise": "line rise"},
    )
    fig = _draw(plot.PairLadderPlot, pair)
    try:
        legend = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
        caption = _flat(" ".join(t.get_text() for t in fig.texts))
    finally:
        plt.close(fig)
    assert legend == ["a: a", "b: rise †", "comparison tau"]
    assert "† line rise: rise." in caption
