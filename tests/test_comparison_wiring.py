"""The comparisons across records wired as nodes, resolved on synthetic carves.

No real record is loaded: the carves are steps over noisy sines, and the 2x2-shaped record
is a synthetic pickle written to a temporary directory.
"""

from __future__ import annotations

import re
import zlib

import numpy as np
import pandas as pd
import pytest

from quebra import comparison_recipe, recipes
from quebra.analyzers import event_table, record_comparison, survival, windows
from quebra.analyzers.windows import SIDE_IN_SPEC, SIDE_OUT_OF_SPEC
from quebra.core.closure import parameter_row
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.core.reference import LocalRef
from quebra.core.runner import _toposort
from quebra.schemas.track912 import track912Schema

pytestmark = pytest.mark.integration

LADDER = [("4", 4.0, True), ("6", 6.0, True)]
RULE = record_comparison.TAU_RULE_MIN_LARGEST_AGE
# The rows the specification asks of each table, spelled here and not read off the module.
PAIR_QUANTITIES = (
    "tau_s",
    "rmst_difference_s",
    "rmst_difference_se_s",
    "logrank_z",
    "logrank_chi2",
    "logrank_p",
    "check_outcome",
)
RECORD_QUANTITIES = ("rmst_s", "rmst_se_s", "n_windows", "n_deaths")
POOLED_QUANTITIES = (
    "k_pooled",
    "pooled_s",
    "pooled_se_s",
    "pooled_lower_s",
    "pooled_upper_s",
    "between_var_s2",
    "between_var_lower_s2",
    "between_var_upper_s2",
    "q_statistic",
    "q_df",
    "i2_fraction",
    "prediction_lower_s",
    "prediction_upper_s",
)


def _carve(
    record: str,
    phase: float,
    period: float = 3.0,
    *,
    ladder: list[tuple[str, float, bool]] = LADDER,
    gap_mult: float = 10.0,
) -> windows.WindowsResult:
    """A noisy sine crossing both rungs many times each way: windows on every side.

    The noise, seeded by the record's name, spreads the windows' ages; a clean sine gives
    windows of one length, so no record has a death before tau and nothing pools.
    """
    t = np.arange(120.0)
    rng = np.random.default_rng(zlib.crc32(record.encode()))
    v = 5.0 + 3.0 * np.sin(t / period + phase) + rng.normal(0.0, 0.8, t.size)
    return windows.run(
        windows.WindowsInputs(
            t_rel_s=t, values=v, thresholds=ladder, dataset_id=record, gap_mult=gap_mult
        )
    )


def _rise(record: str) -> windows.WindowsResult:
    """Out of spec for 20 reads, then in spec to the end, at both rungs.

    One in-spec window, born at read 20 and censored, older than any sine's windows; the
    out-of-spec window starts with the record, so its birth is unseen and it is dropped.
    """
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


def _carves(job: Job, phases: dict[str, float]) -> dict[str, LocalRef]:
    return {
        record: job.step(_carve, name=f"windows_{record}", record=record, phase=phase)
        for record, phase in phases.items()
    }


def _resolve(job: Job, refs: dict[str, LocalRef]) -> dict[str, object]:
    results: dict[str, object] = {}
    for node_id in _toposort(job, [ref.node_id for ref in refs.values()]):
        node = job.dag[node_id]
        inputs = [results[ref.node_id] for ref in node.inputs]
        results[node_id] = node.fn(*inputs, **node.kwargs)
    return {name: results[ref.node_id] for name, ref in refs.items()}


def _rows(table: pd.DataFrame, label: str, side: str) -> pd.DataFrame:
    return table[(table["threshold_label"] == label) & (table["side"] == side)]


def _cell(rows: pd.DataFrame, record: str, quantity: str) -> tuple[object, object]:
    """The one (value, note) of `record`'s `quantity` among `rows`."""
    ((value, note),) = rows[
        (rows["record"] == record) & (rows["quantity"] == quantity)
    ][["value", "note"]].itertuples(index=False)
    return value, note


KEYS = {(label, side) for label, _, _ in LADDER for side in event_table.SIDES}


def test_a_pair_resolves_at_every_rung_and_side_and_its_table_holds_every_number() -> (
    None
):
    """Oracle: specification. One comparison per (label, side), each number on its own row."""
    job = Job("pair_wiring")
    carves = _carves(job, {"r1": 0.0, "r2": 1.0})
    refs = comparison_recipe.wire_record_pair(job, carves, name="R1", tau_rule=RULE)
    out = _resolve(job, refs)
    pairs, table = out["pair"], out["table"]

    assert {sink.name for sink in job.sinks} == {"R1_pair", "R1_table"}
    assert list(table.columns) == list(record_comparison.TABLE_COLUMNS)
    assert set(pairs.items) == KEYS
    assert set(zip(table["threshold_label"], table["side"])) == KEYS
    both = "r1 vs r2"
    for (label, side), pair in pairs.items.items():
        assert pair.not_computed_reason is None, (label, side)
        rows = _rows(table, label, side)
        assert sorted(zip(rows["record"], rows["quantity"])) == sorted(
            [(both, quantity) for quantity in PAIR_QUANTITIES]
            + [(r, quantity) for r in ("r1", "r2") for quantity in RECORD_QUANTITIES]
        )
        assert _cell(rows, both, "tau_s")[0] == pair.tau.rmst_tau_s
        for quantity in PAIR_QUANTITIES[1:-1]:
            assert _cell(rows, both, quantity)[0] == getattr(pair, quantity), quantity
        assert _cell(rows, both, "logrank_p")[1] == record_comparison.UNADJUSTED
        for record, summary in (("r1", pair.summary_a), ("r2", pair.summary_b)):
            for quantity in RECORD_QUANTITIES:
                assert _cell(rows, record, quantity)[0] == getattr(summary, quantity)

    # The rise has no out-of-spec window: that side is not computed, and its tau row says
    # so beside the tau the rule gave over r1 alone.
    lone = record_comparison.pair_table(
        record_comparison.compare_pair_set(
            _carve("r1", 0.0), _rise("rise"), tau_rule=RULE
        ),
        comparison="R1x",
    )
    for label, _, _ in LADDER:
        assert _cell(_rows(lone, label, SIDE_OUT_OF_SPEC), "r1 vs rise", "tau_s")[
            1
        ] == (
            f"{RULE}, set by r1; nothing computed at it "
            f"(rise: {record_comparison.NO_WINDOWS})"
        )


def test_a_forest_carries_its_groups_on_the_parameter_row_and_lists_a_small_group() -> (
    None
):
    """Oracle: `core.closure.parameter_row`, and the pooling floor of the specification.

    A group captured in a closure would change the numbers without changing the run's
    identity or its provenance label.
    """
    job = Job("forest_wiring")
    carves = _carves(job, {f"r{i}": 0.7 * i for i in range(1, 6)})
    groups = {"big": ["r1", "r2", "r3"], "small": ["r4", "r5"]}
    refs = comparison_recipe.wire_record_forest(
        job,
        carves,
        name="R4",
        groups=groups,
        tau_scope=record_comparison.TAU_PER_GROUP,
        tau_rule=RULE,
    )
    line = next(
        r for r in parameter_row(list(job.dag.values())) if r.startswith("R4_forest(")
    )
    for kwarg in (
        "groups={'big': ['r1', 'r2', 'r3'], 'small': ['r4', 'r5']}",
        "tau_scope='per_group'",
        "tau_rule='min_largest_age'",
        "conf_level=0.95",
    ):
        assert re.search(rf"[(,]{re.escape(kwarg)}[,)]", line), (kwarg, line)

    forest = _resolve(job, refs)["forest"]
    for label, side in KEYS:
        small = forest.get("small", label, side)
        assert [row.record for row in small.rows] == ["r4", "r5"]
        assert small.not_pooled_reason is not None and small.pooled_s is None
    assert any(forest.get("big", *key).pooled_s is not None for key in KEYS)


def test_a_forest_takes_tau_per_group_or_across_groups_as_declared() -> None:
    """Oracle: specification; tau computed here as the min over records of the largest age.

    The two groups' sines have different periods, so their windows' largest ages differ
    and the two scopes give different taus. Every row's restricted mean must be the one
    `survival.curve_summaries` gives at the scope's tau, not only the reported tau.
    """
    carves = [_carve(f"a{i}", 0.5 * i, period=3.0) for i in range(3)] + [
        _carve(f"b{i}", 0.5 * i, period=5.0) for i in range(3)
    ]
    groups = {"a": ["a0", "a1", "a2"], "b": ["b0", "b1", "b2"]}
    tables = {
        str(c.meta["dataset_id"]): event_table.event_tables_from_carve(c)
        for c in carves
    }

    def largest(record: str, label: str, side: str) -> float:
        return float(np.max(tables[record].table(label, side).age_s))

    per, across = (
        record_comparison.random_effects_set(
            *carves, groups=groups, tau_scope=scope, tau_rule=RULE
        )
        for scope in (
            record_comparison.TAU_PER_GROUP,
            record_comparison.TAU_ACROSS_GROUPS,
        )
    )
    differ = 0
    for label, side in KEYS:
        own = {
            g: min(largest(r, label, side) for r in names)
            for g, names in groups.items()
        }
        everyone = min(own.values())
        differ += own["a"] != own["b"]
        for group in groups:
            for result, tau in ((per, own[group]), (across, everyone)):
                summary = result.get(group, label, side)
                assert summary.tau.rmst_tau_s == tau
                assert [row.record for row in summary.rows] == groups[group]
                for row in summary.rows:
                    at_tau = survival.curve_summaries(
                        tables[row.record].table(label, side), rmst_tau_s=tau
                    )
                    assert (
                        row.summary.rmst_tau_s,
                        row.summary.rmst_s,
                        row.summary.rmst_se_s,
                    ) == (tau, at_tau.rmst_s, at_tau.rmst_se_s), (group, row.record)
    assert differ > 0, "every rung gave both groups one tau: the scopes went untested"


def test_a_forest_table_holds_every_number_and_every_none_names_its_reason() -> None:
    """Oracle: specification; each value is read back off the summary it flattens.

    The rise, in group big, is listed with no interval on the in-spec side and excluded on
    the out-of-spec side. Group lone holds only another rise, so on the out-of-spec side
    no record has windows and its tau is None.
    """
    job = Job("forest_table")
    carves = _carves(job, {f"r{i}": 0.7 * i for i in range(1, 4)})
    for record in ("rise", "rise2"):
        carves[record] = job.step(_rise, name=f"windows_{record}", record=record)
    refs = comparison_recipe.wire_record_forest(
        job,
        carves,
        name="R4",
        groups={"big": ["r1", "r2", "r3", "rise"], "lone": ["rise2"]},
        tau_scope=record_comparison.TAU_PER_GROUP,
        tau_rule=RULE,
    )
    out = _resolve(job, refs)
    forest, table = out["forest"], out["table"]

    assert list(table.columns) == list(record_comparison.TABLE_COLUMNS)
    assert len(forest.items) == 2 * len(KEYS)
    pooled = "pooled"
    for (group, label, side), summary in forest.items.items():
        tau = summary.tau
        if tau.rmst_tau_s is None:
            why = "; ".join(f"{record}: {reason}" for record, reason in tau.excluded)
            tau_note = f"{RULE}, per_group; no record with windows ({why})"
        else:
            tau_note = f"{RULE}, per_group, set by {tau.tau_set_by}"
        expected: dict[tuple[str, str], tuple[object, object]] = {
            (pooled, "tau_s"): (tau.rmst_tau_s, tau_note),
            (pooled, "conf_level"): (0.95, None),
            (pooled, "check_outcome"): (None, summary.check_summary()),
        }
        for quantity in POOLED_QUANTITIES:
            expected[(pooled, quantity)] = (
                getattr(summary, quantity),
                summary.not_pooled_reason,
            )
        for row in summary.rows:
            for quantity in RECORD_QUANTITIES:
                expected[(row.record, quantity)] = (
                    getattr(row.summary, quantity),
                    row.not_pooled_reason,
                )
            expected[(row.record, "lower_s")] = (row.lower_s, row.not_pooled_reason)
            expected[(row.record, "upper_s")] = (row.upper_s, row.not_pooled_reason)
        for record, reason in tau.excluded:
            expected[(record, "rmst_s")] = (None, reason)
        rows = _rows(table[table["group"] == group], label, side)
        assert sorted(zip(rows["record"], rows["quantity"])) == sorted(expected)
        for (record, quantity), cell in expected.items():
            assert _cell(rows, record, quantity) == cell, (group, label, side, quantity)

    # The loop above met every branch the table has.
    for label, _, _ in LADDER:
        big_in = forest.get("big", label, SIDE_IN_SPEC)
        (rise,) = [row for row in big_in.rows if row.record == "rise"]
        assert rise.lower_s is None
        assert rise.not_pooled_reason == record_comparison.ZERO_VARIANCE
        big_out = forest.get("big", label, SIDE_OUT_OF_SPEC)
        assert big_out.tau.excluded == (("rise", record_comparison.NO_WINDOWS),)
        assert forest.get("lone", label, SIDE_OUT_OF_SPEC).tau.rmst_tau_s is None
    assert any(forest.get("big", *key).pooled_s is not None for key in KEYS)


def test_a_comparison_refuses_halves_doubles_and_carves_that_disagree() -> None:
    """Oracle: specification. A record never meets its halves or itself, one label names
    one threshold, and tau is taken over every record pooled."""
    full, half = "120423_6D2S_qubit1", "120423_6D2S_qubit1_before"
    job = Job("refusals")
    with pytest.raises(ValueError, match="is a half of"):
        comparison_recipe.wire_record_pair(
            job, _carves(job, {full: 0.0, half: 1.0}), name="J", tau_rule=RULE
        )
    with pytest.raises(ValueError, match="is a half of"):
        record_comparison.compare_pair_set(
            _carve(full, 0.0), _carve(half, 1.0), tau_rule=RULE
        )
    with pytest.raises(ValueError, match="is a half of"):
        record_comparison.random_effects_set(
            _carve(full, 0.0),
            _carve(half, 1.0),
            groups={"g": [full], "h": [half]},
            tau_scope=record_comparison.TAU_PER_GROUP,
            tau_rule=RULE,
        )
    # The two halves of one record, without it, are a comparison across a jump.
    record_comparison.compare_pair_set(
        _carve(f"{full}_before", 0.0), _carve(f"{full}_after", 1.0), tau_rule=RULE
    )

    # A forest refuses at graph build, before any carve runs.
    scope = dict(tau_scope=record_comparison.TAU_PER_GROUP, tau_rule=RULE)
    forest_job = Job("forest_refusals")
    carves = _carves(forest_job, {full: 0.0, f"{full}_after": 1.0, "b": 2.0})
    with pytest.raises(ValueError, match="is a half of"):
        comparison_recipe.wire_record_forest(
            forest_job,
            carves,
            name="R4",
            groups={"g": [full, "b"], "h": [f"{full}_after"]},
            **scope,
        )
    with pytest.raises(ValueError, match="appears twice across the groups"):
        comparison_recipe.wire_record_forest(
            forest_job,
            carves,
            name="R4",
            groups={"g": [full, "b"], "h": [full]},
            **scope,
        )
    with pytest.raises(ValueError, match="two carves name one record"):
        record_comparison.random_effects_set(
            _carve("a", 0.0),
            _carve("a", 1.0),
            _carve("b", 2.0),
            groups={"g": ["a", "b"]},
            **scope,
        )
    tables = {
        r: event_table.event_tables_from_carve(_carve(r, 0.5 * i))
        for i, r in enumerate("abc")
    }
    label = LADDER[0][0]
    pooled = [tables[r].table(label, SIDE_IN_SPEC) for r in "abc"]
    with pytest.raises(ValueError, match="another threshold or side"):
        record_comparison.random_effects_summary(
            pooled,
            tau_rule=RULE,
            tau_from=[tables[r].table(label, SIDE_OUT_OF_SPEC) for r in "abc"],
        )
    with pytest.raises(ValueError, match="leave out"):
        record_comparison.random_effects_summary(
            pooled, tau_rule=RULE, tau_from=pooled[:2]
        )

    with pytest.raises(ValueError, match="different ladders"):
        record_comparison.compare_pair_set(
            _carve("a", 0.0), _carve("b", 1.0, ladder=LADDER[:1]), tau_rule=RULE
        )
    with pytest.raises(ValueError, match="two thresholds"):
        record_comparison.compare_pair_set(
            _carve("a", 0.0),
            _carve("b", 1.0, ladder=[("4", 4.5, True), ("6", 6.0, True)]),
            tau_rule=RULE,
        )
    with pytest.raises(ValueError, match="gap rules"):
        record_comparison.compare_pair_set(
            _carve("a", 0.0), _carve("b", 1.0, gap_mult=5.0), tau_rule=RULE
        )
    with pytest.raises(ValueError, match="run_name"):
        comparison_recipe.wire_record_carves(
            Job("renamed"),
            [Dataset(path=f"{full}.pickle", extra={"run_name": "q1_120423"})],
        )


def test_a_2x2_shaped_record_loads_through_track912_and_carves_on_the_t2star_ladder(
    tmp_path,
) -> None:
    """Oracle: `windows.run` on the record's own times and T2* column, called directly.

    The record is synthetic, in the shape the 2x2 records are declared with: a
    `DDMMYY_2x2_qubitN` pickle read through `track912Schema`, device "2x2". Load, filter,
    T2* fit and carve through `wire_record_carves` must give the same windows at every
    rung of the T2* ladder and on both sides.
    """
    seconds = np.arange(300) * 60.0
    t2star_s = 5.5e-6 + 5.0e-6 * np.sin(np.arange(300) / 6.0)
    frame = pd.DataFrame(
        {
            # 2023-07-03 08:00 UTC, the day the file name declares.
            "timestamp": 1688371200.0 + seconds,
            "amp": 1.0,
            "amp error": 0.02,
            "offset": 0.5,
            "offset error": 0.005,
            "T2star": t2star_s,
            "T2star error": 2e-7,
            "frequency": 2.1e6,
            "frequency error": 6e3,
            "phase": 0.1,
            "phase_error": 0.03,
            "normalised chi-square": 0.1,
        }
    )
    path = tmp_path / "030723_2x2_qubit1.pickle"
    frame.to_pickle(path)
    job = Job("carve_2x2")
    dataset = Dataset(path=str(path), schema=track912Schema, qubit=1, device="2x2")
    carves = comparison_recipe.wire_record_carves(job, [dataset])
    with pytest.warns(UserWarning, match="derived from DDMMYY filename prefix"):
        carved = _resolve(job, carves)["030723_2x2_qubit1"]
    direct = windows.run(
        windows.WindowsInputs(
            t_rel_s=seconds,
            values=t2star_s,
            thresholds=recipes.T2STAR_THRESHOLDS,
            gap_mult=recipes.T2STAR_GAP_MULT,
        )
    )

    ours = event_table.event_tables_from_carve(carved)
    theirs = event_table.event_tables_from_carve(direct)
    assert carved.meta["dataset_id"] == "030723_2x2_qubit1"
    assert ours.ladder == tuple(label for label, _, _ in recipes.T2STAR_THRESHOLDS)
    assert set(ours.tables) == set(theirs.tables)
    for key, table in theirs.tables.items():
        assert table.n_windows > 0, key
        np.testing.assert_array_equal(ours.tables[key].age_s, table.age_s)
        np.testing.assert_array_equal(ours.tables[key].event, table.event)
