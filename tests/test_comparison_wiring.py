"""The comparisons across records wired as nodes, resolved on synthetic carves.

No real record is loaded: the carves are steps over noisy sines, and the 2x2-shaped record
is a synthetic pickle written to a temporary directory. The comparison jobs are imported,
which builds their graphs and reads no record.
"""

from __future__ import annotations

import importlib.util
import re
import zlib
from dataclasses import replace
from datetime import date, datetime
from typing import Any
from pathlib import Path
from types import ModuleType

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
from quebra.paper_tables import FLAG_TEXT
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
EVIDENCE = "one synthetic run"
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
    refs = comparison_recipe.wire_record_pair(
        job,
        carves,
        name="R1",
        evidence=EVIDENCE,
        tau_rule=RULE,
        flags={"r2": "flagged synthetic run"},
    )
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
            [(both, quantity) for quantity in (*PAIR_QUANTITIES, "evidence")]
            + [(r, quantity) for r in ("r1", "r2") for quantity in RECORD_QUANTITIES]
            + [("r2", "flag")]
        )
        assert _cell(rows, both, "evidence") == (None, EVIDENCE)
        assert _cell(rows, "r2", "flag") == (None, "flagged synthetic run")
        assert _cell(rows, both, "tau_s")[0] == pair.tau.rmst_tau_s
        for quantity in PAIR_QUANTITIES[1:-1]:
            assert _cell(rows, both, quantity)[0] == getattr(pair, quantity), quantity
        assert _cell(rows, both, "logrank_p")[1] == record_comparison.UNADJUSTED
        # A tested pair's SE carries no note: only a descriptive pair's says why it is None.
        assert _cell(rows, both, "rmst_difference_se_s")[1] is None
        for record, summary in (("r1", pair.summary_a), ("r2", pair.summary_b)):
            for quantity in RECORD_QUANTITIES:
                assert _cell(rows, record, quantity)[0] == getattr(summary, quantity)

    # The rise has no out-of-spec window: that side is not computed, its tau row says so
    # beside the tau the rule gave over r1 alone, and both records' counts are still given,
    # counted here off r1's event table.
    lone = record_comparison.pair_table(
        record_comparison.compare_pair_set(
            _carve("r1", 0.0), _rise("rise"), tau_rule=RULE
        ),
        comparison="R1x",
        evidence=EVIDENCE,
    )
    r1 = event_table.event_tables_from_carve(_carve("r1", 0.0))
    why = f"rise: {record_comparison.NO_WINDOWS}"
    for label, _, _ in LADDER:
        rows = _rows(lone, label, SIDE_OUT_OF_SPEC)
        assert _cell(rows, "r1 vs rise", "tau_s")[1] == (
            f"{RULE}, set by r1; nothing computed at it ({why})"
        )
        seen = r1.table(label, SIDE_OUT_OF_SPEC)
        assert _cell(rows, "r1", "rmst_s") == (None, why)
        assert _cell(rows, "r1", "n_windows") == (len(seen.age_s), None)
        assert _cell(rows, "r1", "n_deaths") == (int(np.sum(seen.event)), None)
        assert _cell(rows, "rise", "n_windows") == (0, None)
        assert _cell(rows, "rise", "n_deaths") == (0, None)
        assert int(np.sum(seen.event)) > 0


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
        evidence={"big": EVIDENCE, "small": EVIDENCE},
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
    evidence = {"big": "four synthetic runs", "lone": "one synthetic run"}
    flags = {"r1": "flag one", "rise": "flag rise", "rise2": "flag rise"}
    refs = comparison_recipe.wire_record_forest(
        job,
        carves,
        name="R4",
        groups={"big": ["r1", "r2", "r3", "rise"], "lone": ["rise2"]},
        evidence=evidence,
        tau_scope=record_comparison.TAU_PER_GROUP,
        tau_rule=RULE,
        flags=flags,
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
            (pooled, "evidence"): (None, evidence[group]),
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
            expected[(record, "n_windows")] = (0, reason)
            expected[(record, "n_deaths")] = (0, reason)
        drawn = [row.record for row in summary.rows] + [r for r, _ in tau.excluded]
        for record in drawn:
            if record in flags:
                expected[(record, "flag")] = (None, flags[record])
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
            job,
            _carves(job, {full: 0.0, half: 1.0}),
            name="J",
            evidence=EVIDENCE,
            tau_rule=RULE,
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
            evidence={"g": EVIDENCE, "h": EVIDENCE},
            **scope,
        )
    with pytest.raises(ValueError, match="appears twice across the groups"):
        comparison_recipe.wire_record_forest(
            forest_job,
            carves,
            name="R4",
            groups={"g": [full, "b"], "h": [full]},
            evidence={"g": EVIDENCE, "h": EVIDENCE},
            **scope,
        )
    # Every group, and every pair, states its evidence: at graph build and in the step.
    with pytest.raises(ValueError, match="evidence needs one line per group"):
        comparison_recipe.wire_record_forest(
            forest_job,
            carves,
            name="R4",
            groups={"g": [full, "b"]},
            evidence={"h": EVIDENCE},
            **scope,
        )
    with pytest.raises(ValueError, match="evidence needs one line per group"):
        comparison_recipe.wire_record_pair(
            forest_job,
            {full: carves[full], "b": carves["b"]},
            name="P",
            evidence=" ",
            tau_rule=RULE,
        )
    three = [_carve(r, 0.5 * i) for i, r in enumerate("abc")]
    with pytest.raises(ValueError, match="evidence needs one line per group"):
        record_comparison.random_effects_table(
            record_comparison.random_effects_set(
                *three, groups={"g": ["a", "b", "c"]}, **scope
            ),
            comparison="R",
            evidence={"g": ""},
        )
    with pytest.raises(ValueError, match="evidence needs one line per group"):
        record_comparison.pair_table(
            record_comparison.compare_pair_set(*three[:2], tau_rule=RULE),
            comparison="P",
            evidence="",
        )
    # Each half of the evidence rule alone: a group with no line, a line with no group,
    # and a line that is not a string.
    forest_three = record_comparison.random_effects_set(
        *three, groups={"g": ["a", "b", "c"]}, **scope
    )
    for bad in ({}, {"g": EVIDENCE, "h": EVIDENCE}, {"g": None}):
        with pytest.raises(ValueError, match="evidence needs one line per group"):
            record_comparison.random_effects_table(
                forest_three, comparison="R", evidence=bad
            )
    # An excluded record's counts are 0 only because it has no windows: a record
    # excluded for any other reason is refused, not reported with 0.
    key = next(iter(forest_three.items))
    summary = forest_three.items[key]
    odd = replace(summary, tau=replace(summary.tau, excluded=(("d", "another"),)))
    with pytest.raises(ValueError, match="known here only when it has no windows"):
        record_comparison.random_effects_table(
            replace(forest_three, items={**forest_three.items, key: odd}),
            comparison="R",
            evidence={"g": EVIDENCE},
        )
    for bad in ({"z": "x"}, {"a": ""}, {"a": None}):
        with pytest.raises(ValueError, match="a flag names a record of the comparison"):
            comparison_recipe.wire_record_forest(
                forest_job,
                carves,
                name="R4",
                groups={"g": [full, "b"]},
                evidence={"g": EVIDENCE},
                flags={(full if k == "a" else k): v for k, v in bad.items()},
                **scope,
            )
    with pytest.raises(ValueError, match="a flag names a record of the comparison"):
        record_comparison.random_effects_table(
            forest_three, comparison="R", evidence={"g": EVIDENCE}, flags={"z": "x"}
        )
    with pytest.raises(ValueError, match="a flag names a record of the comparison"):
        record_comparison.pair_table(
            record_comparison.compare_pair_set(*three[:2], tau_rule=RULE),
            comparison="P",
            evidence=EVIDENCE,
            flags={"c": "x"},
        )
    with pytest.raises(ValueError, match="not a DDMMYY"):
        comparison_recipe.record_dataset("120423_6D2S_q1")
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


JOBS = Path(__file__).resolve().parents[1] / "jobs" / "composite"


def _job_module(stem: str) -> ModuleType:
    """Import one job file by path: its graph is built and no record is read."""
    spec = importlib.util.spec_from_file_location(f"_job_{stem}", JOBS / f"{stem}.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _node(job: Job, name: str) -> Any:
    (node,) = [node for node in job.dag.values() if node.fn_name == name]
    return node


def _carved(job: Job, name: str) -> list[str]:
    """The records whose carves feed the step `name`, in input order."""
    return [
        job.dag[ref.node_id].fn_name.removeprefix("windows_")
        for ref in _node(job, name).inputs
        if job.dag[ref.node_id].fn_name.startswith("windows_")
    ]


# The condition groups as declared at the jobs stop: device, first and last night, records.
SPANS = {
    "6D2S 200 mK April": ("6D2S", date(2023, 4, 4), date(2023, 4, 22), 11),
    "6D2S 300 mK June": ("6D2S", date(2023, 6, 7), date(2023, 6, 9), 6),
    "6D2S 200 mK June-July": ("6D2S", date(2023, 6, 27), date(2023, 7, 11), 7),
    "2x2 11 mK I": ("2x2", date(2023, 6, 9), date(2023, 6, 10), 4),
    "2x2 200 mK I": ("2x2", date(2023, 6, 20), date(2023, 6, 20), 1),
    "2x2 200 mK II": ("2x2", date(2023, 6, 25), date(2023, 7, 12), 13),
    "2x2 11 mK II": ("2x2", date(2023, 8, 22), date(2023, 8, 22), 4),
}


def test_the_comparison_jobs_wire_the_records_and_settings_of_the_specification() -> (
    None
):
    """Oracle: specification, the rows of the paper's dataset tables (Supp. Figs. 10, 11)
    grouped by the nights declared at the jobs stop, and the settings declared there.

    Read off each job's graph, not its constants: the records each step consumes, in
    order, the kwargs it carries, and the Dataset each record loads. R4's device groups must
    hold exactly R3's records, and no record left out at the stop, and no half, appears.
    """
    r1, r2, r3, r4 = (
        _job_module(stem).job
        for stem in (
            "compare_r1_qubits_one_night",
            "compare_r2_nights_one_qubit",
            "compare_r3_records_one_condition",
            "compare_r4_two_devices",
        )
    )
    settings = {"tau_rule": RULE, "conf_level": 0.95}
    for job, name, pair in (
        (r1, "r1", ["210423_6D2S_qubit1", "210423_6D2S_qubit2"]),
        (r2, "r2", ["070423_6D2S_qubit1", "100423_6D2S_qubit1"]),
    ):
        assert _carved(job, f"{name}_pair") == pair
        assert _node(job, f"{name}_pair").kwargs == settings
    assert _node(r1, "r1_table").kwargs["evidence"] == "one run, one voltage setting"

    forests = {}
    for job, name, scope in (
        (r3, "r3", record_comparison.TAU_PER_GROUP),
        (r4, "r4", record_comparison.TAU_ACROSS_GROUPS),
    ):
        kwargs = dict(_node(job, f"{name}_forest").kwargs)
        forests[name] = kwargs.pop("groups")
        assert kwargs == {**settings, "tau_scope": scope}
        assert sorted(_carved(job, f"{name}_forest")) == sorted(
            stem for stems in forests[name].values() for stem in stems
        )
    groups, devices = forests["r3"], forests["r4"]
    assert set(groups) == set(SPANS)
    for group, stems in groups.items():
        device, first, last, count = SPANS[group]
        assert len(set(stems)) == count, group
        for stem in stems:
            night, on, _ = stem.split("_")
            assert on == device, (group, stem)
            assert first <= datetime.strptime(night, "%d%m%y").date() <= last, (
                group,
                stem,
            )
    assert set(devices) == {"6D2S", "2x2"}
    for device in devices:
        assert sorted(devices[device]) == sorted(
            stem
            for group, stems in groups.items()
            if SPANS[group][0] == device
            for stem in stems
        )
    everyone = {stem for stems in groups.values() for stem in stems}
    left_out = {
        "070723_6D2S_qubit4",
        "100723_6D2S_qubit4",
        "070723_6D2S_qubit5",
        "070723_6D2S_qubit6",
    }
    assert len(everyone) == 46 and not everyone & left_out
    assert not any(stem.endswith(record_comparison.HALF_SUFFIXES) for stem in everyone)

    evidence = _node(r3, "r3_table").kwargs["evidence"]
    assert set(evidence) == set(groups)
    lines = {
        device: {evidence[g] for g in groups if SPANS[g][0] == device}
        for device in devices
    }
    assert lines["2x2"] == {
        "temperature and voltage configuration stated (Supp. Fig. S13)"
    }
    (six,) = lines["6D2S"]
    assert six == (
        "temperature period stated; voltage settings can change between nights and are "
        "not listed per night"
    )
    assert _node(r2, "r2_table").kwargs["evidence"] == six
    assert set(_node(r4, "r4_table").kwargs["evidence"]) == {"6D2S", "2x2"}

    # The rows the paper's dataset tables shade, 3 on the 6D2S and 11 on the 2x2, carry
    # one flag on the condition and device tables; no other record does.
    flag = FLAG_TEXT
    shaded = {
        "280623_6D2S_qubit1",
        "280623_6D2S_qubit2",
        "110723_6D2S_qubit3",
        "090623_2x2_qubit2",
        "090623_2x2_qubit3",
        "250623_2x2_qubit1",
        "300623_2x2_qubit1",
        "030723_2x2_qubit1",
        "030723_2x2_qubit2",
        "040723_2x2_qubit2",
        "050723_2x2_qubit1",
        "070723_2x2_qubit1",
        "100723_2x2_qubit1",
        "120723_2x2_qubit2",
    }
    flags = {stem: flag for stem in shaded}
    assert shaded <= everyone
    for job, name in ((r3, "r3"), (r4, "r4")):
        assert _node(job, f"{name}_table").kwargs["flags"] == flags
    for job, name in ((r1, "r1"), (r2, "r2")):
        assert _node(job, f"{name}_table").kwargs["flags"] == {}

    # Each record loads from its file, named by its stem: no run_name renames it.
    for job in (r1, r2, r3, r4):
        for node in job.dag.values():
            dataset = node.kwargs.get("dataset")
            if dataset is None:
                continue
            stem = dataset.path.stem
            _, device, qubit = stem.split("_")
            assert str(dataset.path) == f"data/real_private/{device}/{stem}.pickle"
            assert (dataset.device, dataset.qubit) == (device, int(qubit[5:]))
            assert dataset.schema is track912Schema and "run_name" not in dataset.extra

    # The figures: each pair beside its records' carves, each forest at the one rung and
    # side declared for it, and the pooled ladder; each figure-data step with the evidence
    # its table carries.
    for job, name, pair in (
        (r1, "r1", ["210423_6D2S_qubit1", "210423_6D2S_qubit2"]),
        (r2, "r2", ["070423_6D2S_qubit1", "100423_6D2S_qubit1"]),
    ):
        assert _carved(job, f"{name}_pair_data") == pair
        source, *_ = _node(job, f"{name}_pair_data").inputs
        assert job.dag[source.node_id].fn_name == f"{name}_pair"
        for kwarg in ("evidence", "flags"):
            assert (
                _node(job, f"{name}_pair_data").kwargs[kwarg]
                == _node(job, f"{name}_table").kwargs[kwarg]
            )
        assert f"{name}_pair_figure" in {sink.name for sink in job.sinks}
    for job, name in ((r3, "r3"), (r4, "r4")):
        evidence = _node(job, f"{name}_table").kwargs["evidence"]
        assert _node(job, f"{name}_forest_data").kwargs == {
            "comparison": name,
            "evidence": evidence,
            "threshold_label": "3 µs",
            "side": SIDE_IN_SPEC,
            "flags": flags,
        }
        assert _node(job, f"{name}_ladder_data").kwargs["evidence"] == evidence
        assert {f"{name}_forest_figure", f"{name}_ladder_figure"} <= {
            sink.name for sink in job.sinks
        }

    # J: each split record's halves, before then after, as a pair of their own; no full
    # record is carved beside them, and the shaded record's halves carry a half's flag.
    j = _job_module("compare_j_across_jump").job
    carved = set()
    for record in ("120423_6D2S_qubit1", "270623_6D2S_qubit2", "280623_6D2S_qubit2"):
        name = f"j_{record}"
        halves = [f"{record}_before", f"{record}_after"]
        assert record in everyone
        assert _carved(j, f"{name}_pair") == halves == _carved(j, f"{name}_pair_data")
        assert _node(j, f"{name}_pair").kwargs == settings
        table = _node(j, f"{name}_table").kwargs
        assert table["evidence"] == (
            "two halves of one run, split at a frequency jump: one voltage setting"
        )
        assert table["flags"] == (
            {half: f"half of a record {flag}" for half in halves}
            if record in shaded
            else {}
        )
        assert _node(j, f"{name}_pair_data").kwargs["flags"] == table["flags"]
        carved |= set(halves)
    carves_in_j = {
        node.fn_name.removeprefix("windows_")
        for node in j.dag.values()
        if node.fn_name.startswith("windows_")
    }
    assert carves_in_j == carved
    for node in j.dag.values():
        dataset = node.kwargs.get("dataset")
        if dataset is not None:
            assert (
                str(dataset.path)
                == f"data/real_private/6D2S/{dataset.path.stem}.pickle"
            )
