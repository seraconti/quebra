"""The record comparisons on the 2x2 infidelity windows: their records, their wiring, the
checks each record carries, and what every table and figure says beside a p.

The job graphs are built from the job files and read no record. The comparisons that run do so
on synthetic |x| carves, with ledgers the infidelity ledger step scores on them.
"""

from __future__ import annotations

import warnings
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from quebra.analyzers import fidelity, record_comparison  # noqa: E402
from dataclasses import replace  # noqa: E402

from quebra.analyzers.check_ledger import (  # noqa: E402
    GROUPS,
    VERDICT_FAIL,
    VERDICT_NOT_COMPUTED,
)
from quebra.analyzers.check_outcome import PERMUTATION_KEYS  # noqa: E402
from quebra.analyzers.comparison_checks import (  # noqa: E402
    _refused,
    forest_with_checks,
    pairs_with_checks,
)
from quebra.analyzers.windows import SIDE_IN_SPEC  # noqa: E402
from quebra.cli import _module_from_path  # noqa: E402
from quebra.core.discovery import by_family  # noqa: E402
from quebra.core.job import Job  # noqa: E402
from quebra.core.reference import ArtifactRef  # noqa: E402
from quebra.core.types import RecordFailure  # noqa: E402
from quebra.fidelity_ledger_recipe import fidelity_ledger_step  # noqa: E402
from quebra.ledger_recipe import LEDGER_KNOBS, ledger_step  # noqa: E402
from quebra.paper_tables import SHADED  # noqa: E402
from quebra.plots import record_comparison_plot as plot  # noqa: E402
from quebra.recipes import _fidelity_windows  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
JOBS = REPO / "jobs"
LADDER = list(fidelity.DETUNING_LADDER)
RULE = record_comparison.TAU_RULE_MIN_LARGEST_AGE
PERMUTATION_CHECKS = {check for check, _, _ in PERMUTATION_KEYS}
# Spelled here, not read off the module: the specification of the note beside every p.
P_NOTE = (
    "unadjusted for any other comparison; assumes the windows within each record are "
    "independent"
)
INDEPENDENT = "assumes the windows within each record are independent"


def _job(path: Path) -> Job:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return _module_from_path(path).job


def _nodes(job: Job, name: str) -> list:
    return [node for node in job.dag.values() if node.fn.__name__ == name]


def _ledger_records() -> dict[str, str]:
    """Every infidelity ledger job's records, run name to paper name, off their graphs."""
    out: dict[str, str] = {}
    for declared in by_family(JOBS, "check_ledger_fidelity"):
        for ledger in _nodes(_job(declared.path), "fidelity_ledger_step"):
            out[ledger.node_id.removeprefix("check_ledger_")] = ledger.kwargs[
                "paper_name"
            ]
    return out


def _record(ref: ArtifactRef, suffix: str) -> str:
    assert isinstance(ref, ArtifactRef), ref
    assert ref.included.job.name.startswith("check_ledger_fidelity_2x2_q"), ref
    assert ref.node_name.endswith(suffix), ref.node_name
    return ref.node_name.removesuffix(suffix)


def _carve(record: str, abs_x: np.ndarray) -> object:
    frame = pd.DataFrame(
        {
            "t_rel_s": np.arange(abs_x.size, dtype=float) * 90.0,
            "abs_x": abs_x,
            "sigma_x": np.full(abs_x.size, 2e-3),
        }
    )
    result = fidelity.FidelityResult(frame=frame, meta={"dataset_id": record})
    return _fidelity_windows(result, 10.0, 1.0, True, LADDER)


def _wavy(
    record: str, seed: int, centre: float = 0.06, swing: float = 0.08, n: int = 400
) -> object:
    """An |x| record crossing the nines rungs many times each way; `n` reads."""
    rng = np.random.default_rng(seed)
    x = centre + swing * np.sin(np.arange(n) / 7.0) + rng.normal(0, 0.03, n)
    return _carve(record, np.abs(x))


def _ledger(carve: object, paper: str) -> object:
    return fidelity_ledger_step(
        carve,
        pd.read_csv(REPO / "jobs/bench/results/size_table.csv"),
        paper_name=paper,
        thresholds=LADDER,
        **dict(LEDGER_KNOBS, n_permutations=19),
    )


@pytest.fixture(scope="module")
def scored() -> dict[str, tuple[object, object]]:
    """Three synthetic records, each its carve and its infidelity ledger. `fid_calm` stays
    in spec at 99% from its first read, so it has no window there (its only in-spec window
    opens at the first read and is dropped as an unseen start), and crosses 99.9%."""
    carves = {
        "fid_a": _wavy("fid_a", 7),
        "fid_b": _wavy("fid_b", 8),
        "fid_calm": _wavy("fid_calm", 9, centre=0.055, swing=0.02),
    }
    return {
        r: (c, _ledger(c, f"2x2_q9_{i}")) for i, (r, c) in enumerate(carves.items())
    }


def _by_hand(ledger, label: str, side: str) -> tuple[list[str], list[tuple]]:
    """(the permutation rows' verdicts of that rung and bag on both clocks, every group-row
    verdict of that rung), read straight off the ledger's rows in order. An `(all)` row,
    written where nothing was scored, gives each check of the run-set and each group its
    verdict."""
    rows = ledger.rows[ledger.rows["threshold_label"] == label]
    bag: list[str] = []
    groups: list[tuple] = []
    for r in rows.itertuples():
        if r.check_id == "(all)":
            if r.side == side:
                bag += [r.verdict] * len(PERMUTATION_CHECKS)
            groups += [(g, r.side, r.clock, r.verdict) for g in GROUPS]
        elif r.check_id in GROUPS:
            groups.append((r.check_id, r.side, r.clock, r.verdict))
        elif r.side == side and r.check_id in PERMUTATION_CHECKS:
            bag.append(r.verdict)
    return sorted(bag), groups


def _tally(verdicts: list[str]) -> str:
    """`n verdict` per distinct verdict, sorted by verdict, comma-joined (specification)."""
    return (
        ", ".join(f"{verdicts.count(v)} {v}" for v in sorted(set(verdicts))) or "none"
    )


def _note(groups: list[tuple]) -> str:
    """The group-verdict note, spelled from its specification."""
    return "; ".join(
        f"{g}, {bag} bag, {clock} clock: {v}" for g, bag, clock, v in groups
    )


@pytest.mark.integration
def test_f1_pairs_are_every_pair_of_2x2_records_sharing_a_night(in_repo) -> None:
    """Oracle: the rule over the infidelity ledger jobs' records, and the owner's list of
    12 pairs. F1 compares exactly every pair of records sharing a run date, each on the
    ledger jobs' own carves, a the lower qubit."""
    paper = _ledger_records()
    by_night: dict[str, list[str]] = {}
    for run in sorted(paper, key=lambda r: paper[r]):
        by_night.setdefault(paper[run].rsplit("_", 1)[1], []).append(run)
    rule = {
        pair for runs in by_night.values() for pair in combinations(sorted(runs), 2)
    }
    job = _job(JOBS / "composite/compare_f1_2x2_same_night.py")
    compared = {
        tuple(_record(ref, "_windows") for ref in node.inputs)
        for node in _nodes(job, "compare_pair_set")
    }
    assert compared == rule
    owners = {
        ("2x2_q1_030723", "2x2_q2_030723"),
        ("2x2_q2_040723", "2x2_q3_040723"),
        ("2x2_q1_090623", "2x2_q2_090623"),
        ("2x2_q1_090623", "2x2_q3_090623"),
        ("2x2_q2_090623", "2x2_q3_090623"),
        ("2x2_q1_300623", "2x2_q3_300623"),
        ("2x2_q1_220823", "2x2_q2_220823"),
        ("2x2_q1_220823", "2x2_q3_220823"),
        ("2x2_q1_220823", "2x2_q4_220823"),
        ("2x2_q2_220823", "2x2_q3_220823"),
        ("2x2_q2_220823", "2x2_q4_220823"),
        ("2x2_q3_220823", "2x2_q4_220823"),
    }
    assert {(paper[a], paper[b]) for a, b in compared} == owners


@pytest.mark.integration
def test_f3_groups_are_r3s_2x2_groups(in_repo) -> None:
    """Oracle: R3's graph. F3's groups hold the records of R3's four 2x2 groups, by their
    names in the paper's tables, in R3's order; the 11 mK II evidence names its drive."""
    paper = _ledger_records()
    (r3,) = _nodes(
        _job(JOBS / "composite/compare_r3_records_one_condition.py"),
        "random_effects_set",
    )
    r3_groups = {
        group: [
            "2x2_q{}_{}".format(stem.rsplit("qubit", 1)[1], stem.split("_")[0])
            for stem in stems
        ]
        for group, stems in r3.kwargs["groups"].items()
        if group.startswith("2x2 ")
    }
    job = _job(JOBS / "composite/compare_f3_2x2_conditions.py")
    (f3,) = _nodes(job, "random_effects_set")
    f3_groups = {
        group: [paper[run] for run in runs]
        for group, runs in f3.kwargs["groups"].items()
    }
    assert f3_groups == r3_groups
    (table,) = _nodes(job, "random_effects_table")
    assert "~2.2 MHz drive" in table.kwargs["evidence"]["2x2 11 mK II"]
    assert all(
        "2.2 MHz" not in line
        for group, line in table.kwargs["evidence"].items()
        if group != "2x2 11 mK II"
    )


@pytest.mark.integration
def test_every_comparison_reads_its_records_carve_and_ledger_from_the_ledger_jobs(
    in_repo,
) -> None:
    """Oracle: specification. F1 and F3 include exactly the four infidelity ledger jobs and
    carve and score nothing: each check step reads the comparison, then each record's
    persisted carve, then each record's persisted ledger, the same record each time, with
    the band's run-set; each table reads the filled comparison; F1's table and figure say
    how its pairs were declared."""
    ledger_jobs = sorted(d.job_id for d in by_family(JOBS, "check_ledger_fidelity"))
    f1 = _job(JOBS / "composite/compare_f1_2x2_same_night.py")
    f3 = _job(JOBS / "composite/compare_f3_2x2_conditions.py")
    for job in (f1, f3):
        assert sorted(inc.job.name for inc in job.includes) == ledger_jobs
        assert not _nodes(job, "_fidelity_windows") and not _nodes(
            job, "fidelity_ledger_step"
        )
    declared = _job_module_constant("compare_f1_2x2_same_night", "DECLARED")
    for checks in _nodes(f1, "pairs_with_checks"):
        compare, *refs = checks.inputs
        assert f1.dag[compare.node_id].fn.__name__ == "compare_pair_set"
        carves = [_record(ref, "_windows") for ref in refs[:2]]
        ledgers = [_record(ref, "_check_ledger_data") for ref in refs[2:]]
        assert carves == ledgers
        assert list(f1.dag[compare.node_id].inputs) == refs[:2]
        assert checks.kwargs == {"run_set": PERMUTATION_KEYS}
    assert len(_nodes(f1, "pairs_with_checks")) == 12
    for table in _nodes(f1, "pair_table"):
        (source,) = table.inputs
        assert f1.dag[source.node_id].fn.__name__ == "pairs_with_checks"
        assert table.kwargs["declared"] == declared
    for data in _nodes(f1, "build_pair_ladder_data"):
        assert f1.dag[data.inputs[0].node_id].fn.__name__ == "pairs_with_checks"
        assert data.kwargs["declared"] == declared
    (checks,) = _nodes(f3, "forest_with_checks")
    forest, *refs = checks.inputs
    records = list(checks.kwargs["records"])
    assert len(records) == 22
    assert [_record(ref, "_windows") for ref in refs[:22]] == records
    assert [_record(ref, "_check_ledger_data") for ref in refs[22:]] == records
    assert list(f3.dag[forest.node_id].inputs) == refs[:22]
    (table,) = _nodes(f3, "random_effects_table")
    assert table.inputs[0].node_id == checks.node_id
    # A ledger for a record no group names is refused at graph build, as for a pair.
    from quebra.comparison_recipe import wire_record_forest

    probe = Job("forest_ledgers")
    nodes = {r: probe.step(_wavy, name=f"carve_{r}", record=r, seed=1) for r in "abcd"}
    with pytest.raises(KeyError, match=r"extra \['d'\]"):
        wire_record_forest(
            probe,
            {r: nodes[r] for r in "abc"},
            name="probe",
            groups={"g": ["a", "b", "c"]},
            evidence={"g": "synthetic"},
            tau_scope=record_comparison.TAU_PER_GROUP,
            tau_rule=RULE,
            ledgers=nodes,
        )


def _job_module_constant(stem: str, name: str) -> object:
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        f"_const_{stem}", JOBS / f"composite/{stem}.py"
    )
    module = importlib.util.module_from_spec(spec)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec.loader.exec_module(module)
    return getattr(module, name)


@pytest.mark.unit
def test_a_pairs_checks_are_each_records_own_ledger_rows(scored) -> None:
    """Oracle: the ledgers' rows, read by hand. At every rung and side, each record's check
    outcome holds the five permutation checks, and its cells are that rung and bag's
    permutation rows on both clocks; its group verdicts are the ledger's group rows of that
    rung, 2 groups x 2 bags x 2 clocks; each summary carries its record's outcome; no number
    of the comparison changes."""
    (carve_a, ledger_a), (carve_b, ledger_b) = scored["fid_a"], scored["fid_b"]
    plain = record_comparison.compare_pair_set(carve_a, carve_b, tau_rule=RULE)
    filled = pairs_with_checks(
        plain, carve_a, carve_b, ledger_a, ledger_b, run_set=PERMUTATION_KEYS
    )
    assert set(filled.items) == set(plain.items)
    for (label, side), pair in filled.items.items():
        for checks, groups, summary, ledger in (
            (pair.checks_a, pair.group_verdicts_a, pair.summary_a, ledger_a),
            (pair.checks_b, pair.group_verdicts_b, pair.summary_b, ledger_b),
        ):
            bag, group_rows = _by_hand(ledger, label, side)
            assert len(checks.checks_asked) == len(PERMUTATION_CHECKS)
            assert sorted(v for *_, v in checks.check_verdicts) == bag
            assert list(groups) == group_rows and len(groups) == 8
            assert summary is None or summary.checks == checks
        before = plain.get(label, side)
        for quantity in ("logrank_p", "rmst_difference_s", "rmst_difference_se_s"):
            assert getattr(pair, quantity) == getattr(before, quantity), quantity
    table = record_comparison.pair_table(filled, comparison="F", evidence="synthetic")
    data = record_comparison.build_pair_ladder_data(
        filled, carve_a, carve_b, comparison="F", evidence="synthetic"
    )
    for (label, side), panel in data.panels.items():
        rows = table[(table["threshold_label"] == label) & (table["side"] == side)]
        expected_lines = []
        for tag, record, ledger in (("a", "fid_a", ledger_a), ("b", "fid_b", ledger_b)):
            bag, group_rows = _by_hand(ledger, label, side)
            notes = rows[
                (rows["record"] == record) & (rows["quantity"] == "group_verdicts")
            ]
            assert list(notes["note"]) == [_note(group_rows)], (label, side, record)
            expected_lines.append(
                f"{tag}: checks, {side} bag: {_tally(bag)}; "
                f"groups, both bags: {_tally([v for *_, v in group_rows])}"
            )
        assert list(panel.check_lines) == expected_lines


@pytest.mark.unit
def test_a_forests_rows_carry_each_records_own_ledger_rows(scored) -> None:
    """Oracle: the ledgers' rows, read by hand. Every row of every group, threshold and side
    carries its record's band cells and group rows of that rung; the table prints the group
    verdicts beside the record's numbers; no pooled number changes."""
    records = ("fid_a", "fid_b", "fid_calm")
    carves = [scored[r][0] for r in records]
    ledgers = [scored[r][1] for r in records]
    plain = record_comparison.random_effects_set(
        *carves,
        groups={"one": list(records)},
        tau_scope=record_comparison.TAU_PER_GROUP,
        tau_rule=RULE,
    )
    filled = forest_with_checks(
        plain, *carves, *ledgers, records=records, run_set=PERMUTATION_KEYS
    )
    by_name = dict(zip(records, ledgers, strict=True))
    n_rows = 0
    for (group, label, side), summary in filled.items.items():
        assert summary.pooled_s == plain.get(group, label, side).pooled_s
        for row in summary.rows:
            bag, group_rows = _by_hand(by_name[row.record], label, side)
            assert sorted(v for *_, v in row.summary.checks.check_verdicts) == bag
            assert list(row.group_verdicts) == group_rows
            n_rows += 1
    assert n_rows
    table = record_comparison.random_effects_table(
        filled, comparison="F", evidence={"one": "synthetic"}
    )
    notes = table[table["quantity"] == "group_verdicts"]
    assert len(notes) == n_rows
    assert all(note.count(": ") == 8 for note in notes["note"])
    for (group, label, side), summary in filled.items.items():
        fail, group_fail = [], []
        for row in summary.rows:
            bag, group_rows = _by_hand(by_name[row.record], label, side)
            if VERDICT_FAIL in bag:
                fail.append(row.record)
            if VERDICT_FAIL in [v for *_, v in group_rows]:
                group_fail.append(row.record)
        assert summary.check_line == (
            f"{len(summary.rows)} records checked, every cell in the table; a fail cell: "
            f"{', '.join(fail) or 'none'}; a failing group row: "
            f"{', '.join(group_fail) or 'none'}"
        )
    evidence = {"one": "synthetic"}
    forest_data = record_comparison.build_forest_data(
        filled,
        comparison="F",
        evidence=evidence,
        threshold_label="99.9%",
        side=SIDE_IN_SPEC,
    )
    ladder_data = record_comparison.build_pooled_ladder_data(
        filled, comparison="F", evidence=evidence
    )
    for plot_class, data, line in (
        (
            plot.ForestPlot,
            forest_data,
            filled.get("one", "99.9%", SIDE_IN_SPEC).check_line,
        ),
        (
            plot.PooledLadderPlot,
            ladder_data,
            filled.get("one", "99%", SIDE_IN_SPEC).check_line,
        ),
    ):
        fig = plot_class(name="unit").build_matplotlib(data)
        try:
            text = " ".join(
                " ".join(t.get_text().split())
                for t in [*fig.texts, *(t for ax in fig.axes for t in ax.texts)]
            )
            assert " ".join(line.split()) in text, (plot_class.__name__, text)
            assert "checks asked" not in text, (plot_class.__name__, text)
        finally:
            plt.close(fig)
    with pytest.raises(ValueError, match="differ from the forest's"):
        forest_with_checks(
            plain,
            *carves[:2],
            *ledgers[:2],
            records=records[:2],
            run_set=PERMUTATION_KEYS,
        )


@pytest.mark.unit
def test_a_carve_is_held_to_its_own_ledger_and_a_refused_ledger_reads_not_computed(
    scored,
) -> None:
    """Oracle: specification. A ledger given for another record's carve raises. A ledger
    that recorded a refusal of its record scored no window: it is held to its record only,
    and its record's cells and its 8 group verdicts per rung read `not computed`."""
    (carve_a, ledger_a), (carve_b, ledger_b) = scored["fid_a"], scored["fid_b"]
    plain = record_comparison.compare_pair_set(carve_a, carve_b, tau_rule=RULE)
    with pytest.raises(ValueError, match="the ledger of 'fid_b' given for the carve"):
        pairs_with_checks(
            plain, carve_a, carve_b, ledger_b, ledger_a, run_set=PERMUTATION_KEYS
        )
    refused = ledger_step(
        RecordFailure(dataset_id="fid_b", stage="check_ledger", reason="refused"),
        pd.read_csv(REPO / "jobs/bench/results/size_table.csv"),
        thresholds=LADDER,
        **LEDGER_KNOBS,
    )
    filled = pairs_with_checks(
        plain, carve_a, carve_b, ledger_a, refused, run_set=PERMUTATION_KEYS
    )
    for pair in filled.items.values():
        assert {v for *_, v in pair.checks_b.check_verdicts} == {VERDICT_NOT_COMPUTED}
        assert len(pair.group_verdicts_b) == 8
        assert {v for *_, v in pair.group_verdicts_b} == {VERDICT_NOT_COMPUTED}
    # A carve of the same record that is not the one its ledger scored: fewer reads.
    short_a = _wavy("fid_a", 7, n=300)
    with pytest.raises(ValueError, match="the carve holds"):
        pairs_with_checks(
            plain, short_a, carve_b, ledger_a, ledger_b, run_set=PERMUTATION_KEYS
        )
    forest = record_comparison.random_effects_set(
        short_a,
        carve_b,
        groups={"one": ["fid_a", "fid_b"]},
        tau_scope=record_comparison.TAU_PER_GROUP,
        tau_rule=RULE,
    )
    with pytest.raises(ValueError, match="the carve holds"):
        forest_with_checks(
            forest,
            short_a,
            carve_b,
            ledger_a,
            ledger_b,
            records=("fid_a", "fid_b"),
            run_set=PERMUTATION_KEYS,
        )
    # The refusal's exact shape: one shared `<stage> failed: <reason>` note on `(all)` rows
    # that used no window. A ledger that scored windows, or another note, is not one.
    assert _refused(refused) and not _refused(ledger_a)
    rows = refused.rows.copy()
    rows.loc[rows.index[0], "notes"] = "no segment left"
    assert not _refused(replace(refused, rows=rows))
    rows = refused.rows.copy()
    rows["notes"] = "no segment left"
    assert not _refused(replace(refused, rows=rows))
    rows = refused.rows.copy()
    rows.loc[rows.index[0], "n_windows_used"] = 3
    assert not _refused(replace(refused, rows=rows))


@pytest.mark.unit
def test_a_record_with_no_window_at_a_rung_is_reported_with_its_checks(scored) -> None:
    """Oracle: specification. `fid_calm` has no window in spec at 99%: the comparison there
    is not computed and says why, in its table rows and its figure panel, while the record
    still carries its check outcome and group verdicts; nothing raises."""
    (carve_a, ledger_a), (carve_c, ledger_c) = scored["fid_a"], scored["fid_calm"]
    plain = record_comparison.compare_pair_set(carve_a, carve_c, tau_rule=RULE)
    filled = pairs_with_checks(
        plain, carve_a, carve_c, ledger_a, ledger_c, run_set=PERMUTATION_KEYS
    )
    pair = filled.get("99%", SIDE_IN_SPEC)
    assert pair.n_windows_b == 0
    assert record_comparison.NO_WINDOWS in str(pair.not_computed_reason)
    assert pair.checks_b.checks_asked and len(pair.group_verdicts_b) == 8
    table = record_comparison.pair_table(filled, comparison="F", evidence="synthetic")
    cell = table[
        (table["threshold_label"] == "99%")
        & (table["side"] == SIDE_IN_SPEC)
        & (table["quantity"] == "rmst_difference_s")
    ]
    assert record_comparison.NO_WINDOWS in cell["note"].iloc[0]
    data = record_comparison.build_pair_ladder_data(
        filled, carve_a, carve_c, comparison="F", evidence="synthetic"
    )
    fig = plot.PairLadderPlot(name="unit").build_matplotlib(data)
    try:
        panel = [ax for ax in fig.axes if ax.get_title() == "in spec, 99%"]
        text = " ".join(" ".join(t.get_text().split()) for t in panel[0].texts)
        assert "not computed" in text and "b: checks" in text, text
    finally:
        plt.close(fig)


@pytest.mark.unit
def test_every_log_rank_p_and_the_caption_say_the_windows_are_assumed_independent(
    scored,
) -> None:
    """Oracle: specification. Every computed log-rank p in a pair table carries the note
    that it is unadjusted and assumes independent windows within each record; the figure's
    caption says the same; a declared line appears in the table and the caption only when
    given, and a panel prints check tallies only when the checks are filled."""
    (carve_a, ledger_a), (carve_b, ledger_b) = scored["fid_a"], scored["fid_b"]
    plain = record_comparison.compare_pair_set(carve_a, carve_b, tau_rule=RULE)
    filled = pairs_with_checks(
        plain, carve_a, carve_b, ledger_a, ledger_b, run_set=PERMUTATION_KEYS
    )
    for pairs, declared in ((plain, ""), (filled, "declared after a look")):
        table = record_comparison.pair_table(
            pairs, comparison="F", evidence="synthetic", declared=declared
        )
        p_rows = table[table["quantity"] == "logrank_p"]
        computed = p_rows[p_rows["value"].notna()]
        assert len(computed) and set(computed["note"]) == {P_NOTE}
        declared_rows = table[table["quantity"] == "declared"]
        assert set(declared_rows["note"]) == ({declared} if declared else set())
        data = record_comparison.build_pair_ladder_data(
            pairs,
            carve_a,
            carve_b,
            comparison="F",
            evidence="synthetic",
            declared=declared,
        )
        assert all(
            bool(p.check_lines) == (pairs is filled) for p in data.panels.values()
        )
        fig = plot.PairLadderPlot(name="unit").build_matplotlib(data)
        try:
            text = " ".join(" ".join(t.get_text().split()) for t in fig.texts)
            assert INDEPENDENT in text, text
            assert (declared in text) if declared else ("declared" not in text)
            pointer = "each panel's tallies, every cell in the table"
            assert (pointer in text) == (pairs is filled), text
        finally:
            plt.close(fig)


@pytest.mark.integration
def test_flags_mark_exactly_the_shaded_records_and_f3s_rung_is_the_ledgers(
    in_repo,
) -> None:
    """Oracle: the paper's 2x2 table (`paper_tables.SHADED`) and the infidelity ledger
    jobs' graphs. Every F1 and F3 table and figure flags exactly the records the paper
    shades, 11 of F3's 22; F3's forest is drawn at a rung of the ladder the ledgers scored,
    which F3 declares nowhere."""
    paper = _ledger_records()
    shaded = {run for run, name in paper.items() if name in SHADED}
    assert len(shaded) == 11
    f1 = _job(JOBS / "composite/compare_f1_2x2_same_night.py")
    for table in _nodes(f1, "pair_table"):
        compare = f1.dag[f1.dag[table.inputs[0].node_id].inputs[0].node_id]
        pair = {_record(ref, "_windows") for ref in compare.inputs}
        assert set(table.kwargs["flags"]) == pair & shaded, pair
    f3 = _job(JOBS / "composite/compare_f3_2x2_conditions.py")
    for name in ("random_effects_table", "build_forest_data"):
        (node,) = _nodes(f3, name)
        assert set(node.kwargs["flags"]) == shaded, name
    ladders = {
        tuple(label for label, _, _ in ledger.kwargs["thresholds"])
        for declared in by_family(JOBS, "check_ledger_fidelity")
        for ledger in _nodes(_job(declared.path), "fidelity_ledger_step")
    }
    ((labels),) = ladders
    (forest_data,) = _nodes(f3, "build_forest_data")
    assert forest_data.kwargs["threshold_label"] in labels
