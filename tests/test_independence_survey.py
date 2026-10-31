"""The check-ledger jobs score the windows the panels draw, and the survey re-scores nothing.

Two halves. The graph half imports the `check_ledger` jobs and the device-wide survey and
reads their DAGs: every ledger scores the T2* carve the T2* panel is drawn from, a ledger job
runs checks and nothing else, and the survey stacks every ledger of the family without
carving or scoring. A carve parameter drifting between a ledger job and a T2* job is silent
and produces plausible figures, so these read the graphs rather than the comments. The
analyzer half pins `build_independence_survey`'s reshape on fake ledgers.

Importing a job builds its graph and reads no data, so the graph half needs no data root.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from quebra.analyzers.check_ledger import (
    VERDICT_FAIL,
    VERDICT_INCONCLUSIVE,
    VERDICT_NOT_COMPUTED,
    VERDICT_PASS,
)
from quebra.analyzers.checks.battery import ROW_KEYS
from quebra.analyzers.independence_survey import (
    CHECK_LABELS,
    CHECK_NULL,
    TREND_CHECKS,
    InstrumentGrid,
    build_independence_survey,
    survey_summary,
    trend_group_counts,
)
from quebra.core.dataset import Dataset
from quebra.core.discovery import by_family
from quebra.core.job import Job
from quebra.core.reference import ArtifactRef
from quebra.ledger_recipe import (
    LEDGER_DATA_SUFFIX,
    LEDGER_KNOBS,
    collect_ledgers,
    included_ledger_set,
    ledger_step,
    record_label,
    survey_ledgers,
)
from quebra.panels.check_ledger import CheckLedgerPanel
from quebra.plots.independence_survey_plot import (
    SURVEY_PLOTS,
    IndependenceSurveyOverviewPlot,
)
from quebra.recipes import _windows_run

pytestmark = pytest.mark.unit


REPO = Path(__file__).resolve().parents[1]
JOBS = REPO / "jobs"


def _load_module(path: Path) -> object:
    spec = importlib.util.spec_from_file_location(f"_survey_jobmod_{path.stem}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_job(path: Path) -> Job:
    return _load_module(path).job


def _ledger_jobs() -> dict[str, Job]:
    family = by_family(JOBS, "check_ledger")
    assert len(family) == 6, [j.job_id for j in family]
    return {j.job_id: _load_job(j.path) for j in family}


def _nodes(job: Job, fn: object) -> list[object]:
    return [node for node in job.dag.values() if node.fn is fn]


def _chain(job: Job, node_id: str) -> tuple[Dataset, list[tuple]]:
    """The dataset and the step-by-step signature of the chain that ends at `node_id`.

    Each step is its function's code, the values its closure captured (the filter step
    captures its config there, not as a kwarg) and its kwargs. Node NAMES are left out: a
    ledger job suffixes them per record, and the name decides nothing about the windows.
    """
    steps: list[tuple] = []
    node = job.dag[node_id]
    while node.inputs:
        (parent,) = node.inputs
        cells = tuple(c.cell_contents for c in (node.fn.__closure__ or ()))
        steps.append((node.fn.__code__, cells, sorted(node.kwargs.items())))
        node = job.dag[parent.node_id]
    return node.kwargs["dataset"], steps[::-1]


def _how_loaded(dataset: Dataset) -> tuple:
    """What decides how a record is read: its schema, its loader and the loader's kwargs."""
    return (dataset.schema, Path(dataset.path).suffix, dataset.loader_kwargs)


def _run_name_of(job: Job, ledger_node: object) -> str:
    """The `run_name` of the record a ledger node scores, read off its load node."""
    windows_ref, _bench = ledger_node.inputs
    dataset, _ = _chain(job, windows_ref.node_id)
    return str(dataset.extra["run_name"])


# ------------------------------------------------- the ledger scores the panel's carve


def test_every_ledger_carve_is_the_t2star_jobs_carve_step_for_step():
    """Oracle: the T2* job's own graph. All 28 record carves must equal its chain.

    The T2* job is where the panel's windows come from, so its chain, read off its DAG, is
    the definition. A ledger job reading its record differently, or carving with any other
    function, captured filter config or kwarg, describes different windows from the panel.
    """
    reference = _load_job(JOBS / "active/t2star_q1_070423.py")
    (ref_windows,) = _nodes(reference, _windows_run)
    ref_dataset, ref_steps = _chain(reference, ref_windows.node_id)
    assert len(ref_steps) == 4, "filter, final stage, T2* fit, carve"

    compared = 0
    for job_id, job in _ledger_jobs().items():
        for node in _nodes(job, _windows_run):
            dataset, steps = _chain(job, node.node_id)
            assert _how_loaded(dataset) == _how_loaded(ref_dataset), (
                f"{job_id} {node.node_id} reads its record differently"
            )
            assert steps == ref_steps, f"{job_id} {node.node_id} carves differently"
            compared += 1
    assert compared == 28, f"compared {compared} carves, the family holds 28 records"


def test_every_ledger_is_scored_on_the_shared_knobs_and_its_carves_ladder():
    """Oracle: `LEDGER_KNOBS` and the windows node each ledger reads.

    The ladder a ledger labels its rows with must be the ladder its windows were carved at,
    or every column is mislabelled.
    """
    scored = 0
    for job_id, job in _ledger_jobs().items():
        for node in _nodes(job, ledger_step):
            windows_node, _bench = node.inputs
            carve = job.dag[windows_node.node_id]
            assert carve.fn is _windows_run
            assert node.kwargs["thresholds"] == carve.kwargs["thresholds"], job_id
            knobs = {k: v for k, v in node.kwargs.items() if k != "thresholds"}
            assert knobs == LEDGER_KNOBS, f"{job_id} {node.node_id}: {knobs}"
            scored += 1
    assert scored == 28


def test_the_km_job_scores_on_the_ledger_knobs_with_a_seed_of_its_own():
    """`km_with_checks_6d2s` scores records the ledger jobs also score.

    Its verdict knobs must equal `LEDGER_KNOBS`, or one record could read `pass` in one
    figure and `inconclusive` in the other. Its seed must differ: a shared permutation seed
    makes the two jobs' Monte Carlo error identical, so a p-value near alpha would land the
    same way in both and read as corroboration when it is one draw counted twice.
    """
    km = _load_module(JOBS / "active/km_with_checks_6d2s.py")
    for name in (
        "alpha",
        "min_events_pass",
        "tie_cutoff_distinct",
        "lag_max",
        "n_permutations",
    ):
        assert getattr(km, name.upper()) == LEDGER_KNOBS[name], name
    assert km.SEED != LEDGER_KNOBS["seed"]


def test_c3_is_off_and_its_simulation_count_is_declared():
    """Oracle: issue #19, the ledgers run the permutation rows only, so C3 is off. Its
    count stays declared, not inherited, for a caller that asks for C3.

    The C3 p-value floor is `1/(N+1)`, so a ledger at N = 200 reads a different object from
    one at N = 1000, and the number has to be declared rather than taken from a default.
    """
    assert LEDGER_KNOBS["include_c3"] is False
    assert LEDGER_KNOBS["c3_n_null_sim"] >= 200, (
        "below N = 200 the p-value floor rises above 0.005 and starts to crowd alpha"
    )


# --------------------------------------------------------- a ledger job runs checks only


def test_a_ledger_job_draws_and_writes_only_ledgers_and_their_survey():
    """No T2* panel and no windows artifact: those belong to the T2* jobs."""
    allowed = {CheckLedgerPanel, IndependenceSurveyOverviewPlot, *SURVEY_PLOTS}
    for job_id, job in _ledger_jobs().items():
        for sink in job.sinks:
            if hasattr(sink, "plot_class"):
                assert sink.plot_class in allowed, f"{job_id} draws {sink.plot_class}"
            else:
                assert sink.name.endswith(
                    (LEDGER_DATA_SUFFIX, "_check_ledgers", "_independence_survey_grids")
                ), f"{job_id} materializes {sink.name}"


def test_each_ledger_job_bundles_and_surveys_every_ledger_under_its_own_records_label():
    """Oracle: the record each ledger's chain loads, read off its load node.

    Row i of the survey must be the record ledger i scored. Labels are compared against the
    `run_name` of the dataset behind each ledger, not against node names.
    """
    for job_id, job in _ledger_jobs().items():
        ledgers = _nodes(job, ledger_step)
        (bundle,) = _nodes(job, collect_ledgers)
        assert [r.node_id for r in bundle.inputs] == [n.node_id for n in ledgers], (
            job_id
        )
        expected = tuple(record_label(_run_name_of(job, n)) for n in ledgers)
        assert bundle.kwargs["labels"] == expected, job_id
        (survey,) = _nodes(job, survey_ledgers)
        assert [r.node_id for r in survey.inputs] == [bundle.node_id], job_id
        assert survey.kwargs["dataset_labels"] == expected, job_id


def test_no_record_is_scored_twice():
    """One record, one ledger in the family.

    A `_before` / `_after` file is a row split of a full record, so scoring it beside the
    record scores the same reads twice; and a record in two jobs would be two Monte Carlo
    draws of one question, which a reader would take for corroboration.
    """
    paths = [
        str(node.kwargs["dataset"].path)
        for job in _ledger_jobs().values()
        for node in job.dag.values()
        if node.fn_name == "load"
    ]
    assert len(paths) == 28
    assert len(set(paths)) == len(paths)
    assert not [p for p in paths if p.endswith(("_before.pickle", "_after.pickle"))]


# ------------------------------------------------- the device survey re-scores nothing


def test_the_device_survey_stacks_every_ledger_of_the_family_and_computes_nothing_else():
    """Oracle: the family as `discover` finds it, and the bundle node of each job.

    One node, fed only by refs into included ledger jobs: no load, no carve, no ledger of
    its own. So every cell it draws is a verdict the record's own ledger printed. ONE ref
    per included job, because a composite run without reuse runs an included job once for
    every distinct ref.
    """
    survey_job = _load_job(JOBS / "composite/independence_survey.py")
    family = [j.job_id for j in by_family(JOBS, "check_ledger")]
    assert sorted(inc.job.name for inc in survey_job.includes) == sorted(family)

    assert list(survey_job.dag) == ["independence_survey"]
    (survey,) = survey_job.dag.values()
    assert survey.fn is survey_ledgers
    assert all(isinstance(ref, ArtifactRef) for ref in survey.inputs)

    expected_refs, expected_labels = [], []
    for inc in survey_job.includes:
        (bundle,) = _nodes(inc.job, collect_ledgers)
        (sink,) = [
            s
            for s in inc.job.sinks
            if getattr(s, "node", None) is not None and s.node.node_id == bundle.node_id
        ]
        expected_refs.append((inc.job.name, sink.name))
        expected_labels.extend(bundle.kwargs["labels"])
    assert [(r.included.job.name, r.node_name) for r in survey.inputs] == expected_refs
    assert len(survey.inputs) == len(survey_job.includes) == 6
    labels = survey.kwargs["dataset_labels"]
    assert labels == tuple(expected_labels)
    assert len(set(labels)) == len(labels) == 28
    assert survey.kwargs["alpha"] == LEDGER_KNOBS["alpha"]


def test_including_a_job_with_no_ledger_raises_rather_than_adding_no_rows():
    """A survey row silently missing reads as a record that was never there."""
    composite = Job("survey_probe")
    with pytest.raises(ValueError, match="not a ledger job"):
        included_ledger_set(composite.include("t2star_q1_070423"))


# ------------------------------------------------------- the labels name their records


def test_a_label_attached_to_another_records_ledger_raises():
    a = _fake_ledger("q1_040423", {"1 µs": VERDICT_PASS})
    b = _fake_ledger("q1_050423", {"1 µs": VERDICT_FAIL})
    assert collect_ledgers(a, b, labels=("q1 040423", "q1 050423")).labels == (
        "q1 040423",
        "q1 050423",
    )
    with pytest.raises(ValueError, match="attached to the ledger of"):
        collect_ledgers(a, b, labels=("q1 050423", "q1 040423"))


def test_declared_rows_out_of_step_with_the_bundles_raise():
    a = _fake_ledger("q1_040423", {"1 µs": VERDICT_PASS})
    b = _fake_ledger("q2_210423", {"1 µs": VERDICT_FAIL})
    sets = (
        collect_ledgers(a, labels=("q1 040423",)),
        collect_ledgers(b, labels=("q2 210423",)),
    )
    data = survey_ledgers(*sets, dataset_labels=("q1 040423", "q2 210423"), alpha=0.05)
    assert data.datasets == ["q1 040423", "q2 210423"]
    with pytest.raises(ValueError, match="differ from the ledgers' own"):
        survey_ledgers(*sets, dataset_labels=("q2 210423", "q1 040423"), alpha=0.05)


def test_a_repeated_label_raises_rather_than_drawing_one_record_twice():
    """Oracle: pivot semantics. The grid keeps the first row per label and drops the rest."""
    a = _fake_ledger("A", {"1 µs": VERDICT_PASS})
    b = _fake_ledger("B", {"1 µs": VERDICT_FAIL})
    with pytest.raises(ValueError, match="repeat"):
        build_independence_survey(a, b, dataset_labels=("A", "A"))


def test_c3_is_not_in_the_row_schema_so_it_cannot_be_bench_scored():
    """C3 is not in ROW_KEYS. Nothing may quietly start treating it as a scored check -
    `bench_acceptance_at_n` has no cell for it, so its verdicts carry no power evidence.
    """
    assert not any(key[0] == "c3_serial_copula" for key in ROW_KEYS)


# -------------------------------------------------------------- one figure per check


def test_there_is_exactly_one_figure_per_surveyed_instrument():
    """Adding a check to the battery must add a figure, not silently drop one: one figure
    per entry of SURVEY_KEYS, in its order."""
    from quebra.analyzers.independence_survey import SURVEY_KEYS

    assert len(SURVEY_PLOTS) == len(SURVEY_KEYS)
    assert [p.KEY for p in SURVEY_PLOTS] == list(SURVEY_KEYS)


def test_the_survey_draws_no_c3_grid():
    """Oracle: issue #19, "C3 leaves SURVEY_KEYS". Its code stays; no grid or figure
    draws it."""
    from quebra.analyzers.independence_survey import C3_KEY, SURVEY_KEYS

    assert C3_KEY not in SURVEY_KEYS
    assert C3_KEY not in ROW_KEYS
    assert not any(p.KEY == C3_KEY for p in SURVEY_PLOTS)


def test_every_surveyed_key_has_a_label_and_a_stated_null():
    """A grid whose caption cannot say what a red cell MEANS is not publishable."""
    from quebra.analyzers.independence_survey import SURVEY_KEYS

    for key in SURVEY_KEYS:
        assert key in CHECK_LABELS and CHECK_LABELS[key].strip()
        assert key in CHECK_NULL and "reject" in CHECK_NULL[key]


# --------------------------------------------------------------------- the reshape


def _fake_ledger(dataset_id: str, verdicts: dict[str, str]) -> object:
    rows = []
    for threshold, verdict in verdicts.items():
        for check, calibration, variant in ROW_KEYS:
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "threshold_label": threshold,
                    "clock": "in_spec",
                    "check_id": check,
                    "calibration": calibration,
                    "variant": variant,
                    "p_value": 0.01 if verdict == VERDICT_FAIL else 0.5,
                    "verdict": verdict,
                    "n_events": 100,
                }
            )

    class _L:
        pass

    ledger = _L()
    ledger.rows = pd.DataFrame(rows)
    ledger.dataset_id = dataset_id
    return ledger


def test_the_reshape_keeps_every_dataset_and_threshold():
    a = _fake_ledger("A", {"1 µs": VERDICT_PASS, "2 µs": VERDICT_FAIL})
    b = _fake_ledger("B", {"1 µs": VERDICT_FAIL, "2 µs": VERDICT_PASS})
    data = build_independence_survey(a, b, dataset_labels=("A", "B"))
    assert data.datasets == ["A", "B"]
    assert data.thresholds == ["1 µs", "2 µs"]
    assert len(data.grids) == len(ROW_KEYS)
    for grid in data.grids:
        assert grid.verdicts.shape == (2, 2)


def test_a_label_count_mismatch_raises_rather_than_guessing():
    a = _fake_ledger("A", {"1 µs": VERDICT_PASS})
    with pytest.raises(ValueError, match="will not guess"):
        build_independence_survey(a, dataset_labels=("A", "B"))


def test_a_non_ledger_input_raises():
    with pytest.raises(TypeError, match="not a CheckLedger"):
        build_independence_survey(object())


def test_threshold_order_follows_the_ladder_not_alphabetical_sort():
    """`10 µs` sorts before `2 µs` as a string; the ladder must not be drawn out of order."""
    a = _fake_ledger("A", {"2 µs": VERDICT_PASS, "10 µs": VERDICT_PASS})
    data = build_independence_survey(a, dataset_labels=("A",))
    assert data.thresholds == ["2 µs", "10 µs"]


def test_rejection_share_divides_by_DECIDED_cells_not_by_every_cell():
    """The claim CLAUDE.md's aggregation rule is about: state what the denominator is.

    A grid of 2 rejections and 38 inconclusive cells is "2 of 2 decided cells rejected", not
    "5% rejection". The second reads as reassurance and is not a rate of anything.
    """
    verdicts = pd.DataFrame(
        [[VERDICT_FAIL, VERDICT_INCONCLUSIVE], [VERDICT_INCONCLUSIVE, VERDICT_PASS]]
    )
    grid = InstrumentGrid(
        key=ROW_KEYS[0],
        label="x",
        null_statement="rejects when x",
        clock="in_spec",
        verdicts=verdicts,
        p_values=verdicts,
        n_events=verdicts,
    )
    assert grid.counts[VERDICT_INCONCLUSIVE] == 2
    assert grid.rejection_share_of_decided == pytest.approx(0.5)


def test_a_grid_with_nothing_decided_reports_nan_not_zero():
    """Zero would read as 'nothing rejected', which is not what 'nothing was decided' means."""
    verdicts = pd.DataFrame([[VERDICT_INCONCLUSIVE, VERDICT_INCONCLUSIVE]])
    grid = InstrumentGrid(
        key=ROW_KEYS[0],
        label="x",
        null_statement="rejects when x",
        clock="in_spec",
        verdicts=verdicts,
        p_values=verdicts,
        n_events=verdicts,
    )
    assert np.isnan(grid.rejection_share_of_decided)


def test_the_summary_has_one_row_per_grid_and_one_count_per_verdict():
    a = _fake_ledger("A", {"1 µs": VERDICT_PASS})
    data = build_independence_survey(a, dataset_labels=("A",))
    summary = survey_summary(data)
    assert len(summary) == len(data.grids)
    assert {"pass", "fail", "inconclusive", "not_interpretable", "not_computed"} <= set(
        summary.columns
    )


# ------------------------------------------------------ the trend group across ledgers

_UNDEFINED = "C1/C2 undefined here: tau {v} equals T_N {v}"


def _trend_cell(record, clock, c1, c2, cvm, value="3.5"):
    """One cell's permutation rows of the trend group; a check answers `pass` or not."""
    rows = []
    for check, answers in zip(TREND_CHECKS, (c1, c2, cvm), strict=True):
        rows.append(
            {
                "dataset_id": record,
                "threshold_label": "3 µs",
                "clock": clock,
                "check_id": check,
                "calibration": "permutation",
                "variant": "",
                "verdict": VERDICT_PASS if answers else VERDICT_NOT_COMPUTED,
                "notes": "" if answers else _UNDEFINED.format(v=value),
            }
        )
    return rows


def _trend_ledgers():
    """Three records, hand-built so every count below is read off the rows."""
    from quebra.analyzers.check_ledger import CheckLedger

    asymptotic_answer = {
        **_trend_cell("A", "in_spec", True, False, False)[0],
        "calibration": "asymptotic",
    }
    no_windows = {
        "dataset_id": "B",
        "threshold_label": "3 µs",
        "clock": "in_spec",
        "check_id": "(all)",
        "calibration": "",
        "variant": "",
        "verdict": VERDICT_NOT_COMPUTED,
        "notes": "no windows at this threshold",
    }
    a = _trend_cell("A", "in_spec", False, False, True)
    a += [asymptotic_answer] + _trend_cell("A", "calendar", True, True, True)
    b = [no_windows] + _trend_cell("B", "calendar", False, False, False, value="7")
    c = _trend_cell("C", "in_spec", False, True, True)
    return [
        CheckLedger(rows=pd.DataFrame(rows), dataset_id=record)
        for record, rows in (("A", a), ("B", b), ("C", c))
    ]


def test_the_trend_group_counts_equal_a_hand_table():
    """Oracle: the hand table of `_trend_ledgers`. A in-spec: C1 and C2 silent, CvM
    answers, and an asymptotic C1 answer that must not count. B in-spec has no windows;
    B calendar has nothing at all. C in-spec: C1 silent, C2 answers. Numbers in a message
    read `#`, the digits of `C1/C2` and `T_N` do not."""
    counts = trend_group_counts(_trend_ledgers())
    undefined = "C1/C2 undefined here: tau # equals T_N #"
    assert counts.no_answer.to_dict("records") == [
        {"clock": "calendar", "message": undefined, "n_cells": 1},
        {"clock": "in_spec", "message": undefined, "n_cells": 1},
        {"clock": "in_spec", "message": "no windows at this threshold", "n_cells": 1},
    ]
    per_clock = counts.cvm_alone.set_index("clock")
    assert per_clock.loc[
        "calendar", ["n_cells", "n_trend_answered", "n_cvm_alone"]
    ].tolist() == [2, 1, 0]
    assert per_clock.loc[
        "in_spec", ["n_cells", "n_trend_answered", "n_cvm_alone"]
    ].tolist() == [3, 2, 1]
    assert per_clock.loc["calendar", "share_of_all"] == 0.0
    assert per_clock.loc["in_spec", "share_of_answered"] == 0.5
    assert per_clock.loc["in_spec", "share_of_all"] == pytest.approx(1.0 / 3.0)
    with pytest.raises(ValueError, match="repeat"):
        trend_group_counts(_trend_ledgers() + _trend_ledgers()[:1])


def test_the_trend_count_refuses_a_cell_where_c1_and_c2_were_never_asked():
    """Oracle: the count's contract. A cell with CvM's row but no row of C1's or C2's and
    no `(all)` row never asked them, which is not the same as no answer; a ledger with no
    cell has nothing to count."""
    from quebra.analyzers.check_ledger import CheckLedger

    cvm_only = _trend_cell("D", "in_spec", False, False, True)[2:]
    with pytest.raises(ValueError, match="never asked"):
        trend_group_counts([CheckLedger(rows=pd.DataFrame(cvm_only), dataset_id="D")])
    no_rows = pd.DataFrame(columns=list(cvm_only[0]))
    with pytest.raises(ValueError, match="no cell"):
        trend_group_counts([CheckLedger(rows=no_rows, dataset_id="E")])


def test_one_reason_is_one_message_however_the_ledger_prints_its_numbers():
    """Oracle: the ledger itself. One carve at five time scales fails C1 and C2 on the
    in-spec clock for one reason, `tau == T_N`; the ledger prints its numbers at different
    widths, cuts the note and writes it twice, and the count still reads one message."""
    from quebra.analyzers import check_ledger, windows

    bench = pd.read_csv(
        Path(__file__).resolve().parents[1]
        / "jobs"
        / "bench"
        / "results"
        / "size_table.csv"
    )
    rng = np.random.default_rng(3)
    values: list[float] = []
    for _ in range(60):
        values += [5.0] * int(rng.integers(2, 12))
        values += [1.0] * int(rng.integers(2, 8))
    series = np.asarray(values, dtype=float)
    ladder = [("3 µs", 3.0, True)]
    ledgers = []
    for i, scale in enumerate((1.0, 7.3, 61.7, 1234.5, 98765.4)):
        carved = windows.run(
            windows.WindowsInputs(
                t_rel_s=np.arange(len(series), dtype=float) * 60.0 * scale,
                values=series,
                thresholds=ladder,
                dataset_id=f"s{i}",
            )
        )
        ledgers.append(
            check_ledger.run(
                check_ledger.make_inputs_from_windows(
                    carved,
                    bench,
                    thresholds=ladder,
                    alpha=0.05,
                    min_events_pass=35,
                    tie_cutoff_distinct=5,
                    lag_max=5,
                    n_permutations=19,
                    seed=7,
                )
            )
        )
    raw = set()
    for ledger in ledgers:
        rows = ledger.rows
        c1 = rows[
            (rows["clock"] == "in_spec")
            & (rows["check_id"] == "c1_lewis_robinson")
            & (rows["calibration"] == "permutation")
        ]
        raw.add(str(c1["notes"].iloc[0]))
    assert len(raw) == len(ledgers), "premise: the raw notes differ in their numbers"
    counts = trend_group_counts(ledgers)
    in_spec = counts.no_answer[counts.no_answer["clock"] == "in_spec"]
    assert in_spec["n_cells"].tolist() == [len(ledgers)], in_spec.to_dict("records")
    assert (
        in_spec["message"]
        .iloc[0]
        .startswith("C1/C2 undefined here: tau (#) does not clear T_N (#)")
    )


def test_no_survey_figure_names_c3():
    """Oracle: issue #19, "C3 leaves SURVEY_KEYS". Every survey figure, each instrument's
    and the overview, is rendered, and none of its text names C3."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.text import Text

    a = _fake_ledger("A", {"2 µs": VERDICT_PASS, "3 µs": VERDICT_FAIL})
    b = _fake_ledger("B", {"2 µs": VERDICT_FAIL, "3 µs": VERDICT_PASS})
    data = build_independence_survey(a, b, dataset_labels=("A", "B"))
    for plot in (*SURVEY_PLOTS, IndependenceSurveyOverviewPlot):
        fig = plot(name="c3_probe").build_matplotlib(data)
        try:
            texts = [t.get_text() for t in fig.findobj(Text) if t.get_text()]
            assert texts, f"control: {plot.__name__} draws text"
            named = [t for t in texts if "C3" in t]
            assert not named, (plot.__name__, named)
        finally:
            plt.close(fig)


def _count_script():
    spec = importlib.util.spec_from_file_location(
        "count_trend_answers",
        Path(__file__).resolve().parents[1] / "scripts" / "count_trend_answers.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_count_script_reads_every_ledger_in_a_folder(tmp_path, capsys):
    """Oracle: the hand table again, now through the pickles on disk the script reads."""
    import pickle

    for ledger in _trend_ledgers():
        path = tmp_path / f"q1_{ledger.dataset_id}_check_ledger_data.pkl"
        path.write_bytes(pickle.dumps(ledger))
    assert _count_script().main([str(tmp_path)]) == 0
    printed = capsys.readouterr().out
    assert "3 ledgers from 1 folders" in printed
    assert "no windows at this threshold" in printed
    assert "0.5000" in printed


def test_the_count_script_refuses_a_stale_ledger_and_prints_no_count(tmp_path, capsys):
    """Oracle: `load_artifact`'s refusal. A ledger pickled without one of its fields is
    named, the script prints no count and exits 2: a count over some of the ledgers would
    read as a count over all of them."""
    import pickle

    ledgers = _trend_ledgers()
    del ledgers[0].__dict__["thresholds"]
    for ledger in ledgers:
        path = tmp_path / f"q1_{ledger.dataset_id}_check_ledger_data.pkl"
        path.write_bytes(pickle.dumps(ledger))
    assert _count_script().main([str(tmp_path)]) == 2
    printed = capsys.readouterr().out
    assert "refused 1 of 3 ledgers as stale" in printed
    assert "trend group" not in printed
