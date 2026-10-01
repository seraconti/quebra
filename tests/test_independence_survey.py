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
    VERDICT_PASS,
    VERDICT_UNDERPOWERED,
)
from quebra.analyzers.checks.battery import ROW_KEYS
from quebra.analyzers.independence_survey import (
    CHECK_LABELS,
    CHECK_NULL,
    InstrumentGrid,
    build_independence_survey,
    survey_summary,
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
    figure and `underpowered` in the other. Its seed must differ: a shared permutation seed
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


def test_c3_is_on_and_its_simulation_count_is_declared():
    """C3 is the one instrument with no bench cell. The count is declared, not inherited.

    The C3 p-value floor is `1/(N+1)`, so a ledger at N = 200 reads a different object from
    one at N = 1000, and the number has to be declared rather than taken from a default.
    """
    assert LEDGER_KNOBS["include_c3"] is True
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
    """The consequence of C3 being on: it is in the figures but not in ROW_KEYS.

    Nothing may quietly start treating it as a scored check - `bench_acceptance_at_n` has
    no cell for it, so its verdicts carry no power evidence.
    """
    assert not any(key[0] == "c3_serial_copula" for key in ROW_KEYS)


# -------------------------------------------------------------- one figure per check


def test_there_is_exactly_one_figure_per_surveyed_instrument():
    """Adding a check to the battery must add a figure, not silently drop one.

    The list is SURVEY_KEYS, not ROW_KEYS: C3 is surveyed but deliberately not bench-scored,
    and it was RUN but drawn nowhere until this was separated - the ledger computed C3 rows
    for every dataset and every grid was built from ROW_KEYS, so they went in the bin.
    """
    from quebra.analyzers.independence_survey import SURVEY_KEYS

    assert len(SURVEY_PLOTS) == len(SURVEY_KEYS)
    assert [p.KEY for p in SURVEY_PLOTS] == list(SURVEY_KEYS)


def test_the_survey_draws_c3_even_though_the_bench_cannot_score_it():
    from quebra.analyzers.independence_survey import C3_KEY, SURVEY_KEYS

    assert C3_KEY in SURVEY_KEYS
    assert C3_KEY not in ROW_KEYS
    assert any(p.KEY == C3_KEY for p in SURVEY_PLOTS)


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

    A grid of 2 rejections and 38 powerless cells is "2 of 2 decided cells rejected", not
    "5% rejection". The second reads as reassurance and is not a rate of anything.
    """
    verdicts = pd.DataFrame(
        [[VERDICT_FAIL, VERDICT_UNDERPOWERED], [VERDICT_UNDERPOWERED, VERDICT_PASS]]
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
    assert grid.counts[VERDICT_UNDERPOWERED] == 2
    assert grid.rejection_share_of_decided == pytest.approx(0.5)


def test_a_grid_with_nothing_decided_reports_nan_not_zero():
    """Zero would read as 'nothing rejected', which is not what 'nothing was decided' means."""
    verdicts = pd.DataFrame([[VERDICT_UNDERPOWERED, VERDICT_UNDERPOWERED]])
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


def test_the_summary_has_one_row_per_grid():
    a = _fake_ledger("A", {"1 µs": VERDICT_PASS})
    data = build_independence_survey(a, dataset_labels=("A",))
    assert len(survey_summary(data)) == len(data.grids)
