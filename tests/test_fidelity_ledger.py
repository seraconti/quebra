"""The check ledgers of the 2x2 infidelity windows: their records, carve, failure path, figure
and survey.

The job graphs are built from the job files and read no record. The ledgers that run do so on
synthetic |x| carves.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

import quebra.ledger_recipe as ledger_recipe  # noqa: E402
from quebra.analyzers import fidelity  # noqa: E402
from quebra.analyzers.check_ledger import VERDICT_NOT_COMPUTED, CheckLedger  # noqa: E402
from quebra.analyzers.windows import SIDE_IN_SPEC  # noqa: E402
from quebra.cli import _module_from_path  # noqa: E402
from quebra.core.closure import code_closure  # noqa: E402
from quebra.core.discovery import by_family  # noqa: E402
from quebra.core.job import Job  # noqa: E402
from quebra.core.reference import ArtifactRef  # noqa: E402
from quebra.core.types import RecordFailure  # noqa: E402
from quebra.fidelity_ledger_recipe import (  # noqa: E402
    BENCH_TRANSFER,
    LEDGER_STAGE,
    InfidelityCheckLedger,
    fidelity_ledger_step,
)
from quebra.ledger_recipe import (  # noqa: E402
    LEDGER_KNOBS,
    collect_ledgers,
    ledger_bag,
    record_label,
    survey_ledgers,
)
from quebra.panels.infidelity_check_ledger import InfidelityCheckLedgerPanel  # noqa: E402
from quebra.recipes import _fidelity_windows  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
JOBS = REPO / "jobs"
FAMILY = "check_ledger_fidelity"
LADDER = list(fidelity.DETUNING_LADDER)
SIZE_TABLE = REPO / "jobs/bench/results/size_table.csv"


def _job(path: Path) -> Job:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return _module_from_path(path).job


def _family(family: str) -> dict[str, Job]:
    return {d.job_id: _job(d.path) for d in by_family(JOBS, family)}


def _nodes(job: Job, fn: object) -> list:
    return [node for node in job.dag.values() if node.fn is fn]


def _materialized(job: Job) -> dict[str, str]:
    """Each materialize sink's name, mapped to the node it persists."""
    return {
        sink.name: sink.node.node_id
        for sink in job.sinks
        if getattr(sink, "node", None) is not None
    }


def _carve(record: str, seed: int) -> object:
    """A synthetic |x| record crossing every rung of the nines ladder, carved as the
    infidelity carve is (`recipes._fidelity_windows`)."""
    rng = np.random.default_rng(seed)
    n = 400
    frame = pd.DataFrame(
        {
            "t_rel_s": np.arange(n, dtype=float) * 90.0,
            "abs_x": np.abs(
                0.06 + 0.08 * np.sin(np.arange(n) / 7.0) + rng.normal(0, 0.03, n)
            ),
            "sigma_x": np.full(n, 2e-3),
        }
    )
    result = fidelity.FidelityResult(frame=frame, meta={"dataset_id": record})
    return _fidelity_windows(result, 10.0, 1.0, True, LADDER)


@pytest.mark.integration
def test_the_infidelity_ledgers_score_the_t2star_2x2_ledgers_records_in_order(
    in_repo,
) -> None:
    """Oracle: the T2* 2x2 ledger jobs' graphs. Qubit by qubit, each infidelity ledger job
    includes the Ramsey job with fidelity of exactly the T2* job's records, in its order;
    each record's ledger carries the T2* run name as its name in the paper's tables, and is
    labelled by its own Ramsey job's run name."""
    infidelity = _family(FAMILY)
    assert sorted(infidelity) == [
        f"check_ledger_fidelity_2x2_q{q}" for q in range(1, 5)
    ]
    for q in range(1, 5):
        (t2star_bundle,) = _nodes(
            _job(JOBS / f"active/check_ledger_2x2_q{q}.py"), collect_ledgers
        )
        runs = [label.replace(" ", "_") for label in t2star_bundle.kwargs["labels"]]
        job = infidelity[f"check_ledger_fidelity_2x2_q{q}"]
        assert [inc.job.name for inc in job.includes] == [f"fid_{run}" for run in runs]
        ledgers = _nodes(job, fidelity_ledger_step)
        assert [node.kwargs["paper_name"] for node in ledgers] == runs
        (bundle,) = _nodes(job, collect_ledgers)
        assert bundle.kwargs["labels"] == tuple(
            record_label(f"fid_{run}") for run in runs
        )


@pytest.mark.integration
def test_every_ledger_reads_its_ramsey_jobs_own_carve_with_the_t2star_knobs(
    in_repo,
) -> None:
    """Oracle: each included Ramsey job's own graph, and the T2* ledgers' knobs. Each ledger
    reads `_fidelity_windows` on one ref, to its Ramsey job's persisted fidelity result,
    with that job's own carve kwargs (the nines ladder among them); it is scored on that
    ladder with `LEDGER_KNOBS` and the bench size table; the carve and the ledger are each
    persisted under the record's run name, and each bag gets the infidelity figure."""
    n_records = 0
    for job in _family(FAMILY).values():
        sinks = _materialized(job)
        figures = [s for s in job.sinks if getattr(s, "node", None) is None]
        for ledger in _nodes(job, fidelity_ledger_step):
            carve_ref, bench_ref = ledger.inputs
            carve = job.dag[carve_ref.node_id]
            assert carve.fn is _fidelity_windows
            (source,) = carve.inputs
            assert isinstance(source, ArtifactRef)
            ramsey = source.included.job
            record = ramsey.name
            assert source.node_name == f"{record}_fidelity_raw"
            (own,) = _nodes(ramsey, _fidelity_windows)
            assert carve.kwargs == own.kwargs
            assert own.kwargs["thresholds"] == LADDER
            assert ledger.kwargs["thresholds"] == LADDER
            assert {k: ledger.kwargs[k] for k in LEDGER_KNOBS} == LEDGER_KNOBS
            bench = job.dag[bench_ref.node_id]
            assert Path(bench.kwargs["dataset"].path) == Path(
                "jobs/bench/results/size_table.csv"
            )
            assert sinks[f"{record}_windows"] == carve.node_id
            assert sinks[f"{record}_check_ledger_data"] == ledger.node_id
            drawn = [
                f
                for f in figures
                if job.dag[f.input.node_id].inputs[0].node_id == ledger.node_id
            ]
            assert [f.plot_class for f in drawn] == [InfidelityCheckLedgerPanel] * 2
            n_records += 1
    assert n_records == 22


@pytest.mark.unit
def test_a_ledger_that_refuses_its_record_reads_not_computed_and_the_survey_completes(
    monkeypatch,
) -> None:
    """Oracle: specification, the T2* ledgers' failure contract. A ledger that raises
    ValueError or KeyError for one record becomes that record's `not computed` ledger, with
    the stage and the reason at every threshold, bag and clock; the other record is scored,
    and the survey over all three still builds."""
    real = ledger_recipe.make_inputs_from_windows
    refusals = {"fid_refused": ValueError, "fid_refused_key": KeyError}

    def refusing(window_result, *args, **kwargs):
        error = refusals.get(window_result.meta["dataset_id"])
        if error is not None:
            raise error("synthetic refusal")
        return real(window_result, *args, **kwargs)

    monkeypatch.setattr(ledger_recipe, "make_inputs_from_windows", refusing)
    bench = pd.read_csv(SIZE_TABLE)
    knobs = dict(LEDGER_KNOBS, n_permutations=19)
    sound, refused, refused_key = (
        fidelity_ledger_step(
            _carve(record, seed),
            bench,
            paper_name=paper,
            thresholds=LADDER,
            **knobs,
        )
        for record, seed, paper in (
            ("fid_sound", 7, "2x2_q8_010101"),
            ("fid_refused", 8, "2x2_q9_010101"),
            ("fid_refused_key", 9, "2x2_q7_010101"),
        )
    )
    assert isinstance(refused, InfidelityCheckLedger)
    assert refused.dataset_id == "fid_refused"
    assert (refused.paper_name, refused.bench_transfer) == (
        "2x2_q9_010101",
        BENCH_TRANSFER,
    )
    assert set(refused.rows["verdict"]) == {VERDICT_NOT_COMPUTED}
    assert len(refused.rows) == len(LADDER) * 2 * 2
    assert all(
        f"{LEDGER_STAGE} failed: ValueError: synthetic refusal" in str(note)
        for note in refused.rows["notes"]
    )
    assert refused_key.dataset_id == "fid_refused_key"
    assert set(refused_key.rows["verdict"]) == {VERDICT_NOT_COMPUTED}
    assert all(
        f"{LEDGER_STAGE} failed: KeyError: 'synthetic refusal'" in str(note)
        for note in refused_key.rows["notes"]
    )
    assert not sound.rows["notes"].astype(str).str.contains(LEDGER_STAGE).any()
    labels = tuple(
        record_label(r) for r in ("fid_sound", "fid_refused", "fid_refused_key")
    )
    survey = survey_ledgers(
        collect_ledgers(sound, refused, refused_key, labels=labels),
        dataset_labels=labels,
        alpha=LEDGER_KNOBS["alpha"],
    )
    assert list(survey.datasets) == list(labels)


@pytest.mark.unit
def test_the_figure_says_the_paper_name_and_that_the_bench_is_not_measured_on_x() -> (
    None
):
    """Oracle: specification. Each bag's figure prints the record's name in the paper's
    tables and the bench statement the ledger carries, which says the size cells are not
    measured on an |x| carve and the ties cutoff is not calibrated there; a T2* ledger,
    which carries neither, is refused."""
    for phrase in (
        "size cells",
        "|x| carve",
        "not measured",
        "ties cutoff",
        "not calibrated",
    ):
        assert phrase in BENCH_TRANSFER, phrase
    bench = pd.read_csv(SIZE_TABLE)
    ledger = fidelity_ledger_step(
        _carve("fid_r", 7),
        bench,
        paper_name="2x2_q9_010101",
        thresholds=LADDER,
        **dict(LEDGER_KNOBS, n_permutations=19),
    )
    base = ledger_recipe.ledger_step(
        RecordFailure(dataset_id="fid_r", stage="s", reason="r"),
        bench,
        thresholds=LADDER,
        **LEDGER_KNOBS,
    )
    assert type(base) is CheckLedger
    fig = InfidelityCheckLedgerPanel(name="unit").build_matplotlib(
        ledger_bag(ledger, side=SIDE_IN_SPEC)
    )
    try:
        text = " ".join(" ".join(t.get_text().split()) for t in fig.texts)
        assert "2x2_q9_010101 in the paper's dataset tables" in text, text
        assert " ".join(BENCH_TRANSFER.split()) in text, text
    finally:
        plt.close(fig)
    with pytest.raises(TypeError, match="expects an InfidelityCheckLedger"):
        InfidelityCheckLedgerPanel(name="unit").build_matplotlib(
            ledger_bag(base, side=SIDE_IN_SPEC)
        )


@pytest.mark.unit
def test_the_t2star_ledger_jobs_import_no_infidelity_ledger_module() -> None:
    """Oracle: the identity closure itself. A guard on future imports: none of the T2*
    ledger jobs' step modules may reach the infidelity ledger's modules, or a change there
    would re-key them. It is not the evidence that their keys are unchanged; that is the
    job_code_hash measurement at the checkpoint."""
    new = {
        "quebra.fidelity_ledger_recipe",
        "quebra.panels.infidelity_check_ledger",
        "quebra.infidelity_recipe",
    }
    jobs = _family("check_ledger")
    assert len(jobs) == 10
    for name, job in jobs.items():
        closure = set(code_closure([node.fn for node in job.dag.values()]))
        assert "quebra.ledger_recipe" in closure, f"control: {name} reaches its recipe"
        assert not closure & new, (name, sorted(closure & new))


@pytest.mark.integration
def test_the_survey_reads_one_bundle_per_infidelity_ledger_job_rows_in_job_order(
    in_repo,
) -> None:
    """Oracle: the four ledger jobs' graphs. The survey includes exactly the family, takes
    one ref per job, to the bundle that job persists, and its rows are every job's labels,
    in qubit order: the 22 records of the 2x2."""
    survey_job = _job(JOBS / "composite/independence_survey_fidelity_2x2.py")
    family = sorted(_family(FAMILY))
    assert [inc.job.name for inc in survey_job.includes] == family
    (survey,) = survey_job.dag.values()
    assert survey.fn is survey_ledgers
    expected: list[str] = []
    for ref, job_id in zip(survey.inputs, family, strict=True):
        assert isinstance(ref, ArtifactRef) and ref.included.job.name == job_id
        assert ref.node_name == f"fid_2x2_q{job_id[-1]}_check_ledgers"
        (bundle,) = _nodes(ref.included.job, collect_ledgers)
        expected.extend(bundle.kwargs["labels"])
    assert list(survey.kwargs["dataset_labels"]) == expected
    assert len(expected) == 22
