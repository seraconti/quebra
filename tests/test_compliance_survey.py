"""The compliance survey: the in-spec fraction of every ledger record, with its bounds.

Oracles:
- the T2* tellability step (`tellability_recipe._t2star_tellability`): on the same
  synthetic record, the survey's bounds are its bounds, bit for bit;
- the T2* panel (`recipes._t2star_panel_data`): the survey's fraction is the panel's
  occupancy, bit for bit, and lies between the bounds;
- a hand-built record whose fraction and bounds are hand-computed;
- the job graphs: the composite surveys every record of every ledger job, on the ladder,
  through one carve bundle per job holding the carves the ledgers scored;
- the paper's dataset tables, Supplementary Figs. 10 and 11, transcribed below;
- the artifact: the figure draws its fractions and bounds, rung by rung, in ladder order.

Synthetic data only. Importing a job builds its graph and reads no data.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import matplotlib
import numpy as np
import pytest

from quebra.analyzers import windows
from quebra.analyzers.compliance_survey import (
    PanelCounts,
    ComplianceSurveyData,
    build_compliance_survey,
    record_compliance,
)
from quebra.analyzers.kaplan_meier import CONF_LEVEL
from quebra.analyzers.tellability import COMPUTED
from quebra.core.dataset import Dataset
from quebra.core.discovery import by_family
from quebra.core.job import Job, _load_dataset
from quebra.core.reference import ArtifactRef
from quebra.core.types import RecordFailure
from quebra.ledger_recipe import (
    CARVES_DATA_SUFFIX,
    RecordCarves,
    collect_carves,
    included_carves,
    record_label,
)
from quebra.recipes import (
    _CARVE_STAGES,
    T2STAR_GAP_MULT,
    T2STAR_K,
    T2STAR_THRESHOLDS,
    T2STAR_USE_UNCERTAINTY,
    XI_SEED,
    _event_tables,
    _kaplan_meier_set,
    _t2star_panel_data,
    _windows_run,
    wire_t2star_carve,
)
from quebra.schemas.track912 import track912Schema
from quebra.tellability_recipe import _t2star_tellability

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
JOBS = REPO / "jobs"

_COLUMNS = [
    "timestamp", "amp", "amp error", "offset", "offset error", "T2star", "T2star error",
    "frequency", "frequency error", "phase", "phase_error", "normalised chi-square",
]  # fmt: skip


def _track912_csv(path: Path, *, seed: int, gap: bool, n: int = 260) -> Path:
    """A synthetic record whose T2* wanders across the 2 to 5 us ladder, with errors wide
    enough that some reads are not sure at k = 1 and some not at k = 2; `gap` adds a long
    gap and some failed fits."""
    rng = np.random.default_rng(seed)
    t = np.arange(n) * 30.0 + rng.uniform(-4.0, 4.0, n)
    if gap:
        t[n // 2 :] += 6.0 * 3600.0
    t2star = 3.5e-6 + 1.3e-6 * np.sin(np.arange(n) / 11.0) + rng.normal(0, 2e-7, n)
    sigma = rng.uniform(0.5e-7, 4e-7, n)
    lines = [",".join(_COLUMNS)]
    for i in range(n):
        value = "nan" if gap and i % 37 == 5 else repr(float(t2star[i]))
        row = [repr(float(t[i])), "1.0", "0.02", "0.5", "0.005", value]
        row += [repr(float(sigma[i])), "2.1e6", "6e3", "0.1", "0.03", "0.1"]
        lines.append(",".join(row))
    path.write_text("\n".join(lines) + "\n")
    return path


def _carve(path: Path, name: str) -> tuple[object, windows.WindowsResult]:
    """The ledger jobs' carve, stage by stage, on one synthetic record: (T2* fit, carve)."""
    node: object = _load_dataset(
        Dataset(
            path=str(path),
            schema=track912Schema,
            qubit=1,
            device="probe",
            extra={"run_name": name, "run_start_unix_s": 1.68e9},
        )
    )
    for stage in ("t2star_filter", "t2star_final_filter_stage", "t2star"):
        node = _CARVE_STAGES[stage](node)
    carve = _windows_run(
        node,
        gap_mult=T2STAR_GAP_MULT,
        k=T2STAR_K,
        use_uncertainty=T2STAR_USE_UNCERTAINTY,
        thresholds=T2STAR_THRESHOLDS,
    )
    return node, carve


@pytest.fixture(scope="module")
def carves(tmp_path_factory) -> list[tuple[object, windows.WindowsResult]]:
    root = tmp_path_factory.mktemp("compliance")
    return [
        _carve(_track912_csv(root / f"0{i}0423_r.csv", seed=i, gap=gap), f"q1_0{i}0423")
        for i, gap in ((1, False), (2, True), (3, True))
    ]


def _labels() -> list[str]:
    return [label for label, _, _ in T2STAR_THRESHOLDS]


# ----------------------------------------------- S1: the bounds are tellability's bounds


@pytest.mark.statistical
def test_the_bounds_are_the_t2star_tellability_bounds_bit_for_bit(carves):
    """Oracle: the T2* tellability step on the same record. It reads the T2* fit's frame;
    the survey reads only the carve's read table and settings. Every rung's bounds and
    their status must agree exactly."""
    unsettled = 0
    for result, carve in carves:
        told = _t2star_tellability(
            result,
            carve,
            thresholds=T2STAR_THRESHOLDS,
            gap_mult=T2STAR_GAP_MULT,
            unsure_k=T2STAR_K,
            n_permutations=9,
            seed=1,
        )
        by_label = {row.threshold_label: row for row in told.thresholds}
        survey = record_compliance(carve, T2STAR_THRESHOLDS)
        assert [r.threshold_label for r in survey.rungs] == _labels()
        for rung in survey.rungs:
            row = by_label[rung.threshold_label]
            assert rung.bounds_status == row.status["bounds"]
            assert rung.in_spec_fraction_lower == row.in_spec_fraction_lower
            assert rung.in_spec_fraction_upper == row.in_spec_fraction_upper
            if rung.bounds_status == COMPUTED:
                unsettled += rung.in_spec_fraction_upper > rung.in_spec_fraction_lower
    assert unsettled >= 10, "control: the records must leave intervals unsettled"


# ------------------------------------- S2: the T2* tellability reads its carve's settings


@pytest.mark.policy
@pytest.mark.parametrize("stem", ["t2star_q1_070423", "t2star_q1_100423"])
def test_the_t2star_tellability_reads_its_carves_k_gap_rule_and_ladder(stem):
    """Oracle: the job graph. The survey takes k and the gap rule from the carve's meta,
    and the T2* tellability node must take the same from its job's carve node. With S1
    and `test_every_ledger_carve_is_the_t2star_jobs_carve_step_for_step` (every ledger
    carve is the T2* job's carve), the survey's bounds on the two T2* records are the
    bounds their tellability reports, with no real data loaded."""
    spec = importlib.util.spec_from_file_location(
        f"_compliance_{stem}", JOBS / "active" / f"{stem}.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    dag = module.job.dag
    told, carve = dag["tellability_t2star"], dag["windows"]
    assert any(getattr(ref, "node_id", None) == "windows" for ref in told.inputs)
    assert told.kwargs["unsure_k"] == carve.kwargs["k"]
    assert told.kwargs["gap_mult"] == carve.kwargs["gap_mult"]
    assert list(told.kwargs["thresholds"]) == list(carve.kwargs["thresholds"])


# ------------------------------- S3: the fraction is the panel's occupancy, in its bounds


@pytest.mark.statistical
def test_the_fraction_is_the_panels_occupancy_and_lies_between_its_bounds(carves):
    """Oracle: the T2* panel on the same carve. Its occupancy is the one occupancy
    definition; the survey's fraction must be that number, and between the bounds."""
    for result, carve in carves:
        km = _kaplan_meier_set(_event_tables(carve), CONF_LEVEL)
        panel = _t2star_panel_data(
            result,
            carve,
            km,
            shape_min_reads=5,
            use_uncertainty=T2STAR_USE_UNCERTAINTY,
            xi_seed=XI_SEED,
            thresholds=T2STAR_THRESHOLDS,
        )
        occupancy = panel.reliability.occupancy
        for rung in record_compliance(carve, T2STAR_THRESHOLDS).rungs:
            assert rung.status == COMPUTED
            assert rung.in_spec_fraction == occupancy[rung.threshold_label]
            assert (
                rung.in_spec_fraction_lower
                <= rung.in_spec_fraction
                <= rung.in_spec_fraction_upper
            )


def _hand_carve(
    sigma: np.ndarray,
    *,
    k: float = 1.0,
    gap_mult: float = 10.0,
    dataset_id: str = "hand",
    use_uncertainty: bool = True,
) -> windows.WindowsResult:
    # Spacings 100, 200, 100, 100, then a 1200 s gap (the median is 100 s, the gap rule
    # 10 x 100), then 50, 100. At threshold 3, reads in: T T F F T T T F; read 1 is not
    # sure at k = 1 (|margin| 1 < sigma 2).
    t_s = np.array([0.0, 100.0, 300.0, 400.0, 500.0, 1700.0, 1750.0, 1850.0])
    values = np.array([4.0, 4.0, 2.0, 2.0, 4.0, 4.0, 4.0, 2.0])
    return windows.run(
        windows.WindowsInputs(
            t_rel_s=t_s,
            values=values,
            thresholds=[("3", 3.0, True)],
            sigma=sigma,
            dataset_id=dataset_id,
            gap_mult=gap_mult,
            k=k,
            use_uncertainty=use_uncertainty,
        )
    )


_NO_SIGMA = np.full(8, np.nan)
_SIGMA = np.array([0.1, 2.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])


@pytest.mark.unit
def test_a_hand_built_record_gives_the_hand_computed_fractions():
    """Oracle: a hand table. Observed time is 650 s, the 1200 s gap excluded. Intervals
    held by an in-spec read: 100 + 200 + 50 + 100 = 450, so 9/13. Settled in spec: [0,
    100) and [1700, 1750), 150, so 3/13. Unsettled: the three crossings and nothing else,
    200 + 100 + 100 = 400, so the upper bound is 550/650 = 11/13. Read 1's unsure state
    changes nothing, since its interval holds a crossing; the read-count share, 5/8,
    differs from the time-weighted 9/13."""
    sigma = np.array([0.1, 2.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    (rung,) = record_compliance(_hand_carve(sigma), [("3", 3.0, True)]).rungs
    assert rung.in_spec_fraction == pytest.approx(9 / 13, rel=1e-12)
    assert rung.in_spec_fraction_lower == pytest.approx(3 / 13, rel=1e-12)
    assert rung.in_spec_fraction_upper == pytest.approx(11 / 13, rel=1e-12)

    # One unsure read inside a run: read 0 at margin 1 with sigma 2 unsettles [0, 100).
    sigma[0] = 2.0
    (rung,) = record_compliance(_hand_carve(sigma), [("3", 3.0, True)]).rungs
    assert rung.in_spec_fraction == pytest.approx(9 / 13, rel=1e-12)
    assert rung.in_spec_fraction_lower == pytest.approx(1 / 13, rel=1e-12)
    assert rung.in_spec_fraction_upper == pytest.approx(11 / 13, rel=1e-12)


@pytest.mark.unit
def test_a_rung_or_record_that_cannot_answer_reads_not_computed_never_zero():
    """Oracle: the reasons' contract. No finite sigma leaves the fraction and drops the
    bounds; one read leaves neither; a recorded failure reads its note at every rung."""
    (rung,) = record_compliance(_hand_carve(_NO_SIGMA), [("3", 3.0, True)]).rungs
    assert rung.status == COMPUTED
    assert rung.bounds_status == "not computed: no read with a finite sigma"
    assert rung.in_spec_fraction_lower is None is rung.in_spec_fraction_upper

    one = windows.run(
        windows.WindowsInputs(
            t_rel_s=np.array([0.0]),
            values=np.array([4.0]),
            thresholds=[("3", 3.0, True)],
            sigma=np.array([0.1]),
            dataset_id="one",
            k=1.0,
            use_uncertainty=True,
        )
    )
    (rung,) = record_compliance(one, [("3", 3.0, True)]).rungs
    reason = "not computed: no segment with at least 2 reads"
    assert (rung.status, rung.bounds_status) == (reason, reason)
    assert rung.in_spec_fraction is None

    failed = record_compliance(
        RecordFailure(dataset_id="f", stage="load", reason="E: x"), T2STAR_THRESHOLDS
    )
    assert failed.record_failure == "load failed: E: x"
    assert [r.threshold_label for r in failed.rungs] == _labels()
    assert {r.status for r in failed.rungs} == {"not computed: load failed: E: x"}
    assert {r.in_spec_fraction for r in failed.rungs} == {None}


@pytest.mark.unit
def test_a_fraction_outside_its_bounds_or_a_read_table_out_of_step_raises(monkeypatch):
    """Oracle: the bounds' contract. A fraction that leaves its bounds, or a read table
    that is not the carve's reads, is a defect and raises; neither ships a number."""
    import quebra.analyzers.compliance_survey as compliance

    carve = _hand_carve(np.array([0.1, 2.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]))
    short = windows.WindowsResult(
        windows_in_spec=carve.windows_in_spec,
        windows_out_of_spec=carve.windows_out_of_spec,
        reads=carve.reads.iloc[1:],
        meta=carve.meta,
        diagnostics=carve.diagnostics,
    )
    with pytest.raises(RuntimeError, match="the read table holds 7 reads"):
        record_compliance(short, [("3", 3.0, True)])
    monkeypatch.setattr(
        compliance, "_threshold_in_spec_frac", lambda *a, **k: {"3": 0.95}
    )
    with pytest.raises(RuntimeError, match="outside its bounds"):
        record_compliance(carve, [("3", 3.0, True)])


@pytest.mark.unit
def test_a_carve_cut_on_another_ladder_or_without_its_sigma_is_refused():
    """Oracle: `windows.margin`, the one in-spec definition, and the carve's own settings.
    Another label set, another value or another side under the same label are each
    refused; so is a carve cut with use_uncertainty off, whose read table holds no sigma
    for the bounds the record may declare."""
    carve = _hand_carve(_SIGMA)
    for ladder in ([("4", 3.0, True)], [("3", 3.5, True)], [("3", 3.0, False)]):
        with pytest.raises(ValueError, match="carved"):
            record_compliance(carve, ladder)
    with pytest.raises(ValueError, match="use_uncertainty off"):
        record_compliance(
            _hand_carve(_SIGMA, use_uncertainty=False), [("3", 3.0, True)]
        )


# --------------------------------------------- S4: the composite and the paper's tables


def _job(path: Path) -> object:
    spec = importlib.util.spec_from_file_location(f"_compliance_{path.stem}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _composite():
    return _job(JOBS / "composite" / "compliance_survey.py")


@pytest.mark.policy
def test_the_composite_surveys_every_ledger_record_through_one_bundle_per_job():
    """Oracle: the ledger jobs' graphs. The survey takes one ref per ledger job, to the
    bundle that job persists; each bundle holds that job's carve nodes, the ones its
    ledgers score, in order; the survey's rows are every record of the family, in job
    order, each on its dataset's device, on the T2* ladder the bundles carry."""
    from quebra.recipes import _carve_stage

    family = {j.job_id: j for j in by_family(JOBS, "check_ledger")}
    composite = _composite()
    survey = composite.job.dag["compliance_survey"]
    assert all(isinstance(ref, ArtifactRef) for ref in survey.inputs)
    jobs = [ref.included.job for ref in survey.inputs]
    assert sorted(job.name for job in jobs) == sorted(family)
    assert [job.name for job in jobs[:6]] == [
        f"check_ledger_6d2s_q{q}" for q in range(1, 7)
    ]

    records, devices = [], []
    for ref, job in zip(survey.inputs, jobs):
        (bundle,) = [n for n in job.dag.values() if n.fn is collect_carves]
        persisted = [
            s.name
            for s in job.sinks
            if getattr(s, "node", None) is not None and s.node.node_id == bundle.node_id
        ]
        assert persisted == [ref.node_name] and ref.node_name.endswith(
            CARVES_DATA_SUFFIX
        )
        carves = [
            n
            for n in job.dag.values()
            if n.fn is _carve_stage and n.kwargs["stage"] == "windows"
        ]
        assert [r.node_id for r in bundle.inputs] == [c.node_id for c in carves]
        assert list(bundle.kwargs["records"]) == [c.kwargs["record"] for c in carves]
        assert list(bundle.kwargs["thresholds"]) == list(T2STAR_THRESHOLDS)
        for carve in carves:
            ledgers = [
                n
                for n in job.dag.values()
                if n.node_id.startswith("check_ledger_")
                and n.inputs[0].node_id == carve.node_id
            ]
            assert len(ledgers) == 1, carve.kwargs["record"]
            load = carve
            while "dataset" not in load.kwargs:
                load = job.dag[load.inputs[0].node_id]
            devices.append(load.kwargs["dataset"].device)
        records += list(bundle.kwargs["records"])
    assert list(survey.kwargs["records"]) == records
    assert list(survey.kwargs["devices"]) == devices
    assert len(records) == 28 + 22
    assert "thresholds" not in survey.kwargs, "the ladder is the bundles' own"


def _bundled_job(path: Path, device: str | None) -> SimpleNamespace:
    """An included ledger job in miniature: one recorded carve and its bundle."""
    job = Job("bundle_probe")
    dataset = Dataset(
        path=str(path),
        schema=track912Schema,
        qubit=1,
        device=device,
        extra={"run_name": "q1_010423", "run_start_unix_s": 1.68e9},
    )
    _, carve = wire_t2star_carve(job, dataset=dataset, on_failure="record")
    bundle = job.step(
        collect_carves,
        carve,
        name="record_carves",
        records=("q1_010423",),
        thresholds=T2STAR_THRESHOLDS,
    )
    job.materialize(bundle, name=f"bundle_probe{CARVES_DATA_SUFFIX}")
    return SimpleNamespace(job=job, alias="bundle_probe", ref=lambda name: name)


@pytest.mark.unit
def test_including_a_job_whose_rows_cannot_be_placed_raises(tmp_path):
    """Oracle: the job graph. A T2* job has a carve but no bundle, and a record whose
    dataset declares no device has no panel; the survey refuses both instead of drawing a
    record as absent or under a panel named None."""
    with pytest.raises(ValueError, match="not a ledger job"):
        included_carves(Job("survey_probe").include("t2star_q1_070423"))
    path = _track912_csv(tmp_path / "010423_r.csv", seed=1, gap=False)
    assert included_carves(_bundled_job(path, "probe"))[:2] == (
        ["q1_010423"],
        ["probe"],
    )
    with pytest.raises(ValueError, match="declares no device"):
        included_carves(_bundled_job(path, None))


@pytest.mark.unit
def test_a_bundle_refuses_a_run_name_attached_to_another_records_carve():
    """Oracle: each carve's own id. A run name out of step with its carve raises."""
    failure = RecordFailure(dataset_id="q1_010423", stage="load", reason="E: x")
    ladder = T2STAR_THRESHOLDS
    bundle = collect_carves(failure, records=("q1_010423",), thresholds=ladder)
    assert bundle.records == ("q1_010423",) and list(bundle.thresholds) == ladder
    with pytest.raises(ValueError, match="attached to the carve"):
        collect_carves(failure, records=("q1_020423",), thresholds=ladder)
    with pytest.raises(ValueError, match="run names for"):
        collect_carves(failure, records=(), thresholds=ladder)


# Capannelli et al., Supplementary notes, Supplementary Fig. 10 (6D2S, p. 10) and
# Supplementary Fig. 11 (2x2, p. 11): (qubit, DD/MM, shaded) per row, all in 2023. Shaded
# means the row is red in the table.
_SUPP_FIG_10 = [
    (1, "04/04", 0), (1, "05/04", 0), (1, "07/04", 0), (1, "10/04", 0), (1, "12/04", 0),
    (5, "12/04", 0), (6, "12/04", 0), (1, "21/04", 0), (2, "21/04", 0), (3, "21/04", 0),
    (1, "22/04", 0), (5, "07/06", 0), (6, "07/06", 0), (5, "08/06", 0), (6, "08/06", 0),
    (5, "09/06", 0), (6, "09/06", 0), (2, "27/06", 0), (1, "28/06", 1), (2, "28/06", 1),
    (4, "29/06", 0), (5, "29/06", 0), (6, "29/06", 0), (3, "11/07", 1),
]  # fmt: skip
_SUPP_FIG_11 = [
    (1, "09/06", 0), (2, "09/06", 1), (3, "09/06", 1), (1, "10/06", 0), (1, "20/06", 0),
    (1, "25/06", 1), (3, "26/06", 0), (1, "30/06", 1), (3, "30/06", 0), (1, "03/07", 1),
    (2, "03/07", 1), (2, "04/07", 1), (3, "04/07", 0), (1, "05/07", 1), (1, "07/07", 1),
    (1, "10/07", 1), (3, "11/07", 0), (2, "12/07", 1), (1, "22/08", 0), (2, "22/08", 0),
    (3, "22/08", 0), (4, "22/08", 0),
]  # fmt: skip


def _run_name(device: str, qubit: int, day_month: str) -> str:
    day, month = day_month.split("/")
    return f"{'2x2_' if device == '2x2' else ''}q{qubit}_{day}{month}23"


@pytest.mark.policy
def test_the_declared_table_and_shading_are_the_papers():
    """Oracle: the paper's tables, Supplementary Figs. 10 and 11. The composite declares
    their 24 and 22 rows on their devices, the 3 and 11 shaded rows, and the ledger
    records outside them are the four 6D2S records the tables do not list."""
    composite = _composite()
    table = [("6D2S", _run_name("6D2S", q, dm)) for q, dm, _ in _SUPP_FIG_10] + [
        ("2x2", _run_name("2x2", q, dm)) for q, dm, _ in _SUPP_FIG_11
    ]
    assert sorted(composite.PAPER_RECORDS) == sorted(table)
    shaded = [_run_name("6D2S", q, dm) for q, dm, s in _SUPP_FIG_10 if s] + [
        _run_name("2x2", q, dm) for q, dm, s in _SUPP_FIG_11 if s
    ]
    assert (len(_SUPP_FIG_10), len(_SUPP_FIG_11)) == (24, 22)
    assert sorted(composite.SHADED) == sorted(shaded) and len(shaded) == 3 + 11
    survey = composite.job.dag["compliance_survey"]
    extra = set(survey.kwargs["records"]) - {r for _, r in table}
    assert extra == {"q4_070723", "q4_100723", "q5_070723", "q6_070723"}


def _fake_bundle(records: dict[str, object], ladder=T2STAR_THRESHOLDS) -> RecordCarves:
    return RecordCarves(
        records=tuple(records),
        carves=tuple(records.values()),
        thresholds=tuple(ladder),
    )


@pytest.mark.unit
def test_the_builder_refuses_a_table_or_shading_it_cannot_honour():
    """Oracle: the declarations' contract. A paper record no ledger carved, one on the
    wrong device, a shaded record outside the table, rows out of step with the bundles,
    and bundles cut on two ladders each raise rather than drop or mislabel a row."""
    failure = RecordFailure(dataset_id="q1_010423", stage="load", reason="E: x")
    bundle = _fake_bundle({"q1_010423": failure})
    base = dict(
        records=("q1_010423",),
        devices=("6D2S",),
        paper_records=(("6D2S", "q1_010423"),),
        shaded=(),
        flag_text="f",
        table_source="s",
    )
    assert isinstance(build_compliance_survey(bundle, **base), ComplianceSurveyData)
    for change, match in (
        ({"paper_records": (("6D2S", "q1_020423"),)}, "no ledger job carved"),
        ({"paper_records": (("2x2", "q1_010423"),)}, "carved as a 6D2S"),
        ({"shaded": ("q1_020423",)}, "not in the paper's table"),
        ({"records": ("q1_020423",)}, "differ from the bundles"),
    ):
        with pytest.raises(ValueError, match=match):
            build_compliance_survey(bundle, **{**base, **change})
    other = _fake_bundle(
        {"q1_020423": RecordFailure(dataset_id="q1_020423", stage="load", reason="E")},
        ladder=T2STAR_THRESHOLDS[:-1],
    )
    with pytest.raises(ValueError, match="2 ladders, not one"):
        build_compliance_survey(
            bundle,
            other,
            **{**base, "records": ("q1_010423", "q1_020423"), "devices": ("6D2S",) * 2},
        )
    # Two records carved at another k, or another gap rule: one legend cannot state both.
    hand = [("3", 3.0, True)]
    for setting, match in (({"k": 2.0}, "k "), ({"gap_mult": 20.0}, "gap_mult ")):
        two = _fake_bundle(
            {
                "a": _hand_carve(_SIGMA, dataset_id="a"),
                "b": _hand_carve(_SIGMA, dataset_id="b", **setting),
            },
            ladder=hand,
        )
        with pytest.raises(ValueError, match=f"carved at {match}"):
            build_compliance_survey(
                two, **{**base, "records": ("a", "b"), "devices": ("6D2S",) * 2,
                        "paper_records": ()}
            )  # fmt: skip


# ------------------------------------------------------- S5: the figure draws the artifact


@pytest.mark.unit
def test_the_figure_draws_each_fraction_and_its_bounds_in_ladder_order(carves):
    """Oracle: the artifact. Every drawn record has one bar per computed rung at that
    rung's ladder position, its height the fraction, and a whisker spanning the bounds, or
    a dotted one spanning the row when the bounds were not computed; a record outside the
    paper's table has none and is named under its panel; a failed record's note is
    printed; the panel prints the artifact's counts; the legend names the flag and the
    missing bounds only when a drawn row carries them."""
    from quebra.plots.compliance_survey_plot import (
        ROW_HEIGHT,
        ComplianceSurveyPlot,
        bar_gid,
        bounds_gid,
        no_bounds_gid,
    )

    names = [str(c.meta["dataset_id"]) for _, c in carves]
    first = carves[0][1]
    # The first record again with every sigma dropped: its fractions stand, its bounds go.
    no_sigma = windows.WindowsResult(
        windows_in_spec=first.windows_in_spec,
        windows_out_of_spec=first.windows_out_of_spec,
        reads=first.reads.assign(sigma_v=np.nan),
        meta={**first.meta, "dataset_id": "q3_010423"},
        diagnostics=first.diagnostics,
    )
    failure = RecordFailure(dataset_id="q2_010423", stage="t2star", reason="E: x")
    bundle = RecordCarves(
        records=(*names, "q2_010423", "q3_010423"),
        carves=(*(c for _, c in carves), failure, no_sigma),
        thresholds=tuple(T2STAR_THRESHOLDS),
    )
    paper = (
        ("6D2S", names[0]),
        ("6D2S", names[1]),
        ("6D2S", "q2_010423"),
        ("6D2S", "q3_010423"),
    )
    declared = dict(
        records=bundle.records,
        devices=("6D2S",) * 5,
        paper_records=paper,
        flag_text="FLAG",
        table_source="SOURCE",
    )
    data = build_compliance_survey(bundle, shaded=(names[1],), **declared)
    n = len(T2STAR_THRESHOLDS)
    assert all(
        rung.status == COMPUTED == rung.bounds_status
        for _, carve in carves[:2]
        for rung in record_compliance(carve, T2STAR_THRESHOLDS).rungs
    )
    # Four drawn records: two computed throughout, one failed, one with no sigma.
    assert data.counts == {"6D2S": PanelCounts(4 * n, n, n)}
    fig = ComplianceSurveyPlot(name="probe").build_matplotlib(data)
    try:
        (ax,) = fig.axes
        bars = {p.get_gid(): p for p in ax.patches if p.get_gid()}
        whiskers = {c.get_gid(): c for c in ax.collections if c.get_gid()}
        drawn = 0
        for record in data.records:
            for j, rung in enumerate(record.compliance.rungs):
                key = (record.record, rung.threshold_label)
                if not record.in_paper_table or rung.in_spec_fraction is None:
                    assert not {
                        bar_gid(*key),
                        bounds_gid(*key),
                        no_bounds_gid(*key),
                    } & (set(bars) | set(whiskers))
                    continue
                bar = bars[bar_gid(*key)]
                assert bar.get_x() + bar.get_width() / 2 == pytest.approx(j)
                assert bar.get_height() == pytest.approx(
                    ROW_HEIGHT * rung.in_spec_fraction, rel=1e-12
                )
                base = bar.get_y()
                if rung.bounds_status == COMPUTED:
                    assert no_bounds_gid(*key) not in whiskers
                    span = (
                        ROW_HEIGHT * rung.in_spec_fraction_lower,
                        ROW_HEIGHT * rung.in_spec_fraction_upper,
                    )
                    segment = whiskers[bounds_gid(*key)].get_segments()[0]
                else:
                    assert bounds_gid(*key) not in whiskers
                    span = (0.0, ROW_HEIGHT)
                    segment = whiskers[no_bounds_gid(*key)].get_segments()[0]
                ((x0, y0), (x1, y1)) = segment
                assert x0 == x1 == pytest.approx(j)
                assert (y0 - base, y1 - base) == pytest.approx(
                    span, rel=1e-9, abs=1e-12
                )
                drawn += 1
        assert drawn == 3 * n
        labels = [t.get_text() for t in ax.get_xticklabels()]
        assert labels == [f"{v * 1e6:g}" for _, v, _ in T2STAR_THRESHOLDS]
        texts = " ".join(t.get_text() for t in ax.texts)
        assert record_label(names[2]) in texts and "not in the paper's" in texts
        assert "not computed: t2star failed: E: x" in texts
        assert f"not computed: {n} of {4 * n} rungs" in texts
        assert f"bounds not computed beside a drawn fraction: {n} rungs" in texts
        legend = [t.get_text() for t in fig.legends[0].get_texts()]
        assert "FLAG" in legend and "bounds not computed" in legend
    finally:
        plt.close(fig)

    plain = build_compliance_survey(
        RecordCarves(
            records=(*names, "q2_010423"),
            carves=(*(c for _, c in carves), failure),
            thresholds=tuple(T2STAR_THRESHOLDS),
        ),
        shaded=(),
        **{**declared, "records": (*names, "q2_010423"), "devices": ("6D2S",) * 4,
           "paper_records": paper[:3]},
    )  # fmt: skip
    fig = ComplianceSurveyPlot(name="probe").build_matplotlib(plain)
    try:
        legend = [t.get_text() for t in fig.legends[0].get_texts()]
        assert "FLAG" not in legend and "bounds not computed" not in legend
    finally:
        plt.close(fig)


@pytest.mark.unit
def test_the_label_of_a_row_is_the_independence_surveys():
    """Oracle: `ledger_recipe.record_label`, the row label of every survey grid."""
    failure = RecordFailure(dataset_id="2x2_q1_090623", stage="load", reason="E: x")
    data = build_compliance_survey(
        _fake_bundle({"2x2_q1_090623": failure}),
        records=("2x2_q1_090623",),
        devices=("2x2",),
        paper_records=(),
        shaded=(),
        flag_text="f",
        table_source="s",
    )
    assert data.records[0].label == record_label("2x2_q1_090623")
