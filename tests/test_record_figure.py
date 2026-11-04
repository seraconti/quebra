"""Record figures: the carve held to its ledger, the survival block, and the figure jobs.

Oracles: hand tables for the cross-check; the real `ledger_step` on a synthetic record's
real carve for the run-set, the `(all)` expansion and each band's check outcome; the
artifact's own fields read back off the drawn figure; the job graphs for the one-ref rule
and the record list; and a run directory for the nested-run count.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import to_rgba
from matplotlib.text import Text

from quebra.analyzers.check_ledger import (
    VERDICT_ABSENT,
    VERDICT_NOT_COMPUTED,
    CheckLedger,
)
from quebra.analyzers.check_outcome import PERMUTATION_KEYS
from quebra.analyzers.checks.result import CLOCKS
from quebra.analyzers.event_table import SIDES
from quebra.analyzers.kaplan_meier import CONF_LEVEL
from quebra.analyzers.record_figure import (
    ALL_ROWS,
    build_survival_block,
    cross_check_carve,
    kaplan_meier_of,
    ledger_run_set,
)
from quebra.core.job import Job, _import_job
from quebra.core.reference import ArtifactRef
from quebra.core.runner import run_job
from quebra.core.types import RecordFailure
from quebra.ledger_recipe import LEDGER_KNOBS, ledger_step
from quebra.paper_tables import FLAG_TEXT, PAPER_RECORDS, SHADED
from quebra.plots.record_figure_plot import CELL_GID, SurvivalBlockPlot
from quebra.plots.theme import verdict_color
from quebra.record_figure_recipe import configure_record_figure_job
from quebra.recipes import T2STAR_THRESHOLDS
from tests.test_compliance_survey import _carve, _track912_csv

REPO = Path(__file__).resolve().parents[1]
BENCH = REPO / "jobs" / "bench" / "results" / "size_table.csv"
RECORD = "q1_020423"
FIGURE_JOBS = sorted((REPO / "jobs" / "composite").glob("record_figures_6d2s_q*.py"))


@pytest.fixture(scope="module")
def scored(tmp_path_factory) -> SimpleNamespace:
    """A synthetic record's real carve and the real ledger of it, at 99 permutations."""
    root = tmp_path_factory.mktemp("record_figure")
    _, carve = _carve(
        _track912_csv(root / "020423_r.csv", seed=2, gap=True, n=600), RECORD
    )
    ledger = ledger_step(
        carve,
        pd.read_csv(BENCH),
        thresholds=T2STAR_THRESHOLDS,
        **{**LEDGER_KNOBS, "n_permutations": 99},  # type: ignore[arg-type]
    )
    return SimpleNamespace(carve=carve, ledger=ledger)


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
    `(all)` row has a key of the derived run-set, and the run-set's every key has a row;
    a ledger run with C3 raises."""
    rows = scored.ledger.rows
    answered = rows[rows["check_id"] != ALL_ROWS]
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
    assert not {"q4_070723", "q4_100723", "q5_070723", "q6_070723"} & set(drawn_all)


@pytest.mark.unit
def test_each_block_is_handed_its_own_records_fit_windows_and_carve_kwargs(monkeypatch):
    """Oracle: the job graphs. Every `RecordNodes` a block receives names the fit and
    windows nodes its own record's carve wired (`t2star_<record>`, `windows_<record>`, the
    latter the input the cross-check holds), and carries the carve kwargs of that record's
    windows node in the included ledger job."""
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
    """Oracle: the run directory. A figure job over an invented two-record ledger job
    runs that job once, through its one ref, and renders a survival block per record."""
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
        assert list(run_dir.glob(f"{record}_survival_*.pdf")), sorted(
            p.name for p in run_dir.iterdir()
        )
