"""The infidelity carve wired in a job that includes a Ramsey job with fidelity.

Oracles: the included Ramsey job's own graph (its fidelity carve node, the sink it persists
the fidelity result under, its declared Rabi row and its Dataset), and the names the paper's
dataset tables give their records (`paper_tables.PAPER_RECORDS`).
"""

from __future__ import annotations

import re
from datetime import date

import pytest

from quebra.core.job import Job
from quebra.core.reference import ArtifactRef
from quebra.infidelity_recipe import paper_record_name, wire_infidelity_carve
from quebra.paper_tables import PAPER_RECORDS
from quebra.recipes import _fidelity_windows

pytestmark = pytest.mark.unit

RAMSEY = "ramsey_2x2_q1_030723"


def test_the_carve_is_the_ramsey_jobs_own_on_the_result_it_persists() -> None:
    """Oracle: the included job's graph and the paper's tables. One ref, to the fidelity
    result under the sink the Ramsey job persists it with; the carve is `_fidelity_windows`
    on that ref with the Ramsey carve node's kwargs; the record's run name is the one its job
    declares, and its paper name, from the declared Rabi row, is the 2x2 row the tables
    list."""
    probe = Job("infidelity_probe")
    included = probe.include(RAMSEY)
    ramsey = included.job
    (own,) = [n for n in ramsey.dag.values() if n.fn is _fidelity_windows]
    (persisted,) = [
        sink.name
        for sink in ramsey.sinks
        if getattr(sink, "node", None) is not None
        and sink.node.node_id == own.inputs[0].node_id
    ]
    carve = wire_infidelity_carve(probe, included)
    assert isinstance(carve.result, ArtifactRef) and carve.result.node_name == persisted
    node = probe.dag[carve.windows.node_id]
    assert node.fn is _fidelity_windows and node.inputs == [carve.result]
    assert node.kwargs == own.kwargs == carve.carve_kwargs
    assert node.node_id == "fidelity_windows_2x2_q1_030723"
    # The run name `jobs/active/ramsey_2x2_q1_030723.py` declares on its Dataset.
    assert carve.run_name == "q1_2x2_16h_0307_dataset"
    assert (carve.device, carve.paper_name) == ("2x2", "2x2_q1_030723")
    assert (carve.device, carve.paper_name) in PAPER_RECORDS


def test_a_job_that_is_not_one_ramsey_record_with_fidelity_raises() -> None:
    """Oracle: the wiring rule. A T2* job has no fidelity carve, and a Ramsey job that
    persists its fidelity result under a second sink leaves the ref ambiguous: both raise
    when the job is built."""
    probe = Job("infidelity_probe")
    with pytest.raises(ValueError, match="0 fidelity carves, not one"):
        wire_infidelity_carve(probe, probe.include("t2star_q1_070423"))
    included = probe.include(RAMSEY)
    (own,) = [n for n in included.job.dag.values() if n.fn is _fidelity_windows]
    included.job.materialize(own.inputs[0], name="a_second_copy")
    with pytest.raises(ValueError, match="under 2 sinks, not one"):
        wire_infidelity_carve(probe, included)


def test_a_record_name_is_the_one_the_papers_tables_give_it() -> None:
    """Oracle: the 46 rows of the paper's dataset tables. Each row's device, qubit and date
    (all 2023) give back that row's name: 6D2S rows carry no prefix, 2x2 rows `2x2_`. A
    device the tables do not list raises."""
    for device, name in PAPER_RECORDS:
        match = re.fullmatch(r"(?:2x2_)?q(\d+)_(\d{2})(\d{2})23", name)
        assert match, name
        qubit, day, month = (int(g) for g in match.groups())
        assert paper_record_name(device, qubit, date(2023, month, day)) == name
    assert len(PAPER_RECORDS) == 46
    with pytest.raises(ValueError, match="list no device"):
        paper_record_name("3x3", 1, date(2023, 7, 3))
