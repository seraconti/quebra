"""The serial bench cells refuse a configuration they would mislabel.

A serial cell is labelled by what it asks for: a segment target encoded as a censoring
fraction, and for Arm EL a continuous, uncensored pair of segments. Both ways of asking for
something the generator would not deliver raise, rather than writing a row whose label
describes a layout no replicate had.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from jobs.bench import grid, runner
from jobs.bench.arms import ARM_EL, n_segments_for_censoring
from quebra.analyzers.checks.result import CLOCK_IN_SPEC

pytestmark = pytest.mark.unit


def test_a_segment_target_past_the_arms_cap_raises_in_serial_cells(monkeypatch):
    """Oracle: `n_segments_for_censoring` caps at `n // 5`, so 20 segments at n = 50 realise
    as 10. The shipped grid must round-trip to its own targets and the capped one must raise.
    """
    for cell in grid.serial_cells():
        if cell.arm != ARM_EL:
            assert n_segments_for_censoring(cell.n, cell.censoring) in (
                grid.SERIAL_SEGMENT_GRID
            )
    monkeypatch.setattr(grid, "SERIAL_SEGMENT_GRID", (3, 5, 20))
    with pytest.raises(ValueError, match="20 segments at n = 50 realise as 10"):
        grid.serial_cells()


@pytest.mark.parametrize(
    "asked", [{"quantised": True}, {"censoring": 0.03}], ids=["quantised", "censoring"]
)
def test_an_arm_el_cell_asking_for_what_it_cannot_draw_raises(asked):
    """Contract: Arm EL draws two continuous segments of its own, so a cell asking for
    quantised durations or a censoring target raises out of `run_cell` before its first
    replicate. A raise inside a replicate would be counted as a failure and the run would
    go on.
    """
    cell = grid.Cell(
        arm=ARM_EL,
        n=2 * grid.LEVEL_SHIFT_EVENTS_PER_SEGMENT,
        clock=CLOCK_IN_SPEC,
        kind=grid.KIND_SERIAL,
        shape=grid.SERIAL_SHAPE,
        rho=0.0,
        delta=0.4,
    )
    served = runner.generate(cell, np.random.default_rng(0))
    assert served.n_segments == 2
    with pytest.raises(ValueError, match="Arm EL draws two continuous segments"):
        runner.run_cell(dataclasses.replace(cell, n_replicates=1, **asked), 0)
