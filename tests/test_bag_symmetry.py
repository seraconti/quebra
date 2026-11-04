"""The side clock on both bags: a record and its mirror image about the threshold.

Oracle: reflection symmetry. Reflecting every read about the threshold (`v -> 2T - v`) turns
the record's out-of-spec reads into the mirror's in-spec reads and back, read for read, when
no read sits exactly on the threshold (a read on it is in spec either way, `windows.margin`,
so the mirror is not exact there). The mirror's in-spec bag is then the record's out-of-spec
bag window for window: the same reads, durations, censoring and read gaps. Its ledger rows
must equal the record's bag for bag in every column but `side`, on both clocks: statistics,
p-values (one permutation stream per threshold and clock, `check_ledger.stream_for`, serves
both bags) and bench cells.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from quebra.analyzers import check_ledger, windows
from quebra.analyzers.checks.result import CLOCKS
from quebra.analyzers.event_table import SIDES
from quebra.analyzers.windows import SIDE_IN_SPEC, SIDE_OUT_OF_SPEC

pytestmark = pytest.mark.properties

SIZE_TABLE = (
    Path(__file__).resolve().parents[1]
    / "jobs"
    / "bench"
    / "results"
    / "size_table.csv"
)
THRESHOLD = 3.0


def _record() -> tuple[np.ndarray, np.ndarray]:
    """An AR(1) series around the threshold with two read gaps; no read on the threshold."""
    rng = np.random.default_rng(29)
    n = 1400
    dt = rng.uniform(40.0, 80.0, n)
    dt[[300, 900]] += [6000.0, 15000.0]
    t = np.cumsum(dt) - dt[0]
    z = np.zeros(n)
    for i in range(1, n):
        z[i] = 0.85 * z[i - 1] + rng.normal()
    values = THRESHOLD + 0.2 + 0.9 * z / np.std(z)
    return t, values


def _rows(t: np.ndarray, values: np.ndarray, big_values_good: bool) -> pd.DataFrame:
    ladder = [("3", THRESHOLD, big_values_good)]
    carved = windows.run(
        windows.WindowsInputs(
            t_rel_s=t, values=values, thresholds=ladder, dataset_id="mirror_probe"
        )
    )
    inputs = check_ledger.make_inputs_from_windows(
        carved,
        pd.read_csv(SIZE_TABLE),
        thresholds=ladder,
        alpha=0.05,
        min_events_pass=35,
        tie_cutoff_distinct=5,
        lag_max=5,
        n_permutations=99,
        seed=5,
    )
    return check_ledger.run(inputs).rows


def _bag(rows: pd.DataFrame, side: str) -> pd.DataFrame:
    return rows[rows["side"] == side].drop(columns="side").reset_index(drop=True)


@pytest.mark.parametrize("big_values_good", [True, False])
def test_a_records_bags_are_its_mirrors_bags_swapped(big_values_good):
    """Oracle: reflection symmetry (module docstring). Both directions, every column but
    `side`, exactly. A big-good and a small-good ladder, so the margin's sign is covered."""
    t, values = _record()
    mirror = 2.0 * THRESHOLD - values
    assert not np.any(values == THRESHOLD), "premise: no read on the threshold"
    record_rows = _rows(t, values, big_values_good)
    mirror_rows = _rows(t, mirror, big_values_good)

    for side, other in (
        (SIDE_IN_SPEC, SIDE_OUT_OF_SPEC),
        (SIDE_OUT_OF_SPEC, SIDE_IN_SPEC),
    ):
        pd.testing.assert_frame_equal(
            _bag(record_rows, side), _bag(mirror_rows, other), check_exact=True
        )
    # Premise: each bag answers on each clock, so the equality is not two blank tables.
    answered = record_rows.dropna(subset=["p_value"])
    assert set(zip(answered["side"], answered["clock"])) == {
        (side, clock) for side in SIDES for clock in CLOCKS
    }
    # C1 and C2 answer on both bags and both clocks too, with every side-clock segment that
    # ends at an event read through its first n - 1 events.
    trend = answered[
        answered["check_id"].isin(["c1_lewis_robinson", "c2_anderson_darling"])
    ]
    assert set(zip(trend["side"], trend["clock"], trend["check_id"])) == {
        (side, clock, check)
        for side in SIDES
        for clock in CLOCKS
        for check in ("c1_lewis_robinson", "c2_anderson_darling")
    }
    # And the bags differ from each other, so the swap is a real exchange.
    assert not _bag(record_rows, SIDE_IN_SPEC).equals(
        _bag(record_rows, SIDE_OUT_OF_SPEC)
    )
