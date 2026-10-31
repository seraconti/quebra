"""The battery's row sets: the pipeline carries the permutation rows only, the bench all.

Oracle: issue #19, written out here rather than read back off the code. The thesis reports
permutation calibration only, so the rows a ledger, the survey and a band carry are the five
permutation rows. The bench keeps the asymptotic rows of C1, C2 and CvM, the evidence for
that choice, in the order its tables were written, which must not move.
"""

from __future__ import annotations

import numpy as np
import pytest

from quebra.analyzers.checks import battery
from quebra.analyzers.checks._permutation import block_permutations
from quebra.analyzers.checks.result import Segment

pytestmark = pytest.mark.unit

PERMUTATION_ROWS = [
    ("c1_lewis_robinson", "permutation", ""),
    ("c2_anderson_darling", "permutation", ""),
    ("c5_rank_autocorr", "permutation", "studentized"),
    ("c6_exchangeability", "permutation", ""),
    ("cvm_cramer_von_mises", "permutation", ""),
]
BENCH_ROWS = [
    ("c1_lewis_robinson", "asymptotic", ""),
    ("c1_lewis_robinson", "permutation", ""),
    ("c2_anderson_darling", "asymptotic", ""),
    ("cvm_cramer_von_mises", "asymptotic", ""),
    ("c2_anderson_darling", "permutation", ""),
    ("c5_rank_autocorr", "permutation", "studentized"),
    ("c6_exchangeability", "permutation", ""),
    ("cvm_cramer_von_mises", "permutation", ""),
    ("c5_rank_autocorr", "permutation", "unstudentized"),
]


def _segments(sizes, seed):
    """Exponential durations, censored after the last event so `tau > T_N`."""
    rng = np.random.default_rng(seed)
    segments = []
    for n in sizes:
        x = rng.exponential(1.0, size=n)
        segments.append(Segment(x=x, tau=float(np.sum(x)) + 0.5))
    return segments


def test_the_pipeline_battery_returns_the_five_permutation_rows_only() -> None:
    """The default battery, which every ledger, the survey and a band read, returns the
    five permutation rows in order, on one segment and on two, and no asymptotic row."""
    assert list(battery.ROW_KEYS) == PERMUTATION_ROWS
    for sizes in ([40], [25, 30]):
        segments = _segments(sizes, 3)
        perm = block_permutations(sizes, 99, np.random.default_rng(4))
        rows = [battery.row_key(r) for r in battery.run_battery(segments, perm=perm)]
        assert rows == PERMUTATION_ROWS, sizes


def test_the_bench_path_keeps_every_row_in_its_tables_order() -> None:
    """The bench asks for the asymptotic rows and C5's unstudentized row, and gets the nine
    rows of its tables in their order on one segment."""
    assert list(battery.BENCH_ROW_KEYS) == BENCH_ROWS
    segments = _segments([40], 3)
    perm = block_permutations([40], 99, np.random.default_rng(4))
    rows = battery.run_battery(
        segments, perm=perm, include_asymptotic=True, include_c5_raw=True
    )
    assert [battery.row_key(r) for r in rows] == BENCH_ROWS


@pytest.mark.parametrize("sizes", [[40], [25, 30]])
@pytest.mark.parametrize("tau_checks", [True, False])
def test_leaving_the_asymptotic_rows_out_moves_no_other_row(sizes, tau_checks) -> None:
    """Oracle: the shared permutation set. The asymptotic rows take no permutation set and
    draw nothing, so every row the pipeline keeps is bit-identical with and without them."""
    segments = _segments(sizes, 5)
    perm = block_permutations(sizes, 99, np.random.default_rng(6))
    kept = battery.run_battery(segments, perm=perm, include_tau_checks=tau_checks)
    full = {
        battery.row_key(r): r
        for r in battery.run_battery(
            segments, perm=perm, include_tau_checks=tau_checks, include_asymptotic=True
        )
    }
    assert kept, "control: the battery returned rows"
    for row in kept:
        twin = full[battery.row_key(row)]
        assert (row.statistic, row.p_value) == (twin.statistic, twin.p_value)
