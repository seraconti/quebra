"""The death-time convention, and the non-finite ordering that protects it.

A window dies at the FIRST OUT-OF-SPEC read, so the interval spanning the crossing is
inside the lifetime. This predates the window tables and must not drift: the brief for
the carve port names this trace explicitly as the thing that must not change.

The reference trace ends out-of-spec, so it exercises only down_crossing. A second case
covers scan_end, where the pre-existing `_collect_windows` double-counted the final
interval (an all-in-spec 2-hour span returned 180 minutes, not 120).
"""

from __future__ import annotations

import numpy as np
import pytest

from quebra.analyzers import windows

pytestmark = pytest.mark.statistical


MINUTE_S = 60.0
THRESHOLDS = [("3", 3.0, True)]


def _durations_min(t_rel_s: np.ndarray, values: np.ndarray) -> list[float]:
    result = windows.run(
        windows.WindowsInputs(
            t_rel_s=t_rel_s, values=values, thresholds=THRESHOLDS, dataset_id="unit"
        )
    )
    return [float(d) / MINUTE_S for d in result.windows_in_spec["duration_s"]]


def test_death_at_first_out_of_spec_read() -> None:
    # values crossing down twice; 1-minute spacing; threshold 3.
    t = np.arange(7, dtype=float) * MINUTE_S
    v = np.array([0, 5, 5, 0, 0, 5, 0], dtype=float)
    assert _durations_min(t, v) == [2.0, 1.0]


def test_death_types_of_that_trace() -> None:
    t = np.arange(7, dtype=float) * MINUTE_S
    v = np.array([0, 5, 5, 0, 0, 5, 0], dtype=float)
    result = windows.run(
        windows.WindowsInputs(
            t_rel_s=t, values=v, thresholds=THRESHOLDS, dataset_id="unit"
        )
    )
    assert list(result.windows_in_spec["death_type"]) == [
        "down_crossing",
        "down_crossing",
    ]
    assert list(result.windows_in_spec["birth_type"]) == ["up_crossing", "up_crossing"]
    assert not result.windows_in_spec["censored"].any()


def test_scan_end_window_is_not_inflated() -> None:
    # An all-in-spec 2-hour span is a 120-minute window, censored at scan_end.
    # The pre-existing carve returned 180 here by adding the final interval twice.
    t = np.array([0.0, 3600.0, 7200.0])
    v = np.array([5.0, 5.0, 5.0])
    assert _durations_min(t, v) == [120.0]
    result = windows.run(
        windows.WindowsInputs(
            t_rel_s=t, values=v, thresholds=THRESHOLDS, dataset_id="unit"
        )
    )
    assert list(result.windows_in_spec["death_type"]) == ["scan_end"]
    assert bool(result.windows_in_spec["censored"].iloc[0]) is True


def test_non_finite_values_do_not_fabricate_a_crossing() -> None:
    # A failed Ramsey fit is NaN. NaN >= threshold is False, so leaving it in would
    # look like a down-crossing. It must be dropped before the mask is taken.
    t = np.arange(5, dtype=float) * MINUTE_S
    v = np.array([5.0, 5.0, np.nan, 5.0, 5.0])
    result = windows.run(
        windows.WindowsInputs(
            t_rel_s=t, values=v, thresholds=THRESHOLDS, dataset_id="unit"
        )
    )
    assert len(result.windows_in_spec) == 1
    assert result.diagnostics["n_reads_dropped_nonfinite"] == 1
    assert list(result.windows_in_spec["death_type"]) == ["scan_end"]


def test_backwards_time_raises_rather_than_yielding_negative_lifetimes() -> None:
    # duration = t_death - t_birth, so a descending time array silently produced
    # NEGATIVE durations that flowed onto the survival curve: finite, plausible, wrong.
    t = np.arange(7, dtype=float)[::-1] * MINUTE_S
    v = np.array([5, 5, 0, 5, 5, 5, 0], dtype=float)
    with pytest.raises(ValueError, match="non-decreasing"):
        windows.run(
            windows.WindowsInputs(
                t_rel_s=t, values=v, thresholds=THRESHOLDS, dataset_id="unit"
            )
        )


def test_no_positive_spacing_raises_rather_than_inventing_one() -> None:
    # With every timestamp duplicated there is no observed spacing to scale the gap
    # threshold by. The reference monolith falls back to 1.0; that is a fabricated
    # physical scale in every window boundary downstream.
    t = np.zeros(6, dtype=float)
    v = np.array([5, 5, 0, 5, 5, 5], dtype=float)
    with pytest.raises(ValueError, match="read spacing"):
        windows.run(
            windows.WindowsInputs(
                t_rel_s=t, values=v, thresholds=THRESHOLDS, dataset_id="unit"
            )
        )


def test_duplicate_timestamps_are_still_fine_alongside_real_spacing() -> None:
    # The guard must not reject the case the positive-only median exists FOR.
    t = np.repeat(np.arange(6, dtype=float), 2) * MINUTE_S
    v = np.tile(np.array([5.0, 5.0]), 6)
    result = windows.run(
        windows.WindowsInputs(
            t_rel_s=t, values=v, thresholds=THRESHOLDS, dataset_id="unit"
        )
    )
    assert result.diagnostics["n_nonpositive_steps"] == 6
    assert result.diagnostics["median_spacing_s"] == MINUTE_S


def test_non_finite_drop_happens_before_spacing() -> None:
    # Dropping the middle reads widens the interval past the gap threshold; the gap
    # must be seen on the SURVIVING reads, not on the raw ones.
    # 1-minute spacing -> gap threshold 600 s, so the hole must exceed 10 reads.
    t = np.arange(20, dtype=float) * MINUTE_S
    v = np.full(20, 5.0)
    v[3:15] = np.nan
    result = windows.run(
        windows.WindowsInputs(
            t_rel_s=t, values=v, thresholds=THRESHOLDS, dataset_id="unit"
        )
    )
    assert result.diagnostics["n_gaps"] == 1
    assert list(result.windows_in_spec["death_type"]) == ["gap_start", "scan_end"]
    assert list(result.windows_in_spec["birth_type"]) == ["scan_start", "gap_resume"]
