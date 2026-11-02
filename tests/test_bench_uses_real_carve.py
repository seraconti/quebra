"""The bench must carve with the pipeline's own policy, not a lookalike.

`jobs/bench/carve.py` calls `windows.spacing`, `gap_flags`, `in_spec_mask` and `carve` directly
instead of `windows.run`, to skip a per-read DataFrame the bench never reads and an
unconditional print in an inner loop. That is a performance decision, and it is only sound
if the result is IDENTICAL to what `run` would have produced. This file is where that claim
is enforced rather than asserted in a comment.

The equality has to hold on the awkward inputs, not just a clean one: records with gaps,
records that start and end in spec, records that never leave spec, duplicate timestamps.
Each of those exercises a different branch of the birth/death taxonomy.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from quebra.analyzers import windows
from jobs.bench.carve import BENCH_WINDOW_COLUMNS, carve_windows

pytestmark = pytest.mark.statistical


THRESHOLD = 0.0


def _series(kind: str) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(4242)
    if kind == "plain":
        t = np.arange(200, dtype=float)
        return t, rng.standard_normal(200)
    if kind == "with_gaps":
        # Three blocks separated by jumps far beyond 10x the median spacing.
        t = np.concatenate(
            [np.arange(60), 400 + np.arange(60), 900 + np.arange(60)]
        ).astype(float)
        return t, rng.standard_normal(180)
    if kind == "all_in_spec":
        # No down_crossing anywhere: every window is born at scan_start and dies at
        # scan_end, the case where duration is measured to the last in-spec read.
        t = np.arange(50, dtype=float)
        return t, np.ones(50)
    if kind == "all_out_of_spec":
        t = np.arange(50, dtype=float)
        return t, -np.ones(50)
    if kind == "starts_and_ends_in_spec":
        t = np.arange(40, dtype=float)
        values = rng.standard_normal(40)
        values[0] = values[-1] = 1.0
        return t, values
    if kind == "with_nonfinite":
        # windows.run drops non-finite (t, value) pairs BEFORE measuring spacing, so a
        # single NaN can change the median and therefore the gap threshold.
        t = np.arange(30, dtype=float)
        values = rng.standard_normal(30)
        values[7] = np.nan
        t[19] = np.nan
        return t, values
    if kind == "backwards_time":
        t = np.arange(30, dtype=float)
        t[15] = 4.0
        return t, rng.standard_normal(30)
    if kind == "duplicate_timestamps":
        # Repeated timestamps drive the median spacing to zero unless only positive
        # steps are used; both paths must handle it the same way.
        t = np.repeat(np.arange(40, dtype=float), 2)
        return t, rng.standard_normal(80)
    raise ValueError(kind)


SERIES_KINDS = [
    "plain",
    "with_gaps",
    "all_in_spec",
    "all_out_of_spec",
    "starts_and_ends_in_spec",
    "duplicate_timestamps",
    "with_nonfinite",
]


def _reference_windows(t: np.ndarray, values: np.ndarray) -> pd.DataFrame:
    result = windows.run(
        windows.WindowsInputs(
            t_rel_s=t,
            values=values,
            thresholds=[("thr", THRESHOLD, True)],
        )
    )
    frame = result.windows_in_spec.reset_index(drop=True)
    return frame[BENCH_WINDOW_COLUMNS] if len(frame) else frame


@pytest.mark.parametrize("kind", SERIES_KINDS)
def test_bench_carve_matches_windows_run(kind):
    t, values = _series(kind)
    reference = _reference_windows(t, values)
    produced = carve_windows(t, values, THRESHOLD, big_values_good=True)

    assert len(produced) == len(reference), (
        f"{kind}: bench carved {len(produced)} windows, windows.run carved "
        f"{len(reference)}"
    )
    if not len(reference):
        return
    pd.testing.assert_frame_equal(
        produced.reset_index(drop=True)[BENCH_WINDOW_COLUMNS],
        reference.reset_index(drop=True)[BENCH_WINDOW_COLUMNS],
        check_dtype=False,
        obj=f"bench carve vs windows.run on {kind!r}",
    )


@pytest.mark.parametrize("kind", ["plain", "with_gaps"])
def test_gap_policy_is_not_reimplemented(kind):
    """The bench's gap threshold comes from `windows.spacing`, not from a local rule."""
    t, values = _series(kind)
    produced = carve_windows(t, values, THRESHOLD, big_values_good=True)
    reference = _reference_windows(t, values)
    assert set(produced["death_type"]) == set(reference["death_type"])
    assert set(produced["birth_type"]) == set(reference["birth_type"])


@pytest.mark.parametrize("kind", ["plain", "with_gaps", "starts_and_ends_in_spec"])
def test_a_divergent_bench_carve_is_actually_caught(kind):
    """Positive control, and it must exercise the REAL comparison.

    The previous version tampered with `reference` and compared it against `reference` -
    `carve_windows` was never called, so it was a test of `pandas.testing`. A mutant carve
    reproducing the exact divergence named below passed it. This one mutates the BENCH
    output and asserts the same comparison the real test performs, with the same options.
    """
    t, values = _series(kind)
    reference = _reference_windows(t, values)
    produced = carve_windows(t, values, THRESHOLD, big_values_good=True)
    assert len(reference) > 2, (
        "need a non-trivial reference for the control to mean much"
    )

    # The divergence the module docstring warns about: measuring a down_crossing window to
    # its last IN-SPEC read (the reference monolith's convention) rather than to the first
    # out-of-spec read.
    mutant = produced.copy()
    mutant["duration_s"] = mutant["duration_s"] - 1.0
    with pytest.raises(AssertionError):
        pd.testing.assert_frame_equal(
            mutant.reset_index(drop=True)[BENCH_WINDOW_COLUMNS],
            reference.reset_index(drop=True)[BENCH_WINDOW_COLUMNS],
            check_dtype=False,
        )


def test_the_compared_columns_include_the_ones_that_matter():
    """`BENCH_WINDOW_COLUMNS` defines the test's own coverage, so it needs pinning.

    Dropping `duration_s` from that list would make the equality test stop comparing
    durations and stay green - exactly the failure the control above is named for.
    """
    required = {
        "t_birth_s",
        "t_death_s",
        "duration_s",
        "birth_type",
        "death_type",
        "censored",
    }
    assert required <= set(BENCH_WINDOW_COLUMNS), (
        f"missing from the compared columns: {sorted(required - set(BENCH_WINDOW_COLUMNS))}"
    )


@pytest.mark.parametrize("kind", ["with_nonfinite", "backwards_time"])
def test_bench_carve_matches_windows_run_on_hostile_input(kind):
    """`windows.run` drops non-finite pairs before measuring spacing and raises on
    backwards time. `carve_windows` claims to be identical, so it must do both - the arms
    build times by cumsum of positive gaps and never reach these, but the claim is
    unconditional."""
    t, values = _series(kind)
    if kind == "backwards_time":
        with pytest.raises(ValueError):
            _reference_windows(t, values)
        with pytest.raises(ValueError):
            carve_windows(t, values, THRESHOLD, big_values_good=True)
        return
    pd.testing.assert_frame_equal(
        carve_windows(t, values, THRESHOLD, big_values_good=True).reset_index(
            drop=True
        )[BENCH_WINDOW_COLUMNS],
        _reference_windows(t, values).reset_index(drop=True)[BENCH_WINDOW_COLUMNS],
        check_dtype=False,
    )


@pytest.mark.parametrize("threshold", [float("nan"), float("inf")])
def test_bench_carve_refuses_a_non_finite_threshold_on_a_one_read_record(threshold):
    """Oracle: `windows.run`, which checks the ladder before its two-read early return.
    The bench claims identical output, so it must refuse the same input."""
    t, values = np.zeros(1), np.ones(1)
    with pytest.raises(ValueError, match="non-finite"):
        windows.run(
            windows.WindowsInputs(
                t_rel_s=t, values=values, thresholds=[("thr", threshold, True)]
            )
        )
    with pytest.raises(ValueError, match="non-finite"):
        carve_windows(t, values, threshold, big_values_good=True)


@pytest.mark.parametrize(
    ("clock", "builder_clock"), [("in_spec", "side"), ("calendar", "calendar")]
)
def test_a_carved_arm_segments_through_the_pipelines_clock(clock, builder_clock):
    """Oracle: the pipeline's segment builder. A carved bench arm labels its clocks as the
    bench tables do (`in_spec`, `calendar`) and must reach the builder's side and calendar
    clocks, written out here rather than read back from `BENCH_CLOCK`: its segments are the
    builder's segments of the bench's own carve, the in-spec bag, watched until one read
    interval past the last read."""
    from jobs.bench import arms
    from quebra.analyzers.checks._multiprocess import segments_from_windows

    values = np.random.default_rng(17).standard_normal(400)
    got, _ = arms._carved_segments(values, clock, 2)
    t = np.arange(len(values), dtype=float) * arms.READ_DT_S
    threshold = float(np.quantile(values, 1.0 - arms.IN_SPEC_P))
    want, _ = segments_from_windows(
        carve_windows(t, values, threshold, big_values_good=True),
        clock=builder_clock,
        min_events=2,
        observation_end_s=len(values) * arms.READ_DT_S,
    )
    assert got, "premise: the arm yields segments"
    assert [(list(s.x), s.tau) for s in got] == [(list(s.x), s.tau) for s in want]
