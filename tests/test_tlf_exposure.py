"""`tlf`'s dwell bookkeeping on a gapped record with skewed spacing, against hand arithmetic.

Oracle: a hand-built record whose every run, exposure and exit is counted below, read by
read. Its within-segment spacings are 70% at 10 s and 30% at 40 s, so the median spacing
(10 s) is not the length of a typical run's spacings and a dwell taken as reads times the
median is wrong. It holds one gap, with a state change across it, and cut runs at the
record's start, at the gap from both sides, and at its end. One of them (state 0, reads
10-11) has two reads and a 40 s spacing, so including or excluding it moves the numbers.

The two levels are exact, so the mixture assigns every read correctly: what is tested is
the exposure, exit and run bookkeeping, not mixture-based state assignment
(`test_tlf_dwell_law.py` says the same of its own oracle).
"""

from __future__ import annotations

import dataclasses

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from quebra.analyzers import tlf
from quebra.plots.tlf_plot import TLFPlot

pytestmark = pytest.mark.unit

LOW, HIGH = 0.0, 10.0

# (time in s, state). Segment A: reads 0-11. A 1000 s gap. Segment B: reads 12-21.
_READS = [
    (0, 0), (10, 0), (50, 1), (60, 1), (100, 1), (110, 0), (120, 0), (160, 0),
    (170, 1), (180, 1), (190, 0), (230, 0),
    (1230, 1), (1240, 1), (1280, 1), (1290, 0), (1300, 0), (1340, 1), (1350, 1),
    (1360, 0), (1370, 0), (1380, 0),
]  # fmt: skip

# Runs, by hand. State 0: reads 0-1 (start cut; exits at 50 s; 50 s), 5-7 (complete; 60 s),
# 10-11 (both ends inside segment A, but cut by the gap; 40 s), 15-16 (complete; 50 s),
# 19-21 (cut by the end of the record; 20 s). State 1: 2-4 (complete; 60 s), 8-9 (complete;
# 20 s), 12-14 (born after the gap, so cut; exits at 1290 s; 60 s), 17-18 (complete; 20 s).
# The change from read 11 (state 0) to read 12 (state 1) crosses the gap: no exit.
EXPOSURE_S0, EXITS_S0 = 50.0 + 60.0 + 40.0 + 50.0 + 20.0, 3
EXPOSURE_S1, EXITS_S1 = 60.0 + 20.0 + 60.0 + 20.0, 4
COMPLETE_S0, CUT_S0 = [60.0, 50.0], 3
COMPLETE_S1, CUT_S1 = [60.0, 20.0, 20.0], 1


def _record() -> tuple[np.ndarray, np.ndarray]:
    t = np.array([r[0] for r in _READS], dtype=float)
    values = np.where(np.array([r[1] for r in _READS]) == 1, HIGH, LOW)
    return values, t


@pytest.fixture(scope="module")
def result() -> tlf.TLFResult:
    values, t = _record()
    return tlf.run(values, t, seed=11, gap_mult=10.0)


def test_the_fixture_has_the_spacing_and_gap_it_claims() -> None:
    """Control on the oracle itself: 14 spacings of 10 s, 6 of 40 s and one gap."""
    _, t = _record()
    dt = np.diff(t)
    assert sorted(set(dt)) == [10.0, 40.0, 1000.0]
    assert (int((dt == 10.0).sum()), int((dt == 40.0).sum())) == (14, 6)


def test_exposure_exits_and_mean_dwell_match_hand_arithmetic(result) -> None:
    """Oracle: the run table above. Mean dwell = exposure / exits, per state."""
    assert result.exposure_state0_s == pytest.approx(EXPOSURE_S0, abs=1e-9)
    assert result.exposure_state1_s == pytest.approx(EXPOSURE_S1, abs=1e-9)
    assert (result.n_exits_state0, result.n_exits_state1) == (EXITS_S0, EXITS_S1)
    assert result.mean_dwell_state0_s == pytest.approx(EXPOSURE_S0 / EXITS_S0, abs=1e-9)
    assert result.mean_dwell_state1_s == pytest.approx(EXPOSURE_S1 / EXITS_S1, abs=1e-9)


def test_dwell_cv_is_over_complete_runs_only(result) -> None:
    """Oracle: the complete runs above. State 1 has three, CV = sd(60, 20, 20) / mean
    = 1.2 / sqrt(3); state 0 has two, too few for a CV, and three cut runs."""
    assert (result.n_complete_runs_state0, result.n_cut_runs_state0) == (
        len(COMPLETE_S0),
        CUT_S0,
    )
    assert (result.n_complete_runs_state1, result.n_cut_runs_state1) == (
        len(COMPLETE_S1),
        CUT_S1,
    )
    assert result.dwell_cv_state0 is None
    assert result.dwell_cv_state1 == pytest.approx(1.2 / np.sqrt(3.0), abs=1e-12)


def test_a_change_across_the_gap_is_not_a_transition(result) -> None:
    """Oracle: seven within-segment changes; the rate is over the observed time, 380 s,
    not over the 1380 s record span that includes the gap."""
    assert result.n_transitions == EXITS_S0 + EXITS_S1 == 7
    expected = 7 / ((EXPOSURE_S0 + EXPOSURE_S1) / 3600.0)
    assert result.switching_rate_per_hour == pytest.approx(expected, abs=1e-9)


def test_lobe_separation_is_in_hertz(result) -> None:
    """Oracle: the two levels, 10 Hz apart. A ppm conversion would divide by the mean."""
    assert result.lobe_separation_hz == pytest.approx(HIGH - LOW, abs=1e-6)


def test_gap_mult_reaches_the_gap_rule() -> None:
    """Oracle: the run table above with the gap closed. At gap_mult 200 the threshold is
    2000 s, so the 1000 s step is observed: the change across it is an eighth transition,
    and state 0's run at reads 10-11 now exits at read 12, 1040 s after it began."""
    values, t = _record()
    closed = tlf.run(values, t, seed=11, gap_mult=200.0)
    assert closed.diagnostics["n_gaps"] == 0
    assert closed.n_transitions == 8
    assert closed.exposure_state0_s == pytest.approx(EXPOSURE_S0 + 1000.0, abs=1e-9)
    assert (closed.n_exits_state0, closed.n_cut_runs_state1) == (EXITS_S0 + 1, 0)


def test_a_dropped_read_is_counted_and_changes_nothing_else(result) -> None:
    """Oracle: the same record with one extra read whose value is NaN, at t = 5 s."""
    values, t = _record()
    with_nan = tlf.run(
        np.insert(values, 1, np.nan), np.insert(t, 1, 5.0), seed=11, gap_mult=10.0
    )
    assert with_nan.diagnostics["n_dropped_nonfinite"] == 1
    assert result.diagnostics["n_dropped_nonfinite"] == 0
    assert with_nan.exposure_state0_s == result.exposure_state0_s
    assert with_nan.n_transitions == result.n_transitions


def test_a_repeated_timestamp_inside_a_segment_raises() -> None:
    """Oracle: specification. A state change between two reads at one time would be an
    exit after zero exposure."""
    t = np.array([0.0, 10.0, 10.0, 20.0, 30.0, 40.0])
    with pytest.raises(ValueError, match="repeated timestamp"):
        tlf.run(np.array([0.0, 0.0, 10.0, 10.0, 0.0, 10.0]), t, seed=11)


def _figure_text(result: tlf.TLFResult, values: np.ndarray, t: np.ndarray) -> str:
    payload = {"result": result, "values_hz": values, "t_rel_s": t, "meta": {}}
    figure = TLFPlot("unit").build_matplotlib(payload)
    try:
        return "\n".join(text.get_text() for ax in figure.axes for text in ax.texts)
    finally:
        plt.close(figure)


def test_the_figure_prints_the_counts_from_the_artifact(result) -> None:
    """Oracle: the hand counts above, as the panel prints them."""
    values, t = _record()
    text = _figure_text(result, values, t)
    assert "reads: 22 used, 0 dropped (non-finite time or value)" in text
    assert (
        f"mean dwell, state 1: {EXPOSURE_S1 / EXITS_S1:.2f} s ({EXITS_S1} exits)"
        in text
    )
    assert "dwell CV, state 0: n/a (2 complete, 3 cut runs)" in text


def test_a_failed_fit_reads_as_failed_in_the_figure(result) -> None:
    """Oracle: specification. A fit that did not happen must not print counts of None."""
    values, t = _record()
    failed = dataclasses.replace(
        result,
        gmm2=None,
        fit_failed=True,
        bic_delta=None,
        **dict.fromkeys(tlf._DYNAMICS_FIELDS),
    )
    text = _figure_text(failed, values, t)
    assert "mean dwell, state 0: n/a (fit failed)" in text
    assert "None" not in text


def test_the_histogram_draws_the_reads_the_fit_used(monkeypatch) -> None:
    """Oracle: the fit's own count. A read with a finite value at a NaN time is dropped by
    the fit, so the histogram must drop it too, or its count and the printed one disagree."""
    drawn: list[int] = []
    real_hist = matplotlib.axes.Axes.hist

    def spy(self, x, *args, **kwargs):
        drawn.append(len(np.asarray(x)))
        return real_hist(self, x, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "hist", spy)
    values, t = _record()
    values_bad, t_bad = np.insert(values, 1, 5.0), np.insert(t, 1, np.nan)
    fit = tlf.run(values_bad, t_bad, seed=11, gap_mult=10.0)
    _figure_text(fit, values_bad, t_bad)
    assert fit.diagnostics["n_reads"] == 22
    assert drawn == [22]
