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

import numpy as np
import pytest

from quebra.analyzers import tlf

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
    assert result.exposure_s0 == pytest.approx(EXPOSURE_S0, abs=1e-9)
    assert result.exposure_s1 == pytest.approx(EXPOSURE_S1, abs=1e-9)
    assert (result.n_exits_s0, result.n_exits_s1) == (EXITS_S0, EXITS_S1)
    assert result.mean_dwell_s0 == pytest.approx(EXPOSURE_S0 / EXITS_S0, abs=1e-9)
    assert result.mean_dwell_s1 == pytest.approx(EXPOSURE_S1 / EXITS_S1, abs=1e-9)


def test_dwell_cv_is_over_complete_runs_only(result) -> None:
    """Oracle: the complete runs above. State 1 has three, CV = sd(60, 20, 20) / mean
    = 1.2 / sqrt(3); state 0 has two, too few for a CV, and three cut runs."""
    assert (result.n_complete_runs_s0, result.n_cut_runs_s0) == (
        len(COMPLETE_S0),
        CUT_S0,
    )
    assert (result.n_complete_runs_s1, result.n_cut_runs_s1) == (
        len(COMPLETE_S1),
        CUT_S1,
    )
    assert result.dwell_cv_s0 is None
    assert result.dwell_cv_s1 == pytest.approx(1.2 / np.sqrt(3.0), abs=1e-12)


def test_a_change_across_the_gap_is_not_a_transition(result) -> None:
    """Oracle: seven within-segment changes; the rate is over the observed time, 380 s,
    not over the 1380 s record span that includes the gap."""
    assert result.n_transitions == EXITS_S0 + EXITS_S1 == 7
    expected = 7 / ((EXPOSURE_S0 + EXPOSURE_S1) / 3600.0)
    assert result.switching_rate_per_hour == pytest.approx(expected, abs=1e-9)


def test_lobe_separation_is_in_hertz(result) -> None:
    """Oracle: the two levels, 10 Hz apart. A ppm conversion would divide by the mean."""
    assert result.lobe_separation_hz == pytest.approx(HIGH - LOW, abs=1e-6)
