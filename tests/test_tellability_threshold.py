"""Tellability's per-threshold rows against known answers.

Oracles:
- Q2.2: the definition, enumerated by hand: thin every segment from each of the m offsets in
  turn, count kept pairs and down-crossings, average over the offsets; the exponent against
  an in-test least-squares fit. Independent identically distributed reads give gamma = 1 in
  expectation, as a statistical control.
- Q3.4: a hand classification of every crossing of two built records, at k = 1 and 2.
- Q4.1: hand counts of reads with a successor and transitions, in both forms; and 1 / beta
  for a two-state chain, as a ratio of expectations.
- The in-spec fraction bounds: hand arithmetic, interval by interval; the carve's own
  fraction lies between them.
- rho_k: (1 - alpha - beta)^k for a stationary two-state chain, and hand arithmetic on two
  segments.
- The ladder: L0's, section 4.2.
"""

from __future__ import annotations

import numpy as np
import pytest

from quebra.analyzers import tellability as tl
from quebra.analyzers.windows import observed_reads

# --------------------------------------------------------------------------- Q2.2


def _layout() -> tuple[np.ndarray, np.ndarray]:
    """Segments of 37, 61, 97 and 11 reads, 10 s apart inside, 5000 s gaps between.

    In spec, segment by segment: 37 reads in, out at 20-24, in again to its end; 61 reads
    out at 0-13 (the step from the last read of the first segment, in spec, to this one is
    ACROSS A GAP and is no crossing), in at 14-53, out at 54-60 (a crossing 7 reads from the
    end); 97 reads in, out at 53-62, in again; 11 reads in at 0-5, out at 6-10, a crossing in
    a segment shorter than rungs 13 and 20.
    """
    lengths = (37, 61, 97, 11)
    t, ins, start = [], [], 0.0
    for length in lengths:
        t.append(start + np.arange(length) * 10.0)
        start = t[-1][-1] + 5000.0
    ins = np.r_[
        np.isin(np.arange(37), np.arange(20, 25), invert=True),
        np.isin(np.arange(61), np.arange(14, 54)),
        np.isin(np.arange(97), np.arange(53, 63), invert=True),
        np.arange(11) < 6,
    ]
    return np.concatenate(t), ins


def _by_offsets(lengths, ins, m):
    """n and P by the definition: every offset, every segment, every kept pair."""
    n_total = p_total = 0
    starts = np.r_[0, np.cumsum(lengths)[:-1]]
    for offset in range(m):
        for start, length in zip(starts, lengths):
            kept = start + np.arange(offset, length, m)
            pairs = list(zip(kept[:-1], kept[1:]))
            p_total += len(pairs)
            n_total += sum(1 for i, j in pairs if ins[i] and not ins[j])
    return n_total / m, p_total / m


@pytest.mark.unit
def test_the_count_exponent_on_a_hand_built_layout() -> None:
    """Oracle: the definition enumerated by hand (`_by_offsets`), P_avg(m) = sum_s (L_s -
    min(m, L_s)) / m, count_ratio = [n/P](m) / [n/P](1), and gamma against an in-test fit."""
    t, ins = _layout()
    lengths = (37, 61, 97, 11)
    row = tl.count_exponent_row(observed_reads(t, np.ones(t.size)), ins)
    expected = [_by_offsets(lengths, ins, m) for m in tl.DECIMATION_LADDER]
    for m, (n, p), got_n, got_p in zip(
        tl.DECIMATION_LADDER, expected, row.down_crossings, row.kept_pairs
    ):
        assert got_n == pytest.approx(n, rel=1e-15)
        assert got_p == pytest.approx(p, rel=1e-15)
        assert got_p == pytest.approx(
            sum(L - min(m, L) for L in lengths) / m, rel=1e-15
        )
    base = expected[0][0] / expected[0][1]
    for (n, p), ratio in zip(expected, row.count_ratio):
        assert ratio == pytest.approx(n / p / base, rel=1e-12)
    x = np.log(tl.DECIMATION_LADDER)
    y = np.log([n / p for n, p in expected])
    assert row.count_exponent == pytest.approx(1.0 - np.polyfit(x, y, 1)[0], abs=1e-12)
    assert row.rungs_excluded == ()


@pytest.mark.statistical
def test_independent_reads_give_a_count_exponent_of_one() -> None:
    """Oracle: analytic, as a statistical control. Independent reads in spec with
    probability 0.3 have n / P = 0.21 at every m, so gamma = 1 in expectation on any layout.
    Over seeds 0-299, each drawing the layout and the reads from `default_rng(seed)` as this
    test does, gamma measured 1.0001 with sd 0.0054 and a largest deviation of 0.017; 0.02
    is about four sd."""
    rng = np.random.default_rng(12)
    times, start = [], 0.0
    for length in rng.integers(50, 300, size=40):
        times.append(start + np.arange(length) * 10.0)
        start = times[-1][-1] + 5000.0
    t = np.concatenate(times)
    row = tl.count_exponent_row(
        observed_reads(t, np.ones(t.size)), rng.random(t.size) < 0.3
    )
    assert abs(row.count_exponent - 1.0) < 0.02


@pytest.mark.unit
@pytest.mark.parametrize(
    "t, ins, reason",
    [
        (
            np.arange(30.0),
            np.ones(30, dtype=bool),
            "no down-crossing at this threshold",
        ),
        (
            np.r_[0.0, 1.0, 1000.0, 1001.0],
            np.array([True, False, True, True]),
            "fewer than two rungs with a down-crossing",
        ),
    ],
    ids=["always-in", "one-rung"],
)
def test_an_unanswerable_count_exponent_names_its_reason(t, ins, reason) -> None:
    """Oracle: specification. Two-read segments hold pairs at lag 1 only."""
    with pytest.raises(tl.NotComputed) as caught:
        tl.count_exponent_row(observed_reads(t, np.ones(t.size)), ins)
    assert caught.value.reason == reason


@pytest.mark.unit
def test_rungs_with_no_crossing_or_no_pair_leave_the_fit() -> None:
    """Oracle: hand counts and an in-test least-squares fit. Three segments of twelve reads
    alternating in, out, in, ...: a lag-m pair goes from in to out only for odd m, from an
    even position i with i + m <= 11, so per segment 6, 5 and 4 pairs at m = 1, 3, 5, and
    none at m = 2 or 8, which keep pairs; m = 13 and 20 keep none. With P_avg(m) =
    3 (12 - m) / m, n_avg / P_avg is 18 / 33, 5 / 9 and 2.4 / 4.2 on the fitted rungs."""
    t = np.concatenate([s * 1000.0 + np.arange(12) * 10.0 for s in range(3)])
    ins = np.tile(np.arange(12) % 2 == 0, 3)
    row = tl.count_exponent_row(observed_reads(t, np.ones(t.size)), ins)
    assert row.rungs_excluded == (2, 8, 13, 20)
    assert row.count_ratio[1] == 0.0 and row.count_ratio[4] == 0.0
    assert row.count_ratio[5] is None and row.count_ratio[6] is None
    x = np.log([1.0, 3.0, 5.0])
    y = np.log([18.0 / 33.0, 5.0 / 9.0, 2.4 / 4.2])
    slope = float(np.polyfit(x, y, 1)[0])
    assert row.count_exponent == pytest.approx(1.0 - slope, rel=1e-12)


@pytest.mark.unit
def test_the_ladder_is_l0s() -> None:
    """Oracle: Sera's L0 specification, section 4.2: fitted over m = 1, 2, 3, 5, 8, 13, 20."""
    assert tl.DECIMATION_LADDER == (1, 2, 3, 5, 8, 13, 20)


@pytest.mark.unit
@pytest.mark.parametrize("ladder", [(2, 3, 5), (1, 3, 2), (1, 1, 2), ()])
def test_a_ladder_must_start_at_one_and_increase(ladder) -> None:
    """Oracle: specification. The base rung and every ratio assume rung 1 first."""
    with pytest.raises(ValueError, match="start at 1"):
        tl.count_exponent_row(
            observed_reads(np.arange(30.0), np.ones(30)),
            np.ones(30, dtype=bool),
            ladder,
        )


# --------------------------------------------------------------------------- Q3.4


def _bracket_record():
    """Six segments, 10 s apart inside and 1000 s gaps between, threshold 0, sigma 1 except
    one read whose sigma is unknown. The classes at k = 1 and k = 2, by hand:

    A [3, 1.5, -0.5, -3]: one crossing; its nearest confident reads are 1.5 (k = 1) or 3
       (k = 2), in, and -3, out, with one crossing between: CLEAN at both.
    B [3, 0.5, -0.5, 0.5, -0.5, -3]: two crossings in the bracket 3 .. -3 with three
       crossings in all: MERGED, twice, at both.
    C [1.5, -0.5, -1.5]: CLEAN at k = 1; no confident read at k = 2: UNBRACKETED.
    E [2.5, 3] then D [0.5, -3]: D's crossing has no confident read before it in its own
       segment; E's confident reads are across the gap: UNBRACKETED at both.
    F [3, -0.2 (sigma unknown), 0.3, -3]: the unknown-sigma read is never confident, so both
       crossings sit in the bracket 3 .. -3 with three crossings: MERGED, twice, at both.
    """
    segments = [
        [3, 1.5, -0.5, -3],
        [3, 0.5, -0.5, 0.5, -0.5, -3],
        [1.5, -0.5, -1.5],
        [2.5, 3],
        [0.5, -3],
        [3, -0.2, 0.3, -3],
    ]
    t, v, start = [], [], 0.0
    for seg in segments:
        t.append(start + np.arange(len(seg)) * 10.0)
        v.append(np.asarray(seg, dtype=float))
        start = t[-1][-1] + 1000.0
    v = np.concatenate(v)
    sigma = np.ones(v.size)
    sigma[sum(len(s) for s in segments[:5]) + 1] = np.nan
    return np.concatenate(t), v, sigma


@pytest.mark.unit
@pytest.mark.parametrize("k, expected", [(1, (2, 4, 1)), (2, (1, 4, 2))])
def test_crossings_are_classed_by_their_confident_bracket(k, expected) -> None:
    """Oracle: the hand classification in `_bracket_record`, (clean, merged, unbracketed)."""
    t, v, sigma = _bracket_record()
    counts = tl.crossing_counts(observed_reads(t, v), v, sigma, k)
    assert (counts.clean, counts.merged, counts.unbracketed) == expected


@pytest.mark.unit
@pytest.mark.parametrize("k", [1, 2])
def test_brackets_on_the_wrong_side_and_at_the_threshold(k) -> None:
    """Oracle: hand classification, threshold 0, sigma 1, the same at k = 1 and 2.

    G [-3, 0.5, -0.5, 3]: the bracket runs from -3 (out) to 3 (in): UNBRACKETED.
    H [3, 0.5, -0.5, 3]: the bracket stays in, 3 to 3: UNBRACKETED.
    I [3, 1.5, -0.5, -3, 3]: the bracket ends at -3, one crossing inside it; the up-crossing
       after -3 is outside: CLEAN.
    J [0, -3]: a read at the threshold is in spec, so 0 -> -3 is a crossing, and nothing
       confident comes before it: UNBRACKETED. Read as out of spec it would be no crossing.
    """
    segments = [[-3, 0.5, -0.5, 3], [3, 0.5, -0.5, 3], [3, 1.5, -0.5, -3, 3], [0.0, -3]]
    t = np.concatenate(
        [s * 1000.0 + np.arange(len(seg)) * 10.0 for s, seg in enumerate(segments)]
    )
    v = np.concatenate([np.asarray(seg, dtype=float) for seg in segments])
    counts = tl.crossing_counts(observed_reads(t, v), v, np.ones(v.size), k)
    assert (counts.clean, counts.merged, counts.unbracketed) == (1, 0, 3)


# --------------------------------------------------------------------------- Q4.1


def _windows(t, v, big_values_good=True):
    from quebra.analyzers import windows

    return windows.run(
        windows.WindowsInputs(
            t_rel_s=t, values=v, thresholds=[("u", 0.0, big_values_good)]
        )
    ).windows_in_spec


@pytest.mark.unit
def test_the_mean_run_by_hand_in_both_forms() -> None:
    """Oracle: hand counts. Segment one is in, in, out, in, in, in, out, out; after a gap,
    segment two is in, in, in. Reads with an observed successor in their segment: 2 + 3 + 2
    = 7; in-to-out transitions inside a segment: 2. Mean run 3.5; A = 8 / 11; ratio 3.5 x 3 /
    11. The same mean from the carve's table, sum(n_reads - censored) / sum(1 - censored)."""
    ins = np.array([1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1], dtype=bool)
    t = np.r_[np.arange(8) * 10.0, 1000.0 + np.arange(3) * 10.0]
    v = np.where(ins, 1.0, -1.0)
    row = tl.run_length_row(_windows(t, v), ins)
    assert row.mean_in_spec_run_reads == 3.5
    assert row.read_fraction_in_spec == pytest.approx(8 / 11, rel=1e-15)
    assert row.run_length_ratio == pytest.approx(3.5 * 3 / 11, rel=1e-15)


def _stationary_chain(rng, alpha=0.3, beta=0.1, length=20, n_segments=2000):
    """A two-state chain, P(out -> in) = alpha and P(in -> out) = beta, restarted in every
    segment from its stationary law P(in) = alpha / (alpha + beta)."""
    out = []
    for _ in range(n_segments):
        x = np.empty(length, dtype=bool)
        x[0] = rng.random() < alpha / (alpha + beta)
        u = rng.random(length)
        for i in range(1, length):
            x[i] = (u[i] >= beta) if x[i - 1] else (u[i] < alpha)
        out.append(x)
    t = np.concatenate(
        [s * 1000.0 + np.arange(length) * 10.0 for s in range(n_segments)]
    )
    return t, np.concatenate(out)


@pytest.mark.statistical
def test_the_mean_run_of_a_chain_is_one_over_beta() -> None:
    """Oracle: analytic, as a ratio of expectations. Each in-spec read with a successor
    leaves with probability beta, so E[reads with successor] / E[exits] = 1 / beta = 10,
    whatever the segment length. An inverted margin reads the out-of-spec runs, 1 / alpha
    = 3.3. Over seeds 0-39, each drawing the chain from `default_rng(seed)`, the estimate
    measured 10.005 with sd 0.19 and a largest deviation of 0.40; 0.8 is about four sd."""
    t, x = _stationary_chain(np.random.default_rng(13))
    row = tl.run_length_row(_windows(t, np.where(x, 1.0, -1.0)), x)
    assert abs(row.mean_in_spec_run_reads - 10.0) < 0.8


# --------------------------------------------------------------------------- bounds


def _bounds_record():
    """Spacings 10, 20, 30, 40 s, a gap, then 10 and 20 s; threshold 0, sigma 1."""
    t = np.array([0.0, 10.0, 30.0, 60.0, 100.0, 1000.0, 1010.0, 1030.0])
    v = np.array([3.0, 0.5, -2.0, -3.0, 2.0, 1.5, 2.0, -2.0])
    return t, v, np.ones(t.size)


@pytest.mark.unit
@pytest.mark.parametrize(
    "k, lower, upper",
    [(1.0, 20.0 / 130.0, 100.0 / 130.0), (2.0, 10.0 / 130.0, 100.0 / 130.0)],
)
def test_the_in_spec_fraction_bounds_by_hand(k, lower, upper) -> None:
    """Oracle: hand arithmetic, each interval held by its left read. At k = 1: [0, 10) is
    settled in (3 -> 0.5, its RIGHT read unsure, which does not matter); [10, 30) is unsure
    (0.5); [30, 60) settled out; [60, 100) holds a crossing; the gap is no time; [1000, 1010)
    settled in; [1010, 1030) starts sure in (2) and holds a crossing, so it unsettles.
    Observed 130 s, settled in 20, unsettled 80. At k = 2 the 1.5 read is unsure too, so
    [1000, 1010) unsettles; 2 is still sure, since |margin| = 2 sigma counts."""
    t, v, sigma = _bounds_record()
    b = tl.in_spec_fraction_bounds(observed_reads(t, v), v, sigma, k)
    assert b.in_spec_fraction_lower == pytest.approx(lower, rel=1e-15)
    assert b.in_spec_fraction_upper == pytest.approx(upper, rel=1e-15)
    assert b.unsure_k == k


@pytest.mark.unit
def test_a_read_with_no_sigma_unsettles_the_interval_it_holds() -> None:
    """Oracle: hand arithmetic on `_bounds_record` with the first read's sigma unknown: its
    interval [0, 10), settled in at k = 1, unsettles, so the lower bound loses 10 s and the
    upper keeps them: 10 / 130 and 100 / 130."""
    t, v, sigma = _bounds_record()
    sigma[0] = np.nan
    b = tl.in_spec_fraction_bounds(observed_reads(t, v), v, sigma, 1.0)
    assert b.in_spec_fraction_lower == pytest.approx(10.0 / 130.0, rel=1e-15)
    assert b.in_spec_fraction_upper == pytest.approx(100.0 / 130.0, rel=1e-15)


@pytest.mark.unit
@pytest.mark.parametrize("record", ["bounds", "brackets"])
def test_the_carve_fraction_lies_between_the_bounds(record) -> None:
    """Oracle: the carve's own time-weighted in-spec fraction, its in-spec durations over
    the observed time, which counts every interval by its left read."""
    t, v, sigma = _bounds_record() if record == "bounds" else _bracket_record()
    obs = observed_reads(t, v)
    observed = float(np.sum(np.diff(obs.t_rel_s)[obs.within_segment()]))
    point = float(_windows(t, v)["duration_s"].sum()) / observed
    for k in (1.0, 2.0):
        b = tl.in_spec_fraction_bounds(obs, v, sigma, k)
        assert b.in_spec_fraction_lower <= point <= b.in_spec_fraction_upper


# --------------------------------------------------------------------------- rho_k


@pytest.mark.statistical
def test_rho_k_of_a_stationary_chain_decays_as_a_power() -> None:
    """Oracle: analytic. A stationary two-state chain has rho_k = (1 - alpha - beta)^k, here
    0.6^k. Restarting every 5-read segment from the stationary law keeps it stationary
    inside a segment, and the array stops at lag 4. A pair across a gap pairs independent
    segments, which short segments make frequent: it moves rho_1 to about 0.48. Over seeds
    0-39, each drawing the chain from `default_rng(seed)`, the deviation from 0.6^k stayed
    under 0.012, 0.020 and 0.029 at lags 1-3 (sd up to 0.011), and the cross-gap version's
    above 0.11 at every lag; 0.045 is about four sd at lag 3."""
    t, x = _stationary_chain(np.random.default_rng(14), length=5, n_segments=8000)
    rho = tl.in_spec_autocorr(observed_reads(t, np.ones(t.size)), x)
    assert len(rho) == 4
    for lag in (1, 2, 3):
        assert abs(rho[lag - 1] - 0.6**lag) < 0.045


@pytest.mark.unit
def test_rho_k_by_hand_on_two_segments() -> None:
    """Oracle: hand arithmetic. Segments [1, 1, 0] and, after a gap, [0, 0, 0, 1]: A = 3/7,
    A (1 - A) = 12/49. Centred values 4/7 or -3/7; the within-segment products sum to 10/49
    over 5 pairs at lag 1, -15/49 over 3 at lag 2, -12/49 over 1 at lag 3; no segment holds a
    lag-4 pair. So rho = (1/6, -5/12, -1): one record-wide centre and divisor, not one per
    lag, and the array stops at lag 3."""
    t = np.r_[np.arange(3) * 10.0, 1000.0 + np.arange(4) * 10.0]
    x = np.array([1, 1, 0, 0, 0, 0, 1], dtype=bool)
    rho = tl.in_spec_autocorr(observed_reads(t, np.ones(t.size)), x)
    np.testing.assert_allclose(rho, [1.0 / 6.0, -5.0 / 12.0, -1.0], rtol=1e-12)


@pytest.mark.unit
def test_the_run_length_ratio_and_rho_1_agree_on_a_chain() -> None:
    """CONSISTENCY check, not evidence: for a stationary 0/1 sequence the run-length ratio is
    1 / (1 - rho_1), an identity of the two estimators' targets; both descend from the same
    reads, so this shows they agree, not that either is right."""
    t, x = _stationary_chain(np.random.default_rng(15))
    row = tl.run_length_row(_windows(t, np.where(x, 1.0, -1.0)), x)
    rho_1 = tl.in_spec_autocorr(observed_reads(t, np.ones(t.size)), x)[0]
    assert row.run_length_ratio == pytest.approx(1.0 / (1.0 - rho_1), rel=0.05)


@pytest.mark.unit
@pytest.mark.parametrize(
    "call, reason",
    [
        (
            lambda: tl.in_spec_autocorr(
                observed_reads(np.arange(5.0), np.ones(5)), np.ones(5, dtype=bool)
            ),
            "zero variance",
        ),
        (
            lambda: tl.run_length_row(
                _windows(np.arange(4.0), np.ones(4)), np.ones(4, dtype=bool)
            ),
            "no in-to-out transition within a segment",
        ),
        (
            lambda: tl.in_spec_autocorr(
                observed_reads(np.arange(4.0) * 10.0, np.ones(4), gap_mult=0.5),
                np.array([1, 0, 1, 0], dtype=bool),
            ),
            "no segment with at least 2 reads",
        ),
        (
            lambda: tl.crossing_counts(
                observed_reads(np.arange(4.0), np.ones(4)),
                np.array([1.0, -1.0, 1.0, -1.0]),
                np.full(4, np.nan),
                1,
            ),
            "no read with a finite sigma",
        ),
        (
            lambda: tl.in_spec_fraction_bounds(
                observed_reads(np.arange(4.0), np.ones(4)),
                np.array([1.0, -1.0, 1.0, -1.0]),
                np.full(4, np.nan),
                1.0,
            ),
            "no read with a finite sigma",
        ),
    ],
    ids=[
        "rho-constant",
        "no-exit",
        "rho-one-read-segments",
        "brackets-no-sigma",
        "bounds-no-sigma",
    ],
)
def test_an_unanswerable_threshold_row_names_its_reason(call, reason) -> None:
    """Oracle: specification."""
    with pytest.raises(tl.NotComputed) as caught:
        call()
    assert caught.value.reason == reason
