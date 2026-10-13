"""`tlf`'s dwell statistics against the SAMPLED two-state law, not its continuous limit.

Oracle: a two-state chain observed on a read grid. Sampled at spacing `dt`, a state left
with per-read probability `p` has run lengths exactly Geometric(p) on {1, 2, ...}. `tlf.run`
measures dwells on the ASSIGNED READ SEQUENCE: a state's mean dwell is its exposure (the
forward spacings whose left read is in the state) over its exits (observed changes out of
it), and its CV is over complete runs. On a uniform grid a complete run of L reads lasts
L * dt, so the quantities it reports converge to

    mean dwell -> dt / p            dwell CV -> sqrt(1 - p)

and NOT to the continuous-time `1/lambda` and `CV = 1` of the two-state chain this sampled
one embeds, whose rates satisfy p01 + p10 = 1 - exp(-(lambda0 + lambda1) dt) and split it in
the ratio p01 : p10. They are measurably different here: state 0's sampled mean is 6.67 read
intervals against a continuous 4.37, and its sampled CV is 0.922 against 1.0. An earlier
draft of `spec/spectests06.md` asserted the continuous values; it was wrong, and the
assertions below are written so that version of the test would fail.

**The two rates are deliberately UNEQUAL, and an earlier version of this file was symmetric.**
Under `lambda = mu` the two states have identical dwell laws, so every assertion held whichever
label the mixture happened to assign - inverting `tlf.py`'s label map left the ENTIRE suite at
585 passed, measured. A test that cannot fail when the map it names is inverted is the defect
class this phase exists to remove. With `p01 != p10` the states carry different means (6.67
against 2.22 read intervals) and a swapped map fails immediately.

`switching_rate_per_hour` still has a closed form under asymmetry: the stationary transition
rate is `2 * p01 * p10 / (p01 + p10)` per read. It is asserted, and labelled a consistency
check rather than independent evidence, because it is a function of the same two parameters.

**Restriction, stated rather than hidden.** The draw is noise-free two-level, because
`tlf.run` assigns state by a per-read MAP on the value alone: one misassigned read splits a
true run into three, so the oracle needs misassignment far below `p01`. With measurement
noise there is no closed form and this stops being simulation truth. In that configuration
the Gaussian mixture degenerates to a threshold and does no work, so what this test validates
is the run-length bookkeeping, the label-ordering map and the CV formula - a unit test of
`tlf.py`'s dwell path, not a validation of mixture-based state assignment. Said plainly so
the tier it earns is not overstated.
"""

from __future__ import annotations

import numpy as np
import pytest

from quebra.analyzers import tlf

pytestmark = pytest.mark.statistical

DT_S = 1.0
# UNEQUAL by design: equal rates make the two states indistinguishable and the ordering
# assertion unfalsifiable. Chosen so the sampled and continuous laws are also far apart.
P01 = 0.15  # per-read probability of leaving state 0 (the LOW level) -> longer dwell
P10 = 0.45  # per-read probability of leaving state 1 (the HIGH level) -> shorter dwell
N_READS = 200_000
SEED = 4242
LOW, HIGH = 0.0, 10.0
# Set from the seed-to-seed SPREAD of the statistic, not from how close one draw happens to
# land. Measured over 2000 redraws of the same law: the relative error has sd 0.49% to 0.69%
# depending on the statistic, and a maximum observed 2.40%. A 2% tolerance is therefore only
# about 3 sd and would go red on roughly 1 redraw in 110 - a failure carrying no information
# about the code. 5% is about 7 sd and costs nothing: the discriminating assertions are the
# sampled-versus-continuous comparisons below, whose margins are several times this tolerance.
REL_TOL = 0.05


def _two_state_record() -> tuple[np.ndarray, np.ndarray]:
    """An ASYMMETRIC sampled two-state chain, noise-free, as a metric record."""
    rng = np.random.default_rng(SEED)
    draw = rng.random(N_READS)
    state = np.empty(N_READS, dtype=int)
    state[0] = 0
    for i in range(1, N_READS):
        leaving = P01 if state[i - 1] == 0 else P10
        state[i] = (1 - state[i - 1]) if draw[i] < leaving else state[i - 1]
    values = np.where(state == 1, HIGH, LOW)
    return values, np.arange(N_READS, dtype=float) * DT_S


@pytest.fixture(scope="module")
def result():
    values, timestamps = _two_state_record()
    return tlf.run(values, timestamps, seed=7)


def test_mean_dwell_recovers_the_sampled_law_not_the_continuous_one(result) -> None:
    """Oracle: `dt / p`, the mean of Geometric(p) in read intervals, per state."""
    for state, measured, p in (
        ("s0", result.mean_dwell_state0_s, P01),
        ("s1", result.mean_dwell_state1_s, P10),
    ):
        sampled = DT_S / p
        # the embedded continuous-time chain: the two rates share -ln(1 - p01 - p10) / dt in
        # the ratio of the per-read probabilities
        lam = p / (P01 + P10) * -np.log(1.0 - P01 - P10) / DT_S
        continuous = 1.0 / lam

        assert measured == pytest.approx(sampled, rel=REL_TOL), (
            f"{state}: mean dwell {measured:.4f} against the sampled law {sampled:.4f}"
        )
        assert abs(measured - sampled) < abs(measured - continuous), (
            f"{state}: mean dwell {measured:.4f} must be nearer the SAMPLED law "
            f"{sampled:.4f} than the continuous limit {continuous:.4f}"
        )

    # The states must not be interchangeable: a swapped label map fails here.
    assert result.mean_dwell_state0_s > result.mean_dwell_state1_s * 2.0, (
        "state 0 leaves with probability p01 < p10, so it must dwell markedly longer; "
        "if these are close the label map has been inverted or the draw is symmetric"
    )


def test_dwell_cv_recovers_the_sampled_law_not_unity(result) -> None:
    """Oracle: `sqrt(1 - p)`, the CV of Geometric(p), per state.

    The continuous limit is 1.0 - the CV of an exponential - and is what an implementation
    measuring true sojourn times would give. Here the pair is (0.922, 0.742) against 1.0.
    """
    for state, measured, p in (
        ("s0", result.dwell_cv_state0, P01),
        ("s1", result.dwell_cv_state1, P10),
    ):
        sampled = float(np.sqrt(1.0 - p))
        assert measured == pytest.approx(sampled, rel=REL_TOL), (
            f"{state}: dwell CV {measured:.4f} against the sampled law {sampled:.4f}"
        )
        assert abs(measured - sampled) < abs(measured - 1.0), (
            f"{state}: dwell CV {measured:.4f} must be nearer sqrt(1-p) {sampled:.4f} "
            f"than the exponential's 1.0"
        )


def test_the_switching_rate_agrees_with_the_dwell_statistics(result) -> None:
    """Consistency check, not independent evidence.

    For an asymmetric chain the stationary occupancies are `p10/(p01+p10)` and
    `p01/(p01+p10)`, so transitions per read are `2*p01*p10/(p01+p10)`. Asserted so a change
    that broke the transition count without touching the dwells is visible, and labelled so
    it is not counted as a third oracle: it is a function of the same two parameters.
    """
    stationary_per_read = 2.0 * P01 * P10 / (P01 + P10)
    expected_per_hour = stationary_per_read / DT_S * 3600.0
    assert result.switching_rate_per_hour == pytest.approx(
        expected_per_hour, rel=REL_TOL
    )


def test_the_two_states_are_ordered_low_then_high(result) -> None:
    """Oracle: `tlf.py`'s remap, which sorts components so state 0 is the LOWER mean.

    Without it `s0` and `s1` would swap between runs on the mixture's fitting order, and
    every dwell assertion above would be comparing an arbitrary label to a fixed one. The
    ordering is not exposed as a field, so it is read from the fitted mixture the result
    carries.
    """
    means = np.sort(np.asarray(result.gmm2.means_).ravel())
    assert means[0] == pytest.approx(LOW, abs=1e-6)
    assert means[1] == pytest.approx(HIGH, abs=1e-6)

    # The discriminating half. The record is built so the LOW level is the long-dwell state,
    # so `s0` must carry the longer mean. An earlier symmetric version of this fixture made
    # the two dwell laws identical, and inverting `tlf.py`'s label map then left the whole
    # suite green - the assertion below is what makes that mutation fail.
    assert result.mean_dwell_state0_s == pytest.approx(DT_S / P01, rel=REL_TOL)
    assert result.mean_dwell_state1_s == pytest.approx(DT_S / P10, rel=REL_TOL)
