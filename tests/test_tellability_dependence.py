"""Tellability's sigma rows against known answers: the dependence panel (Q3.0), the
self-consistency of sigma (Q3.2) and sigma against its data-only bound (Q3.3).

Oracles:
- the size of a valid permutation test, computed exactly per replicate by an in-test
  reference that recounts every row's member p by brute force;
- population values: xi(sigma | q) = 1 and xi(q | sigma) = 1/4 for q uniform on [-1, 1] and
  sigma = q^2; sigma = -q is the most extreme pairing for |r| on every row;
- a reference implementation, `dcor`'s naive distance correlation;
- hand arithmetic for the decile spread;
- simulation truth: s^2 / (s^2 + mean sigma^2) for a white signal plus noise of known sigma,
  and a sigma_bar / sigma_hat_D2 of 1 when the reported sigma is the noise's.
"""

from __future__ import annotations

import math

import dcor
import numpy as np
import pytest

from quebra.analyzers import tellability as tl
from quebra.analyzers.min_p import min_p
from quebra.analyzers.windows import observed_reads


def _reference_row_p(table: dict[str, np.ndarray]) -> np.ndarray:
    """Every row's min-p by brute force: #{rows at least as large} / (B + 1) per member."""
    columns = np.column_stack(list(table.values()))
    at_least = columns[None, :, :] >= columns[:, None, :]  # [b, b', j]: row b' >= row b
    return at_least.mean(axis=1).min(axis=1)


@pytest.mark.statistical
@pytest.mark.slow
def test_the_panel_holds_its_size_under_independence() -> None:
    """Oracle: the exact conditional size, against which the rejection rate is checked.

    Under independence every row is equally likely to be the observed one, so the chance
    that the observed row is rejected, given the B + 1 rows, is #{rows b whose own panel p
    would be <= k / (B + 1)} / (B + 1): the in-test reference recounts it by brute force.
    With R = 2000 replicates, B = 99 and alpha = 0.05 (k = 5): (a) the rejection rate is at
    most k / (B + 1) plus 3 binomial SE, and (b) it lies within 3 binomial SE of the mean
    exact conditional size. Ties make the test conservative, so (a) is a bound, not a
    target.
    """
    rng = np.random.default_rng(20261013)
    n_rep, n_perm, n = 2000, 99, 60
    k = math.floor(0.05 * (n_perm + 1))
    rejected = 0
    conditional = []
    for _ in range(n_rep):
        q, s = rng.normal(size=n), rng.normal(size=n)
        perms = np.vstack([rng.permutation(n) for _ in range(n_perm)])
        table = tl.dependence_statistics(q, s, perms)
        rejected += min_p(table).p_value <= k / (n_perm + 1)
        row_min = _reference_row_p(table)
        counts = (row_min[None, :] <= row_min[:, None]).sum(axis=1)
        conditional.append(np.mean(counts <= k))
    rate = rejected / n_rep
    exact = float(np.mean(conditional))
    assert rate <= k / (n_perm + 1) + 3.0 * math.sqrt(0.05 * 0.95 / n_rep)
    assert abs(rate - exact) <= 3.0 * math.sqrt(exact * (1.0 - exact) / n_rep)


@pytest.mark.unit
def test_xi_reads_sigma_as_a_function_of_q() -> None:
    """Oracle: population values. With sigma = q^2, xi(sigma | q) is 1 and xi(q | sigma) is
    1/4, so a member computed in the wrong direction reads near 0.25, not near 1."""
    rng = np.random.default_rng(1)
    q = rng.uniform(-1.0, 1.0, size=2000)
    panel = tl.dependence_panel(q, q**2, n_permutations=19, seed=2)
    assert panel.statistics["chatterjee_xi"] > 0.95


@pytest.mark.unit
def test_a_perfect_negative_pairing_is_the_most_extreme_for_r() -> None:
    """Oracle: |r| is largest when sigma = -q, so on no shuffled row is it as large: the
    Pearson and Spearman member p are the smallest the grid allows, 1 / (B + 1). A signed,
    one-sided r would put them near 1."""
    rng = np.random.default_rng(3)
    q = rng.normal(size=200)
    panel = tl.dependence_panel(
        q, -q + 1e-3 * rng.normal(size=200), n_permutations=19, seed=4
    )
    assert panel.member_p_values["pearson_abs"] == 1 / 20
    assert panel.member_p_values["spearman_abs"] == 1 / 20


@pytest.mark.unit
def test_the_distance_correlation_matches_the_naive_reference() -> None:
    """Oracle: `dcor`'s naive O(n^2) algorithm, independent of the merge-sort one used."""
    rng = np.random.default_rng(5)
    q = rng.normal(size=300)
    s = np.abs(q) + rng.normal(scale=0.5, size=300)
    panel = tl.dependence_panel(q, s, n_permutations=1, seed=6)
    expected = dcor.distance_correlation(
        q, s, method=dcor.DistanceCovarianceMethod.NAIVE
    )
    assert panel.statistics["distance_correlation"] == pytest.approx(
        float(expected), rel=1e-9
    )


@pytest.mark.unit
def test_the_decile_spread_by_hand() -> None:
    """Oracle: hand arithmetic. q = 1..100 and sigma 1 below q = 51, 3 from it: the lower
    five deciles average 1, the upper five 3, the overall mean is 2, so (3 - 1) / 2 = 1."""
    q = np.arange(1.0, 101.0)
    panel = tl.dependence_panel(
        q, np.where(q <= 50, 1.0, 3.0), n_permutations=1, seed=7
    )
    assert panel.statistics["decile_spread"] == pytest.approx(1.0, rel=1e-12)


@pytest.mark.statistical
def test_self_consistency_recovers_the_signal_share_of_the_variance() -> None:
    """Oracle: simulation truth. A white signal of variance s^2 plus noise of known, varying
    sigma_i (uniform on 0.2-1.8, CV about 0.46) has (v - mean sigma^2) / v -> s^2 / (s^2 +
    mean sigma^2), here 0.5. Using mean(sigma)^2 instead reads about 0.59. Over 200 seeds
    the estimate measured sd 0.006 and a largest deviation of 0.018; 0.03 is five sd."""
    rng = np.random.default_rng(9)
    n = 20000
    sigma = rng.uniform(0.2, 1.8, size=n)
    signal_var = 1.2133
    q = rng.normal(scale=math.sqrt(signal_var), size=n) + rng.normal(scale=sigma)
    truth = signal_var / (signal_var + float(np.mean(sigma**2)))
    assert abs(tl.self_consistency(q, sigma) - truth) < 0.03


@pytest.mark.statistical
def test_sigma_meets_its_bound_when_it_is_the_noise() -> None:
    """Oracle: simulation truth. A linear trend plus noise of sd 0.5, reported as 0.5: the
    second-difference scale recovers 0.5 (within 5%, see the Q2.1 test), so the ratio is 1."""
    rng = np.random.default_rng(10)
    n = 6000
    t = np.r_[0.0, np.cumsum(rng.uniform(10.0, 40.0, size=n - 1))]
    obs = observed_reads(t, 0.05 * t + rng.normal(scale=0.5, size=n))
    assert tl.sigma_bound_ratio(obs, np.full(n, 0.5)) == pytest.approx(1.0, rel=0.05)


@pytest.mark.unit
@pytest.mark.parametrize(
    "values, sigma, reason",
    [
        (np.arange(20.0), np.full(20, np.nan), "no read with a finite sigma"),
        (np.arange(20.0), np.ones(20), "zero variance"),
        (np.arange(9.0), np.arange(9.0), "fewer than 10 reads"),
    ],
    ids=["no-sigma", "constant-sigma", "too-few"],
)
def test_an_unanswerable_panel_names_its_reason(values, sigma, reason) -> None:
    """Oracle: specification."""
    with pytest.raises(tl.NotComputed) as caught:
        tl.dependence_panel(values, sigma, n_permutations=9, seed=0)
    assert caught.value.reason == reason
