"""Tier 2: reproduce the numbers Kvaloy and Lindqvist print, from their own data.

This is the tier the suite was missing, and it is the one that catches transcription
errors nothing else catches. The bench measures how a routine BEHAVES - a mis-transcribed
equation that behaves plausibly passes every size and power test we have. Only a published
worked example pins the arithmetic itself.

Source: Kvaloy and Lindqvist, Technometrics 62(1):101-115 (2020). The small bowel record
is in preprint arXiv:1802.08339v1 only (its Section 6.2). Three records, all in
`tests/fixtures/load_haul_dump.py`:

- the single load-haul-dump machine, time censored at 2000 hours: Technometrics Section
  8.1, Tables 1-2 and the Section 8.1 text, p. 112 (arXiv v1 Section 6.1, p. 14);
- the hydraulic systems of six machines, m = 6: Technometrics Section 8.2, Table 4, p. 113.
  This is the external check on the m > 1 path;
- the small bowel motility record, m = 19: arXiv v1 Section 6.2, pp. 15-17, with the data
  from Aalen and Husebye (1991) Table I.

**What is a direct pin and what is a composition.** `gamma_hat`, `gamma_tilde` and `LR`
are returned by shipped functions. `mu_hat`, `sigma_tilde` and the Laplace statistic are
not exposed by the API and are reconstructed here from shipped outputs:

    sigma_tilde = gamma_tilde * tau / N        (definition of gamma_tilde)
    laplace     = LR * gamma_hat               (LR is the Laplace statistic / gamma_hat)

Those two are still genuine checks - they assert that shipped outputs COMPOSE to an
independently published number - but they are weaker than a direct pin and this docstring
says so rather than letting a reader assume otherwise.

**The divisor, stated as what is actually sourced.** Section 3.6 (p. 104) defines
gamma_hat from the "sample mean" and "sample standard deviation" of the complete gaps.
Appendix A.2 (p. 114; arXiv v1 Appendix 1, p. 18) writes that estimator with 1/N(tau), in
a consistency argument where the divisor is asymptotically irrelevant. Every number the
paper prints uses 1/(N-1): Table 2's 48.61 and 0.888, the text's 0.681, and Table 4's
LR p = 0.019 (which 1/N moves to 0.016). The authors' R code computes `sqrt(var(x))`.

`GAMMA_COMPLETE` is the population form, `GAMMA_COMPLETE_SAMPLE` the sample form. The pins
below name which one they use. Consequence of the population form: gamma_hat is smaller
by sqrt((N-1)/N), so the C1 statistic is LARGER by sqrt(N/(N-1)) - 1.4% at N = 36, about
12% on a five-event segment. That direction is anti-conservative.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

import quebra.analyzers.checks.c1_lewis_robinson as c1
from quebra.analyzers.checks._multiprocess import (
    GAMMA_COMPLETE,
    GAMMA_COMPLETE_SAMPLE,
    GAMMA_TRUNCATED,
    gamma_hat,
)
from quebra.analyzers.checks.result import Segment, validate_segment
from tests.fixtures.load_haul_dump import (
    HYDRAULIC_PUBLISHED,
    PUBLISHED,
    SMALL_BOWEL_PUBLISHED,
    SMALL_BOWEL_TABLE_I,
    TAU_H,
    gaps_h,
    hydraulic_segments_data,
)

pytestmark = pytest.mark.statistical


X = gaps_h()
N = len(X)
SEGMENT = Segment(x=X, tau=TAU_H, n_censored_dropped=1)


def test_the_published_record_is_a_well_posed_segment():
    """It must survive our own guards before any number it produces means anything."""
    validate_segment(SEGMENT, require_strict_tau=True)
    assert N == 36
    assert np.all(X > 0.0)
    assert SEGMENT.tau > float(np.cumsum(X)[-1])


# --------------------------------------------------------------- direct pins


def test_gamma_tilde_eq10_matches_the_paper():
    """Eq (10), the estimator that uses the censored time. Direct pin."""
    assert gamma_hat(X, TAU_H, GAMMA_TRUNCATED) == pytest.approx(
        PUBLISHED["gamma_tilde"], abs=0.001
    )


def test_gamma_hat_complete_gaps_is_the_population_form():
    """Oracle: Appendix A.2 (p. 114) 1/N(tau) form, from Table 2's 48.61 rescaled by sqrt(35/36)."""
    assert gamma_hat(X, TAU_H, GAMMA_COMPLETE) == pytest.approx(0.876, abs=0.001)


def test_gamma_hat_sample_form_reproduces_table_2():
    """Oracle: Table 2 row 1 (p. 112), gamma_hat = 0.888, sample divisor. Direct pin."""
    assert gamma_hat(X, TAU_H, GAMMA_COMPLETE_SAMPLE) == pytest.approx(
        PUBLISHED["gamma_hat"], abs=0.0005
    )


def test_lewis_robinson_sample_form_reproduces_the_section_8_1_value():
    """Oracle: Section 8.1 text (p. 112), LR = 0.605/0.888 = 0.681. Direct pin."""
    assert c1.statistic(
        [SEGMENT], gamma_estimator=GAMMA_COMPLETE_SAMPLE
    ) == pytest.approx(PUBLISHED["lr_gamma_hat"], abs=0.0005)


def test_lewis_robinson_matches_the_paper_up_to_the_divisor():
    """LR with the population gamma_hat. The paper's 0.681 uses the sample divisor."""
    assert c1.statistic([SEGMENT], gamma_estimator=GAMMA_COMPLETE) == pytest.approx(
        0.691, abs=0.001
    )


# ------------------------------------------------------- composed quantities


def test_mu_hat_matches_the_paper():
    """Composition: the mean gap is not returned by any shipped function."""
    assert float(np.mean(X)) == pytest.approx(PUBLISHED["mu_hat"], abs=0.01)


def test_sigma_tilde_eq10_matches_the_paper():
    """Composition: `gamma_tilde * tau / N`, since gamma_hat returns only the ratio."""
    sigma_tilde = gamma_hat(X, TAU_H, GAMMA_TRUNCATED) * TAU_H / N
    assert sigma_tilde == pytest.approx(PUBLISHED["sigma_tilde"], abs=0.01)


def test_the_laplace_statistic_matches_the_paper():
    """Composition: `LR * gamma_hat`.

    Not circular. LR and gamma_hat are computed independently by shipped code; asserting
    their product equals a number printed in the paper tests that the eq (4) numerator
    scaling (eq (16) at m = 1) is right, because the Laplace statistic is exactly that
    numerator before the gamma division (Section 3.1, p. 103). A transcription error in the
    `sqrt(12)/(tau*sqrt(N))` factor would show up here and nowhere else in this file.
    """
    lr = c1.statistic([SEGMENT], gamma_estimator=GAMMA_COMPLETE)
    laplace = lr * gamma_hat(X, TAU_H, GAMMA_COMPLETE)
    assert laplace == pytest.approx(PUBLISHED["laplace"], abs=0.001)


# ------------------------------------------------- the documented divergence


def test_the_divergence_is_the_divisor_and_not_a_transcription_error():
    """The 1/N versus 1/(N-1) bridge, at three levels at once.

    The evidential weight is in the SIMULTANEITY, not in any one line: sigma, gamma and LR
    all land on their published values under the same single rescaling. Any one of them
    alone would be a factor of 1.0142 that, inside a 0.01 band at one N, cannot be
    distinguished from a neighbouring constant.
    """
    mu = float(np.mean(X))
    sigma_population = gamma_hat(X, TAU_H, GAMMA_COMPLETE) * mu

    assert sigma_population == pytest.approx(47.93, abs=0.01)
    assert sigma_population * np.sqrt(N / (N - 1)) == pytest.approx(
        PUBLISHED["sigma_hat"], abs=0.01
    )
    # And the same bridge at the level of gamma and of LR.
    assert sigma_population * np.sqrt(N / (N - 1)) / mu == pytest.approx(
        PUBLISHED["gamma_hat"], abs=0.001
    )
    lr_population = c1.statistic([SEGMENT], gamma_estimator=GAMMA_COMPLETE)
    assert lr_population * np.sqrt((N - 1) / N) == pytest.approx(
        PUBLISHED["lr_gamma_hat"], abs=0.001
    )


# ----------------------------------------------------- eq (11), not shipped


def _sigma_star_squared(x: np.ndarray) -> float:
    """Kvaloy and Lindqvist eq (11), p. 104, the successive-difference estimator.

    DELIBERATELY LOCAL TO THIS TEST. The paper evaluates it and declines to use it in its
    own examples "due to apparent less satisfactory significance level properties"
    (Section 3.6, p. 104). It tends to be smaller than sigma_hat under positive dependence
    between neighbouring gaps, which inflates the trend statistic - and neighbour
    dependence is exactly what C5 and C6 exist to detect, so shipping it would entangle the
    trend test with the dependence tests. `_multiprocess.GAMMA_ESTIMATORS` therefore offers
    the complete-gap estimator (in two divisors) and eq (10) only. It is implemented here
    purely to reproduce the published number and confirm the transcription of the
    surrounding machinery.
    """
    return float(np.sum(np.diff(x) ** 2) / (2.0 * (len(x) - 1)))


def test_eq11_reproduces_the_papers_alternative_lr():
    """Oracle: Section 8.1 text (p. 112), 0.605/0.782 = 0.774; the value is not in Table 2."""
    mu = float(np.mean(X))
    gamma_star = np.sqrt(_sigma_star_squared(X)) / mu
    laplace = c1.statistic([SEGMENT], gamma_estimator=GAMMA_COMPLETE) * gamma_hat(
        X, TAU_H, GAMMA_COMPLETE
    )
    assert laplace / gamma_star == pytest.approx(PUBLISHED["lr_sigma_star"], abs=0.001)


def test_eq11_is_not_reachable_through_the_shipped_api():
    """Pins the decision: asking for it must fail, not silently pick another estimator."""
    from quebra.analyzers.checks._multiprocess import GAMMA_ESTIMATORS

    assert "eq11_successive_difference" not in GAMMA_ESTIMATORS
    assert set(GAMMA_ESTIMATORS) == {
        GAMMA_COMPLETE,
        GAMMA_COMPLETE_SAMPLE,
        GAMMA_TRUNCATED,
    }
    with pytest.raises(ValueError, match="unknown gamma estimator"):
        gamma_hat(X, TAU_H, "eq11_successive_difference")


# ------------------------------------------------- m > 1: the hydraulic record


def _hydraulic_segments() -> list[Segment]:
    return [Segment(x=x, tau=tau) for x, tau in hydraulic_segments_data()]


def _two_sided_p(z: float) -> float:
    return float(2.0 * stats.norm.sf(abs(z)))


def test_hydraulic_transcription_reproduces_the_generalised_laplace_p():
    """Oracle: Technometrics Table 4 (p. 113), GL p = 0.062, from the data alone.

    GL = sum_j U_j / sqrt(sum_j U_j^2) with U_j = sum_i T_ij - N_j tau_j / 2 (KL2020
    Section 5.2, eq (17), p. 106; the test is Lawless, Cigsar and Cook's, 2012). GL uses no
    gamma and no weights, so this checks the FIXTURE independently of anything shipped.
    It is a check, not a pin: one p printed to 3 decimals catches most dropped or
    duplicated gaps but no single 1 h error in any of the 152. The gaps themselves were
    read by eye against Kumar and Klefsjo (1992) p. 224 and the authors' R code, and
    `scripts/verify_gold_standard.py` repeats that comparison.
    """
    u = np.array(
        [
            np.cumsum(x).sum() - len(x) * tau / 2.0
            for x, tau in hydraulic_segments_data()
        ]
    )
    gl = float(u.sum() / np.sqrt(np.sum(u**2)))
    assert _two_sided_p(gl) == pytest.approx(HYDRAULIC_PUBLISHED["gl_p"], abs=0.0005)


def test_c1_multiprocess_reproduces_technometrics_table_4():
    """Oracle: Technometrics Table 4 (p. 113), LR p = 0.019, eq (15) weights. Direct pin.

    The external check on the m > 1 path, through the shipped defaults of `c1.run` first
    and then with both choices named. Only the journal weights with the sample divisor
    land on the printed value. The two negative controls have their own oracles: 0.0073
    (preprint weights, sample divisor) is what the authors' R code gives
    (`LRtest_multi`, weights "CVntau", sigma "s": 0.00729); 0.0164 (journal weights, 1/N)
    has no published or R value, since that code has no population-divisor option, and
    is the independent numpy recomputation in `scripts/verify_gold_standard.py`.
    """
    segments = _hydraulic_segments()
    shipped = c1.run(segments)
    assert shipped.p_value == pytest.approx(HYDRAULIC_PUBLISHED["lr_p"], abs=0.0005)
    assert f"weights={c1.WEIGHTS_TECHNOMETRICS}" in shipped.notes.split()

    journal = c1.statistic(
        segments,
        gamma_estimator=GAMMA_COMPLETE_SAMPLE,
        weights=c1.WEIGHTS_TECHNOMETRICS,
    )
    assert _two_sided_p(journal) == pytest.approx(
        HYDRAULIC_PUBLISHED["lr_p"], abs=0.0005
    )

    preprint = c1.statistic(
        segments, gamma_estimator=GAMMA_COMPLETE_SAMPLE, weights=c1.WEIGHTS_ARXIV_V1
    )
    population = c1.statistic(
        segments, gamma_estimator=GAMMA_COMPLETE, weights=c1.WEIGHTS_TECHNOMETRICS
    )
    assert _two_sided_p(preprint) == pytest.approx(0.0073, abs=0.0005)
    assert _two_sided_p(population) == pytest.approx(0.0164, abs=0.0005)


def test_both_weightings_reduce_to_eq4_for_one_segment():
    """Oracle: Section 8.1 (p. 112), LR = 0.681. At m = 1 both weightings are eq (4)."""
    for weights in c1.C1_WEIGHTS:
        assert c1.statistic(
            [SEGMENT], gamma_estimator=GAMMA_COMPLETE_SAMPLE, weights=weights
        ) == pytest.approx(PUBLISHED["lr_gamma_hat"], abs=0.0005)


# ---------------------------------------------- m > 1: the small bowel record


def test_small_bowel_transcription_reproduces_the_papers_summaries():
    """Oracle: arXiv v1 Section 6.2 (p. 16) and Table 4 (p. 17), from Table I data alone.

    The paper pools one gamma over the 80 complete periods (sample divisor) and uses the
    preprint eq (16) with that common gamma, so LR^m = Laplace^m / gamma_hat. Computed
    here with numpy because the shipped C1 estimates gamma per segment; this pins the
    fixture and the paper's arithmetic, not the shipped path.
    """
    periods = np.concatenate([np.asarray(x, float) for x, _ in SMALL_BOWEL_TABLE_I])
    assert len(SMALL_BOWEL_TABLE_I) == SMALL_BOWEL_PUBLISHED["n_subjects"]
    assert len(periods) == SMALL_BOWEL_PUBLISHED["n_complete_periods"]
    mu, sd = float(periods.mean()), float(periods.std(ddof=1))
    assert mu == pytest.approx(SMALL_BOWEL_PUBLISHED["mu_hat"], abs=0.005)
    assert sd == pytest.approx(SMALL_BOWEL_PUBLISHED["sigma_hat"], abs=0.005)
    assert sd / mu == pytest.approx(SMALL_BOWEL_PUBLISHED["gamma_hat"], abs=0.0005)

    u, taus, ns = [], [], []
    for x, censored in SMALL_BOWEL_TABLE_I:
        t = np.cumsum(np.asarray(x, float))
        tau = float(t[-1] + censored)
        u.append(t.sum() - len(t) * tau / 2.0)
        taus.append(tau)
        ns.append(len(t))
    u, taus, ns = np.array(u), np.array(taus), np.array(ns)
    laplace = float(np.sqrt(12.0) * u.sum() / np.sqrt(np.sum(taus**2 * ns)))
    assert laplace == pytest.approx(SMALL_BOWEL_PUBLISHED["laplace"], abs=0.005)
    assert laplace / (sd / mu) == pytest.approx(
        SMALL_BOWEL_PUBLISHED["lr_multiprocess"], abs=0.005
    )
    gl = float(u.sum() / np.sqrt(np.sum(u**2)))
    assert _two_sided_p(gl) == pytest.approx(SMALL_BOWEL_PUBLISHED["gl_p"], abs=0.0005)


def test_small_bowel_motility_shipped_path():
    """arXiv v1 Section 6.2, m = 19, through shipped C1. SKIPPED - no pooled-gamma mode.

    The data are in the fixture and the paper's arithmetic is reproduced by
    `test_small_bowel_transcription_reproduces_the_papers_summaries`. The shipped C1 cannot
    reproduce LR^m = 3.67 because the paper tests the stronger null of ONE common gap
    distribution and pools gamma across subjects (p. 15), while C1 estimates gamma per
    segment and refuses a segment with fewer than two complete gaps (subject 5 has one).
    Per-segment gamma on the 18 usable subjects gives LR^m = 4.49 under the shipped
    defaults (journal weights, sample divisor), as `scripts/verify_gold_standard.py`
    computes. Activating this needs a pooled-gamma estimator, which is new functionality
    and a design decision, not a fix.
    """
    pytest.skip(
        "C1 has no pooled-gamma estimator; the arXiv v1 Section 6.2 LR^m = 3.67 pools gamma "
        "over all 19 subjects. The m > 1 path is checked externally by "
        "test_c1_multiprocess_reproduces_technometrics_table_4."
    )
