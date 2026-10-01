"""Published recurrent-event records analysed by Kvaloy and Lindqvist, as test fixtures.

Source paper: Kvaloy and Lindqvist, "A class of tests for trend in time censored recurrent
event data", Technometrics 62(1):101-115 (2020), doi:10.1080/00401706.2019.1605936.
Preprint: arXiv:1802.08339v1 (2018), which differs from Technometrics. The two versions
number their sections differently and use different second examples; every locator below
names its version.

Three records:

- `FAILURE_TIMES_H` - one load-haul-dump machine, single process. Original source Kumar,
  Klefsjo and Granholm (1989), Reliability Engineering and System Safety 26:341-361
  (KL2020's reference list, p. 115, and the authors' R code print volume 24).
  Technometrics Section 8.1, Table 1, p. 112 (arXiv v1 Section 6.1, Table 1, p. 14).
- `HYDRAULIC_GAPS_H` - hydraulic systems of six load-haul-dump machines, m = 6. Original
  source Kumar and Klefsjo (1992), Reliability Engineering and System Safety 35:217-224.
  Technometrics Section 8.2, Table 4, pp. 112-113. Not in arXiv v1. The gaps are not
  printed in the Technometrics paper. Kumar and Klefsjo print them in their Appendix,
  p. 224, and they equal the authors' `example_LHDdata_multi.R`
  (github.com/jtkgithub/trendtests, identical to the journal's supplementary material).
- `SMALL_BOWEL_TABLE_I` - migrating motor complex periods of 19 subjects, m = 19. Original
  source Aalen and Husebye (1991), Statistics in Medicine 10:1227-1240, Table I, p. 1229.
  arXiv v1 Section 6.2, pp. 15-17. Not in the Technometrics version.

The single-process record: TIMES are cumulative failure times in hours from the start of
observation, NOT gaps. The record is TIME censored at tau = 2000 hours (Table 1 note), so
the trailing 30 hours after the last failure at 1970 is an incomplete gap - which is the
whole point of the time-censored formulation and the reason this fixture exercises
`Segment` rather than a plain array.
"""

from __future__ import annotations

import numpy as np

# 36 cumulative failure times, hours.
FAILURE_TIMES_H = np.array(
    [
        16,
        39,
        71,
        95,
        98,
        110,
        114,
        226,
        294,
        344,
        555,
        599,
        757,
        822,
        963,
        1077,
        1167,
        1202,
        1257,
        1317,
        1345,
        1372,
        1402,
        1536,
        1625,
        1643,
        1675,
        1726,
        1736,
        1772,
        1796,
        1799,
        1814,
        1868,
        1894,
        1970,
    ],
    dtype=float,
)

# Time censoring point, hours. tau > T_N by 30 hours.
TAU_H = 2000.0


def gaps_h() -> np.ndarray:
    """Inter-failure gaps, hours. `T_1` is itself a gap from the origin at t = 0."""
    return np.diff(np.concatenate([[0.0], FAILURE_TIMES_H]))


# Kvaloy and Lindqvist's printed values for the single-process record. Locators are for
# Technometrics (arXiv v1 in brackets); the numbers are identical in both versions.
#
#   Table 2, p. 112 [p. 14]:     mu_hat, sigma_hat, gamma_hat   (complete gaps)
#                                sigma_tilde, gamma_tilde       (eq (10), with censored time)
#   Section 8.1 text, p. 112 [Section 6.1, p. 14]:
#                                laplace 0.605, lr_gamma_hat 0.681 (= 0.605/0.888),
#                                sigma_star 42.77, gamma_star 0.782, lr_sigma_star 0.774
#
# Every complete-gap number uses the sample 1/(N-1) divisor. `tests/test_checks_published_
# values.py` asserts both the agreement and the size of the population-divisor divergence.
PUBLISHED = {
    "mu_hat": 54.72,  # mean of the complete gaps, Table 2 row 1
    "sigma_tilde": 47.23,  # eq (10), uses the censored time, Table 2 row 2
    "gamma_tilde": 0.850,  # eq (10), Table 2 row 2
    "laplace": 0.605,  # eq (4) with the 1/gamma_hat factor dropped (Section 3.1), text
    "lr_sigma_star": 0.774,  # eq (4) with gamma* = sigma*/mu_hat, eq (11); text, NOT Table 2
    "sigma_hat": 48.61,  # complete gaps, Table 2 row 1 (1/(N-1) divisor)
    "gamma_hat": 0.888,  # complete gaps, Table 2 row 1
    "lr_gamma_hat": 0.681,  # follows from gamma_hat; text
}


# ---------------------------------------------------------------------------
# Hydraulic systems of six load-haul-dump machines, m = 6 (Technometrics Section 8.2)
# ---------------------------------------------------------------------------
#
# Times BETWEEN events for machines 1, 3, 9, 11, 17 and 20, as in the authors' R example.
# "for the purpose of this example we will follow their analysis and assume that for each
# machine the data are time censored at the last recorded event time" (p. 112): the last
# gap of each list ends at tau_j and is not itself an event. `hydraulic_segments_data()`
# applies that convention; it is the authors' convention, not one chosen here.
HYDRAULIC_GAPS_H = {
    1: [327, 125, 7, 6, 107, 277, 54, 332, 510, 110, 10, 9, 85,
        27, 59, 16, 8, 34, 21, 152, 158, 44, 18],
    3: [637, 40, 397, 36, 54, 53, 97, 63, 216, 118, 125, 25, 4,
        101, 184, 167, 81, 46, 18, 32, 219, 405, 20, 248, 140],
    9: [278, 261, 990, 191, 107, 32, 51, 10, 132, 176, 247, 165, 454, 142,
        38, 249, 212, 204, 182, 116, 30, 24, 32, 38, 10, 311, 61],
    11: [353, 96, 49, 211, 82, 175, 79, 117, 26, 4, 5, 60, 39, 35, 258,
         97, 59, 3, 37, 8, 245, 79, 49, 31, 259, 283, 150, 24],
    17: [401, 36, 18, 159, 341, 171, 24, 350, 72, 303, 34, 45, 324, 2, 70, 57,
         103, 11, 5, 3, 144, 80, 53, 84, 218, 122],
    20: [231, 20, 361, 260, 176, 16, 101, 293, 5, 119, 9, 80, 112, 10,
         162, 90, 176, 370, 90, 15, 315, 32, 266],
}  # fmt: skip


def hydraulic_segments_data() -> list[tuple[np.ndarray, float]]:
    """`(complete gaps, tau)` per machine: tau is the last recorded event time."""
    out = []
    for gaps in HYDRAULIC_GAPS_H.values():
        g = np.asarray(gaps, dtype=float)
        out.append((g[:-1], float(g.sum())))
    return out


# Technometrics Table 4, p. 113: p-values, with no tail stated. LR and GL reproduce as
# two-sided normal p-values. LR and ELR use the eq (15) weights; CvM and AD are upper tails
# of sums weighted by tau_j, with a null simulated from the limiting law (so their third
# decimal carries Monte Carlo error); GL is Lawless, Cigsar and Cook's (2012) generalised
# Laplace test, KL2020 Section 5.2, eq (17). The ELR value matches the authors' R code
# (`avec = 1 - 0.5*tau_j/max tau`), not the text's `a_j tau_j = 4743/2`, for which
# `scripts/verify_gold_standard.py` gives 0.0079. Tests assert only `lr_p` and `gl_p`.
HYDRAULIC_PUBLISHED = {
    "lr_p": 0.019,
    "cvm_p": 0.064,
    "ad_p": 0.076,
    "elr_p": 0.003,
    "gl_p": 0.062,
}


# ---------------------------------------------------------------------------
# Small bowel motility, m = 19 (arXiv v1 Section 6.2 only)
# ---------------------------------------------------------------------------
#
# Aalen and Husebye (1991) Table I, p. 1229: per subject, the complete MMC periods in
# minutes, then the censored period that ends the record. Each record starts at a phase III
# (an event), and is time censored at the end of measurement.
SMALL_BOWEL_TABLE_I = [
    ([112, 145, 39, 52, 21, 34, 33, 51], 54),
    ([206, 147], 30),
    ([284, 59, 186], 4),
    ([94, 98, 84], 87),
    ([67], 131),
    ([124, 34, 87, 75, 43, 38, 58, 142, 75], 23),
    ([116, 71, 83, 68, 125], 111),
    ([111, 59, 47, 95], 110),
    ([98, 161, 154, 55], 44),
    ([166, 56], 122),
    ([63, 90, 63, 103, 51], 85),
    ([47, 86, 68, 144], 72),
    ([120, 106, 176], 6),
    ([112, 25, 57, 166], 85),
    ([132, 267, 89], 86),
    ([120, 47, 165, 64, 113], 12),
    ([162, 141, 107, 69], 39),
    ([106, 56, 158, 41, 41, 168], 13),
    ([147, 134, 78, 66, 100], 4),
]

# arXiv v1 Section 6.2, p. 16, and Table 4, p. 17. The paper tests the STRONGER null that
# all 19 processes share one gap distribution and so estimates ONE gamma from the 80
# complete periods pooled (p. 15), sample divisor. Its Table 4 prints LR p < 0.0001, which
# contradicts the p = 0.00024 of its own text; the text value is the one recorded here. No
# test asserts it: it is 2 Phi(-3.67) of the printed statistic, which the tests pin.
SMALL_BOWEL_PUBLISHED = {
    "n_subjects": 19,
    "n_complete_periods": 80,
    "mu_hat": 98.76,
    "sigma_hat": 52.62,
    "gamma_hat": 0.533,
    "laplace": 1.95,
    "lr_multiprocess": 3.67,
    "p_value": 0.00024,
    "gl_p": 0.007,
}
