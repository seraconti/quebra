#!/usr/bin/env python3
"""Recompute the "recomputed" columns of GOLD_STANDARD.md sections 3, 4 and 5.

Contract:

- numpy and scipy only. Nothing is imported from quebra, and no formula is copied from
  `analyzers/checks/` except `adinf_repo`, a copy of the shipped `adinf` approximation,
  whose error against the exact limit is what it is here to measure. Every statistic is
  computed twice where that is possible: once by integrating the tied-down counting
  process piece by piece (derived here), once from the equation as printed in the paper
  (eq (6) with the `/tau^2` its derivation gives; the printed `/tau` is evaluated too).
  The two must agree before a value is reported.
- The data are transcribed by eye from the printed pages:
  - load-haul-dump, one machine: Kvaloy and Lindqvist, Technometrics 62(1) 2020 (KL2020),
    Table 1, p. 112;
  - hydraulic systems, six machines: Kumar and Klefsjo, RESS 35 (1992), Appendix table,
    p. 224. The same gaps are parsed from the authors' `example_LHDdata_multi.R` inside
    `8040713.zip` when that archive is given, and the two must be equal;
  - small bowel motility: Aalen and Husebye, Stat Med 10 (1991), Table I, p. 1229.
  `--fixture` optionally compares these transcriptions with the repository fixture by
  parsing it as text; the fixture is never the data source.
- Limiting laws are evaluated by numerical inversion (Imhof) of the closed-form
  characteristic functions: A^2 from Anderson and Darling (1954) eq (7), p. 768, and W^2
  from the product sin(sqrt(s))/sqrt(s). They are checked against two series: Anderson and
  Darling (1954) eq (8) for A^2, and the Bessel-K series for W^2 (GOLD_STANDARD attributes
  it to Anderson and Darling 1952 eq (4.35); that paper is not in spec/sources).
- The two simulated nulls (AD^m and CvM^m with weights proportional to tau_j) are drawn from
  the limit with a fixed seed and reported with their Monte Carlo standard errors, next to
  the exact value from the same Imhof inversion.

Every row prints this recomputation, the value GOLD_STANDARD prints for it, and the
published value. Those GOLD strings are held in this file, so each one is also looked up in
`docs/GOLD_STANDARD.md` (`--gold`), and a string the doc no longer contains fails the run:
an edit to the doc cannot pass unseen. The REPO rows are held to the docstring of
`c2_anderson_darling.ad_limiting_cdf` the same way, read as text and never imported.

Exit status 1 if a GOLD or REPO row differs from its recomputation by more than half a
unit in its last printed digit, if a GOLD or REPO string is missing from its file, or if
the `--fixture` comparison fails. A `--klcode-zip` mismatch or a failed internal
cross-check raises. The printed column is informational and never sets the exit status.
"""

from __future__ import annotations

import argparse
import ast
import math
import re
import sys
import zipfile
from pathlib import Path

import numpy as np
from scipy import integrate, optimize, special, stats

SQRT12 = math.sqrt(12.0)

# --------------------------------------------------------------------------- data

# KL2020 Table 1, p. 112: failure times, hours. "The data are time censored at 2000 hr."
LHD_TIMES_H = [
    16, 39, 71, 95, 98, 110, 114, 226, 294, 344, 555, 599,
    757, 822, 963, 1077, 1167, 1202, 1257, 1317, 1345, 1372, 1402, 1536,
    1625, 1643, 1675, 1726, 1736, 1772, 1796, 1799, 1814, 1868, 1894, 1970,
]  # fmt: skip
LHD_TAU_H = 2000.0

# Kumar and Klefsjo (1992), Appendix, p. 224: time between successive failures of the
# hydraulic systems, hours. KL2020 Section 8.2, p. 112, follows their analysis: each
# machine is time censored at its last recorded event, so the last gap is the censored one.
HYDRAULIC_KK1992 = {
    1: [327, 125, 7, 6, 107, 277, 54, 332, 510, 110, 10, 9, 85, 27, 59, 16, 8, 34, 21,
        152, 158, 44, 18],
    3: [637, 40, 397, 36, 54, 53, 97, 63, 216, 118, 125, 25, 4, 101, 184, 167, 81, 46,
        18, 32, 219, 405, 20, 248, 140],
    9: [278, 261, 990, 191, 107, 32, 51, 10, 132, 176, 247, 165, 454, 142, 38, 249, 212,
        204, 182, 116, 30, 24, 32, 38, 10, 311, 61],
    11: [353, 96, 49, 211, 82, 175, 79, 117, 26, 4, 5, 60, 39, 35, 258, 97, 59, 3, 37, 8,
         245, 79, 49, 31, 259, 283, 150, 24],
    17: [401, 36, 18, 159, 341, 171, 24, 350, 72, 303, 34, 45, 324, 2, 70, 57, 103, 11, 5,
         3, 144, 80, 53, 84, 218, 122],
    20: [231, 20, 361, 260, 176, 16, 101, 293, 5, 119, 9, 80, 112, 10, 162, 90, 176, 370,
         90, 15, 315, 32, 266],
}  # fmt: skip

# Aalen and Husebye (1991) Table I, p. 1229: complete MMC periods, then the censored one,
# minutes. Each record starts at a phase III.
BOWEL_AH1991 = [
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

# --------------------------------------------------------------------------- reporting

ROWS: list[tuple[str, str, float, str | None, str | None, str, float | None]] = []


def _decimals(text: str) -> int:
    text = text.strip().lstrip("+-")
    if "e" in text.lower():
        mant, exp = text.lower().split("e")
        return _decimals(mant) - int(exp)
    return len(text.split(".")[1]) if "." in text else 0


def row(
    section: str,
    label: str,
    value: float,
    gold: str | None,
    printed: str | None = None,
    note: str = "",
    mc_tol: float | None = None,
) -> None:
    """Record one comparison. `gold` and `printed` are the strings as printed.

    `mc_tol` is set only on a Monte Carlo row: GOLD's own value there is a simulation with a
    stated standard error, so the row agrees when the difference is within half a unit in
    the last digit OR within `mc_tol` (three combined Monte Carlo standard errors). The
    half-unit verdict is printed as well, so neither reading is hidden.
    """
    ROWS.append((section, label, float(value), gold, printed, note, mc_tol))


def verdict(value: float, ref: str | None, mc_tol: float | None = None) -> str:
    if ref is None:
        return "-"
    half = 0.5 * 10.0 ** (-_decimals(ref))
    diff = value - float(ref)
    if abs(diff) <= half * (1 + 1e-9):
        return "agree"
    if mc_tol is not None and abs(diff) <= mc_tol:
        return f"agree within 3 MC s.e. ({diff:+.2g})"
    return f"DISAGREE ({diff:+.3g})"


def fmt(value: float) -> str:
    if value == 0 or 1e-3 <= abs(value) < 1e5:
        return f"{value:.10g}"
    return f"{value:.6e}"


def print_rows() -> tuple[int, int]:
    """Print every row; return the disagreement counts for GOLD rows and for REPO rows."""
    n_bad = n_repo = 0
    current = None
    for section, label, value, gold, printed, note, mc_tol in ROWS:
        if section != current:
            print(f"\n== {section}")
            print(
                f"   {'quantity':<46} {'recomputed':>16} | {'GOLD':>10} {'vs GOLD':<22} |"
                f" {'printed':>9} {'vs printed':<20}"
            )
            current = section
        vg = verdict(value, gold, mc_tol)
        vp = verdict(value, printed)
        if section.startswith("GOLD"):
            n_bad += vg.startswith("DISAGREE")
        else:
            n_repo += vg.startswith("DISAGREE")
        print(
            f"   {label:<46} {fmt(value):>16} | {gold or '':>10} {vg:<22} |"
            f" {printed or '':>9} {vp:<20}{('  ' + note) if note else ''}"
        )
    return n_bad, n_repo


# --------------------------------------------------------------------------- estimators


def gaps(times: np.ndarray) -> np.ndarray:
    """Complete gaps: T_1 is a gap from the origin."""
    return np.diff(np.concatenate([[0.0], times]))


def moments(x: np.ndarray, ddof: int) -> tuple[float, float, float]:
    mu = float(np.mean(x))
    sd = float(np.sqrt(np.sum((x - mu) ** 2) / (len(x) - ddof)))
    return mu, sd, sd / mu


def eq10(x: np.ndarray, tau: float) -> tuple[float, float, float]:
    """KL2020 eq (10), p. 104: estimators that use the censored time tau - T_N."""
    n = len(x)
    mu = tau / n
    var = (np.sum(x**2) + (tau - np.sum(x)) ** 2) / n - mu**2
    return mu, math.sqrt(var), math.sqrt(var) / mu


def eq11(x: np.ndarray) -> float:
    """KL2020 eq (11), p. 104: sigma* from successive differences."""
    n = len(x)
    return math.sqrt(np.sum(np.diff(x) ** 2) / (2 * (n - 1)))


def weibull_mle(complete: np.ndarray, censored: np.ndarray) -> dict[str, float]:
    """Weibull MLE with right-censored observations, by the profile score for the shape.

    Checked against a direct two-parameter maximisation of the log-likelihood.
    """
    t = np.concatenate([complete, censored]).astype(float)
    d = len(complete)
    scale0 = t.max()
    ts = t / scale0
    mean_log_x = float(np.mean(np.log(complete / scale0)))

    def score(k: float) -> float:
        tk = ts**k
        return float(np.sum(tk * np.log(ts)) / np.sum(tk) - 1.0 / k - mean_log_x)

    k = optimize.brentq(score, 0.05, 50.0, xtol=1e-15, rtol=1e-15)
    lam = scale0 * (np.sum(ts**k) / d) ** (1.0 / k)

    def nll(p: np.ndarray) -> float:
        kk, ll = math.exp(p[0]), math.exp(p[1])
        logf = (
            np.log(kk / ll) + (kk - 1) * np.log(complete / ll) - (complete / ll) ** kk
        )
        return float(-(np.sum(logf) - np.sum((censored / ll) ** kk)))

    res = optimize.minimize(
        nll,
        x0=[0.0, math.log(np.mean(t))],
        method="Nelder-Mead",
        options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 20000},
    )
    k2, lam2 = math.exp(res.x[0]), math.exp(res.x[1])
    g1 = special.gamma(1 + 1 / k)
    g2 = special.gamma(1 + 2 / k)
    mu = lam * g1
    sd = lam * math.sqrt(g2 - g1**2)
    return {
        "shape": k,
        "scale": lam,
        "mu": mu,
        "sigma": sd,
        "gamma": sd / mu,
        "shape_direct": k2,
        "scale_direct": lam2,
    }


# --------------------------------------------------------------------------- statistics
# The tied-down process: V(s) = (N(s tau) - s N) / (gamma sqrt(N)), s in [0, 1]. N(s tau)
# equals i on [T_i/tau, T_{i+1}/tau) with T_0 = 0 and T_{N+1} = tau, so V is linear on
# each piece and every functional below is an exact sum over pieces.


def _pieces(
    times: np.ndarray, tau: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    s = np.concatenate([[0.0], times / tau, [1.0]])
    return s[:-1], s[1:], np.arange(len(s) - 1, dtype=float), len(times)


def _int_nu(i: np.ndarray, lo: np.ndarray, hi: np.ndarray, n: int) -> np.ndarray:
    """Integral of (i - n s) over [lo, hi]."""
    return i * (hi - lo) - n * (hi**2 - lo**2) / 2.0


def lr_direct(times, tau, g):
    a, b, i, n = _pieces(times, tau)
    return -SQRT12 * np.sum(_int_nu(i, a, b, n)) / (g * math.sqrt(n))


def ks_direct(times, tau, g):
    a, b, i, n = _pieces(times, tau)
    return max(np.max(np.abs(i - n * a)), np.max(np.abs(i - n * b))) / (
        g * math.sqrt(n)
    )


def cvm_direct(times, tau, g):
    a, b, i, n = _pieces(times, tau)
    return np.sum(((n * b - i) ** 3 - (n * a - i) ** 3) / (3.0 * n)) / (g**2 * n)


def ad_direct(times, tau, g):
    """(i - n s)^2 / (s(1-s)) = -n^2 + i^2/s + (n-i)^2/(1-s), integrated piece by piece."""
    a, b, i, n = _pieces(times, tau)
    tot = -(n**2) * np.sum(b - a)
    m = i > 0
    tot += np.sum(i[m] ** 2 * np.log(b[m] / a[m]))
    m = i < n
    tot += np.sum((n - i[m]) ** 2 * np.log((1 - a[m]) / (1 - b[m])))
    return tot / (g**2 * n)


def elr_direct(times, tau, g, aa):
    """(int_0^a V - int_a^1 V) / sqrt(1/12 - a^2 (1-a)^2): KL2020 eq (8), scaled."""
    a, b, i, n = _pieces(times, tau)
    lo_hi = np.minimum(b, aa)
    below = np.where(a < aa, _int_nu(i, a, lo_hi, n), 0.0)
    hi_lo = np.maximum(a, aa)
    above = np.where(b > aa, _int_nu(i, hi_lo, b, n), 0.0)
    return (np.sum(below) - np.sum(above)) / (
        g * math.sqrt(n) * math.sqrt(1 / 12 - aa**2 * (1 - aa) ** 2)
    )


def quad_check(times, tau, g) -> tuple[float, float]:
    """CvM and AD by adaptive quadrature of V^2 and V^2/(s(1-s)) on each piece."""
    a, b, i, n = _pieces(times, tau)
    cvm = ad = 0.0
    for lo, hi, ii in zip(a, b, i):
        cvm += integrate.quad(
            lambda s: (ii - n * s) ** 2, lo, hi, epsabs=0, epsrel=1e-13
        )[0]
        ad += integrate.quad(
            lambda s: (ii - n * s) ** 2 / (s * (1 - s)),
            lo,
            hi,
            epsabs=0,
            epsrel=1e-13,
            limit=200,
        )[0]
    return cvm / (g**2 * n), ad / (g**2 * n)


# The same statistics from the equations as printed in KL2020 (pp. 103-104).


def lr_eq4(times, tau, g):
    n = len(times)
    return SQRT12 * (np.sum(times) - n * tau / 2) / (g * tau * math.sqrt(n))


def ks_eq5(times, tau, g):
    n = len(times)
    i = np.arange(1, n + 1)
    u = n * times / tau
    return np.max(np.maximum(np.abs(i - u), np.abs(i - 1 - u))) / (g * math.sqrt(n))


def cvm_eq6(times, tau, g, middle_power: int = 2):
    """eq (6). As printed the middle term divides by tau; middle_power=2 is the integral."""
    n = len(times)
    t = np.concatenate([[0.0], times])
    i = np.arange(0, n)
    x_next = t[1:] - t[:-1]
    s = np.sum(
        i**2 * x_next / tau - i * n * (t[1:] ** 2 - t[:-1] ** 2) / tau**middle_power
    )
    s += n**2 * (times[-1] ** 2 / tau**2 - times[-1] / tau + 1 / 3)
    return s / (g**2 * n)


def ad_eq7(times, tau, g):
    n = len(times)
    i = np.arange(1, n)
    ti, tip = times[:-1], times[1:]
    s = np.sum(
        (n - i) ** 2 * np.log((tau - ti) / (tau - tip)) + i**2 * np.log(tip / ti)
    )
    s += n**2 * (math.log(tau / (tau - times[0])) + math.log(tau / times[-1]) - 1)
    return s / (g**2 * n)


def elr_eq9(times, tau, g, aa):
    n = len(times)
    brace = np.sum(np.abs(times - aa * tau)) - (0.5 - aa * (1 - aa)) * tau * n
    return brace / (g * tau * math.sqrt(n) * math.sqrt(1 / 12 - aa**2 * (1 - aa) ** 2))


def agreed(*values: float, tol: float = 1e-10) -> float:
    ref = values[0]
    for v in values[1:]:
        if not abs(v - ref) <= tol * max(1.0, abs(ref)):
            raise AssertionError(f"implementations disagree: {values}")
    return ref


def reversal_selftest(times: np.ndarray, tau: float, a_values) -> None:
    """Time reversal t -> tau - t maps V(s) to -V(1 - s).

    So LR changes sign, KS, CvM and AD are unchanged, and ELR(a) becomes ELR(1 - a). Both
    implementations of each statistic must obey this. It is the only check here that sees a
    one-sided KS: on the LHD record the supremum is a negative excursion, attained at the
    right end of a piece, so a KS that skipped left ends would still agree with eq (5); on
    the reversed record the supremum moves to a left end.
    """
    rev = np.sort(tau - times)
    for f in (lr_eq4, lr_direct):
        agreed(f(rev, tau, 1.0), -f(times, tau, 1.0))
    for f in (ks_eq5, ks_direct, cvm_eq6, cvm_direct, ad_eq7, ad_direct):
        agreed(f(rev, tau, 1.0), f(times, tau, 1.0))
    agreed(ks_eq5(rev, tau, 1.0), ks_direct(rev, tau, 1.0))
    agreed(cvm_eq6(rev, tau, 1.0), cvm_direct(rev, tau, 1.0))
    agreed(ad_eq7(rev, tau, 1.0), ad_direct(rev, tau, 1.0))
    for aa in a_values:
        for f in (elr_eq9, elr_direct):
            agreed(f(rev, tau, 1.0, 1.0 - aa), f(times, tau, 1.0, aa))


# --------------------------------------------------------------------------- limit laws
# A quadratic form Q = sum_k c_k Z_k^2 with c_k = w * lambda_k. Its survival function by
# Imhof (1961): P(Q > x) = 1/2 + (1/pi) int_0^inf sin(A(u) - x u/2) / (u rho(u)) du, with
# A(u) = (1/2) sum arctan(c_k u) = (1/2) arg P(w u), log rho(u) = (1/2) log|P(w u)|, and
# P(u) = prod_k (1 + i u lambda_k) in closed form:
#   A^2:  lambda_k = 1/(k(k+1)),  P = -cos(pi sqrt(1/4 + s)) / (pi s),  s = -i u
#   W^2:  lambda_k = 1/(k pi)^2,  P = sin(sqrt(s)) / sqrt(s),           s = -i u
# The first is AD1954 eq (7) with s = 2it; both follow from the Gamma reflection formula.

LAMBDA_SUM = {"ad": 1.0, "cvm": 1.0 / 6.0}
LAMBDA_SUMSQ = {"ad": (math.pi**2 - 9.0) / 3.0, "cvm": 1.0 / 90.0}


def log_prod(family: str, u: np.ndarray) -> np.ndarray:
    """log prod_k (1 + i u lambda_k), imaginary part modulo 2 pi."""
    u = np.asarray(u, dtype=float)
    s = -1j * u
    out = np.empty(u.shape, dtype=complex)
    small = u < 1.0
    if family == "ad":
        c = np.sqrt(0.25 + s[small])
        out[small] = np.log(-np.cos(np.pi * c) / (np.pi * s[small]))
        z = np.pi * np.sqrt(0.25 + s[~small])  # Im z < 0
        logcos = 1j * z - math.log(2.0) + np.log1p(np.exp(-2j * z))
        out[~small] = 1j * math.pi + logcos - math.log(math.pi) - np.log(s[~small])
    elif family == "cvm":
        z = np.sqrt(s[small])
        out[small] = np.log(np.sin(z) / z)
        z = np.sqrt(s[~small])  # Im z < 0
        logsin = 1j * z - np.log(2j) + np.log1p(-np.exp(-2j * z))
        out[~small] = logsin - np.log(z)
    else:
        raise ValueError(family)
    return out


def log_prod_truncated(family: str, u: float, k_max: int = 200000) -> complex:
    """The same product by direct summation plus a tail series: an independent check."""
    k = np.arange(1, k_max + 1, dtype=float)
    lam = 1 / (k * (k + 1)) if family == "ad" else 1 / (k * math.pi) ** 2
    s1 = LAMBDA_SUM[family] - lam.sum()
    s2 = LAMBDA_SUMSQ[family] - (lam**2).sum()
    arg = np.sum(np.arctan(u * lam)) + u * s1
    mod = 0.5 * np.sum(np.log1p((u * lam) ** 2)) + 0.5 * u**2 * s2
    return complex(mod, arg)


class LimitSum:
    """P(sum_j w_j X_j > x) for independent X_j from one limiting law."""

    def __init__(
        self,
        family: str,
        weights,
        panel: float = 0.5,
        nodes: int = 24,
        log_rho_cut: float = 60.0,
    ):
        self.family = family
        self.w = np.atleast_1d(np.asarray(weights, dtype=float))
        # Find U where rho(U) > exp(log_rho_cut): the integrand beyond is negligible.
        upper = 1.0
        while self._log_rho(np.array([upper]))[0] < log_rho_cut:
            upper *= 1.5
        n_panels = int(math.ceil(upper / panel))
        x, wt = np.polynomial.legendre.leggauss(nodes)
        left = np.arange(n_panels) * panel
        self.u = (left[:, None] + (x[None, :] + 1) * panel / 2).ravel()
        self.wt = np.tile(wt * panel / 2, n_panels)
        a = np.zeros_like(self.u)
        lr = np.zeros_like(self.u)
        for wj in self.w:
            lp = log_prod(self.family, wj * self.u)
            a += 0.5 * np.unwrap(lp.imag)
            lr += 0.5 * lp.real
        self.a = a
        self.kernel = self.wt * np.exp(-lr) / self.u
        self.upper = upper

    def _log_rho(self, u):
        return sum(0.5 * log_prod(self.family, wj * u).real for wj in self.w)

    def sf(self, x: float) -> float:
        return (
            0.5 + float(np.sum(np.sin(self.a - x * self.u / 2) * self.kernel)) / math.pi
        )

    def cdf(self, x: float) -> float:
        return 1.0 - self.sf(x)

    def isf(self, alpha: float) -> float:
        lo, hi = 1e-3, 1.0
        while self.sf(hi) > alpha:
            hi *= 2
        return optimize.brentq(
            lambda x: self.sf(x) - alpha, lo, hi, xtol=1e-14, rtol=1e-15
        )


def ad_cdf_eq8(z: float) -> float:
    """AD1954 eq (8), p. 768, term by term with each integral by adaptive quadrature."""
    total = 0.0
    for j in range(200):
        k = 4 * j + 1
        coef = (
            (-1) ** j * math.exp(special.gammaln(j + 0.5) - special.gammaln(j + 1)) * k
        )

        def f(w, k=k):
            return math.exp(
                z / (8 * (w * w + 1)) - k * k * math.pi**2 * (1 + w * w) / (8 * z)
            )

        val = integrate.quad(f, 0, np.inf, epsabs=1e-300, epsrel=1e-13, limit=400)[0]
        term = math.sqrt(2.0) / z * coef * val
        total += term
        if j > 2 and abs(term) < 1e-18:
            break
    return total


def cvm_cdf_bessel(x: float) -> float:
    """W^2 limit, sum_k Gamma(k+1/2)/(Gamma(1/2) k!) sqrt(4k+1) e^{-q} K_{1/4}(q)/(pi sqrt x)."""
    total = 0.0
    for k in range(200):
        y = 4 * k + 1
        q = y * y / (16 * x)
        c = math.exp(special.gammaln(k + 0.5) - special.gammaln(k + 1)) / math.pi**1.5
        term = c * math.sqrt(y) * math.exp(-2 * q) * special.kve(0.25, q) / math.sqrt(x)
        total += term
        if k > 2 and term < 1e-18:
            break
    return total


def kolmogorov_sf(x: float) -> float:
    """2 sum (-1)^(k-1) exp(-2 k^2 x^2), checked against scipy.special.kolmogorov."""
    k = np.arange(1, 200)
    own = float(2 * np.sum((-1.0) ** (k - 1) * np.exp(-2 * k**2 * x**2)))
    return agreed(own, float(special.kolmogorov(x)), tol=1e-12)


def sample_limit(
    family: str,
    rng: np.random.Generator,
    n: int,
    k_terms: int = 100,
    chunk: int = 20000,
) -> np.ndarray:
    """Draws of sum_k lambda_k Z_k^2: k_terms exact, the tail as a normal with its moments."""
    k = np.arange(1, k_terms + 1, dtype=float)
    lam = 1 / (k * (k + 1)) if family == "ad" else 1 / (k * math.pi) ** 2
    s1 = LAMBDA_SUM[family] - lam.sum()
    s2 = LAMBDA_SUMSQ[family] - (lam**2).sum()
    out = np.empty(n)
    for start in range(0, n, chunk):
        m = min(chunk, n - start)
        z = rng.standard_normal((m, k_terms))
        out[start : start + m] = (
            (z * z) @ lam + s1 + math.sqrt(2 * s2) * rng.standard_normal(m)
        )
    return out


def adinf_repo(z: float) -> float:
    """The 13-coefficient approximation shipped in c2_anderson_darling.ad_limiting_cdf.

    Coefficients transcribed from the repository source (GOLD_STANDARD attributes them to
    Marsaglia and Marsaglia 2004, not in spec/sources). Used only to measure its error.
    """
    if z < 2:
        return (
            math.exp(-1.2337141 / z)
            / math.sqrt(z)
            * (
                2.00012
                + (
                    0.247105
                    - (0.0649821 - (0.0347962 - (0.011672 - 0.00168691 * z) * z) * z)
                    * z
                )
                * z
            )
        )
    return math.exp(
        -math.exp(
            1.0776
            - (
                2.30695
                - (0.43424 - (0.082433 - (0.008056 - 0.0003146 * z) * z) * z) * z
            )
            * z
        )
    )


# --------------------------------------------------------------------------- multi-process


def segments_hydraulic(data: dict[int, list[int]]):
    out = []
    for machine, g in data.items():
        g = np.asarray(g, dtype=float)
        x = g[:-1]
        out.append({"id": machine, "x": x, "t": np.cumsum(x), "tau": float(g.sum())})
    return out


def segments_bowel(data):
    out = []
    for j, (complete, cens) in enumerate(data, start=1):
        x = np.asarray(complete, dtype=float)
        out.append(
            {
                "id": j,
                "x": x,
                "t": np.cumsum(x),
                "tau": float(x.sum() + cens),
                "cens": float(cens),
            }
        )
    return out


def lr_m(segs, weights: str, ddof: int) -> float:
    """sum_j A_j LR_j with A_j normalised to unit sum of squares (KL2020 eq (13))."""
    lr = np.array([lr_eq4(s["t"], s["tau"], moments(s["x"], ddof)[2]) for s in segs])
    g = np.array([moments(s["x"], ddof)[2] for s in segs])
    n = np.array([len(s["t"]) for s in segs], dtype=float)
    tau = np.array([s["tau"] for s in segs])
    if weights == "KL2020":  # eq (14), p. 105
        a = np.sqrt(n) / g
    elif weights == "KLv1":  # arXiv v1 eq (14), p. 8
        a = g * tau * np.sqrt(n)
    else:
        raise ValueError(weights)
    a = a / np.sqrt(np.sum(a**2))
    value = float(np.sum(a * lr))
    u = np.array([np.sum(s["t"]) - len(s["t"]) * s["tau"] / 2 for s in segs])
    if weights == "KL2020":  # eq (16), p. 105
        closed = SQRT12 / math.sqrt(np.sum(n / g**2)) * np.sum(u / tau / g**2)
    else:  # arXiv v1 eq (16), p. 8
        closed = SQRT12 / math.sqrt(np.sum(g**2 * tau**2 * n)) * np.sum(u)
    return agreed(value, float(closed))


def lr_m_mixed(segs, ddof_lr: int, ddof_w: int) -> float:
    """KL2020 weights with gamma_j from `ddof_w` in A_j and from `ddof_lr` in LR_j.

    KLcode's `findwvec("sqrtNCV")` computes the weight CV from `var(xvec)` whatever `sigma`
    says (unless sigma = "l"), so its weights are always 1/(N-1).
    """
    n = np.array([len(s["t"]) for s in segs], dtype=float)
    g_w = np.array([moments(s["x"], ddof_w)[2] for s in segs])
    a = np.sqrt(n) / g_w
    a = a / np.sqrt(np.sum(a**2))
    lr = np.array([lr_eq4(s["t"], s["tau"], moments(s["x"], ddof_lr)[2]) for s in segs])
    return float(np.sum(a * lr))


def linrank_klcode(times: np.ndarray, use_order: bool) -> float:
    """KLcode `LinRanktestobs`, which assigns scores through R's `order()`, not `rank()`.

    `use_order=False` substitutes the ranks. The printed LinR p-value decides which one the
    authors ran; LCC2012, the method's source, is not in spec/sources.
    """
    x = gaps(times)
    n = len(x)
    if use_order:
        r = (
            np.argsort(x, kind="stable") + 1
        )  # R order(): 1-based indices of the sorted gaps
    else:
        r = stats.rankdata(x, method="ordinal").astype(int)
    e = np.array([np.sum(1.0 / (n + 1 - np.arange(1, rj + 1))) for rj in r])
    j = np.arange(1, n + 1) - (n + 1) / 2
    u = float(np.sum(e * j))
    var = float(np.sum(j**2) * np.sum((e - e.mean()) ** 2 / (n - 1)))
    return u / math.sqrt(var)


def gl(segs) -> float:
    """KL2020 eq (17), p. 106."""
    u = np.array([np.sum(s["t"]) - len(s["t"]) * s["tau"] / 2 for s in segs])
    return float(np.sum(u) / math.sqrt(np.sum(u**2)))


def p2(z: float) -> float:
    return float(2 * stats.norm.sf(abs(z)))


# --------------------------------------------------------------------------- inputs


def parse_klcode(zip_path: Path) -> dict[int, list[int]]:
    with zipfile.ZipFile(zip_path) as zf:
        name = next(n for n in zf.namelist() if n.endswith("example_LHDdata_multi.R"))
        text = zf.read(name).decode("latin-1")
    out = {}
    for m in re.finditer(r"LHD(\d+)_xtimes\s*<-\s*c\(([^)]*)\)", text):
        out[int(m.group(1))] = [int(v) for v in re.findall(r"\d+", m.group(2))]
    return out


def parse_fixture(path: Path) -> dict[str, object]:
    """Read literal assignments from the fixture module as text. Nothing is executed."""
    tree = ast.parse(path.read_text())
    out: dict[str, object] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            name = getattr(node.targets[0], "id", None)
            if name in {"HYDRAULIC_GAPS_H", "SMALL_BOWEL_TABLE_I"}:
                out[name] = ast.literal_eval(node.value)
            if name == "FAILURE_TIMES_H":
                out[name] = ast.literal_eval(node.value.args[0])
    return out


# --------------------------------------------------------------------------- main

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SEED = 20260929
DEFAULT_DRAWS = 400_000
DEFAULT_BIG_DRAWS = 4_000_000


def _docstring_of(path: Path, function: str) -> str:
    """The docstring of one top-level function, read from the file as text."""
    for node in ast.parse(path.read_text()).body:
        if isinstance(node, ast.FunctionDef) and node.name == function:
            return ast.get_docstring(node) or ""
    raise SystemExit(f"{function} not found in {path}")


def strings_missing_from(
    gold_path: Path, repo_path: Path
) -> list[tuple[str, str, str]]:
    """Every GOLD or REPO string that its file no longer contains.

    The GOLD strings must occur in the gold-standard doc, the REPO strings in the docstring
    they quote. A string that is absent means the file changed and this script did not.
    """
    gold_text = gold_path.read_text()
    repo_text = " ".join(_docstring_of(repo_path, "ad_limiting_cdf").split())
    missing = []
    for section, label, _value, gold, _printed, _note, _tol in ROWS:
        if gold is None:
            continue
        if section.startswith("GOLD") and gold not in gold_text:
            missing.append((str(gold_path.name), label, gold))
        elif section.startswith("REPO") and gold not in repo_text:
            missing.append((str(repo_path.name), label, gold))
    return missing


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--klcode-zip", type=Path, default=None)
    ap.add_argument("--fixture", type=Path, default=None)
    ap.add_argument("--draws", type=int, default=DEFAULT_DRAWS)
    ap.add_argument("--big-draws", type=int, default=DEFAULT_BIG_DRAWS)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--gold", type=Path, default=ROOT / "docs" / "GOLD_STANDARD.md")
    ap.add_argument(
        "--repo-docstring",
        type=Path,
        default=ROOT / "src/quebra/analyzers/checks/c2_anderson_darling.py",
    )
    ap.add_argument(
        "--json", type=Path, default=None, help="also write every row as JSON"
    )
    args = ap.parse_args(argv)
    if args.draws < 400_000:
        ap.error("--draws must be at least 4e5")
    # The CvM^m draws follow AD^m's on one generator, so all three settle the values.
    default_mc = (args.seed, args.draws, args.big_draws) == (
        DEFAULT_SEED,
        DEFAULT_DRAWS,
        DEFAULT_BIG_DRAWS,
    )

    print(
        "verify_gold_standard: numpy",
        np.__version__,
        "| scipy",
        __import__("scipy").__version__,
    )

    # ---- transcription checks
    print("\n== transcription")
    if args.klcode_zip is not None:
        kl = parse_klcode(args.klcode_zip)
        same = kl == HYDRAULIC_KK1992
        print(
            f"   KLcode example_LHDdata_multi.R == Kumar and Klefsjo (1992) p. 224: {same}"
        )
        if not same:
            raise SystemExit("hydraulic transcription differs from the authors' R code")
    n_fixture = 0
    if args.fixture is not None:
        fx = parse_fixture(args.fixture)
        expected = {
            "FAILURE_TIMES_H": ("KL2020 Table 1", [float(v) for v in LHD_TIMES_H]),
            "HYDRAULIC_GAPS_H": ("KK1992 p. 224", HYDRAULIC_KK1992),
            "SMALL_BOWEL_TABLE_I": (
                "AH1991 Table I",
                [(list(a), b) for a, b in BOWEL_AH1991],
            ),
        }
        for name, (source, ours) in expected.items():
            if name not in fx:
                print(f"   fixture {name}: absent from {args.fixture}")
                n_fixture += 1
                continue
            theirs = fx[name]
            if name == "FAILURE_TIMES_H":
                theirs = [float(v) for v in theirs]
            elif name == "SMALL_BOWEL_TABLE_I":
                theirs = [(list(a), b) for a, b in theirs]
            print(f"   fixture {name} == {source}:", theirs == ours)
            n_fixture += theirs != ours

    # ---- limit laws first: the p-values below use them
    ad1 = LimitSum("ad", [1.0])
    cvm1 = LimitSum("cvm", [1.0])
    for fam, lim, xs in (
        ("ad", ad1, [0.5, 1.0, 2.0, 3.0, 6.0]),
        ("cvm", cvm1, [0.1, 0.3, 0.6, 1.0]),
    ):
        for x in xs:
            other = ad_cdf_eq8(x) if fam == "ad" else cvm_cdf_bessel(x)
            agreed(lim.cdf(x), other, tol=2e-12)
        for u in (0.3, 7.0, 150.0):
            a = log_prod(fam, np.array([u]))[0]
            b = log_prod_truncated(fam, u)
            # imaginary parts agree modulo 2 pi; the real parts must agree outright
            d_im = (a.imag - b.imag + math.pi) % (2 * math.pi) - math.pi
            assert abs(a.real - b.real) < 1e-9 and abs(d_im) < 1e-9, (fam, u, a, b)
    print("\n== limit-law self-checks")
    print(
        f"   Imhof on eq (7) vs AD1954 eq (8) series: agree to 2e-12 at z in"
        f" {{0.5,1,2,3,6}}; integration upper limit U = {ad1.upper:.0f}"
    )
    print(
        f"   Imhof on sin(sqrt s)/sqrt s vs Bessel-K series: agree to 2e-12 at x in"
        f" {{0.1,0.3,0.6,1}}; U = {cvm1.upper:.0f}"
    )
    print("   closed-form products vs direct sums (2e5 terms + tail): agree to 1e-9")
    for fam, lim, mean_ref, var_ref in (
        ("A^2", ad1, 1.0, 2 * (math.pi**2 - 9) / 3),
        ("W^2", cvm1, 1 / 6, 1 / 45),
    ):
        m1 = integrate.quad(lim.sf, 0, 60, limit=400, epsabs=1e-13)[0]
        m2 = (
            2
            * integrate.quad(lambda v: v * lim.sf(v), 0, 60, limit=400, epsabs=1e-13)[0]
        )
        agreed(m1, mean_ref, tol=1e-9)
        agreed(m2 - m1**2, var_ref, tol=1e-8)
        print(
            f"   {fam}: mean and variance from the inverted CDF: {m1:.10f}, {m2 - m1**2:.10f}"
            f" (analytic {mean_ref:.10f}, {var_ref:.10f})"
        )

    def ad_sf(z):
        return ad1.sf(z)

    def cvm_sf(x):
        return cvm1.sf(x)

    # ================================================================= section 3
    sec = "GOLD 3: LHD single machine (KL2020 Section 8.1, Tables 1-3)"
    t = np.asarray(LHD_TIMES_H, dtype=float)
    tau = LHD_TAU_H
    x = gaps(t)
    n = len(t)
    print(
        f"\n   N = {n}, sum of gaps = {x.sum():.0f}, censored gap = {tau - t[-1]:.0f}"
    )
    reversal_selftest(t, tau, (0.3, 0.5))
    hyd_rev = segments_hydraulic(HYDRAULIC_KK1992)
    tau_max_rev = max(s["tau"] for s in hyd_rev)
    for s in hyd_rev:
        reversal_selftest(
            s["t"],
            s["tau"],
            (0.5, 1 - 0.5 * s["tau"] / tau_max_rev, 0.5 * tau_max_rev / s["tau"]),
        )
    print(
        "   time reversal (LR -> -LR; KS, CvM, AD unchanged; ELR(a) -> ELR(1-a)) holds for"
        " both implementations on LHD and the six machines"
    )
    row(sec, "N (failure times in Table 1)", n, "36", None)
    row(sec, "sum of the 36 gaps, h", x.sum(), "1970", None)
    row(sec, "censored gap tau - T_N, h", tau - t[-1], "30", None)
    mu1, sd1, g1 = moments(x, 1)
    mu0, sd0, g0 = moments(x, 0)
    row(sec, "mu_hat", mu1, "54.722", "54.72")
    row(sec, "sigma_hat, 1/(N-1)", sd1, "48.611", "48.61")
    row(sec, "sigma_hat, 1/N", sd0, "47.931", "48.61")
    row(sec, "gamma_hat, 1/(N-1)", g1, "0.8883", "0.888")
    row(sec, "gamma_hat, 1/N", g0, "0.8759", "0.888")
    mt, st, gt = eq10(x, tau)
    row(sec, "mu_tilde, eq (10)", mt, "55.556", "55.56")
    row(sec, "sigma_tilde, eq (10)", st, "47.228", "47.23")
    row(sec, "gamma_tilde, eq (10)", gt, "0.8501", "0.850")
    wb = weibull_mle(x, np.array([tau - t[-1]]))
    agreed(wb["shape"], wb["shape_direct"], tol=1e-6)
    agreed(wb["scale"], wb["scale_direct"], tol=1e-6)
    row(sec, "Weibull mu (30 h gap right-censored)", wb["mu"], "55.459", "55.46")
    row(sec, "Weibull sigma", wb["sigma"], "47.219", "47.22")
    row(sec, "Weibull gamma", wb["gamma"], "0.8514", "0.851")
    row(sec, "Weibull shape", wb["shape"], "1.179", None, f"scale {wb['scale']:.4f}")
    lap = agreed(lr_eq4(t, tau, 1.0), lr_direct(t, tau, 1.0))
    row(sec, "Laplace (eq (4) without 1/gamma)", lap, "0.6051", "0.605")
    lr1 = agreed(lr_eq4(t, tau, g1), lr_direct(t, tau, g1))
    lr0 = agreed(lr_eq4(t, tau, g0), lr_direct(t, tau, g0))
    row(sec, "LR, 1/(N-1)", lr1, "0.6811", "0.681")
    row(sec, "LR p, 1/(N-1)", p2(lr1), "0.4958", "0.50")
    row(sec, "LR, 1/N", lr0, "0.6908", "0.681")
    ss = eq11(x)
    gs = ss / mu1
    row(sec, "sigma*, eq (11)", ss, "42.771", "42.77")
    row(sec, "gamma* = sigma*/mu_hat", gs, "0.7816", "0.782")
    lrs = lap / gs
    row(sec, "LR with gamma*", lrs, "0.7741", "0.774")
    row(sec, "LR with gamma*, p", p2(lrs), "0.4389", "0.44")

    sec = "GOLD 3: KL2020 Table 3, p. 113 (two-sided p-values)"
    gold_t3 = {  # test: (p 1/(N-1), p 1/N, statistic 1/(N-1), printed p)
        "LR": ("0.4958", "0.4897", "0.6811", "0.50"),
        "KS": ("0.2864", "0.2711", "0.9850", "0.29"),
        "CvM": ("0.1312", "0.1240", "0.30462", "0.13"),
        "AD": ("0.0856", "0.0795", "2.05555", "0.086"),
        "ELR": ("0.0115", "0.0103", "2.5283", "0.011"),
    }
    printed_eq6_as_is = None
    t3_p: dict[int, dict[str, float]] = {}
    for ddof, gg in ((1, g1), (0, g0)):
        cq, aq = quad_check(t, tau, gg)
        st_ = {
            "LR": agreed(lr_eq4(t, tau, gg), lr_direct(t, tau, gg)),
            "KS": agreed(ks_eq5(t, tau, gg), ks_direct(t, tau, gg)),
            "CvM": agreed(cvm_eq6(t, tau, gg), cvm_direct(t, tau, gg), cq),
            "AD": agreed(ad_eq7(t, tau, gg), ad_direct(t, tau, gg), aq),
            "ELR": agreed(elr_eq9(t, tau, gg, 0.5), elr_direct(t, tau, gg, 0.5)),
        }
        rel_exact = max(
            abs(a_ - b_) / abs(b_)
            for a_, b_ in (
                (lr_eq4(t, tau, gg), lr_direct(t, tau, gg)),
                (ks_eq5(t, tau, gg), ks_direct(t, tau, gg)),
                (cvm_eq6(t, tau, gg), cvm_direct(t, tau, gg)),
                (ad_eq7(t, tau, gg), ad_direct(t, tau, gg)),
                (elr_eq9(t, tau, gg, 0.5), elr_direct(t, tau, gg, 0.5)),
            )
        )
        rel_quad = max(
            abs(cvm_eq6(t, tau, gg) - cq) / cq, abs(ad_eq7(t, tau, gg) - aq) / aq
        )
        print(
            f"   LHD ddof {ddof}: printed-equation closed forms vs exact piecewise integral of"
            f" the bridge, largest relative difference {rel_exact:.1e}; CvM and AD vs adaptive"
            f" quadrature (epsrel 1e-13) {rel_quad:.1e}"
        )
        if ddof == 1:
            row(
                sec,
                "INFO closed forms vs exact piecewise integral, max rel diff",
                rel_exact,
                None,
                None,
                note="GOLD: closed forms agree with numerical integration to 1e-14",
            )
            row(
                sec,
                "INFO CvM, AD closed forms vs adaptive quadrature, max rel diff",
                rel_quad,
                None,
                None,
                note="quadrature requested at epsrel 1e-13",
            )
        pv = {
            "LR": p2(st_["LR"]),
            "KS": kolmogorov_sf(st_["KS"]),
            "CvM": cvm_sf(st_["CvM"]),
            "AD": ad_sf(st_["AD"]),
            "ELR": p2(st_["ELR"]),
        }
        t3_p[ddof] = pv
        tag = "1/(N-1)" if ddof == 1 else "1/N"
        for k in gold_t3:
            gp1, gp0, gs1, pr = gold_t3[k]
            row(sec, f"{k} p, {tag}", pv[k], gp1 if ddof == 1 else gp0, pr)
            row(sec, f"{k} statistic, {tag}", st_[k], gs1 if ddof == 1 else None, None)
        if ddof == 1:
            printed_eq6_as_is = cvm_eq6(t, tau, gg, middle_power=1)
    for use_order in (True, False):
        z = linrank_klcode(t, use_order)
        row(
            sec,
            f"INFO LinR p, scores via {'order()' if use_order else 'ordinal ranks'}",
            p2(z),
            None,
            "0.76",
            note=f"LinR = {z:.4f}; not recomputed in GOLD",
        )
    print(
        "\n   eq (6) evaluated exactly as printed (middle term / tau):"
        f" {printed_eq6_as_is:.6g}; with / tau^2 it equals the bridge integral."
    )

    # ================================================================= section 4.1
    sec = "GOLD 4.1: hydraulic, six LHD machines (KL2020 Section 8.2, Table 4)"
    hyd = segments_hydraulic(HYDRAULIC_KK1992)
    gold_m = {
        1: ("22", "2496", "1.1952", "1.6538"),
        3: ("24", "3526", "1.0864", "0.6848"),
        9: ("26", "4743", "1.1075", "1.9106"),
        11: ("27", "2913", "0.9277", "0.2606"),
        17: ("25", "3230", "1.0102", "1.5297"),
        20: ("22", "3309", "0.8672", "0.1102"),
    }
    for s in hyd:
        gN, gtau, gg, glr = gold_m[s["id"]]
        g_1 = moments(s["x"], 1)[2]
        row(sec, f"machine {s['id']}: N_j", len(s["t"]), gN)
        row(sec, f"machine {s['id']}: tau_j", s["tau"], gtau)
        row(
            sec,
            f"machine {s['id']}: gamma_j, 1/(N-1)",
            g_1,
            gg,
            note=f"1/N: {moments(s['x'], 0)[2]:.4f}",
        )
        row(
            sec,
            f"machine {s['id']}: LR_j, 1/(N-1)",
            lr_eq4(s["t"], s["tau"], g_1),
            glr,
            note=f"1/N: {lr_eq4(s['t'], s['tau'], moments(s['x'], 0)[2]):.4f}",
        )
    lrm_p: dict[tuple[str, int], float] = {}
    for wname, ddof, gstat, gp in (
        ("KL2020", 1, "2.349", "0.0188"),
        ("KL2020", 0, "2.399", "0.0164"),
        ("KLv1", 1, "2.683", "0.0073"),
        ("KLv1", 0, "2.739", "0.0062"),
    ):
        v = lr_m(hyd, wname, ddof)
        lrm_p[(wname, ddof)] = p2(v)
        tag = "1/(N-1)" if ddof else "1/N"
        row(sec, f"LR^m, {wname} weights, {tag}", v, gstat)
        row(sec, f"LR^m p, {wname} weights, {tag}", p2(v), gp, "0.019")
    for ddof_lr, ddof_w in ((0, 1), (1, 0)):
        v = lr_m_mixed(hyd, ddof_lr, ddof_w)
        row(
            sec,
            f"INFO LR^m KL2020 weights, LR_j ddof {ddof_lr}, A_j ddof {ddof_w}",
            v,
            None,
            None,
            note=f"p {p2(v):.4f}; ddof 1 = 1/(N-1). KLcode's weights are always ddof 1",
        )
    glv = gl(hyd)
    row(sec, "GL", glv, "1.864")
    row(sec, "GL p", p2(glv), "0.0623", "0.062")

    # ELR^m: eq (15) weights; a_j by the authors' code or by the Section 8.2 text.
    tau_all = np.array([s["tau"] for s in hyd])
    elrm_p: dict[tuple[bool, int], float] = {}
    for ddof in (1, 0):
        g_all = np.array([moments(s["x"], ddof)[2] for s in hyd])
        n_all = np.array([len(s["t"]) for s in hyd], dtype=float)
        a15 = np.sqrt(n_all) / g_all
        a15 /= np.sqrt(np.sum(a15**2))
        tag = "1/(N-1)" if ddof else "1/N"
        for conv, avec, gold_p in (
            (
                "a_j = 1 - tau_j/(2 max tau) (KLcode)",
                1 - 0.5 * tau_all / tau_all.max(),
                "0.0027",
            ),
            (
                "a_j tau_j = 4743/2 (Section 8.2 text)",
                (tau_all.max() / 2) / tau_all,
                "0.0079",
            ),
        ):
            e = [
                agreed(
                    elr_eq9(s["t"], s["tau"], g, aj),
                    elr_direct(s["t"], s["tau"], g, aj),
                )
                for s, g, aj in zip(hyd, g_all, avec)
            ]
            v = float(np.sum(a15 * np.array(e)))
            elrm_p[("KLcode" in conv, ddof)] = p2(v)
            row(
                sec,
                f"ELR^m p, {conv}, {tag}",
                p2(v),
                gold_p if ddof else None,
                "0.003",
                note=f"ELR^m = {v:.4f}",
            )

    # AD^m, CvM^m: weights prop. tau_j, per-machine gamma_j (KLcode sigma = "s"). The null
    # is sum_j w_j X_j with X_j iid from the limit; it does not depend on the divisor, so one
    # set of draws serves both divisors.
    w_tau = tau_all / np.sqrt(np.sum(tau_all**2))
    ad_lim = LimitSum("ad", w_tau)
    cvm_lim = LimitSum("cvm", w_tau)
    rng = np.random.default_rng(args.seed)
    stats_m = {}
    for ddof in (1, 0):
        g_all = [moments(s["x"], ddof)[2] for s in hyd]
        stats_m[("AD^m", ddof)] = float(
            sum(
                w * agreed(ad_eq7(s["t"], s["tau"], g), ad_direct(s["t"], s["tau"], g))
                for w, s, g in zip(w_tau, hyd, g_all)
            )
        )
        stats_m[("CvM^m", ddof)] = float(
            sum(
                w
                * agreed(cvm_eq6(s["t"], s["tau"], g), cvm_direct(s["t"], s["tau"], g))
                for w, s, g in zip(w_tau, hyd, g_all)
            )
        )
    # GOLD 4.1 prints the exact value and one 4e5-draw simulation with its s.e.
    gold_mc = {
        "AD^m": ("0.0752", "0.0750", "0.076", 0.0004),
        "CvM^m": ("0.0651", "0.0650", "0.064", 0.0004),
    }
    mc_exact: dict[tuple[str, int], float] = {}
    for name, lim, fam in (("AD^m", ad_lim, "ad"), ("CvM^m", cvm_lim, "cvm")):
        gold_exact, gold_p, pr, gold_se = gold_mc[name]
        for ddof in (1, 0):
            tag = "1/(N-1)" if ddof else "1/N"
            stat = stats_m[(name, ddof)]
            exact = lim.sf(stat)
            mc_exact[(name, ddof)] = exact
            row(sec, f"{name} statistic, {tag}", stat, None, None)
            row(
                sec,
                f"{name} p exact (Imhof), {tag}",
                exact,
                gold_exact if ddof else None,
                pr,
            )
        for k_run, n_draws in enumerate((args.draws, args.big_draws)):
            null = sum(w * sample_limit(fam, rng, n_draws) for w in w_tau)
            for ddof in (1, 0):
                tag = "1/(N-1)" if ddof else "1/N"
                stat = stats_m[(name, ddof)]
                exact = lim.sf(stat)
                p_mc = float(np.mean(null > stat))
                se = math.sqrt(p_mc * (1 - p_mc) / n_draws)
                row(
                    sec,
                    f"{name} p MC {n_draws:.0e} draws, {tag}",
                    p_mc,
                    gold_p if ddof and k_run == 0 else None,
                    pr,
                    note=f"MC s.e. {se:.5f}; (MC - exact)/se = {(p_mc - exact) / se:+.2f}",
                    # At the default seed and draws the run reproduces GOLD's value, so
                    # the Monte Carlo allowance applies only when either is changed.
                    mc_tol=(
                        3 * math.hypot(se, gold_se)
                        if ddof and k_run == 0 and not default_mc
                        else None
                    ),
                )
            # the simulated null against the exact law: a check on the sampler itself
            for q in (0.5, 0.9, 0.95):
                x_q = float(np.quantile(null, q))
                agreed_mc = abs(lim.cdf(x_q) - q) <= 4 * math.sqrt(
                    q * (1 - q) / n_draws
                )
                if not agreed_mc:
                    raise AssertionError(
                        f"{name} sampler off at q = {q}: {lim.cdf(x_q)}"
                    )
        del null
    for name, exact in (
        ("AD^m", ad_lim.sf(stats_m[("AD^m", 1)])),
        ("CvM^m", cvm_lim.sf(stats_m[("CvM^m", 1)])),
    ):
        row(
            sec,
            f"{name}: MC s.e. of p at KLcode's Npsim = 1e4",
            math.sqrt(exact * (1 - exact) / 1e4),
            "0.0025" if name == "CvM^m" else None,
            None,
            note="GOLD: 's.e. about 0.0025'",
        )
        row(
            sec,
            f"{name}: MC s.e. of p at 4e5 draws",
            math.sqrt(exact * (1 - exact) / 4e5),
            "0.0004",
            None,
        )

    # ================================================================= section 4.2
    sec = "GOLD 4.2: small bowel motility (KLv1 Section 6.2)"
    bow = segments_bowel(BOWEL_AH1991)
    pooled = np.concatenate([s["x"] for s in bow])
    cens = np.array([s["cens"] for s in bow])
    print(
        f"\n   bowel: {len(bow)} subjects, {len(pooled)} complete periods, {len(cens)} censored"
    )
    bm1, bs1, bg1 = moments(pooled, 1)
    bm0, bs0, bg0 = moments(pooled, 0)
    row(sec, "mu_hat (80 complete periods)", bm1, "98.76", "98.76")
    row(sec, "sigma_hat, 1/(N-1)", bs1, "52.62", "52.62")
    row(sec, "sigma_hat, 1/N", bs0, "52.29", "52.62")
    row(sec, "gamma_hat, 1/(N-1)", bg1, "0.5328", "0.533")
    row(sec, "gamma_hat, 1/N", bg0, "0.5294", "0.533")
    u = np.array([np.sum(s["t"]) - len(s["t"]) * s["tau"] / 2 for s in bow])
    nb = np.array([len(s["t"]) for s in bow], dtype=float)
    tb = np.array([s["tau"] for s in bow])
    lapm = SQRT12 * np.sum(u) / math.sqrt(np.sum(tb**2 * nb))  # KLv1 eq (16), gamma = 1
    row(sec, "Laplace^m (KLv1 eq (16), gamma = 1)", lapm, "1.9530", "1.95")
    row(sec, "Laplace^m p", p2(lapm), "0.0508", "0.051")
    row(sec, "LR^m = Laplace^m / gamma_hat, 1/(N-1)", lapm / bg1, "3.6658", "3.67")
    row(sec, "LR^m p, 1/(N-1)", p2(lapm / bg1), "0.000247", "0.00024")
    row(sec, "LR^m, 1/N", lapm / bg0, "3.6890", "3.67")
    row(sec, "2 Phi(-3.67), the printed rounded statistic", p2(3.67), None, "0.00024")
    lap_kl2020 = SQRT12 * np.sum(u / tb) / math.sqrt(np.sum(nb))
    row(
        sec,
        "INFO: KL2020 eq (16), common gamma, 1/(N-1)",
        lap_kl2020 / bg1,
        None,
        None,
        note="not in GOLD; the journal weighting with one pooled gamma",
    )
    glb = gl(bow)
    row(sec, "GL", glb, "2.685")
    row(sec, "GL p", p2(glb), "0.0073", "0.007")
    wbb = weibull_mle(pooled, cens)
    agreed(wbb["shape"], wbb["shape_direct"], tol=1e-6)
    agreed(wbb["scale"], wbb["scale_direct"], tol=1e-6)
    row(sec, "Weibull mu (19 censored periods)", wbb["mu"], "104.50", "104.49")
    row(sec, "Weibull sigma", wbb["sigma"], "52.46", "52.45")
    row(sec, "Weibull gamma", wbb["gamma"], "0.502", "0.502")
    row(sec, "Weibull shape", wbb["shape"], "2.092", None, f"scale {wbb['scale']:.4f}")
    no5 = [s for s in bow if len(s["t"]) >= 2]
    for wname, ddof in (("KL2020", 1), ("KL2020", 0), ("KLv1", 1), ("KLv1", 0)):
        tag = "1/(N-1)" if ddof else "1/N"
        row(
            sec,
            f"per-segment gamma, no subject 5, {wname}, {tag}",
            lr_m(no5, wname, ddof),
            {("KLv1", 0): "4.37", ("KL2020", 1): "4.49"}.get((wname, ddof)),
            None,
            note=f"{len(no5)} subjects"
            + {
                ("KLv1", 0): "; KLv1 weights with 1/N",
                ("KL2020", 1): "; the shipped defaults",
            }.get((wname, ddof), ""),
        )

    # ================================================================= section 3, "The divisor"
    # GOLD: "Every printed number (Table 2, Section 8.1, Table 3, Table 4, and KLv1 Section
    # 6.2) lands on 1/(N-1) and misses on 1/N." Only the printed values that depend on the
    # divisor are tallied. A deterministic value lands when it rounds to the printed digits;
    # the two Table 4 values the authors simulated (Npsim = 1e4) land when within 3 of their
    # Monte Carlo standard errors; 0.00024 is judged from the statistic rounded to two
    # decimals, which is how GOLD 4.2 says the authors computed it.
    def lands(value: float, printed: str, npsim: float | None) -> bool:
        if npsim is None:
            return abs(value - float(printed)) <= 0.5 * 10.0 ** (
                -_decimals(printed)
            ) * (1 + 1e-9)
        return abs(value - float(printed)) <= 3 * math.sqrt(value * (1 - value) / npsim)

    divisor_entries = [
        ("Table 2 sigma_hat", sd1, sd0, "48.61", None),
        ("Table 2 gamma_hat", g1, g0, "0.888", None),
        ("Section 8.1 LR", lr1, lr0, "0.681", None),
        *[
            (f"Table 3 {k} p", t3_p[1][k], t3_p[0][k], gold_t3[k][3], None)
            for k in gold_t3
        ],
        (
            "Table 4 LR^m p, eq (15) weights",
            lrm_p[("KL2020", 1)],
            lrm_p[("KL2020", 0)],
            "0.019",
            None,
        ),
        (
            "Table 4 ELR^m p, KLcode a_j",
            elrm_p[(True, 1)],
            elrm_p[(True, 0)],
            "0.003",
            None,
        ),
        (
            "Table 4 AD^m p, simulated by the authors",
            mc_exact[("AD^m", 1)],
            mc_exact[("AD^m", 0)],
            "0.076",
            1e4,
        ),
        (
            "Table 4 CvM^m p, simulated by the authors",
            mc_exact[("CvM^m", 1)],
            mc_exact[("CvM^m", 0)],
            "0.064",
            1e4,
        ),
        ("KLv1 6.2 sigma_hat", bs1, bs0, "52.62", None),
        ("KLv1 6.2 gamma_hat", bg1, bg0, "0.533", None),
        ("KLv1 6.2 LR^m", lapm / bg1, lapm / bg0, "3.67", None),
        (
            "KLv1 6.2 LR^m p, from the statistic to 2 dp",
            p2(round(lapm / bg1, 2)),
            p2(round(lapm / bg0, 2)),
            "0.00024",
            None,
        ),
    ]
    sec = "GOLD 3: the divisor claim"
    n_ok = 0
    for label, v1, v0, pr, npsim in divisor_entries:
        l1, l0 = lands(v1, pr, npsim), lands(v0, pr, npsim)
        n_ok += l1 and not l0
        row(
            sec,
            f"INFO {label}",
            v1,
            None,
            None,
            note=(
                f"printed {pr}: 1/(N-1) {'lands' if l1 else 'MISSES'}; 1/N {fmt(v0)}"
                f" {'LANDS' if l0 else 'misses'}"
                + ("; within 3 MC s.e. at Npsim 1e4" if npsim else "")
            ),
        )
    row(
        sec,
        "divisor-dependent printed values: 1/(N-1) lands, 1/N misses",
        n_ok,
        str(len(divisor_entries)),
        None,
        note="GOLD column = the count 'every' implies",
    )

    # ================================================================= section 5
    sec = "GOLD 5.1: limiting A^2"
    s1974 = {0.15: "1.610", 0.10: "1.933", 0.05: "2.492", 0.025: "3.070", 0.01: "3.857"}
    gold_q = {
        0.15: "1.6212385",
        0.10: "1.9329578",
        0.05: "2.4923672",
        0.025: "3.0774642",
        0.01: "3.8781250",
    }
    gold_tail = {
        0.15: "0.1523",
        0.10: "0.1000",
        0.05: "0.0500",
        0.025: "0.0252",
        0.01: "0.0102",
    }
    mm2004 = {
        0.10: "1.9329578327415937304",
        0.05: "2.4923671600494096176",
        0.01: "3.8781250216053948842",
    }
    for alpha in (0.15, 0.10, 0.05, 0.025, 0.01):
        q = ad1.isf(alpha)
        row(
            sec,
            f"{alpha:.3g} upper point, exact limit",
            q,
            gold_q[alpha],
            s1974[alpha],
            note=(
                f"MM2004 as quoted by GOLD: {mm2004[alpha]}; diff {q - float(mm2004[alpha]):+.1e}"
                if alpha in mm2004
                else ""
            ),
        )
        row(
            sec,
            f"upper tail at S1974 {s1974[alpha]}",
            ad_sf(float(s1974[alpha])),
            gold_tail[alpha],
            f"{alpha:g}",
        )
    grid = np.concatenate([np.linspace(0.05, 2.0, 400), np.linspace(2.0, 8.0, 241)])
    errs = np.array([adinf_repo(z) - ad1.cdf(z) for z in grid])
    for z0, gold_err in (
        (1.0, "-1.9e-5"),
        (1.933, "-1.1e-5"),
        (2.492, "+8.1e-6"),
        (3.878, "-2.6e-6"),
        (5.0, "+8.4e-6"),
    ):
        row(
            sec,
            f"adinf(z) - exact CDF at z = {z0}",
            adinf_repo(z0) - ad1.cdf(z0),
            gold_err,
        )
    j = int(np.argmax(np.abs(errs[grid < 2])))
    row(
        sec,
        "max |adinf error| on 0.05 <= z < 2",
        float(np.abs(errs[grid < 2]).max()),
        None,
        "2e-6",
        note=f"at z = {grid[grid < 2][j]:.3f}; printed = MM2004 p. 3 as quoted by GOLD",
    )
    j2 = int(np.argmax(np.abs(errs[grid >= 2])))
    row(
        sec,
        "max |adinf error| on 2 <= z <= 8",
        float(np.abs(errs[grid >= 2]).max()),
        None,
        "8e-7",
        note=f"at z = {grid[grid >= 2][j2]:.3f}",
    )
    for lo_z, hi_z in ((0.6, 1.6), (2.0, 4.0), (4.0, 7.0)):
        m = (grid >= lo_z) & (grid <= hi_z)
        jj = int(np.argmax(np.abs(errs[m])))
        zz = optimize.minimize_scalar(
            lambda z: -abs(adinf_repo(z) - ad1.cdf(z)),
            bounds=(max(lo_z, grid[m][jj] - 0.05), min(hi_z, grid[m][jj] + 0.05)),
            method="bounded",
            options={"xatol": 1e-6},
        ).x
        row(
            sec,
            f"INFO extreme adinf error on [{lo_z}, {hi_z}]",
            adinf_repo(zz) - ad1.cdf(zz),
            None,
            None,
            note=f"at z = {zz:.4f}; GOLD quotes point values, not extrema",
        )
    # Not GOLD: the numbers in the docstring of c2_anderson_darling.ad_limiting_cdf, in the
    # GOLD column so the same half-unit rule applies.
    sec = "REPO (not GOLD): c2_anderson_darling.ad_limiting_cdf docstring"
    row(sec, "adinf upper tail at 1.933", 1 - adinf_repo(1.933), "0.10001")
    row(sec, "adinf upper tail at 2.492", 1 - adinf_repo(2.492), "0.05001")
    row(sec, "adinf upper tail at 3.070", 1 - adinf_repo(3.070), "0.0252")
    row(sec, "adinf upper tail at 3.857", 1 - adinf_repo(3.857), "0.0102")
    row(sec, "exact 2.5% point", ad1.isf(0.025), "3.0775")
    row(sec, "exact 1% point", ad1.isf(0.01), "3.8781")
    zz = optimize.minimize_scalar(
        lambda z: -abs(adinf_repo(z) - ad1.cdf(z)),
        bounds=(0.9, 1.05),
        method="bounded",
        options={"xatol": 1e-6},
    ).x
    row(
        sec,
        "largest |error| below z = 2",
        abs(adinf_repo(zz) - ad1.cdf(zz)),
        "1.95e-5",
    )
    row(sec, "where it occurs, z", zz, "0.97")
    zz_hi = optimize.minimize_scalar(
        lambda z: -abs(adinf_repo(z) - ad1.cdf(z)),
        bounds=(2.3, 2.9),
        method="bounded",
        options={"xatol": 1e-6},
    ).x
    row(
        sec,
        "largest |error| on z >= 2 (grid 2..8)",
        float(np.abs(errs[grid >= 2]).max()),
        "9.1e-6",
        note=f"refined extremum {abs(adinf_repo(zz_hi) - ad1.cdf(zz_hi)):.4e}",
    )
    row(sec, "where it occurs, z", zz_hi, "2.59")

    sec = "GOLD 5.2: limiting W^2"
    ad1952 = {
        0.10: "0.34730",
        0.05: "0.46136",
        0.025: None,
        0.01: "0.74346",
        0.001: "1.16786",
    }
    s_w2 = {0.10: "0.347", 0.05: "0.461", 0.025: "0.581", 0.01: "0.743", 0.001: None}
    gold_series = {
        0.10: "0.3473049",
        0.05: "0.4613613",
        0.025: "0.5806147",
        0.01: "0.7434593",
        0.001: "1.1678583",
    }
    for alpha in (0.10, 0.05, 0.025, 0.01, 0.001):
        q = cvm1.isf(alpha)
        agreed(1 - cvm_cdf_bessel(q), alpha, tol=1e-11)
        row(
            sec,
            f"{alpha:g} upper point, exact limit",
            q,
            gold_series[alpha],
            ad1952[alpha] or s_w2[alpha],
            note=(
                "printed = AD1952 Table 1 as quoted by GOLD (not in spec/sources); S1974 "
                + (s_w2[alpha] or "-")
            )
            if ad1952[alpha]
            else "printed = S1974 modified W^2",
        )
        if alpha == 0.025:
            row(
                sec,
                "2.5% point to five decimals ('0.58061 is correct')",
                round(q, 5),
                "0.58061",
                None,
                note=f"unrounded {q:.10f}",
            )

    n_bad, n_repo = print_rows()
    print(
        f"\n{n_bad} GOLD row(s) disagree with GOLD_STANDARD beyond half a unit in its last"
        f" digit; {n_repo} REPO docstring row(s) disagree with the docstring."
    )
    missing = strings_missing_from(args.gold, args.repo_docstring)
    for where, label, text in missing:
        print(f"   MISSING from {where}: {text!r} ({label})")
    print(f"{len(missing)} GOLD or REPO string(s) missing from the file they quote.")
    if args.fixture is not None:
        print(f"{n_fixture} fixture comparison(s) failed.")
    if args.json is not None:
        import json

        args.json.write_text(
            json.dumps(
                [
                    {
                        "section": r[0],
                        "label": r[1],
                        "value": r[2],
                        "gold": r[3],
                        "printed": r[4],
                        "note": r[5],
                        "vs_gold": verdict(r[2], r[3], r[6]),
                        "vs_printed": verdict(r[2], r[4]),
                    }
                    for r in ROWS
                ],
                indent=1,
            )
        )
    return 1 if n_bad or n_repo or missing or n_fixture else 0


if __name__ == "__main__":
    sys.exit(main())
