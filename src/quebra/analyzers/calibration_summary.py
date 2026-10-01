"""Turning the bench tables into the four typed artifacts the calibration figures draw.

These are steps: pure compute, no disk, no matplotlib. The bench tables arrive as
DataFrames that the job loaded as datasets, so nothing here imports `jobs/bench/` - the
dependency runs through provenance (their sha256 enters the run identity) rather than
around it.

Four questions, four artifacts, one per figure:

    SizeVsN            does a check hold its nominal level at the event counts we have?
    PowerVsDependence  how does its power grow with the dependence it is shown?
    ValidationCurve    is the transcription right, and what does estimating gamma cost?
    ReadDependence     does read-level correlation reach the durations at all?

`ValidationCurve` is the one that is not a reshape of a table: it re-simulates, because a
P-P plot needs the statistics themselves and the tables only carry rejection rates. It is
a deterministic function of `(seed, n_values, n_replicates)`, and the seed travels as a
step kwarg so it lands on the provenance label.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

import quebra.analyzers.checks.c2_anderson_darling as c2
from quebra.analyzers.checks.result import CALIB_ASYMPTOTIC, CALIB_PERMUTATION, Segment
from quebra.core._artifact_guard import StaleArtifactGuard

# Nominal significance every size curve is read against.
NOMINAL_ALPHA = 0.05


def _label(check: str, calibration: str, variant: str) -> str:
    suffix = f"/{variant}" if isinstance(variant, str) and variant else ""
    return f"{check} [{calibration}{suffix}]"


def _with_label(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["row"] = [
        _label(c, k, v)
        for c, k, v in zip(frame["check"], frame["calibration"], frame["variant"])
    ]
    return frame


@dataclass
class SizeVsN(StaleArtifactGuard):
    """Empirical size against event count, on the primary null configuration.

    `n_values` is shared; `rate_by_row` and `se_by_row` are keyed by the row label. The
    Monte Carlo SE travels WITH the rate because a size table without its own error bar
    invites exactly the over-reading the bench exists to prevent.
    """

    n_values: list[int]
    rate_by_row: dict[str, list[float]] = field(default_factory=dict)
    se_by_row: dict[str, list[float]] = field(default_factory=dict)
    alpha: float = NOMINAL_ALPHA
    cell_description: str = ""
    # `rate_by_row` is a MEAN over the cells that remain after the primary-cell filter -
    # in practice the two Weibull shapes, which are not interchangeable: at n = 20 C2
    # asymptotic is 0.0575 at shape 0.75 and 0.0790 at shape 1.5, a spread wider than the
    # effect the large-n end of the curve shows. So the spread travels with the mean and
    # the figure draws it; a single line labelled "the primary null cell" would have been
    # the kind of silent averaging this project keeps finding.
    n_cells_pooled: int = 0
    pooled_over: list[str] = field(default_factory=list)
    rate_lo_by_row: dict[str, list[float]] = field(default_factory=dict)
    rate_hi_by_row: dict[str, list[float]] = field(default_factory=dict)


@dataclass
class PowerVsDependence(StaleArtifactGuard):
    """Power against the dependence actually induced, not against the latent knob.

    `induced_lag1` is the x-axis on purpose: `rho` is a Gaussian-copula parameter no reader
    can interpret, while the induced duration-level lag-1 is the dependence a check sees.
    """

    induced_lag1_by_n: dict[int, list[float]] = field(default_factory=dict)
    power_by_row_and_n: dict[str, dict[int, list[float]]] = field(default_factory=dict)
    alpha: float = NOMINAL_ALPHA
    # One representative row per check, chosen HERE rather than at draw time. The artifact
    # keeps every row (completeness); this says which the figure shows and why.
    rows_drawn: list[str] = field(default_factory=list)
    selection_reason: str = ""


@dataclass
class ValidationCurve(StaleArtifactGuard):
    """P-P data for eq (7): theoretical vs empirical CDF, both gamma paths.

    Eq (7) is Kvaloy and Lindqvist's Anderson-Darling statistic, Technometrics 62(1) 2020,
    p. 104. Two curves per n. `gamma = 1` isolates the transcription; `gamma_hat` is what ships.
    The vertical gap between them at small n IS the finite-N cost of estimating gamma, and
    it is the reason the transcription test must not be read as covering production.
    """

    n_values: list[int] = field(default_factory=list)
    theoretical: list[float] = field(default_factory=list)
    empirical_gamma_one: dict[int, list[float]] = field(default_factory=dict)
    empirical_gamma_hat: dict[int, list[float]] = field(default_factory=dict)
    n_replicates: int = 0
    seed: int = 0


@dataclass
class ReadDependence(StaleArtifactGuard):
    """Induced duration-level dependence against read-level correlation.

    The finding this artifact exists to carry is a NULL one: the curve stays flat near zero
    while `y = x` climbs away from it. Read correlation does not reach the durations.
    """

    read_rho: list[float] = field(default_factory=list)
    induced_lag1: list[float] = field(default_factory=list)
    # Standard error of the MEAN ACROSS the pooled cells (the n and clock levels at this
    # rho), not the sampling error of any single lag-1 estimate. It answers "how much do
    # the cells disagree", which is the right question for a curve drawn through their
    # mean - but it is not a confidence interval on the statistic, and naming it so here
    # keeps the figure's error bars from being read as one.
    induced_lag1_se: list[float] = field(default_factory=list)
    arm: str = ""


def size_vs_n(
    size_table: pd.DataFrame,
    *,
    arm: str = "A_iid_weibull",
    alpha: float = NOMINAL_ALPHA,
) -> SizeVsN:
    """The primary null configuration: one arm, continuous durations, no extra censoring.

    Narrowed rather than pooled over the whole grid, because the figure's claim is about a
    check's behaviour at an event count and averaging over censoring would flatten the very
    effect the asymptotic curves exist to show. What the filter does NOT remove is the
    shape factor - `pooled_over` names whatever survives it and `rate_lo/hi_by_row` carry
    the spread, so the figure can draw the pooling instead of hiding it.
    """
    frame = _with_label(size_table)
    cell = frame[
        (frame["arm"] == arm)
        & (~frame["quantised"].astype(bool))
        & (frame["censoring_target"] == 0.0)
    ]
    if not len(cell):
        raise ValueError(f"no {arm} continuous, uncensored size cells in the table")
    n_values = sorted(int(n) for n in cell["n_target"].unique())
    # Which factors still vary after the filter. Anything listed here is averaged over,
    # and saying so is the whole difference between "the primary null cell" and a mean.
    pooled_over = [
        column
        for column in ("shape", "quantised", "censoring_target")
        if cell[column].nunique() > 1
    ]
    rate_by_row: dict[str, list[float]] = {}
    se_by_row: dict[str, list[float]] = {}
    lo_by_row: dict[str, list[float]] = {}
    hi_by_row: dict[str, list[float]] = {}
    for row_label, group in cell.groupby("row"):
        by_n = group.groupby("n_target")
        rate_by_row[row_label] = [
            float(by_n.get_group(n)["rejection_rate"].mean()) for n in n_values
        ]
        lo_by_row[row_label] = [
            float(by_n.get_group(n)["rejection_rate"].min()) for n in n_values
        ]
        hi_by_row[row_label] = [
            float(by_n.get_group(n)["rejection_rate"].max()) for n in n_values
        ]
        # The NULL standard error, not the observed rate's: the observed one collapses to
        # zero on a cell that never rejected, drawing a conservative cell as infinitely
        # precise. `n_used` is averaged over the pooled cells, exact when they share a
        # replicate count, which they do here.
        se_by_row[row_label] = [
            float(np.sqrt(alpha * (1.0 - alpha) / by_n.get_group(n)["n_used"].mean()))
            for n in n_values
        ]
    n_pooled = int(round(len(cell) / max(len(rate_by_row), 1) / max(len(n_values), 1)))
    description = f"{arm}, continuous durations, censoring 0"
    if pooled_over:
        description += f" - mean over {n_pooled} cells varying {', '.join(pooled_over)}"
    return SizeVsN(
        n_values=n_values,
        rate_by_row=rate_by_row,
        se_by_row=se_by_row,
        alpha=alpha,
        cell_description=description,
        n_cells_pooled=n_pooled,
        pooled_over=pooled_over,
        rate_lo_by_row=lo_by_row,
        rate_hi_by_row=hi_by_row,
    )


def power_vs_dependence(
    power_table: pd.DataFrame,
    *,
    arm: str = "E_copula_ar1_durations",
    alpha: float = NOMINAL_ALPHA,
) -> PowerVsDependence:
    frame = _with_label(power_table)
    cell = frame[frame["arm"] == arm]
    if not len(cell):
        raise ValueError(f"no {arm} power cells in the table")
    induced: dict[int, list[float]] = {}
    power: dict[str, dict[int, list[float]]] = {}
    for row_label, group in cell.groupby("row"):
        per_n: dict[int, list[float]] = {}
        for n, sub in group.groupby("n_target"):
            ordered = sub.groupby("rho").mean(numeric_only=True).sort_index()
            per_n[int(n)] = [float(v) for v in ordered["rejection_rate"]]
            lag1 = [float(v) for v in ordered["mean_induced_lag1"]]
            # `induced` is shared by every row at this n, which is only sound because
            # `mean_induced_lag1` is a CELL-level scalar stamped onto each of the cell's
            # rows. Assert it rather than rely on it: if a future table made it row-level,
            # the last row written would win and every other row would be plotted against
            # the wrong x, silently.
            if int(n) in induced and not np.allclose(induced[int(n)], lag1, atol=1e-12):
                raise ValueError(
                    f"induced lag-1 differs between rows at n={n}; it is assumed "
                    "cell-level and shared as the x-axis of every curve"
                )
            induced[int(n)] = lag1
        power[row_label] = per_n
    return PowerVsDependence(
        induced_lag1_by_n=induced,
        power_by_row_and_n=power,
        alpha=alpha,
        rows_drawn=_representative_rows(power),
        selection_reason=(
            "permutation calibration, studentized where a variant exists - the asymptotic "
            "curves are within Monte Carlo error of them and the bench rejected both"
        ),
    )


def _representative_rows(power: dict[str, dict[int, list[float]]]) -> list[str]:
    """One row per check: the calibration a reader should actually use.

    Drawing all nine put C1, C2 and CvM on the page twice in visually identical panels - the
    asymptotic and permutation power curves differ by less than the Monte Carlo error,
    while their SIZE (the other figure) differs enormously. Showing both here implies a
    distinction this figure cannot resolve.
    """
    chosen: list[str] = []
    for check in sorted({_check_of(row) for row in power}):
        candidates = [row for row in power if _check_of(row) == check]
        preferred = [r for r in candidates if "permutation" in r] or candidates
        studentized = [r for r in preferred if "unstudentized" not in r]
        chosen.append(sorted(studentized or preferred)[0])
    return chosen


def _check_of(row_label: str) -> str:
    return row_label.split(" [", 1)[0]


def read_dependence(
    power_table: pd.DataFrame,
    *,
    arm: str = "C_ar1_reads_carved",
    clock: str = "in_spec",
) -> ReadDependence:
    frame = power_table[(power_table["arm"] == arm) & (power_table["clock"] == clock)]
    if not len(frame):
        raise ValueError(f"no {arm} power cells on the {clock} clock")
    grouped = frame.groupby("rho")
    rho = sorted(float(r) for r in frame["rho"].unique())
    induced = [float(grouped.get_group(r)["mean_induced_lag1"].mean()) for r in rho]
    se = [
        float(
            grouped.get_group(r)["mean_induced_lag1"].std(ddof=1)
            / np.sqrt(max(len(grouped.get_group(r)), 1))
        )
        for r in rho
    ]
    return ReadDependence(
        read_rho=rho, induced_lag1=induced, induced_lag1_se=se, arm=arm
    )


def _exponential_segment(tau: float, rng: np.random.Generator) -> Segment:
    """A time-truncated exponential segment: tau fixed, event count random.

    `tau` is a truncation TIME, not an event count. Gaps are unit-mean, so about `tau` of
    them land; the draw is generously oversized so the cumulative sum reaches `tau` with
    overwhelming probability.
    """
    tau = float(tau)
    gaps = rng.exponential(size=int(tau + 10 * np.sqrt(tau) + 50))
    kept = gaps[: int(np.searchsorted(np.cumsum(gaps), tau, side="left"))]
    return Segment(x=kept, tau=tau, n_censored_dropped=1)


def validation_curve(
    n_values: tuple[int, ...] = (20, 50, 355),
    n_replicates: int = 4000,
    seed: int = 20260811,
    n_points: int = 200,
) -> ValidationCurve:
    """Simulate the P-P curve for C2's eq (7) under both gamma paths.

    Deterministic in `seed`, which is a step kwarg, so the figure is reproducible and the
    seed appears on the provenance label. This is the one calibration artifact that is not
    a reshape of the bench tables: a P-P plot needs the statistics, and the tables carry
    only rejection rates.
    """
    rng = np.random.default_rng(seed)
    theoretical = list(np.linspace(1.0 / n_points, 1.0 - 1.0 / n_points, n_points))
    gamma_one: dict[int, list[float]] = {}
    gamma_hat: dict[int, list[float]] = {}
    for n in n_values:
        p_one = np.empty(n_replicates, dtype=float)
        p_hat = np.empty(n_replicates, dtype=float)
        for i in range(n_replicates):
            segment = _exponential_segment(n, rng)
            p_one[i] = 1.0 - float(
                c2.ad_limiting_cdf(c2._eq7(segment.x, segment.tau, 1.0))
            )
            p_hat[i] = float(c2.run([segment], calibration=CALIB_ASYMPTOTIC).p_value)
        # Empirical CDF of the p-values, read at the theoretical quantiles. Under a correct
        # calibration p is uniform, so this lands on the diagonal.
        gamma_one[int(n)] = [float(np.mean(p_one <= q)) for q in theoretical]
        gamma_hat[int(n)] = [float(np.mean(p_hat <= q)) for q in theoretical]
    return ValidationCurve(
        n_values=[int(n) for n in n_values],
        theoretical=theoretical,
        empirical_gamma_one=gamma_one,
        empirical_gamma_hat=gamma_hat,
        n_replicates=n_replicates,
        seed=seed,
    )


__all__ = [
    "CALIB_ASYMPTOTIC",
    "CALIB_PERMUTATION",
    "NOMINAL_ALPHA",
    "PowerVsDependence",
    "ReadDependence",
    "SizeVsN",
    "ValidationCurve",
    "power_vs_dependence",
    "read_dependence",
    "size_vs_n",
    "validation_curve",
]


# --------------------------------------------------------------------- bench acceptance
#
# "Is this check calibrated at this event count?" has ONE arithmetic definition, and it
# lives here rather than in `jobs/bench/report.py` so that the ledger (a pipeline step) and the
# report (a study) cannot drift apart. The dependency direction allows it: `jobs/bench/` may
# import `analyzers/`, never the reverse.
#
# What is shared is the per-cell z-test and its multiplicity correction. What is NOT shared
# is the aggregation: the report takes the worst cell over all n >= 35 to score a check
# overall, while the ledger asks about one n at a time. `jobs/bench/report.py` imports
# `null_se`, `bonferroni_z_crit` and the three constants from here, so the difference
# between them is the grouping they apply and nothing else.
#
# They can therefore disagree, by design and in both directions. Two live examples:
# C5-unstudentized is flagged by the ledger at n = 20 (z = -2.97 against 2.955) but reads
# calibrated in the report, which only scores n >= 35; C2-asymptotic reads REJECT overall
# in the report but is accepted by the ledger on the in-spec clock at n = 50, 75, 100 and
# 355, which is the more useful statement for a record with 355 windows.

# Censoring the real data exhibits (0.000-0.026 wherever n >= 20). Cells beyond it are a
# corner the data never reaches and must not condemn a check that works where it lives.
ENVELOPE_MAX_CENSORING = 0.03

# A rate is only read as a rate if it rests on at least this fraction of its replicates.
MIN_SUPPORT_FRACTION = 0.5

# Family-wise error rate for the calibration decision.
FAMILYWISE_ALPHA = 0.05


def null_se(alpha: float, n_used: float) -> float:
    """Standard error of a rejection rate UNDER THE NULL.

    Deliberately not the observed rate's SE, which collapses to zero on a cell that never
    rejected and would make a perfectly conservative cell look infinitely miscalibrated.
    """
    return float(np.sqrt(alpha * (1.0 - alpha) / max(float(n_used), 1.0)))


def bonferroni_z_crit(n_cells: int, familywise: float = FAMILYWISE_ALPHA) -> float:
    """Two-sided z threshold for the worst of `n_cells` comparisons.

    Without this the rule is a coin flip: the largest of 56-90 deviations is ~2.5-3 MC SE
    by chance alone, so a flat tolerance flags a correct check as often as a broken one.
    """
    from scipy import stats as _stats

    return float(_stats.norm.ppf(1.0 - familywise / (2.0 * max(int(n_cells), 1))))


def bench_acceptance_at_n(
    size_table: pd.DataFrame,
    *,
    alpha: float = NOMINAL_ALPHA,
    null_arms: tuple[str, ...] = (
        "A_iid_weibull",
        "B_iid_reads_carved",
        "C_ar1_reads_carved",
        "D_trp_power_law",
        "E_copula_ar1_durations",
    ),
) -> pd.DataFrame:
    """Per (check, calibration, variant, clock, n_target): is the check calibrated there?

    `bench_size_at_n` is the WORST adequately-supported null cell inside the censoring
    envelope at that event count - the conservative reading, so a pass means "calibrated
    even in the least favourable configuration tested". It is judged against a Bonferroni
    threshold for the number of cells that worst was taken over, which is what makes a
    flag mean miscalibrated rather than unlucky.
    """
    frame = _with_label(size_table)
    cells = frame[
        frame["arm"].isin(null_arms)
        & (frame["censoring_target"] <= ENVELOPE_MAX_CENSORING)
        & (frame["n_used"] >= MIN_SUPPORT_FRACTION * frame["n_replicates"])
    ].copy()
    if not len(cells):
        raise ValueError("no adequately supported null cells inside the envelope")
    cells["se"] = [null_se(alpha, u) for u in cells["n_used"]]
    cells["z"] = (cells["rejection_rate"] - alpha) / cells["se"]

    out = []
    for (row_label, clock, n_target), group in cells.groupby(
        ["row", "clock", "n_target"]
    ):
        worst = group.loc[group["z"].abs().idxmax()]
        z_crit = bonferroni_z_crit(len(group))
        out.append(
            {
                "row": row_label,
                "check": worst["check"],
                "calibration": worst["calibration"],
                "variant": worst["variant"]
                if isinstance(worst["variant"], str)
                else "",
                "clock": clock,
                "n_target": int(n_target),
                "bench_size_at_n": float(worst["rejection_rate"]),
                "bench_size_se": float(worst["se"]),
                "bench_size_z": float(worst["z"]),
                "bench_z_crit": float(z_crit),
                "bench_accepted": bool(abs(float(worst["z"])) <= z_crit),
                "bench_n_cells": int(len(group)),
                "bench_size_cell": (
                    f"{worst['arm']} q={bool(worst['quantised'])} "
                    f"c={worst['censoring_target']:g}"
                ),
            }
        )
    return pd.DataFrame(out)


def nearest_bracketing_n(available: list[int], n_events: int) -> int:
    """The grid point a real event count is judged at. Never interpolate.

    The bench measured six event counts; a record with 356 windows is judged at the
    nearest one and the table says which. Interpolating between grid points would invent
    a size that was never measured.
    """
    if not available:
        raise ValueError("no bench grid points to bracket against")
    return int(min(available, key=lambda g: (abs(g - n_events), g)))
