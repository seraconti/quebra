"""Turning the bench tables into the five typed artifacts the calibration figures draw.

These are steps: pure compute, no disk, no matplotlib. The bench tables arrive as
DataFrames that the job loaded as datasets, so nothing here imports `jobs/bench/` - the
dependency runs through provenance (their sha256 enters the run identity) rather than
around it.

Five questions, five artifacts, one per figure:

    SizeVsN            does a check hold its nominal level at the event counts we have?
    PowerVsDependence  how does its power grow with the serial dependence it is shown?
    PowerVsTrend       how does its power grow with a trend in the event rate?
    ValidationCurve    is the transcription right, and what does estimating gamma cost?
    ReadDependence     does read-level correlation reach the durations at all?

`ValidationCurve` is the one that is not a reshape of a table: it re-simulates, because a
P-P plot needs the statistics themselves and the tables only carry rejection rates. It is
a deterministic function of `(seed, n_values, n_replicates)`, and the seed travels as a
step kwarg so it lands on the provenance label.

Every artifact carries the rows of `battery.ROW_KEYS` only. C5's unstudentized variant is in
the bench tables as the measured comparison for what studentizing buys, and it reaches no
figure. A dependence axis is a quantity that is 0 with no dependence: the generator's own
rank correlation for Arm E, and the measured lag-1 with its exchangeable mean added back for
Arm C.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

import quebra.analyzers.checks.c2_anderson_darling as c2
from quebra.analyzers.checks.battery import BENCH_ROW_KEYS, ROW_KEYS
from quebra.analyzers.checks.result import CALIB_ASYMPTOTIC, CALIB_PERMUTATION, Segment
from quebra.core._artifact_guard import StaleArtifactGuard

# Nominal significance every size curve is read against.
NOMINAL_ALPHA = 0.05

# The bench arms the two power figures read: Arm E's latent correlation is the serial
# dependence axis, Arm D's power-law exponent the trend axis.
ARM_E = "E_copula_ar1_durations"
ARM_D = "D_trp_power_law"

# Why the power figures show one row per check, the permutation one.
REPRESENTATIVE_REASON = (
    "permutation calibration - the asymptotic curves are within Monte Carlo error of them "
    "and the bench rejected both"
)


def _label(check: str, calibration: str, variant: str) -> str:
    suffix = f"/{variant}" if isinstance(variant, str) and variant else ""
    return f"{check} [{calibration}{suffix}]"


def _refuse_failed_cells(selection: pd.DataFrame, what: str) -> None:
    """Raise if a cell in `selection` failed every replicate.

    `jobs/bench/runner.run_cell` writes such a cell as one row naming no check. Dropping it
    would leave the point it belongs to pooling fewer cells while still labelled a mean over
    all of them. Callers narrow `selection` to the cells they draw first, so a failed cell
    elsewhere in the table does not stop them.
    """
    failed = selection[selection["check"].isna()]
    if len(failed):
        first = failed.iloc[0]
        raise ValueError(
            f"{what}: {len(failed)} cell(s) it draws failed every replicate, the first "
            f"{first['arm']} n={first['n_target']} shape={first.get('shape')} "
            f"rho={first.get('rho')}: {first.get('fail_reason')}"
        )


def _require_full_pooling(counts: pd.Series, n_pooled: int, where: str) -> None:
    """Raise unless every point pools all `n_pooled` cells, so a mean labelled as over the
    pooled factors is one."""
    short = counts[counts != n_pooled]
    if len(short):
        raise ValueError(
            f"{where}: a point pools {int(short.iloc[0])} of {n_pooled} cells "
            f"(at {short.index[0]}); a cell is missing from the table"
        )


def _drawn_rows(selection: pd.DataFrame, what: str) -> pd.DataFrame:
    """The rows of `battery.ROW_KEYS` among the cells an artifact draws.

    A cell that failed every replicate raises (`_refuse_failed_cells`). C5's unstudentized
    row is dropped: it is in `BENCH_ROW_KEYS` only. Any other row in neither list raises, so
    a table carrying a check this code does not know is not drawn by omission.
    """
    _refuse_failed_cells(selection, what)
    frame = selection
    keys = list(
        zip(
            frame["check"],
            frame["calibration"],
            frame["variant"].fillna(""),
            strict=True,
        )
    )
    unknown = sorted(set(keys) - set(BENCH_ROW_KEYS), key=str)
    if unknown:
        raise ValueError(f"rows outside battery.BENCH_ROW_KEYS: {unknown}")
    return frame[[key in ROW_KEYS for key in keys]]


def neighbour_rank_correlation(rho: float) -> float:
    """Spearman's correlation of two neighbouring durations in Arm E: `(6/pi) arcsin(rho/2)`.

    Arm E's durations are a monotone map of a stationary Gaussian AR(1) whose lag-1
    correlation is `rho`, so their rank correlation is the Gaussian copula's. It is exactly
    0 at `rho = 0` and does not depend on the event count.
    """
    return float(6.0 / np.pi * np.arcsin(rho / 2.0))


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
    """Power against the rank correlation of neighbouring durations the generator induces.

    The x-axis is `neighbour_rank_correlation(rho)`, one value per `rho` shared by every n.
    It is 0 with no dependence and exact, so the curves for different n line up by
    construction. The measured `mean_induced_lag1` is not used: it is C5's uncentred `r_1`,
    about `-1/N` with no dependence and pulled further toward 0 at small N under dependence,
    so curves drawn against it shift with n.
    """

    rho: list[float] = field(default_factory=list)
    neighbour_rank_corr: list[float] = field(default_factory=list)
    power_by_row_and_n: dict[str, dict[int, list[float]]] = field(default_factory=dict)
    # Each power value is a MEAN over the cells sharing (row, n, rho) - in practice the two
    # Weibull shapes. `pooled_over` names what varies and the lo/hi dicts carry the spread,
    # as `SizeVsN` does, so a reader is not shown a pooled mean as "the cell".
    pooled_over: list[str] = field(default_factory=list)
    power_lo_by_row_and_n: dict[str, dict[int, list[float]]] = field(
        default_factory=dict
    )
    power_hi_by_row_and_n: dict[str, dict[int, list[float]]] = field(
        default_factory=dict
    )
    alpha: float = NOMINAL_ALPHA
    # One representative row per check, chosen HERE rather than at draw time. The artifact
    # keeps every row (completeness); this says which the figure shows and why.
    rows_drawn: list[str] = field(default_factory=list)
    selection_reason: str = ""


@dataclass
class PowerVsTrend(StaleArtifactGuard):
    """Rejection rate against the trend exponent `b` of Arm D, a Weibull trend-renewal
    process whose event intensity is `b t^(b-1)`.

    `b = 1` is no trend; `b < 1` is a falling event rate and `b > 1` a rising one. The
    `b = 1` point is a size, so the curves start from the nominal level there. As in
    `PowerVsDependence`, each value is a mean over the cells sharing (row, n, b), named in
    `pooled_over`, with its range in the lo/hi dicts.
    """

    b: list[float] = field(default_factory=list)
    power_by_row_and_n: dict[str, dict[int, list[float]]] = field(default_factory=dict)
    pooled_over: list[str] = field(default_factory=list)
    power_lo_by_row_and_n: dict[str, dict[int, list[float]]] = field(
        default_factory=dict
    )
    power_hi_by_row_and_n: dict[str, dict[int, list[float]]] = field(
        default_factory=dict
    )
    alpha: float = NOMINAL_ALPHA
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

    `centred_lag1` is each cell's `mean_induced_lag1` plus `1/mean_n_events`, then averaged
    over the event counts at one `rho`. The column is C5's uncentred `r_1`, whose mean with
    no dependence is `-1/N` on one segment, so without the correction each point would mix
    offsets from about -1/20 at n = 20 to about -1/1000 at n = 1000. `1/mean_n_events`
    stands in for the mean of `1/N` over replicates.
    """

    read_rho: list[float] = field(default_factory=list)
    centred_lag1: list[float] = field(default_factory=list)
    # Standard error of the MEAN ACROSS the pooled cells (one per event count at this rho),
    # not the sampling error of any single lag-1 estimate. It answers "how much do the
    # cells disagree", which is the right question for a curve drawn through their mean -
    # but it is not a confidence interval on the statistic, and naming it so here keeps the
    # figure's error bars from being read as one.
    centred_lag1_se: list[float] = field(default_factory=list)
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
    selected = size_table[
        (size_table["arm"] == arm)
        & (~size_table["quantised"].astype(bool))
        & (size_table["censoring_target"] == 0.0)
    ]
    if not len(selected):
        raise ValueError(f"no {arm} continuous, uncensored size cells in the table")
    cell = _with_label(_drawn_rows(selected, "size_vs_n"))
    n_values = sorted(int(n) for n in cell["n_target"].unique())
    # Which factors still vary after the filter. Anything listed here is averaged over,
    # and saying so is the whole difference between "the primary null cell" and a mean.
    pooled_over = [
        column
        for column in ("shape", "quantised", "censoring_target")
        if cell[column].nunique() > 1
    ]
    n_pooled = len(cell[pooled_over].drop_duplicates()) if pooled_over else 1
    rate_by_row: dict[str, list[float]] = {}
    se_by_row: dict[str, list[float]] = {}
    lo_by_row: dict[str, list[float]] = {}
    hi_by_row: dict[str, list[float]] = {}
    for row_label, group in cell.groupby("row"):
        by_n = group.groupby("n_target")
        _require_full_pooling(
            by_n.size().reindex(n_values, fill_value=0),
            n_pooled,
            f"size_vs_n {row_label}",
        )
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
    arm: str = ARM_E,
    alpha: float = NOMINAL_ALPHA,
) -> PowerVsDependence:
    """Power per row and event count against `neighbour_rank_correlation(rho)`.

    Validity assumptions

    - Assumption: the arm's durations are a monotone map of a stationary Gaussian AR(1)
      latent with lag-1 correlation `rho`, so `(6/pi) arcsin(rho/2)` is their neighbours'
      rank correlation.
      Diagnostic: `tests/test_calibration_artifacts.py` compares the centred measured `r_1`
      at the largest event count with it; any other arm raises here.
      Consequence of violation: every curve is drawn at the wrong x, by the same amount at
      every n.
      Reference: `jobs/bench/arms.py`, `arm_e`; no source located for the identity itself.
    """
    if arm != ARM_E:
        raise ValueError(
            f"the x-axis is Arm E's Gaussian-copula rank correlation; {arm!r} is not Arm E"
        )
    selected = power_table[power_table["arm"] == arm]
    if not len(selected):
        raise ValueError(f"no {arm} power cells in the table")
    cell = _with_label(_drawn_rows(selected, "power_vs_dependence"))
    grid = _power_grid(cell, "rho", "power_vs_dependence")
    return PowerVsDependence(
        rho=grid.values,
        neighbour_rank_corr=[neighbour_rank_correlation(r) for r in grid.values],
        power_by_row_and_n=grid.mean,
        pooled_over=grid.pooled_over,
        power_lo_by_row_and_n=grid.lo,
        power_hi_by_row_and_n=grid.hi,
        alpha=alpha,
        rows_drawn=_representative_rows(grid.mean),
        selection_reason=REPRESENTATIVE_REASON,
    )


@dataclass
class _PowerGrid:
    values: list[float]
    mean: dict[str, dict[int, list[float]]]
    lo: dict[str, dict[int, list[float]]]
    hi: dict[str, dict[int, list[float]]]
    pooled_over: list[str]


def _power_grid(cell: pd.DataFrame, parameter: str, what: str) -> _PowerGrid:
    """Rejection rate per row and event count at each value of `parameter`: the mean over
    the cells that share them, with its range, and the factors that vary among those cells.

    Every curve must cover the same `parameter` values, since the x-axis is shared, and
    every point must pool the same number of cells (`_require_full_pooling`).
    """
    values = sorted(float(v) for v in cell[parameter].unique())
    pooled_over = [
        column
        for column in ("shape", "clock", "quantised", "censoring_target")
        if cell[column].nunique() > 1
    ]
    n_pooled = len(cell[pooled_over].drop_duplicates()) if pooled_over else 1
    grid = _PowerGrid(values=values, mean={}, lo={}, hi={}, pooled_over=pooled_over)
    for row_label, group in cell.groupby("row"):
        mean_n: dict[int, list[float]] = {}
        lo_n: dict[int, list[float]] = {}
        hi_n: dict[int, list[float]] = {}
        for n, sub in group.groupby("n_target"):
            by_value = sub.groupby(parameter)["rejection_rate"]
            ordered = by_value.mean().sort_index()
            if [float(v) for v in ordered.index] != values:
                raise ValueError(
                    f"{row_label} at n={n} covers {parameter} {list(ordered.index)}, "
                    f"not {values}"
                )
            _require_full_pooling(
                by_value.size(), n_pooled, f"{what} {row_label} at n={n}"
            )
            mean_n[int(n)] = [float(v) for v in ordered]
            lo_n[int(n)] = [float(v) for v in by_value.min().sort_index()]
            hi_n[int(n)] = [float(v) for v in by_value.max().sort_index()]
        grid.mean[row_label] = mean_n
        grid.lo[row_label] = lo_n
        grid.hi[row_label] = hi_n
    return grid


def power_vs_trend(
    power_table: pd.DataFrame,
    size_table: pd.DataFrame,
    *,
    alpha: float = NOMINAL_ALPHA,
) -> PowerVsTrend:
    """Rejection rate per row and event count against Arm D's trend exponent `b`.

    The `b != 1` cells are Arm D's power cells. `b = 1` is no trend, Arm D's uncensored size
    cells, so each curve passes through the check's size there. Only uncensored cells are
    drawn, which is what every power cell is.
    """
    trend = power_table[power_table["arm"] == ARM_D]
    null = size_table[
        (size_table["arm"] == ARM_D)
        & (size_table["b"] == 1.0)
        & (size_table["censoring_target"] == 0.0)
    ]
    if not len(trend):
        raise ValueError(f"no {ARM_D} power cells in the table")
    if not len(null):
        raise ValueError(f"no {ARM_D} b = 1 uncensored size cells in the table")
    selected = pd.concat([trend, null], ignore_index=True)
    cell = _with_label(_drawn_rows(selected, "power_vs_trend"))
    grid = _power_grid(cell, "b", "power_vs_trend")
    return PowerVsTrend(
        b=grid.values,
        power_by_row_and_n=grid.mean,
        pooled_over=grid.pooled_over,
        power_lo_by_row_and_n=grid.lo,
        power_hi_by_row_and_n=grid.hi,
        alpha=alpha,
        rows_drawn=_representative_rows(grid.mean),
        selection_reason=REPRESENTATIVE_REASON,
    )


def _representative_rows(power: dict[str, dict[int, list[float]]]) -> list[str]:
    """One row per check: the calibration a reader should actually use.

    Drawing every row put C1, C2 and CvM on the page twice in visually identical panels -
    the asymptotic and permutation power curves differ by less than the Monte Carlo error,
    while their SIZE (the other figure) differs enormously. Showing both here implies a
    distinction this figure cannot resolve.
    """
    chosen: list[str] = []
    for check in sorted({_check_of(row) for row in power}):
        candidates = [row for row in power if _check_of(row) == check]
        preferred = [r for r in candidates if "permutation" in r] or candidates
        chosen.append(sorted(preferred)[0])
    return chosen


def _check_of(row_label: str) -> str:
    return row_label.split(" [", 1)[0]


def read_dependence(
    power_table: pd.DataFrame,
    *,
    arm: str = "C_ar1_reads_carved",
    clock: str = "in_spec",
) -> ReadDependence:
    """`centred_lag1` per read-level `rho`, one cell per event count; see `ReadDependence`.

    Raises on a cell with more than one segment, where `-1/N` is not the exchangeable mean;
    on a cell whose rows disagree on `mean_induced_lag1`, which is a cell-level value; and
    on a cell whose every replicate failed, which would otherwise enter a mean as NaN.

    Validity assumptions

    - Assumption: every cell is one segment, so the exchangeable mean of `r_1` is `-1/N`.
      Diagnostic: the raise above, on `mean_n_segments`.
      Consequence of violation: the centring misses the level effect of the layout.
      Reference: Dufour and Roy (1985), Eq. (2.2), p. 260.
    - Assumption: `1/mean_n_events` stands in for the mean of `1/N` over replicates.
      Diagnostic: none here; the table carries no per-replicate event count.
      Consequence of violation: a small offset is left in, largest at the smallest n.
      Reference: no source located.
    """
    frame = power_table[(power_table["arm"] == arm) & (power_table["clock"] == clock)]
    if not len(frame):
        raise ValueError(f"no {arm} power cells on the {clock} clock")
    _refuse_failed_cells(frame, "read_dependence")
    by_cell = frame.groupby(["n_target", "rho"])
    if (by_cell["mean_induced_lag1"].nunique() > 1).any():
        raise ValueError("mean_induced_lag1 differs between rows of one cell")
    cells = by_cell.first().reset_index()
    several = cells[cells["mean_n_segments"] != 1.0]
    if len(several):
        raise ValueError(
            f"{len(several)} {arm} cells have more than one segment, where -1/N is not "
            "the exchangeable mean of r_1; centring them by 1/N would be wrong"
        )
    cells = cells.assign(
        centred=cells["mean_induced_lag1"] + 1.0 / cells["mean_n_events"]
    )
    grouped = cells.groupby("rho")
    rho = sorted(float(r) for r in cells["rho"].unique())
    centred = [float(grouped.get_group(r)["centred"].mean()) for r in rho]
    se = [
        float(
            grouped.get_group(r)["centred"].std(ddof=1)
            / np.sqrt(len(grouped.get_group(r)))
        )
        for r in rho
    ]
    return ReadDependence(
        read_rho=rho, centred_lag1=centred, centred_lag1_se=se, arm=arm
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
# They can therefore disagree, by design and in both directions. A live example:
# C2-asymptotic reads REJECT overall in the report but is accepted by the ledger on the
# in-spec clock at every n from 50 up, which is the more useful statement for a record of a
# few hundred windows. The other direction comes from the multiplicity: the ledger's per-n
# threshold corrects for fewer cells than the report's, so a cell can exceed it while
# staying inside the report's band.

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

    Without this the rule is a coin flip: the largest of 86-144 deviations is ~2.5-3 MC SE
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

    The bench measured a fixed set of event counts (`jobs/bench/grid.py:N_GRID`); a record
    is judged at the nearest one and the table says which. Interpolating between grid
    points would invent a size that was never measured.
    """
    if not available:
        raise ValueError("no bench grid points to bracket against")
    return int(min(available, key=lambda g: (abs(g - n_events), g)))
