"""The instrument report: what licenses each statistical routine for use in the pipeline.

The bench measures how a routine BEHAVES; it cannot catch a
mis-transcribed equation that behaves plausibly. This module assembles the four tiers of
evidence that can, into one typed artifact a figure draws:

  Tier 1  the routine computes what its source's equation says          (code review)
  Tier 2  it reproduces the numbers the source PRINTS, from the source's own data
  Tier 3  it holds its nominal level under a null it is entitled to - a p-value's
          size, or, for an estimator with no p-value, a band's coverage
  Tier 4  it agrees with an independent implementation, usually the authors' own R package

No tier subsumes another. A routine can pass Tier 2 and fail Tier 3 (right arithmetic,
wrong null), or pass Tier 3 and fail Tier 4 (well-calibrated, but calibrated for a
different statistic than the one claimed).

**Pure compute.** Every input arrives as a DataFrame the job loaded; nothing here touches
disk and nothing imports matplotlib. The published tables, the R reference values and the
tie experiment are all declared as Datasets by `jobs/active/instrument_validation.py`, so
a figure's dependence on them runs through provenance rather than around it.

**What this report is NOT.** It is not a promotion decision and it does not score power.
`jobs/bench/results/promotion_report.md` does that, from the calibration bench, where every
check in the battery has cells and C3 has none. A routine can be perfectly validated here and
still be unusable on a short window.

**Cited, not copied.** A number in the report is either computed here, from the published
record, the R reference values or this run's own Monte Carlo, or it is not quoted: a claim
resting on the test suite names the test, and one resting on the bench names the table and
the cells. A copied number cannot follow its source.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

# Verdict vocabulary. Deliberately three-valued in the same spirit as the ledger's: an
# absent tier is NOT a failing tier, and the figure must not let the two look alike.
TIER_PASS = "pass"
TIER_FAIL = "fail"
TIER_PARTIAL = "partial"
TIER_ABSENT = "absent"

# Agreement at or below this counts as "the same number, computed twice". Chosen far below
# any plausible transcription error and far above cross-language float noise.
AGREEMENT_TOL = 1e-9


@dataclass(frozen=True)
class TierRow:
    """One instrument's evidence at one tier."""

    instrument: str
    tier: int
    verdict: str
    detail: str


@dataclass(frozen=True)
class PublishedComparison:
    """Tier 2: one quantity, ours against the number the paper prints."""

    quantity: str
    ours: float
    published: float
    agrees: bool
    note: str

    @property
    def relative_difference(self) -> float:
        if self.published == 0.0:
            return float("nan")
        return abs(self.ours - self.published) / abs(self.published)


@dataclass(frozen=True)
class CrossImplementation:
    """Tier 4: one statistic on one case, ours against R."""

    case: str
    statistic: str
    ours: float
    reference: float
    abs_difference: float
    reference_is_random: bool
    reference_spread: float


@dataclass(frozen=True)
class InstrumentValidationData:
    """The complete typed artifact. The renderer reads fields and draws; it computes nothing."""

    published: list[PublishedComparison]
    cross_implementation: list[CrossImplementation]
    tiers: list[TierRow]
    tie_experiment: pd.DataFrame
    divergence_tie_fraction: dict[int, float]
    # The same crossing in the form a reader can act on: how few distinct values the
    # response must take before the closed form stops agreeing with the permutation.
    divergence_levels: dict[int, float]
    mc_floor: float
    dropped: dict[str, int] = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    @property
    def instruments(self) -> list[str]:
        seen: list[str] = []
        for row in self.tiers:
            if row.instrument not in seen:
                seen.append(row.instrument)
        return seen

    def tier_verdict(self, instrument: str, tier: int) -> TierRow | None:
        for row in self.tiers:
            if row.instrument == instrument and row.tier == tier:
                return row
        return None


# A truncation time, not an event count. Gaps are unit-mean, so about this many events land,
# and how many is random.
ASYMPTOTIC_SIZE_TAU = 20.0
ASYMPTOTIC_SIZE_REPLICATES = 1200
ASYMPTOTIC_SIZE_SEED = 777


def measure_asymptotic_size(
    module,
    *,
    tau: float = ASYMPTOTIC_SIZE_TAU,
    seed: int = ASYMPTOTIC_SIZE_SEED,
    alpha: float = 0.05,
    replicates: int = ASYMPTOTIC_SIZE_REPLICATES,
) -> float:
    """Rejection rate of an asymptotic check on iid gaps - its actual size at this tau.

    LIVES HERE, not in the test module, and that is the point. A previous version of this
    measurement sat in `tests/test_tier3_calibration.py` while three of its docstrings and
    the report's tier-3 verdicts claimed the report imported it. It did not - `grep` matched
    that one file - so "asymptotic size measured" was a verdict resting on no number in the
    artifact. The pipeline may not import `tests/`, so the only way to make the claim true
    is for the measurement to live on this side and the TEST to import it. It now does.

    The null simulated here has to be the one the theory assumes: the truncation time is
    chosen without reference to the events, which `checks/result` states as a requirement.
    Nothing downstream can check it for you - `checks/battery` gates the asymptotic rows on
    `tau == T_N`, which an event-derived `T_N(1 + 1/n)` clears.

    Deriving tau from the draw as `g.sum() + g.mean()` gives exactly that, so `T_N/tau` is
    pinned at `n/(n+1)` where the null has it Beta(n, 1), and the count is fixed where it
    should be random. The leftover window then carries no variance, which for C2 is eq (7)'s
    tail term. `_exponential_segment` fixes tau and lets the count fall where it falls.

    Pure compute: seeded generator in, float out, no disk and no matplotlib.
    """
    import numpy as _np

    from quebra.analyzers.calibration_summary import _exponential_segment
    from quebra.analyzers.checks.result import CALIB_ASYMPTOTIC

    rng = _np.random.default_rng(seed)
    rejected = 0
    for _ in range(replicates):
        segment = _exponential_segment(tau, rng)
        p = module.run([segment], calibration=CALIB_ASYMPTOTIC).p_value
        # A `None` p-value counts as a non-rejection, which is only defensible because the
        # event count is now random: at N < 2 the checks raise rather than declining, so a
        # None here would be a new failure mode rather than a quiet zero.
        rejected += int(p is not None and p <= alpha)
    return rejected / replicates


BAND_COVERAGE_SEED = 20260907
BAND_COVERAGE_REPLICATES = 5000
BAND_COVERAGE_N = 60
BAND_COVERAGE_CENSOR_AT = 1.5
BAND_COVERAGE_EVAL_AT = 0.5


# Deterministic in its arguments: same seed and settings, same float. Cached because
# `build_instrument_validation` is constructed once per test in the instrument suite and
# 5000 Kaplan-Meier fits per construction is 1 s each time, for a number that cannot change.
@functools.lru_cache(maxsize=8)
def measure_band_coverage(
    *,
    seed: int = BAND_COVERAGE_SEED,
    replicates: int = BAND_COVERAGE_REPLICATES,
    n: int = BAND_COVERAGE_N,
    censor_at: float = BAND_COVERAGE_CENSOR_AT,
    eval_at: float = BAND_COVERAGE_EVAL_AT,
) -> float:
    """Coverage of Kaplan-Meier's 95% log-log band against a known exponential survival.

    The tier-3 analogue for an estimator that produces a BAND rather than a p-value: the
    question "does it hold its nominal level under a null it is entitled to" becomes "does
    the 95% band contain the truth 95% of the time". Same question, different instrument.

    Lives here rather than in the test module for the reason `measure_asymptotic_size`
    gives above: a tier-3 verdict must quote a number the artifact produced, so the
    measurement is on this side and the test imports it.

    The truth is `S(t) = exp(-t)` and censoring is administrative at a FIXED time, so it is
    non-informative by construction and the estimator is entitled to the data. Coverage is
    the mean over replicates at ONE evaluation time; it is not a max or a min over a grid,
    which would need a multiplicity correction before any verdict could be read off it.

    Pure compute: seeded generator in, float out, no disk and no matplotlib.
    """
    import contextlib
    import io

    import numpy as _np

    from quebra.analyzers import kaplan_meier as _km

    rng = _np.random.default_rng(seed)
    truth = float(_np.exp(-eval_at))
    at = _np.array([eval_at])
    covered = 0
    # `kaplan_meier.run` logs one line per call by design; 5000 of them would bury the job's
    # own output. Suppressed here rather than made conditional in the estimator.
    with contextlib.redirect_stdout(io.StringIO()):
        for _ in range(replicates):
            draw = rng.exponential(size=n)
            curve = _km.run(
                _km.KaplanMeierInputs(
                    duration_min=_np.minimum(draw, censor_at),
                    death_observed=draw <= censor_at,
                )
            )
            lo = _km._step_eval(curve.time_min, curve.band_lower, at)[0]
            hi = _km._step_eval(curve.time_min, curve.band_upper, at)[0]
            covered += int(_np.isfinite(lo) and _np.isfinite(hi) and lo <= truth <= hi)
    return covered / replicates


def asymptotic_size_se(
    size: float, replicates: int = ASYMPTOTIC_SIZE_REPLICATES
) -> float:
    """Monte Carlo standard error of a measured size.

    Reported beside the size because at these replicate counts it is around 0.007, so the
    third and fourth decimals of the rate are noise and a bare `0.0708` reads as far more
    precise than the measurement is. The bench tables already carry an SE column.
    """
    return float((size * (1.0 - size) / replicates) ** 0.5)


def measure_all_asymptotic_sizes(
    *,
    tau: float = ASYMPTOTIC_SIZE_TAU,
    seed: int = ASYMPTOTIC_SIZE_SEED,
    replicates: int = ASYMPTOTIC_SIZE_REPLICATES,
) -> dict[str, float]:
    """The three asymptotic instruments' measured size, keyed by the report's own names."""
    import quebra.analyzers.checks.c1_lewis_robinson as _c1
    import quebra.analyzers.checks.c2_anderson_darling as _c2
    import quebra.analyzers.checks.cvm_cramer_von_mises as _cvm

    return {
        name: measure_asymptotic_size(mod, tau=tau, seed=seed, replicates=replicates)
        for name, mod in (
            ("C1 Lewis-Robinson", _c1),
            ("C2 Anderson-Darling", _c2),
            ("CvM", _cvm),
        )
    }


def _published_lookup(published_df: pd.DataFrame) -> dict[str, float]:
    return {
        str(r.quantity): float(r.value) for r in published_df.itertuples(index=False)
    }


def load_haul_dump_record(
    gaps_df: pd.DataFrame, published_df: pd.DataFrame
) -> tuple[int, float, int]:
    """(complete gaps, censoring time in hours, censored gaps) of the published record.

    The record is time censored at `tau_h`: observation stops there, so whatever time is left
    after the last failure is one censored gap, which the complete-gap estimators drop. The
    counts are derived from the record rather than written down, so a caption quoting them
    cannot outlive it. `tau_h` is compared with `last_event_time` under `TAU_MARGIN`, as every
    guard on a segment is, so summation order cannot decide whether a gap is censored.
    """
    from quebra.analyzers.checks.result import TAU_MARGIN, last_event_time

    x = gaps_df["gap_h"].to_numpy(dtype=float)
    tau_h = _published_lookup(published_df)["tau_h"]
    t_n = last_event_time(x)
    if tau_h < t_n * (1.0 - TAU_MARGIN):
        raise ValueError(
            f"the complete gaps end at {t_n:.10g} h, past the censoring time {tau_h:.10g} h"
        )
    return len(x), tau_h, int(tau_h > t_n * (1.0 + TAU_MARGIN))


def build_published_comparisons(
    gaps_df: pd.DataFrame, published_df: pd.DataFrame
) -> list[PublishedComparison]:
    """Tier 2 against Kvaloy and Lindqvist's load-haul-dump record, Section 6.1.

    The DIVISOR divergence is carried, not hidden. Their Table 2 prints sigma_hat, gamma_hat
    and LR on the sample `1/(N-1)` divisor; `_multiprocess.gamma_hat` uses the population
    `1/N` form deliberately, to match eq (10)'s divisor so the two estimators differ ONLY by
    the residual term. All three of their numbers land simultaneously under one rescaling by
    `sqrt(N/(N-1))`, and it is that SIMULTANEITY that identifies the cause as the divisor -
    a single number rescaled by 1.014 at one N cannot be told from any nearby constant.

    Consequence, stated rather than buried: our gamma_hat is smaller by `sqrt((N-1)/N)`, so
    our C1 statistic is LARGER by `sqrt(N/(N-1))` - 1.4% at N = 36, about 12% on a five-event
    segment. That direction is anti-conservative.
    """
    import quebra.analyzers.checks.c1_lewis_robinson as c1
    from quebra.analyzers.checks._multiprocess import (
        GAMMA_COMPLETE,
        GAMMA_TRUNCATED,
        gamma_hat,
    )
    from quebra.analyzers.checks.result import Segment

    x = gaps_df["gap_h"].to_numpy(dtype=float)
    pub = _published_lookup(published_df)
    tau = pub["tau_h"]
    n = len(x)
    segment = Segment(
        x=x, tau=tau, n_censored_dropped=load_haul_dump_record(gaps_df, published_df)[2]
    )

    gamma_complete = gamma_hat(x, tau, GAMMA_COMPLETE)
    gamma_truncated = gamma_hat(x, tau, GAMMA_TRUNCATED)
    lr = c1.statistic([segment], gamma_estimator=GAMMA_COMPLETE)
    mu = float(np.mean(x))
    bridge = float(np.sqrt(n / (n - 1)))

    rows = [
        PublishedComparison("mu_hat", mu, pub["mu_hat"], True, "mean gap (h)"),
        PublishedComparison(
            "gamma_tilde (eq 10)",
            gamma_truncated,
            pub["gamma_tilde"],
            True,
            "uses the censored time",
        ),
        PublishedComparison(
            "sigma_tilde (eq 10)",
            gamma_truncated * tau / n,
            pub["sigma_tilde"],
            True,
            "composed: gamma_tilde * tau / N",
        ),
        PublishedComparison(
            "Laplace",
            lr * gamma_complete,
            pub["laplace"],
            True,
            "composed: LR * gamma_hat",
        ),
        PublishedComparison(
            "sigma_hat",
            gamma_complete * mu,
            pub["sigma_hat"],
            False,
            f"1/N vs their 1/(N-1); x{bridge:.4f} reconciles",
        ),
        PublishedComparison(
            "gamma_hat",
            gamma_complete,
            pub["gamma_hat"],
            False,
            f"1/N vs their 1/(N-1); x{bridge:.4f} reconciles",
        ),
        PublishedComparison(
            "LR",
            lr,
            pub["lr_gamma_hat"],
            False,
            f"1/N vs their 1/(N-1); /{bridge:.4f} reconciles",
        ),
    ]
    for row in rows:
        if row.agrees and row.relative_difference > 1e-3:
            raise ValueError(
                f"{row.quantity} was marked as agreeing with the published value but "
                f"differs by {row.relative_difference:.2%}. The artifact must not assert "
                "an agreement the numbers do not show."
            )
    return rows


def build_cross_implementation(
    inputs_df: pd.DataFrame, values_df: pd.DataFrame
) -> list[CrossImplementation]:
    """Tier 4 against the R packages, several of them by the methods' own authors.

    Handles the case the experiment turned up: `XICOR::xicor` is a RANDOM VARIABLE when x
    has ties, because eq (8) breaks those ties uniformly at random. There, "the reference
    value" is a distribution and the honest comparison is membership, so the row carries
    `reference_is_random` and the observed spread instead of pretending to an equality.
    """
    from quebra.analyzers.shape_stats import chatterjee_xi, dcor

    def value(case: str, quantity: str) -> float | None:
        rows = values_df[
            (values_df["case"] == case) & (values_df["quantity"] == quantity)
        ]
        return float(rows["value"].iloc[0]) if len(rows) == 1 else None

    def inputs(case: str) -> tuple[np.ndarray, np.ndarray]:
        rows = inputs_df[inputs_df["case"] == case].sort_values("i")
        return (
            rows["x"].to_numpy(dtype=float),
            rows["y"].to_numpy(dtype=float),
        )

    out: list[CrossImplementation] = []
    for case, statistic, fn in [
        ("xi_tie_free", "xi", chatterjee_xi),
        ("xi_tied_y_only", "xi", chatterjee_xi),
        ("dcor_nonmonotone", "xi", chatterjee_xi),
        ("xi_tie_free", "dcor", dcor),
        ("xi_tied", "dcor", dcor),
        ("dcor_nonmonotone", "dcor", dcor),
    ]:
        reference = value(case, "xicor" if statistic == "xi" else "dcor")
        if reference is None:
            continue
        x, y = inputs(case)
        ours = float(fn(x, y))
        out.append(
            CrossImplementation(
                case=case,
                statistic=statistic,
                ours=ours,
                reference=reference,
                abs_difference=abs(ours - reference),
                reference_is_random=False,
                reference_spread=0.0,
            )
        )

    # The tied case, where the reference is a distribution rather than a number.
    mean, sd = value("xi_tied", "xicor_mean"), value("xi_tied", "xicor_sd")
    if mean is not None and sd is not None:
        x, y = inputs("xi_tied")
        ours = float(chatterjee_xi(x, y))
        out.append(
            CrossImplementation(
                case="xi_tied",
                statistic="xi",
                ours=ours,
                reference=mean,
                abs_difference=abs(ours - mean),
                reference_is_random=True,
                reference_spread=sd,
            )
        )
    return out


def _divergence_levels(tie_df: pd.DataFrame, threshold: float) -> dict[int, float]:
    """Same crossing as `_divergence_points`, reported in DISTINCT RESPONSE LEVELS.

    Levels is the reportable form. Tie fraction saturates - quantising to 20 levels already
    ties 83 to 99 per cent of a sample - so "diverges at tie fraction 1.000" is true of
    several different cells and tells a reader nothing they can act on. "Diverges once the
    response has only 2 distinct values" is the same fact in a form they can check against
    their own metric.
    """
    out: dict[int, float] = {}
    for n, sub in tie_df.groupby("n"):
        # Coarsest first: y_levels 0 means continuous, which is the finest, so it sorts last.
        ordered = sub.assign(_order=sub["y_levels"].replace(0, 10**9)).sort_values(
            "_order", ascending=False
        )
        crossed = ordered[ordered["p_signed_diff_mean"].abs() > threshold]
        # `iloc[0]`, the FINEST crossing: the reportable number is where divergence STARTS
        # as the response coarsens, not the coarsest cell that also happens to cross. With
        # both 3 and 2 levels crossing, the answer is 3.
        out[int(n)] = (
            float(crossed.iloc[0]["y_levels"]) if len(crossed) else float("nan")
        )
    return out


def _divergence_points(tie_df: pd.DataFrame, threshold: float) -> dict[int, float]:
    """The deliverable: tie fraction at which the two p-values SYSTEMATICALLY diverge.

    Read off the signed mean, not the absolute one. At B = 999 the permutation p-value
    carries Monte Carlo error near 0.016, which is most of a 0.02 threshold, so the
    absolute difference reports a divergence that is mostly its own noise. Noise averages
    out of the signed mean; systematic disagreement does not.

    A tie fraction of NaN means the grid never crossed the threshold, which is a result and
    is drawn as one - not as a missing value.
    """
    out: dict[int, float] = {}
    for n, sub in tie_df.groupby("n"):
        ordered = sub.sort_values("tie_fraction_y")
        crossed = ordered[ordered["p_signed_diff_mean"].abs() > threshold]
        out[int(n)] = (
            float(crossed.iloc[0]["tie_fraction_y"]) if len(crossed) else float("nan")
        )
    return out


SIZE_TABLE_PATH = "jobs/bench/results/size_table.csv"
PROMOTION_REPORT_PATH = "jobs/bench/results/promotion_report.md"


def _bench_size_note(check: str) -> str:
    """What the bench's size cells show for one check, cited by table and cell key."""
    return (
        "On Weibull gaps its asymptotic cells at n = 20 have point estimates above nominal "
        f"for both shapes the bench runs ({SIZE_TABLE_PATH}: check={check}, calibration=asymptotic, "
        "arm=A_iid_weibull, clock=in_spec, quantised=False, censoring_target=0.0, "
        "n_target=20). Its permutation cells hold size at n >= 35 and censoring <= 0.03 "
        f"({SIZE_TABLE_PATH}: check={check}, calibration=permutation), as "
        f"{PROMOTION_REPORT_PATH} judges them: Bonferroni within each row, none across rows"
    )


def build_instrument_validation(
    gaps_df: pd.DataFrame,
    published_df: pd.DataFrame,
    r_inputs_df: pd.DataFrame,
    r_values_df: pd.DataFrame,
    tie_df: pd.DataFrame,
    *,
    divergence_threshold: float = 0.02,
    n_perm_in_tie_study: int = 999,
    asymptotic_size_seed: int = ASYMPTOTIC_SIZE_SEED,
    asymptotic_size_tau: float = ASYMPTOTIC_SIZE_TAU,
    asymptotic_size_replicates: int = ASYMPTOTIC_SIZE_REPLICATES,
    band_coverage_seed: int = BAND_COVERAGE_SEED,
    band_coverage_replicates: int = BAND_COVERAGE_REPLICATES,
) -> InstrumentValidationData:
    """Assemble the whole artifact. This is the step the job calls."""
    published = build_published_comparisons(gaps_df, published_df)
    band_coverage = measure_band_coverage(
        seed=band_coverage_seed, replicates=band_coverage_replicates
    )
    # Measured HERE, into the artifact, so the tier-3 rows below quote a number this run
    # produced rather than one a docstring remembers.
    sizes = measure_all_asymptotic_sizes(
        replicates=asymptotic_size_replicates,
        tau=asymptotic_size_tau,
        seed=asymptotic_size_seed,
    )
    cross = build_cross_implementation(r_inputs_df, r_values_df)

    def measured(instrument: str) -> str:
        size = sizes[instrument]
        se = asymptotic_size_se(size, asymptotic_size_replicates)
        return (
            f"size measured at tau={asymptotic_size_tau:g}: {size:.4f} +/- {se:.4f} on "
            "exponential gaps (nominal 0.05). "
        )

    n_exact = sum(1 for row in published if row.agrees)
    xi_exact = [r for r in cross if r.statistic == "xi" and not r.reference_is_random]
    dcor_exact = [
        r for r in cross if r.statistic == "dcor" and not r.reference_is_random
    ]
    r_ref = {
        (str(r.case), str(r.quantity)): float(r.value)
        for r in r_values_df.itertuples(index=False)
    }
    c3_seed = int(r_ref[("durations_iid", "serial_indep_sim_seed")])
    c3_n_sim = int(r_ref[("durations_iid", "serial_indep_N")])
    c3_lag_max = int(r_ref[("durations_iid", "serial_indep_lag_max")])
    exactness = (
        "tests/test_tier3_calibration.py::"
        "test_the_permutation_p_value_is_exact_under_exchangeability"
    )

    # Tier table. Every verdict here is a statement about evidence that EXISTS in this
    # repo, not a judgement of the underlying method.
    tiers = [
        TierRow(
            "C1 Lewis-Robinson",
            2,
            TIER_PARTIAL,
            f"{n_exact} of {len(published)} published quantities exact; "
            f"{len(published) - n_exact} differ by the 1/N vs 1/(N-1) divisor",
        ),
        TierRow(
            "C1 Lewis-Robinson",
            3,
            TIER_PARTIAL,
            measured("C1 Lewis-Robinson") + _bench_size_note("c1_lewis_robinson"),
        ),
        TierRow(
            "C1 Lewis-Robinson", 4, TIER_ABSENT, "no independent implementation in hand"
        ),
        TierRow(
            "C2 Anderson-Darling",
            2,
            TIER_PASS,
            "the published limiting critical values reproduced "
            "(tests/test_checks_statistics.py::"
            "test_ad_limiting_cdf_reproduces_published_critical_values)",
        ),
        TierRow(
            "C2 Anderson-Darling",
            3,
            TIER_PARTIAL,
            measured("C2 Anderson-Darling") + _bench_size_note("c2_anderson_darling"),
        ),
        TierRow(
            "C2 Anderson-Darling",
            4,
            TIER_PARTIAL,
            "the gamma=1 path equals a textbook Anderson-Darling written in our own test "
            "file (tests/test_checks_statistics.py::"
            "test_eq7_is_the_classical_anderson_darling) - same language and same author, "
            "so it is tier-1 evidence wearing a tier-4 label. No independent "
            "implementation of eq (7) is in hand",
        ),
        TierRow(
            "C3 serial copula", 2, TIER_ABSENT, "no published worked example in hand"
        ),
        TierRow(
            "C3 serial copula",
            3,
            TIER_ABSENT,
            "no bench cell; run only on the iid reference input "
            "(tests/test_r_cross_implementation.py::"
            "test_the_bridge_reproduces_the_reference_serial_indep_values, which needs R)",
        ),
        TierRow(
            "C3 serial copula",
            4,
            TIER_PARTIAL,
            "BRIDGE FIDELITY, not cross-implementation agreement. "
            "tests/test_r_cross_implementation.py::"
            "test_the_bridge_reproduces_the_reference_serial_indep_values reproduces "
            "jobs/reference/r_reference_values.csv on BOTH the statistic and the p-value at "
            f"seed {c3_seed}, N {c3_n_sim}, lag.max {c3_lag_max}, which pins argument "
            "marshalling and the CSV round trip. It is not tier-4 evidence: the fixture and "
            "the bridge both call copula::serialIndepTest, so this is R against itself. An "
            "independent implementation of the statistic is still absent",
        ),
        TierRow(
            "C5 rank autocorr",
            3,
            TIER_PASS,
            f"permutation exactness asserted ({exactness})",
        ),
        TierRow(
            "C5 rank autocorr",
            4,
            TIER_ABSENT,
            "no independent implementation of C5's statistic is in hand. R's lag-1 rank "
            "autocorrelation in the reference fixture is a Pearson correlation of the lagged "
            "pairs; C5 divides every lag by the full centred sum of squares, so the two "
            "differ by construction and one cannot vouch for the other",
        ),
        TierRow(
            "C6 exchangeability",
            3,
            TIER_PASS,
            f"permutation exactness asserted ({exactness})",
        ),
        TierRow(
            "C6 exchangeability",
            4,
            TIER_ABSENT,
            "no independent implementation in hand",
        ),
        TierRow(
            "CvM",
            2,
            TIER_PASS,
            "the published limiting critical values reproduced (tests/test_checks_cvm.py::"
            "test_the_limiting_cdf_reproduces_the_published_critical_values)",
        ),
        TierRow(
            "CvM",
            3,
            TIER_PASS,
            measured("CvM") + _bench_size_note("cvm_cramer_von_mises"),
        ),
        TierRow(
            "CvM",
            4,
            TIER_PASS,
            "the gamma=1 path equals scipy.stats.cramervonmises (tests/test_checks_cvm.py::"
            "test_at_gamma_one_it_is_the_classical_cramer_von_mises)",
        ),
        TierRow(
            "Chatterjee xi",
            3,
            TIER_PARTIAL,
            "closed-form null holds at nominal on tie-free data; under ties it is not "
            "entitled (the xi tie experiment, below)",
        ),
        TierRow(
            "Chatterjee xi",
            4,
            TIER_PASS,
            f"matches XICOR on the {len(xi_exact)} exact cases of the tier-4 detail (max "
            f"abs difference {max(r.abs_difference for r in xi_exact):.1e}) and lies inside "
            "its tie-break distribution on tied data (tests/test_r_cross_implementation.py::"
            "test_xi_lies_inside_the_XICOR_tie_break_distribution_on_TIED_data); matches "
            "scipy.stats.chatterjeexi too (tests/test_r_cross_implementation.py::"
            "test_xi_matches_scipys_chatterjeexi). Two independent implementations, one of "
            "which runs in the suite without R",
        ),
        TierRow(
            "distance correlation",
            4,
            TIER_PASS,
            f"matches the energy package on the {len(dcor_exact)} cases of the tier-4 "
            f"detail (max abs difference {max(r.abs_difference for r in dcor_exact):.1e})",
        ),
        TierRow(
            "Kaplan-Meier",
            2,
            TIER_ABSENT,
            "no published worked example with printed intermediates is in hand, as there "
            "is for C1 and C2. The derivation-from-definition evidence that would be tier 1 "
            "exists (a censored case hand-derived independently, product-limit, risk sets "
            "and Greenwood, in tests/test_kaplan_meier.py) but no renderer draws tier 1, so "
            "it is not claimed as a row here",
        ),
        TierRow(
            "Kaplan-Meier",
            3,
            TIER_PASS,
            f"band coverage measured at t={BAND_COVERAGE_EVAL_AT:g}: "
            f"{band_coverage:.4f} +/- {asymptotic_size_se(band_coverage, band_coverage_replicates):.4f} "
            f"against S(t)=exp(-t), n={BAND_COVERAGE_N}, administrative censoring at "
            f"{BAND_COVERAGE_CENSOR_AT:g} (nominal 0.95). Coverage at ONE evaluation time, "
            "not a max over a grid, so no multiplicity correction is owed. No direction is "
            "claimed",
        ),
        TierRow(
            "Kaplan-Meier",
            4,
            TIER_PARTIAL,
            "agrees with scipy.stats.ecdf on CensoredData, point estimate and log-log band "
            "(tests/test_kaplan_meier.py::"
            "test_the_product_limit_estimate_agrees_with_scipys_ecdf, "
            "tests/test_kaplan_meier.py::"
            "test_the_log_log_band_agrees_with_scipys_log_log_interval). The point estimate "
            "is compared as step functions on the union grid; the band at Kaplan-Meier's own "
            "jump points, where both are finite. PARTIAL and not pass because that "
            "comparison runs in the suite and produces no CrossImplementation row here",
        ),
    ]

    divergence = _divergence_points(tie_df, divergence_threshold)
    divergence_levels = _divergence_levels(tie_df, divergence_threshold)
    n_gaps, tau_h, n_censored = load_haul_dump_record(gaps_df, published_df)
    return InstrumentValidationData(
        published=published,
        cross_implementation=cross,
        tiers=tiers,
        tie_experiment=tie_df.copy(),
        divergence_tie_fraction=divergence,
        divergence_levels=divergence_levels,
        mc_floor=0.5 / float(np.sqrt(n_perm_in_tie_study)),
        dropped={"load_haul_dump_censored_gaps": n_censored},
        meta={
            "load_haul_dump_complete_gaps": n_gaps,
            "load_haul_dump_tau_h": tau_h,
            "divergence_threshold": divergence_threshold,
            "n_published_quantities": len(published),
            "n_cross_checks": len(cross),
            "agreement_tol": AGREEMENT_TOL,
            # The seed and n DECIDE the three tier-3 numbers above, so they belong on the
            # artifact. Left as bare defaults they reached neither meta nor the provenance
            # label - the same defect class as an un-declared seed.
            "asymptotic_size_seed": asymptotic_size_seed,
            "asymptotic_size_tau": asymptotic_size_tau,
            "asymptotic_size_replicates": asymptotic_size_replicates,
            # The seed and replicate count DECIDE the tier-3 band row, so they
            # belong on the artifact for the same reason the asymptotic ones do.
            "band_coverage_seed": band_coverage_seed,
            "band_coverage_replicates": band_coverage_replicates,
            "band_coverage_n": BAND_COVERAGE_N,
            "band_coverage_eval_at": BAND_COVERAGE_EVAL_AT,
            "band_coverage_censor_at": BAND_COVERAGE_CENSOR_AT,
        },
    )


def _span(values: pd.Series, fmt: str) -> str:
    """`lo to hi` of a column, formatted; one value when they coincide."""
    lo, hi = format(float(values.min()), fmt), format(float(values.max()), fmt)
    return lo if lo == hi else f"{lo} to {hi}"


def render_tier_table_markdown(data: InstrumentValidationData) -> str:
    """The instrument x tier table as markdown. Pure: returns a string, writes nothing.

    Markdown rather than a rendered panel, deliberately. CLAUDE.md bans status tables in
    agent-facing docs because they go stale silently; the answer here is not to draw the
    table into a PNG - which goes stale just as silently and cannot be diffed - but to
    GENERATE it, so a stale table is a visible diff on the next run.
    """
    lines = [
        "# Instrument report - what licenses each routine",
        "",
        "GENERATED by `jobs/bench/instrument_report.py` from the same typed artifact the four",
        "`instrument_validation` figures draw. Do not hand-edit: the next run overwrites it.",
        "It lives beside `promotion_report.md` because that is this repo's one tracked home",
        "for generated markdown, NOT because it is a bench product - the bench does not",
        "score any of the tiers below.",
        "",
        "Tiers are independent kinds of evidence, and none subsumes another:",
        "",
        "| tier | question it answers |",
        "|---|---|",
        "| 1 | does the code compute what its source's equation says (review) |",
        "| 2 | does it reproduce the numbers the source PRINTS, on the source's data |",
        "| 3 | does it hold its nominal level under a null it is entitled to - a p-value's size, or a band's coverage |",
        "| 4 | does it agree with an independent implementation of the same statistic |",
        "",
        "`absent` means no evidence was gathered. It is NOT a failing grade, and it is not",
        "a pass either - read it as a gap in this repo's evidence, nothing more.",
        "",
        "**None of this is power.** An instrument can be fully validated here and still be",
        "unable to detect anything on a window this project actually carves. Power is in",
        "`promotion_report.md`, which scores every check in the battery; C3 has no bench",
        "cell.",
        "",
        "## The table",
        "",
        "| instrument | tier 2 | tier 3 | tier 4 |",
        "|---|---|---|---|",
    ]
    for instrument in data.instruments:
        cells = []
        for tier in (2, 3, 4):
            row = data.tier_verdict(instrument, tier)
            cells.append(f"{row.verdict} - {row.detail}" if row else TIER_ABSENT)
        lines.append(f"| {instrument} | " + " | ".join(cells) + " |")

    lines += [
        "",
        "## Tier 2 detail - published values",
        "",
        "Kvaloy and Lindqvist Section 6.1, load-haul-dump record, "
        f"{data.meta['load_haul_dump_complete_gaps']} complete gaps,",
        (
            f"time censored at {data.meta['load_haul_dump_tau_h']:.10g} h."
            if data.dropped["load_haul_dump_censored_gaps"]
            else f"observed to {data.meta['load_haul_dump_tau_h']:.10g} h, ending on a failure."
        ),
        "",
        "| quantity | ours | published | agrees | note |",
        "|---|---|---|---|---|",
    ]
    for row in data.published:
        lines.append(
            f"| {row.quantity} | {row.ours:.4f} | {row.published:.4f} | "
            f"{'yes' if row.agrees else 'NO'} | {row.note} |"
        )

    lines += [
        "",
        "## Tier 4 detail - cross-implementation",
        "",
        "| case | statistic | ours | reference | difference | note |",
        "|---|---|---|---|---|---|",
    ]
    for row in data.cross_implementation:
        note = (
            f"reference is RANDOM (sd {row.reference_spread:.4f}); membership, not equality"
            if row.reference_is_random
            else "exact reference value"
        )
        lines.append(
            f"| {row.case} | {row.statistic} | {row.ours:.10f} | {row.reference:.10f} | "
            f"{row.abs_difference:.2e} | {note} |"
        )

    threshold = data.meta["divergence_threshold"]
    tie = data.tie_experiment
    # Each sentence selects by the property it states: the control is tie-FREE, and the
    # saturation claim is about quantising to a number of levels.
    free = tie.loc[tie["tie_fraction_y"] == 0.0, "p_abs_diff_mean"]
    control = (
        "On the tie-FREE control cells, where the true difference is ~0, `|dp|` still "
        f"reads {_span(free, '.4f')} across n."
        if len(free)
        else "The table holds no tie-free control cell."
    )
    quantised = tie[tie["y_levels"] > 0]
    if len(quantised):
        finest = int(quantised["y_levels"].max())
        tied = quantised.loc[quantised["y_levels"] == finest, "tie_fraction_y"]
        saturation = (
            f"saturates: quantising to {finest} levels already ties {_span(tied, '.1%')} of "
            "a sample across n, so"
        )
    else:
        saturation = "saturates as the response coarsens, so"
    lines += [
        "",
        "## The xi tie experiment",
        "",
        "How coarse the response must get before the closed-form p-value stops agreeing",
        f"with the permutation - divergence beyond {threshold} in the SIGNED mean.",
        "",
        "Read off the signed mean, not `|dp|`: the permutation p-value carries Monte Carlo",
        f"error near {data.mc_floor:.4f}, which is most of the threshold, so `|dp|` reports",
        "its own noise as a divergence.",
        control,
        "",
        "Reported in DISTINCT RESPONSE VALUES rather than tie fraction, because tie fraction",
        saturation,
        "'diverges at tie fraction 1.000' is true of several different cells and cannot be",
        "checked against a real metric.",
        "",
        "| n | diverges once the response has | (tie fraction there) |",
        "|---|---|---|",
    ]
    for n, levels in sorted(data.divergence_levels.items()):
        tie = data.divergence_tie_fraction.get(n, float("nan"))
        if not np.isfinite(levels):
            lines.append(f"| {n} | never diverges | - |")
        else:
            tie_shown = "-" if not np.isfinite(tie) else f"{tie:.3f}"
            lines.append(
                f"| {n} | {int(levels)} distinct values or fewer | {tie_shown} |"
            )
    lines += [
        "",
        "The cost of ties shows up as LEVEL, not as p-value disagreement. See the centre",
        "panel of the tie-experiment figure: the permutation holds its nominal rate across",
        "every cell, while the closed form inflates as the response coarsens.",
    ]
    lines.append("")
    return "\n".join(lines)
