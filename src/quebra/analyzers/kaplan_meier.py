"""Kaplan-Meier survival of the windows on one side of a threshold.

This is a STEP: pure compute, no I/O, no matplotlib. It reads an event table
(`analyzers/event_table.py`) built from the window table `analyzers/windows.py`
produces, so the gap policy, the censoring and the window identity are the same facts
here as in every other consumer of that table. The product-limit arithmetic lives in
`product_limit` and the quantile rule in `survival_quantile`, both unit-agnostic, so
every estimator that needs them calls the same code.

Two decisions shape the estimate:

- **Right-censored windows are kept, not dropped.** A window that died at a read gap or
  at the end of the scan is not a completed lifetime, but it IS evidence that the window
  survived at least that long. Discarding it biases S(t) toward short lifetimes exactly
  where the record is thinnest; Kaplan-Meier uses it.

- **Windows whose birth was not observed are excluded.** Such a window
  (`birth_observed` False) was already on its side when observation began or
  resumed, so its age at first sight is unknown: its recorded duration is a residual
  lifetime, not a lifetime. Treating it as a lifetime understates survival; treating it
  as censored at that duration is also wrong. The honest handling is left truncation,
  which needs an entry-age this record does not carry, so these windows are excluded and
  the count is carried on the artifact - `n_unobserved_birth_dropped` - for the figure to
  state. FIGURE_STANDARD requires the exclusion be visible in the panel, not just here.

`KaplanMeierSet` holds one curve per threshold and side; it is a job's one Kaplan-Meier
node, and the within-calibration panel draws from it.

Curve durations are MINUTES (`_min`), because windows on this record run from seconds to a
few hours and hours would put every interesting feature below 0.1. Event tables, and every
other estimator, are in seconds.

Validity assumptions

- Assumption: the windows on the declared side are draws from one lifetime law, and a
  window's censoring is independent of its future.
  Diagnostic: the attached check outcome (`assumptions.A1_RENEWAL_DURATIONS`); censoring
  here happens only at read gaps and at the end of the scan.
  Consequence of violation: the curve estimates a mixture whose weights shift as windows
  leave the risk set, and the band need not hold its nominal level.
  Reference: R `survival` 3.8.6 `survfit(..., conf.type = "log-log")`, whose curve and
  band `tests/test_survival_r_reference.py` pins; no source located for the original
  papers.
- Assumption: the band is pointwise and asymptotic: the Greenwood variance, normal on the
  log-log scale.
  Diagnostic: `instrument_validation.measure_band_coverage`, coverage at one time on a
  simulated exponential with fixed censoring.
  Consequence of violation: read as simultaneous, or at small risk sets, the band claims
  more than it holds.
  Reference: as above, and `scipy.stats.ecdf(...).sf.confidence_interval(method="log-log")`,
  cross-checked in `tests/test_kaplan_meier.py`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.stats import norm

import icontract

from quebra.analyzers import event_table as _event_table
from quebra.analyzers.check_attachment import CheckAttachment, summarise
from quebra.analyzers.windows import (
    BIRTH_DOWN_CROSSING,
    BIRTH_GAP_RESUME,
    BIRTH_SCAN_START,
    BIRTH_UP_CROSSING,
    DEATH_DOWN_CROSSING,
    DEATH_GAP_START,
    DEATH_SCAN_END,
    DEATH_UP_CROSSING,
)

# The codes `windows.run` emits, on either side. Selection reads the `birth_observed` and
# `censored` booleans; `make_inputs_from_windows` refuses unrecognised codes and booleans
# that disagree with the codes, because either means the table was not produced by the
# carve. Both taxonomies are guarded: guarding one and not its twin is what AGENTS.md
# section 4 calls fixing one site of a class.
KNOWN_DEATH_TYPES = frozenset(
    {DEATH_DOWN_CROSSING, DEATH_UP_CROSSING, DEATH_GAP_START, DEATH_SCAN_END}
)
KNOWN_BIRTH_TYPES = frozenset(
    {BIRTH_UP_CROSSING, BIRTH_DOWN_CROSSING, BIRTH_SCAN_START, BIRTH_GAP_RESUME}
)
# A birth is observed exactly when it is a crossing; a window is censored exactly when it
# ended at a gap or at the end of the scan.
_OBSERVED_BIRTH_TYPES = frozenset({BIRTH_UP_CROSSING, BIRTH_DOWN_CROSSING})
_CENSORING_DEATH_TYPES = frozenset({DEATH_GAP_START, DEATH_SCAN_END})
# The event table's columns plus the two code columns checked here: one list, so a table
# cannot clear this function's check and fail the event table's.
_KM_WINDOW_COLUMNS = (*_event_table.WINDOW_COLUMNS_NEEDED, "birth_type", "death_type")

# The default confidence level of the band. The level is a parameter (`conf_level`), so
# the normal quantile is computed from it (`z_two_sided`) rather than stored beside it.
CONF_LEVEL = 0.95

# R's survival::quantile.survfit tolerance for "S sits exactly at p", sqrt(machine eps).
QUANTILE_TOL = float(np.sqrt(np.finfo(float).eps))


def z_two_sided(conf_level: float) -> float:
    """Two-sided normal quantile for a pointwise interval at `conf_level`."""
    if not 0.0 < conf_level < 1.0:
        raise ValueError(f"conf_level must lie in (0, 1); got {conf_level!r}")
    return float(norm.ppf(0.5 + conf_level / 2.0))


@dataclass(slots=True)
class KaplanMeierInputs:
    """Lifetimes and their censoring indicator, already reduced to observed births.

    death_observed[i] is True when window i died of an observed crossing back to the other
    side, False when it was right-censored (gap start or scan end).
    """

    duration_min: np.ndarray
    death_observed: np.ndarray
    label: str = ""
    dataset_id: str = ""
    threshold_label: str = ""
    # The side of the threshold the windows were on; empty for inputs built by hand.
    side: str = ""
    n_windows_carved: int = 0
    n_unobserved_birth_dropped: int = 0
    conf_level: float = CONF_LEVEL


@dataclass
class KaplanMeierCurve:
    """A complete Kaplan-Meier estimate: the step function, its band, and its support.

    `time_min` / `survival` are the left ends of the step function's segments, starting
    at (0, 1). Draw with `drawstyle="steps-post"`; the arrays are not resampled onto a
    grid, so a renderer never has to guess where a step fell.

    `band_lower` / `band_upper` are the log-log-transformed pointwise interval at
    `conf_level`. The transform is used rather than Greenwood-on-S directly because the
    plain interval leaves [0, 1] in both tails, which on a survival axis draws a band
    the estimator cannot mean. Entries are NaN where the interval is undefined
    (S = 1 before the first death, S = 0 after the last, or a risk set fully consumed).
    """

    time_min: np.ndarray
    survival: np.ndarray
    band_lower: np.ndarray
    band_upper: np.ndarray
    n_at_risk: np.ndarray

    # Censoring marks, at the S(t) the curve holds when each censored window leaves.
    censor_time_min: np.ndarray
    censor_survival: np.ndarray

    label: str = ""
    dataset_id: str = ""
    threshold_label: str = ""
    # The side of the threshold the curve estimates, so a curve opened alone says which.
    side: str = ""
    conf_level: float = CONF_LEVEL

    n_windows: int = 0
    n_deaths: int = 0
    n_censored: int = 0
    n_windows_carved: int = 0
    n_unobserved_birth_dropped: int = 0
    n_zero_duration: int = 0

    # `survival_quantile` on this minute curve. `survival.curve_summaries` applies the same
    # rule in seconds, so the two medians agree up to the rounding of the division by 60.
    median_survival_min: float | None = None
    # Largest duration in the risk set. S(t) is undefined beyond it; a figure that
    # extends the curve past this is drawing an extrapolation.
    max_observed_min: float = 0.0


@dataclass
class KaplanMeierComparison:
    """Several curves at ONE threshold, plus which two are furthest apart.

    `ranking` is every unordered pair, sorted by `separation` descending, so the pair the
    figure draws is a value in the artifact rather than an eyeball judgement made at draw
    time.

    ONE statistic decides the choice: `separation` = `log_time_separation`, the area between
    the two step curves integrated against d(log10 t) - the vertical gap the eye reads off
    the figure's own log-time axis, summed over the decades it spans. Units are
    survival-fraction x decades.

    The obvious alternative, the vertical supremum sup|S_a - S_b|, is NOT used and must
    not be reintroduced: it is capped at 1 and saturates whenever one curve reaches zero
    before the other starts falling, which is exactly the regime this comparison lives in.
    On the shipped five-dataset carve it returned 0.976 for three different pairs and
    could not rank them at all.

    `separation` is a DISTANCE, not a test: no p-value is attached and none is implied,
    because the curves being ranked were selected on the T2* mean of the same records.

    The check fields below are what makes `AGENTS.md` section 5 true structurally rather
    than editorially. A band and its check outcome in two sibling artifacts are "attached"
    only in the sense that both files land in the same directory; a reader who opens the
    band alone learns nothing about what was checked. They are plain strings and tuples of
    strings so they enter the run identity through `core/closure.py` and so this module
    keeps importing nothing from the check layer.

    `assumption_id` names the record in `analyzers/assumptions.py` the band rests on.
    `checks_asked` and `checks_unanswered` are the two halves section 5 requires be named:
    a grid of `not computed` cells satisfies "attached" and says nothing, so the outcome has
    to say what was asked for as well as what came back. Empty means NOT ASSESSED, which is
    a different claim from "assessed and nothing rejected".
    """

    curves: list[KaplanMeierCurve]
    threshold_label: str
    ranking: list[tuple[str, str, float]] = field(default_factory=list)
    pair: tuple[str, str] | None = None
    assumption_id: str = ""
    checks_asked: tuple[str, ...] = ()
    checks_unanswered: tuple[str, ...] = ()
    check_verdicts: tuple[tuple[str, str, str], ...] = ()

    def check_summary(self) -> str:
        """One line naming what was asked and what came back, for a caption to state.

        The renderer states this; it does not derive it (`AGENTS.md` section 3).
        """
        return summarise(self.checks_asked, self.checks_unanswered, self.check_verdicts)

    def curve(self, label: str) -> KaplanMeierCurve:
        for c in self.curves:
            if c.label == label:
                return c
        raise KeyError(
            f"no Kaplan-Meier curve labelled {label!r}. Have: "
            f"{[c.label for c in self.curves]}"
        )

    def pair_curves(self) -> tuple[KaplanMeierCurve, KaplanMeierCurve]:
        if self.pair is None:
            raise ValueError(
                "KaplanMeierComparison has no selected pair - it was built from fewer "
                "than two curves. Nothing to draw."
            )
        return self.curve(self.pair[0]), self.curve(self.pair[1])


def make_inputs_from_windows(
    windows: pd.DataFrame,
    *,
    threshold_label: str,
    label: str,
    dataset_id: str = "",
) -> KaplanMeierInputs:
    """Select one threshold's windows, drop unobserved births, carry the dropped count.

    The table must hold exactly one side (`windows_in_spec` or `windows_out_of_spec`), and
    that side is the one estimated. Raises on a mixed or empty side column, and on an
    unknown threshold label rather than returning an empty estimate: a silently empty
    survival curve is the wrong-but-plausible result this repo raises to avoid.
    """
    missing = [c for c in _KM_WINDOW_COLUMNS if c not in windows.columns]
    if missing:
        raise KeyError(
            f"Kaplan-Meier requires columns {missing} in the window table. "
            f"Columns: {list(windows.columns)}"
        )
    known = set(windows["threshold_label"].unique())
    if threshold_label not in known:
        raise KeyError(
            f"threshold {threshold_label!r} is not in the window table. "
            f"Carved thresholds: {sorted(known)}"
        )
    sides = sorted(set(windows["side"].unique()))
    if len(sides) != 1:
        raise ValueError(
            f"a window table for Kaplan-Meier must hold one side; this one holds {sides}. "
            "Pass windows_in_spec or windows_out_of_spec, not a concatenation."
        )
    # NOT a contract: this reads a frame column, and the useful message names the offending
    # values. A precondition would report the whole Series. Same invariant class, different
    # tool, and the boundary is the point rather than an inconsistency.
    unknown_births = set(windows["birth_type"].unique()) - KNOWN_BIRTH_TYPES
    if unknown_births:
        raise ValueError(
            f"window table carries unknown birth_type(s) {sorted(unknown_births)}. Known: "
            f"{sorted(KNOWN_BIRTH_TYPES)}. A table with unrecognised codes was not produced "
            f"by the carve, so its birth_observed column cannot be trusted either."
        )
    unknown = set(windows["death_type"].unique()) - KNOWN_DEATH_TYPES
    if unknown:
        raise ValueError(
            f"window table carries unknown death_type(s) {sorted(unknown)}. Known: "
            f"{sorted(KNOWN_DEATH_TYPES)}. A table with unrecognised codes was not produced "
            f"by the carve, so its censored column cannot be trusted either."
        )
    # Built first so its bool-dtype check on the two flags runs before they are compared.
    table = _event_table.from_windows(
        windows, threshold_label=threshold_label, side=sides[0], dataset_id=dataset_id
    )
    _require_flags_match_codes(windows)
    return make_inputs_from_event_table(table, label=label)


def _require_flags_match_codes(windows: pd.DataFrame) -> None:
    births = windows["birth_type"].isin(_OBSERVED_BIRTH_TYPES).to_numpy()
    deaths = windows["death_type"].isin(_CENSORING_DEATH_TYPES).to_numpy()
    for flag, expected, codes in (
        ("birth_observed", births, "birth_type"),
        ("censored", deaths, "death_type"),
    ):
        n_disagree = int(np.count_nonzero(windows[flag].to_numpy() != expected))
        if n_disagree:
            raise ValueError(
                f"{n_disagree} window(s) have a {flag} flag that disagrees with their "
                f"{codes} code. The carve sets both from one fact, so this table was not "
                "produced by it."
            )


def make_inputs_from_event_table(
    table: _event_table.EventTable, *, label: str, conf_level: float = CONF_LEVEL
) -> KaplanMeierInputs:
    """An event table (seconds) as Kaplan-Meier inputs (minutes)."""
    return KaplanMeierInputs(
        duration_min=table.age_s / 60.0,
        death_observed=table.event.copy(),
        label=label,
        dataset_id=table.dataset_id,
        threshold_label=table.threshold_label,
        side=table.side,
        n_windows_carved=table.n_windows + table.n_unobserved_birth_dropped,
        n_unobserved_birth_dropped=table.n_unobserved_birth_dropped,
        conf_level=conf_level,
    )


def product_limit(
    event_age: np.ndarray, n_events: np.ndarray, n_at_risk: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Kaplan-Meier steps from a risk table, in whatever unit the ages carry.

    Returns (time, survival, n_at_risk, greenwood), starting at (0, 1) and holding one
    entry per age with at least one death. `greenwood` is the running sum
    d / (n (n - d)); it is inf from a risk set fully consumed onward, where the variance
    is not defined.

    Reference: R `survival` 3.8.6 `survfit`: `surv`, and `std.err` squared, which is this
    sum at each death time; no source located for the original papers.
    """
    times, surv, at_risk, green = (
        [0.0],
        [1.0],
        [int(n_at_risk[0]) if len(n_at_risk) else 0],
        [0.0],
    )
    s, g = 1.0, 0.0
    for age, d, n in zip(event_age, n_events, n_at_risk):
        if d == 0:
            continue
        s *= 1.0 - d / n
        denom = n * (n - d)
        g += np.inf if denom == 0 else d / denom
        times.append(float(age))
        surv.append(float(s))
        at_risk.append(int(n))
        green.append(float(g))
    return (
        np.asarray(times, dtype=float),
        np.asarray(surv, dtype=float),
        np.asarray(at_risk, dtype=int),
        np.asarray(green, dtype=float),
    )


def survival_quantile(
    time: np.ndarray, survival: np.ndarray, p: float, max_time: float
) -> float | None:
    """The p-quantile of a right-continuous step curve, by R's survival convention.

    The first step time where S <= 1 - p. When S sits exactly at 1 - p there (within
    `QUANTILE_TOL`), the curve is flat at the quantile until its next step, and the
    quantile is that flat stretch's midpoint; a stretch with no later step ends at
    `max_time`, the largest observed age. None when S never reaches 1 - p. Without
    censoring this is the ordinary sample quantile: ages 10, 20, 30, 40 give a median of 25.

    `p` must lie in (0, 1). Outside it the rule returns plausible numbers for a quantile
    that does not exist (p = 1.5 reads as "never reached"), so it raises instead.

    Reference: R `survival` 3.8.6 `survival:::findq`, which `quantile.survfit` reaches
    through `doquant`.
    """
    if not 0.0 < p < 1.0:
        raise ValueError(f"quantile probability p must lie in (0, 1); got {p!r}")
    target = 1.0 - p
    reached = np.flatnonzero(survival <= target + QUANTILE_TOL)
    if not len(reached):
        return None
    j = int(reached[0])
    if abs(survival[j] - target) > QUANTILE_TOL:
        return float(time[j])
    end = float(time[j + 1]) if j + 1 < len(time) else float(max_time)
    return 0.5 * (float(time[j]) + end)


@icontract.require(
    lambda inputs: np.asarray(inputs.death_observed).dtype == np.bool_,
    description=(
        "death_observed must already be a boolean array. np.asarray(..., dtype=bool) "
        "coerces silently, so any nonzero float would be read as an OBSERVED DEATH and a "
        "censored window would enter the curve as a failure"
    ),
)
def run(inputs: KaplanMeierInputs) -> KaplanMeierCurve:
    """Kaplan-Meier product-limit estimate with a log-log pointwise band.

    The precondition is a hard computational invariant in the sense
    `spec/quebraplan.md` 4.2 means: violating it does not make the estimate uncertain, it
    makes it an estimate of something else. It is a contract rather than a bare raise
    because the condition is one expression and `icontract` reports the offending dtype
    without a hand-written message. `make_inputs_from_windows` below shows the other half
    of that boundary: a check needing frame inspection stays a bare raise.
    """
    t = np.asarray(inputs.duration_min, dtype=float)
    observed = np.asarray(inputs.death_observed, dtype=bool)
    if len(t) != len(observed):
        raise ValueError(
            f"duration_min ({len(t)}) and death_observed ({len(observed)}) must have "
            "the same length"
        )
    if len(t) == 0:
        raise ValueError(
            f"Kaplan-Meier has no windows to estimate from ({inputs.label!r}, "
            f"{inputs.threshold_label!r}). The carve produced "
            f"{inputs.n_windows_carved} window(s), all with unobserved births."
        )
    if not np.all(np.isfinite(t)):
        raise ValueError("duration_min contains non-finite values")
    if np.any(t < 0.0):
        raise ValueError("duration_min contains negative durations")

    # Ties: a death and a censoring at the same recorded time are ordered death-first, so
    # the censored window is still counted in that time's risk set (`risk_table`).
    n_total = len(t)
    time_min, survival, n_at_risk, greenwood = product_limit(
        *_risk_columns(t, observed)
    )
    # A risk set fully consumed makes the Greenwood term infinite; NaN propagates into the
    # band and the renderer draws no band there, the honest rendering of "the variance is
    # not defined here".
    lower, upper = loglog_band(survival, greenwood, z_two_sided(inputs.conf_level))

    censored_t = np.sort(t[~observed])
    censor_survival = _step_eval(time_min, survival, censored_t)

    median_min = survival_quantile(time_min, survival, 0.5, float(np.max(t)))

    n_deaths = int(np.count_nonzero(observed))
    curve = KaplanMeierCurve(
        time_min=time_min,
        survival=survival,
        band_lower=lower,
        band_upper=upper,
        n_at_risk=n_at_risk,
        censor_time_min=censored_t,
        censor_survival=censor_survival,
        label=inputs.label,
        dataset_id=inputs.dataset_id,
        threshold_label=inputs.threshold_label,
        side=inputs.side,
        conf_level=inputs.conf_level,
        n_windows=n_total,
        n_deaths=n_deaths,
        n_censored=n_total - n_deaths,
        n_windows_carved=inputs.n_windows_carved,
        n_unobserved_birth_dropped=inputs.n_unobserved_birth_dropped,
        n_zero_duration=int(np.count_nonzero(t == 0.0)),
        median_survival_min=median_min,
        max_observed_min=float(np.max(t)),
    )
    return curve


def _risk_columns(
    age: np.ndarray, event: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ages, n_events, _n_censored, n_at_risk = _event_table.risk_table(age, event)
    return ages, n_events, n_at_risk


def loglog_band(
    survival: np.ndarray, greenwood: np.ndarray, z: float
) -> tuple[np.ndarray, np.ndarray]:
    """Pointwise log-log interval: S^exp(+/- z * sqrt(V) / |log S|).

    `z` has no default: it is `z_two_sided(conf_level)` for the level the caller carries,
    so a band cannot silently be built at a level its result does not name.

    Undefined where log S is 0 (S = 1) or S is 0; NaN there, so a renderer cannot draw a
    band across a region the estimator says nothing about.

    Reference: R `survival` 3.8.6 `survfit(..., conf.type = "log-log")`; no source located
    for the original paper.
    """
    lower = np.full_like(survival, np.nan, dtype=float)
    upper = np.full_like(survival, np.nan, dtype=float)
    ok = (survival > 0.0) & (survival < 1.0) & np.isfinite(greenwood)
    if not np.any(ok):
        return lower, upper
    s_ok = survival[ok]
    se = np.sqrt(greenwood[ok]) / np.abs(np.log(s_ok))
    lower[ok] = np.clip(s_ok ** np.exp(z * se), 0.0, 1.0)
    upper[ok] = np.clip(s_ok ** np.exp(-z * se), 0.0, 1.0)
    return lower, upper


def _step_eval(
    step_times: np.ndarray, step_values: np.ndarray, at: np.ndarray
) -> np.ndarray:
    """Right-continuous step function S evaluated at `at`: the value AFTER the last step <= t."""
    if len(at) == 0:
        return np.asarray([], dtype=float)
    idx = np.searchsorted(step_times, np.asarray(at, dtype=float), side="right") - 1
    return step_values[np.clip(idx, 0, len(step_values) - 1)]


def _support_max_min(curve: KaplanMeierCurve) -> float:
    """Largest t at which S is defined.

    A curve whose last recorded window died of an observed crossing has S = 0 there, and
    S = 0 for every later t is a statement the estimator makes, not an extrapolation - so
    its support is unbounded. A curve still above zero at its longest window ended in a
    censored observation: S beyond that point is genuinely unknown, and its support stops.
    """
    return np.inf if curve.survival[-1] == 0.0 else curve.max_observed_min


def _require_estimated(curve: KaplanMeierCurve) -> None:
    if curve.n_windows == 0:
        raise ValueError(
            f"curve {curve.label!r} ({curve.threshold_label!r}, {curve.side!r}) holds no "
            "windows: its S = 1 is the typed empty curve, not an estimate, and has no "
            "distance to another curve."
        )


def log_time_separation(a: KaplanMeierCurve, b: KaplanMeierCurve) -> float:
    """Area between the two step curves against d(log10 t), in fraction x decades.

    Integrated over [t_lo, t_hi]: t_lo is the earlier of the two first step times (below
    it both curves are 1 and contribute nothing), t_hi the end of the shorter support,
    capped at the longer of the two records - past that both curves are flat and the
    integral would run forever.

    Both curves are step functions, so this is an EXACT rectangle sum on the union of
    their step times, not a quadrature approximation. Raises on a curve with no windows
    (an empty side): its support ends at 0, so it would score 0 against any curve and rank
    as identical to it.
    """
    _require_estimated(a)
    _require_estimated(b)
    t_hi = min(_support_max_min(a), _support_max_min(b))
    t_hi = min(t_hi, max(a.max_observed_min, b.max_observed_min))
    positive_steps = np.concatenate([a.time_min, b.time_min])
    positive_steps = positive_steps[positive_steps > 0.0]
    if len(positive_steps) == 0 or t_hi <= 0.0:
        return 0.0
    t_lo = float(np.min(positive_steps))
    if t_hi <= t_lo:
        return 0.0

    edges = np.unique(np.concatenate([positive_steps, [t_lo, t_hi]]))
    edges = edges[(edges >= t_lo) & (edges <= t_hi)]
    # The curves hold their value on [edges[i], edges[i+1]); evaluate at the left edge.
    left = edges[:-1]
    gap = np.abs(
        _step_eval(a.time_min, a.survival, left)
        - _step_eval(b.time_min, b.survival, left)
    )
    return float(np.sum(gap * np.diff(np.log10(edges))))


def compare(
    curves: list[KaplanMeierCurve], threshold_label: str
) -> KaplanMeierComparison:
    """Rank every pair by `log_time_separation` and record the widest-apart pair.

    All curves must share `threshold_label`: survival at different thresholds measures
    different events, and ranking distances across them compares nothing. Every curve
    must hold windows (`log_time_separation`).
    """
    for curve in curves:
        _require_estimated(curve)
    mismatched = [c.label for c in curves if c.threshold_label != threshold_label]
    if mismatched:
        raise ValueError(
            f"compare() was given curves carved at a different threshold than "
            f"{threshold_label!r}: {mismatched}. Survival at two thresholds is survival "
            "of two different events."
        )
    ranking = [
        (curves[i].label, curves[j].label, log_time_separation(curves[i], curves[j]))
        for i in range(len(curves))
        for j in range(i + 1, len(curves))
    ]
    ranking.sort(key=lambda row: row[2], reverse=True)
    pair = (ranking[0][0], ranking[0][1]) if ranking else None
    return KaplanMeierComparison(
        curves=list(curves),
        threshold_label=threshold_label,
        ranking=ranking,
        pair=pair,
    )


def curve_from_event_table(
    table: _event_table.EventTable, *, label: str = "", conf_level: float = CONF_LEVEL
) -> KaplanMeierCurve:
    """Kaplan-Meier on an event table, with a typed empty curve for an empty table.

    `run` refuses an empty estimate; an empty SIDE is data (a record always in spec has no
    out-of-spec windows), so it gets a curve that holds S = 1, says n = 0, and carries the
    dropped count.
    """
    if table.n_windows > 0:
        return run(
            make_inputs_from_event_table(table, label=label, conf_level=conf_level)
        )
    nan = np.array([np.nan])
    return KaplanMeierCurve(
        time_min=np.array([0.0]),
        survival=np.array([1.0]),
        band_lower=nan,
        band_upper=nan.copy(),
        n_at_risk=np.array([0]),
        censor_time_min=np.array([], dtype=float),
        censor_survival=np.array([], dtype=float),
        label=label,
        dataset_id=table.dataset_id,
        threshold_label=table.threshold_label,
        side=table.side,
        conf_level=conf_level,
        n_windows_carved=table.n_unobserved_birth_dropped,
        n_unobserved_birth_dropped=table.n_unobserved_birth_dropped,
    )


@dataclass
class KaplanMeierSet:
    """One Kaplan-Meier curve per (threshold label, side), from one set of event tables.

    The one Kaplan-Meier node of a job: the panel draws from it and every other consumer
    reads it, so a figure and a table built from the same windows cannot show two bands.
    `checks` is the check outcome its bands rest on; empty means NOT ASSESSED.
    """

    ladder: tuple[str, ...]
    conf_level: float
    curves: dict[tuple[str, str], KaplanMeierCurve]
    checks: CheckAttachment = field(default_factory=CheckAttachment)

    def curve(self, threshold_label: str, side: str) -> KaplanMeierCurve:
        key = (threshold_label, side)
        if key not in self.curves:
            raise KeyError(
                f"no Kaplan-Meier curve for {key}; ladder {list(self.ladder)}"
            )
        return self.curves[key]


def kaplan_meier_set(
    tables: _event_table.EventTables, *, conf_level: float = CONF_LEVEL
) -> KaplanMeierSet:
    """One curve per (threshold label, side) of `tables`, keyed as the tables are.

    `checks` is left empty, which means NOT ASSESSED.
    """
    return KaplanMeierSet(
        ladder=tables.ladder,
        conf_level=conf_level,
        curves={
            key: curve_from_event_table(table, label=key[0], conf_level=conf_level)
            for key, table in tables.tables.items()
        },
    )
