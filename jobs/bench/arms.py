"""The four data-generating arms.

Each arm returns `list[Segment]` - whatever it does internally, the checks see the same
object, so a size number from Arm A and one from Arm B are comparable.

    A  iid Weibull durations, direct                  null for C1, C2, C5, C6
    B  iid Gaussian reads, thresholded and CARVED     null, through the real carve
    C  AR(1) READ-level dependence, carved            null (see below) - not a power arm
    D  Weibull TRP, power-law trend, direct           power for C1, C2 (b != 1)
    E  copula AR(1) DURATION-level dependence         power for C2, C5, C6 (rho > 0)

A, D and E produce durations directly; B and C produce READS and put them through
`bench/carve.py`, so they measure the checks against the pipeline's own window definition
rather than against an idealisation of it. The pair (A, B) at the same n is the bench's own
consistency check: if carving an iid read series changes a size, the difference is the
carve's doing and not the check's.

**Arm C turned out to be a null arm, and Arm E was added because of it.** The plan had
Arm C as the power arm for C5/C6, on the assumption that correlated reads yield correlated
durations. They do not: the induced duration-level lag-1 never exceeds +0.002 over the
whole rho grid (`promotion_report.md`, "Two findings the grid was not designed to
produce"), because successive level crossings of a stationary Gaussian process are very
nearly a renewal process. Arm C is kept exactly as briefed since that is a real and
reportable finding, but Arm E is what actually varies duration dependence, and without it
criterion 2 could not be scored for C5 or C6 at all.

**tau is chosen BEFORE the events, and the record is truncated at it.** Gaps are drawn,
cumulated, and every event with `T_i < tau` is kept; the leftover `tau - T_N` is a
partial gap. This is what time censoring means, and two earlier constructions failed on
it, both caught by measurement:

- `tau = sum(x)` over exactly `n` gaps makes the `i = N-1` summand of Kvaloy & Lindqvist's
  eq (7), `ln((tau - T_{N-1})/(tau - T_N))`, `+inf`,
  so the briefed "n durations right-censored at tau" would have failed every replicate.
- Drawing `n + 1` gaps and setting `tau = sum` of all of them fixes that but breaks eq
  (10): the residual is then a COMPLETE gap rather than a partial one, so `mu = tau/N`
  overstates the mean, `sigma^2 = (1/N)[sum(x^2) + residual^2] - mu^2` goes NEGATIVE, and
  `gamma_hat` collapses to zero. At n = 50 in 17 segments that killed every replicate.

Truncating at a pre-chosen `tau` makes the residual partial by construction and `tau > T_N`
strict with probability one. It also means the event count is RANDOM around its target -
which is how Kvaloy & Lindqvist index their own Figure 1 (Technometrics 62(1) 2020, p. 107,
caption: "expected number of events"), so
it makes the comparison against their numbers more direct rather than less. `n` is a target
throughout; `mean_n_events` in the tables is what was realised.

The exact exponential oracle survives: with `tau` fixed and exponential gaps, conditional
on `N` the ratios `T_i/tau` are exactly uniform order statistics. That property is checked
by `tests/test_checks_statistics.py::test_eq7_at_gamma_one_matches_the_limiting_ad_null`,
which re-implements the construction rather than importing this module, so it is evidence
about the property and not a guard on this generator: no test exercises `arm_a`.

**Censoring is realised as SEGMENTS.** A censored unit in this pipeline is a window that
died at a read gap or at the end of the scan, and each of those terminates a segment. So a
censoring fraction of `c` means `m ~ c*n/(1-c)` segments, each contributing one incomplete
trailing gap - which is exactly the multi-process structure eqs (13)-(16) exist for. The
elevated `c = 0.25` arm is therefore the one that actually exercises what C1 and C2 buy
over a test that ignores censoring.

The realised fraction is reported, never assumed, and it does not always reach the target:
`m` is capped so each segment expects at least `MIN_EXPECTED_EVENTS_PER_SEGMENT` events.
Without that cap, `c = 0.25` at `n = 20` asks for 7 segments of expected size 2.9, and
since the count within a segment is random, about 86% of replicates would contain a segment
with fewer than the two events every check needs. A capped, honestly-labelled 0.17 is worth
more than a nominal 0.25 that fails most of the time.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy import stats
from scipy.signal import lfilter
from scipy.special import gamma as gamma_fn

from jobs.bench.carve import carve_windows
from quebra.analyzers.checks._multiprocess import segments_from_windows
from quebra.analyzers.checks.result import CLOCK_IN_SPEC, Segment

ARM_A = "A_iid_weibull"
ARM_B = "B_iid_reads_carved"
ARM_C = "C_ar1_reads_carved"
ARM_D = "D_trp_power_law"
ARM_E = "E_copula_ar1_durations"
ARMS = (ARM_A, ARM_B, ARM_C, ARM_D, ARM_E)

# Reads per unit of mean duration in the quantised variant of Arms A and D. Chosen to sit
# in the regime the real data occupies at the 4-7 us thresholds, where 79% of windows are
# a single read long and ties dominate the rank statistics; the continuous variant stands
# in for the 3 us threshold, which spans k = 1..72 with 43 distinct values.
QUANTUM_MEAN_READS = 3.0

# Arm B/C read grid. Uniform spacing, so `spacing()` sees a clean median and `gap_flags`
# finds nothing - these arms have no injected gaps, which is why their censoring is
# emergent (one terminal window) rather than a factor that can be dialled.
READ_DT_S = 1.0

# In-spec probability for the thresholded arms. p = 0.5 maximises the number of runs per
# read (N*p*(1-p)), so it is the cheapest place to buy a given event count, and it makes
# in-spec and out-of-spec runs identically distributed.
IN_SPEC_P = 0.5

# Discarded leading draws in Arm C, so the retained AR(1) series is stationary. At the
# largest rho on the grid (0.5) the transient decays as 0.5^k, so 200 is many orders of
# magnitude beyond where it matters and costs nothing.
AR_BURN_IN = 200

# Floor on the expected events per segment when censoring is dialled up. Two is what the
# checks need; five leaves enough headroom that the random count rarely falls below it.
MIN_EXPECTED_EVENTS_PER_SEGMENT = 5

# Segments too short to test are dropped rather than failing the whole replicate, exactly
# as `segments_from_windows` does for carved records. The realised censoring is recomputed
# from what survived, so a dropped segment shows up in the table instead of being hidden.
MIN_EVENTS_PER_SEGMENT = 2


def _weibull_unit_mean(shape: float, size: int, rng: np.random.Generator) -> np.ndarray:
    """Weibull draws scaled to mean 1, so `shape` varies dispersion and nothing else."""

    return rng.weibull(shape, size=size) / float(gamma_fn(1.0 + 1.0 / shape))


def _quantise(x: np.ndarray, mean_reads: float = QUANTUM_MEAN_READS) -> np.ndarray:
    """Round durations up onto a read grid.

    `ceil`, not `round`: a duration is a whole number of read intervals and the smallest
    observable one is a single interval, so the grid has no zero cell. That matters -
    a zero first duration puts `T_1 = 0` into eq (7)'s `ln(T_2/T_1)` and makes it singular,
    and the validator raises on any zero duration.
    """
    quantum = 1.0 / float(mean_reads)
    return np.ceil(x / quantum) * quantum


def _ar1_innovation_scale(rho: float) -> float:
    """Innovation SD giving a stationary AR(1) of unit variance."""
    if not -1.0 < rho < 1.0:
        raise ValueError(f"rho must be in (-1, 1); got {rho}")
    return float(np.sqrt(1.0 - rho**2))


def _stationary_ar1(size: int, rho: float, rng: np.random.Generator) -> np.ndarray:
    noise = rng.standard_normal(size + AR_BURN_IN)
    return lfilter([_ar1_innovation_scale(rho)], [1.0, -rho], noise)[AR_BURN_IN:]


def _segment_sizes(n: int, n_segments: int) -> list[int]:
    """Split `n` expected events into `n_segments` blocks as evenly as possible."""
    base, extra = divmod(n, n_segments)
    return [base + (1 if i < extra else 0) for i in range(n_segments)]


def _direct_segments(
    n: int,
    censoring: float,
    build: Callable[[int], tuple[np.ndarray, float]],
) -> list[Segment]:
    """Build time-truncated segments from a per-segment generator.

    `build(expected)` returns `(gaps, tau)` for one segment whose expected event count is
    `expected`. Each segment gets its own truncation time, so dialling censoring up splits
    the same total expected event count into more, shorter, independently truncated
    processes - which is exactly the multi-process structure eqs (13)-(16) address.
    """
    m = n_segments_for_censoring(n, censoring)
    segments: list[Segment] = []
    for expected in _segment_sizes(n, m):
        gaps, tau = build(expected)
        if len(gaps) < MIN_EVENTS_PER_SEGMENT:
            continue
        segments.append(Segment(x=gaps, tau=tau, n_censored_dropped=1))
    if not segments:
        raise ValueError(
            f"no segment reached {MIN_EVENTS_PER_SEGMENT} events at n={n}, "
            f"censoring={censoring}"
        )
    return segments


def n_segments_for_censoring(n: int, censoring: float) -> int:
    """How many segments realise a target censoring fraction at this `n`.

    One incomplete trailing gap per segment, so `c = m/(n+m)` and `m = c*n/(1-c)`. Capped
    so each segment EXPECTS at least `MIN_EXPECTED_EVENTS_PER_SEGMENT`; see the module
    docstring for why a smaller expected size makes most replicates unusable.
    """
    if not 0.0 <= censoring < 1.0:
        raise ValueError(f"censoring must be in [0, 1); got {censoring}")
    if censoring == 0.0:
        return 1
    target = int(round(censoring * n / (1.0 - censoring)))
    cap = max(1, n // MIN_EXPECTED_EVENTS_PER_SEGMENT)
    return int(max(1, min(cap, target)))


def _overdraw(expected: float) -> int:
    """How many gaps to draw so cumulating them is sure to pass `expected`.

    Ten standard deviations of slack at the widest dispersion on the grid, plus a floor.
    If it is ever not enough the caller raises rather than redrawing: a redraw conditioned
    on the first draw being too short would bias the retained sample toward short gaps.
    """
    return int(expected + 10.0 * np.sqrt(max(expected, 1.0)) + 50)


def _truncate_at(gaps: np.ndarray, tau: float) -> np.ndarray:
    """Keep the gaps whose cumulative sum stays strictly below `tau`."""
    event_times = np.cumsum(gaps)
    if event_times[-1] <= tau:
        raise ValueError(
            f"drew {len(gaps)} gaps totalling {event_times[-1]:.6g}, which does not reach "
            f"tau={tau:.6g}; _overdraw needs widening for this dispersion"
        )
    return gaps[: int(np.searchsorted(event_times, tau, side="left"))]


def _tau_for(expected_events: float, quantised: bool) -> float:
    """Truncation time for a unit-mean gap process expecting `expected_events` events.

    Offset by half a grid cell when quantised. On a grid every event time is a multiple of
    the quantum, so an integer `tau` can land exactly on one, leaving a zero residual and
    eq (7) singular; half a cell off the grid makes that impossible while shifting the
    truncation time by a sixth of a mean duration.
    """
    tau = float(expected_events)
    if quantised:
        tau += 0.5 / QUANTUM_MEAN_READS
    return tau


def arm_a(
    n: int,
    rng: np.random.Generator,
    *,
    shape: float = 1.0,
    quantised: bool = False,
    censoring: float = 0.0,
) -> list[Segment]:
    """iid Weibull durations. The null for every check in the battery."""

    def build(expected: int) -> tuple[np.ndarray, float]:
        tau = _tau_for(expected, quantised)
        gaps = _weibull_unit_mean(shape, _overdraw(expected), rng)
        if quantised:
            gaps = _quantise(gaps)
        return _truncate_at(gaps, tau), tau

    return _direct_segments(n, censoring, build)


def arm_d(
    n: int,
    rng: np.random.Generator,
    *,
    shape: float = 1.0,
    b: float = 1.0,
    quantised: bool = False,
    censoring: float = 0.0,
) -> list[Segment]:
    """Weibull trend-renewal process with a power-law trend.

    The trend intensity is `lambda(t) = b * t^(b-1)`, so the cumulative intensity is
    `Lambda(t) = t^b`. A TRP is defined by `Lambda(T_i)` forming a renewal process, so the
    generator is: draw unit-mean renewal gaps `V_i`, cumulate to `S_i`, and invert -
    `T_i = S_i^(1/b)`. Durations are the differences.

    `b = 1` is the identity transform, i.e. no trend at all, so it reproduces Arm A exactly
    and doubles as a size row. `b < 1` stretches later durations (improving); `b > 1`
    compresses them (degrading). This is the alternative C1 and C2 are actually directed
    at - without it, criterion 2 would be unevaluable for both.
    """
    if b <= 0.0:
        raise ValueError(f"the power-law exponent b must be positive; got {b}")
    if quantised:
        raise ValueError(
            "Arm D has no quantised variant: quantising the durations would change "
            "their sum and so break the truncation time the trend transform defines. "
            "Arm A carries the quantisation factor."
        )

    def build(expected: int) -> tuple[np.ndarray, float]:
        # Truncation happens in TRANSFORMED time, where the process is a plain renewal
        # process, and both the events and tau are then mapped back. Truncating in real
        # time instead would make the retained count depend on b through the transform as
        # well as through the trend, confounding the two.
        tau_transformed = float(expected)
        renewal = np.cumsum(_weibull_unit_mean(shape, _overdraw(expected), rng))
        kept = renewal[renewal < tau_transformed]
        if len(kept) < MIN_EVENTS_PER_SEGMENT:
            return np.zeros(0, dtype=float), 1.0
        event_times = kept ** (1.0 / b)
        gaps = np.diff(np.concatenate([[0.0], event_times]))
        return gaps, float(tau_transformed ** (1.0 / b))

    return _direct_segments(n, censoring, build)


def _carved_segments(
    values: np.ndarray, clock: str, min_events: int
) -> tuple[list[Segment], int]:
    """Carve a read series and map it onto segments.

    The observation interval is `[0, n_reads*dt)` - one full read interval PAST the last
    read - not `[0, t_last]`. That distinction is the difference between a truncation time
    chosen in advance and one read off the data: ending at the last read makes `tau` an
    event-determined boundary, and when the final window is a single read it collapses to
    `tau == T_N`, where eq (7) is singular. Measured: that hits 25% of Arm B
    replicates on the calendar clock.
    """
    n_reads = len(values)
    t = np.arange(n_reads, dtype=float) * READ_DT_S
    threshold = float(np.quantile(values, 1.0 - IN_SPEC_P))
    windows = carve_windows(t, values, threshold, big_values_good=True)
    return segments_from_windows(
        windows,
        clock=clock,
        min_events=min_events,
        observation_end_s=n_reads * READ_DT_S,
    )


def n_reads_for_events(n: int) -> int:
    """Reads needed to expect about `n` complete windows.

    An iid read series in spec with probability `p` starts a new in-spec run at each
    out-to-in transition, so the expected count over `N` reads is `N*p*(1-p)`; at p = 0.5
    that is `N/4`. The `+ 8` is slack so a short record still clears the minimum event
    count after the trailing censored window is removed.
    """
    return int(n / (IN_SPEC_P * (1.0 - IN_SPEC_P))) + 8


def arm_b(
    n: int,
    rng: np.random.Generator,
    *,
    clock: str = CLOCK_IN_SPEC,
    min_events: int = 2,
) -> tuple[list[Segment], int]:
    """iid Gaussian reads, thresholded at the occupancy-matched quantile, then carved.

    Under this generator the carved durations are iid geometric run lengths, so the null
    holds on BOTH clocks: successive in-spec runs are iid, and so are successive
    birth-to-birth intervals (an in-spec run plus an out-of-spec run). Being intrinsically
    quantised - every duration is a whole number of read intervals, most of them one - it
    is also the arm where tie handling in the permutation p-value has to work.
    """
    values = rng.standard_normal(n_reads_for_events(n))
    return _carved_segments(values, clock, min_events)


def arm_c(
    n: int,
    rng: np.random.Generator,
    *,
    rho: float = 0.0,
    clock: str = CLOCK_IN_SPEC,
    min_events: int = 2,
) -> tuple[list[Segment], int]:
    """AR(1) reads at read-level correlation `rho`, thresholded and carved identically.

    `z_t = rho*z_{t-1} + sqrt(1 - rho^2)*e_t`, stationary with unit variance, so `rho` moves
    the dependence and nothing else - the marginal distribution and hence the threshold and
    the occupancy are unchanged. At `rho = 0` this is Arm B with a different random path,
    which is what makes the two comparable as a self-consistency check.

    Read-level `rho` is NOT duration-level dependence. So the runner reports the INDUCED
    duration-level lag-1 for every cell, and the power curve is read against that, not
    against `rho`.
    """
    if not -1.0 < rho < 1.0:
        raise ValueError(f"rho must be in (-1, 1); got {rho}")
    values = _stationary_ar1(n_reads_for_events(n), rho, rng)
    return _carved_segments(values, clock, min_events)


def arm_e(
    n: int,
    rng: np.random.Generator,
    *,
    shape: float = 1.0,
    rho: float = 0.0,
    quantised: bool = False,
    censoring: float = 0.0,
) -> list[Segment]:
    """Weibull durations with a Gaussian-copula AR(1) dependence, direct.

    Draw a stationary latent AR(1) `w`, map it through the normal CDF to uniforms, and
    push those through the Weibull quantile function. The result has EXACTLY Arm A's
    marginal distribution - so at `rho = 0` this arm is Arm A - while `rho` moves the rank
    dependence and nothing else. For a Gaussian copula the induced Spearman correlation is
    `(6/pi)*arcsin(rho/2)`, which is within 4% of `rho` itself over the grid used here, so
    the swept parameter is directly comparable to a duration-level lag-1.

    **Why this arm exists.** It was not in the plan; Arm C was to be the power arm for
    C5 and C6. Measurement killed that: read-level AR(1) induces NO positive duration-level
    dependence at any rho on the grid, because the level crossings of a stationary
    Gaussian process regenerate and successive excursion lengths are very nearly a
    renewal process. Arm C is kept exactly as specified because that null result is worth
    reporting, but without a generator that actually varies duration dependence,
    criterion 2 would be unevaluable for half the assessed checks.
    """

    def build(expected: int) -> tuple[np.ndarray, float]:
        tau = _tau_for(expected, quantised)
        latent = _stationary_ar1(_overdraw(expected), rho, rng)
        uniforms = stats.norm.cdf(latent)
        # Weibull(shape) quantile function, then scaled to unit mean exactly as
        # _weibull_unit_mean does, so the marginal matches Arm A term for term.
        gaps = (-np.log1p(-uniforms)) ** (1.0 / shape) / float(
            gamma_fn(1.0 + 1.0 / shape)
        )
        if quantised:
            gaps = _quantise(gaps)
        return _truncate_at(gaps, tau), tau

    return _direct_segments(n, censoring, build)
