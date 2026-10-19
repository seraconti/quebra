"""Telegraph (two-level fluctuator) diagnostics: a two-lobe fit and the dwell dynamics.

The lobes are a 1- versus 2-component Gaussian mixture on the values. The dynamics assign
each read to the lower or the upper lobe and measure the two states' dwells on the reads
the record observed, under the carve's gap rule (`windows.observed_reads`):

- **exposure** of a state is the sum of within-segment forward spacings whose LEFT read is
  in that state, the convention the Kaplan-Meier carve uses. A spacing across a gap is
  charged to no state;
- **exits** of a state are its observed changes to the other state within a segment. A state
  change across a gap was not observed and is neither an exit nor a transition;
- **mean dwell** is exposure / exits;
- **dwell CV** is over COMPLETE runs only, those whose start and end were both observed. A
  run cut by a gap or by either end of the record is counted in `n_cut_runs_state*`;
- **switching rate** is the within-segment state changes over the total exposure.

Two reads at one time inside a segment raise: a state change between them would be an exit
after zero exposure.

`values` are in hertz (the Ramsey recipe passes `raw_frequency_hz` or `delta_hz`), so the
lobe separation and spread are in hertz. Fractional units would need a declared operating
frequency, which this analyzer does not take.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from sklearn.mixture import GaussianMixture

from quebra.analyzers.windows import (
    BIRTH_UP_CROSSING,
    DEATH_DOWN_CROSSING,
    DEFAULT_GAP_MULT,
    carve,
    observed_reads,
)


@dataclass(slots=True)
class TLFResult:
    # Group 1 - histogram-based
    is_bimodal: bool
    # None when the 2-component fit raised. A number here always means the comparison happened.
    bic_delta: float | None
    normalized_bic_delta: float | None
    lobe_separation_hz: float | None
    # weighted sqrt of per-lobe variances - reflects within-state spread, NOT a clean noise floor; confounded by drift and unresolved sub-fluctuators.
    within_lobe_spread_hz: float | None
    lobe_snr: float | None
    gmm1: Any
    gmm2: Any

    # Group 2 - dynamics. State 0 is the lower lobe, state 1 the upper.
    n_transitions: int | None
    switching_rate_per_hour: float | None
    mean_dwell_state0_s: float | None
    mean_dwell_state1_s: float | None
    dwell_cv_state0: float | None
    dwell_cv_state1: float | None
    exposure_state0_s: float | None
    exposure_state1_s: float | None
    n_exits_state0: int | None
    n_exits_state1: int | None
    n_complete_runs_state0: int | None
    n_complete_runs_state1: int | None
    n_cut_runs_state0: int | None
    n_cut_runs_state1: int | None

    diagnostics: dict[str, object] = field(default_factory=dict)
    # Read this before `is_bimodal`: False there means no comparison was made, not that one
    # lobe was established.
    fit_failed: bool = False


# Every dynamics field, None together when the 2-component fit failed.
_DYNAMICS_FIELDS = (
    "n_transitions",
    "switching_rate_per_hour",
    *(
        f"{name}_state{k}{unit}"
        for name, unit in (
            ("mean_dwell", "_s"),
            ("dwell_cv", ""),
            ("exposure", "_s"),
            ("n_exits", ""),
            ("n_complete_runs", ""),
            ("n_cut_runs", ""),
        )
        for k in (0, 1)
    ),
)


@dataclass(slots=True)
class _Dwells:
    exposure_s: float
    exits: int
    complete_s: list[float]
    n_cut: int


def _state_dwells(t: np.ndarray, in_state: np.ndarray, is_gap: np.ndarray) -> _Dwells:
    """Exposure, exits and complete-run durations of one state, from the carve's runs."""
    exposure = 0.0
    exits = 0
    complete: list[float] = []
    n_runs = 0
    for run in carve(t, in_state, is_gap):
        s, e = int(run["s"]), int(run["e"])
        n_runs += 1
        exit_observed = run["death_type"] == DEATH_DOWN_CROSSING
        # An observed exit closes the run at the first read of the other state, inside the
        # segment; a cut run's last forward spacing crosses a gap or does not exist.
        duration = float(t[e] - t[s]) if exit_observed else float(t[e - 1] - t[s])
        exposure += duration
        if exit_observed:
            exits += 1
            if run["birth_type"] == BIRTH_UP_CROSSING:
                complete.append(duration)
    return _Dwells(
        exposure_s=exposure,
        exits=exits,
        complete_s=complete,
        n_cut=n_runs - len(complete),
    )


def _dwell_cv(durations: list[float]) -> float | None:
    # Defined only from three complete runs: the sample sd of two is too unstable to quote.
    if len(durations) < 3:
        return None
    return float(np.std(durations, ddof=1) / np.mean(durations))


def run(
    values: np.ndarray,
    timestamps: np.ndarray,
    *,
    seed: int,
    gap_mult: float = DEFAULT_GAP_MULT,
) -> TLFResult:
    """Fit 1- and 2-component GMMs and return the TLF diagnostics.

    `seed` is required, not defaulted. `GaussianMixture` initialises by k-means, so an
    unseeded fit makes `bic_delta`, `is_bimodal`, the state assignment and every dwell
    statistic a fresh random variable per call while the run identity stays byte-identical.
    It arrives as a step argument, like `xi_seed`, so it reaches the provenance label.
    `gap_mult` sets the gap rule, as in `windows.run`; the Ramsey recipe passes it as a step
    argument too.

    Non-finite (time, value) pairs are dropped first and counted in `diagnostics`. A
    2-component fit that raises ValueError or LinAlgError sets `fit_failed`, records the
    error in `diagnostics["fit_error"]` and leaves `bic_delta`, `gmm2` and every dynamics
    field as None rather than reporting a comparison that did not happen. A repeated
    timestamp inside a segment raises ValueError.

    Validity assumptions:

    - Assumption: the state process is a two-state Markov chain in continuous time, leaving
      state 0 at rate lambda0 and state 1 at rate lambda1.
      Diagnostic: the dwell CV over complete runs, which is 1 for exponential dwells when
      the reads are dense and the segments long compared with the dwells. Complete runs
      are biased toward short dwells otherwise, so a CV below 1 is not by itself evidence
      against exponential dwells.
      Consequence of violation: exposure / exits stays the ratio estimate of a mean dwell
      under a stationary renewal process. Even under the chain it is the maximum-likelihood
      estimate of 1 / lambda0 only for a continuously observed path; with discrete reads it
      is the dense-read approximation of the next assumption.
      Reference: the censored-exponential maximum-likelihood estimate, total time at risk
      over the number of events; no source located.
    - Assumption: the reads are dense on the scale of the dwells, the chain is stationary
      (the probability that a read is in state 0 is the same at every read), and the
      spacings do not depend on the state.
      Diagnostic: `diagnostics["effective_spacing_s"]`, E[delta^2]/E[delta] over the
      within-segment spacings, times (lambda0 + lambda1).
      Consequence of violation: switches missed between reads bias the mean dwell UP and
      the switching rate DOWN. For the chain above, P(the right read is in state 1 | the
      left read is in state 0) = lambda0 delta - lambda0 (lambda0 + lambda1) delta^2 / 2 +
      O(delta^3). So, as ratios of expectations with the moments taken over the
      within-segment spacings, mean dwell = (1 / lambda0) (1 + (lambda0 + lambda1)
      E[delta^2] / (2 E[delta])) to first order, symmetrically for state 1, and the
      switching rate = (2 lambda0 lambda1 / (lambda0 + lambda1)) (1 - (lambda0 + lambda1)
      E[delta^2] / (2 E[delta])).
      Reference: derived here from the two-state transition probability; no source located.
    - Assumption: each read is assigned its true state.
      Diagnostic: `lobe_snr`.
      Consequence of violation: one misassigned read splits a run in three, adding two exits
      and biasing the mean dwell DOWN.
      Reference: no source located.
    """
    if timestamps is None:
        raise ValueError(
            "timestamps must be provided to TLF.run for dynamics computation"
        )

    if np.asarray(values).size == 0:
        raise ValueError("No values provided to TLF analysis")
    if np.asarray(timestamps).shape[0] != np.asarray(values).shape[0]:
        raise ValueError("timestamps must have the same length as values")
    obs = observed_reads(timestamps, values, gap_mult=gap_mult)
    vals = obs.values.reshape(-1, 1)
    if vals.size == 0:
        raise ValueError("No finite (time, value) pair provided to TLF analysis")
    within = np.diff(obs.t_rel_s)[obs.within_segment()]
    if np.any(within == 0.0):
        raise ValueError(
            f"{int(np.sum(within == 0.0))} repeated timestamp(s) inside a segment; a state "
            f"change between two reads at one time would be an exit after zero exposure"
        )
    diagnostics: dict[str, object] = {
        "n_reads": int(vals.shape[0]),
        "n_dropped_nonfinite": obs.n_dropped_nonfinite,
        "n_nonpositive_steps": obs.n_nonpositive_steps,
        "median_spacing_s": obs.median_spacing_s,
        "effective_spacing_s": obs.effective_spacing_s(),
        "gap_threshold_s": obs.gap_threshold_s,
        "n_gaps": obs.n_gaps,
        "gap_mult": float(gap_mult),
    }

    gmm1 = GaussianMixture(n_components=1, covariance_type="full", random_state=seed)
    gmm1.fit(vals)

    try:
        gmm2 = GaussianMixture(
            n_components=2, covariance_type="full", random_state=seed
        )
        gmm2.fit(vals)
    except (ValueError, np.linalg.LinAlgError) as exc:
        diagnostics["fit_error"] = f"{type(exc).__name__}: {exc}"
        return TLFResult(
            is_bimodal=False,
            bic_delta=None,
            normalized_bic_delta=None,
            lobe_separation_hz=None,
            within_lobe_spread_hz=None,
            lobe_snr=None,
            gmm1=gmm1,
            gmm2=None,
            **dict.fromkeys(_DYNAMICS_FIELDS),
            diagnostics=diagnostics,
            fit_failed=True,
        )

    bic_delta = float(gmm1.bic(vals) - gmm2.bic(vals))
    normalized_bic_delta = float(bic_delta / len(vals))
    is_bimodal = bic_delta > 6.0

    means = np.sort(gmm2.means_.ravel())
    weights = gmm2.weights_.ravel()
    covs = gmm2.covariances_.ravel()
    lobe_separation_hz = float(abs(means[-1] - means[0]))
    # approximate within-lobe spread as sqrt(weighted average of variances)
    within_lobe_spread_hz = float(np.sqrt(np.sum(weights * covs)))
    lobe_snr = (
        float(lobe_separation_hz / within_lobe_spread_hz)
        if within_lobe_spread_hz > 0.0
        else None
    )

    # State 0 is the component with the LOWER mean, whatever order the fit returned.
    lower_idx = int(np.argmin(gmm2.means_.ravel()))
    state = np.where(gmm2.predict(vals).ravel() == lower_idx, 0, 1)
    t = obs.t_rel_s
    low = _state_dwells(t, state == 0, obs.is_gap)
    high = _state_dwells(t, state == 1, obs.is_gap)

    n_transitions = low.exits + high.exits
    total_exposure_s = low.exposure_s + high.exposure_s
    switching_rate_per_hour = (
        float(n_transitions / (total_exposure_s / 3600.0))
        if total_exposure_s > 0.0
        else None
    )

    return TLFResult(
        is_bimodal=is_bimodal,
        bic_delta=bic_delta,
        normalized_bic_delta=normalized_bic_delta,
        lobe_separation_hz=lobe_separation_hz,
        within_lobe_spread_hz=within_lobe_spread_hz,
        lobe_snr=lobe_snr,
        gmm1=gmm1,
        gmm2=gmm2,
        n_transitions=n_transitions,
        switching_rate_per_hour=switching_rate_per_hour,
        mean_dwell_state0_s=low.exposure_s / low.exits if low.exits else None,
        mean_dwell_state1_s=high.exposure_s / high.exits if high.exits else None,
        dwell_cv_state0=_dwell_cv(low.complete_s),
        dwell_cv_state1=_dwell_cv(high.complete_s),
        exposure_state0_s=low.exposure_s,
        exposure_state1_s=high.exposure_s,
        n_exits_state0=low.exits,
        n_exits_state1=high.exits,
        n_complete_runs_state0=len(low.complete_s),
        n_complete_runs_state1=len(high.complete_s),
        n_cut_runs_state0=low.n_cut,
        n_cut_runs_state1=high.n_cut,
        diagnostics=diagnostics,
    )
