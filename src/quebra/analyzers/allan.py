"""Overlapping Allan deviation of a gapped Ramsey record, pooled over its segments.

The overlapping estimator assumes equally spaced reads. A gapped record is cut into
segments under the carve's gap rule (`windows.observed_reads`), and every segment is
treated as sampled at ONE spacing, tau0, the median positive spacing of the whole record.
At each averaging factor m (tau = m * tau0) the second differences of every segment long
enough to hold one are pooled, so each segment weighs in by its number of terms. No term
spans a gap.

Two exponents are easy to confuse and are named apart here. The local slope of
log ADEV against log tau is mu / 2 in the convention AVAR ∝ tau^mu. The spectral exponent
alpha of S_y(f) ∝ f^alpha is a different number (mu = -alpha - 1 for -3 < alpha < 1).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from quebra.analyzers.windows import DEFAULT_GAP_MULT, observed_reads

TAUS_MODES = ("all", "octave")


@dataclass(slots=True)
class AllanResult:
    modes: dict[str, pd.DataFrame]
    meta: dict[str, object]
    fractional_adev: np.ndarray | None = None
    # The declared operating frequency the fractional series divides by; None when no
    # fractional series was asked for.
    nu0_hz: float | None = None
    diagnostics: dict[str, object] = field(default_factory=dict)


def check_nu0_hz(nu0_hz: object) -> float:
    """Return a declared operating frequency as a float, or raise.

    A non-positive or non-finite value has no fractional meaning, and dividing by it would
    give a series of the right shape and the wrong size.
    """
    try:
        value = float(nu0_hz)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise TypeError(f"nu0_hz must be a number of hertz, got {nu0_hz!r}") from exc
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"nu0_hz must be finite and positive, got {value!r}")
    return value


def _segment_terms(y: np.ndarray, m: int) -> tuple[float, int]:
    """Sum of squared second differences of one segment's phase at factor m, and their
    count. In units of y squared times tau0 squared, so tau0 is applied once, at pooling.

    Phase in units of tau0 is the running sum of y, starting at 0; the second difference
    x[j + 2m] - 2 x[j + m] + x[j] is m times the difference of two adjacent m-read means.
    """
    n_terms = len(y) - 2 * m + 1
    if n_terms < 1:
        return 0.0, 0
    x = np.concatenate(([0.0], np.cumsum(y)))
    d2 = x[2 * m :] - 2.0 * x[m:-m] + x[: -2 * m]
    return float(np.sum(d2 * d2)), n_terms


def pooled_overlapping_adev(
    segments: list[np.ndarray], tau0_s: float, min_terms: int, taus_mode: str
) -> pd.DataFrame:
    """Overlapping ADEV at tau = m * tau0, pooled over segments by term count.

    AVAR(m) = sum of squared second differences over every segment /
    (2 m^2 tau0^2 * total terms), the overlapping estimator of each segment weighted by its
    number of terms. A tau is kept when its pooled term count is at least `min_terms`.
    `adev_err` is adev / sqrt(terms), the term-count standard error; the overlapping terms
    are correlated, so it understates the uncertainty, and it is not a confidence interval.
    """
    if taus_mode not in TAUS_MODES:
        raise ValueError(f"taus_mode must be one of {TAUS_MODES}, got {taus_mode!r}")
    longest = max((len(seg) for seg in segments), default=0)
    m_max = longest // 2
    if taus_mode == "all":
        factors = list(range(1, m_max + 1))
    else:
        factors = [2**k for k in range(int(math.log2(m_max)) + 1)] if m_max else []

    rows = []
    for m in factors:
        total, terms, used = 0.0, 0, 0
        for seg in segments:
            s, n = _segment_terms(seg, m)
            if n:
                total, terms, used = total + s, terms + n, used + 1
        if terms < min_terms:
            continue
        # tau0 cancels between the phase (x = tau0 * cumsum(y)) and the 1 / (m tau0)^2.
        adev = math.sqrt(total / (2.0 * m * m * terms))
        rows.append(
            {
                "m": m,
                "tau_s": m * tau0_s,
                "adev": adev,
                "adev_err": adev / math.sqrt(terms),
                "n_terms": terms,
                "n_segments": used,
            }
        )
    return pd.DataFrame(
        rows, columns=["m", "tau_s", "adev", "adev_err", "n_terms", "n_segments"]
    )


def _transform_signal(freq_hz: np.ndarray, mode: str) -> np.ndarray:
    """The series a mode computes on. One mode exists: `raw`, the detuning as read."""
    if mode == "raw":
        return np.asarray(freq_hz, dtype=float)
    raise ValueError(f"Unknown Allan mode '{mode}'; the only mode is 'raw'.")


def run(
    norm: Mapping[str, object],
    config: Mapping[str, object],
    *,
    gap_mult: float = DEFAULT_GAP_MULT,
    fractional: bool = False,
    nu0_hz: float | None = None,
) -> AllanResult:
    """Overlapping Allan deviation of `delta_hz`, per segment and pooled.

    Time in seconds, frequency in hertz. `fractional=True` also returns ADEV / `nu0_hz`,
    with `nu0_hz` the declared operating frequency: one fixed value, never a read's own
    frequency. `fractional` and `nu0_hz` must be given together.

    Raises on a repeated timestamp inside a segment (two reads at one time have no spacing
    for the equal-spacing estimator), on `fractional` without a valid `nu0_hz` or the
    other way round, and when no tau keeps `min_points_for_allan` pooled terms.

    Validity assumptions:

    - Assumption: reads within a segment are equally spaced at tau0.
      Diagnostic: `diagnostics["spacing_cv_within_segments"]`.
      Consequence of violation: tau = m * tau0 is the nominal averaging time only; uneven
      spacing mixes averaging times and smears the slope.
      Reference: NIST SP 1065 (Riley), the overlapping Allan variance, which
      `allantools.oadev` cites as eqn (11), p. 16.
    - Assumption: the noise process is the same in every segment.
      Diagnostic: none shipped; `n_segments` per tau says how many segments a point pools.
      Consequence of violation: the pooled deviation is a term-weighted average, led by
      the longest segments.
      Reference: no source located.
    - Assumption: the terms behind `adev_err` are independent.
      Diagnostic: none; overlapping terms share reads, so they are not.
      Consequence of violation: `adev_err` understates the uncertainty, more at large m.
      Reference: `allantools` reports the same form (`calc_adev_phase`); confidence
      intervals with equivalent degrees of freedom are not computed here.
    """
    if not isinstance(norm, Mapping):
        raise TypeError("allan.run expects normalized dataset mapping input.")
    if "t_rel_s" not in norm:
        raise KeyError(
            "allan.run requires 't_rel_s' (relative seconds) in normalized mapping"
        )
    if "delta_hz" not in norm:
        raise KeyError("allan.run requires 'delta_hz' (hertz) in normalized mapping")
    if fractional and nu0_hz is None:
        raise ValueError(
            "fractional=True needs nu0_hz, the declared operating frequency in hertz; "
            "the records carry no absolute frequency to divide by"
        )
    if nu0_hz is not None and not fractional:
        raise ValueError("nu0_hz was given but fractional=False; nothing would read it")
    nu0 = check_nu0_hz(nu0_hz) if fractional else None

    allan_cfg = config.get("allan", {})
    if not isinstance(allan_cfg, Mapping):
        raise TypeError("config['allan'] must be a mapping")
    modes = list(allan_cfg.get("modes", []))
    if len(modes) == 0:
        raise ValueError("config['allan']['modes'] must be a non-empty list")
    min_terms = int(allan_cfg.get("min_points_for_allan", 8))
    taus_mode = str(allan_cfg.get("taus_mode", "all"))

    obs = observed_reads(norm["t_rel_s"], norm["delta_hz"], gap_mult=gap_mult)
    segments = obs.segments()
    within = np.diff(obs.t_rel_s)[obs.within_segment()]
    if np.any(within == 0.0):
        raise ValueError(
            f"{int(np.sum(within == 0.0))} repeated timestamp(s) inside a segment; the "
            f"overlapping Allan estimator needs every within-segment spacing positive"
        )
    tau0_s = obs.median_spacing_s
    meta_in = (
        dict(norm.get("meta", {})) if isinstance(norm.get("meta", {}), Mapping) else {}
    )

    summary: dict[str, pd.DataFrame] = {}
    for mode in modes:
        signal = _transform_signal(obs.values, mode=str(mode))
        frame = pooled_overlapping_adev(
            [signal[s:e] for s, e in segments], tau0_s, min_terms, taus_mode
        )
        if frame.empty:
            raise ValueError(
                f"no tau keeps {min_terms} pooled terms: {len(segments)} segment(s), the "
                f"longest of {max(e - s for s, e in segments)} reads"
            )
        summary[str(mode)] = frame
        print(
            f"[allan] mode={mode} tau_points={len(frame)} tau0={tau0_s:.6g}s "
            f"segments={len(segments)}",
            flush=True,
        )

    fractional_adev: np.ndarray | None = None
    if nu0 is not None:
        first = summary[str(modes[0])]
        fractional_adev = first["adev"].to_numpy(dtype=float) / nu0

    return AllanResult(
        modes=summary,
        meta={
            "tau0_s": float(tau0_s),
            "dataset_id": meta_in.get("dataset_id"),
            "run_name": meta_in.get("run_name"),
            "qubit": meta_in.get("qubit"),
        },
        fractional_adev=fractional_adev,
        nu0_hz=nu0,
        diagnostics={
            "gap_mult": float(gap_mult),
            "gap_threshold_s": obs.gap_threshold_s,
            "n_segments": len(segments),
            # A segment of one read holds no difference at any tau.
            "n_segments_too_short": sum(1 for s, e in segments if e - s < 2),
            "n_dropped_nonfinite": obs.n_dropped_nonfinite,
            "min_terms_per_tau": min_terms,
            "spacing_cv_within_segments": float(np.std(within) / np.mean(within))
            if len(within)
            else float("nan"),
        },
    )
