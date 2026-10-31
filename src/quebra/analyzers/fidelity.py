from __future__ import annotations

from dataclasses import dataclass, field
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from quebra.analyzers.kaplan_meier import KaplanMeierSet
from quebra.analyzers.windows import check_k
from quebra.core.types import LOAD_COUNT_KEYS


@dataclass(slots=True)
class FidelityResult:
    frame: pd.DataFrame
    meta: dict[str, object]
    diagnostics: dict[str, object] = field(default_factory=dict)


@dataclass(slots=True)
class FidelityInputs:
    t_rel_s: np.ndarray
    raw_frequency_hz: np.ndarray
    frequency_error_hz: np.ndarray
    rabi_hz: np.ndarray
    profile: str
    use_angular_frequency: bool = False
    dataset_id: str = ""
    rabi_source: str = ""
    rabi_declared_key: str = ""
    # Rows the load dropped before this step (LOAD_COUNT_KEYS), copied onto the result.
    load_counts: dict[str, int] = field(default_factory=dict)


# The three frequencies x is built from, each in cyclic hertz, the unit their `_hz` keys
# declare (AGENTS.md section 3). x = Delta / Omega is a ratio, so all three must share it.
FREQUENCY_KEYS = ("raw_frequency_hz", "frequency_error_hz", "rabi_hz")


def make_inputs_from_norm(
    norm: Mapping[str, object], config: Mapping[str, object]
) -> FidelityInputs:
    """Extract FidelityInputs from a normalized mapping and config dict.

    Raises ValueError when the norm carries a Rabi frequency or a frequency under an
    angular key (`*_rad_s`), naming it: x divides one frequency by another, and a
    factor 2 pi between them would move every label. Units are declared by key names in
    this repo, so a value given in another unit under a cyclic-hertz key is not caught here.
    """
    if "t_rel_s" not in norm:
        raise KeyError(
            "Fidelity requires 't_rel_s' (relative seconds) in the normalized mapping."
        )
    angular = sorted(
        key
        for key in norm
        if isinstance(key, str)
        and key.endswith("_rad_s")
        and ("rabi" in key or "frequency" in key)
    )
    if angular:
        raise ValueError(
            f"Fidelity reads its frequencies in cyclic hertz ({', '.join(FREQUENCY_KEYS)}), "
            f"but the normalized mapping carries {angular} in radians per second; x = "
            f"Delta / Omega needs one unit for all three"
        )
    if "rabi_hz" not in norm:
        raise KeyError(
            "Fidelity requires 'rabi_hz' in the normalized mapping. Declare the record's "
            "Rabi frequency (configure_ramsey_job(..., rabi_table=Dataset(...))) "
            "or disable fidelity for this job."
        )
    for key in ("raw_frequency_hz", "frequency_error_hz"):
        if key not in norm:
            raise KeyError(
                f"Fidelity requires {key!r} in the normalized mapping: x is the detuning "
                f"from the first read over the Rabi frequency, and its sigma the read's "
                f"frequency error over the same Rabi frequency."
            )
    fidelity_cfg = config.get("fidelity", {})
    if not isinstance(fidelity_cfg, Mapping):
        raise TypeError("config['fidelity'] must be a mapping")
    meta = norm.get("meta", {}) if isinstance(norm.get("meta", {}), Mapping) else {}
    dataset_id = str(meta.get("dataset_id") or "").strip()
    if not dataset_id:
        raise KeyError(
            "Fidelity requires a non-empty meta['dataset_id'] in the normalized mapping: "
            "the panel caption names the record it draws."
        )
    profile = str(meta.get("profile", config.get("dataset_profile", "")))
    # The artifact must say where its Rabi frequency came from, so a value with no named
    # source is refused rather than recorded as "" or "missing".
    rabi_source = str(meta.get("rabi_source") or "").strip()
    if rabi_source in ("", "missing"):
        raise KeyError(
            "Fidelity requires meta['rabi_source'] naming where 'rabi_hz' came from; the "
            "declared join sets it to 'declared'."
        )
    return FidelityInputs(
        t_rel_s=np.asarray(norm["t_rel_s"], dtype=float),
        raw_frequency_hz=np.asarray(norm["raw_frequency_hz"], dtype=float),
        frequency_error_hz=np.asarray(norm["frequency_error_hz"], dtype=float),
        rabi_hz=np.asarray(norm["rabi_hz"], dtype=float),
        profile=profile,
        use_angular_frequency=bool(fidelity_cfg.get("use_angular_frequency", False)),
        dataset_id=dataset_id,
        rabi_source=rabi_source,
        rabi_declared_key=str(meta.get("rabi_declared_key", "")),
        load_counts={key: int(meta[key]) for key in LOAD_COUNT_KEYS if key in meta},
    )


def _to_angular_frequency(
    freq_hz: np.ndarray, use_angular_frequency: bool
) -> np.ndarray:
    frequency_hz = np.asarray(freq_hz, dtype=float)
    if use_angular_frequency:
        return 2.0 * np.pi * frequency_hz
    return frequency_hz


def gate_fidelity(delta: np.ndarray, rabi_drive: np.ndarray) -> np.ndarray:
    """Average gate fidelity of an X_pi/2 rotation under a quasi-static detuning `delta`.

    F = (1/3) [1 + Omega^2/lambda^2 + (Delta^2/lambda^2) cos^2(lambda pi / 4 Omega)
               + (Omega/lambda) sin(lambda pi / 2 Omega)],  lambda = sqrt(Delta^2 + Omega^2),

    transcribed from Capannelli et al., "Tracking spin qubit frequency variations over 912
    days", Supplementary notes, section V, Eq. (9) (lambda from Eq. (6)). It depends on
    `delta` and `rabi_drive` only through their ratio, so both must share one unit, hertz or
    radians per second alike.

    Validity assumptions:
    - Assumption: the detuning is constant over the rotation (quasi-static) and the rotation
      is otherwise ideal (Supplementary section V). Diagnostic: none computed here. Consequence
      of violation: the fidelity is optimistic. Reference: Supplementary notes, section V.
    """
    delta = np.asarray(delta, dtype=float)
    rabi_drive = np.asarray(rabi_drive, dtype=float)
    lam = np.sqrt(delta**2 + rabi_drive**2)
    term_1 = 1.0
    term_2 = (rabi_drive**2) / (lam**2)
    term_3 = (delta**2) / (lam**2) * np.cos((lam * np.pi) / (4.0 * rabi_drive)) ** 2
    term_4 = (rabi_drive / lam) * np.sin((lam * np.pi) / (2.0 * rabi_drive))
    return (1.0 / 3.0) * (term_1 + term_2 + term_3 + term_4)


# ---------------------------------------------------------------------------
# The infidelity as a function of x = Delta / Omega
#
# Eq. (9) rearranges, with u = sqrt(1 + x^2), to F = 1/3 + (u cos(u pi/4) + sin(u pi/4))^2
# / (3 u^2). So I = 1 - F <= 2/3 everywhere, with equality exactly where
# tan(u pi/4) = -u, and I rises strictly from 0 at x = 0 to 2/3 at the first such point,
# PEAK_X. Past it I stays above 0.315 (Cauchy-Schwarz gives I >= 1/3 - 1/(3 u^2), above
# 0.315 for x >= 4.2, and a grid covers the rest; both in tests/test_infidelity_eq9.py).
# For any rung T below MAX_RUNG_INFIDELITY, "I <= T" is therefore exactly "|x| <= x_T".
# ---------------------------------------------------------------------------

# sqrt(u*^2 - 1), u* the root of u cos(u pi/4) + sin(u pi/4) on (2, 4).
PEAK_X = 2.2768823885676457
# The largest rung for which "I <= T" is an interval in |x|: below the smallest I past PEAK_X.
MAX_RUNG_INFIDELITY = 0.315
# The smallest rung x_threshold accepts. 1 - F cancels as I falls, so x_T loses relative
# precision: it agrees with a cancellation-free form of Eq. (9) to 1e-12 or better on the
# shipped rungs (1e-2 to 1e-4) and to 1e-9 at 1e-7, both pinned by tests.
MIN_RUNG_INFIDELITY = 1e-7


def infidelity_at(x: np.ndarray | float) -> np.ndarray:
    """I(x) = 1 - F(x, 1), clipped at 0: the one infidelity of x the helpers below share.

    `x` is Delta / Omega with both in one unit; Eq. (9) depends on nothing else.
    """
    x = np.asarray(x, dtype=float)
    return np.clip(1.0 - gate_fidelity(x, np.ones_like(x)), 0.0, None)


def x_threshold(threshold: float) -> float:
    """x_T >= 0 with I(x_T) = threshold: a rung T on the infidelity becomes x_T on |x|.

    Raises ValueError unless MIN_RUNG_INFIDELITY <= threshold < MAX_RUNG_INFIDELITY. Above
    that range "I <= T" stops being an interval in |x|; below it, the cancellation in
    1 - F makes x_T imprecise. The root is bracketed on [0, PEAK_X], where I is strictly
    increasing, and found to a relative tolerance of 4 machine epsilons.
    """
    threshold = float(threshold)
    if not MIN_RUNG_INFIDELITY <= threshold < MAX_RUNG_INFIDELITY:
        raise ValueError(
            f"x_threshold needs {MIN_RUNG_INFIDELITY} <= T < {MAX_RUNG_INFIDELITY}, got "
            f"{threshold!r}: at or above {MAX_RUNG_INFIDELITY} 'I <= T' is no longer "
            f"'|x| <= x_T', and below {MIN_RUNG_INFIDELITY} 1 - F cancels and x_T "
            f"loses precision"
        )
    return float(
        brentq(
            lambda x: float(infidelity_at(x)) - threshold,
            0.0,
            PEAK_X,
            xtol=np.finfo(float).tiny,
            rtol=4.0 * np.finfo(float).eps,
        )
    )


def infidelity_bounds(
    x: np.ndarray, sigma_x: np.ndarray, k: float
) -> tuple[np.ndarray, np.ndarray]:
    """(I(max(x - k s, 0)), I(x + k s)): the infidelity over a read's k-sigma interval in x.

    The read's error is symmetric in x, and I is monotone on [0, PEAK_X], so these are the
    extremes of I over the interval; the spread of I itself is lopsided, which is why the
    infidelity carries no sigma. For plots only, never stored by a step here: its caller is
    the thesis record figure's output builder (#21), which stores the bounds in its
    artifact so the renderer only draws them (AGENTS.md section 3).

    Raises ValueError on a bad k (`windows.check_k`), a negative x or sigma, and where
    x + k s passes PEAK_X, where I turns over and the endpoints stop being the extremes.
    A read with a non-finite sigma gets NaN bounds.

    Validity assumptions:

    - Assumption: the read's error is symmetric in x with sd sigma_x, Gaussian as `run`
      assumes, so x +/- k sigma_x is the read's k-sigma interval.
      Diagnostic: tellability's sigma rows (self-consistency, sigma bound) on the signed x.
      Consequence of violation: the interval in x has the wrong width and both bounds move
      with it; their order, and I(x) lying between them, hold for any error shape.
      Reference: issue #18; no equation locator. The other assumptions of `run` (the drive
      reading, Omega exact, one unit) carry over unchanged.
    """
    k = check_k(k)
    x = np.asarray(x, dtype=float)
    s = np.asarray(sigma_x, dtype=float)
    s = np.where(np.isfinite(s), s, np.nan)
    if np.any(x < 0.0):
        raise ValueError("x is |Delta| / Omega and cannot be negative")
    if np.any(s[np.isfinite(s)] < 0.0):
        raise ValueError("a reported sigma is negative")
    top = x + k * s
    if np.any(top > PEAK_X):
        raise ValueError(
            f"{int(np.sum(top > PEAK_X))} read(s) have x + k sigma above PEAK_X = "
            f"{PEAK_X}, where the infidelity turns over"
        )
    return infidelity_at(np.maximum(x - k * s, 0.0)), infidelity_at(top)


def run(inputs: FidelityInputs) -> FidelityResult:
    """The infidelity of each read, and its x = Delta / Omega with sigma, from typed inputs.

    The reference is the first read: Delta_i = f_i - f_0, x_i = Delta_i / Omega_i and
    sigma_x,i = sigma_f,i / Omega_i. The reads are the ones that survived the load, so f_0
    is the first read kept after the schema's row drops, the non-finite drops and the
    filter's backstop, not the file's first row. The chi-square check never moves it: a
    read at or above its threshold refuses the whole record (`transforms/filter.py`).

    The frame holds, per read: `delta_from_first_hz` (Delta), `frequency_error_hz`,
    `rabi_frequency_hz`, `x` (signed), `abs_x`, `sigma_x`, `fidelity` and `infidelity`
    (Eq. (9)). The infidelity carries no sigma: its spread is lopsided, so the carve
    classifies |x| with sigma_x instead (`recipes._fidelity_windows`).

    Raises ValueError on an empty record, on arrays of different lengths, on a Rabi
    frequency that is not finite and positive, and on a non-finite reference frequency.

    Validity assumptions:

    - Assumption: the read's error is Gaussian in x, with sd sigma_f / Omega.
      Diagnostic: tellability's sigma rows (self-consistency, sigma bound) on the signed x.
      Consequence of violation: a read's unsure label is too wide or too narrow, so the
      unsure count at each rung is off; the in/out side, decided by |x| alone, is not.
      Reference: issue #18; no equation locator.
    - Assumption: the drive reading. The X_pi/2 rotation is driven at f_0, the fitted
      frequency of the first read that survives the load's drops, a set number, so the
      detuning of read i is f_true,i - f_0 and its only error is sigma_f,i (read 0 is
      0 +/- sigma_f,0).
      Diagnostic: none from the data; it is a modelling choice.
      Consequence of violation: if the rotation instead tracked the true frequency at
      read 0, every Delta_i also carries read 0's error, so its sd is
      sqrt(sigma_f,i^2 + sigma_f,0^2) and every label is too confident by the factor
      sqrt(1 + (sigma_f,0 / sigma_f,i)^2): sqrt(2) when the two errors are equal, and with
      no upper bound when read 0's error is the larger.
      Reference: Capannelli et al., Supplementary notes, section V, Eq. (9); issue #18.
    - Assumption: Omega is exact. It is the declared per-qubit Rabi frequency
      (`transforms/declared_rabi.py`), never estimated here.
      Diagnostic: none; its error is a per-qubit systematic, and a rerun at Omega (1 +/- e)
      is backlog work.
      Consequence of violation: x scales by 1 / (1 + e) on every read, moving all reads
      together relative to every rung.
      Reference: issue #18 ("Not here").
    - Assumption: x is a ratio of two frequencies in one unit: Delta, sigma_f and Omega all
      in cyclic hertz.
      Diagnostic: `make_inputs_from_norm` and the declared join refuse a frequency declared
      in radians per second; a mislabelled value shows only on a real run, as the
      self-consistency and sigma-bound rows off by orders of magnitude.
      Consequence of violation: a factor 2 pi in x, so every label is wrong by that factor.
      Reference: Eq. (9) depends on Delta and Omega only through Delta / Omega.
    """
    t_rel_s = np.asarray(inputs.t_rel_s, dtype=float)
    frequency_hz = np.asarray(inputs.raw_frequency_hz, dtype=float)
    sigma_f_hz = np.asarray(inputs.frequency_error_hz, dtype=float)
    rabi_hz = np.asarray(inputs.rabi_hz, dtype=float)

    if len(t_rel_s) == 0:
        raise ValueError("Cannot run fidelity on empty dataset.")
    lengths = {len(a) for a in (t_rel_s, frequency_hz, sigma_f_hz, rabi_hz)}
    if len(lengths) != 1:
        raise ValueError(
            f"t_rel_s, raw_frequency_hz, frequency_error_hz and rabi_hz must have one "
            f"length; got {len(t_rel_s)}, {len(frequency_hz)}, {len(sigma_f_hz)} and "
            f"{len(rabi_hz)}"
        )
    if not np.all(np.isfinite(rabi_hz) & (rabi_hz > 0.0)):
        raise ValueError("the Rabi frequency must be finite and positive on every read")
    if not np.isfinite(frequency_hz[0]):
        raise ValueError(
            "the first read's frequency is not finite, so there is no reference to take "
            "the detuning from"
        )

    delta_hz = frequency_hz - float(frequency_hz[0])
    x = delta_hz / rabi_hz
    fidelity = gate_fidelity(
        _to_angular_frequency(delta_hz, inputs.use_angular_frequency),
        _to_angular_frequency(rabi_hz, inputs.use_angular_frequency),
    )
    infidelity = np.clip(1.0 - fidelity, 0.0, None)
    rabi_base_hz = float(np.median(rabi_hz))

    frame = pd.DataFrame(
        {
            "t_rel_s": t_rel_s,
            "delta_from_first_hz": delta_hz,
            "frequency_error_hz": sigma_f_hz,
            "rabi_frequency_hz": rabi_hz,
            "x": x,
            "abs_x": np.abs(x),
            "sigma_x": sigma_f_hz / rabi_hz,
            "fidelity": fidelity,
            "infidelity": infidelity,
            "rabi_base_hz": np.full_like(t_rel_s, rabi_base_hz, dtype=float),
        }
    )
    diag = {
        "mean_infidelity": float(np.mean(infidelity)),
        "max_infidelity": float(np.max(infidelity)),
        "rabi_base_hz": rabi_base_hz,
    }
    print(
        f"[fidelity] points={len(frame)} mean_infidelity={diag['mean_infidelity']:.6e} "
        f"max_infidelity={diag['max_infidelity']:.6e} rabi_base_hz={rabi_base_hz:.6g}",
        flush=True,
    )
    return FidelityResult(
        frame=frame,
        meta={
            "profile": inputs.profile,
            "rabi_base_hz": rabi_base_hz,
            "dataset_id": inputs.dataset_id,
            "rabi_source": inputs.rabi_source,
            "rabi_declared_key": inputs.rabi_declared_key,
            **inputs.load_counts,
        },
        diagnostics=diag,
    )


# ---------------------------------------------------------------------------
# Panel-data factory
#
# Kept out of the render layer: building panel data inside build_matplotlib happens at draw
# time, where no DAG node can supply the window and read tables. The adapter belongs beside
# its analyzer, as t2star's does.
# ---------------------------------------------------------------------------


# The nines ladder, in infidelity: 1e-2 (99%) to 1e-4 (99.99%), the thesis rungs, small
# values good. DECLARED, not chosen from the record: every rung is carved whether or not
# the record crosses it, so two records are carved against the same rungs and a rung the
# data never reaches gives one-sided windows rather than disappearing.
FIDELITY_LADDER: tuple[tuple[str, float, bool], ...] = tuple(
    ("99" + ("." + "9" * n if n > 0 else "") + "%", 10.0 ** (-(n + 2)), False)
    for n in range(3)
)


def detuning_ladder(
    ladder: Sequence[tuple[str, float, bool]],
) -> tuple[tuple[str, float, bool], ...]:
    """The same rungs on |x|: (label, x_threshold(T), False), labels unchanged.

    Raises ValueError on a big-values-good rung: infidelity is small-good, and so is |x|.
    """
    out = []
    for label, value, big_values_good in ladder:
        if big_values_good:
            raise ValueError(f"rung {label!r} is big-values-good; infidelity is not")
        out.append((label, x_threshold(value), False))
    return tuple(out)


# `detuning_ladder(FIDELITY_LADDER)` written out, so the carve's step kwargs (and with them
# the run identity) cannot move with the last bit of a root finder. Pinned to the helper by
# tests/test_infidelity_eq9.py.
DETUNING_LADDER: tuple[tuple[str, float, bool], ...] = (
    ("99%", 0.17365994622903097, False),
    ("99.9%", 0.05478655089906184, False),
    ("99.99%", 0.017320959849165627, False),
)

# The floor |x| is clipped at for a log axis, so the reference read (x = 0) can be drawn.
# I(ABS_X_FLOOR) is 1e-16 to leading order. Far below every rung, so the carve's
# classification is the same with or without it.
ABS_X_FLOOR = float(np.sqrt(3.0 * 1e-16))


def panel_series(result: FidelityResult) -> np.ndarray:
    """|x| clipped at ABS_X_FLOOR: the one series the carve, the panel and tellability's
    threshold rows read, so they classify and draw the same values."""
    return np.clip(result.frame["abs_x"].to_numpy(dtype=float), ABS_X_FLOOR, None)


# The caption's note on the damage axis, which integrates |x| above each rung.
DAMAGE_NOTE = (
    "cumulative damage: detuning excess above each rung, in units of Ω, times hours"
)


def make_panel_data(
    result: FidelityResult,
    windows: pd.DataFrame,
    reads: pd.DataFrame,
    gap_spans_s: list[tuple[float, float]],
    kaplan_meier: KaplanMeierSet,
    thresholds: Sequence[tuple[str, float, bool]],
    shape_min_reads: int = 5,
    xi_seed: int = 0,
    k: float = 1.0,
    use_uncertainty: bool = False,
):
    """Convert FidelityResult + the window tables to WithinCalibrationPanelData.

    The panel draws what the carve classified: |x| = |Delta| / Omega on a log axis
    (`panel_series`), sigma_x as its error bars, and `thresholds`, the ladder on |x|
    (`DETUNING_LADDER`), whose labels name the fidelity each rung keeps. The relative
    error sigma / value is computed but not drawn: x has no natural zero, and the
    reference read is 0 by construction. The caption's dataset is
    `result.meta["dataset_id"]`, set by `run` from the norm.
    """
    from quebra.analyzers.within_calibration_compute import (
        build_within_calibration_panel_data,
    )
    from quebra.plots.theme import qubit_color

    frame = result.frame
    t_h = frame["t_rel_s"].to_numpy(dtype=float) / 3600.0

    dataset_id = str(result.meta["dataset_id"])
    meta: dict[str, object] = {"dataset": dataset_id, "caption_note": DAMAGE_NOTE}
    rabi_base_hz = result.meta.get("rabi_base_hz")
    if rabi_base_hz is not None:
        meta["rabi_base_hz"] = f"{rabi_base_hz:.4g} Hz"
    profile = result.meta.get("profile")
    if profile:
        meta["profile"] = str(profile)

    if rabi_base_hz is not None:
        rabi_mhz = float(rabi_base_hz) / 1e6
        primary_label = f"|Δ|/Ω  ({rabi_mhz:.4g} MHz Rabi)"
    else:
        primary_label = "|Δ|/Ω"

    return build_within_calibration_panel_data(
        t_h=t_h,
        primary_series=panel_series(result),
        primary_label=primary_label,
        thresholds=list(thresholds),
        meta=meta,
        windows=windows,
        reads=reads,
        gap_spans_s=gap_spans_s,
        kaplan_meier=kaplan_meier,
        primary_sigma=frame["sigma_x"].to_numpy(dtype=float),
        shape_min_reads=shape_min_reads,
        xi_seed=xi_seed,
        k=k,
        use_uncertainty=use_uncertainty,
        use_log_scale=True,
        include_relative_error=False,
        color=qubit_color(dataset_id=dataset_id),
    )
