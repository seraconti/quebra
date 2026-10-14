from __future__ import annotations

from dataclasses import dataclass, field
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from quebra.analyzers.kaplan_meier import KaplanMeierSet


@dataclass(slots=True)
class FidelityResult:
    frame: pd.DataFrame
    meta: dict[str, object]
    diagnostics: dict[str, object] = field(default_factory=dict)


@dataclass(slots=True)
class FidelityInputs:
    t_rel_s: np.ndarray
    delta_hz: np.ndarray
    rabi_hz: np.ndarray
    profile: str
    raw_frequency_hz: np.ndarray | None = None
    use_angular_frequency: bool = False
    dataset_id: str = ""
    rabi_source: str = ""
    rabi_declared_key: str = ""


def make_inputs_from_norm(
    norm: Mapping[str, object], config: Mapping[str, object]
) -> FidelityInputs:
    """Extract FidelityInputs from a normalized mapping and config dict."""
    if "t_rel_s" not in norm:
        raise KeyError(
            "Fidelity requires 't_rel_s' (relative seconds) in the normalized mapping."
        )
    if "rabi_hz" not in norm:
        raise KeyError(
            "Fidelity requires 'rabi_hz' in the normalized mapping. Declare the record's "
            "Rabi frequency (configure_ramsey_job(..., rabi_table=Dataset(...))) "
            "or disable fidelity for this job."
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
    raw_frequency_hz = (
        np.asarray(norm["raw_frequency_hz"], dtype=float)
        if "raw_frequency_hz" in norm
        else None
    )
    return FidelityInputs(
        t_rel_s=np.asarray(norm["t_rel_s"], dtype=float),
        delta_hz=np.asarray(norm["delta_hz"], dtype=float),
        rabi_hz=np.asarray(norm["rabi_hz"], dtype=float),
        profile=profile,
        raw_frequency_hz=raw_frequency_hz,
        use_angular_frequency=bool(fidelity_cfg.get("use_angular_frequency", False)),
        dataset_id=dataset_id,
        rabi_source=rabi_source,
        rabi_declared_key=str(meta.get("rabi_declared_key", "")),
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


def run(inputs: FidelityInputs) -> FidelityResult:
    """Compute gate fidelity from typed FidelityInputs.

    Frequencies are in hertz unless inputs.use_angular_frequency is True.
    """
    t_rel_s = inputs.t_rel_s
    delta_hz = inputs.delta_hz
    rabi_hz = inputs.rabi_hz

    if len(t_rel_s) == 0:
        raise ValueError("Cannot run fidelity on empty dataset.")

    rabi_drive = _to_angular_frequency(rabi_hz, inputs.use_angular_frequency)
    rabi_base_hz = float(np.median(rabi_hz))
    result_meta = {
        "profile": inputs.profile,
        "rabi_base_hz": rabi_base_hz,
        "dataset_id": inputs.dataset_id,
        "rabi_source": inputs.rabi_source,
        "rabi_declared_key": inputs.rabi_declared_key,
    }

    if inputs.profile == "longrun":
        if inputs.raw_frequency_hz is None:
            raise KeyError(
                "Longrun fidelity requires 'raw_frequency_hz' in FidelityInputs."
            )
        raw_frequency_hz = inputs.raw_frequency_hz
        if len(raw_frequency_hz) != len(t_rel_s):
            raise ValueError(
                "raw_frequency_hz length must match t_rel_s for longrun fidelity."
            )

        delta_f0_hz = raw_frequency_hz - float(raw_frequency_hz[0])
        delta_fmean_hz = raw_frequency_hz - float(np.mean(raw_frequency_hz))
        delta_f0 = _to_angular_frequency(delta_f0_hz, inputs.use_angular_frequency)
        delta_fmean = _to_angular_frequency(
            delta_fmean_hz, inputs.use_angular_frequency
        )

        fidelity_f0 = gate_fidelity(delta_f0, rabi_drive)
        fidelity_fmean = gate_fidelity(delta_fmean, rabi_drive)
        infidelity_f0 = np.clip(1.0 - fidelity_f0, 0.0, None)
        infidelity_fmean = np.clip(1.0 - fidelity_fmean, 0.0, None)

        frame = pd.DataFrame(
            {
                "t_rel_s": t_rel_s,
                "delta_hz": delta_fmean_hz,
                "delta_hz_f_minus_f0": delta_f0_hz,
                "delta_hz_f_minus_fmean": delta_fmean_hz,
                "delta_hz_f0": delta_f0_hz,
                "delta_hz_fmean": delta_fmean_hz,
                "rabi_frequency_hz": rabi_hz,
                "fidelity": fidelity_fmean,
                "fidelity_f_minus_f0": fidelity_f0,
                "fidelity_f_minus_fmean": fidelity_fmean,
                "fidelity_f0": fidelity_f0,
                "fidelity_fmean": fidelity_fmean,
                "infidelity": infidelity_fmean,
                "infidelity_f_minus_f0": infidelity_f0,
                "infidelity_f_minus_fmean": infidelity_fmean,
                "infidelity_f0": infidelity_f0,
                "infidelity_fmean": infidelity_fmean,
                "rabi_base_hz": np.full_like(t_rel_s, rabi_base_hz, dtype=float),
            }
        )
        diag = {
            "mean_infidelity_fmean": float(np.mean(infidelity_fmean)),
            "mean_infidelity_f0": float(np.mean(infidelity_f0)),
            "rabi_base_hz": rabi_base_hz,
        }
        print(
            f"[fidelity] points={len(frame)} mean_infidelity_fmean={diag['mean_infidelity_fmean']:.6e} mean_infidelity_f0={diag['mean_infidelity_f0']:.6e} rabi_base_hz={rabi_base_hz:.6g}",
            flush=True,
        )
        return FidelityResult(
            frame=frame,
            meta=dict(result_meta),
            diagnostics=diag,
        )

    delta = _to_angular_frequency(delta_hz, inputs.use_angular_frequency)
    fidelity = gate_fidelity(delta, rabi_drive)
    infidelity = np.clip(1.0 - fidelity, 0.0, None)

    frame = pd.DataFrame(
        {
            "t_rel_s": t_rel_s,
            "delta_hz": delta_hz,
            "rabi_frequency_hz": rabi_hz,
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
        f"[fidelity] points={len(frame)} mean_infidelity={diag['mean_infidelity']:.6e} max_infidelity={diag['max_infidelity']:.6e} rabi_base_hz={rabi_base_hz:.6g}",
        flush=True,
    )
    return FidelityResult(
        frame=frame,
        meta=dict(result_meta),
        diagnostics=diag,
    )


# ---------------------------------------------------------------------------
# Panel-data factory
#
# Kept out of the render layer: building panel data inside build_matplotlib happens at draw
# time, where no DAG node can supply the window and read tables. The adapter belongs beside
# its analyzer, as t2star's does.
# ---------------------------------------------------------------------------


# The nines ladder, in infidelity: 1e-2 (99%) to 1e-9 (99.9999999%), small values good.
# DECLARED, not chosen from the record: every rung is carved whether or not the record
# crosses it, so two records are carved against the same rungs and a rung the data never
# reaches gives one-sided windows rather than disappearing. The recipe passes this one tuple
# to both the carve and the panel-data step.
FIDELITY_LADDER: tuple[tuple[str, float, bool], ...] = tuple(
    ("99" + ("." + "9" * n if n > 0 else "") + "%", 10.0 ** (-(n + 2)), False)
    for n in range(8)
)


# The floor infidelity is clipped at, so that its logarithm stays finite. One value for the
# panel, the carve and tellability's clipped fraction.
INFIDELITY_FLOOR = 1e-16


def panel_series(result: FidelityResult) -> np.ndarray:
    """The infidelity series the panel and the carve both use (clipped identically)."""
    return np.clip(
        result.frame["infidelity"].to_numpy(dtype=float), INFIDELITY_FLOOR, None
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

    The caption's dataset is `result.meta["dataset_id"]`, set by `run` from the norm.
    """
    from quebra.analyzers.within_calibration_compute import (
        build_within_calibration_panel_data,
    )
    from quebra.plots.theme import qubit_color

    frame = result.frame
    t_h = frame["t_rel_s"].to_numpy(dtype=float) / 3600.0
    infidelity = panel_series(result)

    traces: list[tuple[str, np.ndarray]] | None = None
    if "infidelity_f0" in frame.columns:
        traces = [
            ("f−fₘₑₐₙ", infidelity),
            (
                "f−f₀",
                np.clip(
                    frame["infidelity_f0"].to_numpy(dtype=float), INFIDELITY_FLOOR, None
                ),
            ),
        ]

    dataset_id = str(result.meta["dataset_id"])
    meta: dict[str, object] = {"dataset": dataset_id}
    rabi_base_hz = result.meta.get("rabi_base_hz")
    if rabi_base_hz is not None:
        meta["rabi_base_hz"] = f"{rabi_base_hz:.4g} Hz"
    profile = result.meta.get("profile")
    if profile:
        meta["profile"] = str(profile)

    if rabi_base_hz is not None:
        rabi_mhz = float(rabi_base_hz) / 1e6
        primary_label = f"Infidelity  ({rabi_mhz:.4g} MHz Rabi)"
    else:
        primary_label = "Infidelity"

    return build_within_calibration_panel_data(
        t_h=t_h,
        primary_series=infidelity,
        primary_label=primary_label,
        thresholds=list(thresholds),
        meta=meta,
        windows=windows,
        reads=reads,
        gap_spans_s=gap_spans_s,
        kaplan_meier=kaplan_meier,
        shape_min_reads=shape_min_reads,
        xi_seed=xi_seed,
        k=k,
        use_uncertainty=use_uncertainty,
        traces=traces,
        use_log_scale=True,
        color=qubit_color(dataset_id=dataset_id),
    )
