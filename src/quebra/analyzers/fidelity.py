from __future__ import annotations

from dataclasses import dataclass, field
from collections.abc import Mapping

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
            "Fidelity requires 'rabi_hz' in the normalized mapping. "
            "Add it via lookup_prior(..., fields=['Rabi_frequency'], aliases={'Rabi_frequency': 'rabi_hz'}) "
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
    )


def _to_angular_frequency(
    freq_hz: np.ndarray, use_angular_frequency: bool
) -> np.ndarray:
    frequency_hz = np.asarray(freq_hz, dtype=float)
    if use_angular_frequency:
        return 2.0 * np.pi * frequency_hz
    return frequency_hz


def gate_fidelity(delta: np.ndarray, rabi_drive: np.ndarray) -> np.ndarray:
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
            meta={
                "profile": inputs.profile,
                "rabi_base_hz": rabi_base_hz,
                "dataset_id": inputs.dataset_id,
            },
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
        meta={
            "profile": inputs.profile,
            "rabi_base_hz": rabi_base_hz,
            "dataset_id": inputs.dataset_id,
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


def panel_thresholds(infidelity: np.ndarray) -> list[tuple[str, float, bool]]:
    """Nines thresholds (infidelity units, big_values_good=False) that the data crosses.

    Data-derived, so the windows step and the panel-data step must both call this on the
    SAME series to carve and to render the same ladder. It is pure, so they agree.
    """
    inf = np.asarray(infidelity, dtype=float)
    inf_finite = inf[np.isfinite(inf) & (inf > 0)]
    if len(inf_finite) == 0:
        return []
    inf_min, inf_max = float(np.min(inf_finite)), float(np.max(inf_finite))
    thresholds: list[tuple[str, float, bool]] = []
    for n in range(0, 8):
        # n=0 → infidelity<0.01 = 99% fidelity; n=1 → infidelity<0.001 = 99.9%; ...
        inf_thr = 10.0 ** (-(n + 2))
        label = "99" + ("." + "9" * n if n > 0 else "") + "%"
        if inf_min < inf_thr < inf_max:
            # big_values_good=False: above this infidelity threshold = out-of-spec
            thresholds.append((label, inf_thr, False))
    return thresholds


def panel_series(result: FidelityResult) -> np.ndarray:
    """The infidelity series the panel and the carve both use (clipped identically)."""
    return np.clip(result.frame["infidelity"].to_numpy(dtype=float), 1e-16, None)


def make_panel_data(
    result: FidelityResult,
    windows: pd.DataFrame,
    reads: pd.DataFrame,
    gap_spans_s: list[tuple[float, float]],
    kaplan_meier: KaplanMeierSet,
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
                np.clip(frame["infidelity_f0"].to_numpy(dtype=float), 1e-16, None),
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
        thresholds=panel_thresholds(infidelity),
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
