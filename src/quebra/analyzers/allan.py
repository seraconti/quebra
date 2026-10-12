from __future__ import annotations

from dataclasses import dataclass, field
from collections.abc import Mapping

import allantools
import numpy as np
import pandas as pd


@dataclass(slots=True)
class AllanResult:
    modes: dict[str, pd.DataFrame]
    meta: dict[str, object]
    fractional_adev: np.ndarray | None = None
    carrier_hz: float | None = None
    diagnostics: dict[str, object] = field(default_factory=dict)


def infer_sample_period_s(timestamps: np.ndarray) -> float:
    dt_s = np.diff(np.asarray(timestamps, dtype=float))
    dt_s = dt_s[dt_s > 0]
    if len(dt_s) == 0:
        raise ValueError("Cannot infer sample period from timestamp data.")
    return float(np.median(dt_s))


def _allan_deviation_overlapping(
    values: np.ndarray, dt_s: float, min_points: int, taus_mode: object
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    signal = np.asarray(values, dtype=float)
    signal = signal[np.isfinite(signal)]
    if len(signal) < min_points:
        raise ValueError(
            f"Not enough points for Allan deviation: got {len(signal)}, need at least {min_points}."
        )

    fs_hz = 1.0 / float(dt_s)
    tau_s, adev, adev_err, n_pairs = allantools.oadev(
        signal,
        rate=fs_hz,
        data_type="freq",
        taus=taus_mode,
    )
    if len(tau_s) == 0:
        raise ValueError("No valid tau points were produced by allantools.oadev.")
    return (
        np.asarray(tau_s, dtype=float),
        np.asarray(adev, dtype=float),
        np.asarray(adev_err, dtype=float),
        np.asarray(n_pairs, dtype=int),
    )


def _transform_signal(freq_hz: np.ndarray, mode: str) -> np.ndarray:
    """The series a mode computes on. One mode exists: `raw`, the detuning as read."""
    if mode == "raw":
        return np.asarray(freq_hz, dtype=float)
    raise ValueError(f"Unknown Allan mode '{mode}'; the only mode is 'raw'.")


def run(
    norm: Mapping[str, object],
    config: Mapping[str, object],
    fractional: bool = False,
    carrier_col: str = "qubit_frequency_hz",
) -> AllanResult:
    """Compute Allan deviation for a normalized Ramsey dataset.

    Inputs are expected to use seconds for time and hertz for frequency.
    """
    if not isinstance(norm, Mapping):
        raise TypeError("allan.run expects normalized dataset mapping input.")

    if "t_rel_s" not in norm:
        raise KeyError(
            "allan.run requires 't_rel_s' (relative seconds) in normalized mapping"
        )
    t_rel_s = np.asarray(norm["t_rel_s"], dtype=float)
    delta_hz = np.asarray(norm["delta_hz"], dtype=float)
    dt_s = infer_sample_period_s(t_rel_s)
    fs_hz = 1.0 / float(dt_s)
    meta_in = (
        dict(norm.get("meta", {})) if isinstance(norm.get("meta", {}), Mapping) else {}
    )

    allan_cfg = config.get("allan", {})
    if not isinstance(allan_cfg, Mapping):
        raise TypeError("config['allan'] must be a mapping")
    modes = list(allan_cfg.get("modes", []))
    if len(modes) == 0:
        raise ValueError("config['allan']['modes'] must be a non-empty list")
    min_points = int(allan_cfg.get("min_points_for_allan", 8))
    taus_mode = allan_cfg.get("taus_mode", "all")

    summary: dict[str, pd.DataFrame] = {}
    for mode in modes:
        signal = _transform_signal(delta_hz, mode=str(mode))
        tau_s, adev, adev_err, n_pairs = _allan_deviation_overlapping(
            signal, dt_s, min_points, taus_mode
        )
        summary[str(mode)] = pd.DataFrame(
            {
                "tau_s": tau_s,
                "adev": adev,
                "adev_err": adev_err,
                "n_pairs": n_pairs,
                "dt_s": np.full_like(tau_s, dt_s, dtype=float),
            }
        )
        print(
            f"[allan] mode={mode} tau_points={len(tau_s)} avg_f_acq={fs_hz:.6g}Hz",
            flush=True,
        )

    fractional_adev: np.ndarray | None = None
    carrier_hz: float | None = None
    if fractional:
        # carrier column must exist in the normalized mapping
        if carrier_col not in norm:
            raise ValueError(
                f"fractional=True requires column '{carrier_col}' in the DataFrame. Found: {list(norm.keys())}. "
                f"Pass carrier_col= naming a column present in the norm, or run with "
                f"fractional=False."
            )
        carrier_vals = np.asarray(norm[carrier_col], dtype=float)
        # mean carrier (Hz)
        carrier_hz = float(np.mean(carrier_vals))
        # compute fractional ADEV using the same adev values as the first mode (if present)
        # choose first available mode's adev as base reference for fractional conversion
        first_mode_vals = None
        for mode_name in modes:
            dfm = summary[str(mode_name)]
            first_mode_vals = np.asarray(dfm["adev"], dtype=float)
            break
        if first_mode_vals is not None and carrier_hz != 0.0:
            fractional_adev = first_mode_vals / carrier_hz

    return AllanResult(
        modes=summary,
        meta={
            "sample_period_s": float(dt_s),
            "acquisition_frequency_hz": float(fs_hz),
            "dataset_id": meta_in.get("dataset_id"),
            "run_name": meta_in.get("run_name"),
            "qubit": meta_in.get("qubit"),
        },
        fractional_adev=fractional_adev,
        carrier_hz=carrier_hz,
    )
