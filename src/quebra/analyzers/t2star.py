"""T2* coherence-time analyzer.

T2STAR_DEFAULT_LADDER - analyzer-supplied starting point for 912-day silicon
quantum dot data (T2* range 1–10 µs).  Production jobs MUST supply their own
thresholds via the `thresholds=` argument to make_panel_data:

    - Devices with T2* ranges outside 1–10 µs (e.g. superconducting qubits
      with T2* > 100 µs) MUST supply their own thresholds; the default ladder
      will not be informative.
    - Values near 1 µs may indicate Ramsey fit failure rather than genuinely
      short T2*.  Reader caveat, not auto-filtered.

Values in the ladder are in SI seconds; make_panel_data converts to µs for
display.  Labels are in µs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from collections.abc import Mapping

import numpy as np
import pandas as pd

from quebra.analyzers.kaplan_meier import KaplanMeierSet
from quebra.analyzers.within_calibration_compute import (
    build_within_calibration_panel_data,
)
from quebra.analyzers.within_calibration_data import WithinCalibrationPanelData


@dataclass(slots=True)
class T2StarResult:
    frame: pd.DataFrame  # columns: t_rel_s, t2star_s, t2star_error_s (optional)
    meta: dict[str, object]
    diagnostics: dict[str, object] = field(default_factory=dict)


@dataclass(slots=True)
class T2StarInputs:
    t_rel_s: np.ndarray
    t2star_s: np.ndarray
    t2star_error_s: np.ndarray | None = None
    dataset_id: str = ""


def make_inputs_from_norm(
    norm: Mapping[str, object], config: Mapping[str, object] | None = None
) -> T2StarInputs:
    """Extract T2StarInputs from a normalized mapping."""
    if "t_rel_s" not in norm:
        raise KeyError("T2* requires 't_rel_s' in the normalized mapping.")
    if "T2star_s" not in norm:
        raise KeyError(
            "T2* requires 'T2star_s' in the normalized mapping. "
            "The dataset must contain a 'T2star' column."
        )
    meta = norm.get("meta", {}) if isinstance(norm.get("meta", {}), Mapping) else {}
    dataset_id = str(meta.get("dataset_id") or "").strip()
    if not dataset_id:
        raise KeyError(
            "T2* requires a non-empty meta['dataset_id'] in the normalized mapping: the "
            "panel caption names the record it draws."
        )
    t2star_error_s = (
        np.asarray(norm["T2star_error_s"], dtype=float)
        if "T2star_error_s" in norm
        else None
    )
    return T2StarInputs(
        t_rel_s=np.asarray(norm["t_rel_s"], dtype=float),
        t2star_s=np.asarray(norm["T2star_s"], dtype=float),
        t2star_error_s=t2star_error_s,
        dataset_id=dataset_id,
    )


def run(inputs: T2StarInputs) -> T2StarResult:
    """Extract and summarise T2* time series from pre-measured inputs."""
    t_rel_s = inputs.t_rel_s
    t2star_s = inputs.t2star_s
    if len(t_rel_s) == 0:
        raise ValueError("Cannot run T2* analysis on empty dataset.")
    if len(t_rel_s) != len(t2star_s):
        raise ValueError("t_rel_s and t2star_s must have the same length.")

    valid = np.isfinite(t2star_s)
    n_valid = int(np.sum(valid))
    n_total = len(t_rel_s)

    cols: dict[str, np.ndarray] = {"t_rel_s": t_rel_s, "t2star_s": t2star_s}
    if inputs.t2star_error_s is not None:
        cols["t2star_error_s"] = inputs.t2star_error_s

    t2star_valid = t2star_s[valid]
    # NaN, not 0.0, where undefined: a sample std needs two reads, a mean one.
    mean_us = float(np.mean(t2star_valid)) * 1e6 if n_valid else float("nan")
    std_us = float(np.std(t2star_valid, ddof=1)) * 1e6 if n_valid > 1 else float("nan")

    diag: dict[str, object] = {
        "n_valid": n_valid,
        "n_total": n_total,
        "n_nan": n_total - n_valid,
        "mean_us": mean_us,
        "std_us": std_us,
    }
    print(
        f"[t2star] points={n_valid}/{n_total} mean={mean_us:.3f}µs std={std_us:.3f}µs",
        flush=True,
    )
    return T2StarResult(
        frame=pd.DataFrame(cols),
        meta={"dataset_id": inputs.dataset_id, "mean_us": mean_us, "std_us": std_us},
        diagnostics=diag,
    )


# ---------------------------------------------------------------------------
# Analyzer-supplied default ladder (912-day silicon quantum dot context)
# Values in SI seconds; labels in µs.  Jobs should declare their own
# thresholds - see module docstring.
# ---------------------------------------------------------------------------

T2STAR_DEFAULT_LADDER: list[tuple[str, float, bool]] = [
    ("1 µs", 1e-6, True),
    ("2 µs", 2e-6, True),
    ("3 µs", 3e-6, True),
    ("4 µs", 4e-6, True),
    ("5 µs", 5e-6, True),
    ("6 µs", 6e-6, True),
    ("7 µs", 7e-6, True),
    ("8 µs", 8e-6, True),
    ("9 µs", 9e-6, True),
    ("10 µs", 10e-6, True),
]


# ---------------------------------------------------------------------------
# Panel-data factory
# ---------------------------------------------------------------------------


def make_panel_data(
    result: T2StarResult,
    windows: pd.DataFrame,
    reads: pd.DataFrame,
    gap_spans_s: list[tuple[float, float]],
    kaplan_meier: KaplanMeierSet,
    shape_min_reads: int = 5,
    xi_seed: int = 0,
    k: float = 1.0,
    use_uncertainty: bool = False,
    thresholds: list[tuple[str, float, bool]] | None = None,
    primary_label: str | None = None,
) -> WithinCalibrationPanelData:
    """Convert T2StarResult to WithinCalibrationPanelData for WithinCalibrationPanel.

    `thresholds`: job-supplied list of (label, value_s, big_values_good) triples
    where values are in SI seconds.  If None, falls back to T2STAR_DEFAULT_LADDER
    (appropriate for the 912-day context; not universal - see module docstring).
    If [], no derived threshold views render.
    """
    frame = result.frame
    t_h = frame["t_rel_s"].to_numpy(dtype=float) / 3600.0
    t2star_us = frame["t2star_s"].to_numpy(dtype=float) * 1e6

    resolved = thresholds if thresholds is not None else T2STAR_DEFAULT_LADDER
    # Convert threshold values from SI seconds to µs to match primary_series units.
    panel_thresholds = [(lbl, val * 1e6, bvg) for lbl, val, bvg in resolved]
    # The carve ran in SI seconds on the SI ladder; the panel plots µs. The panel takes
    # every in/out classification from the carve's tables and never reclassifies the
    # rescaled values, so a read that sits on a threshold cannot change side here.
    sigma_us = (
        frame["t2star_error_s"].to_numpy(dtype=float) * 1e6
        if "t2star_error_s" in frame.columns
        else None
    )

    label = primary_label if primary_label is not None else "T2* (µs)"
    # The caption's mean and std are the signal band's (`value_mean`, `value_std`).
    meta: dict[str, object] = {"dataset": str(result.meta["dataset_id"])}

    return build_within_calibration_panel_data(
        t_h=t_h,
        primary_series=t2star_us,
        primary_label=label,
        thresholds=panel_thresholds,
        meta=meta,
        windows=windows,
        reads=reads,
        gap_spans_s=gap_spans_s,
        kaplan_meier=kaplan_meier,
        shape_min_reads=shape_min_reads,
        xi_seed=xi_seed,
        k=k,
        use_uncertainty=use_uncertainty,
        primary_sigma=sigma_us,
        use_log_scale=False,
    )
