from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np

Norm = dict[str, object]

# The filter config keys `run` reads; any other key is refused.
FILTER_KEYS = frozenset({"apply_chi_squared", "chi_squared_threshold"})


@dataclass(slots=True)
class FilterResult:
    stages: dict[str, Norm]
    stage_order: list[str]
    final_stage: str
    counts: list[dict[str, int | str]]
    meta: dict[str, object]
    diagnostics: dict[str, object] = field(default_factory=dict)

    @property
    def final_norm(self) -> Norm:
        return self.stages[self.final_stage]


def _copy_norm(norm: Mapping[str, object]) -> Norm:
    copied = dict(norm)
    meta = copied.get("meta")
    if isinstance(meta, dict):
        copied["meta"] = dict(meta)
    return copied


def _subset_norm(norm: Mapping[str, object], mask: np.ndarray) -> Norm:
    """Apply a row mask to every row-aligned column, or raise naming the one that cannot be.

    Every array a Norm carries is one value per read, so a length disagreeing with the mask
    is a defect rather than a column with nothing to filter. Left unfiltered it would be
    index-misaligned against `t_rel_s` for every step downstream, and no consumer raises on
    it. Scalars pass through; an array that cannot be boxed at all says so.
    """
    out = _copy_norm(norm)
    mask = np.asarray(mask, dtype=bool)
    for key, value in list(out.items()):
        if key == "meta":
            continue
        try:
            arr = np.asarray(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"filter cannot mask '{key}': it is neither a scalar nor an array "
                f"({type(value).__name__}). A Norm value has to be one or the other, or it "
                f"cannot be kept aligned with 't_rel_s'."
            ) from exc
        if arr.ndim == 0:
            continue
        if len(arr) != len(mask):
            raise ValueError(
                f"filter cannot mask '{key}': length {len(arr)} against a mask of "
                f"{len(mask)}. Every array in a Norm is one value per read, so leaving it "
                f"unfiltered would misalign it against 't_rel_s'."
            )
        out[key] = arr[mask]

    if "t_rel_s" not in out:
        raise KeyError(
            "filter requires 't_rel_s' (relative seconds) in normalized mapping"
        )
    t_rel_s = np.asarray(out["t_rel_s"], dtype=float)
    if len(t_rel_s) > 0:
        t_rel_s = t_rel_s - float(t_rel_s[0])
    out["t_rel_s"] = t_rel_s
    if "meta" not in out or not isinstance(out["meta"], dict):
        out["meta"] = {}
    out["meta"]["n_points"] = int(len(t_rel_s))
    out["meta"]["duration_h"] = (
        float((np.max(t_rel_s) - np.min(t_rel_s)) / 3600.0) if len(t_rel_s) > 1 else 0.0
    )
    return out


def run(
    norm: Mapping[str, object], config: Mapping[str, object], output_dir: object
) -> FilterResult:
    """Filter a normalized Ramsey dataset.

    Two stages: a chi-square GUARD that refuses a record containing a read at or above the
    published fit-quality cut, and a stage that drops a read with no finite detuning or time.
    No read is dropped for its value: an outlier is an observation. Time values are in seconds
    and frequency values are in hertz.
    """
    del output_dir

    if not isinstance(norm, Mapping):
        raise TypeError("filter.run expects normalized dataset mapping input.")

    dataset_id = norm["meta"]["dataset_id"]
    if "t_rel_s" not in norm:
        raise KeyError(
            "filter.run requires 't_rel_s' (relative seconds) in normalized mapping"
        )
    t_rel_s = np.asarray(norm["t_rel_s"], dtype=float)
    if len(t_rel_s) == 0:
        raise ValueError("Cannot run filter on empty dataset.")

    print(f"[filter] dataset={dataset_id} start_points={len(t_rel_s)}", flush=True)

    current = _copy_norm(norm)
    stages: dict[str, Norm] = {"raw": _copy_norm(norm)}
    stage_order = ["raw"]
    rows: list[dict[str, int | str]] = [{"stage": "raw", "n_points": int(len(t_rel_s))}]

    filter_cfg = (
        config.get("filter", {})
        if isinstance(config.get("filter", {}), Mapping)
        else {}
    )
    # A key this filter does not read is refused, not ignored: a caller passing a cut that no
    # longer exists would otherwise get none and be told nothing.
    unknown = sorted(set(filter_cfg) - FILTER_KEYS)
    if unknown:
        raise KeyError(
            f"unknown filter config key(s) {unknown}; the filter reads only "
            f"{sorted(FILTER_KEYS)}"
        )
    apply_chi = bool(filter_cfg.get("apply_chi_squared", True))

    if apply_chi:
        if "chi_squared" not in current:
            raise KeyError(
                "the chi-square guard requires a 'chi_squared' key in the normalized mapping."
            )
        if "chi_squared_threshold" not in filter_cfg:
            raise KeyError(
                "the chi-square guard requires config.filter.chi_squared_threshold."
            )
        chi_threshold = float(filter_cfg["chi_squared_threshold"])
        # A GUARD, not a cut. Capannelli et al., "Tracking spin qubit frequency variations
        # over 912 days", Supplementary notes, section III.C, p. 9, filter out every drift
        # dataset whose normalised chi-squared score is 0.8 or more, and the released records
        # hold the points left after that filter (Supplementary Figs. 10 and 11, pp. 10-11,
        # fitted minus filtered points, for the linear and the 2x2 device). A read at or above the threshold therefore means the record is not
        # the published one; dropping it here would hide that, so the record is refused. A
        # missing chi-square is refused too, though a Norm built through `track912Schema` never
        # has one: its validation drops such rows on load.
        chi = np.asarray(current["chi_squared"], dtype=float)
        failing = ~(chi < chi_threshold)
        if np.any(failing):
            worst = (
                np.nanmax(chi[failing])
                if np.any(np.isfinite(chi[failing]))
                else chi[failing][0]
            )
            raise ValueError(
                f"[{dataset_id}] {int(failing.sum())} read(s) have a normalised chi-square "
                f">= {chi_threshold} or none at all (largest {worst!r}). The records are "
                f"expected already filtered at this threshold; filter them upstream, or set "
                f"filter.apply_chi_squared=False deliberately (that choice is not a step "
                f"kwarg, so it does not appear on the provenance label)."
            )
        rows.append({"stage": "chi_squared", "n_points": int(len(current["t_rel_s"]))})
        stages["chi_squared"] = _copy_norm(current)
        stage_order.append("chi_squared")
        print(
            f"[filter] chi-square guard: all {len(current['t_rel_s'])} reads below "
            f"{chi_threshold}",
            flush=True,
        )

    # A BACKSTOP. A read with no finite detuning or time carries no observation, and is
    # dropped and counted here. On the 6D2S and 2x2 records it drops nothing: on load the
    # schemas already remove rows with a missing (NaN) T2*, phase or chi-square
    # (`track912Schema.validate`) and rows with a non-finite time or frequency
    # (`RamseySeriesSchema.to_norm`), and those removals are not counted here. It also masks every row-aligned column, which is what checks their
    # alignment.
    if "delta_hz" not in current:
        raise KeyError("filter.run requires 'delta_hz' (hertz) in normalized mapping")
    delta_hz = np.asarray(current["delta_hz"], dtype=float)
    finite = np.isfinite(delta_hz) & np.isfinite(
        np.asarray(current["t_rel_s"], dtype=float)
    )
    if not np.any(finite):
        raise ValueError(
            f"filter has no read with both a finite 'delta_hz' and a finite 't_rel_s' "
            f"({delta_hz.size} read(s): {int((~np.isfinite(delta_hz)).sum())} non-finite "
            f"delta_hz, "
            f"{int((~np.isfinite(np.asarray(current['t_rel_s'], dtype=float))).sum())} "
            f"non-finite t_rel_s)."
        )
    current = _subset_norm(current, finite)
    rows.append(
        {
            "stage": "non_finite",
            "n_points": int(len(current["t_rel_s"])),
            "n_dropped": int((~finite).sum()),
        }
    )
    stages["non_finite"] = _copy_norm(current)
    stage_order.append("non_finite")
    print(
        f"[filter] non_finite: kept {len(current['t_rel_s'])}, dropped {int((~finite).sum())}",
        flush=True,
    )

    diagnostics: dict[str, object] = {"stage_counts": rows}
    return FilterResult(
        stages=stages,
        stage_order=stage_order,
        final_stage=stage_order[-1],
        counts=rows,
        meta=dict(current["meta"]),
        diagnostics=diagnostics,
    )
