"""Render-local helpers for WithinCalibrationPanel: functions of axes and theme only.

Nothing here derives a SCIENTIFIC quantity or writes into the artifact. These are
view-only functions of already-computed arrays: decade guides, bin geometry, the
observed-stretch slicing, and the claim-discipline annotation. Every verdict the
annotation states (whether a peak is symmetric, whether pooling is too thin) is a flag
computed in the analyzer band and handed here already decided. They do arithmetic - `adaptive_ylim` takes
`observed_slices` searches for cut indices - but only to decide where to put
ink, which CLAUDE.md places on the renderer's side of the line. Every value they read
was computed by analyzers.within_calibration_compute and lives in the artifact.

Kept out of within_calibration.py so the renderer file is the drawing sequence and this
file is the arithmetic that drawing needs.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from quebra.analyzers.fidelity import ABS_X_FLOOR
from quebra.analyzers.within_calibration_data import WithinCalibrationPanelData


def draw_decade_guides(ax: plt.Axes, values: np.ndarray) -> None:
    clipped = np.clip(np.asarray(values, dtype=float), 1e-16, None)
    lo, hi = float(np.min(clipped)), float(np.max(clipped))
    if not (np.isfinite(lo) and np.isfinite(hi) and lo > 0.0 and hi > 0.0):
        return
    for power in range(int(np.floor(np.log10(lo))), int(np.ceil(np.log10(hi))) + 1):
        ax.axhline(
            10.0**power,
            color="gray",
            linestyle="--",
            linewidth=0.7,
            alpha=0.3,
            zorder=0,
        )


def floor_reads(pd_: WithinCalibrationPanelData) -> np.ndarray:
    """Reads a log axis cannot draw as themselves: those at the |x| floor.

    The fidelity panel lifts |x| = 0, read 0 by construction, to `fidelity.ABS_X_FLOOR` so
    a log axis keeps it. That lift is a drawing floor, not a measurement, so the axis takes
    its range from the reads above it and marks the floor reads at its bottom edge. Only a
    log-scale series can hold one; every other panel gets an all-False mask.
    """
    values = np.asarray(pd_.signal.values, dtype=float)
    if not pd_.use_log_scale:
        return np.zeros(len(values), dtype=bool)
    return values == ABS_X_FLOOR


def lower_limit_above_floor(values: np.ndarray, at_floor: np.ndarray) -> float | None:
    """The decade boundary at or below the smallest finite read above the floor."""
    above = np.asarray(values, dtype=float)[~at_floor]
    above = above[np.isfinite(above) & (above > 0.0)]
    if not len(above):
        return None
    return float(10.0 ** np.floor(np.log10(float(np.min(above)))))


def observed_slices(pd_: WithinCalibrationPanelData) -> list[tuple[int, int]]:
    """[start, stop) index ranges of contiguously observed reads.

    Pure lookup of the gap boundaries the carve recorded - the renderer decides
    nothing about where a gap is, only that it does not draw through one.
    """
    if not pd_.signal.gap_spans_h:
        return [(0, len(pd_.signal.t_h))]
    cut_points = []
    for _t_before_h, t_after_h in pd_.signal.gap_spans_h:
        idx = int(np.searchsorted(pd_.signal.t_h, t_after_h, side="left"))
        if 0 < idx < len(pd_.signal.t_h):
            cut_points.append(idx)
    bounds = [0, *sorted(set(cut_points)), len(pd_.signal.t_h)]
    return [(lo, hi) for lo, hi in zip(bounds[:-1], bounds[1:]) if hi > lo]


def bin_centers_and_widths(
    edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    widths = np.diff(edges)
    return edges[:-1] + widths / 2.0, widths


# A symmetric excursion (peak_frac ~ 0.50) drives Spearman toward zero BY CONSTRUCTION,
# whatever the underlying predictability. Anything within this of 0.50 gets the note.
SYMMETRIC_PEAK_TOLERANCE = 0.08


def shape_annotation(
    medians: dict[str, float],
    defined: dict[str, int],
    peak_is_symmetric: bool = False,
) -> str:
    """The only way shape statistics reach a figure.

    Claim discipline, enforced rather than conventional: rho2 is NEVER emitted alone.
    Chatterjee's xi and the falling-limb rho always accompany it, because a near-zero
    rho2 on a symmetric excursion says nothing about predictability. `peak_is_symmetric`
    is decided in the analyzer band and handed in already resolved - the renderer states
    the condition, it does not classify.

    Every median carries the count it rests on: `limb_rho` in particular is undefined
    whenever the peak falls in the last three reads, which selects for left-peaked
    windows, so a bare median would overstate its own support.
    """
    if not medians:
        return "no complete windows"

    def cell(key: str, name: str) -> str:
        value = medians.get(key, float("nan"))
        n = defined.get(key, 0)
        return (
            f"{name}={value:.3g} (n={n})" if np.isfinite(value) else f"{name}=- (n={n})"
        )

    # Two lines at most: this sits under a figure axis, and a block tall enough to
    # collide with the next row is a block nobody reads. rho2_excess stays in the CSV.
    lines = [
        # rho2 first, but never without its two companions on the same line
        f"{cell('rho2', 'rho2')}  {cell('xi_m_to_f', 'xi')}  "
        f"{cell('limb_rho', 'limb rho')}  {cell('peak_frac', 'peak')}"
    ]
    if peak_is_symmetric:
        # States the condition; does not tell the reader what to conclude from it.
        lines.append("peak at mid-window: rho2 suppressed by symmetry")
    return "\n".join(lines)
