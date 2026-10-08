"""Carving a synthetic read series with the REAL window policy.

Arms B and C generate reads and need windows out of them. The one thing they must not do
is re-derive how a read series becomes windows: the gap threshold, the strict `>` test,
gap-wins-ties, death-at-first-out-of-spec and censored-dies-at-last-in-spec are decisions
that live in `analyzers/windows.py`, and a bench that reimplemented any of them would be
measuring a test against a carve the tool does not use.

So this calls `windows.spacing`, `windows.gap_flags`, `windows.in_spec_mask` and
`windows.carve` directly, plus the same duration/censoring derivation the analyzer applies
to their output. `tests/test_bench_uses_real_carve.py` asserts the result is column-for-
column identical to `windows.run`'s window table on the same series - that equality is
where the no-reimplementation claim is actually enforced, rather than in this comment.

Why not just call `windows.run`? Two reasons, both about the hot path: it builds a
12-column per-read DataFrame that the bench never looks at (at 4n reads per replicate over
~324k replicates that is the single largest cost in the run), and it prints unconditionally
on every call. Neither is a defect in `run` - the pipeline wants both - but neither belongs
in an inner loop.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from quebra.analyzers.windows import (
    DEATH_DOWN_CROSSING,
    DEFAULT_GAP_MULT,
    carve,
    gap_flags,
    in_spec_mask,
    spacing,
)

# The subset of `windows.WINDOW_COLUMNS` the checks consume. Deliberately a subset: the
# bench has no dataset_id or threshold ladder, and inventing values for those columns
# would make the equality test in tests/ pass on fabricated agreement.
BENCH_WINDOW_COLUMNS = [
    "window_index",
    "t_birth_s",
    "t_death_s",
    "duration_s",
    "birth_type",
    "death_type",
    "censored",
    "n_reads",
]


def carve_windows(
    t_rel_s: np.ndarray,
    values: np.ndarray,
    threshold_value: float,
    *,
    big_values_good: bool = True,
    gap_mult: float = DEFAULT_GAP_MULT,
) -> pd.DataFrame:
    """One threshold's window table, built from the real carve primitives."""
    t = np.asarray(t_rel_s, dtype=float)
    v = np.asarray(values, dtype=float)
    if len(t) != len(v):
        raise ValueError(f"t and values must match; got {len(t)} and {len(v)}")
    if not np.isfinite(threshold_value):
        raise ValueError(f"threshold has non-finite value {threshold_value!r}")
    # The next two mirror `windows.run` exactly, and are not optional decoration: the
    # claim this module makes is that its output is IDENTICAL to `run`'s, and `run` drops
    # non-finite pairs BEFORE measuring spacing (so one NaN moves the median and hence the
    # gap threshold) and refuses backwards time. Without them the bench carved 8 windows
    # where `run` carved 7, and a backwards step produced a negative duration_s.
    finite = np.isfinite(t) & np.isfinite(v)
    t, v = t[finite], v[finite]
    if np.any(np.diff(t) < 0.0):
        raise ValueError(
            "t_rel_s must be non-decreasing; the carve assumes observation order"
        )
    if len(t) < 2:
        return pd.DataFrame(columns=BENCH_WINDOW_COLUMNS)

    _median_s, gap_threshold_s, _n_gaps, _n_nonpositive = spacing(t, gap_mult)
    is_gap = gap_flags(t, gap_threshold_s)
    ins = in_spec_mask(v, threshold_value, big_values_good)
    windows = carve(t, ins, is_gap)

    rows: list[dict[str, object]] = []
    for index, w in enumerate(windows):
        s, e = int(w["s"]), int(w["e"])
        death_type = str(w["death_type"])
        t_birth_s = float(t[s])
        # Identical to analyzers/windows.py: a down_crossing dies at the first
        # out-of-spec read; anything else dies at its last in-spec read.
        t_death_s = (
            float(t[e]) if death_type == DEATH_DOWN_CROSSING else float(t[e - 1])
        )
        rows.append(
            {
                "window_index": index,
                "t_birth_s": t_birth_s,
                "t_death_s": t_death_s,
                "duration_s": t_death_s - t_birth_s,
                "birth_type": str(w["birth_type"]),
                "death_type": death_type,
                "censored": death_type != DEATH_DOWN_CROSSING,
                "n_reads": e - s,
            }
        )
    return pd.DataFrame(rows, columns=BENCH_WINDOW_COLUMNS)
