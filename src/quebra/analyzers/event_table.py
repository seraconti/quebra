"""The event table: the one input every survival estimator reads.

One table per record, threshold and side. Per window it holds the age the estimators use,
whether the death was observed, and the interval the true length is known to lie in. From
those it derives the risk table: each distinct age, the deaths and censorings there, and
the number at risk.

The estimators read this table and nothing else, so none re-reads the reads or re-carves.
The same table serves both sides of a threshold, and the constructor takes any other
right-censored durations directly and validates them.

Ages are seconds (`_s`). The risk table follows the Kaplan-Meier tie convention: at a tied
age deaths come first, so a window censored at age a is still at risk at a.

Validity assumptions

- Assumption: a window's KM age is `t_death_s - t_birth_s`: first in-window read to the
  first read past the death crossing, or to the last read for a censored window.
  Diagnostic: `windows.run` raises unless both sides and the gaps tile the record.
  Consequence of violation: every estimator reads ages on a different clock from the carve.
  Reference: `docs/WINDOW_SEMANTICS.md`.
- Assumption: windows whose birth was not observed carry no age, so they are excluded and
  counted (`n_unobserved_birth_dropped`), not entered as lifetimes.
  Diagnostic: the count travels on every result built from this table.
  Consequence of violation: a residual lifetime enters as a lifetime and survival is
  understated. Reference: `analyzers/kaplan_meier.py` module docstring.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import cast

import numpy as np
import pandas as pd

from quebra.analyzers.windows import SIDE_IN_SPEC, SIDE_OUT_OF_SPEC, WindowsResult

SIDES = (SIDE_IN_SPEC, SIDE_OUT_OF_SPEC)

# The window-table columns `from_windows` reads.
WINDOW_COLUMNS_NEEDED = (
    "threshold_label",
    "side",
    "t_before_birth_s",
    "t_birth_s",
    "t_last_s",
    "t_death_s",
    "duration_s",
    "birth_observed",
    "censored",
)


def risk_table(
    age: np.ndarray, event: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Distinct ages with their deaths, censorings and number at risk. Any unit.

    At risk at age a means age >= a, so a window censored at a is counted at a: the
    deaths-first tie convention.

    Two ages tie only when they are the same float. That is R `survival`'s
    `survfit(..., timefix = FALSE)`; R's default, `timefix = TRUE` (`aeqSurv`), also
    merges ages that differ by rounding noise, so near-ties a few ulps apart are two ages
    here and one there.

    `event` must be a numpy boolean array, as in `EventTable`: a bool cast would read any
    nonzero flag as a death.
    """
    if not isinstance(event, np.ndarray) or event.dtype != np.bool_:
        raise TypeError("event must be a numpy boolean array")
    age = np.asarray(age, dtype=float)
    if len(age) == 0:
        empty_f = np.array([], dtype=float)
        empty_i = np.array([], dtype=int)
        return empty_f, empty_i, empty_i.copy(), empty_i.copy()
    ages, inverse = np.unique(age, return_inverse=True)
    n_events = np.bincount(inverse, weights=event.astype(float), minlength=len(ages))
    n_total = np.bincount(inverse, minlength=len(ages))
    n_at_risk = len(age) - np.concatenate([[0], np.cumsum(n_total)[:-1]])
    n_events = n_events.astype(int)
    return ages, n_events, (n_total - n_events).astype(int), n_at_risk.astype(int)


@dataclass
class EventTable:
    """Right-censored durations with the interval each true length lies in.

    `age_s` is the KM age, the duration the estimators use. The window's true length lies
    in the CLOSED interval [age_lo_s, age_hi_s], with age_hi_s = inf for a censored window;
    duplicate timestamps can make it a single point (lo == hi), an exact length. `age_s`
    lies in the same interval, which the constructor checks. Durations known exactly may
    leave both bounds at None: lo is then `age_s`, and hi is `age_s` for a death and inf
    for a censored window.
    """

    dataset_id: str
    threshold_label: str
    side: str
    age_s: np.ndarray
    event: np.ndarray
    age_lo_s: np.ndarray | None = None
    age_hi_s: np.ndarray | None = None
    n_unobserved_birth_dropped: int = 0

    event_age_s: np.ndarray = field(init=False)
    n_events: np.ndarray = field(init=False)
    n_censored: np.ndarray = field(init=False)
    n_at_risk: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        if self.side not in SIDES:
            raise ValueError(f"side must be one of {SIDES}; got {self.side!r}")
        if not isinstance(self.event, np.ndarray) or self.event.dtype != np.bool_:
            # np.asarray(..., dtype=bool) coerces silently, so a float flag of 0.5 would
            # enter as an observed death. Require the caller to have meant a boolean.
            raise TypeError("event must be a numpy boolean array")
        # Copies, not views: the risk table is computed once, below, from these arrays.
        self.event = self.event.copy()
        self.age_s = np.array(self.age_s, dtype=float)
        n = len(self.age_s)
        if len(self.event) != n:
            raise ValueError(
                f"age_s ({n}) and event ({len(self.event)}) differ in length"
            )
        if not np.all(np.isfinite(self.age_s)) or np.any(self.age_s < 0.0):
            raise ValueError("age_s must be finite and non-negative")
        if self.age_lo_s is None:
            self.age_lo_s = self.age_s.copy()
        if self.age_hi_s is None:
            self.age_hi_s = np.where(self.event, self.age_s, np.inf)
        self.age_lo_s = np.array(self.age_lo_s, dtype=float)
        self.age_hi_s = np.array(self.age_hi_s, dtype=float)
        if len(self.age_lo_s) != n or len(self.age_hi_s) != n:
            raise ValueError("age_lo_s and age_hi_s must match age_s in length")
        # Every comparison with NaN is False, so a NaN bound would pass the checks below.
        if not np.all(np.isfinite(self.age_lo_s)) or np.any(np.isnan(self.age_hi_s)):
            raise ValueError("age_lo_s must be finite and age_hi_s must not be NaN")
        if np.any(self.age_lo_s > self.age_s) or np.any(self.age_s > self.age_hi_s):
            raise ValueError("every age must lie in its interval: lo <= age <= hi")
        if np.any(np.isinf(self.age_hi_s) & self.event):
            raise ValueError("an observed death needs a finite age_hi_s")
        if np.any(np.isfinite(self.age_hi_s) & ~self.event):
            raise ValueError("a censored window has age_hi_s = inf")
        (
            self.event_age_s,
            self.n_events,
            self.n_censored,
            self.n_at_risk,
        ) = risk_table(self.age_s, self.event)

    @property
    def n_windows(self) -> int:
        return int(len(self.age_s))

    @property
    def n_deaths(self) -> int:
        return int(np.count_nonzero(self.event))


def from_windows(
    windows: pd.DataFrame,
    *,
    threshold_label: str,
    side: str,
    dataset_id: str = "",
    known_labels: set[str] | None = None,
) -> EventTable:
    """One threshold's windows on one side, births observed, as an event table.

    A birth lies in [t_before_birth_s, t_birth_s] and a death in [t_last_s, t_death_s], so
    a death's length lies in [t_last_s - t_birth_s, t_death_s - t_before_birth_s]; a
    censored window's in [t_last_s - t_birth_s, inf]. The KM age is `duration_s`.

    `side` must match the table's own `side` column. `known_labels` is the carve's ladder:
    with it, a label on the ladder with no windows on this side gives an empty table, and a
    label off the ladder raises. Without it the label must appear in the table, because a
    typo and an empty side would otherwise look the same.
    """
    missing = [c for c in WINDOW_COLUMNS_NEEDED if c not in windows.columns]
    if missing:
        raise KeyError(f"the window table lacks {missing}; build it with windows.run")
    if side not in SIDES:
        raise ValueError(f"side must be one of {SIDES}; got {side!r}")
    present = set(windows["threshold_label"].unique())
    ladder = known_labels if known_labels is not None else present
    if threshold_label not in ladder:
        raise KeyError(
            f"threshold {threshold_label!r} is not on the ladder: {sorted(ladder)}"
        )
    for flag in ("birth_observed", "censored"):
        # to_numpy(dtype=bool) would read NaN and the string "False" as True. An empty
        # frame has no dtype to check (pandas gives an empty column object dtype).
        if len(windows) and windows[flag].dtype != np.bool_:
            raise TypeError(
                f"window column {flag!r} must be bool; got {windows[flag].dtype}"
            )
    rows = windows[windows["threshold_label"] == threshold_label]
    wrong_side = set(rows["side"].unique()) - {side}
    if wrong_side:
        raise ValueError(
            f"asked for the {side} side but the table holds {sorted(wrong_side)} windows"
        )
    seen = rows[rows["birth_observed"].to_numpy(dtype=bool)]
    event = ~seen["censored"].to_numpy(dtype=bool)
    t_birth = seen["t_birth_s"].to_numpy(dtype=float)
    age_lo = seen["t_last_s"].to_numpy(dtype=float) - t_birth
    age_hi = np.where(
        event,
        seen["t_death_s"].to_numpy(dtype=float)
        - seen["t_before_birth_s"].to_numpy(dtype=float),
        np.inf,
    )
    return EventTable(
        dataset_id=dataset_id,
        threshold_label=threshold_label,
        side=side,
        age_s=seen["duration_s"].to_numpy(dtype=float),
        event=event,
        age_lo_s=age_lo,
        age_hi_s=age_hi,
        n_unobserved_birth_dropped=int(len(rows) - len(seen)),
    )


def from_carve(carved: WindowsResult, *, threshold_label: str, side: str) -> EventTable:
    """`from_windows` on a `windows.WindowsResult`, with the carve's own ladder."""
    return _from_carve_on_ladder(carved, threshold_label, side, set(_ladder(carved)))


@dataclass
class EventTables:
    """One event table per (threshold label, side), from one carve, on its whole ladder."""

    ladder: tuple[str, ...]
    tables: dict[tuple[str, str], EventTable]

    def table(self, threshold_label: str, side: str) -> EventTable:
        key = (threshold_label, side)
        if key not in self.tables:
            raise KeyError(f"no event table for {key}; ladder {list(self.ladder)}")
        return self.tables[key]


def event_tables_from_carve(carved: WindowsResult) -> EventTables:
    """Every (label, side) of a `windows.WindowsResult`, empty sides included.

    The ladder is the carve's own (`diagnostics["per_threshold"]`), so a threshold the
    record had too few reads to carve is still on it, with two empty tables.
    """
    ladder = _ladder(carved)
    tables = {
        (label, side): _from_carve_on_ladder(carved, label, side, set(ladder))
        for label in ladder
        for side in SIDES
    }
    return EventTables(ladder=ladder, tables=tables)


def _ladder(carved: WindowsResult) -> tuple[str, ...]:
    # `diagnostics` is typed dict[str, object]; `windows.run` keys this entry by label.
    return tuple(cast("dict[str, object]", carved.diagnostics["per_threshold"]))


def _from_carve_on_ladder(
    carved: WindowsResult, label: str, side: str, ladder: set[str]
) -> EventTable:
    table = (
        carved.windows_in_spec if side == SIDE_IN_SPEC else carved.windows_out_of_spec
    )
    return from_windows(
        table,
        threshold_label=label,
        side=side,
        dataset_id=str(carved.meta.get("dataset_id", "")),
        known_labels=ladder,
    )
