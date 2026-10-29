"""The design grid: every cell the bench evaluates, and what it means.

A CELL is one data-generating configuration at one event count on one clock. Each cell is
replicated many times and yields one ROW per (check, calibration, variant) - so the tables
are keyed by cell x row, and a check's size at n = 20 is a single number with a stated
Monte Carlo error rather than an impression.

The event-count levels cover the 6D2S check-ledger records; the censoring and shape levels
are set from what the real data showed, measured before the bench was written; the
dependence grid is a sweep, not a calibration to a record:

- `N_GRID` starts at 20 because below that `lag_layout` runs out of pairs, and ends at 1000
  because the 6D2S check-ledger records (`jobs/active/check_ledger_6d2s_q*.py`) reach 936
  events; above 355 the points are spaced so every one of their event counts is within 150
  of a grid point. 35 is included because it is roughly the count at the thresholds that
  matter most. A cell's seed is its index in `all_cells()`, so a new n goes at the END:
  the size cells before it keep their seeds, and every power cell, which follows all the
  size cells, is re-seeded.
- `CENSORING_GRID` is `{0, 0.03, 0.25}`. The real data sits at 0.000-0.026 wherever
  n >= 20, so 0 and 0.03 bracket it; 0.25 is deliberately outside it, and exists only to
  show what C1 and C2's censoring machinery buys when there is enough censoring to matter.
- `RHO_GRID` is dense at small `rho`, so the weak dependence a check is least likely to see
  is resolved finely. The power curve is read against the rank correlation of neighbouring
  durations each `rho` induces, `(6/pi) arcsin(rho/2)`; `mean_induced_lag1` records C5's
  uncentred `r_1` per cell, which reads about -1/N with no dependence, N the cell's event
  count.
- `B_GRID` brackets 1 on both sides, so C1's power is measured against a decreasing
  intensity and an increasing one rather than only the direction that happens to be easier.
- `SHAPE_GRID` is `{0.75, 1.5}` - over- and under-dispersed relative to exponential
  respectively (coefficient of variation about 1.35 and 0.68) - because Kvaloy & Lindqvist
  (Technometrics 62(1) 2020, Section 6.1, Figure 1, p. 107) use exactly these two shapes and
  show the asymptotic calibrations degrade differently on the two sides.

Censoring is an Arm A/D factor only in the size and power cells. Arms B and C carve a
uniformly spaced read series with no injected gaps, which can produce exactly one
incomplete window (the terminal one), so their censoring is ~1/n and emergent. Sweeping a
censoring level there would silently mean something different from what it means in Arm A,
so it is not swept.

The SERIAL cells (`serial_cells`, kind `serial`, their own table) put the rank checks on the
layouts the size and power cells leave out:

- Arm E split into a target of 3, 5 and 10 segments (`SERIAL_SEGMENT_GRID`) at n = 50 and
  60, through the censoring knob, over `RHO_GRID` and `rho = 0`;
- Arm E2, dependence at lag 2 only, on the same layouts over `RHO_GRID`;
- Arm EL, two segments of `LEVEL_SHIFT_EVENTS_PER_SEGMENT` expected events with the second
  shifted by `delta`, over `LEVEL_SHIFT_DELTA_GRID` x `LEVEL_SHIFT_PHI_GRID`, the latent
  AR(1) coefficient carried in the `rho` column. `phi = 0` is a null cell at every `delta`.

They use one Weibull shape, `SERIAL_SHAPE` = 0.75. C5 and C6 read ranks, and the copula
map to durations is monotone, so the shape reaches them only through truncation: how many
events a segment keeps before its `tau`. The over-dispersed shape makes that count the more
variable of the two (the count variance of a renewal process scales with the squared
coefficient of variation, about 1.8 against 0.46), so it is the harder layout to read.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from jobs.bench.arms import (
    ARM_A,
    ARM_B,
    ARM_C,
    ARM_D,
    ARM_E,
    ARM_E2,
    ARM_EL,
    n_segments_for_censoring,
)
from quebra.analyzers.checks.result import CLOCK_CALENDAR, CLOCK_IN_SPEC

KIND_SIZE = "size"
KIND_POWER = "power"
KIND_SERIAL = "serial"

ALPHA = 0.05
N_PERM = 999
N_REPLICATES_SIZE = 2000
N_REPLICATES_POWER = 1000

N_GRID = (20, 35, 50, 75, 100, 355, 500, 700, 1000)
SHAPE_GRID = (0.75, 1.5)
CENSORING_GRID = (0.0, 0.03, 0.25)
RHO_GRID = (0.05, 0.1, 0.15, 0.2, 0.3, 0.5)
B_GRID = (0.7, 0.85, 1.2, 1.5)

SERIAL_SHAPE = 0.75
SERIAL_N_GRID = (50, 60)
SERIAL_SEGMENT_GRID = (3, 5, 10)
LEVEL_SHIFT_EVENTS_PER_SEGMENT = 150
LEVEL_SHIFT_DELTA_GRID = (0.0, 0.4, 0.8)
LEVEL_SHIFT_PHI_GRID = (-0.12, 0.0, 0.12)

MASTER_SEED = 20260806


@dataclass(frozen=True)
class Cell:
    arm: str
    n: int
    clock: str
    kind: str
    shape: float | None = None
    quantised: bool = False
    censoring: float = 0.0
    rho: float | None = None
    b: float | None = None
    delta: float | None = None
    # False on the in-spec clock of a carved arm, where tau is not well posed - see
    # `analyzers.checks.battery.run_battery`.
    include_tau_checks: bool = True
    n_replicates: int = N_REPLICATES_SIZE

    def key(self) -> tuple:
        return (
            self.arm,
            self.n,
            self.clock,
            self.kind,
            self.shape,
            self.quantised,
            self.censoring,
            self.rho,
            self.b,
            self.delta,
        )

    def label(self) -> str:
        parts = [self.arm, f"n={self.n}", self.clock, self.kind]
        if self.shape is not None:
            parts.append(f"shape={self.shape}")
        if self.quantised:
            parts.append("quantised")
        if self.censoring:
            parts.append(f"c={self.censoring}")
        if self.rho is not None:
            parts.append(f"rho={self.rho}")
        if self.b is not None:
            parts.append(f"b={self.b}")
        if self.delta is not None:
            parts.append(f"delta={self.delta}")
        return " ".join(parts)


def _carved_clocks(arm: str, kind: str, n: int, **kwargs) -> list[Cell]:
    """One cell per clock for a carved arm.

    The calendar clock gets the full battery; the in-spec clock gets the rank checks and CvM
    by permutation, the rows that stay defined when `tau == T_N`.
    That asymmetry is not a preference - on the in-spec clock of a carved record, time
    stops accumulating when the record ends out of spec, so `tau == T_N` on ~73% of
    replicates and C2 is genuinely undefined there. Reporting C1/C2 from the surviving
    replicates would report a size conditioned on how the record happened to end.
    """
    reps = N_REPLICATES_SIZE if kind == KIND_SIZE else N_REPLICATES_POWER
    return [
        Cell(
            arm=arm,
            n=n,
            clock=CLOCK_CALENDAR,
            kind=kind,
            include_tau_checks=True,
            n_replicates=reps,
            **kwargs,
        ),
        Cell(
            arm=arm,
            n=n,
            clock=CLOCK_IN_SPEC,
            kind=kind,
            include_tau_checks=False,
            n_replicates=reps,
            **kwargs,
        ),
    ]


def size_cells() -> list[Cell]:
    """Every cell whose data-generating process satisfies the checks' null."""
    cells: list[Cell] = []
    for n in N_GRID:
        for shape in SHAPE_GRID:
            for quantised in (False, True):
                for censoring in CENSORING_GRID:
                    cells.append(
                        Cell(
                            arm=ARM_A,
                            n=n,
                            clock=CLOCK_IN_SPEC,
                            kind=KIND_SIZE,
                            shape=shape,
                            quantised=quantised,
                            censoring=censoring,
                        )
                    )
            # b = 1 is the identity trend, so Arm D reduces to a renewal process and this
            # is a second, independently coded route to the same null - a check on the
            # TRP generator itself as much as on the tests.
            for censoring in CENSORING_GRID:
                cells.append(
                    Cell(
                        arm=ARM_D,
                        n=n,
                        clock=CLOCK_IN_SPEC,
                        kind=KIND_SIZE,
                        shape=shape,
                        censoring=censoring,
                        b=1.0,
                    )
                )
            # rho = 0 is the null anchor the Arm E power curve is read against.
            cells.append(
                Cell(
                    arm=ARM_E,
                    n=n,
                    clock=CLOCK_IN_SPEC,
                    kind=KIND_SIZE,
                    shape=shape,
                    rho=0.0,
                )
            )
        cells.extend(_carved_clocks(ARM_B, KIND_SIZE, n))
        cells.extend(_carved_clocks(ARM_C, KIND_SIZE, n, rho=0.0))
    return cells


def power_cells() -> list[Cell]:
    """Every cell whose data-generating process violates a null the checks test."""
    cells: list[Cell] = []
    for n in N_GRID:
        for shape in SHAPE_GRID:
            for b in B_GRID:
                cells.append(
                    Cell(
                        arm=ARM_D,
                        n=n,
                        clock=CLOCK_IN_SPEC,
                        kind=KIND_POWER,
                        shape=shape,
                        b=b,
                        n_replicates=N_REPLICATES_POWER,
                    )
                )
            for rho in RHO_GRID:
                cells.append(
                    Cell(
                        arm=ARM_E,
                        n=n,
                        clock=CLOCK_IN_SPEC,
                        kind=KIND_POWER,
                        shape=shape,
                        rho=rho,
                        n_replicates=N_REPLICATES_POWER,
                    )
                )
        for rho in RHO_GRID:
            # Expected to be FLAT at the nominal level: read-level dependence does not
            # propagate to durations. Kept at full resolution because that flatness is the
            # finding, and a finding needs the same evidence a positive result would.
            cells.extend(_carved_clocks(ARM_C, KIND_POWER, n, rho=rho))
    return cells


def censoring_for_segments(n: int, n_segments: int) -> float:
    """The censoring target at which `arms.n_segments_for_censoring(n, .)` is `n_segments`:
    one incomplete gap per segment, so `c = m / (n + m)`. Raises when that function caps
    the count below `n_segments`, which it does once a segment would expect fewer than
    `arms.MIN_EXPECTED_EVENTS_PER_SEGMENT` events."""
    censoring = n_segments / (n + n_segments)
    realised = n_segments_for_censoring(n, censoring)
    if realised != n_segments:
        raise ValueError(
            f"{n_segments} segments at n = {n} realise as {realised}; "
            "the cell would be labelled with a layout it does not have"
        )
    return censoring


def serial_cells() -> list[Cell]:
    """Gapped and lag-2 layouts for the rank checks, and the level shift; see the module
    docstring. A `rho = 0` (or `phi = 0`) cell is a null cell and gets the size replicates."""

    def replicates(rho: float) -> int:
        return N_REPLICATES_SIZE if rho == 0.0 else N_REPLICATES_POWER

    cells: list[Cell] = []
    for n in SERIAL_N_GRID:
        for n_segments in SERIAL_SEGMENT_GRID:
            censoring = censoring_for_segments(n, n_segments)
            for arm, rhos in ((ARM_E, (0.0, *RHO_GRID)), (ARM_E2, RHO_GRID)):
                for rho in rhos:
                    cells.append(
                        Cell(
                            arm=arm,
                            n=n,
                            clock=CLOCK_IN_SPEC,
                            kind=KIND_SERIAL,
                            shape=SERIAL_SHAPE,
                            censoring=censoring,
                            rho=rho,
                            n_replicates=replicates(rho),
                        )
                    )
    for delta in LEVEL_SHIFT_DELTA_GRID:
        for phi in LEVEL_SHIFT_PHI_GRID:
            cells.append(
                Cell(
                    arm=ARM_EL,
                    n=2 * LEVEL_SHIFT_EVENTS_PER_SEGMENT,
                    clock=CLOCK_IN_SPEC,
                    kind=KIND_SERIAL,
                    shape=SERIAL_SHAPE,
                    rho=phi,
                    delta=delta,
                    n_replicates=replicates(phi),
                )
            )
    return cells


def all_cells() -> list[Cell]:
    # Serial cells go last, so every size and power cell keeps its index and its seed.
    return [*size_cells(), *power_cells(), *serial_cells()]


@dataclass(frozen=True)
class GridSummary:
    n_size_cells: int
    n_power_cells: int
    n_replicates: int
    notes: list[str] = field(default_factory=list)


def summarise() -> GridSummary:
    size, power = size_cells(), power_cells()
    return GridSummary(
        n_size_cells=len(size),
        n_power_cells=len(power),
        n_replicates=sum(c.n_replicates for c in [*size, *power]),
    )
