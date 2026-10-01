"""The design grid: every cell the bench evaluates, and what it means.

A CELL is one data-generating configuration at one event count on one clock. Each cell is
replicated many times and yields one ROW per (check, calibration, variant) - so the tables
are keyed by cell x row, and a check's size at n = 20 is a single number with a stated
Monte Carlo error rather than an impression.

The event-count, censoring and shape levels are set from what the real data shows,
measured before the bench was written; the dependence grid is a sweep, not a calibration to
a record:

- `N_GRID` ends at 355 because that is the largest usable event count in the record, and
  starts at 20 because below that `lag_layout` runs out of pairs. 35 is included because it
  is roughly the count at the thresholds that matter most.
- `CENSORING_GRID` is `{0, 0.03, 0.25}`. The real data sits at 0.000-0.026 wherever
  n >= 20, so 0 and 0.03 bracket it; 0.25 is deliberately outside it, and exists only to
  show what C1 and C2's censoring machinery buys when there is enough censoring to matter.
- `RHO_GRID` is dense at small `rho`, so the weak dependence a check is least likely to see
  is resolved finely. The duration-level lag-1 each `rho` induces is recorded per
  cell (`mean_induced_lag1`), and the power curve is read against that.
- `B_GRID` brackets 1 on both sides, so C1's power is measured against a decreasing
  intensity and an increasing one rather than only the direction that happens to be easier.
- `SHAPE_GRID` is `{0.75, 1.5}` - over- and under-dispersed relative to exponential
  respectively (coefficient of variation about 1.35 and 0.68) - because Kvaloy & Lindqvist
  (Technometrics 62(1) 2020, Section 6.1, Figure 1, p. 107) use exactly these two shapes and
  show the asymptotic calibrations degrade differently on the two sides.

Censoring is an Arm A/D factor only. Arms B and C carve a uniformly spaced read series
with no injected gaps, which can produce exactly one incomplete window (the terminal one),
so their censoring is ~1/n and emergent. Sweeping a censoring level there would silently
mean something different from what it means in Arm A, so it is not swept.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from jobs.bench.arms import ARM_A, ARM_B, ARM_C, ARM_D, ARM_E
from quebra.analyzers.checks.result import CLOCK_CALENDAR, CLOCK_IN_SPEC

KIND_SIZE = "size"
KIND_POWER = "power"

ALPHA = 0.05
N_PERM = 999
N_REPLICATES_SIZE = 2000
N_REPLICATES_POWER = 1000

N_GRID = (20, 35, 50, 75, 100, 355)
SHAPE_GRID = (0.75, 1.5)
CENSORING_GRID = (0.0, 0.03, 0.25)
RHO_GRID = (0.05, 0.1, 0.15, 0.2, 0.3, 0.5)
B_GRID = (0.7, 0.85, 1.2, 1.5)

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


def all_cells() -> list[Cell]:
    return [*size_cells(), *power_cells()]


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
