"""The in-spec fraction of every record on the T2* ladder, with its bounds.

The compliance survey. For every record a ledger job carved and every rung of the ladder,
three numbers read off that carve alone:

- the time-weighted in-spec fraction, the panel's occupancy
  (`within_calibration_compute._threshold_in_spec_frac`, the one occupancy definition), so
  a bar here is the number the record's T2* panel reports;
- its lower and upper bounds (`tellability.in_spec_fraction_bounds`), at the gap rule and
  the confidence level k the carve recorded in `WindowsResult.meta`, so they are the bounds
  the record's T2* tellability reports.

The fraction lies between its bounds by construction. Each interval between two reads of a
segment is held by its left read: the fraction counts it in spec when that read is in
spec, the lower bound only when it is also settled, the upper bound whenever it is in spec
or unsettled. A fraction outside its bounds is a defect, and raises.

A record whose load or carve failed arrives as a `RecordFailure` and reads `not computed`
with its note at every rung; a rung the record cannot answer reads `not computed:
<reason>`. Neither is reported as 0, and no record is dropped.

The job declares which records the paper's dataset table holds and which of its rows the
paper shades; both are step kwargs, so they reach provenance. The builder checks them
against the records the ledger jobs carved. A carved record outside the table stays on the
artifact, marked, and the figure lists it instead of drawing it.

Pure compute: carves in, typed artifact out.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from quebra.analyzers.tellability import (
    COMPUTED,
    NotComputed,
    in_spec_fraction_bounds,
)
from quebra.analyzers.windows import WindowsResult, margin, observed_reads
from quebra.analyzers.within_calibration_compute import _threshold_in_spec_frac
from quebra.core._artifact_guard import StaleArtifactGuard
from quebra.core.types import RecordFailure

# How far a fraction may sit outside its bounds before it is a defect: summation order
# alone moves either side by a few ulps.
_BOUNDS_TOL = 1e-12


@dataclass(frozen=True)
class RungCompliance:
    """One rung of one record: the time-weighted in-spec fraction and its bounds.

    `status` belongs to the fraction and `bounds_status` to the two bounds; each is
    `computed` or `not computed: <reason>`, and its numbers are None unless computed.
    """

    threshold_label: str
    threshold_value: float
    big_values_good: bool
    in_spec_fraction: float | None
    in_spec_fraction_lower: float | None
    in_spec_fraction_upper: float | None
    status: str
    bounds_status: str


@dataclass
class RecordCompliance(StaleArtifactGuard):
    """Every rung of one record, in ladder order, with the carve settings it was read at.

    `gap_mult` and `unsure_k` are the carve's own, None when the record failed;
    `record_failure` is the failure's note, empty unless it failed.
    """

    dataset_id: str
    rungs: tuple[RungCompliance, ...]
    gap_mult: float | None
    unsure_k: float | None
    record_failure: str = ""


@dataclass(frozen=True)
class SurveyRecord:
    """One record of the survey: its run name, row label, device and the paper's marks."""

    record: str
    label: str
    device: str
    in_paper_table: bool
    # The flag text when the paper's table shades this record, else empty.
    flag: str
    compliance: RecordCompliance


@dataclass(frozen=True)
class PanelCounts:
    """The rungs of one device's drawn records, and how many of them read `not computed`:
    the fraction (`fraction_not_computed`), or the bounds alone beside a computed fraction
    (`bounds_not_computed`)."""

    n_rungs: int
    fraction_not_computed: int
    bounds_not_computed: int


@dataclass
class ComplianceSurveyData(StaleArtifactGuard):
    """Every record the ledger jobs carved, in their order, on one ladder.

    `devices` is the panel order and `counts` each panel's tally of `not computed` rungs.
    `gap_mult` and `unsure_k` are the one gap rule and the one confidence level every
    computed record was read at (None when no record was).
    """

    thresholds: tuple[tuple[str, float, bool], ...]
    devices: tuple[str, ...]
    records: tuple[SurveyRecord, ...]
    counts: dict[str, PanelCounts]
    gap_mult: float | None
    unsure_k: float | None
    flag_text: str
    table_source: str

    def drawn(self, device: str) -> tuple[SurveyRecord, ...]:
        """The records of `device` in the paper's dataset table, in survey order."""
        return tuple(r for r in self.records if r.device == device and r.in_paper_table)

    def not_drawn(self, device: str) -> tuple[SurveyRecord, ...]:
        """The records of `device` the ledger jobs carved outside the paper's table."""
        return tuple(
            r for r in self.records if r.device == device and not r.in_paper_table
        )


def _ladder(
    thresholds: Sequence[Sequence[object]],
) -> list[tuple[str, float, bool]]:
    return [(str(r[0]), float(r[1]), bool(r[2])) for r in thresholds]  # type: ignore[arg-type]


def _check_carved_ladder(
    carve: WindowsResult, ladder: list[tuple[str, float, bool]]
) -> None:
    """The carve must be cut on exactly this ladder: the same labels in order, and on
    every label with reads the margins `windows.margin` gives at this value and side."""
    carved = list(carve.diagnostics.get("per_threshold", {}))
    labels = [label for label, _, _ in ladder]
    if carved != labels:
        raise ValueError(
            f"{carve.meta['dataset_id']!r} was carved on {carved}, not on the survey's "
            f"ladder {labels}"
        )
    reads = carve.reads
    for label, value, big_values_good in ladder:
        sub = reads[reads["threshold_label"] == label]
        if len(sub) and not np.array_equal(
            margin(sub["value"].to_numpy(dtype=float), value, big_values_good),
            sub["margin"].to_numpy(dtype=float),
        ):
            raise ValueError(
                f"{carve.meta['dataset_id']!r}: threshold {label!r} was carved at another "
                f"value or side than ({value}, {big_values_good})"
            )


def record_compliance(
    carve: WindowsResult | RecordFailure, thresholds: Sequence[Sequence[object]]
) -> RecordCompliance:
    """The in-spec fraction and its bounds at every rung of `thresholds`, off one carve.

    Raises ValueError when the carve was cut on another ladder or without its sigma
    (`use_uncertainty` off, so its read table holds no sigma for the bounds to read), and
    RuntimeError when its read table does not hold the carve's own reads or a fraction
    falls outside its bounds: all are defects, not properties of the record.

    Validity assumptions:

    - Assumption: no excursion to the other side happens between two reads on the same
      side.
      Diagnostic: the count exponent and rho_k rows of the record's tellability.
      Consequence of violation: both the fraction and its bounds describe the time held
      by reads in spec, not the true time in spec.
      Reference: `tellability.in_spec_fraction_bounds`; no source located.
    - Assumption: the carve's sigma is a Gaussian 1-sigma, so k reads as a confidence
      level.
      Diagnostic: the self-consistency row (Q3.2) of the record's tellability.
      Consequence of violation: the bounds are at a nominal level only; the fraction is
      unaffected, since it reads no sigma.
      Reference: `tellability.in_spec_fraction_bounds`; no source located.
    """
    ladder = _ladder(thresholds)
    if isinstance(carve, RecordFailure):
        status = f"not computed: {carve.note()}"
        return RecordCompliance(
            dataset_id=carve.dataset_id,
            rungs=tuple(
                RungCompliance(label, value, big, None, None, None, status, status)
                for label, value, big in ladder
            ),
            gap_mult=None,
            unsure_k=None,
            record_failure=carve.note(),
        )
    if not isinstance(carve, WindowsResult):
        raise TypeError(
            f"{type(carve).__name__} is neither a carve nor a recorded failure"
        )
    _check_carved_ladder(carve, ladder)
    dataset_id = str(carve.meta["dataset_id"])
    if not carve.meta["use_uncertainty"]:
        raise ValueError(
            f"{dataset_id!r} was carved with use_uncertainty off: its read table holds no "
            "sigma, so its bounds would read none where the record may declare one"
        )
    gap_mult = float(carve.meta["gap_mult"])
    k = float(carve.meta["k"])
    n_reads = int(carve.diagnostics["n_reads_finite"])
    gap_spans_h = [
        (lo / 3600.0, hi / 3600.0) for lo, hi in carve.diagnostics["gap_spans_s"]
    ]
    reads = carve.reads
    rungs = []
    for label, value, big_values_good in ladder:
        sub = reads[reads["threshold_label"] == label]
        if len(sub) != (n_reads if n_reads >= 2 else 0):
            raise RuntimeError(
                f"{dataset_id!r}: the read table holds {len(sub)} reads at {label!r}, "
                f"the carve read {n_reads}"
            )
        t_s = sub["t_read_s"].to_numpy(dtype=float)
        values = sub["value"].to_numpy(dtype=float)
        fraction = _threshold_in_spec_frac(
            t_s / 3600.0, values, reads, [(label, value, big_values_good)], gap_spans_h
        )[label]
        status = COMPUTED
        if not math.isfinite(fraction):
            status = NotComputed("no segment with at least 2 reads").status
        obs = observed_reads(
            t_s,
            values,
            sigma=pd.to_numeric(sub["sigma_v"], errors="coerce").to_numpy(dtype=float),
            gap_mult=gap_mult,
        )
        lower = upper = None
        try:
            bounds = in_spec_fraction_bounds(
                obs, sub["margin"].to_numpy(dtype=float), obs.sigma, k
            )
        except NotComputed as exc:
            bounds_status = exc.status
        else:
            bounds_status = COMPUTED
            lower = bounds.in_spec_fraction_lower
            upper = bounds.in_spec_fraction_upper
            if status == COMPUTED and not (
                lower - _BOUNDS_TOL <= fraction <= upper + _BOUNDS_TOL
            ):
                raise RuntimeError(
                    f"{dataset_id!r} at {label!r}: the in-spec fraction {fraction!r} lies "
                    f"outside its bounds [{lower!r}, {upper!r}]"
                )
        rungs.append(
            RungCompliance(
                threshold_label=label,
                threshold_value=value,
                big_values_good=big_values_good,
                in_spec_fraction=fraction if status == COMPUTED else None,
                in_spec_fraction_lower=lower,
                in_spec_fraction_upper=upper,
                status=status,
                bounds_status=bounds_status,
            )
        )
    return RecordCompliance(
        dataset_id=dataset_id, rungs=tuple(rungs), gap_mult=gap_mult, unsure_k=k
    )


def _panel_counts(rows: Sequence[SurveyRecord]) -> PanelCounts:
    rungs = [rung for row in rows for rung in row.compliance.rungs]
    return PanelCounts(
        n_rungs=len(rungs),
        fraction_not_computed=sum(rung.status != COMPUTED for rung in rungs),
        bounds_not_computed=sum(
            rung.status == COMPUTED and rung.bounds_status != COMPUTED for rung in rungs
        ),
    )


def _one(values: list[float], what: str) -> float | None:
    distinct = sorted(set(values))
    if len(distinct) > 1:
        raise ValueError(
            f"the records were carved at {what} {distinct}; one survey, one"
        )
    return distinct[0] if distinct else None


def build_compliance_survey(
    *carve_sets: object,
    records: tuple[str, ...],
    devices: tuple[str, ...],
    paper_records: tuple[tuple[str, str], ...],
    shaded: tuple[str, ...],
    flag_text: str,
    table_source: str,
) -> ComplianceSurveyData:
    """The survey of every carve in `carve_sets`, records in `records` order.

    `carve_sets` are the ledger jobs' carve bundles (`ledger_recipe.RecordCarves`: their
    `records` and `carves`, in order, and the `thresholds` they were cut on). The survey
    reads its ladder from the bundles, which must share one. `records` and `devices` are
    declared when the graph
    is built, so the row order and each record's device reach the provenance label;
    `records` must equal the bundles' own run names in order. `paper_records` holds the
    paper's table as (device, run name) pairs, and `shaded` the run names of its shaded
    rows. Raises when a record of the paper's table was not carved, was carved on another
    device, or a shaded record is not in the table: a survey silently missing a row of the
    paper's table is the failure this guards.
    """
    bundled = tuple(r for s in carve_sets for r in s.records)  # type: ignore[attr-defined]
    if bundled != tuple(records):
        raise ValueError(f"declared rows {records} differ from the bundles' {bundled}")
    carves = [c for s in carve_sets for c in s.carves]  # type: ignore[attr-defined]
    if len(devices) != len(records) or len(carves) != len(records):
        raise ValueError(
            f"{len(records)} records, {len(devices)} devices, {len(carves)} carves"
        )
    if len(set(records)) != len(records):
        raise ValueError(f"a record is repeated: {records}")
    device_of = dict(zip(records, devices))
    table = {}
    for device, record in paper_records:
        if record not in device_of:
            raise ValueError(
                f"the paper's table holds {record!r} ({device}), which no ledger job carved"
            )
        if device_of[record] != device:
            raise ValueError(
                f"the paper's table lists {record!r} on {device}, but it was carved as a "
                f"{device_of[record]} record"
            )
        table[record] = device
    outside = sorted(set(shaded) - set(table))
    if outside:
        raise ValueError(f"shaded records {outside} are not in the paper's table")
    ladders = {tuple(_ladder(s.thresholds)) for s in carve_sets}  # type: ignore[attr-defined]
    if len(ladders) != 1:
        raise ValueError(f"the bundles were cut on {len(ladders)} ladders, not one")
    (ladder,) = (list(lad) for lad in ladders)
    rows = []
    for record, device, carve in zip(records, devices, carves):
        compliance = record_compliance(carve, ladder)
        if compliance.dataset_id != record:
            raise ValueError(
                f"row {record!r} holds the carve of {compliance.dataset_id!r}"
            )
        rows.append(
            SurveyRecord(
                record=record,
                label=record.replace("_", " "),
                device=device,
                in_paper_table=record in table,
                flag=flag_text if record in shaded else "",
                compliance=compliance,
            )
        )
    computed = [r.compliance for r in rows if not r.compliance.record_failure]
    panels = tuple(dict.fromkeys(devices))
    return ComplianceSurveyData(
        thresholds=tuple(ladder),
        devices=panels,
        records=tuple(rows),
        counts={
            device: _panel_counts(
                [r for r in rows if r.device == device and r.in_paper_table]
            )
            for device in panels
        },
        gap_mult=_one([c.gap_mult for c in computed], "gap_mult"),  # type: ignore[misc]
        unsure_k=_one([c.unsure_k for c in computed], "k"),  # type: ignore[misc]
        flag_text=flag_text,
        table_source=table_source,
    )
