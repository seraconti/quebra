"""Two devices' calibration cycles compared as records, descriptively.

`compare_cycle_pair` reads each device's cycles as an event table (`cycles.cycle_event_table`)
and hands the two to `record_comparison.compare_pair` with `test=False`: one comparison tau
by the pair's rule, each device's restricted mean to it with its own interval, and the
point difference, with no log-rank and no SE for the difference. Each device's check
outcome is the one its cycles figure draws beside the MTBC and the band, the within-stretch
standing (`CyclesPanelData.checks`), attached to its side of the pair and to its summary.

The restricted means run to the shared tau, the smaller of the two devices' largest observed
cycles. A device's MTBC runs to its own largest cycle (`cycles.build_cycles_panel_data`), so
the device that sets the shared tau is restricted where its MTBC is, and its restricted mean
and SE are its MTBC's. The other device is restricted short of its own largest cycle, unless
the two are equal, and its curve is above zero up to that cycle, so its restricted mean is
smaller than its MTBC.

The result is a `ByThresholdSide` of one entry, the cycles' empty threshold label and side
`event_table.CYCLES`, so `record_comparison.pair_table` reads it as it reads a ladder.

A module of its own so that neither `cycles` nor `record_comparison` imports the other:
whatever a module imports joins the identity closure of its steps, and each stays out of the
other's jobs.

Validity assumptions

- Assumption: each device's cycles are exchangeable draws from one law (a2).
  Diagnostic: the device's check outcome, attached to its side of the pair.
  Consequence of violation: that device's interval does not hold; the restricted means and
  the point difference are descriptive either way.
  Reference: `assumptions.A2_RENEWAL_CYCLES`.
- The assumptions of `record_comparison`, but the independence of the two records: no
  number here rests on it, since nothing is tested between the devices.
"""

from __future__ import annotations

from dataclasses import replace

from quebra.analyzers.cycles import CycleTable, CyclesPanelData, cycle_event_table
from quebra.analyzers.event_table import CYCLES
from quebra.analyzers.kaplan_meier import CONF_LEVEL
from quebra.analyzers.record_comparison import PairComparison, compare_pair
from quebra.analyzers.survival import ByThresholdSide, CurveSummaries


def compare_cycle_pair(
    cycles_a: CycleTable,
    cycles_b: CycleTable,
    panel_a: CyclesPanelData,
    panel_b: CyclesPanelData,
    *,
    tau_rule: str,
    conf_level: float = CONF_LEVEL,
) -> ByThresholdSide[PairComparison]:
    """Device a against device b on their cycles, descriptive only (module docstring).

    `panel_a` and `panel_b` are each device's cycles figure data, read for its check
    outcome. Raises when a panel is another record's.
    """
    for table, panel in ((cycles_a, panel_a), (cycles_b, panel_b)):
        if panel.dataset_id != table.dataset_id:
            raise ValueError(
                f"the cycles figure of {panel.dataset_id!r} given for the cycles of "
                f"{table.dataset_id!r}"
            )
    pair = compare_pair(
        cycle_event_table(cycles_a),
        cycle_event_table(cycles_b),
        tau_rule=tau_rule,
        conf_level=conf_level,
        test=False,
    )
    pair = replace(
        pair,
        summary_a=_with_checks(pair.summary_a, panel_a),
        summary_b=_with_checks(pair.summary_b, panel_b),
        checks_a=panel_a.checks,
        checks_b=panel_b.checks,
    )
    return ByThresholdSide(ladder=("",), items={("", CYCLES): pair})


def _with_checks(
    summary: CurveSummaries | None, panel: CyclesPanelData
) -> CurveSummaries | None:
    """The summary carrying the device's check outcome; None, for a device with no cycle."""
    return None if summary is None else replace(summary, checks=panel.checks)
