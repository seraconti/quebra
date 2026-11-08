"""The calibration-cycle jobs: one Quantum Inspire status history, its cycles and figure.

Per job: load the status history (`schemas/qi_status_history`), turn its samples into
spells, the spells into cycles and status shares, run the check battery on the cycles
against the bench's size table, build `CyclesPanelData`, materialize it and draw
`CalibrationCyclesPanel`. The cycles are materialized too, so a composite can compare two
devices on them (`analyzers/cycle_comparison.py`). The spells also give the two context
tables, materialized and not drawn: one row per cycle (`cycle_context`) and one per period
of observed exposure (`cycle_periods`). Every step is a function of `analyzers/cycles.py`
called as it stands, so a node's identity is that module's code and the knobs on its label.

A module of its own, not part of `quebra.recipes`, for the reason `ledger_recipe.py` gives:
whatever a module imports joins the identity closure of every step function it defines, and
the cycles must not re-key the T2* jobs, nor the T2* code the cycle jobs.
"""

from __future__ import annotations

from quebra.analyzers.check_outcome import PERMUTATION_KEYS
from quebra.analyzers.cycles import (
    N_PERIODS,
    build_cycles_panel_data,
    cycle_checks,
    cycle_context,
    cycle_periods,
    cycles,
    spells_from_samples,
    status_shares,
)
from quebra.analyzers.kaplan_meier import CONF_LEVEL
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.panels.calibration_cycles import CalibrationCyclesPanel

# The knobs that turn a p-value into a verdict, with the values `ledger_recipe.LEDGER_KNOBS`
# gives the windows, so a row on cycles is scored on the rule a row on windows is. The seed
# is each job's own. Each knob reaches the provenance label as a step kwarg.
CYCLE_KNOBS: dict[str, object] = {
    # The permutation rows: the a2 record's diagnostic (`assumptions.A2_RENEWAL_CYCLES`).
    "run_set": PERMUTATION_KEYS,
    "alpha": 0.05,
    "min_events_pass": 35,
    "tie_cutoff_distinct": 5,
    "lag_max": 5,
    "n_permutations": 999,
    # `check_ledger.CheckLedgerInputs`' default for windows.
    "min_events_per_segment": 2,
    "conf_level": CONF_LEVEL,
}

TARGETS = ["static", "academic"]


def configure_qi_cycles_job(
    job: Job,
    *,
    prefix: str,
    dataset: Dataset,
    bench_size_table: Dataset,
    spacing_s: float,
    seed: int,
    run_set: tuple[tuple[str, str, str], ...],
    alpha: float,
    min_events_pass: int,
    tie_cutoff_distinct: int,
    lag_max: int,
    n_permutations: int,
    min_events_per_segment: int,
    conf_level: float,
) -> None:
    """Wire one status history from load to figure. `prefix` names the sinks:
    `{prefix}_cycles`, `{prefix}_panel_data`, `{prefix}_cycle_context` and
    `{prefix}_cycle_periods` (materialized) and the figure titled `{prefix}`."""
    samples = job.load(dataset)
    spells = job.step(spells_from_samples, samples, name="spells", spacing_s=spacing_s)
    table = job.step(cycles, spells, name="cycles")
    job.materialize(table, name=f"{prefix}_cycles")
    shares = job.step(status_shares, spells, name="status_shares")
    bench = job.load_df(bench_size_table)
    checks = job.step(
        cycle_checks,
        table,
        bench,
        name="cycle_checks",
        run_set=run_set,
        alpha=alpha,
        min_events_pass=min_events_pass,
        tie_cutoff_distinct=tie_cutoff_distinct,
        lag_max=lag_max,
        n_permutations=n_permutations,
        seed=seed,
        min_events_per_segment=min_events_per_segment,
    )
    panel = job.step(
        build_cycles_panel_data,
        table,
        checks,
        shares,
        name="cycles_panel_data",
        conf_level=conf_level,
    )
    # Not `{prefix}`: that is the figure's title, and two sinks may not share a
    # provenance name.
    job.materialize(panel, name=f"{prefix}_panel_data")
    context = job.step(cycle_context, spells, name="cycle_context")
    periods = job.step(cycle_periods, spells, name="cycle_periods", n_periods=N_PERIODS)
    job.materialize(context, name=f"{prefix}_cycle_context")
    job.materialize(periods, name=f"{prefix}_cycle_periods")
    job.figure(CalibrationCyclesPanel, panel, targets=TARGETS, title=prefix)
