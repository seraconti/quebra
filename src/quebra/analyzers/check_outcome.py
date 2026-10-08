"""The check-outcome artifact, and the two selections that decide what is in it.

Ledger rows in, one typed artifact out, no disk and no matplotlib - the same shape as
`independence_survey.build_independence_survey`, read by one `BasePlot` subclass.

**Two independent selections over the same vocabulary.** The RUN-SET is which check rows are
computed; the DISPLAY-SET is which reach a figure. Running all six and showing one is a normal
thing to want, so they are separate tuples rather than one setting with a display flag. Both
address the full `(check_id, calibration, variant)` key, because C1's asymptotic and
permutation rows can disagree and the asymptotic ones are measured as oversized at the event
counts this project has (`docs/iid_checks/C1_lewis_robinson.md`).

Both are STEP KWARGS. `core/closure.py`'s `_render` takes scalars and containers of scalars,
so a tuple of string triples enters the run identity and reaches the Mermaid label; a typed
object would raise there instead, which is why the selection is data and not a class.

**No verdict stops anything.** A `fail`, an `underpowered`, an absent row and an empty ledger
all produce a grid that renders. It raises only on wiring mistakes - a display-set outside
the run-set, a threshold that is not in the ledger - because those are programming errors
and a silent one draws a figure about the wrong thing.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import pandas as pd

from quebra.analyzers.check_ledger import VERDICT_ABSENT
from quebra.analyzers.checks.result import CALIB_ASYMPTOTIC, CALIB_PERMUTATION
from quebra.analyzers.independence_survey import (
    C3_KEY,
    CHECK_LABELS,
    CHECK_NULL,
    SURVEY_KEYS,
)
from quebra.core._artifact_guard import StaleArtifactGuard


RowKey = tuple[str, str, str]

# Everything selectable: the nine rows the bench scores, plus C3.
ALL_KEYS: tuple[RowKey, ...] = tuple(SURVEY_KEYS)

# The rows `include_tau_checks=False` drops. MEASURED by running the battery both ways on a
# single well-posed segment, not transcribed from the flag's docstring: the three asymptotic
# rows plus C1's and C2's permutation rows, five of nine. CvM's PERMUTATION row survives,
# which is easy to get wrong from prose alone because its asymptotic sibling does not.
TAU_DEPENDENT_KEYS: frozenset[RowKey] = frozenset(
    {
        ("c1_lewis_robinson", CALIB_ASYMPTOTIC, ""),
        ("c1_lewis_robinson", CALIB_PERMUTATION, ""),
        ("c2_anderson_darling", CALIB_ASYMPTOTIC, ""),
        ("c2_anderson_darling", CALIB_PERMUTATION, ""),
        ("cvm_cramer_von_mises", CALIB_ASYMPTOTIC, ""),
    }
)

C2_ASYMPTOTIC_KEY: RowKey = ("c2_anderson_darling", CALIB_ASYMPTOTIC, "")

# A useful default with a reason behind it: at the event counts this project has, the
# asymptotic calibrations are measured as oversized, so a reader shown only the permutation
# rows is shown the rows entitled to the data. Not the default anywhere - a job asks for it.
PERMUTATION_KEYS: tuple[RowKey, ...] = tuple(
    key for key in ALL_KEYS if key[1] == CALIB_PERMUTATION
)


def normalise(keys: Iterable[Sequence[str]], *, field: str) -> tuple[RowKey, ...]:
    """Validate a selection and return it as row-key tuples in `ALL_KEYS` order.

    Order is imposed rather than preserved so that two jobs naming the same rows produce the
    same step kwarg and therefore the same run identity. Duplicates collapse for the same
    reason.
    """
    seen: set[RowKey] = set()
    for raw in keys:
        key = tuple(str(part) for part in raw)
        if len(key) != 3:
            raise ValueError(
                f"{field} entry {raw!r} is not a (check_id, calibration, variant) triple. "
                f"Rows are addressed by the full key because a check's asymptotic and "
                f"permutation rows can disagree."
            )
        if key not in ALL_KEYS:
            raise ValueError(
                f"{field} names {key!r}, which is not a selectable row. "
                f"Known rows: {[list(k) for k in ALL_KEYS]}"
            )
        seen.add(key)  # type: ignore[arg-type]
    return tuple(key for key in ALL_KEYS if key in seen)


def resolve(
    run_set: Iterable[Sequence[str]],
    display_set: Iterable[Sequence[str]],
) -> tuple[tuple[RowKey, ...], tuple[RowKey, ...]]:
    """Normalise both selections and enforce that the display-set is a subset of the run-set.

    Showing a row that was never computed is not a rendering choice, it is a claim about a
    check that did not happen, so it is refused here rather than rendered as a blank.
    """
    run = normalise(run_set, field="run_set")
    display = normalise(display_set, field="display_set")
    if not run:
        raise ValueError("run_set is empty: there would be no check to run at all")
    extra = [key for key in display if key not in run]
    if extra:
        raise ValueError(
            f"display_set names rows outside the run_set: {[list(k) for k in extra]}. "
            f"A row that is not run cannot be shown - it would render as a blank that "
            f"looks like 'no answer' when the truth is 'never asked'."
        )
    return run, display


def battery_flags(run_set: Sequence[RowKey]) -> dict[str, bool]:
    """The battery's own switches, derived from the run-set.

    These are a COMPUTE optimisation and not the contract. They are coarse: asking for C1's
    asymptotic row alone still turns `include_tau_checks` on, which computes C2's rows too.
    `select_rows` is what makes the result match the run-set exactly.
    """
    keys = set(run_set)
    return {
        "include_tau_checks": bool(keys & TAU_DEPENDENT_KEYS),
        "include_c2_asymptotic": C2_ASYMPTOTIC_KEY in keys,
        "include_c3": C3_KEY in keys,
    }


def missing_rows(
    rows: Iterable[object], run_set: Sequence[RowKey], key_of
) -> tuple[RowKey, ...]:
    """Rows the run-set asked for that produced nothing at all.

    A row can vanish legitimately at runtime: C1 and C2 are dropped when `tau` is singular,
    and C2's asymptotic row for m > 1 segments. The caller renders these as `not computed`
    with the reason, because a row that was asked for and answered nothing must not look the
    same as a row that was never asked for.
    """
    present = {key_of(row) for row in rows}
    return tuple(key for key in ALL_KEYS if key in set(run_set) and key not in present)


# `VERDICT_ABSENT` lives in `check_ledger` beside the five the ledger emits, so the verdict
# vocabulary has one home. It is deliberately NOT `not computed`, which means the check ran and
# declined to answer: a reader who cannot tell those apart cannot tell "we asked and got
# nothing" from "we never asked".


@dataclass(frozen=True, slots=True)
class OutcomeGrid:
    """One displayed check, across the datasets drawn and the clocks scored.

    `verdicts[i][j]` is the verdict for `datasets[i]` on `clocks[j]`, `VERDICT_ABSENT` where the
    ledger has no such row. `p_values` and `notes` are aligned with it so the renderer can
    annotate without going back to a frame.

    `datasets` holds the ledger's own `dataset_id`, which is what rows are keyed by;
    `dataset_labels` holds what a reader should see. Looking a row up by its label finds
    nothing and renders as an honest grid of absent cells, so the two are separate fields.
    """

    key: RowKey
    label: str
    null_statement: str
    datasets: tuple[str, ...]
    dataset_labels: tuple[str, ...]
    clocks: tuple[str, ...]
    verdicts: tuple[tuple[str, ...], ...]
    p_values: tuple[tuple[float | None, ...], ...]
    notes: tuple[tuple[str, ...], ...]

    @property
    def counts(self) -> dict[str, int]:
        """How many cells carry each verdict. Named to match `InstrumentGrid.counts`."""
        tally: dict[str, int] = {}
        for row in self.verdicts:
            for verdict in row:
                tally[verdict] = tally.get(verdict, 0) + 1
        return tally


# NOT `slots=True`. `StaleArtifactGuard.__setstate__` restores through `self.__dict__`, which
# a slots dataclass does not have, and every other guard user is a plain dataclass for the
# same reason. Declaring slots here would make the staleness guard raise on every valid
# artifact instead of only on a stale one.
@dataclass
class CheckOutcome(StaleArtifactGuard):
    """Everything the check-outcome figure draws, and what was asked for to produce it.

    `run_set` and `display_set` are ON the artifact, not only in the job file, because a
    figure read six months from now has to say which checks were asked for as well as which
    are shown. `asked_but_unanswered` is the difference that matters most: a row that was run
    and produced nothing is not the same as a row nobody asked for, and only the artifact can
    tell them apart.
    """

    grids: list[OutcomeGrid]
    threshold_label: str = ""
    run_set: tuple[RowKey, ...] = ()
    display_set: tuple[RowKey, ...] = ()
    asked_but_unanswered: tuple[RowKey, ...] = ()
    datasets: tuple[str, ...] = ()
    clocks: tuple[str, ...] = ()
    alpha: float = 0.05

    def check_thresholds(self, labels: list[str]) -> None:
        """Match the artifact-guard contract the other panel artifacts use."""
        if self.threshold_label and labels and self.threshold_label not in labels:
            raise ValueError(
                f"CheckOutcome was built at {self.threshold_label!r}, which is not on the "
                f"requested ladder {labels}"
            )


def _row_key_of(row: pd.Series) -> RowKey:
    return (str(row["check_id"]), str(row["calibration"]), str(row["variant"]))


def build_check_outcome(
    ledger_rows: pd.DataFrame,
    *,
    threshold_label: str,
    run_set: tuple[RowKey, ...],
    display_set: tuple[RowKey, ...],
    datasets: tuple[str, ...],
    clocks: tuple[str, ...],
    dataset_labels: tuple[str, ...] | None = None,
    alpha: float,
) -> CheckOutcome:
    """Turn ledger rows into one grid per displayed key.

    `ledger_rows` is the concatenation of every drawn dataset's `CheckLedger.rows`. Passing
    the frame rather than the ledgers keeps this a pure reshape: it reads columns and writes
    a typed artifact, and computes no statistic of its own.
    """
    run, display = resolve(run_set, display_set)
    labels = tuple(dataset_labels) if dataset_labels is not None else tuple(datasets)
    if len(labels) != len(datasets):
        raise ValueError(
            f"dataset_labels has {len(labels)} entries for {len(datasets)} datasets; they "
            f"are positional, so a mismatch would label a row with another record's name."
        )

    if len(ledger_rows):
        known = set(ledger_rows["threshold_label"].astype(str))
        if threshold_label not in known:
            # A wiring mistake, not a result. Silently yielding an all-absent grid renders
            # as an honest figure about a threshold nobody scored - the same mistake
            # `kaplan_meier.make_inputs_from_windows` raises on, for the same reason.
            raise KeyError(
                f"threshold {threshold_label!r} is not in the ledger. Scored thresholds: "
                f"{sorted(known)}"
            )
        at_threshold = ledger_rows[
            ledger_rows["threshold_label"].astype(str) == threshold_label
        ]
    else:
        at_threshold = ledger_rows

    # Indexed once so a missing cell is a lookup miss rather than a repeated frame scan.
    cells: dict[tuple[str, str, RowKey], pd.Series] = {}
    for _, row in at_threshold.iterrows():
        cells[(str(row["dataset_id"]), str(row["clock"]), _row_key_of(row))] = row

    grids: list[OutcomeGrid] = []
    for key in display:
        verdicts: list[tuple[str, ...]] = []
        p_values: list[tuple[float | None, ...]] = []
        notes: list[tuple[str, ...]] = []
        for dataset in datasets:
            v_row: list[str] = []
            p_row: list[float | None] = []
            n_row: list[str] = []
            for clock in clocks:
                cell = cells.get((dataset, clock, key))
                if cell is None:
                    v_row.append(VERDICT_ABSENT)
                    p_row.append(None)
                    n_row.append("no row in the ledger for this cell")
                    continue
                v_row.append(str(cell["verdict"]))
                raw_p = cell["p_value"]
                p_row.append(None if pd.isna(raw_p) else float(raw_p))
                note = cell["notes"]
                n_row.append("" if pd.isna(note) else str(note))
            verdicts.append(tuple(v_row))
            p_values.append(tuple(p_row))
            notes.append(tuple(n_row))

        grids.append(
            OutcomeGrid(
                key=key,
                label=CHECK_LABELS.get(key, " ".join(p for p in key if p)),
                null_statement=CHECK_NULL.get(key, ""),
                datasets=tuple(datasets),
                dataset_labels=labels,
                clocks=tuple(clocks),
                verdicts=tuple(verdicts),
                p_values=tuple(p_values),
                notes=tuple(notes),
            )
        )

    present = list(at_threshold.itertuples(index=False)) if len(at_threshold) else []
    unanswered = missing_rows(
        present,
        run,
        lambda r: (str(r.check_id), str(r.calibration), str(r.variant)),
    )

    return CheckOutcome(
        grids=grids,
        threshold_label=threshold_label,
        run_set=run,
        display_set=display,
        asked_but_unanswered=unanswered,
        datasets=tuple(datasets),
        clocks=tuple(clocks),
        alpha=alpha,
    )
