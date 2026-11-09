"""Each compared record's check outcome and group-row verdicts, read off its own ledger.

`pairs_with_checks` fills a pair set and `forest_with_checks` a random-effects set. For each
record, at each threshold and side, they fill two things.
- The band outcome, on the record's check-outcome fields: `checks_a` / `checks_b` and each
  summary of a pair, and each row's summary in a forest. It is the ledger's permutation rows
  of that rung and bag on both clocks, as the record figures read a band's outcome
  (`record_figure._band_outcomes`), naming a1.
- Its group-row verdicts at that rung, (group, bag, clock, verdict), over every bag and
  clock. A ledger's `(all)` row, written where nothing was scored, gives both groups its
  verdict, as `record_figure.expand_all_rows` gives every check of a run-set.

Each carve is held to the ledger that scored it (`record_figure.cross_check_carve`): the
same record, and as many windows at every rung and bag as the ledger rows used. A ledger
that recorded a refusal of its record (`_refused`) scored no window; it is held to its
record only, and its cells read `not computed`.

A forest summary also gets the one line its figures print (`check_line`): how many records
carry checks and which have a `fail` cell or a failing group row, every cell being in the
table. Listing every record's outcome would run to a line per record.

Nothing here computes a number the ledger did not print, and no number of the comparison
changes: the checks annotate it. A forest record with no window at a threshold and side has
no row there, so nothing is attached to it; its table row names the reason.

A module of its own so the check layer stays out of the identity closure of every
comparison that reads no ledger.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import replace

from quebra.analyzers.check_attachment import A1_RENEWAL_DURATIONS_ID, CheckAttachment
from quebra.analyzers.check_ledger import (
    GROUPS,
    VERDICT_FAIL,
    VERDICT_NOT_COMPUTED,
    CheckLedger,
)
from quebra.analyzers.check_outcome import RowKey, check_attachment_of
from quebra.analyzers.record_comparison import (
    GroupVerdicts,
    PairComparison,
    RandomEffectsSet,
)
from quebra.analyzers.record_figure import ALL_ROWS, _band_outcomes, cross_check_carve
from quebra.analyzers.survival import ByThresholdSide
from quebra.analyzers.windows import WindowsResult


def _refused(ledger: CheckLedger) -> bool:
    """Whether `ledger` is a record's recorded refusal, `ledger_recipe.ledger_step` on a
    `RecordFailure`: every row an `(all)` row, `not computed`, using no window, with one
    note shared by every row, `<stage> failed: <reason>` (`core.types.RecordFailure.note`).
    A ledger whose checks failed one by one scored windows, so it is not one."""
    rows = ledger.rows
    if not len(rows):
        return False
    notes = set(rows["notes"].astype(str))
    return (
        bool((rows["check_id"] == ALL_ROWS).all())
        and bool((rows["verdict"] == VERDICT_NOT_COMPUTED).all())
        and bool((rows["n_windows_used"] == 0).all())
        and len(notes) == 1
        and re.match(r"^\S+ failed: ", next(iter(notes))) is not None
    )


def hold_carve_to_ledger(carve: WindowsResult, ledger: CheckLedger) -> None:
    """Raise unless `ledger` scored `carve` (module docstring)."""
    record = str(carve.meta["dataset_id"])
    if record != ledger.dataset_id:
        raise ValueError(
            f"the ledger of {ledger.dataset_id!r} given for the carve of {record!r}"
        )
    if not _refused(ledger):
        cross_check_carve(carve, ledger)


def group_verdicts(ledger: CheckLedger, threshold_label: str) -> GroupVerdicts:
    """Every group-row verdict of `ledger` at `threshold_label`, in row order: (group, bag,
    clock, verdict), an `(all)` row giving every group its verdict. Raises when the ledger
    has no row at that threshold."""
    rows = ledger.rows[ledger.rows["threshold_label"] == threshold_label]
    if not len(rows):
        raise ValueError(
            f"the ledger of {ledger.dataset_id!r} has no row at {threshold_label!r}"
        )
    verdicts: list[tuple[str, str, str, str]] = []
    for row in rows.itertuples(index=False):
        if row.check_id in GROUPS:
            verdicts.append((row.check_id, row.side, row.clock, row.verdict))
        elif row.check_id == ALL_ROWS:
            verdicts += [(group, row.side, row.clock, row.verdict) for group in GROUPS]
    return tuple(verdicts)


def _attachments(
    ledger: CheckLedger, ladder: Sequence[str], run_set: Sequence[RowKey]
) -> dict[tuple[str, str], CheckAttachment]:
    run = tuple(tuple(key) for key in run_set)
    return {
        key: check_attachment_of(outcome, A1_RENEWAL_DURATIONS_ID)
        for key, outcome in _band_outcomes(ledger, tuple(ladder), run).items()  # type: ignore[arg-type]
    }


def pairs_with_checks(
    pairs: ByThresholdSide[PairComparison],
    carve_a: WindowsResult,
    carve_b: WindowsResult,
    ledger_a: CheckLedger,
    ledger_b: CheckLedger,
    *,
    run_set: Sequence[RowKey],
) -> ByThresholdSide[PairComparison]:
    """`pairs` with each record's check outcome and group-row verdicts filled from its own
    ledger, at every threshold and side (module docstring). `run_set` is the checks a band
    reads: the ledger's permutation rows (`check_outcome.PERMUTATION_KEYS`)."""
    for carve, ledger in ((carve_a, ledger_a), (carve_b, ledger_b)):
        hold_carve_to_ledger(carve, ledger)
    checks_a = _attachments(ledger_a, pairs.ladder, run_set)
    checks_b = _attachments(ledger_b, pairs.ladder, run_set)
    items: dict[tuple[str, str], PairComparison] = {}
    for key, pair in pairs.items.items():
        if (pair.record_a, pair.record_b) != (ledger_a.dataset_id, ledger_b.dataset_id):
            raise ValueError(
                f"ledgers of {ledger_a.dataset_id!r} and {ledger_b.dataset_id!r} given for "
                f"the pair {pair.record_a} vs {pair.record_b}"
            )
        label = key[0]
        items[key] = replace(
            pair,
            summary_a=(
                None
                if pair.summary_a is None
                else replace(pair.summary_a, checks=checks_a[key])
            ),
            summary_b=(
                None
                if pair.summary_b is None
                else replace(pair.summary_b, checks=checks_b[key])
            ),
            checks_a=checks_a[key],
            checks_b=checks_b[key],
            group_verdicts_a=group_verdicts(ledger_a, label),
            group_verdicts_b=group_verdicts(ledger_b, label),
        )
    return ByThresholdSide(ladder=pairs.ladder, items=items)


def forest_with_checks(
    forest: RandomEffectsSet,
    *carves_then_ledgers: WindowsResult | CheckLedger,
    records: Sequence[str],
    run_set: Sequence[RowKey],
) -> RandomEffectsSet:
    """`forest` with each record row's check outcome and group-row verdicts filled from its
    own ledger (module docstring).

    `carves_then_ledgers` is every record's carve, in `records` order, then every record's
    ledger in the same order; `records` must name every record of the forest's groups.
    """
    n = len(records)
    if len(carves_then_ledgers) != 2 * n:
        raise ValueError(
            f"{len(carves_then_ledgers)} inputs for {n} records: a carve and a ledger each"
        )
    named = sorted(r for names in forest.groups.values() for r in names)
    if sorted(records) != named:
        raise ValueError(f"records {sorted(records)} differ from the forest's {named}")
    carves = carves_then_ledgers[:n]
    ledgers = carves_then_ledgers[n:]
    by_record: dict[
        str, tuple[dict[tuple[str, str], CheckAttachment], CheckLedger]
    ] = {}
    for record, carve, ledger in zip(records, carves, ledgers, strict=True):
        if ledger.dataset_id != record:  # type: ignore[union-attr]
            raise ValueError(
                f"the ledger of {ledger.dataset_id!r} given for {record!r}"  # type: ignore[union-attr]
            )
        hold_carve_to_ledger(carve, ledger)  # type: ignore[arg-type]
        by_record[record] = (
            _attachments(ledger, forest.ladder, run_set),  # type: ignore[arg-type]
            ledger,  # type: ignore[assignment]
        )
    items = {}
    for (group, label, side), summary in forest.items.items():
        rows = []
        for row in summary.rows:
            attachments, ledger = by_record[row.record]
            rows.append(
                replace(
                    row,
                    summary=replace(row.summary, checks=attachments[(label, side)]),
                    group_verdicts=group_verdicts(ledger, label),
                )
            )
        items[(group, label, side)] = replace(
            summary, rows=tuple(rows), check_line=_check_line(rows)
        )
    return replace(forest, items=items)


def _check_line(rows: Sequence) -> str:
    """The forest summary's figure line (module docstring)."""
    failing = [
        row.record
        for row in rows
        if any(v == VERDICT_FAIL for *_, v in row.summary.checks.check_verdicts)
    ]
    group_failing = [
        row.record
        for row in rows
        if any(v == VERDICT_FAIL for *_, v in row.group_verdicts)
    ]
    return (
        f"{len(rows)} records checked, every cell in the table; a fail cell: "
        f"{', '.join(failing) or 'none'}; a failing group row: "
        f"{', '.join(group_failing) or 'none'}"
    )
