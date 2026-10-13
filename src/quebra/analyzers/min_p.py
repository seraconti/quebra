"""Single-step min-P: one calibrated p-value for a panel of statistics on shared permutations.

Each member of a panel is a statistic computed on the observed data and on the SAME B
permutations of it, oriented so that a larger value is more extreme: a two-sided statistic
enters as its absolute value, as the checks' `two_sided_p_value` takes it, and a signed one
entered as it is gives a one-sided member p. A member's own p is the
permutation p the checks report (`checks._permutation.permutation_p_value`): the share of the
B + 1 rows, the observed one included, whose statistic is at least the observed one, ties
counted against. Applying that rule to every row in turn gives each row a p per member; the
row's smallest is its min-p. The panel's p is the share of rows whose min-p is at most the
observed row's.

Because every member is computed on the same permutation rows, the null distribution of the
min-p keeps whatever dependence the members have on each other. Members that rank every row
the same way (an increasing copy) cost nothing over one member; independent members cost
about 1 - (1 - p)^K for K members; no dependence costs more than K p, the Bonferroni cap.

Nothing here decides anything: the result carries every member's p, the panel's p and the
members that attain the observed min-p, for a reader to read.

A leaf: it imports nothing from quebra.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class MinPResult:
    member_p_values: dict[str, float]
    p_value: float
    # Every member whose own p equals the observed row's min-p; ties keep them all.
    argmin_labels: tuple[str, ...]
    n_permutations: int


def row_p_values(statistics: np.ndarray) -> np.ndarray:
    """p of every row of every column: #{rows with a statistic >= this one} / n_rows.

    `statistics` has one column per member and B + 1 rows, row 0 the observed one. Row 0's
    p is the checks' `(1 + #{null >= observed}) / (1 + B)`.
    """
    n_rows = statistics.shape[0]
    p = np.empty_like(statistics, dtype=float)
    for j in range(statistics.shape[1]):
        column = statistics[:, j]
        ordered = np.sort(column)
        # count of values >= x is n_rows minus the count strictly below x
        p[:, j] = (n_rows - np.searchsorted(ordered, column, side="left")) / n_rows
    return p


def min_p(statistics: Mapping[str, Sequence[float] | np.ndarray]) -> MinPResult:
    """The single-step min-P p-value of a panel, with every member's own p.

    `statistics` maps each member's label to its B + 1 values, the observed one first and
    then one per shared permutation, larger meaning more extreme. Raises ValueError on an
    empty panel, a member that is not one-dimensional, unequal lengths, fewer than one
    permutation, or a non-finite value (NaN or infinite): a statistic that a permutation
    makes undefined is a bug in the statistic.

    The panel's p always lies between the smallest member p and the member count times it.
    Ties are counted against, as in the checks. Each member's p lives on the grid
    {1, ..., B + 1} / (B + 1), so the row min-ps tie across rows, and the panel's null size
    at level alpha is at most floor(alpha (B + 1)) / (B + 1) and can be well below it: the
    test is conservative, not exact. The panel's p is also never below d / (B + 1), with d,
    at most the member count, the number of rows that hold some member's largest value
    alone, since each such row has min-p 1 / (B + 1). With K members, B + 1 must be well
    above K / alpha for the panel to reach alpha at all.

    Validity assumptions:

    - Assumption: under the null the B + 1 rows are exchangeable, which holds when every
      member is computed on the same permutations of data whose permutation leaves its
      distribution unchanged.
      Diagnostic: none here; the members' own construction carries it.
      Consequence of violation: neither the members' p-values nor the panel's has its
      stated size. Members computed on DIFFERENT permutations also break the joint row
      structure, so the panel's p no longer reflects how the members depend on each other.
      Reference: Westfall, P. H. and Young, S. S. (1993), Resampling-Based Multiple Testing,
      Wiley, the single-step min-P adjustment; no section locator located.
    """
    if not statistics:
        raise ValueError("min_p needs at least one member statistic")
    labels = list(statistics)
    columns = [np.asarray(statistics[label], dtype=float) for label in labels]
    lengths = {label: col.shape for label, col in zip(labels, columns)}
    if len({shape for shape in lengths.values()}) != 1 or columns[0].ndim != 1:
        raise ValueError(
            f"every member needs one value per row, the same number of rows: {lengths}"
        )
    n_rows = columns[0].size
    if n_rows < 2:
        raise ValueError(
            f"min_p needs the observed row and at least one permutation; got {n_rows} row(s)"
        )
    table = np.column_stack(columns)
    bad = [label for label, col in zip(labels, columns) if not np.all(np.isfinite(col))]
    if bad:
        raise ValueError(f"non-finite statistic(s) in member(s) {bad}")

    p = row_p_values(table)
    row_min = p.min(axis=1)
    observed_min = float(row_min[0])
    return MinPResult(
        member_p_values={label: float(p[0, j]) for j, label in enumerate(labels)},
        p_value=float(np.count_nonzero(row_min <= observed_min) / n_rows),
        argmin_labels=tuple(
            label for j, label in enumerate(labels) if p[0, j] == observed_min
        ),
        n_permutations=n_rows - 1,
    )
