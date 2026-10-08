"""What an estimate rests on, and what was checked: the fields every result carries.

A leaf: it imports nothing from quebra (`tests/test_check_attachment.py` holds it to
that). Estimators and bands carry these fields, and importing the check layer from them
would put every check into the identity of every job that draws a curve. So the
assumption id is a string literal (the same test file pins it to
`analyzers/assumptions.py`), and the verdicts are plain tuples of strings.

Every number a survival step reports carries one standing:

- `SURVIVES`: a true statement about this record whatever the checks say (the in-spec
  fraction, the mean window length, the placement bracket).
- `DESCRIPTIVE`: an exact summary of this record when no window is censored. With censored
  windows it leans on the multiplicative intensity model, one hazard shared by every window
  at risk and independent censoring, so `n_censored` travels beside it (Kaplan-Meier point
  values, quantile points, the restricted mean).
- `RESTS_ON_EXCHANGEABILITY`: wrong if the windows are not draws from one law; it carries
  the check outcome (every band, CI and SE, the Nelson-Aalen shape, Turnbull).

An empty attachment means NOT ASSESSED, which is a different claim from "assessed and
nothing rejected". A check annotates; it never changes what is computed or drawn.
"""

from __future__ import annotations

from dataclasses import dataclass

A1_RENEWAL_DURATIONS_ID = "a1_renewal_durations"

SURVIVES = "survives a rejection"
DESCRIPTIVE = "descriptive"
RESTS_ON_EXCHANGEABILITY = "rests on exchangeability"


def summarise(
    checks_asked: tuple[str, ...],
    checks_unanswered: tuple[str, ...],
    check_verdicts: tuple[tuple[str, str, str], ...],
) -> str:
    """One line naming what was asked and what came back, for a caption to state."""
    if not checks_asked:
        return "independence checks: NOT ASSESSED for this band"
    tally: dict[str, int] = {}
    for _label, _dataset, verdict in check_verdicts:
        tally[verdict] = tally.get(verdict, 0) + 1
    shown = ", ".join(f"{n} {v}" for v, n in sorted(tally.items()))
    line = f"{len(checks_asked)} checks asked; cells: {shown or 'none'}"
    if checks_unanswered:
        line += f"; no answer anywhere from {', '.join(checks_unanswered)}"
    return line


@dataclass(frozen=True)
class CheckAttachment:
    """The assumption a result rests on, the checks asked for it, and what came back.

    `check_verdicts` holds (check label, dataset label, verdict) triples.
    """

    assumption_id: str = A1_RENEWAL_DURATIONS_ID
    checks_asked: tuple[str, ...] = ()
    checks_unanswered: tuple[str, ...] = ()
    check_verdicts: tuple[tuple[str, str, str], ...] = ()

    def summary(self) -> str:
        return summarise(self.checks_asked, self.checks_unanswered, self.check_verdicts)
