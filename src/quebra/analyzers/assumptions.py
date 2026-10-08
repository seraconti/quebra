"""The validity conditions the reliability arithmetic runs on, as data.

An assumption record carries the four fields `AGENTS.md` section 10 already mandates for a
statistical docstring - assumption, diagnostic, consequence of violation, reference - plus an
id, the disposition, and the checks that diagnose it. Written as records rather than as prose
so a figure can name them, a test can cite them, and one grep can find every place an
assumption is relied on.

Ids follow the repo's existing lowercase convention (`a1_renewal_durations`), the same shape
as `CHECK_NAME` strings and the `docs/iid_checks/C1_*.md` pages, so the trace query spans
code, tests and docs. This declines `spec/quebraplan.md` 4.4's uppercase `test_A3_*` style on
purpose.

**Adding a record is one literal and one line.** Write the `Assumption`, add it to `_RECORDS`.
Everything else - the registry, the id-shape check, the trace guard in
`tests/test_assumptions.py` - follows from that, and `diagnostic_checks` lives ON the record
so there is no second registry to keep in sync. That is deliberate: more of these are
expected, and one that is painful to add is one that does not get added.

**What this module is not.** It records validity conditions; it does not test them and does
not enforce them. No band is suppressed here: the band is always computed and always drawn,
with its check outcome attached (`AGENTS.md` section 5).

Pure data. The only import is the checks' own names, so a record cannot name a diagnostic
that does not ship.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from quebra.analyzers.checks import c1_lewis_robinson as _c1
from quebra.analyzers.checks import c2_anderson_darling as _c2
from quebra.analyzers.checks import c3_serial_copula as _c3
from quebra.analyzers.checks import c5_rank_autocorr as _c5
from quebra.analyzers.checks import c6_exchangeability as _c6
from quebra.analyzers.checks import cvm_cramer_von_mises as _cvm

SHIPPED_CHECK_NAMES: frozenset[str] = frozenset(
    m.CHECK_NAME for m in (_c1, _c2, _c3, _c5, _c6, _cvm)
)

# `a<number>_<snake_case>`. Enforced so a fifth record cannot arrive as `A5-Foo` and quietly
# fall outside the grep the whole traceability scheme rests on.
_ID_SHAPE = re.compile(r"^a[0-9]+_[a-z0-9]+(?:_[a-z0-9]+)*$")

# What the code does when the assumption fails. Closed: a record whose disposition is not one
# of these is claiming a behaviour the package does not have.
#
# `undetected` has no members yet and is kept on purpose: it is the honest value for an
# assumption the package cannot see fail, so adding such a record does not also require
# widening this vocabulary.
#
# `warns` is absent for the opposite reason: `warnings.warn` occurs once in the package
# (`schemas/ramsey_series.py`, clock provenance) and never for a statistical assumption, so it
# describes a behaviour that does not exist here. Add it when something needs it.
DISPOSITION_RAISES = "raises"
DISPOSITION_REPORTS = "reports"
DISPOSITION_UNDETECTED = "undetected"
DISPOSITIONS = (DISPOSITION_RAISES, DISPOSITION_REPORTS, DISPOSITION_UNDETECTED)


@dataclass(frozen=True, slots=True)
class Assumption:
    """One validity condition, and what happens when it does not hold.

    Frozen because a record is a statement of fact about the package, not configuration: a
    caller that could edit one could make a figure claim a diagnostic that never ran.

    `diagnostic_checks` is the machine-readable half of `diagnostic`, by `CHECK_NAME`. An
    EMPTY tuple is a claim, not an omission: it says no shipped check diagnoses this, and the
    disposition must then be `undetected`.
    """

    id: str
    statement: str
    diagnostic: str
    consequence: str
    reference: str
    disposition: str
    diagnostic_checks: tuple[str, ...] = field(default=())

    def __post_init__(self) -> None:
        for name in ("id", "statement", "diagnostic", "consequence", "reference"):
            if not str(getattr(self, name)).strip():
                raise ValueError(
                    f"assumption {self.id!r} has an empty {name!r}. Every field is "
                    f"load-bearing: an empty diagnostic reads as 'not checked' and an empty "
                    f"reference as 'no source', and neither should be inferred from a blank."
                )
        if not _ID_SHAPE.match(self.id):
            raise ValueError(
                f"assumption id {self.id!r} is not of the form a<number>_<snake_case>. "
                f"The traceability scheme is a grep over this shape, so an id outside it "
                f"is invisible to it."
            )
        if self.disposition not in DISPOSITIONS:
            raise ValueError(
                f"assumption {self.id!r} has disposition {self.disposition!r}; "
                f"known: {list(DISPOSITIONS)}"
            )
        unknown = set(self.diagnostic_checks) - SHIPPED_CHECK_NAMES
        if unknown:
            raise ValueError(
                f"assumption {self.id!r} names checks that do not ship: {sorted(unknown)}. "
                f"Shipped: {sorted(SHIPPED_CHECK_NAMES)}"
            )
        if not self.diagnostic_checks and self.disposition != DISPOSITION_UNDETECTED:
            raise ValueError(
                f"assumption {self.id!r} names no diagnosing check but claims disposition "
                f"{self.disposition!r}. With nothing to run, the package cannot report or "
                f"raise anything, so the honest value is {DISPOSITION_UNDETECTED!r}."
            )


A1_RENEWAL_DURATIONS = Assumption(
    id="a1_renewal_durations",
    statement=(
        "The window durations on the declared side behave like a renewal process: "
        "independent, identically distributed, no trend. Reading a survival curve as the "
        "law of a window's duration rests on it, as do its intervals and any reading of "
        "occupancy or mean time between failures beyond this record."
    ),
    diagnostic=(
        "The six checks in `analyzers/checks/`, scored into verdicts by "
        "`analyzers/check_ledger.py`. Read the diagnostic's SAMPLE before transferring its "
        "verdict: `segments_from_windows` on the in-spec clock uses complete-window "
        "durations and folds the censored window into tau, and on the calendar clock the "
        "events are window BIRTHS, while `kaplan_meier.make_inputs_from_windows` keeps "
        "observed-birth windows and retains the censored ones as censored observations. A "
        "verdict about one of those samples is not automatically a verdict about the other, "
        "and whether the transfer is valid is UNRESOLVED."
    ),
    consequence=(
        "The arithmetic still produces numbers and the numbers are wrong in ways nothing "
        "else catches. A trend or serial dependence in the durations biases the survival "
        "curve and every summary read off it, without any symptom at the call site."
    ),
    reference=(
        "Kvaloy & Lindqvist, Technometrics 62(1):101-115 (2020), for the "
        "trend-against-renewal tests this "
        "package implements. Ascher & Feingold for the renewal-process framing itself; no "
        "locator located."
    ),
    disposition=DISPOSITION_REPORTS,
    diagnostic_checks=(
        _c1.CHECK_NAME,
        _c2.CHECK_NAME,
        _c3.CHECK_NAME,
        _c5.CHECK_NAME,
        _c6.CHECK_NAME,
        _cvm.CHECK_NAME,
    ),
)


# Add a record here and it is registered, id-checked and covered by the trace guard. Order is
# the order a reader should meet them.
_RECORDS: tuple[Assumption, ...] = (A1_RENEWAL_DURATIONS,)

ASSUMPTIONS: dict[str, Assumption] = {record.id: record for record in _RECORDS}


__all__ = [
    "ASSUMPTIONS",
    "DISPOSITIONS",
    "DISPOSITION_RAISES",
    "DISPOSITION_REPORTS",
    "DISPOSITION_UNDETECTED",
    "SHIPPED_CHECK_NAMES",
    "A1_RENEWAL_DURATIONS",
    "Assumption",
]
