"""Statistical checks on carved window durations.

Each check answers "is this duration sequence consistent with <null>?" and returns a
`CheckResult` carrying the calibration and the event counts alongside the p-value, because
a p-value from 20 events calibrated asymptotically is a different object from the same
number from 355 events calibrated by permutation.

    C1  Lewis-Robinson trend, time-censored + multi-process   (Kvaloy & Lindqvist eq 4/16)
    C2  Anderson-Darling renewal, time-censored               (Kvaloy & Lindqvist eq 7)
    C3  copula serial independence, via R                     (UNASSESSED - no bench cell)
    C5  studentized max-lag rank autocorrelation, permutation
    C6  portmanteau rank exchangeability, permutation
    CvM Cramer-von Mises renewal, time-censored               (promoted after C1-C6)

Equation numbers are those of Kvaloy & Lindqvist, Technometrics 62(1):101-115 (2020).

Nothing here reads disk (except C3's optional subprocess), imports matplotlib, or prints.
"""

from quebra.analyzers.checks.result import (  # noqa: F401
    CALIB_ASYMPTOTIC,
    CALIB_PERMUTATION,
    CALIB_R_COPULA,
    CLOCK_CALENDAR,
    CLOCK_IN_SPEC,
    CLOCKS,
    CheckResult,
    Segment,
)
