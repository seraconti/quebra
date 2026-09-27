"""Run every check on one record off ONE permutation set.

This is the bench's entry point and the reason `_permutation.PermutationSet` is a shared
object rather than each check's private business. Two things follow from computing the
permuted matrix once here:

- **Speed.** Each check that needs permutations gathers a `(B, N)` matrix; at N = 355,
  B = 999 that is a 355k-element copy, and C1, C2, C5-studentized, C5-raw and C6 all want
  the identical one. Sharing it roughly halves a full replicate at N = 355 and the n sweep
  with it.
- **Pairing.** Every check sees the same permutations, so any two rows from one call differ
  only in the statistic, not in the Monte Carlo noise of the reference set. The bench
  leans on that when it compares checks against each other and when it compares Arm B
  against Arm C at rho = 0: a paired difference has a smaller Monte Carlo error than the
  difference of two independently calibrated runs.

C3 is not included. It is an out-of-process R call with its own simulated null, so it
shares nothing with the permutation set and would serialise the whole battery behind a
subprocess. The bench does not call it at all: C3 has no arm, no `ROW_KEYS` entry and no
row in either results table, so its size and power are unmeasured whether or not R is
installed. A caller that wants C3 sets `include_c3` on the `check_ledger` inputs; the
ledger invokes the bridge directly, never through this function.

Nothing is caught. A degenerate record - every duration tied, a zero-length window, a
segment whose last event lands exactly on its truncation time - raises out of here, and the
BENCH decides what to do with it, counting the failure into a `n_failed` column. Swallowing
it at this level would turn a degenerate cell into a reported size.
"""

from __future__ import annotations

import numpy as np

import quebra.analyzers.checks.c1_lewis_robinson as c1
import quebra.analyzers.checks.c2_anderson_darling as c2
import quebra.analyzers.checks.c5_rank_autocorr as c5
import quebra.analyzers.checks.c6_exchangeability as c6
import quebra.analyzers.checks.cvm_cramer_von_mises as cvm
from quebra.analyzers.checks._multiprocess import GAMMA_COMPLETE
from quebra.analyzers.checks._permutation import (
    resolve_perm,
    DEFAULT_N_PERM,
    PermutationSet,
)
from quebra.analyzers.checks._rank_serial import (
    MAX_LAG_CAP,
    autocorrelations,
    global_ranks,
    lag_layout,
)
from quebra.analyzers.checks.result import (
    CALIB_ASYMPTOTIC,
    CALIB_PERMUTATION,
    CLOCK_IN_SPEC,
    CheckResult,
    Segment,
    concatenated_gaps,
    segment_sizes,
)

# One row per (check, calibration, variant). This tuple IS the row schema of the bench
# tables, so adding a check here adds a row everywhere without touching the runner.
ROW_KEYS = (
    ("c1_lewis_robinson", CALIB_ASYMPTOTIC, ""),
    ("c1_lewis_robinson", CALIB_PERMUTATION, ""),
    ("c2_anderson_darling", CALIB_ASYMPTOTIC, ""),
    ("cvm_cramer_von_mises", CALIB_ASYMPTOTIC, ""),
    ("c2_anderson_darling", CALIB_PERMUTATION, ""),
    ("c5_rank_autocorr", CALIB_PERMUTATION, c5.VARIANT_STUDENTIZED),
    ("c5_rank_autocorr", CALIB_PERMUTATION, c5.VARIANT_RAW),
    ("c6_exchangeability", CALIB_PERMUTATION, ""),
    # CvM is the fourth functional of the same Brownian bridge as
    # C1 and C2, and it is here for a reason the other two cannot cover: eq (7) carries a
    # `1/(s(1-s))` weight, so C2 is singular when the last event lands on the truncation
    # time, and on the IN-SPEC clock of a carved record that is the common case - measured,
    # 311 of 340 survey cells have no C1/C2 answer for exactly that reason. CvM's integrand
    # has no such weight and is finite there. It is also the source paper's own preference
    # for m > 1, which is the gapped case.
    ("cvm_cramer_von_mises", CALIB_PERMUTATION, ""),
)


def row_key(result: CheckResult) -> tuple[str, str, str]:
    """The `ROW_KEYS` entry a result belongs to, read back off the result itself."""
    variant = ""
    for token in result.notes.split():
        if token.startswith("variant="):
            variant = token.split("=", 1)[1]
    return (result.check, result.calibration, variant)


def run_battery(
    segments: list[Segment],
    *,
    clock: str = CLOCK_IN_SPEC,
    perm: PermutationSet | None = None,
    n_perm: int = DEFAULT_N_PERM,
    rng: np.random.Generator | None = None,
    gamma_estimator: str = GAMMA_COMPLETE,
    max_lag: int = MAX_LAG_CAP,
    include_c2_asymptotic: bool = True,
    include_tau_checks: bool = True,
) -> list[CheckResult]:
    """Up to nine results in `ROW_KEYS` order.

    `include_c2_asymptotic` exists because C2's asymptotic calibration is defined only for
    a single segment (Kvaloy & Lindqvist Section 4.2 reject the normal approximation for
    the summed statistic). The battery drops that row for m > 1 rather than raising, since
    a gapped record legitimately has m > 1 and the other eight rows are still wanted.

    `include_tau_checks=False` drops five rows - the three asymptotic ones and C1's and C2's
    permutation rows - leaving four. It is for records where `tau` is not well posed - specifically the IN-SPEC clock of a
    carved series, where in-spec time stops accumulating the moment the record ends out of
    spec, so `tau == T_N` and eq (7) is singular. Measured on an iid read series that is
    73% of replicates, and the event triggering it depends on the data, so running C1/C2
    on the surviving 27% would report a size conditioned on how the record happened to
    end. C5 and C6 never touch `tau` and stay valid there, which is why they are kept
    rather than the whole clock being dropped.
    """
    if not segments:
        raise ValueError("the battery needs at least one segment")
    sizes = segment_sizes(segments)

    if perm is None and rng is None:
        # Raised here rather than three frames down, because this is the call a pipeline
        # step makes and the fix belongs at the job level: declare an integer seed beside
        # alpha and the ladder, pass it as a step kwarg so it reaches the provenance
        # label, and build the generator from it.
        raise ValueError(
            "run_battery needs either a prebuilt `perm` or an `rng`. Six of its nine "
            "rows are permutation-calibrated, so without one the p-values are a fresh "
            "random draw on every call while the run identity stays unchanged."
        )
    perm = resolve_perm(sizes, perm, n_perm, rng)

    # The two shared intermediates. Everything below is an aggregation of one of these.
    # The gap gather is skipped when C1/C2 are off, since nothing else reads it.
    gaps = concatenated_gaps(segments)
    # CvM wants the permuted matrix as well, and unlike C1/C2 it runs even when
    # `include_tau_checks` is False - so this is not gated on that flag alone.
    permuted = perm.apply(gaps)
    layout = lag_layout(sizes, max_lag)
    ranks = global_ranks(gaps)
    observed_r = autocorrelations(ranks, layout)[0]
    null_r = autocorrelations(perm.apply(ranks), layout)
    autocorr = (observed_r, null_r, layout)

    results: list[CheckResult] = []
    if include_tau_checks:
        results.extend(
            [
                c1.run(
                    segments,
                    calibration=CALIB_ASYMPTOTIC,
                    clock=clock,
                    gamma_estimator=gamma_estimator,
                ),
                c1.run(
                    segments,
                    calibration=CALIB_PERMUTATION,
                    clock=clock,
                    gamma_estimator=gamma_estimator,
                    perm=perm,
                    permuted=permuted,
                ),
            ]
        )
        if include_c2_asymptotic and len(segments) == 1:
            results.append(
                c2.run(
                    segments,
                    calibration=CALIB_ASYMPTOTIC,
                    clock=clock,
                    gamma_estimator=gamma_estimator,
                )
            )
        # CvM's ASYMPTOTIC row belongs inside this gate with C1 and C2. Its STATISTIC is
        # finite when `tau == T_N` - there is no `1/(s(1-s))` weight to blow up - but its
        # limiting null still assumes a truncation time chosen INDEPENDENTLY of the events,
        # and an event-determined `tau` breaks the tied-down bridge for CvM exactly as it
        # does for the other two. Only the permutation row below is entitled to that case.
        if len(segments) == 1:
            results.append(
                cvm.run(
                    segments,
                    calibration=CALIB_ASYMPTOTIC,
                    clock=clock,
                    gamma_estimator=gamma_estimator,
                )
            )
        results.append(
            c2.run(
                segments,
                calibration=CALIB_PERMUTATION,
                clock=clock,
                gamma_estimator=gamma_estimator,
                perm=perm,
                permuted=permuted,
            )
        )
    results.extend(
        [
            c5.run(
                segments,
                variant=c5.VARIANT_STUDENTIZED,
                clock=clock,
                autocorr=autocorr,
            ),
            c5.run(segments, variant=c5.VARIANT_RAW, clock=clock, autocorr=autocorr),
            c6.run(segments, clock=clock, autocorr=autocorr),
        ]
    )
    # CvM by PERMUTATION runs always, including when the tau checks are off. That is the
    # point of promoting it: on the in-spec clock of a carved record `tau == T_N` silences
    # C1 and C2 - measured, 311 of 340 survey cells - and the permutation calibration
    # conditions on the observed data, so it is entitled to that case while their
    # asymptotic routes are not.
    results.append(
        cvm.run(
            segments,
            calibration=CALIB_PERMUTATION,
            clock=clock,
            gamma_estimator=gamma_estimator,
            perm=perm,
            permuted=permuted,
        )
    )
    return results
