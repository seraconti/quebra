"""Tier 4: our estimators against independent R implementations.

Tier 2 pins arithmetic against a published table. Tier 3 measures whether a p-value holds
its level. Neither catches a routine that is SELF-CONSISTENTLY WRONG - one that computes
some well-behaved statistic which is not the statistic its source defines. A second
implementation, written by other people from the same paper, is what catches that.

For two of these the R package is by authors of the method: `XICOR` is by Holmes and
Chatterjee, `energy` by Rizzo and Szekely. That is as close to a definitional reference as exists.

**Most of this file never runs R.** `rscripts/reference_values.R` writes two committed CSV fixtures
and this file reads them. `pytest` therefore works on a machine with no R, which was the
may be the state of a reviewer's machine.

**The fixture carries the INPUT DATA, not just the answers.** Reproducing R's RNG stream
from Python would be fragile, and a mismatch in the data generator would surface here as a
statistical disagreement - the most misleading possible failure. Python reads the exact
numbers R saw.

**The headline case is `xi_tied`.** `chatterjee_xi` uses the tie-corrected form (Chatterjee 2021), not the
tie-free reduction, and no shipped window exercises it: all 309 window-rows on the T2* ladder
are tie-free. This file is the ONLY external evidence that the tie-corrected path is right,
and `XICOR::xicor(ties = TRUE)` is the reference
implementation of the same equation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.stats import rankdata, spearmanr

from quebra.analyzers.shape_stats import chatterjee_xi, dcor
from tests.fixtures import R_REFERENCE_INPUTS, R_REFERENCE_VALUES

pytestmark = pytest.mark.statistical


_INPUTS = pd.read_csv(R_REFERENCE_INPUTS)
_VALUES = pd.read_csv(R_REFERENCE_VALUES)

# R and numpy both use IEEE doubles, but the two implementations sum in different orders and
# `energy::dcor` goes through a C routine. Agreement to 1e-10 is far tighter than any
# transcription error could hide in, and looser than bit-equality, which is not achievable
# across languages.
TOL = 1e-10


def r_value(case: str, quantity: str) -> float:
    rows = _VALUES[(_VALUES["case"] == case) & (_VALUES["quantity"] == quantity)]
    if len(rows) != 1:
        raise AssertionError(
            f"expected exactly one fixture row for {case}/{quantity}, got {len(rows)}. "
            "Regenerate with `Rscript rscripts/reference_values.R`."
        )
    return float(rows["value"].iloc[0])


def r_inputs(case: str) -> tuple[np.ndarray, np.ndarray]:
    rows = _INPUTS[_INPUTS["case"] == case].sort_values("i")
    if rows.empty:
        raise AssertionError(f"no fixture inputs for case {case!r}")
    return rows["x"].to_numpy(dtype=float), rows["y"].to_numpy(dtype=float)


def test_the_fixture_records_which_r_produced_it():
    """A reference value with no version behind it cannot be re-derived or challenged."""
    assert r_value("meta", "r_major") >= 4
    for package in ("XICOR", "energy", "randtests", "copula"):
        assert r_value("meta", f"{package}_major") >= 0
        assert r_value("meta", f"{package}_minor") >= 0


# ------------------------------------------------------------------ Chatterjee's xi


def test_xi_matches_XICOR_on_tie_free_data():
    """Where the tie-free reduction and the tie-corrected form coincide, so a regression is unambiguous."""
    x, y = r_inputs("xi_tie_free")
    assert chatterjee_xi(x, y) == pytest.approx(
        r_value("xi_tie_free", "xicor"), abs=TOL
    )


def test_xi_lies_inside_the_XICOR_tie_break_distribution_on_TIED_data():
    """The headline check, and the only external evidence for the tie-corrected form.

    EQUALITY IS NOT THE RIGHT TEST HERE, and finding that out was the point of running it.
    Chatterjee (2021) breaks ties in x uniformly AT RANDOM, so `XICOR::xicor` is a random variable on
    this input: measured, 7 distinct values in 8 calls on identical data, spanning 0.425 to
    0.563. A single R draw is not a reference constant, and asserting our value equals it
    would pin R's RNG state. (The first version of this test did exactly that and failed at
    0.4649 against 0.4976 - a disagreement that looked like a transcription error and was
    not one.)

    What is well defined is the DISTRIBUTION over tie-breaks. We take a deterministic
    stable-sort tie-break, which this repo chose on purpose - "a seeded shuffle inside a
    reproducibility tool is its own problem" - so the correct statement is that our value is
    one admissible member of that distribution, not an outlier from it. That is a weaker
    claim than equality, and it is the strongest claim the estimator actually supports.

    The COST of the deterministic choice is measured, not waved at: `jobs/bench/xi_ties.py`
    quantifies it against a randomised break across tie fractions.
    """
    x, y = r_inputs("xi_tied")
    ours = chatterjee_xi(x, y)
    lo, hi = r_value("xi_tied", "xicor_q001"), r_value("xi_tied", "xicor_q999")
    assert lo <= ours <= hi, (
        f"our tie-break gives {ours:.6f}, outside the 0.1-99.9 percentile band "
        f"[{lo:.6f}, {hi:.6f}] of {int(r_value('xi_tied', 'xicor_draws_n'))} XICOR draws. "
        "That is a transcription difference, not a tie-breaking difference."
    )
    # And it must be a TYPICAL member, not a boundary one: within 3 SD of the mean.
    mean, sd = r_value("xi_tied", "xicor_mean"), r_value("xi_tied", "xicor_sd")
    assert abs(ours - mean) <= 3.0 * sd


def test_the_tied_case_really_does_make_XICOR_random():
    """A positive control on the reasoning above, not on the estimator.

    If a regenerated fixture ever made this case deterministic, the membership test would
    silently become a much weaker statement than its docstring claims.
    """
    assert r_value("xi_tied", "xicor_ties_on_x") > 0
    assert r_value("xi_tied", "xicor_sd") > 0.01


def test_XICOR_is_deterministic_when_only_the_response_is_tied():
    """The other half: with continuous x there is nothing to break, so equality DOES hold.

    This is what localises the disagreement above to X-tie-breaking rather than to the tie-corrected
    form itself - the y-only case ties heavily and still matches to 1e-10.
    """
    assert r_value("xi_tied_y_only", "xicor_ties_on_x") == 0
    assert r_value("xi_tied_y_only", "xicor_distinct_draws") == 1


def test_xi_matches_XICOR_when_only_the_response_is_tied():
    """Ties on y only - the case Chatterjee's Thm 2.1 excludes.

    Our closed-form null is not entitled to this data, and `xi_p_value` routes it to the
    permutation branch. The ESTIMATOR still has to be right, which is what this pins.
    """
    x, y = r_inputs("xi_tied_y_only")
    assert chatterjee_xi(x, y) == pytest.approx(
        r_value("xi_tied_y_only", "xicor"), abs=TOL
    )


def test_the_tie_free_case_is_genuinely_tie_free_and_the_tied_case_genuinely_tied():
    """A positive control on the FIXTURE: the two cases must actually differ in ties.

    Without this, a fixture regenerated with a different seed could quietly make the tied
    case tie-free, and `test_xi_matches_XICOR_on_TIED_data` would keep passing while
    testing nothing - the exact defect CLAUDE.md's positive-control rule names.
    """
    x_free, y_free = r_inputs("xi_tie_free")
    assert len(np.unique(x_free)) == len(x_free)
    assert len(np.unique(y_free)) == len(y_free)

    x_tied, y_tied = r_inputs("xi_tied")
    assert len(np.unique(y_tied)) < len(y_tied), "the tied case has no ties on y"
    # And the two estimator forms must actually disagree here, or the case is not a test.
    n = len(x_tied)
    order = np.argsort(x_tied, kind="stable")
    reduction = 1.0 - 3.0 * np.abs(
        np.diff(rankdata(y_tied, method="max")[order])
    ).sum() / (n**2 - 1)
    assert abs(reduction - chatterjee_xi(x_tied, y_tied)) > 0.01, (
        "the tie-free reduction and the tie-corrected form agree on this case, so it cannot detect the "
        "tie-corrected form"
    )


# ------------------------------------------------------- distance correlation, Spearman


@pytest.mark.parametrize("case", ["xi_tie_free", "xi_tied", "dcor_nonmonotone"])
def test_dcor_matches_the_energy_package(case):
    """`energy` is Szekely and Rizzo's own package for their own statistic."""
    x, y = r_inputs(case)
    assert dcor(x, y) == pytest.approx(r_value(case, "dcor"), abs=1e-9)


@pytest.mark.parametrize(
    "case", ["xi_tie_free", "xi_tied", "xi_tied_y_only", "dcor_nonmonotone"]
)
def test_spearman_matches_r(case):
    """Cheap, and it localises a failure: if xi and Spearman BOTH miss on one case, the
    fixture or the input handling is wrong, not the estimator."""
    x, y = r_inputs(case)
    assert float(spearmanr(x, y).statistic) == pytest.approx(
        r_value(case, "spearman"), abs=TOL
    )


def test_the_three_statistics_disagree_where_they_should():
    """Not a cross-check - a check that the instruments are not redundant.

    On `dcor_nonmonotone` the relation is strong but not monotone. Spearman is near zero,
    while xi and dcor both see it. If all three ever agreed on such a case, one of them
    would not be measuring what its docstring claims and the panel would be carrying three
    columns of the same number.
    """
    rho = abs(r_value("dcor_nonmonotone", "spearman"))
    xi = r_value("dcor_nonmonotone", "xicor")
    dc = r_value("dcor_nonmonotone", "dcor")
    assert xi > rho and dc > rho


# ------------------------------------------------------------- the duration-series checks


@pytest.mark.parametrize("case", ["durations_iid", "durations_trend"])
def test_lag1_rank_autocorrelation_matches_r(case):
    """The quantity C5 is built on, isolated from C5's max-over-lags aggregation.

    Cross-checking C5's own statistic against `randtests::bartels.rank.test` would be an
    apples-to-oranges comparison and is NOT done here: Bartels is the rank von Neumann
    RATIO, a different functional of the same ranks, and our C5 is a studentized maximum
    over lags. What the two genuinely share is the lag-1 rank autocorrelation, so that is
    what is pinned. The Bartels values are in the fixture as context, not as a target.
    """
    durations, _ = r_inputs(case)
    ranks = rankdata(durations)
    ours = float(np.corrcoef(ranks[:-1], ranks[1:])[0, 1])
    assert ours == pytest.approx(r_value(case, "lag1_rank_autocorr"), abs=1e-9)


def test_the_c3_reference_is_pinned_with_its_simulation_seed():
    """C3's p-value is simulated, so a value without its seed is not reproducible.

    This does not re-run R. It pins that the fixture carries the seed, N and lag.max
    alongside the number, so the value can be regenerated and challenged. C3's own
    agreement with this reference is a separate exercise, and it now exists:
    `test_the_bridge_reproduces_the_reference_serial_indep_values` below, behind the `r`
    marker so it skips where there is no interpreter.
    """
    assert r_value("durations_iid", "serial_indep_sim_seed") == 707
    assert r_value("durations_iid", "serial_indep_N") == 1000
    assert r_value("durations_iid", "serial_indep_lag_max") == 5
    p = r_value("durations_iid", "serial_indep_global_p_value")
    assert 0.0 < p <= 1.0
    # iid input: not rejecting is the expected outcome, and a smoke test, not calibration.
    assert p > 0.05


# ------------------------------------- a THIRD implementation, and it needs no fixture


@pytest.mark.parametrize(
    "case", ["xi_tie_free", "xi_tied", "xi_tied_y_only", "dcor_nonmonotone"]
)
def test_xi_matches_scipys_chatterjeexi(case):
    """`scipy.stats.chatterjeexi` (scipy >= 1.15) - independent, Python, no R needed.

    This is STRONGER tier-4 evidence than the R fixture, for a reason that has nothing to do
    with which is more authoritative: it runs unconditionally in the suite. The R comparison
    depends on a committed fixture that only regenerates on a machine with R, so a future
    change to `chatterjee_xi` is checked against a frozen file; this one re-derives the
    reference from a maintained library on every run.

    It also covers the TIED case, which the R comparison cannot pin by equality - XICOR
    breaks x-ties at random. scipy, like us, breaks them deterministically, so here equality
    holds and the tied path gets an exact external check for the first time.
    """
    scipy_xi = pytest.importorskip(
        "scipy.stats", reason="scipy.stats is required to compare against it"
    )
    if not hasattr(scipy_xi, "chatterjeexi"):
        pytest.skip(
            "chatterjeexi needs scipy >= 1.15; the R fixture still covers tier 4"
        )
    x, y = r_inputs(case)
    assert chatterjee_xi(x, y) == pytest.approx(
        float(scipy_xi.chatterjeexi(x, y).statistic), abs=TOL
    )


def test_the_three_implementations_agree_where_all_three_are_defined():
    """Ours, scipy's and R's on one tie-free case - three routes to one number.

    Agreement among three independent implementations is what makes a transcription claim
    hard to argue with. It is stated here as one assertion so a reader does not have to
    assemble it from three parametrised cases.
    """
    scipy_xi = pytest.importorskip("scipy.stats")
    if not hasattr(scipy_xi, "chatterjeexi"):
        pytest.skip("this scipy has no chatterjeexi")
    x, y = r_inputs("xi_tie_free")
    ours = chatterjee_xi(x, y)
    assert ours == pytest.approx(float(scipy_xi.chatterjeexi(x, y).statistic), abs=TOL)
    assert ours == pytest.approx(r_value("xi_tie_free", "xicor"), abs=TOL)


# ------------------------------------------------------- the R bridge, against this fixture
#
# These three need a real interpreter and carry the `r` marker. They skip where `Rscript` is
# absent; they never pass with a mocked one. The tier marker stays `statistical` from the
# module-level `pytestmark`, because `r` is a COST marker on the orthogonal axis.


@pytest.mark.r
def test_the_bridge_reproduces_the_reference_serial_indep_values(requires_rscript):
    """Oracle: `copula::serialIndepTest` itself, via jobs/reference/r_reference_values.csv.

    BOTH quantities are asserted, and that is the point of the test. The statistic is
    invariant to the simulation seed and to N, so a statistic-only assertion cannot detect a
    swapped `seed`/`n_sim` argument order in the bridge. The p-value is a function of the
    simulated null and therefore of both, so it is what pins the call.

    `analyzers/instrument_validation.py` reported this comparison as tier-4 ABSENT and named
    it exactly, on the grounds that it "needs R at test time, which the suite refuses to
    require". The `r` tier is what makes it requirable.
    """
    from quebra.analyzers.checks import c3_serial_copula as c3
    from quebra.analyzers.checks.result import CLOCK_IN_SPEC, Segment

    x, _unused = r_inputs("durations_iid")
    seed = int(r_value("durations_iid", "serial_indep_sim_seed"))
    n_sim = int(r_value("durations_iid", "serial_indep_N"))
    lag_max = int(r_value("durations_iid", "serial_indep_lag_max"))

    # `run`, not `_invoke_rscript`. The ledger calls `run`, and the argument plumbing
    # between the two is what a helper-level test leaves uncovered: swapping run's
    # positional arguments to the helper is a mutation the whole suite survives, and the
    # shipped row would then carry a seed that does not reproduce its own p-value.
    result = c3.run(
        [Segment(x=x, tau=float(x.sum()), n_censored_dropped=0)],
        clock=CLOCK_IN_SPEC,
        max_lag=lag_max,
        seed=seed,
        n_null_sim=n_sim,
    )
    statistic, p_value = result.statistic, result.p_value

    assert f"seed={seed}" in result.notes and f"N={n_sim}" in result.notes, (
        f"the row must record the settings that produced it: {result.notes}"
    )
    assert statistic == pytest.approx(
        r_value("durations_iid", "serial_indep_global_statistic"), abs=TOL
    )
    assert p_value == pytest.approx(
        r_value("durations_iid", "serial_indep_global_p_value"), abs=TOL
    )


@pytest.mark.r
def test_the_fixture_meta_rows_match_the_local_r(requires_rscript):
    """Oracle: the installed R, asked directly.

    `test_the_fixture_records_which_r_produced_it` asserts these rows are PRESENT, which is
    a real property and runs everywhere. This asserts they are TRUE of the machine, which
    only an interpreter can answer. A fixture whose recorded versions have drifted from the
    R that can be run is a fixture whose values nobody can re-derive.
    """
    import subprocess

    probe = (
        "cat(R.version$major, R.version$minor, "
        + ", ".join(
            f'as.character(packageVersion("{p}"))'
            for p in ("XICOR", "energy", "randtests", "copula")
        )
        + ")"
    )
    out = subprocess.run(
        [requires_rscript, "-e", probe],
        capture_output=True,
        text=True,
        errors="replace",
        timeout=120,
        check=True,
    ).stdout.split()

    r_major, r_minor_full = out[0], out[1]
    assert int(r_major) == int(r_value("meta", "r_major"))
    assert int(r_minor_full.split(".")[0]) == int(r_value("meta", "r_minor"))

    packages = ("XICOR", "energy", "randtests", "copula")
    for package, reported in zip(packages, out[2:], strict=True):
        parts = reported.replace("-", ".").split(".")
        assert int(parts[0]) == int(r_value("meta", f"{package}_major")), package
        assert int(parts[1]) == int(r_value("meta", f"{package}_minor")), package


@pytest.mark.r
def test_the_bridge_times_out_rather_than_hanging(requires_rscript):
    """Oracle: `subprocess.run`'s own timeout contract.

    The bridge's caller in `check_ledger` catches `TimeoutExpired` and writes a `declined:`
    row, so a bridge that swallowed the timeout instead of raising would turn a hung
    interpreter into a silent `not computed`. This pins the raise, not the ledger's handling
    of it. The ledger does not pass `timeout_s`, so a ledger row runs under the 900 s
    default: that branch is reachable in production, and simply not cheap to reach here.
    """
    import subprocess

    from quebra.analyzers.checks import c3_serial_copula as c3
    from quebra.analyzers.checks.result import CLOCK_IN_SPEC, Segment

    x, _unused = r_inputs("durations_iid")
    segment = Segment(x=x, tau=float(x.sum()), n_censored_dropped=0)

    with pytest.raises(subprocess.TimeoutExpired) as excinfo:
        c3.run([segment], clock=CLOCK_IN_SPEC, seed=707, timeout_s=0.5, n_null_sim=1000)
    # `_invoke_rscript` also calls `r_library_paths`, whose own 60 s ceiling raises the same
    # type, so the timeout value is what says which subprocess this was.
    assert excinfo.value.timeout == 0.5
