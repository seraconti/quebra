"""The carve recovers a known duration law from a metric record: record -> carve -> durations.

Oracle: analytic, from the generating process. iid reads thresholded at the occupancy-matched
quantile are in spec with probability `q` independently at each read, so an in-spec run of
length `k` has probability `(1 - q) * q**(k-1)` and a complete window's duration is its run
length times the read interval. No small-`dt` limit is involved and no Monte Carlo truth is
needed: the law is exact.

**`q` is deliberately NOT 1/2, and an earlier version of this file used 1/2.** At `q = 1/2`
the law is symmetric under an inversion of the in-spec convention: swapping
`windows.in_spec_mask`'s two branches leaves the run-length distribution Geometric(1/2) and
this file passed, measured (chi2 7.30, p = 0.199). `AGENTS.md` section 5 calls inverting the
in-spec test a scientific error rather than a style problem, so the one end-to-end
record-to-carve-to-duration-law test must not be blind to it. At `q = 0.3` an inversion makes
the reads in spec with probability 0.7, the runs three times longer, and the chi-square
rejects outright.

This construction is Arm B's in `jobs/bench/arms.py`, reproduced in four lines rather than
imported: `arm_b` routes through `jobs/bench/carve.py` and returns `Segment` objects already
filtered by `min_events` and stripped of their trailing censored window, so it would test the
bench's carve - already pinned column-for-column by `tests/test_bench_uses_real_carve.py` -
and would lose the birth and death taxonomy this test selects on. `windows.run` is the
subject; the bench is where the same construction is scored for size and power.

This is the only test in the suite that closes metric record to carve to duration law end to
end. Arms A, D and E emit `Segment` objects directly and never touch a threshold, so they
exercise the estimators but not the carve that feeds them.

**No committed fixture.** `spec/spectests06.md` R6.4 records why an earlier draft's plan to
commit a CSV was wrong: the geometric law is a statement about a DISTRIBUTION and a file is
one realisation of it, whose byte-identity is exactly what stops the seed being varied. A
single draw cannot separate a correct carve from a subtly wrong one.

**Seeds are written out, and the draws are POOLED into one test.** Asserting per seed would
make one weak test per seed and a multiplicity problem - at alpha = 0.05 about one in twenty
fails by chance - which is the flakiness `AGENTS.md` section 4 exists to prevent. Pooling
gives one verdict with the whole sample behind it and nothing to correct for. The seed tuple
follows the idiom at `tests/test_tier3_calibration.py`: written out rather than derived, so
a rerun reproduces exactly.
"""

from __future__ import annotations

import numpy as np
import pytest

from quebra.analyzers import windows

pytestmark = pytest.mark.statistical

# One fixed seed per draw, written out rather than derived, so a rerun reproduces exactly.
SEEDS = (11, 23, 37, 41)
READS_PER_DRAW = 2000
READ_DT_S = 1.0
# NOT 1/2: at 1/2 the run-length law is invariant under an inversion of the in-spec
# convention, and this test would not notice one. See the module docstring.
IN_SPEC_P = 0.3
# Bins for the goodness-of-fit: run lengths 1..4 and a 5-or-more tail. Chosen so every
# expected count clears 5 at the sample size these seeds produce.
BIN_EDGES = (1, 2, 3, 4)


def _complete_run_lengths(seed: int) -> np.ndarray:
    """Carve one iid draw with the real analyzer and return complete windows' run lengths.

    Only `up_crossing` births with `down_crossing` deaths are complete lifetimes. The first
    and last windows of a record are an endurance bag and a censored window respectively,
    and including either would bias the law being tested.
    """
    rng = np.random.default_rng(seed)
    values = rng.standard_normal(READS_PER_DRAW)
    threshold = float(np.quantile(values, 1.0 - IN_SPEC_P))
    t_rel_s = np.arange(READS_PER_DRAW, dtype=float) * READ_DT_S

    carved = windows.run(
        windows.WindowsInputs(
            t_rel_s=t_rel_s,
            values=values,
            thresholds=[("q70", threshold, True)],
            dataset_id="carve-law",
        )
    )
    w = carved.windows_in_spec
    complete = w[
        (w["birth_type"] == "up_crossing") & (w["death_type"] == "down_crossing")
    ]
    return (complete["duration_s"].to_numpy() / READ_DT_S).round().astype(int)


def _binned(run_lengths: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Observed and expected bin counts under Geometric(1 - q) on {1, 2, ...}.

    The pmf is written as a function of `IN_SPEC_P` rather than as a literal, because a bare
    constant swap would leave the expected frequencies summing to something other than the
    observed total and `scipy.stats.chisquare` would RAISE rather than assert - a change in
    the constant would then break the test for the wrong reason.
    """
    q = IN_SPEC_P
    n = len(run_lengths)
    observed = np.array(
        [np.count_nonzero(run_lengths == k) for k in BIN_EDGES]
        + [np.count_nonzero(run_lengths > BIN_EDGES[-1])],
        dtype=float,
    )
    probs = np.array(
        [(1.0 - q) * q ** (k - 1) for k in BIN_EDGES] + [q ** BIN_EDGES[-1]]
    )
    assert probs.sum() == pytest.approx(1.0, abs=1e-12), (
        "the bin scheme must be complete"
    )
    return observed, probs * n


def test_carved_durations_follow_the_geometric_law_of_the_generating_process() -> None:
    """Oracle: Geometric(1 - q) on {1, 2, ...}, exact for this construction.

    One chi-square on the POOLED sample, with the aggregation and the bin scheme stated
    before the verdict is read off it. Bins are run lengths 1..4 and a 5-or-more tail;
    `q` is known rather than estimated, so the null has 4 degrees of freedom.

    A duration is `(index of the first out-of-spec read) - (index of the birth read)`, so
    an off-by-one in the carve's death index shifts the whole distribution by one and this
    fails on the very first bin.
    """
    from scipy.stats import chisquare

    per_seed = {seed: _complete_run_lengths(seed) for seed in SEEDS}
    pooled = np.concatenate(list(per_seed.values()))
    assert len(pooled) > 1000, f"only {len(pooled)} complete windows; the test is thin"

    observed, expected = _binned(pooled)
    result = chisquare(observed, expected)
    means = ", ".join(f"seed {s}: {v.mean():.3f}" for s, v in per_seed.items())
    assert result.pvalue > 0.001, (
        f"pooled run lengths depart from Geometric(1 - q) at q={IN_SPEC_P}: "
        f"chi2={result.statistic:.2f}, p={result.pvalue:.4g} over {len(pooled)} complete "
        f"windows from {len(SEEDS)} seeds. Per-seed mean run length ({means}); the law's "
        f"mean is {1 / (1 - IN_SPEC_P):.2f}"
    )


def test_the_law_is_a_property_of_the_carve_and_not_of_the_draw() -> None:
    """Negative control: a shifted duration must fail the same chi-square.

    Without this the test above could pass on any distribution whose first six bins happen
    to look geometric, and there would be no evidence the comparison discriminates. Adding
    one read interval to every duration - exactly what an off-by-one death index would do -
    must reject.
    """
    from scipy.stats import chisquare

    pooled = np.concatenate([_complete_run_lengths(seed) for seed in SEEDS]) + 1

    observed, expected = _binned(pooled)
    assert chisquare(observed, expected).pvalue < 1e-6, (
        "a uniformly shifted duration must reject; if it does not, the chi-square in the "
        "test above is not discriminating and its pass means nothing"
    )
