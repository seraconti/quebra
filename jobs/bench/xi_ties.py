"""How far can Chatterjee's xi be trusted as the response ties?

Three questions, one experiment, because they share a generator and separating them would
mean three chances for the data to differ.

**Q1 - where does the closed-form p-value stop agreeing with the permutation p-value?**
`shape_stats.xi_p_value` branches on `XI_TIE_CUTOFF = 0.0`: any tie at all routes to the
permutation. That cutoff was set to zero because zero is where the two forms provably
coincide, not because anything was measured. This measures the cost of being wrong about it
- the DELIVERABLE is the tie fraction at which the two p-values first differ by more than
0.02, as a number.

**Q2 - what does our deterministic X-tie break cost against a randomised one?**
Chatterjee (2021) breaks ties in x uniformly at random. We use a deterministic stable sort,
on purpose: "a seeded shuffle inside a reproducibility tool is its own problem". Tier 4
established that this makes `XICOR::xicor` a random variable while ours is a constant -
7 distinct values in 8 calls on the same tied input. So the honest question is not whether
we match a draw, but where our fixed choice sits in the distribution of admissible ones,
and how wide that distribution gets.

**Q3 - does the asymptotic null survive ties at all?**
Chatterjee's Thm 2.1 gives `sqrt(n) xi -> N(0, 2/5)` for CONTINUOUS y. Under ties the
variance is not 2/5 and the closed form is not entitled to the data. This measures the
resulting type-I error directly, which is what "not entitled" costs in practice.

**This is a study, not a pipeline layer.** Nothing outside `jobs/bench/` may import it
(`tests/test_bench_isolation.py` enforces that). It writes
`jobs/bench/results/xi_tie_experiment.csv`; a figure that wants these numbers declares that CSV
as a Dataset and loads it through the pipeline, so the dependency runs through provenance
rather than around it.

Run:  python jobs/bench/xi_ties.py
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata


from quebra.analyzers.permutation import paired_permutation_test
from quebra.analyzers.shape_stats import (
    chatterjee_xi,
    tie_fraction,
    xi_p_value_asymptotic,
)

RESULTS = Path(__file__).resolve().parent / "results" / "xi_tie_experiment.csv"

# n values chosen to bracket what this project actually runs on. Reference 10.3 caveat 2
# flags 35-355 as the disputed range; a shipped window is 5 to a few dozen reads, and the
# pooled statistic sees a few hundred.
N_VALUES = (35, 75, 150, 355)

# Ties are induced by ROUNDING y to a grid, which is how ties arise here for real: the
# instrument reads a quantised metric. A coarser grid means a higher tie fraction. The grid
# is expressed as the number of distinct levels y is collapsed onto.
Y_LEVELS = (0, 20, 10, 5, 3, 2)  # 0 means "no rounding", the tie-free control

N_REPLICATES = 300
N_PERM = 999
N_TIE_BREAKS = 200  # draws of a RANDOM x-tie break, for Q2
ALPHA = 0.05

# The deliverable's threshold. 0.02 on a p-value is the point at which a reader comparing
# against 0.05 or 0.01 could be led to a different verdict.
DIVERGENCE_THRESHOLD = 0.02

# THE THRESHOLD MUST BE READ AGAINST THIS. The permutation p-value is itself estimated from
# B resamples, so it carries Monte Carlo error of about `sqrt(p(1-p)/B)` - at B = 999 that is
# 0.0158 near p = 0.5, which is most of the 0.02 threshold. So `mean |p_closed - p_perm|`
# conflates real disagreement with the permutation's own noise and would report a divergence
# that is not there.
# The experiment therefore records the SIGNED mean difference as well. Noise averages out of
# the signed mean and systematic disagreement does not, which makes the signed quantity the
# one the deliverable is read off. `p_abs_diff_mean` is kept because it bounds the per-window
# error a reader would actually see, but it is labelled, not used as the answer.
MC_FLOOR = 0.5 / np.sqrt(N_PERM)


@dataclass
class Cell:
    n: int
    y_levels: int
    tie_fraction_y: float = 0.0
    tie_fraction_x: float = 0.0
    xi_mean: float = 0.0
    p_closed_mean: float = 0.0
    p_perm_mean: float = 0.0
    p_signed_diff_mean: float = 0.0
    p_signed_diff_se: float = 0.0
    p_abs_diff_mean: float = 0.0
    p_abs_diff_p90: float = 0.0
    type_i_closed: float = 0.0
    type_i_perm: float = 0.0
    tie_fraction_x_q2: float = 0.0
    xi_tie_break_sd: float = 0.0
    xi_tie_break_range: float = 0.0
    xi_deterministic_percentile: float = float("nan")
    details: dict = field(default_factory=dict)


def _quantise(v: np.ndarray, levels: int) -> np.ndarray:
    """Collapse `v` onto `levels` equal-width bins spanning its observed range."""
    lo, hi = float(v.min()), float(v.max())
    if levels <= 0 or hi <= lo:
        return v
    edges = np.linspace(lo, hi, levels + 1)
    return np.clip(np.digitize(v, edges[1:-1]), 0, levels - 1).astype(float)


def _make_pair(
    rng: np.random.Generator,
    n: int,
    y_levels: int,
    dependent: bool,
    x_levels: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """An (x, y) pair, with y and optionally x quantised onto a fixed number of levels.

    `dependent=False` is the NULL used for type-I error; `dependent=True` carries a real
    non-monotone relation, which is the case xi exists to detect.

    `x_levels` is for Q2 ONLY, and defaults to off. Q1 and Q3 leave x continuous because
    that is the shipped situation - x is the read time and y is the quantised metric, so
    ties arise on y. Q2 asks what the deterministic X-tie break costs, and a first version
    of this experiment reported a tie-break spread of exactly 0.0000 in every cell because
    it never generated an x-tie for the tie break to act on. A measurement of something
    that cannot happen is not a measurement, and a column of zeros would have read as
    "the deterministic choice costs nothing".
    """
    x = rng.normal(size=n)
    y = (x**2 if dependent else rng.normal(size=n)) + rng.normal(size=n, scale=0.3)
    if x_levels > 0:
        x = _quantise(x, x_levels)
    if y_levels > 0:
        y = _quantise(y, y_levels)
    return x, y


def _xi_random_tie_break(
    x: np.ndarray, y: np.ndarray, rng: np.random.Generator
) -> float:
    """The tie-corrected xi with x-ties broken uniformly at random, as Chatterjee (2021) says.

    Implemented here rather than in `shape_stats` deliberately: the shipped estimator is
    deterministic by decision, and this exists only to measure what that decision costs.
    """
    n = len(x)
    jitter_order = np.lexsort((rng.random(n), x))
    r = rankdata(y, method="max")[jitter_order]
    ell = rankdata(-y, method="max")
    denominator = 2.0 * float(np.sum(ell * (n - ell)))
    if denominator == 0.0:
        return float("nan")
    return float(1.0 - n * float(np.abs(np.diff(r)).sum()) / denominator)


def run_cell(n: int, y_levels: int, seed: int) -> Cell:
    rng = np.random.default_rng(seed)
    cell = Cell(n=n, y_levels=y_levels)

    xis, p_closed, p_perm, tf_y, tf_x = [], [], [], [], []
    reject_closed = reject_perm = 0

    for _ in range(N_REPLICATES):
        # Q1 and Q3 use the INDEPENDENT pair: that is the null both p-values claim to
        # calibrate against, so a disagreement there is a calibration disagreement.
        x, y = _make_pair(rng, n, y_levels, dependent=False)
        if len(np.unique(y)) < 2:
            continue
        xi = chatterjee_xi(x, y)
        pc = xi_p_value_asymptotic(xi, n)
        pp = paired_permutation_test(
            x, y, chatterjee_xi, seed=int(rng.integers(1 << 31)), n_resamples=N_PERM
        ).p_value
        xis.append(xi)
        p_closed.append(pc)
        p_perm.append(pp)
        tf_y.append(tie_fraction(y))
        tf_x.append(tie_fraction(x))
        reject_closed += int(pc <= ALPHA)
        reject_perm += int(pp <= ALPHA)

    signed = np.asarray(p_closed) - np.asarray(p_perm)
    diffs = np.abs(signed)
    cell.p_signed_diff_mean = float(np.mean(signed))
    cell.p_signed_diff_se = float(np.std(signed) / np.sqrt(max(len(signed), 1)))
    cell.tie_fraction_y = float(np.mean(tf_y))
    cell.tie_fraction_x = float(np.mean(tf_x))
    cell.xi_mean = float(np.mean(xis))
    cell.p_closed_mean = float(np.mean(p_closed))
    cell.p_perm_mean = float(np.mean(p_perm))
    cell.p_abs_diff_mean = float(np.mean(diffs))
    cell.p_abs_diff_p90 = float(np.percentile(diffs, 90))
    cell.type_i_closed = reject_closed / max(len(p_closed), 1)
    cell.type_i_perm = reject_perm / max(len(p_perm), 1)

    # Q2: the spread over random x-tie breaks, on ONE representative dependent pair, with
    # our deterministic value located inside it as a percentile.
    # X IS QUANTISED HERE AND NOWHERE ELSE - see `_make_pair`. The x grid is tied to the y
    # grid so the table has one "how coarse is this cell" axis rather than a second sweep.
    x, y = _make_pair(rng, n, y_levels, dependent=True, x_levels=y_levels)
    if len(np.unique(y)) >= 2:
        draws = np.array([_xi_random_tie_break(x, y, rng) for _ in range(N_TIE_BREAKS)])
        draws = draws[np.isfinite(draws)]
        ours = chatterjee_xi(x, y)
        if draws.size:
            cell.tie_fraction_x_q2 = float(tie_fraction(x))
            cell.xi_tie_break_sd = float(np.std(draws))
            cell.xi_tie_break_range = float(draws.max() - draws.min())
            cell.xi_deterministic_percentile = float(100.0 * np.mean(draws <= ours))
    return cell


def main() -> None:
    rows = []
    for n in N_VALUES:
        for y_levels in Y_LEVELS:
            cell = run_cell(n, y_levels, seed=hash((n, y_levels)) % (1 << 31))
            rows.append({k: v for k, v in cell.__dict__.items() if k != "details"})
            print(
                f"n={n:<4} levels={y_levels:<3} tie_y={cell.tie_fraction_y:.3f}  "
                f"signed={cell.p_signed_diff_mean:+.4f}+/-{cell.p_signed_diff_se:.4f} "
                f"|dp|={cell.p_abs_diff_mean:.4f}  "
                f"typeI closed={cell.type_i_closed:.3f} perm={cell.type_i_perm:.3f}  "
                f"tiebreak(x_tie={cell.tie_fraction_x_q2:.2f}) "
                f"sd={cell.xi_tie_break_sd:.4f} "
                f"pct={cell.xi_deterministic_percentile:.0f}",
                flush=True,
            )
    frame = pd.DataFrame(rows)
    RESULTS.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(RESULTS, index=False)

    print(f"\nwrote {RESULTS}")
    print("\n=== DELIVERABLE ===")
    print(
        f"Tie fraction at which the closed form and the permutation p-value first diverge "
        f"SYSTEMATICALLY by more than {DIVERGENCE_THRESHOLD}.\n"
        f"Read off the SIGNED mean over {N_REPLICATES} replicates, NOT |dp|: at B="
        f"{N_PERM} the permutation p carries Monte Carlo error near {MC_FLOOR:.4f}, which "
        f"is most of the threshold, so |dp| reports its own noise as a divergence. "
        f"(Measured: the tie-FREE control cell, where the true difference is ~0, still "
        f"gives |dp| ~ 0.010.)"
    )
    for n in N_VALUES:
        sub = frame[frame["n"] == n].sort_values("tie_fraction_y")
        crossed = sub[sub["p_signed_diff_mean"].abs() > DIVERGENCE_THRESHOLD]
        if crossed.empty:
            worst = sub.loc[sub["p_signed_diff_mean"].abs().idxmax()]
            print(
                f"  n={n:<4} NEVER diverges systematically; worst is "
                f"{worst['p_signed_diff_mean']:+.4f} +/-{worst['p_signed_diff_se']:.4f} "
                f"at tie fraction {worst['tie_fraction_y']:.3f} "
                f"({int(worst['y_levels'])} levels)"
            )
        else:
            first = crossed.iloc[0]
            print(
                f"  n={n:<4} tie fraction {first['tie_fraction_y']:.3f} "
                f"({int(first['y_levels'])} levels), signed mean = "
                f"{first['p_signed_diff_mean']:+.4f} "
                f"+/-{first['p_signed_diff_se']:.4f}"
            )

    print(
        "\nType-I error of the CLOSED FORM at alpha = 0.05 - the operational consequence, "
        "and the one that actually moves. The permutation column is the reference."
    )
    for _, r in frame.sort_values(["n", "tie_fraction_y"]).iterrows():
        flag = " <-- inflated" if r["type_i_closed"] > 0.10 else ""
        print(
            f"  n={int(r['n']):<4} levels={int(r['y_levels']):<3} "
            f"tie_y={r['tie_fraction_y']:.3f}  closed={r['type_i_closed']:.3f}  "
            f"perm={r['type_i_perm']:.3f}{flag}"
        )


if __name__ == "__main__":
    main()
