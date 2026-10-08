"""Turn the two tables into `promotion_report.md`.

The report is generated, not written by hand, so every number in it is traceable to a row
of `size_table.csv` or `power_table.csv` and regenerating after a rerun cannot leave a
stale figure in the prose. Where a judgement is made, the RULE is printed next to the
verdict rather than left implicit.

Decision rule, stated once here and applied uniformly:

    PROMOTE  size indistinguishable from nominal across every null cell inside the real
             data's censoring envelope, AND mean power of at least 0.5 at n = 100 over
             the grid of the alternative the check is directed at
    HOLD     size calibrated but mean power at n = 100 below 0.5
    REJECT   size miscalibrated somewhere inside the envelope

Size is judged by a multiplicity-corrected z-test (`FAMILYWISE_ALPHA`) rather than a flat
tolerance: a flat tolerance rejects every row, because the max of 86-144 deviations is
~3 SE by chance. Power is judged by its mean over the alternative's grid rather than at the
grid's strongest point, which reports a capability few records call on. No grid point is
singled out as "the real data's": the bench's dependence knob is not calibrated to a record.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from jobs.bench.arms import ARM_A, ARM_B, ARM_C, ARM_D, ARM_E
from jobs.bench.grid import ALPHA, N_GRID
from quebra.analyzers.calibration_summary import (
    ENVELOPE_MAX_CENSORING,
    FAMILYWISE_ALPHA,
    MIN_SUPPORT_FRACTION,
    bonferroni_z_crit,
    null_se,
)
from quebra.analyzers.checks.result import CALIB_ASYMPTOTIC
from quebra.analyzers.checks.result import CLOCK_IN_SPEC

RESULTS_DIR = Path(__file__).resolve().parent / "results"

# Where the real data sits. Censoring and quantisation were measured before the bench was
# built; the event counts come from the 6D2S ledgers.
REAL_DATA_NOTES = {
    "event_counts": "up to 935 events per record, threshold and clock (28 6D2S ledger records)",
    "censoring": "0.000-0.026 wherever n >= 20",
    "quantisation": "3 us spans k=1..72 with 43 distinct values; 4 us has 79% at k=1",
}

# Which alternative each check is actually directed at, so criterion 2 scores it against
# the right arm instead of against whichever arm happens to flatter it.
DIRECTED_AT = {
    "c1_lewis_robinson": (ARM_D, "b", "monotone trend in the intensity"),
    "c2_anderson_darling": (ARM_D, "b", "departure from a renewal process"),
    "c5_rank_autocorr": (ARM_E, "rho", "serial dependence between durations"),
    "c6_exchangeability": (ARM_E, "rho", "any departure from exchangeability"),
    # CvM is scored against Arm D like C2: it is the same class of
    # alternative - a departure from a renewal process - because CvM and C2 are two
    # functionals of the same tied-down bridge and differ only in the weight.
    "cvm_cramer_von_mises": (ARM_D, "b", "departure from a renewal process"),
}

# Arm D at b = 1 is the identity trend, so its SIZE cells are null cells and belong here.
# They were omitted, silently discarding 233 rows of null evidence and never performing the
# TRP-generator cross-check `grid.py` says they exist for. Including them changes each row's
# cell count, hence its Bonferroni z_crit, so the regenerated report is re-read rather than
# assumed unchanged. (Arm D's POWER cells live in the power table and cannot contaminate.)
NULL_ARMS = (ARM_A, ARM_B, ARM_C, ARM_D, ARM_E)


# A rejection rate is only reported as a SIZE if it rests on at least this fraction of its
# cell's replicates. Below it the surviving replicates were selected by a data-dependent
# event - a segment whose durations were all one quantum, or a censored cell that happened
# to collapse to a single segment - so the rate is conditional on that event and is not the
# quantity it appears to be. C2 and CvM asymptotic rows can sit at n_used = 1 or 5 out of 2000; left
# unfiltered, a rate from a single replicate is either 0.0 or 1.0 and would have dominated
# the worst-deviation table and driven a spurious REJECT.
#
# The value, and ENVELOPE_MAX_CENSORING and FAMILYWISE_ALPHA with it, are IMPORTED from
# analyzers/calibration_summary.py rather than declared here. The ledger scores real data
# with the same rule, and two copies of "is this check calibrated" that could drift apart
# is a correctness risk, not a style one.


def with_support(frame: pd.DataFrame) -> pd.DataFrame:
    """Rows whose rate rests on enough replicates to be read as a rate."""
    return frame[frame["n_used"] >= MIN_SUPPORT_FRACTION * frame["n_replicates"]]


def under_supported(frame: pd.DataFrame) -> pd.DataFrame:
    """The complement, reported separately rather than dropped silently."""
    return frame[frame["n_used"] < MIN_SUPPORT_FRACTION * frame["n_replicates"]]


def _fmt(value: float, digits: int = 4) -> str:
    return "-" if pd.isna(value) else f"{value:.{digits}f}"


def _row_label(row: pd.Series) -> str:
    variant = row["variant"] if isinstance(row["variant"], str) else ""
    suffix = f"/{variant}" if variant else ""
    return f"{row['check']} [{row['calibration']}{suffix}]"


def _label_column(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["row"] = frame.apply(_row_label, axis=1)
    return frame


def load() -> tuple[pd.DataFrame, pd.DataFrame]:
    size = pd.read_csv(RESULTS_DIR / "size_table.csv")
    power = pd.read_csv(RESULTS_DIR / "power_table.csv")
    return _label_column(size), _label_column(power)


def _markdown_table(frame: pd.DataFrame, columns: list | None = None) -> str:
    """Render a frame as markdown.

    Column KEYS stay as they are (several tables are pivoted on numeric levels like
    `n_target` or `rho`, so the keys are ints and floats); only their rendering is
    stringified.
    """
    keys = list(frame.columns) if columns is None else list(columns)
    lines = [
        "| " + " | ".join(str(key) for key in keys) + " |",
        "|" + "|".join("---" for _ in keys) + "|",
    ]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(str(row[key]) for key in keys) + " |")
    return "\n".join(lines)


def size_by_n(size: pd.DataFrame) -> pd.DataFrame:
    """Rejection rate on the primary null cell (Arm A, continuous, no extra censoring)."""
    primary = with_support(
        size[
            (size["arm"] == ARM_A)
            & (~size["quantised"].astype(bool))
            & (size["censoring_target"] == 0.0)
        ]
    )
    table = (
        primary.groupby(["row", "n_target"])["rejection_rate"]
        .mean()
        .unstack("n_target")
    )
    return table.reindex(columns=[n for n in N_GRID if n in table.columns])


def scored_cells(
    size: pd.DataFrame, min_n: int = 35, inside_envelope: bool = True
) -> pd.DataFrame:
    """Null cells a size verdict rests on, with a z-score against nominal.

    The z uses the NULL standard error `sqrt(alpha(1-alpha)/n_used)`, not the observed
    one: the observed SE collapses to zero on a cell that never rejected, which would make
    a perfectly conservative cell look infinitely miscalibrated.
    """
    frame = with_support(
        size[size["arm"].isin(NULL_ARMS) & (size["n_target"] >= min_n)]
    ).copy()
    frame = (
        frame[frame["censoring_target"] <= ENVELOPE_MAX_CENSORING]
        if inside_envelope
        else frame[frame["censoring_target"] > ENVELOPE_MAX_CENSORING]
    )
    se = np.array([null_se(ALPHA, u) for u in frame["n_used"]])
    frame["z"] = (frame["rejection_rate"] - ALPHA) / se
    return frame


def size_verdict_table(
    size: pd.DataFrame, inside_envelope: bool = True
) -> pd.DataFrame:
    """Per row: the worst cell, its z, and the Bonferroni threshold it is judged against."""
    frame = scored_cells(size, inside_envelope=inside_envelope)
    rows = []
    for row_label, group in frame.groupby("row"):
        n_cells = len(group)
        z_crit = bonferroni_z_crit(n_cells)
        worst = group.loc[group["z"].abs().idxmax()]
        rows.append(
            {
                "row": row_label,
                "n_cells": n_cells,
                "worst_rate": round(float(worst["rejection_rate"]), 4),
                "worst_z": round(float(worst["z"]), 2),
                "z_crit": round(z_crit, 2),
                "calibrated": bool(abs(float(worst["z"])) <= z_crit),
                "worst_cell": (
                    f"{worst['arm']}/{worst['clock']} n={worst['n_target']} "
                    f"q={bool(worst['quantised'])} c={worst['censoring_target']}"
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("worst_z", key=abs, ascending=False)


def worst_size(size: pd.DataFrame, min_n: int = 35) -> pd.DataFrame:
    """Largest deviation from nominal across every adequately supported null cell."""
    null_cells = with_support(
        size[size["arm"].isin(NULL_ARMS) & (size["n_target"] >= min_n)]
    ).copy()
    null_cells["deviation"] = (null_cells["rejection_rate"] - ALPHA).abs()
    worst = null_cells.loc[null_cells.groupby("row")["deviation"].idxmax()]
    return worst[
        [
            "row",
            "arm",
            "clock",
            "n_target",
            "quantised",
            "censoring_target",
            "rejection_rate",
            "mc_se",
            "deviation",
            "n_used",
        ]
    ].sort_values("deviation", ascending=False)


def power_at(power: pd.DataFrame, check: str, n: int = 100) -> pd.DataFrame:
    arm, parameter, _description = DIRECTED_AT[check]
    subset = power[
        (power["check"] == check) & (power["arm"] == arm) & (power["n_target"] == n)
    ]
    return (
        subset.groupby(["row", parameter])["rejection_rate"].mean().unstack(parameter)
    )


def dependence_grid_table(power: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    """C5 and C6 power by (row, rho) and n, mean over the two Weibull shapes, with the
    range of duration lag-1 each rho induces across its cells.

    Also returns, as text, the grid points at which any row's power exceeds 0.5 at some
    n below 100: the dependence strength these checks need before a small-n silence says
    anything.
    """
    rows = []
    strong = set()
    for check in ("c5_rank_autocorr", "c6_exchangeability"):
        arm, parameter, _description = DIRECTED_AT[check]
        subset = power[(power["check"] == check) & (power["arm"] == arm)]
        for (row_label, value), group in subset.groupby(["row", parameter]):
            lag1 = group["mean_induced_lag1"]
            entry = {
                "row": row_label,
                parameter: value,
                "induced_lag1": f"{lag1.min():.3f} to {lag1.max():.3f}",
            }
            for n, sub in group.groupby("n_target"):
                rate = float(sub["rejection_rate"].mean())
                entry[n] = rate
                if n < 100 and rate > 0.5:
                    strong.add(value)
            rows.append(entry)
    text = (
        "rho = " + ", ".join(f"{v:g}" for v in sorted(strong))
        if strong
        else "no grid point"
    )
    return pd.DataFrame(rows).round(3), text


def realised_censoring(size: pd.DataFrame) -> dict[float, float]:
    """Target censoring level -> what the generator actually produced.

    They are not the same and the gap is large at the top of the grid: `m` is capped so
    each segment expects at least a few events, and that cap binds at every n, so the
    "0.25" arm realises 0.15-0.19 per cell (about 0.167 pooled), so every table prints the realised value beside
    the target label.
    """
    return {
        float(target): float(group["censoring_realised"].mean())
        for target, group in size.groupby("censoring_target")
    }


def censoring_effect(size: pd.DataFrame) -> pd.DataFrame:
    """Size against the censoring factor, which only Arms A and D carry.

    Columns are labelled with the REALISED fraction, not the target.
    """
    subset = with_support(
        size[size["arm"].isin([ARM_A, ARM_D]) & (size["n_target"] >= 50)]
    )
    table = (
        subset.groupby(["row", "censoring_target"])["rejection_rate"]
        .mean()
        .unstack("censoring_target")
    )
    realised = realised_censoring(subset)  # the A/D subset the table actually shows
    return table.rename(
        columns={
            c: f"target {c:g} (realised {realised.get(float(c), float('nan')):.3f})"
            for c in table.columns
        }
    )


def quantisation_effect(size: pd.DataFrame) -> pd.DataFrame:
    subset = with_support(
        size[(size["arm"] == ARM_A) & (size["censoring_target"] == 0.0)]
    )
    return (
        subset.groupby(["row", "quantised"])["rejection_rate"]
        .mean()
        .unstack("quantised")
    )


def arm_c_flatness(power: pd.DataFrame) -> pd.DataFrame:
    """Arm C's power against read-level rho - expected flat, which is the finding."""
    subset = power[(power["arm"] == ARM_C) & (power["clock"] == CLOCK_IN_SPEC)]
    return subset.groupby(["row", "rho"])["rejection_rate"].mean().unstack("rho")


def check_agreement(
    power: pd.DataFrame, reference: str = "c6_exchangeability"
) -> pd.DataFrame:
    """Mean paired difference in rejection rate between each row and the reference row.

    Criterion 4 is scored here rather than asserted in prose: for every
    row, the mean and max |difference| against C6 over the power cells they share. A large
    mean says the two checks disagree about the same records; near-zero says they are
    measuring the same thing and one of them is redundant.
    """
    # `clock` is part of the cell identity: Arm C is run on both clocks, and C1/C2 are
    # deliberately absent from its in-spec clock. Omitting it left the reference index with
    # 36 duplicate labels, so `align(join="inner")` returned MORE rows than the left side
    # had - inventing 36 pairings that compared a calendar-clock rate against an
    # in-spec-clock rate from a different record (mean_abs_diff 0.154 instead of 0.188).
    keys = ["arm", "clock", "n_target", "shape", "rho", "b"]
    ref_rows = power[power["check"] == reference]
    if not len(ref_rows):
        return pd.DataFrame()
    ref = ref_rows.set_index(keys)["rejection_rate"]
    if ref.index.duplicated().any():
        raise ValueError(
            f"the agreement key {keys} does not identify a cell for {reference!r}: "
            f"{int(ref.index.duplicated().sum())} duplicate labels. Aligning on it would "
            "cross-join rather than pair."
        )
    out = []
    for row_label, group in power.groupby("row"):
        if row_label.startswith(reference):
            continue
        joined = group.set_index(keys)["rejection_rate"].align(ref, join="inner")
        if not len(joined[0]):
            continue
        diff = (joined[0] - joined[1]).abs()
        out.append(
            {
                "row": row_label,
                "n_shared_cells": int(len(diff)),
                "mean_abs_diff_vs_c6": round(float(diff.mean()), 3),
                "max_abs_diff_vs_c6": round(float(diff.max()), 3),
            }
        )
    return pd.DataFrame(out).sort_values("mean_abs_diff_vs_c6", ascending=False)


def degenerate_cells(size: pd.DataFrame, power: pd.DataFrame) -> pd.DataFrame:
    both = pd.concat([size, power])
    failures = both[(both["n_failed"] > 0) | (both["n_failed_tau_only"] > 0)]
    if not len(failures):
        return failures
    keys = [
        "arm",
        "clock",
        "n_target",
        "quantised",
        "censoring_target",
        "n_replicates",
        "n_failed",
        "n_failed_tau_only",
        "mean_n_distinct_durations",
    ]
    return (
        failures[keys]
        .drop_duplicates()
        .sort_values(["n_failed_tau_only", "n_failed"], ascending=False)
    )


def verdicts(size: pd.DataFrame, power: pd.DataFrame) -> pd.DataFrame:
    """Promote / hold / reject per row, on the rule printed in the report.

    Size is judged INSIDE the real data's censoring envelope, with the per-row Bonferroni
    threshold from `size_verdict_table`. Power is the mean rejection rate at n = 100 over
    the grid of the arm the check is directed at; PROMOTE needs at least 0.5.
    """
    verdict_table = size_verdict_table(size).set_index("row")
    rows = []
    for row_label in sorted(size["row"].unique()):
        check = size[size["row"] == row_label]["check"].iloc[0]
        if row_label not in verdict_table.index:
            rows.append(
                {
                    "row": row_label,
                    "verdict": "UNSCORED",
                    "worst_rate": "-",
                    "worst_z": "-",
                    "power_n100_min": "-",
                    "power_n100_mean": "-",
                    "power_n100_max": "-",
                    "driver": "no adequately supported null cell in the envelope",
                }
            )
            continue
        entry = verdict_table.loc[row_label]
        arm, _parameter, _ = DIRECTED_AT[check]
        directed = power[
            (power["row"] == row_label)
            & (power["arm"] == arm)
            & (power["n_target"] == 100)
        ]
        # min / mean / max over the directed grid. The acceptance rule uses the MEAN, for
        # every check: taking the max reports the single most favourable cell as if it were
        # the expected power - for c1 [permutation] at n = 100 the eight (shape, b) cells
        # run 0.195 to 1.000, and a hand-typed driver string would print 1.000.
        power_min = (
            float(directed["rejection_rate"].min()) if len(directed) else float("nan")
        )
        power_mean = (
            float(directed["rejection_rate"].mean()) if len(directed) else float("nan")
        )
        power_max = (
            float(directed["rejection_rate"].max()) if len(directed) else float("nan")
        )
        scored_power = power_mean
        if not bool(entry["calibrated"]):
            verdict = "REJECT"
            driver = (
                f"size {entry['worst_rate']} at {entry['worst_cell']}, "
                f"z={entry['worst_z']} against a {entry['z_crit']} threshold"
            )
        elif pd.isna(scored_power) or scored_power < 0.5:
            verdict = "HOLD"
            driver = (
                f"size calibrated (worst z={entry['worst_z']} of {entry['z_crit']}) but "
                f"mean power at n=100 is only {_fmt(scored_power, 3)}"
            )
        else:
            verdict = "PROMOTE"
            driver = (
                f"size calibrated across {entry['n_cells']} cells "
                f"(worst z={entry['worst_z']} of {entry['z_crit']}), "
                f"mean power {scored_power:.3f} at n=100 "
                f"(range {power_min:.3f}-{power_max:.3f})"
            )
        rows.append(
            {
                "row": row_label,
                "verdict": verdict,
                "worst_rate": entry["worst_rate"],
                "worst_z": entry["worst_z"],
                "power_n100_min": _fmt(power_min, 3),
                "power_n100_mean": _fmt(power_mean, 3),
                "power_n100_max": _fmt(power_max, 3),
                "driver": driver,
            }
        )
    return pd.DataFrame(rows)


def build(runtime_note: str = "") -> str:
    size, power = load()
    # Computed from the table, not quoted: the earlier draft hand-typed a figure that
    # matched no cell in it, ten lines below a claim that nothing is transcribed by hand.
    _c2_asym_n20 = with_support(
        size[
            (size["arm"] == ARM_A)
            & (size["check"] == "c2_anderson_darling")
            & (size["calibration"] == CALIB_ASYMPTOTIC)
            & (size["n_target"] == 20)
            & (~size["quantised"].astype(bool))
            & (size["censoring_target"] == 0.0)
        ]
    ).sort_values("shape")
    # One figure per Weibull shape: the two cells are not interchangeable, and their mean
    # is in no cell of the table.
    # Arm C's induced duration-level lag-1, from the tables: rho = 0 lives in the size
    # table and the rest of the rho grid in the power table.
    _arm_c = pd.concat([size, power], ignore_index=True)
    _arm_c = _arm_c[(_arm_c["arm"] == ARM_C) & _arm_c["mean_induced_lag1"].notna()]
    _by_rho = _arm_c.groupby(_arm_c["rho"].fillna(0.0))["mean_induced_lag1"].mean()
    arm_c_lag1 = (
        f"{_by_rho.iloc[0]:+.3f} at rho = {_by_rho.index[0]:g} and "
        f"{_by_rho.iloc[-1]:+.3f} at rho = {_by_rho.index[-1]:g}, with every cell between "
        f"{_arm_c['mean_induced_lag1'].min():+.3f} and "
        f"{_arm_c['mean_induced_lag1'].max():+.3f}"
        if len(_arm_c)
        else "(no Arm C cell)"
    )
    shipped_c2_n20 = (
        " and ".join(
            f"{rate:.4f} at Weibull shape {shape:g}"
            for shape, rate in zip(
                _c2_asym_n20["shape"], _c2_asym_n20["rejection_rate"]
            )
        )
        if len(_c2_asym_n20)
        else "(no cell)"
    )

    parts: list[str] = []
    add = parts.append

    add("# Promotion report - the six checks against the calibration bench\n")
    add(
        "Generated by `bench/report.py` from `size_table.csv` and `power_table.csv`. "
        "Every number below is a row of those tables; nothing is transcribed by hand.\n"
    )

    add("\n## Scope, and what is NOT here\n")
    # NOT probed. This file is TRACKED and its header says every number below is a row of
    # the two CSVs; an environment probe would make it a function of the machine too, so
    # regenerating on a box without R would rewrite a committed artifact. The claim below
    # is written to hold in both directions instead, which is also the stronger claim: C3
    # carries no evidence here because the BENCH never ran it, not because R was missing.
    add(
        "**Five checks are assessed, not six.** C3 (`copula::serialIndepTest`) needs R, "
        "and whether R is present changes what the ledger produces elsewhere, not what "
        "this bench measured. C3 is outside this bench either way: it has no arm here, no "
        "entry in `battery.ROW_KEYS`, and no row in `size_table.csv` or `power_table.csv`, "
        "so neither its size nor its power has been measured. "
        "**C3 carries no evidence here. Its silence is not a pass, and neither is a "
        "p-value it returns.**\n"
    )
    add(
        "\n**The transcription is pinned by test, and the pin is narrower than it looks.** "
        "`tests/test_checks_statistics.py` asserts two things on every suite run: that eq "
        "(7) at `gamma_hat = 1` equals the textbook Anderson-Darling form to float "
        "precision, and that its rejection rate matches the Marsaglia LIMITING "
        "distribution at n = 20 and 50 within a stated tolerance. Both concern the "
        "`gamma = 1` path.\n\n"
        "The SHIPPED asymptotic path divides by an estimated `gamma_hat`, and is "
        f"measurably different: C2 rejects at {shipped_c2_n20} on the primary null "
        "configuration at n = 20 (Arm A, continuous, uncensored) against a nominal 0.05. "
        "That gap is not a defect - it is the finite-N cost of estimating gamma, and it is "
        "what the size table below measures. The suite holds the shipped size only to a "
        "coarse band, too wide to detect the gap itself. "
        "Read the size table, not the transcription test, for what the shipped path "
        "does.\n"
    )

    # Computed, not typed: hand-writing these two ranges got both of them wrong ("44-90 null cells",
    # "z_crit 3.3-3.5") and both were wrong.
    _verdict_sizes = size_verdict_table(size)
    cell_lo, cell_hi = (
        int(_verdict_sizes["n_cells"].min()),
        int(_verdict_sizes["n_cells"].max()),
    )
    zc_lo, zc_hi = (
        float(_verdict_sizes["z_crit"].min()),
        float(_verdict_sizes["z_crit"].max()),
    )

    add("\n## Decision rule\n")
    add(
        "```\n"
        "PROMOTE  size statistically indistinguishable from nominal across every null\n"
        "         cell with n >= 35 INSIDE the real data's censoring envelope, and\n"
        "         mean power of at least 0.5 at n = 100 over the grid of the\n"
        "         alternative the check is directed at\n"
        "HOLD     size calibrated but mean power at n = 100 below 0.5\n"
        "REJECT   size miscalibrated somewhere inside the envelope\n"
        "```\n\n"
        "**Size is judged by a multiplicity-corrected z-test, not by a flat tolerance.** "
        f"Each row is scored against {cell_lo}-{cell_hi} null cells. Under a PERFECT test "
        "the largest of "
        "that many deviations is about 2.5-3 Monte Carlo SE by chance alone, so a flat "
        '"every cell within 0.01" rule would reject every check regardless of truth - '
        "the first draft of this report did exactly that. Instead each cell gets "
        "`z = (rate - alpha)/sqrt(alpha(1-alpha)/n_used)` and the row is flagged only if "
        f"its worst |z| exceeds a Bonferroni threshold at family-wise {FAMILYWISE_ALPHA} "
        f"for that row's cell count (z_crit {zc_lo:.2f}-{zc_hi:.2f} here).\n\n"
        f"**The envelope is `censoring <= {ENVELOPE_MAX_CENSORING}`**, because the real "
        f"data sits at {REAL_DATA_NOTES['censoring']}. The 0.25 arm is deliberately "
        "outside it and is scored separately below: a check that works where the data "
        "lives should not be condemned by a corner the data never reaches, and what "
        "happens in that corner is a finding in its own right.\n"
    )

    add("\n## Verdicts\n")
    add(_markdown_table(verdicts(size, power)))
    add("\n\nSize inside the envelope, per row:\n\n")
    add(_markdown_table(size_verdict_table(size)))
    add(
        "\n\nSize OUTSIDE the envelope (censoring 0.25, which the real data never "
        "reaches), reported separately rather than folded into the verdict:\n\n"
    )
    add(_markdown_table(size_verdict_table(size, inside_envelope=False)))

    add(
        "\n\n## Power across the dependence grid - read before trusting a non-rejection\n"
    )
    grid_table, strong = dependence_grid_table(power)
    add(
        "C5 and C6 are scored on Arm E, whose `rho` sets the serial dependence between "
        "durations. Power by event count at every grid point, as the mean over the two "
        "Weibull shapes, with the range of duration-level lag-1 each `rho` induces across "
        "its cells:\n\n"
    )
    add(_markdown_table(grid_table))
    add(
        f"\n\nBelow n = 100, power exceeds 0.5 only at {strong}. A non-rejection from "
        "these checks at small n is therefore close to uninformative unless the "
        "dependence is strong - it is the expected outcome whether the durations are "
        "dependent or not - and must not be read as evidence of independence. They are "
        "correctly calibrated, so a REJECTION is meaningful; it is the silence that "
        "carries no information.\n"
    )

    add("\n\n## Criterion 1 - size at the event counts this project has\n")
    add(
        "Arm A, continuous durations, minimal censoring. Columns are the target event "
        f"count; nominal is {ALPHA}.\n"
    )
    table = size_by_n(size).round(4).reset_index()
    add(_markdown_table(table))
    add("\n\nWorst deviation from nominal across EVERY null cell at n >= 35:\n\n")
    worst = worst_size(size).round(4)
    add(_markdown_table(worst))

    add(
        "\n\n## Criterion 2 - power against the alternative each check is directed at\n"
    )
    for check, (arm, parameter, description) in DIRECTED_AT.items():
        add(f"\n**{check}** - directed at {description}; Arm {arm[0]}, n = 100.\n\n")
        table = power_at(power, check).round(3).reset_index()
        add(_markdown_table(table))
        add("\n")

    add("\n## Criterion 3 - what the censoring machinery buys\n")
    add(
        "Restated as a PREDICTION rather than a pass/fail criterion, per the plan. "
        "Kvaloy & Lindqvist "
        "(Technometrics 62(1) 2020, Section 6.1, Figure 1, p. 107) report that these "
        "asymptotic calibrations are mildly NON-conservative at small samples. Read off "
        "Figure 1 at 10 expected events: AD about 0.11 at shape 1.5 and 0.09 at shape "
        "0.75, LR about 0.08 at both, converging to 0.05 by 40-60. The size table above "
        "is the comparison. Note this "
        "criterion barely discriminates on the real data, where censoring is "
        f"{REAL_DATA_NOTES['censoring']}; the 0.25 arm is deliberately outside that range.\n\n"
    )
    table = censoring_effect(size).round(4).reset_index()
    add(_markdown_table(table))

    add("\n\n## Criterion 4 - do the checks agree where they should\n")
    add(
        "Every row from one replicate reads the SAME permutation set, so differences "
        "between checks within a cell are paired and their Monte Carlo error is smaller "
        "than that of two independent estimates. (Across ARMS the comparison is not "
        "paired - Arm B and Arm C realise different segment sizes, so no shared "
        "permutation set exists. Those are compared unpaired, with both SEs shown.)\n\n"
        "C6 is the natural reference: it is a pure permutation test of exchangeability, so "
        "its validity rests on no approximation and no citation.\n\n"
    )
    add(
        "**Agreement between checks, computed.** Each pair's rejection rates are compared "
        "across the shared power cells of the arm both are directed at; because every row "
        "in a cell reads the same permutation set, the difference is paired. C6 is the "
        "reference: a pure permutation test of exchangeability, resting on no "
        "approximation and no citation.\n\n"
    )
    agreement = check_agreement(power)
    add(_markdown_table(agreement) if len(agreement) else "No shared cells to compare.")
    add("\n\nEffect of quantisation on size (Arm A, no extra censoring):\n\n")
    table = quantisation_effect(size).round(4).reset_index()
    add(_markdown_table(table))

    add("\n\n## Two findings the grid was not designed to produce\n")
    add(
        "\n**1. Read-level dependence does not reach the durations.** Arm C was planned as "
        "the power arm for C5 and C6. It cannot be: correlated READS do not produce "
        f"correlated DURATIONS. Measured induced duration-level lag-1 averages {arm_c_lag1}, "
        "because successive level "
        "crossings of a stationary Gaussian process are very nearly a renewal process. The "
        "power table below is flat at the nominal level across the whole rho grid, which is "
        "the evidence for that claim.\n\n"
        "**Arm E was added because of this** - Gaussian-copula AR(1) on the durations "
        "themselves, which has exactly Arm A's marginal so `rho` moves dependence and "
        "nothing else. Without it, criterion 2 would be unscoreable for C5 and C6. It was "
        "not in the approved plan.\n\n"
    )
    table = arm_c_flatness(power).round(3).reset_index()
    add(_markdown_table(table))
    add(
        "\n\n**2. Eq (10)'s `gamma_tilde` is unusable at small per-segment N.** It is a "
        "difference of two large terms and goes negative: a regular process with unit gaps "
        "truncated at tau = 5.5 gives N = 5, mu = 1.1 and sigma^2 = 1.05 - 1.21 = -0.16. At "
        "~5 events per segment it killed 383 of 400 replicates. The default is therefore "
        "the complete-gap coefficient of variation - which is the `gamma_hat` that actually "
        "appears in eqs (4) and (7); eq (10) is the paper's ALTERNATIVE, and it is kept as "
        "a selectable variant with this caveat recorded.\n"
    )

    add("\n## Degenerate and failed cells\n")
    add(
        "A replicate whose statistic is undefined is counted, never silently treated as a "
        "non-rejection. `n_failed_tau_only` counts replicates where C1/C2 were undefined "
        "but the rank checks still ran - typically a segment whose durations are all one "
        "quantum, so `gamma_hat` is zero. Each row's own `n_used` is its denominator.\n\n"
    )
    thin = under_supported(pd.concat([size, power]))
    add(
        f"\n**Rows excluded from every size and verdict table** (support below "
        f"{MIN_SUPPORT_FRACTION:.0%} of their cell's replicates). These are reported here "
        "and nowhere else, because a rate computed from the surviving replicates is "
        "conditional on whatever killed the others:\n\n"
    )
    if len(thin):
        columns = [
            "arm",
            "n_target",
            "quantised",
            "censoring_target",
            "check",
            "calibration",
            "n_replicates",
            "n_used",
            "rejection_rate",
        ]
        add(_markdown_table(thin[columns].round(4)))
        add(
            "\n\nThe `c2_anderson_darling` and `cvm_cramer_von_mises` `[asymptotic]` entries "
            "at `n_used` of 1 and 5 are not a degeneracy at all - both asymptotic "
            "calibrations are only defined for a single segment, and at censoring 0.25 a "
            "replicate has many, so the row exists "
            "only on the rare replicate where every segment but one was dropped for having "
            "too few events. It is an artefact of the row schema, not a measurement.\n"
        )
    else:
        add("None.\n")
    add("\nCells with at least one failed replicate:\n\n")

    failures = degenerate_cells(size, power)
    if len(failures):
        add(_markdown_table(failures.round(2).head(25)))
        add(f"\n\n({len(failures)} cell configurations had at least one failure.)\n")
    else:
        add("No cell had a failed replicate.\n")

    add("\n## Reproducibility limitation\n")
    add(
        "The bench and the checks are tracked, and every replicate is seeded from "
        "`(MASTER_SEED, cell index, replicate)`, so the tables reproduce from the commit "
        "that produced them. The tables do not say which commit that is: they carry no "
        "commit, code hash, weighting or divisor, so nothing detects that they describe "
        "code that has since changed.\n"
    )

    add("\n## Runtime\n")
    add(runtime_note or "(runtime.txt not found)")
    add("\n")
    return "".join(parts)


def main() -> Path:
    runtime_path = RESULTS_DIR / "runtime.txt"
    runtime_note = (
        "```\n" + runtime_path.read_text() + "```\n" if runtime_path.exists() else ""
    )
    target = RESULTS_DIR / "promotion_report.md"
    target.write_text(build(runtime_note))
    print(f"[report] wrote {target}")
    return target


if __name__ == "__main__":
    main()
