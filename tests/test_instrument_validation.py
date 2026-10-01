"""The instrument report must not overstate the evidence it summarises.

That is the whole risk of this module. It is a REPORT: its job is to say what is known, and
the failure mode is not a wrong number but a verdict that claims more than the tier behind
it supports. So most of these tests are about honesty rather than arithmetic.

The arithmetic that IS here is the Tier 2 recomputation, which must agree with
`test_checks_published_values.py`. Two files computing the same published comparison from
the same record is deliberate: if they ever disagree, one of them is reading the reference
data wrongly and the report is the one that reaches a figure.
"""

from __future__ import annotations

import functools

import numpy as np
import pandas as pd
import pytest

from quebra.analyzers.instrument_validation import (
    TIER_ABSENT,
    TIER_FAIL,
    TIER_PARTIAL,
    TIER_PASS,
    build_instrument_validation,
    build_published_comparisons,
    render_tier_table_markdown,
)
from tests.fixtures import R_REFERENCE_INPUTS, R_REFERENCE_VALUES

pytestmark = pytest.mark.statistical


REFERENCE = R_REFERENCE_VALUES.parent
GAPS = pd.read_csv(REFERENCE / "load_haul_dump.csv")
PUBLISHED = pd.read_csv(REFERENCE / "load_haul_dump_published.csv")
R_INPUTS = pd.read_csv(R_REFERENCE_INPUTS)
R_VALUES = pd.read_csv(R_REFERENCE_VALUES)


def _tie_frame() -> pd.DataFrame:
    """A minimal tie table. Built here rather than read, so these tests never wait on the
    ~1 hour bench study and never silently pass because its CSV is stale."""
    rows = []
    for n in (35, 355):
        for tie, signed in [(0.0, 0.001), (0.5, 0.004), (0.9, 0.030), (0.99, 0.045)]:
            rows.append(
                {
                    "n": n,
                    "y_levels": 0,
                    "tie_fraction_y": tie,
                    "p_signed_diff_mean": signed,
                    "p_signed_diff_se": 0.001,
                    "p_abs_diff_mean": abs(signed) + 0.01,
                    "type_i_closed": 0.05 + signed,
                    "type_i_perm": 0.05,
                    "xi_tie_break_sd": tie * 0.07,
                }
            )
    return pd.DataFrame(rows)


@functools.lru_cache(maxsize=1)
def _build():
    """Memoised: about 2 s of Monte Carlo per construction, and every test wants the same one."""
    return build_instrument_validation(
        GAPS, PUBLISHED, R_INPUTS, R_VALUES, _tie_frame()
    )


# ------------------------------------------------------------------- Tier 2 arithmetic


def test_the_published_comparison_agrees_with_the_tier_2_test_module():
    """Two independent readers of the same published record must not disagree."""
    from tests.fixtures.load_haul_dump import PUBLISHED as PUB_DICT

    rows = {r.quantity: r for r in build_published_comparisons(GAPS, PUBLISHED)}
    assert rows["mu_hat"].ours == pytest.approx(PUB_DICT["mu_hat"], abs=0.01)
    assert rows["gamma_tilde (eq 10)"].ours == pytest.approx(
        PUB_DICT["gamma_tilde"], abs=0.001
    )
    assert rows["Laplace"].ours == pytest.approx(PUB_DICT["laplace"], abs=0.001)
    # The complete-gap rows use the shipped default; it is the sample divisor, so they are
    # direct pins against Table 2 and the Section 8.1 text rather than a rescaled bridge.
    assert rows["gamma_hat"].ours == pytest.approx(PUB_DICT["gamma_hat"], abs=0.0005)
    assert rows["LR"].ours == pytest.approx(PUB_DICT["lr_gamma_hat"], abs=0.0005)
    assert rows["gamma_hat"].agrees and rows["LR"].agrees


def test_the_reference_csv_holds_the_published_record_itself():
    """The figure reads `reference/`, the suite reads the module. They must be one record."""
    from tests.fixtures.load_haul_dump import gaps_h

    assert np.allclose(GAPS["gap_h"].to_numpy(dtype=float), gaps_h())


def test_a_quantity_cannot_be_marked_as_agreeing_when_it_does_not():
    """The positive control on the artifact's own honesty.

    `build_published_comparisons` raises rather than emitting `agrees=True` beside numbers
    that differ. Without this guard the figure would print a green bar over a disagreement,
    which is the one output this whole pass exists to prevent.
    """
    corrupted = PUBLISHED.copy()
    corrupted.loc[corrupted["quantity"] == "mu_hat", "value"] = 99.0
    with pytest.raises(ValueError, match="marked as agreeing"):
        build_published_comparisons(GAPS, corrupted)


def test_the_c1_tier_2_verdict_follows_the_published_rows(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Oracle: the tier rule, `pass` only when every published quantity agrees.

    The published record agrees on every row and must read `pass`; the same build with one
    row flipped must read `partial`. A verdict typed in, or computed from anything but the
    rows, fails one of the two.
    """
    import dataclasses

    import quebra.analyzers.instrument_validation as iv

    assert _build().tier_verdict("C1 Lewis-Robinson", 2).verdict == TIER_PASS
    rows = build_published_comparisons(GAPS, PUBLISHED)
    rows[0] = dataclasses.replace(rows[0], agrees=False)
    monkeypatch.setattr(iv, "build_published_comparisons", lambda _g, _p: rows)
    flipped = build_instrument_validation(
        GAPS, PUBLISHED, R_INPUTS, R_VALUES, _tie_frame()
    )
    assert flipped.tier_verdict("C1 Lewis-Robinson", 2).verdict == TIER_PARTIAL


def test_a_population_default_is_reconciled_to_table_2_by_the_divisor_bridge(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Oracle: Table 2's printed values, reached from the 1/N estimates by sqrt(N/(N-1)).

    With `GAMMA_DEFAULT` switched to the population form, sigma_hat, gamma_hat and LR must
    disagree, each note must name the bridge, and the bridge must carry our value onto the
    printed one. The figure's caption must then say that rows differ.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import quebra.analyzers.checks._multiprocess as mp
    from quebra.plots.instrument_validation_plot import PublishedValuesPlot

    monkeypatch.setattr(mp, "GAMMA_DEFAULT", mp.GAMMA_COMPLETE)
    rows = {r.quantity: r for r in build_published_comparisons(GAPS, PUBLISHED)}
    n = len(GAPS)
    bridge = (n / (n - 1)) ** 0.5
    for quantity, carried in [
        ("sigma_hat", rows["sigma_hat"].ours * bridge),
        ("gamma_hat", rows["gamma_hat"].ours * bridge),
        ("LR", rows["LR"].ours / bridge),
    ]:
        row = rows[quantity]
        assert not row.agrees, quantity
        assert "1/N vs their 1/(N-1)" in row.note and "reconciles" in row.note, row.note
        assert carried == pytest.approx(row.published, rel=1e-3), quantity
    data = build_instrument_validation(
        GAPS, PUBLISHED, R_INPUTS, R_VALUES, _tie_frame()
    )
    fig = PublishedValuesPlot(name="unit").build_matplotlib(data)
    try:
        caption = " ".join(text.get_text() for text in fig.texts)
    finally:
        plt.close(fig)
    assert "3 differ, each row's note says why." in caption, caption


# ------------------------------------------------------------------------- honesty


def test_every_tier_verdict_is_from_the_known_vocabulary():
    known = {TIER_PASS, TIER_FAIL, TIER_PARTIAL, TIER_ABSENT}
    for row in _build().tiers:
        assert row.verdict in known
        assert row.detail.strip(), (
            f"{row.instrument} tier {row.tier} has no justification"
        )


def test_c3_is_not_credited_with_calibration_it_does_not_have():
    """C3 runs and has NO bench cell.

    A smoke test on iid input is not calibration, and this pass has already had to correct
    five documents that implied otherwise. The report must not become the sixth.
    """
    data = _build()
    tier3 = data.tier_verdict("C3 serial copula", 3)
    assert tier3 is not None and tier3.verdict == TIER_ABSENT
    assert "no bench cell" in tier3.detail


def test_xi_tier_3_is_partial_because_the_closed_form_is_not_entitled_to_ties():
    data = _build()
    row = data.tier_verdict("Chatterjee xi", 3)
    assert row is not None and row.verdict == TIER_PARTIAL


def test_cvm_is_promoted_and_its_tier_3_no_longer_says_it_has_no_bench_cell():
    """Once the bench carries CvM rows, "no bench cell exists for CvM" is false, and a
    report repeating it understates the evidence it holds.
    """
    from quebra.analyzers.checks.battery import ROW_KEYS

    assert any(key[0] == "cvm_cramer_von_mises" for key in ROW_KEYS)
    row = _build().tier_verdict("CvM", 3)
    assert row is not None
    assert "no bench cell exists for CvM" not in row.detail


def test_the_markdown_says_it_is_generated_and_disclaims_power():
    text = render_tier_table_markdown(_build())
    assert "GENERATED by" in text
    assert "Do not hand-edit" in text
    assert "None of this is power" in text
    # Every instrument reaches the table.
    for instrument in _build().instruments:
        assert instrument in text


def test_the_markdown_reports_a_never_crossing_row_as_a_result_not_a_blank():
    """A tie fraction that never crosses the threshold is a FINDING - the closed form held.

    Rendering it as an empty cell would read as missing data, which is the opposite of what
    it means.
    """
    frame = _tie_frame()
    frame["p_signed_diff_mean"] = 0.0001  # nothing crosses
    data = build_instrument_validation(GAPS, PUBLISHED, R_INPUTS, R_VALUES, frame)
    assert all(not np.isfinite(v) for v in data.divergence_tie_fraction.values())
    assert all(not np.isfinite(v) for v in data.divergence_levels.values())
    text = render_tier_table_markdown(data)
    assert "never diverges" in text
    # And it must not render as an empty cell, which would read as missing data.
    assert "|  |" not in text


# --------------------------------------------------------------------- tier 4 rows


def test_the_random_reference_row_is_flagged_as_random():
    """The tied-x case compares against a DISTRIBUTION; the artifact must carry that fact,
    or the figure would draw an equality that does not exist."""
    rows = _build().cross_implementation
    xi_rows = [r for r in rows if r.case == "xi_tied" and r.statistic == "xi"]
    assert len(xi_rows) == 1
    assert xi_rows[0].reference_is_random
    assert xi_rows[0].reference_spread > 0.0

    # dcor on the SAME tied input is NOT random and must not be flagged as such: it has no
    # tie-breaking step at all, so it agrees exactly. Having both statistics on one case is
    # what makes this a test of the flag rather than of the case.
    dcor_rows = [r for r in rows if r.case == "xi_tied" and r.statistic == "dcor"]
    assert len(dcor_rows) == 1
    assert not dcor_rows[0].reference_is_random
    assert dcor_rows[0].abs_difference < 1e-9


def test_the_exact_reference_rows_actually_agree():
    for row in _build().cross_implementation:
        if not row.reference_is_random:
            assert row.abs_difference < 1e-9, f"{row.case}/{row.statistic} disagrees"


def test_the_divergence_is_reported_in_levels_not_only_tie_fraction():
    """Tie fraction saturates; levels is the actionable form.

    Quantising to 20 levels already ties 83-99% of a sample, so "diverges at tie fraction
    1.000" is true of several different cells and a reader cannot check it against their
    own metric. "Diverges once the response has 2 distinct values" they can.
    """
    frame = _tie_frame()
    frame["y_levels"] = [20, 10, 3, 2] * 2
    frame["p_signed_diff_mean"] = [0.001, 0.002, 0.003, 0.050] * 2
    data = build_instrument_validation(GAPS, PUBLISHED, R_INPUTS, R_VALUES, frame)
    assert data.divergence_levels == {35: 2.0, 355: 2.0}


def test_divergence_levels_picks_the_FINEST_crossing_not_the_coarsest():
    """If 2 and 3 levels both cross, the answer is 3 - that is where it starts."""
    frame = _tie_frame()
    frame["y_levels"] = [20, 10, 3, 2] * 2
    frame["p_signed_diff_mean"] = [0.001, 0.002, 0.030, 0.050] * 2
    data = build_instrument_validation(GAPS, PUBLISHED, R_INPUTS, R_VALUES, frame)
    assert data.divergence_levels == {35: 3.0, 355: 3.0}


def test_no_tier_verdict_claims_evidence_that_does_not_exist():
    """The audit, kept as a test.

    Four verdicts claimed more than their tier supported: two tier-4 `pass` cells with no
    independent implementation behind them (C2's reference is five hand-written lines in
    our own test file; C3's was 'it IS the R implementation', which compares R to itself),
    and tier-3 cells asserting a direction that the measurement contradicts.

    C3's tier-4 cell is now `pass` on real evidence: an `r`-marked test puts the bridge
    against a committed fixture on both the statistic and the p-value. The assertion below
    moved with it, and what it guards moved too. It no longer asks whether the cell is
    `absent`; it asks that a `pass` there NAME the artifact it rests on, so the circular
    claim this audit removed cannot return under the same verdict.
    """
    data = _build()
    c2_t4 = data.tier_verdict("C2 Anderson-Darling", 4)
    assert c2_t4 is not None and c2_t4.verdict != TIER_PASS, (
        "a reimplementation in the same language by the same author is tier-1 evidence"
    )
    c5_t4 = data.tier_verdict("C5 rank autocorr", 4)
    assert c5_t4 is not None and c5_t4.verdict != TIER_PASS
    assert "matches R" not in c5_t4.detail, (
        "no test compares C5's own statistic with R: the lag-1 check computes a Pearson "
        "correlation inline, and C5 divides by the full centred sum of squares instead"
    )
    c3_t4 = data.tier_verdict("C3 serial copula", 4)
    assert c3_t4 is not None
    assert c3_t4.verdict != TIER_PASS, (
        "tier 4 is agreement with an INDEPENDENT implementation of the same statistic. The "
        "fixture and the bridge both call copula::serialIndepTest, so no evidence here can "
        "be tier-4 pass; bridge fidelity is a weaker claim and belongs at partial"
    )
    if c3_t4.verdict != TIER_ABSENT:
        assert "r_reference_values.csv" in c3_t4.detail, (
            "a non-absent cell must name the artifact it rests on; the claim this audit "
            "removed was 'it IS the R implementation', which named nothing"
        )
        assert "test_r_cross_implementation.py" in c3_t4.detail, (
            "and it must name the test, so the claim cannot outlive its evidence"
        )

    # No tier-3 row may assert a DIRECTION: measured, it flips between Weibull and
    # exponential gaps, so any single direction is wrong for one of them.
    for row in data.tiers:
        if row.tier == 3:
            assert (
                "anti-conservative at small n as the source reports" not in row.detail
            )


def test_the_tier_3_rows_quote_a_number_this_run_produced():
    """`measure_all_asymptotic_sizes` must actually be called, not merely cited.

    The defect this guards: the measurement lived in `tests/`, three docstrings claimed the
    report imported it, and `grep` matched only the test file - so 'asymptotic size measured'
    rested on no number in the artifact at all.
    """
    import re

    data = _build()
    # Driven from the artifact, not from a literal instrument tuple: a hardcoded list is
    # outrun by the next instrument added, which is what happened when Kaplan-Meier landed
    # with a tier-3 row this loop could not see.
    #
    # Every tier-3 row, not only the passing ones: C1 and C2 are `partial` and still quote
    # a measured size, and a number that decides a partial verdict has to be as real as one
    # that decides a pass.
    #
    # Not every row quotes a number - C5 and C6 pass on permutation exactness, which is an
    # identity rather than a measurement. The rule is therefore: any row that CLAIMS a
    # measured number must be re-derivable here, and a row that says "measured" in a shape
    # this guard cannot parse is itself a failure.
    from quebra.analyzers.instrument_validation import (
        BAND_COVERAGE_REPLICATES,
        measure_all_asymptotic_sizes,
        measure_band_coverage,
    )

    remeasured: set[str] = set()
    for row in [r for r in data.tiers if r.tier == 3]:
        size_match = re.search(r"size measured at tau=([\d.]+): (0\.\d{4})", row.detail)
        band_match = re.search(
            r"band coverage measured at t=[\d.]+: (0\.\d{4})", row.detail
        )
        if size_match:
            expected = measure_all_asymptotic_sizes(tau=float(size_match.group(1)))[
                row.instrument
            ]
            assert float(size_match.group(2)) == pytest.approx(expected, abs=5e-5), (
                f"{row.instrument} quotes {size_match.group(2)} but measures {expected:.4f}"
            )
            remeasured.add(row.instrument)
        elif band_match:
            expected_cov = measure_band_coverage(replicates=BAND_COVERAGE_REPLICATES)
            assert float(band_match.group(1)) == pytest.approx(
                expected_cov, abs=5e-5
            ), (
                f"{row.instrument} quotes {band_match.group(1)} but "
                f"measures {expected_cov:.4f}"
            )
            remeasured.add(row.instrument)
        else:
            assert "measured" not in row.detail, (
                f"{row.instrument} tier 3 says 'measured' but quotes no number this guard "
                f"can re-derive: {row.detail!r}"
            )

    assert "Kaplan-Meier" in remeasured, "the KM tier-3 row must be re-measured here"
    assert remeasured == {
        "C1 Lewis-Robinson",
        "C2 Anderson-Darling",
        "CvM",
        "Kaplan-Meier",
    }, f"re-measured {sorted(remeasured)}; a row stopped quoting its number"


def test_the_caption_and_report_carry_the_record_they_were_built_from(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Oracle: a record with known counts, made by keeping 30 published gaps and ending at them.

    AGENTS.md section 4: a measured number in prose comes from the artifact it cites, in the
    state it ships. The load-haul-dump caption and the report's tier-2 heading quote the
    record's gap count, censoring time and censored-gap count. The published record must give
    36, 2000 h and 1; a shortened one, 30 gaps observed to their own sum, must give 30, that
    sum, and none. A number fixed in the text reads the same on both, so it fails one of them.

    Tier 2 compares against the published values, and its agreement flags are fixed per
    quantity, so a shortened record cannot pass it: its rows are taken from the published
    record, and every other part of the build runs on the shortened one.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import quebra.analyzers.instrument_validation as iv
    from quebra.plots.instrument_validation_plot import PublishedValuesPlot

    import dataclasses

    # One agreeing row flipped, so tier 2's counts differ from the published record's.
    published_rows = iv.build_published_comparisons(GAPS, PUBLISHED)
    first = next(i for i, row in enumerate(published_rows) if row.agrees)
    published_rows[first] = dataclasses.replace(published_rows[first], agrees=False)
    n_rows, n_exact = len(published_rows), sum(r.agrees for r in published_rows)
    monkeypatch.setattr(
        iv, "build_published_comparisons", lambda _g, _p: published_rows
    )
    # The C3 bridge's seed moved too, so a seed typed into its tier-4 row would disagree.
    r_values = R_VALUES.copy()
    r_values.loc[r_values["quantity"] == "serial_indep_sim_seed", "value"] = 708
    fewer = GAPS.iloc[:30]
    tau_h = float(fewer["gap_h"].sum())
    ending_on_a_failure = PUBLISHED.copy()
    ending_on_a_failure.loc[ending_on_a_failure["quantity"] == "tau_h", "value"] = tau_h
    data = build_instrument_validation(
        fewer, ending_on_a_failure, R_INPUTS, r_values, _tie_frame()
    )
    tier2 = data.tier_verdict("C1 Lewis-Robinson", 2).detail
    assert tier2.startswith(
        f"{n_exact} of {n_rows} published quantities exact; {n_rows - n_exact} differ"
    ), tier2
    assert "seed 708," in data.tier_verdict("C3 serial copula", 4).detail
    assert data.meta["load_haul_dump_complete_gaps"] == 30
    assert data.meta["load_haul_dump_tau_h"] == tau_h
    assert data.dropped == {"load_haul_dump_censored_gaps": 0}

    fig = PublishedValuesPlot(name="unit").build_matplotlib(data)
    try:
        caption = " ".join(text.get_text() for text in fig.texts)
    finally:
        plt.close(fig)
    expected = (
        f"30 complete gaps, observed to {tau_h:.10g} h, ending on a failure, "
        "so no gap is censored."
    )
    assert expected in caption, caption

    report = render_tier_table_markdown(data)
    assert "load-haul-dump record, 30 complete gaps," in report
    assert f"observed to {tau_h:.10g} h, ending on a failure." in report

    real = _build()
    assert real.meta["load_haul_dump_complete_gaps"] == len(GAPS) == 36
    assert real.meta["load_haul_dump_tau_h"] == 2000.0
    assert real.dropped == {"load_haul_dump_censored_gaps": 1}
    fig = PublishedValuesPlot(name="unit").build_matplotlib(real)
    try:
        caption = " ".join(text.get_text() for text in fig.texts)
    finally:
        plt.close(fig)
    assert (
        "36 complete gaps, time censored at 2000 h, 1 censored gap dropped." in caption
    )

    # A record whose gaps run past its own censoring time is not a record: refuse it.
    past = PUBLISHED.copy()
    past.loc[past["quantity"] == "tau_h", "value"] = 1000.0
    with pytest.raises(ValueError, match="past the censoring time"):
        iv.load_haul_dump_record(GAPS, past)


def test_every_citation_in_the_report_resolves_and_the_bench_claim_holds():
    """Oracle: the test files on disk, and the committed size table the report cites.

    The report cites instead of copying. A pointer naming a test that no longer exists, or a
    bench claim its own cells no longer bear out, is a copied number's failure in a new form:
    text that outlives its evidence.
    """
    import re
    from pathlib import Path

    from quebra.analyzers.instrument_validation import SIZE_TABLE_PATH

    repo = Path(__file__).resolve().parents[1]
    text = render_tier_table_markdown(_build())

    # The tie paragraph is computed from the tie table it was given. `_tie_frame()` has a
    # tie-free |dp| of 0.0110 at every n and no quantised row, unlike the bench's table.
    assert "`|dp|` still reads 0.0110 across n." in text, text
    assert "saturates as the response coarsens" in text

    pointers = set(re.findall(r"(tests/[\w/]+\.py)::(\w+)", text))
    assert len(pointers) >= 10, (
        f"the report cites too few tests to be the new report: {pointers}"
    )
    for path, name in sorted(pointers):
        source = (repo / path).read_text(encoding="utf-8")
        assert re.search(rf"^def {name}\(", source, re.MULTILINE), (
            f"{path}::{name} is gone"
        )

    table = pd.read_csv(repo / SIZE_TABLE_PATH)
    key = re.compile(
        r"check=(\w+), calibration=asymptotic, arm=(\w+), clock=(\w+), "
        r"quantised=(\w+), censoring_target=([\d.]+), n_target=(\d+)"
    )
    cited = key.findall(text)
    assert sorted({c[0] for c in cited}) == [
        "c1_lewis_robinson",
        "c2_anderson_darling",
        "cvm_cramer_von_mises",
    ], cited
    for check, arm, clock, quantised, censoring, n_target in set(cited):
        cells = table[
            (table["kind"] == "size")
            & (table["check"] == check)
            & (table["calibration"] == "asymptotic")
            & (table["arm"] == arm)
            & (table["clock"] == clock)
            & (table["quantised"].astype(str) == quantised)
            & (table["censoring_target"] == float(censoring))
            & (table["n_target"] == int(n_target))
        ]
        assert sorted(cells["shape"]) == [0.75, 1.5], (check, cells["shape"].tolist())
        assert (cells["rejection_rate"] > 0.05).all(), (
            f"the report says {check}'s cited asymptotic cells have point estimates above "
            f"nominal; the table says {cells['rejection_rate'].tolist()}"
        )

    # The permutation half is the promotion report's own verdict, inside its envelope, and the
    # envelope the report states is the one the promotion report defines.
    promotion_text = (repo / "jobs/bench/results/promotion_report.md").read_text()
    assert "**The envelope is `censoring <= 0.03`**" in promotion_text
    assert "cell with n >= 35 INSIDE" in promotion_text
    assert "hold size at n >= 35 and censoring <= 0.03" in text
    promotion = promotion_text.splitlines()
    start = promotion.index("Size inside the envelope, per row:")
    envelope = {}
    for line in promotion[start + 3 :]:
        if not line.startswith("|"):
            break
        cells = [c.strip() for c in line.strip("|").split("|")]
        envelope[cells[0]] = cells[5]
    for check in {c[0] for c in cited}:
        assert f"check={check}, calibration=permutation" in text
        assert envelope[f"{check} [permutation]"] == "True", (check, envelope)

    from quebra.analyzers.instrument_validation import _span

    assert _span(pd.Series([0.25, 0.5]), ".2f") == "0.25 to 0.50"
    assert _span(pd.Series([0.25, 0.25]), ".2f") == "0.25"
