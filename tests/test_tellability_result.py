"""The tellability step as a whole: its status vocabulary and its wiring.

Oracles:
- specification: every status is `computed`, `not declared` or `not computed: <reason>` with
  a reason from the fixed set; a group's fields are all set when it is computed and all None
  when it is not, and no computed number, inside a tuple or a panel included, is NaN;
- the row functions themselves: `run` must report what each row's function gives on the same
  observed reads, threshold by threshold, on both sides of a threshold, at the gap rule and
  confidence level it is given;
- specification for sigma: a metric whose sigma is not declared, known but not usable, or
  declared with no finite value carries that status on every row that reads sigma;
- the load's own counts, which the result repeats.
"""

from __future__ import annotations

import dataclasses
import io
import math
import pickle

import numpy as np
import pytest

from quebra.analyzers import tellability as tl
from quebra.analyzers import windows
from quebra.core._artifact_guard import load_artifact

# One threshold is small-is-good, so every run reads both sides of the in-spec definition.
THRESHOLDS = (("u0", 0.0, True), ("u1", 0.8, True), ("d0", 0.3, False))


def _record(seed: int = 21, n: int = 1500):
    """A drifting signal read at irregular times, with one gap; sigma known and varying."""
    rng = np.random.default_rng(seed)
    steps = rng.uniform(10.0, 40.0, size=n - 1)
    steps[n // 2] += 20000.0
    t = np.r_[0.0, np.cumsum(steps)]
    sigma = rng.uniform(0.2, 0.6, size=n)
    q = np.sin(2.0 * np.pi * t / 4000.0) + rng.normal(scale=sigma)
    return t, q, sigma


def _carve(t, q, thresholds, gap_mult=10.0, k=1.0, dataset_id="synthetic"):
    return windows.run(
        windows.WindowsInputs(
            t_rel_s=t,
            values=q,
            thresholds=list(thresholds),
            gap_mult=gap_mult,
            k=k,
            dataset_id=dataset_id,
        )
    )


def _inputs(
    t,
    q,
    sigma,
    *,
    sigma_status="declared",
    thresholds=THRESHOLDS,
    gap_mult=10.0,
    carve_gap_mult=None,
    unsure_k=1.0,
    clip_floor=None,
    windows_result="carve",
    n_permutations=19,
    seed=5,
    load_counts=None,
):
    if windows_result == "carve":
        windows_result = (
            _carve(
                t,
                q,
                thresholds,
                gap_mult if carve_gap_mult is None else carve_gap_mult,
                unsure_k,
            )
            if thresholds
            else None
        )
    return tl.TellabilityInputs(
        dataset_id="synthetic",
        metric_key="metric",
        value_key="q",
        t_rel_s=t,
        values=q,
        sigma=sigma,
        sigma_status=sigma_status,
        thresholds=tuple(thresholds),
        windows=windows_result,
        unsure_k=unsure_k,
        gap_mult=gap_mult,
        clip_floor=clip_floor,
        n_permutations=n_permutations,
        seed=seed,
        load_counts={} if load_counts is None else load_counts,
    )


def _all_finite(value) -> bool:
    """Every number inside `value` is finite; None inside a tuple is an absent entry."""
    if isinstance(value, (str, bool)):
        return True
    if isinstance(value, (int, float, np.integer, np.floating)):
        return math.isfinite(float(value))
    if isinstance(value, np.ndarray):
        return bool(np.all(np.isfinite(value)))
    if isinstance(value, tuple):
        return all(v is None or _all_finite(v) for v in value)
    if isinstance(value, dict):
        return all(_all_finite(v) for v in value.values())
    if dataclasses.is_dataclass(value):
        return all(
            _all_finite(getattr(value, f.name)) for f in dataclasses.fields(value)
        )
    raise TypeError(type(value))


def _check_groups(row, groups) -> None:
    valid_status = {tl.COMPUTED, tl.NOT_DECLARED}
    assert set(row.status) == set(groups)
    for group, fields in groups.items():
        status = row.status[group]
        if status not in valid_status:
            assert status.startswith("not computed: "), status
            tl.check_reason(status[len("not computed: ") :])
        values = [getattr(row, name) for name in fields]
        if status == tl.COMPUTED:
            for name, value in zip(fields, values):
                assert value is not None, (group, name)
                assert _all_finite(value), (group, name, value)
        else:
            assert all(value is None for value in values), (group, status)


@pytest.mark.unit
@pytest.mark.parametrize("record", ["full", "degenerate"])
def test_every_status_comes_from_the_vocabulary_and_matches_its_fields(record) -> None:
    """Oracle: specification. On a full record every group is computed, at every threshold;
    on a constant, short one most are not, each with its reason, and none ships NaN as
    computed."""
    if record == "full":
        t, q, sigma = _record()
        result = tl.run(_inputs(t, q, sigma))
        assert set(result.metric.status.values()) - {tl.COMPUTED} == {
            "not computed: no floor for this metric"
        }
        for row in result.thresholds:
            assert set(row.status.values()) == {tl.COMPUTED}, row.threshold_label
    else:
        t = np.arange(8) * 10.0
        q = np.full(8, 0.5)
        result = tl.run(_inputs(t, q, np.full(8, 0.1), thresholds=(("u0", 0.0, True),)))
        assert (
            result.metric.status["second_difference"] == "not computed: zero variance"
        )
        assert result.thresholds[0].status["count_exponent"] == (
            "not computed: no down-crossing at this threshold"
        )
    _check_groups(result.metric, tl.METRIC_GROUPS)
    for row in result.thresholds:
        _check_groups(row, tl.THRESHOLD_GROUPS)


def _degenerate(shape: str):
    rng = np.random.default_rng(40)
    n = 400
    t = np.r_[0.0, np.cumsum(rng.uniform(10.0, 40.0, size=n - 1))]
    if shape == "gap":
        t[n // 2 :] += 20000.0
    sigma = (
        np.full(n, 0.3) if shape == "constant-sigma" else rng.uniform(0.2, 0.6, size=n)
    )
    if shape == "means-equal":
        # Every time read twice, at 1 + d and 1 - d (d dyadic, so each mean is exactly 1):
        # the reads vary, the averaged reads do not.
        times = t[: n // 2]
        d = rng.integers(1, 500, size=times.size) / 1024.0
        return np.repeat(times, 2), np.ravel(np.column_stack([1.0 + d, 1.0 - d])), sigma
    q = {
        "constant": np.full(n, 0.5),
        "linear": 0.003 * t,
        "constant-inexact": np.full(n, 0.3),
        "constant-sigma": rng.normal(size=n),
        "half-constant": np.r_[np.zeros(n // 2), rng.normal(size=n - n // 2)],
        "all-in": 5.0 + rng.normal(scale=0.1, size=n),
        "all-out": -5.0 + rng.normal(scale=0.1, size=n),
        "gap": rng.normal(size=n),
        "no-gap": rng.normal(size=n),
    }[shape]
    return t, q, sigma


@pytest.mark.unit
@pytest.mark.parametrize(
    "shape",
    [
        "constant",
        "constant-inexact",
        "constant-sigma",
        "half-constant",
        "all-in",
        "all-out",
        "gap",
        "no-gap",
        "two",
        "one",
        "forty-one",
        "linear",
        "means-equal",
    ],
)
def test_a_legal_record_of_any_shape_ships_a_result(shape) -> None:
    """Oracle: specification. A record the rows cannot fully answer is not an error: `run`
    returns, and every group carries a status from the vocabulary. A row function that let
    a plain error escape on a legal record would raise here."""
    if shape == "two":
        t, q, sigma = np.array([0.0, 10.0]), np.array([0.1, 0.9]), np.array([0.2, 0.2])
    elif shape == "one":
        t, q, sigma = np.array([0.0]), np.array([0.1]), np.array([0.2])
    elif shape == "forty-one":
        # Halves of 21 and 20 reads: the first holds exactly ten grid frequencies, and at
        # this spacing rounding puts f_top / 10 a bit below f_0.
        rng = np.random.default_rng(41)
        t = np.arange(41) * 42.8899947017729
        q, sigma = rng.normal(size=41), rng.uniform(0.2, 0.6, size=41)
    else:
        t, q, sigma = _degenerate(shape)
    result = tl.run(_inputs(t, q, sigma))
    _check_groups(result.metric, tl.METRIC_GROUPS)
    for row in result.thresholds:
        _check_groups(row, tl.THRESHOLD_GROUPS)
    zero = "not computed: zero variance"
    if shape == "constant-inexact":
        # 400 reads of 0.3 have an sd of rounding, about 6e-17: zero range is the test.
        for group in (
            "spectrum",
            "second_difference",
            "self_consistency",
            "sigma_bound",
        ):
            assert result.metric.status[group] == zero, group
    if shape == "constant-sigma":
        assert result.metric.status["dependence"] == zero
    if shape == "linear":
        # Exactly linear values: sigma_hat_D2 is rounding, so the sigma bound is zero.
        assert result.metric.status["sigma_bound"] == zero
    if shape == "means-equal":
        assert result.metric.status["spectrum"] == zero
        assert result.metric.status["g10_halves"] == zero


def _fields(groups, status) -> set[str]:
    return {
        name
        for group, names in groups.items()
        if status[group] == tl.COMPUTED
        for name in names
    }


def _assert_fields(row, expected) -> None:
    for name, value in expected.items():
        got = getattr(row, name)
        if isinstance(value, np.ndarray):
            assert np.array_equal(got, value), name
        else:
            assert got == value, name


@pytest.mark.unit
def test_run_reports_what_the_row_functions_give() -> None:
    """Oracle: the row functions, called directly on the same observed reads, for EVERY
    field of every computed group, threshold by threshold and in order, on both sides of
    the in-spec definition. It sees what `run` adds (which value lands in which field,
    order, side, settings); a bug inside a row function or `observed_reads` is for the
    row functions' own tests."""
    t, q, sigma = _record()
    result = tl.run(_inputs(t, q, sigma))
    obs = windows.observed_reads(t, q, sigma=sigma, gap_mult=10.0)
    metric = result.metric
    assert (metric.n_reads, metric.n_segments) == (obs.t_rel_s.size, 2)
    assert len(obs.segments()) == 2
    spacing = tl.spacing_row(obs)
    first_alias, edge = tl.band_edges(obs)
    rows = tl.spectrum_rows(obs)
    floor, ratio = tl.noise_floor_ratio(obs, rows, obs.sigma)
    second = tl.second_difference_row(obs)
    first_half, second_half = tl.g10_halves(obs)
    expected = {
        "median_spacing_s": obs.median_spacing_s,
        "gap_threshold_s": obs.gap_threshold_s,
        "effective_spacing_s": spacing.effective_spacing_s,
        "mean_spacing_s": spacing.mean_spacing_s,
        "spacing_cv": spacing.spacing_cv,
        "first_alias_hz": first_alias,
        "band_edge_hz": edge,
        "spectral_window_at_first_alias": tl.window_at_first_alias(obs),
        "spectrum_frequency_hz": rows.frequency_hz,
        "spectrum_per_hz": rows.spectrum_per_hz,
        "edge_level_per_hz": rows.edge_level_per_hz,
        "noise_floor_per_hz": floor,
        "edge_noise_ratio": ratio,
        "second_difference_scale": second.second_difference_scale,
        "value_sd": second.value_sd,
        "second_difference_ratio": second.second_difference_ratio,
        "g10": tl.g10_row(obs, rows).g10,
        "g10_first_half": first_half.g10,
        "g10_second_half": second_half.g10,
        "dependence": tl.dependence_panel(
            obs.values, obs.sigma, n_permutations=19, seed=5
        ),
        "self_consistency": tl.self_consistency(obs.values, obs.sigma),
        "sigma_bound_ratio": tl.sigma_bound_ratio(obs, obs.sigma),
    }
    assert set(expected) == _fields(tl.METRIC_GROUPS, metric.status)
    _assert_fields(metric, expected)
    assert (result.settings["seed"], result.settings["n_permutations"]) == (5, 19)
    assert [row.threshold_label for row in result.thresholds] == [
        label for label, _, _ in THRESHOLDS
    ]
    table = _carve(t, q, THRESHOLDS).windows_in_spec
    for (label, value, big), row in zip(THRESHOLDS, result.thresholds, strict=True):
        assert (row.threshold_label, row.threshold_value, row.big_values_good) == (
            label,
            value,
            big,
        )
        margin_v = windows.margin(obs.values, value, big)
        inside = margin_v >= 0
        exponent = tl.count_exponent_row(obs, inside)
        counts = tuple(
            tl.crossing_counts(obs, margin_v, obs.sigma, k)
            for k in tl.CONFIDENCE_LEVELS
        )
        k1 = counts[tl.CONFIDENCE_LEVELS.index(1)]
        runs = tl.run_length_row(table[table["threshold_label"] == label], inside)
        bounds = tl.in_spec_fraction_bounds(obs, margin_v, obs.sigma, 1.0)
        expected = {
            "decimation_ladder": exponent.decimation_ladder,
            "down_crossings": exponent.down_crossings,
            "kept_pairs": exponent.kept_pairs,
            "count_ratio": exponent.count_ratio,
            "count_exponent": exponent.count_exponent,
            "rungs_excluded": exponent.rungs_excluded,
            "crossing_counts": counts,
            "clean_fraction_k1": k1.clean / (k1.clean + k1.merged + k1.unbracketed),
            "read_fraction_in_spec": runs.read_fraction_in_spec,
            "mean_in_spec_run_reads": runs.mean_in_spec_run_reads,
            "run_length_ratio": runs.run_length_ratio,
            "in_spec_fraction_lower": bounds.in_spec_fraction_lower,
            "in_spec_fraction_upper": bounds.in_spec_fraction_upper,
            "unsure_k": bounds.unsure_k,
            "in_spec_autocorr": tl.in_spec_autocorr(obs, inside),
        }
        assert set(expected) == _fields(tl.THRESHOLD_GROUPS, row.status), label
        _assert_fields(row, expected)
        assert runs.read_fraction_in_spec == float(inside.mean())


@pytest.mark.unit
def test_run_reads_at_the_gap_rule_and_confidence_level_it_is_given() -> None:
    """Oracle: the row functions at gap_mult 20 and k = 2, neither of them a default."""
    t, q, sigma = _record()
    result = tl.run(_inputs(t, q, sigma, gap_mult=20.0, unsure_k=2.0))
    obs = windows.observed_reads(t, q, sigma=sigma, gap_mult=20.0)
    assert result.metric.gap_threshold_s == obs.gap_threshold_s
    for (_, value, big), row in zip(THRESHOLDS, result.thresholds, strict=True):
        bounds = tl.in_spec_fraction_bounds(
            obs, windows.margin(obs.values, value, big), obs.sigma, 2.0
        )
        assert row.unsure_k == 2.0
        assert row.in_spec_fraction_lower == bounds.in_spec_fraction_lower
        assert row.in_spec_fraction_upper == bounds.in_spec_fraction_upper


@pytest.mark.unit
@pytest.mark.parametrize(
    "sigma_status",
    ["not computed: sigma propagation not derived", "not declared"],
    ids=["not-derived", "not-declared"],
)
def test_the_sigma_rows_carry_the_sigma_status(sigma_status) -> None:
    """Oracle: specification. Without a usable sigma, every row that reads sigma carries the
    metric's sigma status and no number, at the metric and at every threshold; every other
    row is still computed."""
    t, q, _ = _record()
    result = tl.run(_inputs(t, q, None, sigma_status=sigma_status, clip_floor=-5.0))
    # On a record too short for a spectrum, the noise floor still names the sigma first.
    short = tl.run(
        _inputs(t[:8], q[:8], None, sigma_status=sigma_status, thresholds=())
    )
    assert short.metric.status["spectrum"].startswith("not computed: span shorter")
    assert short.metric.status["noise_floor"] == sigma_status
    for group in tl.METRIC_GROUPS:
        expected = sigma_status if group in tl.SIGMA_GROUPS else tl.COMPUTED
        assert result.metric.status[group] == expected, group
    for row in result.thresholds:
        for group in tl.THRESHOLD_GROUPS:
            expected = sigma_status if group in tl.SIGMA_GROUPS else tl.COMPUTED
            assert row.status[group] == expected, group
        _check_groups(row, tl.THRESHOLD_GROUPS)
    _check_groups(result.metric, tl.METRIC_GROUPS)


@pytest.mark.unit
def test_a_declared_sigma_with_no_finite_value_is_no_sigma_everywhere() -> None:
    """Oracle: specification. Declared but never finite is no sigma: every sigma row, at the
    metric and at every threshold, says so, instead of a crossing count or bounds computed
    as if every read were sure."""
    t, q, _ = _record(n=300)
    result = tl.run(_inputs(t, q, np.full(t.size, np.nan)))
    no_sigma = "not computed: no read with a finite sigma"
    for row, groups in [(result.metric, tl.METRIC_GROUPS)] + [
        (row, tl.THRESHOLD_GROUPS) for row in result.thresholds
    ]:
        for group in groups:
            if group in tl.SIGMA_GROUPS:
                assert row.status[group] == no_sigma, group
        _check_groups(row, groups)
    assert result.metric.sigma_status == no_sigma
    assert result.metric.n_reads_sigma_unknown == t.size


@pytest.mark.unit
def test_the_result_counts_what_was_dropped_before_and_inside_the_step() -> None:
    """Oracle: the load's own counts, repeated; one read with a NaN value dropped by the
    step, one kept with no sigma. Counts the load did not record are None, not 0."""
    t, q, sigma = _record(n=300)
    # Read 10 is dropped for its value, its sigma with it; read 20 is kept with no sigma.
    q[10], sigma[10], sigma[20] = np.nan, np.nan, np.nan
    counts = {"n_rows_dropped_by_schema": 3, "n_rows_dropped_invalid": 0}
    metric = tl.run(_inputs(t, q, sigma, load_counts=counts)).metric
    assert (metric.n_rows_dropped_by_schema, metric.n_rows_dropped_invalid) == (3, 0)
    assert (metric.n_reads_dropped_nonfinite, metric.n_reads_sigma_unknown) == (1, 1)
    bare = tl.run(_inputs(t, q, sigma)).metric
    assert (bare.n_rows_dropped_by_schema, bare.n_rows_dropped_invalid) == (None, None)


@pytest.mark.unit
@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"sigma_status": "not declared"}, "sigma is given"),
        ({"carve_gap_mult": 20.0}, "one gap rule"),
        ({"windows_result": None}, "without the carve's windows"),
        ({"sigma_status": "unknown"}, "unknown sigma status"),
        ({"sigma_status": "not computed: zero variance"}, "unknown sigma status"),
        ({"load_counts": {"n_rows_lost": 1}}, "unknown load counts"),
    ],
    ids=[
        "sigma-given-not-declared",
        "other-gap-rule",
        "no-windows",
        "bad-status",
        "reason-as-status",
        "unknown-count",
    ],
)
def test_contradictory_inputs_raise(kwargs, match) -> None:
    """Oracle: specification. These describe no record that exists."""
    t, q, sigma = _record(n=300)
    with pytest.raises(ValueError, match=match):
        tl.run(_inputs(t, q, sigma, **kwargs))


@pytest.mark.unit
@pytest.mark.parametrize(
    "case, match",
    [
        ("declared-without-sigma", "sigma is absent"),
        ("negative-sigma", "negative"),
        ("zero-k", "unsure_k must be finite and positive"),
        ("repeated-label", "repeated"),
        ("other-record", "one set of reads"),
        ("other-id", "carved on record"),
        ("prefix", "one set of reads"),
        ("no-permutations", "n_permutations must be"),
        ("negative-seed", "seed must be"),
        ("no-finite-sigma-as-input", "set by the step"),
    ],
)
def test_inputs_that_contradict_themselves_raise_whatever_the_sigma_status(
    case, match
) -> None:
    """Oracle: specification. The checks do not wait for a row that happens to read the
    input: a negative sigma, a non-positive k, no permutations, a negative seed or a
    repeated label is refused up front, with sigma not declared where the panel would never
    run; and windows carved on other reads, of another record (by its id, or by its reads)
    or of part of this one, are refused."""
    t, q, sigma = _record(n=300)
    if case == "declared-without-sigma":
        inputs = _inputs(t, q, None)
    elif case == "negative-sigma":
        inputs = _inputs(t, q, -sigma)
    elif case == "zero-k":
        # The carve refuses k = 0 itself, so the windows are carved at 1 and the step is
        # told 0.
        inputs = dataclasses.replace(
            _inputs(t, q, None, sigma_status="not declared"), unsure_k=0.0
        )
    elif case == "repeated-label":
        inputs = dataclasses.replace(
            _inputs(t, q, sigma), thresholds=(("u0", 0.0, True), ("u0", 0.0, True))
        )
    elif case == "other-record":
        t2, q2, _ = _record(seed=22, n=300)
        inputs = _inputs(t, q, sigma, windows_result=_carve(t2, q2, THRESHOLDS))
    elif case == "other-id":
        other = _carve(t, q, THRESHOLDS, dataset_id="elsewhere")
        inputs = _inputs(t, q, sigma, windows_result=other)
    elif case == "prefix":
        # Evenly spaced, so the prefix has the record's median spacing and only the read
        # count tells it apart.
        even = np.arange(300) * 20.0
        inputs = _inputs(
            even, q, sigma, windows_result=_carve(even[:200], q[:200], THRESHOLDS)
        )
    elif case == "no-finite-sigma-as-input":
        inputs = _inputs(t, q, None, sigma_status=tl.NO_FINITE_SIGMA)
    else:
        bad = {"no-permutations": {"n_permutations": 0}, "negative-seed": {"seed": -1}}
        inputs = dataclasses.replace(
            _inputs(t, q, None, sigma_status="not declared"), **bad[case]
        )
    with pytest.raises(ValueError, match=match):
        tl.run(inputs)


@pytest.mark.unit
@pytest.mark.parametrize(
    "thresholds, k, match",
    [
        ([("elsewhere", 0.5, True)], 1.0, "not carved at threshold"),
        ([("u0", 0.3, True), ("u1", 0.8, True)], 1.0, "was carved at"),
        ([("u0", 0.0, False), ("u1", 0.8, True)], 1.0, "was carved at"),
        (list(THRESHOLDS), 2.0, "one confidence level"),
    ],
    ids=["other-label", "same-label-other-value", "same-label-other-side", "other-k"],
)
def test_windows_carved_at_other_thresholds_are_refused(thresholds, k, match) -> None:
    """Oracle: specification. A threshold's rows must come from windows carved at that
    threshold: its label, its value, its side, and the carve's own confidence level."""
    t, q, sigma = _record(n=300)
    other = _carve(t, q, thresholds, k=k)
    with pytest.raises(ValueError, match=match):
        tl.run(_inputs(t, q, sigma, thresholds=THRESHOLDS[:2], windows_result=other))


NAN = float("nan")


@pytest.mark.unit
@pytest.mark.parametrize(
    "change, match",
    [
        (lambda m: {"g10": NAN}, "missing or not finite"),
        (lambda m: {"g10": None}, "missing or not finite"),
        (lambda m: {"clipped_fraction": 0.5}, "carries a value"),
        (lambda m: {"status": {"spacing": tl.COMPUTED}}, "do not name the groups"),
        (
            lambda m: {"status": {**m.status, "clipped": "not computed: banana"}},
            "not a tellability reason",
        ),
        (lambda m: {"sigma_status": "banana"}, "unknown sigma status"),
        (
            lambda m: {"dependence": dataclasses.replace(m.dependence, p_value=NAN)},
            "missing or not finite",
        ),
        (
            lambda m: {
                "dependence": dataclasses.replace(
                    m.dependence,
                    statistics={**m.dependence.statistics, "pearson_abs": NAN},
                )
            },
            "missing or not finite",
        ),
    ],
    ids=[
        "nan-computed",
        "none-computed",
        "value-not-computed",
        "missing-group",
        "unknown-reason",
        "unknown-sigma-status",
        "nan-panel-p",
        "nan-panel-statistic",
    ],
)
def test_a_row_that_breaks_its_statuses_cannot_be_built(change, match) -> None:
    """Oracle: specification, enforced at construction: a row whose numbers or statuses
    contradict the vocabulary is refused when it is made, a NaN inside the dependence panel
    included, not only when a test inspects it."""
    t, q, sigma = _record(n=300)
    metric = tl.run(_inputs(t, q, sigma)).metric
    with pytest.raises(ValueError, match=match):
        dataclasses.replace(metric, **change(metric))


@pytest.mark.unit
def test_a_nan_inside_a_tuple_or_a_record_cannot_be_shipped_as_computed() -> None:
    """Oracle: specification. The finiteness rule reaches inside tuples and the records
    they hold."""
    t, q, sigma = _record(n=300)
    row = tl.run(_inputs(t, q, sigma)).thresholds[0]
    with pytest.raises(ValueError, match="missing or not finite"):
        dataclasses.replace(row, in_spec_autocorr=(NAN,))
    bad_count = dataclasses.replace(row.crossing_counts[0], k=NAN)
    with pytest.raises(ValueError, match="missing or not finite"):
        dataclasses.replace(row, crossing_counts=(bad_count, *row.crossing_counts[1:]))


@pytest.mark.unit
def test_an_error_a_row_did_not_name_is_raised_not_shipped(monkeypatch) -> None:
    """Oracle: specification (`_attempt`'s contract). NotComputed becomes a status; any
    other error is a bug, and `run` raises it rather than ship a plausible status."""
    t, q, sigma = _record(n=300)

    def broken(obs):
        raise ValueError("a bug in a row")

    monkeypatch.setattr(tl, "spacing_row", broken)
    with pytest.raises(ValueError, match="a bug in a row") as excinfo:
        tl.run(_inputs(t, q, sigma))
    assert not isinstance(excinfo.value, tl.NotComputed)


@pytest.mark.unit
def test_a_reason_outside_the_set_is_refused() -> None:
    """Oracle: specification. The vocabulary cannot drift one call site at a time."""
    assert tl.NotComputed("fewer than 3 reads").reason == "fewer than 3 reads"
    for bad in (
        "fewer than many reads",
        "skipped",
        "zero variance!",
        "zero variance\n",
    ):
        with pytest.raises(ValueError, match="not a tellability reason"):
            tl.NotComputed(bad)


@pytest.mark.unit
def test_the_result_loads_back_and_a_stale_one_is_refused() -> None:
    """Oracle: the artifact guard's contract. A result missing a field it now has is stale."""
    t, q, sigma = _record(n=300)
    result = tl.run(_inputs(t, q, sigma))
    loaded = load_artifact(io.BytesIO(pickle.dumps(result)))
    assert loaded.metric.g10 == result.metric.g10
    stale = dataclasses.replace(result)
    del stale.__dict__["settings"]
    with pytest.raises(ValueError, match="stale TellabilityResult"):
        load_artifact(io.BytesIO(pickle.dumps(stale)))
