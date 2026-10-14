"""The tellability step as a whole: its status vocabulary and its wiring.

Oracles:
- specification: every status is `computed`, `not declared` or `not computed: <reason>` with
  a reason from the fixed set; a group's fields are all set when it is computed and all None
  when it is not, and no computed field is NaN;
- the row functions themselves: `run` must report what each row's function gives on the same
  observed reads;
- specification for sigma: a metric whose sigma is not declared, or known but not usable,
  carries that status on every row that reads sigma.
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

THRESHOLDS = (("u0", 0.0, True), ("u1", 0.8, True))


def _record(seed: int = 21, n: int = 1500):
    """A drifting signal read at irregular times, with one gap; sigma known and varying."""
    rng = np.random.default_rng(seed)
    steps = rng.uniform(10.0, 40.0, size=n - 1)
    steps[n // 2] += 20000.0
    t = np.r_[0.0, np.cumsum(steps)]
    sigma = rng.uniform(0.2, 0.6, size=n)
    q = np.sin(2.0 * np.pi * t / 4000.0) + rng.normal(scale=sigma)
    return t, q, sigma


def _inputs(
    t,
    q,
    sigma,
    *,
    sigma_status="declared",
    thresholds=THRESHOLDS,
    gap_mult=10.0,
    carve_gap_mult=None,
    clip_floor=None,
    windows_result="carve",
):
    if windows_result == "carve":
        windows_result = (
            windows.run(
                windows.WindowsInputs(
                    t_rel_s=t,
                    values=q,
                    thresholds=list(thresholds),
                    gap_mult=gap_mult if carve_gap_mult is None else carve_gap_mult,
                )
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
        unsure_k=1.0,
        gap_mult=gap_mult,
        clip_floor=clip_floor,
        n_permutations=19,
        seed=5,
    )


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
            for value in values:
                assert value is not None, (group, fields)
                if isinstance(value, float):
                    assert math.isfinite(value), (group, value)
                if isinstance(value, np.ndarray):
                    assert np.all(np.isfinite(value)), group
        else:
            assert all(value is None for value in values), (group, status)


@pytest.mark.unit
@pytest.mark.parametrize("record", ["full", "degenerate"])
def test_every_status_comes_from_the_vocabulary_and_matches_its_fields(record) -> None:
    """Oracle: specification. On a full record every group is computed; on a constant,
    short one most are not, each with its reason, and none ships NaN as computed."""
    if record == "full":
        t, q, sigma = _record()
        result = tl.run(_inputs(t, q, sigma))
        assert set(result.metric.status.values()) - {tl.COMPUTED} == {
            "not computed: no floor for this metric"
        }
    else:
        t = np.arange(8.0) * 10.0
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


@pytest.mark.unit
def test_run_reports_what_the_row_functions_give() -> None:
    """Oracle: the row functions, called directly on the same observed reads."""
    t, q, sigma = _record()
    result = tl.run(_inputs(t, q, sigma))
    obs = windows.observed_reads(t, q, sigma=sigma, gap_mult=10.0)
    rows = tl.spectrum_rows(obs)
    assert result.metric.g10 == tl.g10_row(obs, rows).g10
    assert result.metric.second_difference_ratio == (
        tl.second_difference_row(obs).second_difference_ratio
    )
    assert result.metric.self_consistency == tl.self_consistency(obs.values, obs.sigma)
    for (label, value, big), row in zip(THRESHOLDS, result.thresholds):
        margin_v = windows.margin(obs.values, value, big)
        assert (
            row.count_exponent
            == tl.count_exponent_row(obs, margin_v >= 0).count_exponent
        )
        bounds = tl.in_spec_fraction_bounds(obs, margin_v, obs.sigma, 1.0)
        assert row.in_spec_fraction_lower == bounds.in_spec_fraction_lower
        assert row.threshold_label == label


@pytest.mark.unit
@pytest.mark.parametrize(
    "sigma_status",
    ["not computed: sigma propagation not derived", "not declared"],
    ids=["not-derived", "not-declared"],
)
def test_the_sigma_rows_carry_the_sigma_status(sigma_status) -> None:
    """Oracle: specification. Without a usable sigma, every row that reads sigma carries the
    metric's sigma status and no number; every other row is still computed."""
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
    _check_groups(result.metric, tl.METRIC_GROUPS)


@pytest.mark.unit
@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"sigma_status": "not declared"}, "sigma is given"),
        ({"carve_gap_mult": 20.0}, "one gap rule"),
        ({"windows_result": None}, "without the carve's windows"),
        ({"sigma_status": "unknown"}, "unknown sigma status"),
    ],
    ids=["sigma-given-not-declared", "other-gap-rule", "no-windows", "bad-status"],
)
def test_contradictory_inputs_raise(kwargs, match) -> None:
    """Oracle: specification. These describe no record that exists."""
    t, q, sigma = _record(n=300)
    with pytest.raises(ValueError, match=match):
        tl.run(_inputs(t, q, sigma, **kwargs))


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
    other = windows.run(
        windows.WindowsInputs(t_rel_s=t, values=q, thresholds=thresholds, k=k)
    )
    with pytest.raises(ValueError, match=match):
        tl.run(_inputs(t, q, sigma, windows_result=other))


@pytest.mark.unit
def test_a_reason_outside_the_set_is_refused() -> None:
    """Oracle: specification. The vocabulary cannot drift one call site at a time."""
    assert tl.NotComputed("fewer than 3 reads").reason == "fewer than 3 reads"
    for bad in ("fewer than many reads", "skipped", "zero variance!"):
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
