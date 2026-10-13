"""The pooled overlapping Allan deviation of a gapped record, against three oracles.

- A reference implementation: on one gapless segment the pooled estimator must equal
  `allantools.oadev`, which implements the overlapping Allan deviation independently.
- Hand arithmetic: on a two-segment record the pooled variance is the sum of every
  segment's squared differences of adjacent m-read means over twice their total count,
  computed below with explicit loops; no term spans the gap. A second layout, whose
  whole-record and within-segment median spacings differ, shows tau0 is the former.
- An analytic value: white frequency noise of variance sigma^2 has AVAR(m tau0) =
  sigma^2 / m in expectation, on any segment layout.
"""

from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from quebra.analyzers import allan
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.plots.allan_plot import AllanPlot
from quebra.recipes import RAMSEY_GAP_MULT, configure_ramsey_job


def _config(min_points: int = 1, taus_mode: str = "all") -> dict[str, object]:
    return {
        "allan": {
            "modes": ["raw"],
            "taus_mode": taus_mode,
            "min_points_for_allan": min_points,
        }
    }


def _norm(t: np.ndarray, y: np.ndarray, **extra: object) -> dict[str, object]:
    return {"t_rel_s": t, "delta_hz": y, "meta": {"dataset_id": "synthetic"}, **extra}


# Segment A: five reads 10 s apart. A 1000 s gap. Segment B: three reads 12 s apart. The
# positive spacings are 10 x4, 1000, 12 x2, so tau0 is 10 s, not segment B's own 12 s.
_T = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 1040.0, 1052.0, 1064.0])
_Y = np.array([1.0, 3.0, 2.0, 5.0, 4.0, 0.0, 2.0, 1.0])
_SEGMENTS = [_Y[:5], _Y[5:]]


def _hand_terms(y: np.ndarray, m: int) -> list[float]:
    """Squared differences of adjacent m-read means, one per overlapping position."""
    return [
        (float(np.mean(y[j + m : j + 2 * m])) - float(np.mean(y[j : j + m]))) ** 2
        for j in range(len(y) - 2 * m + 1)
    ]


@pytest.mark.unit
def test_one_segment_matches_the_allantools_reference() -> None:
    """Oracle: `allantools.oadev`, an independent implementation of the same estimator.
    allantools drops a tau holding a single term, so only shared taus are compared."""
    import allantools  # a dev dependency: a missing reference fails, it does not skip

    rng = np.random.default_rng(3)
    y = rng.normal(size=200)
    t = np.arange(200, dtype=float) * 24.0
    result = allan.run(_norm(t, y), _config())
    frame = result.modes["raw"]
    taus, adev, _, ns = allantools.oadev(
        y, rate=1.0 / 24.0, data_type="freq", taus=frame["tau_s"].to_numpy()
    )
    shared = frame.set_index("tau_s").loc[taus]
    np.testing.assert_allclose(shared["adev_hz"].to_numpy(), adev, rtol=1e-12)
    np.testing.assert_array_equal(shared["n_terms"].to_numpy(), ns)


@pytest.mark.unit
def test_two_segments_pool_by_term_count_with_one_tau0() -> None:
    """Oracle: the explicit loops in `_hand_terms`. At m = 1 segment A holds four terms
    and B two; at m = 2 only A holds any. adev_err_hz is the pooled adev over the square
    root of the pooled count."""
    result = allan.run(_norm(_T, _Y), _config())
    frame = result.modes["raw"].set_index("m")
    assert result.meta["tau0_s"] == 10.0
    for m, n_segments in ((1, 2), (2, 1)):
        terms = [x for seg in _SEGMENTS for x in _hand_terms(seg, m)]
        expected = math.sqrt(sum(terms) / (2.0 * len(terms)))
        row = frame.loc[m]
        assert row["tau_s"] == m * 10.0
        assert row["adev_hz"] == pytest.approx(expected, rel=1e-12)
        assert row["n_terms"] == len(terms)
        assert row["n_segments"] == n_segments
        assert row["adev_err_hz"] == pytest.approx(
            expected / math.sqrt(len(terms)), rel=1e-12
        )
    assert list(frame.index) == [1, 2]


@pytest.mark.statistical
def test_white_frequency_noise_gives_sigma_over_root_m_on_a_gapped_record() -> None:
    """Oracle: AVAR = sigma^2 / m for white frequency noise, in expectation.

    Twenty segments of 300 to 700 reads, with gaps between them. A statistical control on
    one seeded draw, with a 10% tolerance on the deviation that was not set from a measured
    spread; the hand-arithmetic test above is the evidence for the pooling.
    """
    rng = np.random.default_rng(20261013)
    sigma = 2.5
    times, values, start = [], [], 0.0
    for length in rng.integers(300, 701, size=20):
        times.append(start + np.arange(length) * 24.0)
        values.append(rng.normal(scale=sigma, size=length))
        start = times[-1][-1] + 5000.0
    result = allan.run(
        _norm(np.concatenate(times), np.concatenate(values)),
        _config(taus_mode="octave"),
    )
    frame = result.modes["raw"].set_index("m")
    assert result.diagnostics["n_segments"] == 20
    for m in (1, 2, 4, 8, 16):
        ratio = frame.loc[m, "adev_hz"] * math.sqrt(m) / sigma
        assert abs(ratio - 1.0) < 0.10, f"m={m}: adev*sqrt(m)/sigma = {ratio:.4f}"


@pytest.mark.unit
def test_the_fractional_series_divides_by_the_declared_nu0_only() -> None:
    """Oracle: specification. ADEV / nu0 with the kwarg's value, which the artifact
    carries; the record's own frequency column, about 2 MHz here, is never read."""
    rng = np.random.default_rng(5)
    t = np.arange(100, dtype=float) * 24.0
    y = rng.normal(size=100)
    norm = _norm(t, y, raw_frequency_hz=2.0e6 + y)
    result = allan.run(norm, _config(), fractional=True, nu0_hz=5.5e9)
    assert result.nu0_hz == 5.5e9
    np.testing.assert_allclose(
        result.fractional_adev,
        result.modes["raw"]["adev_hz"].to_numpy() / 5.5e9,
        rtol=1e-12,
    )


@pytest.mark.unit
@pytest.mark.parametrize("nu0_hz", [0.0, -1.0, float("nan")])
def test_a_non_positive_or_non_finite_nu0_raises(nu0_hz: float) -> None:
    """Oracle: specification. Such a nu0 has no fractional meaning."""
    t = np.arange(20, dtype=float)
    with pytest.raises(ValueError, match="nu0_hz must be finite and positive"):
        allan.run(_norm(t, np.ones(20)), _config(), fractional=True, nu0_hz=nu0_hz)


def _ramsey_job(**allan_kwargs: object) -> Job:
    job = Job("synthetic_allan_wiring")
    configure_ramsey_job(
        job,
        Dataset(path="data/none.pickle", qubit=1, device="6D2S"),
        profile="overnight",
        include_fidelity=False,
        include_tlf=True,
        **allan_kwargs,  # type: ignore[arg-type]
    )
    return job


@pytest.mark.unit
@pytest.mark.parametrize(
    "call",
    [
        lambda: allan.run(
            _norm(np.arange(20.0), np.ones(20)), _config(), fractional=True
        ),
        lambda: allan.run(_norm(np.arange(20.0), np.ones(20)), _config(), nu0_hz=5.5e9),
        lambda: _ramsey_job(allan_fractional=True),
        lambda: _ramsey_job(allan_nu0_hz=5.5e9),
        lambda: _ramsey_job(allan_fractional=True, allan_nu0_hz=-5.5e9),
    ],
    ids=["run-no-nu0", "run-nu0-unread", "job-no-nu0", "job-nu0-unread", "job-bad-nu0"],
)
def test_fractional_and_nu0_go_together_at_run_and_at_graph_build(call) -> None:
    """Oracle: specification. A fractional series without a declared nu0, or a nu0 that
    nothing reads, is refused; the job refuses it before any node runs."""
    with pytest.raises(ValueError, match="nu0"):
        call()


@pytest.mark.unit
def test_a_repeated_timestamp_inside_a_segment_raises() -> None:
    """Oracle: specification. The equal-spacing estimator has no spacing for two reads at
    one time."""
    t = np.array([0.0, 10.0, 10.0, 20.0, 30.0, 40.0])
    with pytest.raises(ValueError, match="repeated timestamp"):
        allan.run(_norm(t, np.arange(6.0)), _config())


@pytest.mark.unit
def test_a_tau_is_kept_by_its_pooled_term_count() -> None:
    """Oracle: the hand counts above. At m = 1 the two segments pool 4 + 2 = 6 terms and at
    m = 2 only segment A's 2, so a minimum of 5 keeps m = 1, though A alone holds 4, and
    drops m = 2; a minimum of 7 keeps nothing, which raises."""
    frame = allan.run(_norm(_T, _Y), _config(min_points=5)).modes["raw"]
    assert frame["m"].tolist() == [1]
    assert frame["n_terms"].tolist() == [6]
    with pytest.raises(ValueError, match="no tau keeps 7 pooled terms"):
        allan.run(_norm(_T, _Y), _config(min_points=7))


@pytest.mark.unit
def test_gap_mult_reaches_the_segmentation() -> None:
    """Oracle: the explicit loops on the whole record. At gap_mult 200 the 1000 s step is
    no gap, so the record is one segment and its m = 1 terms include the step from 4 to 0."""
    result = allan.run(_norm(_T, _Y), _config(), gap_mult=200.0)
    terms = _hand_terms(_Y, 1)
    assert result.diagnostics["n_segments"] == 1
    assert result.modes["raw"].set_index("m").loc[1, "adev_hz"] == pytest.approx(
        math.sqrt(sum(terms) / (2.0 * len(terms))), rel=1e-12
    )


@pytest.mark.unit
def test_octave_mode_reaches_the_largest_power_of_two() -> None:
    """Oracle: arithmetic. Sixteen reads allow m up to 8, so octave mode gives 1, 2, 4, 8."""
    t = np.arange(16, dtype=float)
    frame = allan.run(_norm(t, np.sin(t)), _config(taus_mode="octave")).modes["raw"]
    assert frame["m"].tolist() == [1, 2, 4, 8]


@pytest.mark.unit
def test_tau0_and_the_spacing_diagnostics_on_a_layout_with_a_lone_read() -> None:
    """Oracle: arithmetic. Spacings 10, 10, 1000, 12, 1000 cut the record into segments of 3,
    2 and 1 reads. The whole record's median positive spacing is 12 s while the
    within-segment median is 10 s, and tau0 is the former. Only the lone read is too short
    for a difference; the two-read segment holds one. Within segments the spacings are 10,
    10, 12: CV = sqrt(2) / 16 with the 1/n variance."""
    t = np.array([0.0, 10.0, 20.0, 1020.0, 1032.0, 2032.0])
    result = allan.run(_norm(t, np.arange(6.0)), _config())
    assert result.meta["tau0_s"] == 12.0
    assert result.diagnostics["n_segments"] == 3
    assert result.diagnostics["n_segments_too_short"] == 1
    assert result.diagnostics["spacing_cv_within_segments"] == pytest.approx(
        math.sqrt(2.0) / 16.0, rel=1e-12
    )


@pytest.mark.unit
@pytest.mark.parametrize("nu0_hz", [True, "5.5e9"], ids=["bool", "string"])
def test_a_nu0_that_is_not_a_number_is_refused(nu0_hz: object) -> None:
    """Oracle: specification. `True` would read as 1 Hz and a string as a number."""
    with pytest.raises(TypeError, match="number of hertz"):
        allan.check_nu0_hz(nu0_hz)


@pytest.mark.unit
def test_allan_and_the_telegraph_read_one_ramsey_gap_rule() -> None:
    """Oracle: specification. One gap rule for every node that segments the filtered reads,
    passed as a step kwarg from one constant."""
    dag = _ramsey_job().dag
    gap_mults = {
        node_id: node.kwargs["gap_mult"]
        for node_id, node in dag.items()
        if node_id in {"allan", "tlf"}
    }
    assert gap_mults == {"allan": RAMSEY_GAP_MULT, "tlf": RAMSEY_GAP_MULT}


@pytest.mark.unit
def test_the_figure_prints_the_spacing_and_what_it_dropped() -> None:
    """Oracle: the layout of the tau0 test, as the caption prints it: three segments, one
    too short, spacing CV sqrt(2) / 16 = 0.088, nothing dropped."""
    t = np.array([0.0, 10.0, 20.0, 1020.0, 1032.0, 2032.0])
    result = allan.run(_norm(t, np.arange(6.0) ** 2), _config())
    figure = AllanPlot("unit").build_matplotlib(result)
    try:
        text = "\n".join(x.get_text() for ax in figure.axes for x in ax.texts)
    finally:
        plt.close(figure)
    assert "3 segments, 1 too short; spacing CV within segments 0.09" in text
    assert (
        "dropped: 0 non-finite reads, 0 non-positive or non-finite ADEV points" in text
    )
