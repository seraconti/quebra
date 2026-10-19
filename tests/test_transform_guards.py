"""The filter's and TLF's guards, each shown to fail when its guard is removed.

What they assert is one property: a step hands downstream either a value the instrument
produced or nothing, never a value invented to stand in for one, and never a column that has
lost its alignment with `t_rel_s`.

Scope is deliberately the GUARDS, not the numbers. The chi-square threshold is a citation, not
a tuned value; whether a failed fit can silently empty a record, a read can be dropped without
being counted, or a column can come out misaligned, is answerable here.

Every assertion here was checked by mutating the guard it covers and confirming the test goes
red. A guard whose test passes either way is a claim, not a check.
"""

from __future__ import annotations

import numpy as np
import pytest

from quebra.analyzers.tlf import run as run_tlf
from quebra.transforms.filter import run as run_filter

pytestmark = pytest.mark.unit


def _norm(**overrides) -> dict[str, object]:
    """A minimal well-formed Norm: 10 reads, a clean linear detuning."""
    n = 10
    norm: dict[str, object] = {
        "t_rel_s": np.arange(n, dtype=float) * 10.0,
        "delta_hz": np.linspace(100.0, 200.0, n),
        "meta": {"dataset_id": "unit"},
    }
    norm.update(overrides)
    return norm


_FILTER_CFG = {
    "filter": {
        "apply_chi_squared": False,
    },
    "dataset_profile": "unit",
}


# --- filter: the chi-square guard refuses a record that is not the published one ----


@pytest.mark.parametrize(
    "bad", [0.8, 0.85, float("nan")], ids=["at_cut", "above", "nan"]
)
def test_the_chi_square_guard_refuses_a_read_at_or_above_the_published_cut(bad):
    """Oracle: Capannelli et al., Supplementary notes, section III.C, p. 9 (drift datasets
    scoring a normalised chi-squared of 0.8 or more are filtered out; the released records hold
    what is left, Supplementary Figs. 10 and 11, pp. 10-11).

    A record inside the cut passes untouched; the same record with one read at the cut,
    above it, or without a chi-square is refused rather than silently trimmed.
    """
    cfg = {
        "filter": {
            "apply_chi_squared": True,
            "chi_squared_threshold": 0.8,
        }
    }
    clean = np.linspace(0.1, 0.79, 10)
    kept = run_filter(_norm(chi_squared=clean), cfg, None).final_norm
    assert len(np.asarray(kept["t_rel_s"])) == 10, (
        "a record inside the cut loses no read"
    )
    dirty = clean.copy()
    dirty[3] = bad
    with pytest.raises(ValueError, match="normalised chi-square >= 0.8"):
        run_filter(_norm(chi_squared=dirty), cfg, None)


@pytest.mark.parametrize("key", ["apply_sigma", "apply_frequency_window", "typo"])
def test_the_filter_refuses_a_config_key_it_does_not_read(key):
    """Oracle: specification. A cut that no longer exists is refused, not silently skipped."""
    cfg = {"filter": {"apply_chi_squared": False, key: True}}
    with pytest.raises(KeyError, match="unknown filter config key"):
        run_filter(_norm(), cfg, None)


# --- filter: every column keeps its alignment with t_rel_s ------------------


@pytest.mark.parametrize("column", ["T2star_s", "chi_squared"])
def test_filter_raises_on_a_column_it_cannot_keep_aligned(column):
    """A column left at full length while its siblings shrink is index-misaligned for every
    step downstream, and no consumer raises on it. Checked on a value column and on the
    chi-square column, with the guard off."""
    with pytest.raises(ValueError, match="filter cannot mask"):
        run_filter(_norm(**{column: np.ones(4)}), _FILTER_CFG, None)


def test_filter_keeps_every_row_aligned_column_the_same_length():
    """The negative control: a Norm that loses a row keeps every column in step.

    The dropped row is a failed fit (a non-finite detuning), one of the two kinds of read the
    filter removes (a non-finite detuning or a non-finite time).
    """
    values = np.linspace(100.0, 200.0, 10)
    values[4] = np.nan
    result = run_filter(
        _norm(delta_hz=values, T2star_s=np.linspace(1e-6, 2e-6, 10)),
        _FILTER_CFG,
        None,
    )
    lengths = {
        key: len(np.asarray(value))
        for key, value in result.final_norm.items()
        if key != "meta" and np.asarray(value).ndim
    }
    assert len(set(lengths.values())) == 1, lengths
    assert lengths["t_rel_s"] == 9, "the non-finite read should have been dropped"


def test_a_read_without_a_finite_time_is_dropped_and_counted():
    """Oracle: specification. A read with a finite detuning but no finite time cannot be
    placed: it is dropped, and the drop is counted in the stage table rather than lost."""
    times = np.arange(10, dtype=float) * 10.0
    times[4] = np.nan
    result = run_filter(_norm(t_rel_s=times), _FILTER_CFG, None)
    assert len(np.asarray(result.final_norm["t_rel_s"])) == 9
    assert result.counts[-1] == {"stage": "non_finite", "n_points": 9, "n_dropped": 1}


def test_filter_rejects_a_ragged_column():
    """A column that is not one value per read cannot be masked with the reads, so the
    filter raises rather than passing it through misaligned."""
    ragged = [[1, 2], [3]]
    with pytest.raises(ValueError, match="filter cannot mask"):
        run_filter(_norm(odd=ragged), _FILTER_CFG, None)


def test_an_all_non_finite_series_raises_naming_its_own_cause():
    with pytest.raises(ValueError, match="no read with both a finite 'delta_hz'"):
        run_filter(_norm(delta_hz=np.full(10, np.nan)), _FILTER_CFG, None)


# --- tlf: reproducible under a fixed identity ------------------------------


def _bimodal(n: int = 400) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(11)
    values = np.concatenate(
        [rng.normal(1.0e6, 300.0, n), rng.normal(1.002e6, 300.0, n)]
    )
    return values, np.arange(values.size, dtype=float) * 24.0


def test_the_seed_reaches_both_estimators():
    """Asserted on `random_state`, not on two runs agreeing.

    Agreement is the wrong oracle: `GaussianMixture` initialises by k-means, and on
    well-separated lobes it converges to the same solution from any start - so two unseeded
    runs agree too, and a test of agreement passes with `random_state=None`. Checked by
    removing the seed: the agreement form stayed green. What has to hold is that the seed
    reaches sklearn, which is observable directly.
    """
    values, timestamps = _bimodal()
    result = run_tlf(values, timestamps, seed=4242)
    assert result.gmm1.random_state == 4242
    assert result.gmm2.random_state == 4242, (
        "the 2-component fit is the one whose init decides is_bimodal"
    )


def test_the_same_seed_gives_the_same_verdict():
    """Weaker than the check above and kept as a sanity floor: whatever the seed does, it
    must not make one input give two answers."""
    values, timestamps = _bimodal()
    first = run_tlf(values, timestamps, seed=20260902)
    second = run_tlf(values, timestamps, seed=20260902)
    assert first.bic_delta == second.bic_delta
    assert first.is_bimodal == second.is_bimodal
    assert first.mean_dwell_state0_s == second.mean_dwell_state0_s


def test_the_seed_is_required():
    """Not defaulted, so a caller cannot omit it and get an irreproducible answer."""
    values, timestamps = _bimodal()
    with pytest.raises(TypeError):
        run_tlf(values, timestamps)


def _two_component_fit_raises(monkeypatch, error: Exception) -> None:
    import quebra.analyzers.tlf as tlf_module

    real = tlf_module.GaussianMixture

    class _FailsOnTwo:
        def __init__(self, *args, n_components: int = 1, **kwargs):
            self._n = n_components
            self._inner = real(*args, n_components=n_components, **kwargs)

        def fit(self, *args, **kwargs):
            if self._n == 2:
                raise error
            return self._inner.fit(*args, **kwargs)

        def __getattr__(self, name):
            return getattr(self._inner, name)

    monkeypatch.setattr(tlf_module, "GaussianMixture", _FailsOnTwo)


def test_a_failed_fit_reports_fit_failed_rather_than_a_verdict(monkeypatch):
    """`bic_delta = 0.0` with `gmm2` aliased to `gmm1` let a consumer read "no evidence of
    bimodality" off a comparison that never happened."""
    _two_component_fit_raises(monkeypatch, ValueError("singular covariance"))
    values, timestamps = _bimodal()
    result = run_tlf(values, timestamps, seed=1)

    assert result.fit_failed is True
    assert result.diagnostics["fit_error"] == "ValueError: singular covariance"
    assert result.bic_delta is None, (
        "a comparison that did not happen reported a number"
    )
    assert result.gmm2 is None, "gmm2 aliased to the 1-component fit"
    assert result.is_bimodal is False


def test_a_programming_error_in_the_fit_is_not_read_as_a_failed_fit(monkeypatch):
    """Only the fit's own failure types become `fit_failed`; anything else propagates."""
    _two_component_fit_raises(
        monkeypatch, TypeError("a bug, not a singular covariance")
    )
    values, timestamps = _bimodal()
    with pytest.raises(TypeError, match="a bug"):
        run_tlf(values, timestamps, seed=1)
