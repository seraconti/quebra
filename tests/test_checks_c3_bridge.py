"""C3's R bridge must degrade to `None`, never to a wrong number and never to a crash.

R is an OPTIONAL out-of-process dependency, so the decline paths are the ones that have to
hold everywhere: a missing interpreter, a multi-process input and a series too short for
the lag all return a result whose `p_value` is None and whose `notes` say why, and nothing
raises. Three tests monkeypatch the probe, so those assert the same thing on a machine with
R and on one without. The other three are machine-independent for their own reasons: two
never reach the probe (the data guard at the top of `run` fires first, and one checks only
that the R file ships beside its wrapper), and
`test_records_whether_rscript_is_absent_here` reads the real machine and skips where R is
present. The numeric path past `_invoke_rscript` is NOT pinned here.
"""

from __future__ import annotations

import numpy as np
import pytest

import quebra.analyzers.checks.c3_serial_copula as c3
from quebra.analyzers.checks.result import Segment

pytestmark = pytest.mark.unit


def _segment(n: int = 30, seed: int = 3) -> Segment:
    rng = np.random.default_rng(seed)
    gaps = rng.exponential(size=n + 40)
    tau = float(n)
    kept = gaps[: int(np.searchsorted(np.cumsum(gaps), tau, side="left"))]
    return Segment(x=kept, tau=tau)


def test_missing_r_returns_none_without_raising(monkeypatch):
    monkeypatch.setattr(c3, "rscript_path", lambda: None)
    result = c3.run([_segment()])
    assert result.p_value is None
    assert "R unavailable" in result.notes
    assert result.n_events > 0


def test_multiple_segments_are_declined_with_a_reason(monkeypatch):
    monkeypatch.setattr(c3, "rscript_path", lambda: "/usr/bin/Rscript")
    result = c3.run([_segment(20, 1), _segment(20, 2)])
    assert result.p_value is None
    assert "multi-process" in result.notes


def test_a_series_too_short_for_the_lag_is_declined(monkeypatch):
    monkeypatch.setattr(c3, "rscript_path", lambda: "/usr/bin/Rscript")
    rng = np.random.default_rng(0)
    short = Segment(x=rng.exponential(size=3), tau=99.0)
    result = c3.run([short], max_lag=5)
    assert result.p_value is None
    assert "too short" in result.notes


def test_invalid_segments_still_raise():
    """Graceful R handling must not soften the ordinary data guards."""
    with pytest.raises(ValueError, match="strictly positive"):
        c3.run([Segment(x=np.array([1.0, 0.0, 2.0]), tau=99.0)])


def test_records_whether_rscript_is_absent_here():
    """Records whether `Rscript` is absent on the machine this suite runs on.

    Skips rather than fails where R exists: a machine with R is a healthier machine, and a
    red test there would punish the fix. The skip message is the actionable part.
    """
    if c3.rscript_path() is not None:
        pytest.skip(
            "Rscript is on PATH here, so C3's numeric path is exercisable and this "
            "test has nothing to assert. The promotion report's five-check scope is "
            "unaffected: C3 has no bench cell, so its size and power stay unmeasured."
        )
    assert c3.rscript_path() is None


def test_the_r_script_ships_next_to_its_wrapper():
    assert c3.R_SCRIPT.exists(), f"missing {c3.R_SCRIPT}"
    body = c3.R_SCRIPT.read_text()
    assert "serialIndepTest" in body and "commandArgs" in body


def test_a_probe_that_answers_nothing_is_distinguishable_from_one_that_never_ran(
    monkeypatch, tmp_path
):
    """Oracle: the CheckLedger field defaults these two probes feed.

    Both probes report into fields whose `not asked` and `absent` states already carry an
    empty value, so a probe that RAN and said nothing must not return that same empty value.
    Neither sentinel was guarded when it was added: deleting either left the whole suite
    green, which is why this test names both in one place.
    """
    stub = tmp_path / "Rscript"
    stub.write_text("#!/bin/sh\nexit 0\n")
    stub.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    assert c3.rscript_path() is not None, "the stub interpreter must be visible"

    assert c3.r_version() == "probe returned nothing"

    # The two probes answer differently ON PURPOSE. `r_library_paths` returns `()`, because
    # its other consumer joins the tuple into `R_LIBS` for the real bridge call and a
    # diagnostic string there would become a filesystem search path. Distinguishing an empty
    # answer from a probe that never ran is the ledger's job, so the sentinel lives there.
    assert c3.r_library_paths() == ()

    from quebra.analyzers.check_ledger import _r_provenance

    assert _r_provenance(True)[2] == ("probe returned nothing",)
