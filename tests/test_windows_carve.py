"""Carving semantics, ported from monoliths/v13fig/test_carve.py.

The monolith file could not be used as-is: it has no test_ functions, its `check()`
helper prints instead of raising, and a module-scope `raise SystemExit` makes pytest
abort the whole session with an INTERNALERROR. Its 19 checks are reproduced here as
asserts against analyzers/windows.py.

ONE CHECK DID NOT PORT: the monolith's `min_reads=3` filter test. `min_reads` is
deliberately not implemented (see the analyzers/windows.py docstring) - this repo's carve
has never had that filter, so omitting it is what preserves parity. The check is gone
rather than faked.

The reference fixture: reads at t = 0..9, then a 100-unit gap, then 10 more reads.
  idx 0,1      in   -> window opens at index 0        -> birth scan_start
  idx 2,3      out
  idx 4,5,6    in   -> real up-crossing and down      -> COMPLETE
  idx 7        out
  idx 8,9      in   -> open when the gap hits         -> death gap_start
  idx 10,11,12 in   -> first reads after the gap      -> birth gap_resume
  idx 13       out
  idx 14..19   in   -> runs to the end                -> death scan_end
"""

from __future__ import annotations

import math
import re

import numpy as np
import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from quebra.analyzers import windows
from quebra.core.paths import repo_root


T = np.r_[np.arange(10) * 1.0, np.arange(10) * 1.0 + 110.0]
Y = np.array(
    [2, 2, 0, 0, 2, 2, 2, 0, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2],
    dtype=float,
)


def _carve(t: np.ndarray, y: np.ndarray, u: float, gap_mult: float = 10.0):
    _, gap_threshold_s, _, _ = windows.spacing(t, gap_mult)
    return windows.carve(
        t,
        windows.in_spec_mask(y, u, True),
        windows.gap_flags(t, gap_threshold_s),
    )


def _bounds(wins: list[dict[str, object]]) -> list[tuple[int, int]]:
    return [(int(w["s"]), int(w["e"])) for w in wins]


def _complete(wins: list[dict[str, object]]) -> list[dict[str, object]]:
    """Windows with an observed birth AND an observed death."""
    return [
        w
        for w in wins
        if w["birth_type"] == windows.BIRTH_UP_CROSSING
        and w["death_type"] == windows.DEATH_DOWN_CROSSING
    ]


@pytest.mark.unit
def test_spacing() -> None:
    median_s, gap_threshold_s, n_gaps, n_nonpositive = windows.spacing(T, gap_mult=10.0)
    assert median_s == 1.0
    assert gap_threshold_s == 10.0
    assert n_gaps == 1
    assert n_nonpositive == 0


@pytest.mark.statistical
def test_carving_bounds_births_and_deaths() -> None:
    wins = _carve(T, Y, u=1.0)
    assert len(wins) == 5
    assert _bounds(wins) == [(0, 2), (4, 7), (8, 10), (10, 13), (14, 20)]
    assert [w["birth_type"] for w in wins] == [
        windows.BIRTH_SCAN_START,
        windows.BIRTH_UP_CROSSING,
        windows.BIRTH_UP_CROSSING,
        windows.BIRTH_GAP_RESUME,
        windows.BIRTH_UP_CROSSING,
    ]
    assert [w["death_type"] for w in wins] == [
        windows.DEATH_DOWN_CROSSING,
        windows.DEATH_DOWN_CROSSING,
        windows.DEATH_GAP_START,
        windows.DEATH_DOWN_CROSSING,
        windows.DEATH_SCAN_END,
    ]


@pytest.mark.statistical
def test_complete_excludes_every_unobserved_endpoint() -> None:
    wins = _carve(T, Y, u=1.0)
    complete = _complete(wins)
    assert _bounds(complete) == [(4, 7)]
    assert all(w["birth_type"] != windows.BIRTH_SCAN_START for w in complete)
    assert all(w["birth_type"] != windows.BIRTH_GAP_RESUME for w in complete)
    assert all(w["death_type"] != windows.DEATH_GAP_START for w in complete)
    assert all(w["death_type"] != windows.DEATH_SCAN_END for w in complete)


@pytest.mark.unit
def test_gap_mult_sensitivity() -> None:
    _, gap_threshold_s, n_gaps, _ = windows.spacing(T, gap_mult=200.0)
    assert (gap_threshold_s, n_gaps) == (200.0, 0)
    # Without a gap, idx 8..12 merge into a single window.
    wins = _carve(T, Y, u=1.0, gap_mult=200.0)
    assert _bounds(wins) == [(0, 2), (4, 7), (8, 13), (14, 20)]


@pytest.mark.unit
def test_duplicate_timestamps_do_not_collapse_the_threshold() -> None:
    t_dup = np.repeat(np.arange(10) * 2.0, 2)
    median_s, gap_threshold_s, _, n_nonpositive = windows.spacing(t_dup, gap_mult=10.0)
    assert median_s == 2.0  # median over POSITIVE steps only
    assert n_nonpositive == 10
    assert gap_threshold_s == 20.0


@pytest.mark.unit
def test_gap_beats_a_simultaneous_down_crossing() -> None:
    # The read after the gap is out-of-spec: the open window must die gap_start, not
    # down_crossing. Gap wins ties because the gap check runs first.
    t = np.r_[np.arange(4) * 1.0, np.array([100.0])]
    y = np.array([2, 2, 2, 2, 0], dtype=float)
    wins = _carve(t, y, u=1.0)
    assert len(wins) == 1
    assert wins[0]["death_type"] == windows.DEATH_GAP_START


@pytest.mark.unit
def test_a_reading_exactly_at_the_threshold_is_in_spec_in_both_directions() -> None:
    """Oracle: `AGENTS.md` section 5, in spec means `margin >= 0`, so equality is in spec.

    Every carve fixture elsewhere sits well clear of its threshold, so without this test a
    `>=` mutated to `>` in `windows.in_spec_mask` leaves the suite green. Both directions
    are pinned, and the neighbours too, so a shifted comparison fails as well as a flipped
    one.
    """
    at = np.array([5.0])
    assert windows.in_spec_mask(at, 5.0, big_values_good=True).tolist() == [True]
    assert windows.in_spec_mask(at, 5.0, big_values_good=False).tolist() == [True]
    # and the neighbours, so the test fails on a shifted comparison rather than only a
    # flipped one
    assert windows.in_spec_mask(np.array([4.99, 5.01]), 5.0, True).tolist() == [
        False,
        True,
    ]
    assert windows.in_spec_mask(np.array([4.99, 5.01]), 5.0, False).tolist() == [
        True,
        False,
    ]


def _run_short(thresholds):
    return windows.run(
        windows.WindowsInputs(
            t_rel_s=np.arange(4.0), values=np.ones(4), thresholds=thresholds
        )
    )


@pytest.mark.unit
@pytest.mark.parametrize(
    ("call", "match"),
    [
        (lambda: windows.margin(np.ones(2), np.nan, True), "threshold_value"),
        (lambda: windows.in_spec_mask(np.r_[np.nan, 1.0], 0.5, True), "finite values"),
        (
            lambda: windows.run(
                windows.WindowsInputs(
                    t_rel_s=np.zeros(1),
                    values=np.ones(1),
                    thresholds=[("a", np.inf, True)],
                )
            ),
            "non-finite value",
        ),
        (lambda: _run_short([("a", 0.5, True), ("a", 2.0, True)]), "unique"),
    ],
    ids=["nan_threshold", "nan_value", "inf_threshold_one_read", "repeated_label"],
)
def test_a_read_with_no_side_or_a_repeated_label_raises(call, match) -> None:
    """Oracle: AGENTS.md section 3, errors are raised, not swallowed.

    A NaN compared with zero reads as out of spec, so a NaN threshold or value would
    classify silently. A one-read record never reaches `margin`, so `run` checks the
    ladder itself. Two thresholds sharing a label would merge in every keyed table.
    """
    with pytest.raises(ValueError, match=match):
        call()


@pytest.mark.unit
def test_equality_is_not_a_gap() -> None:
    # spacing exactly == gap_threshold must not count; the test is strict `>`.
    t = np.array([0.0, 1.0, 2.0, 12.0])  # median 1.0 -> threshold 10.0, last step 10.0
    _, gap_threshold_s, n_gaps, _ = windows.spacing(t, gap_mult=10.0)
    assert gap_threshold_s == 10.0
    assert n_gaps == 0


# -------------------- what does and does not move a boundary (8.4a)


PRIVATE_ROOT = repo_root() / "data" / "real_private"


# READ from the jobs, not transcribed. A hardcoded copy is not a guard: it agreed with
# itself while either job drifted, which is the dead-control defect AGENTS.md section 4
# records. `GAP_MULT` is the parameter that actually moves a boundary.
def _job_gap_mult(stem: str) -> float:
    path = repo_root() / "jobs" / "active" / f"{stem}.py"
    if not path.is_file():
        pytest.skip(f"{path} absent")
    text = path.read_text(encoding="utf-8")
    found = re.search(r"gap_mult=([0-9.]+)", text) or re.search(
        r"^GAP_MULT\s*=\s*([0-9.]+)", text, re.M
    )
    assert found, f"no gap_mult literal found in {stem}.py"
    return float(found.group(1))


def _poster_carve() -> dict:
    return {
        "gap_mult": _job_gap_mult("km_poster_6d2s"),
        "k": 1.0,
        "use_uncertainty": False,
    }


def _survey_carve() -> dict:
    return {
        "gap_mult": _job_gap_mult("km_with_checks_6d2s"),
        "k": 1.0,
        "use_uncertainty": True,
    }


THRESHOLD_LABEL = "3.0 µs"
THRESHOLD = [(THRESHOLD_LABEL, 3.0e-6, True)]


def _synthetic_reads(n: int = 400):
    """Reads that cross the rung many times, with a sigma large enough to make most of them
    'uncertain' if uncertainty were allowed to move anything."""
    rng = np.random.default_rng(11)
    t = np.arange(n, dtype=float) * 30.0
    values = 3.0e-6 + 1.2e-6 * np.sin(np.arange(n) / 7.0) + rng.normal(0, 3e-7, n)
    sigma = np.full(n, 8e-7)
    return t, values, sigma


def _carve_annotated(*, use_uncertainty: bool):
    t, values, sigma = _synthetic_reads()
    return windows.run(
        windows.WindowsInputs(
            t_rel_s=t,
            values=values,
            thresholds=THRESHOLD,
            sigma=sigma if use_uncertainty else None,
            k=1.0,
            gap_mult=10.0,
            use_uncertainty=use_uncertainty,
        )
    )


@pytest.mark.unit
def test_uncertainty_annotates_reads_and_never_moves_a_window_boundary():
    """The contract, on synthetic reads, in CI. Both halves are asserted: the window table
    must be identical AND the read states must actually differ, or the test would pass
    vacuously on a sigma too small to annotate anything."""
    plain = _carve_annotated(use_uncertainty=False)
    annotated = _carve_annotated(use_uncertainty=True)

    wp = plain.windows_in_spec.reset_index(drop=True)
    wa = annotated.windows_in_spec.reset_index(drop=True)
    assert len(wp) > 5, "the fixture must produce several windows to be worth comparing"
    assert wp.equals(wa), (
        "use_uncertainty moved a window boundary; it must only annotate"
    )

    states_plain = set(plain.reads["state"])
    states_annotated = set(annotated.reads["state"])
    assert states_plain != states_annotated, (
        "the fixture's sigma is too small to annotate anything, so the test above proved "
        f"nothing: {states_plain} vs {states_annotated}"
    )
    assert any("uncertain" in s for s in states_annotated)


@pytest.mark.integration
@pytest.mark.real
@pytest.mark.parametrize(
    "stem,qubit",
    [
        ("280623_6D2S_qubit2", 2),
        ("040423_6D2S_qubit1", 1),
        ("220423_6D2S_qubit1", 1),
        ("090623_6D2S_qubit6", 6),
        ("070723_6D2S_qubit4", 4),
    ],
)
def test_the_two_job_configurations_agree_on_the_real_records(stem, qubit, in_repo):
    """The differential CHECKPOINT 8.4a owes: real records, the real parameter sets.

    Skips without the private tree, following `tests/test_data_manifest.py`. If this ever
    fails, the two jobs have diverged on something that DOES move boundaries and the new job's
    check outcome would describe a different carve from the band beside it.
    """
    path = PRIVATE_ROOT / "6D2S" / f"{stem}.pickle"
    if not path.is_file():
        pytest.skip(f"{path} absent (expected without the private data)")

    from quebra.analyzers import t2star
    from quebra.core.dataset import Dataset
    from quebra.core.job import _load_dataset
    from quebra.recipes import RAMSEY_CONFIG, _filter_step, _final_stage
    from quebra.schemas.track912 import track912Schema

    dataset = Dataset(
        path=f"data/real_private/6D2S/{stem}.pickle",
        schema=track912Schema,
        qubit=qubit,
        device="6D2S",
        extra={"run_name": stem},
    )
    norm = _final_stage(_filter_step(RAMSEY_CONFIG)(_load_dataset(dataset)))
    frame = t2star.run(t2star.make_inputs_from_norm(norm)).frame

    def carve(config):
        return windows.run(
            windows.make_inputs_from_frame(
                frame,
                time_col="t_rel_s",
                value_col="t2star_s",
                thresholds=THRESHOLD,
                dataset_id=stem,
                sigma_col="t2star_error_s" if config["use_uncertainty"] else None,
                **config,
            )
        ).windows_in_spec

    poster = carve(_poster_carve()).reset_index(drop=True)
    survey = carve(_survey_carve()).reset_index(drop=True)

    assert len(poster), f"{stem} produced no windows at {THRESHOLD_LABEL}"
    assert poster.equals(survey), (
        f"{stem}: the poster and survey carve configurations disagree on the window table. "
        f"That is a divergence in something that moves boundaries, not an annotation."
    )


_STEPS = st.sampled_from([0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 50.0])
_LEVELS = st.sampled_from([1.0, 4.0, 5.0, 6.0, 9.0])


@pytest.mark.properties
@given(
    steps=st.lists(_STEPS, min_size=1, max_size=40),
    levels=st.lists(_LEVELS, min_size=41, max_size=41),
    big_values_good=st.booleans(),
)
def test_both_sides_and_the_gaps_tile_the_record(steps, levels, big_values_good):
    """Oracle: arithmetic on the record itself, computed here without the carve's help.

    Duplicate timestamps (step 0), reads exactly at the threshold (level 5) and gaps (step
    50 against a median spacing near 1) are all drawn. The expected windows are the
    maximal runs of reads on one side with no gap between neighbours, so their [s, e)
    ranges partition the reads. Each carved row must match its run read for read, the
    read table must point in-spec reads at their run and no other read anywhere, the
    windows and gap spans must cover the record end to end, and each side's durations
    must equal the observed time charged to that side's reads, each interval charged to
    the read that opens it.
    """
    assume(any(step > 0 for step in steps))
    t = np.r_[0.0, np.cumsum(steps)]
    v = np.asarray(levels[: len(t)])
    n = len(t)
    result = windows.run(
        windows.WindowsInputs(
            t_rel_s=t, values=v, thresholds=[("5", 5.0, big_values_good)]
        )
    )
    w_in, w_out = result.windows_in_spec, result.windows_out_of_spec
    assert set(w_in["side"]) <= {windows.SIDE_IN_SPEC}
    assert set(w_out["side"]) <= {windows.SIDE_OUT_OF_SPEC}

    dt = np.diff(t)
    gap_threshold_s = windows.DEFAULT_GAP_MULT * float(np.median(dt[dt > 0]))
    is_gap = np.r_[False, dt > gap_threshold_s]
    margin_v = v - 5.0 if big_values_good else 5.0 - v
    in_spec = margin_v >= 0.0
    starts = [
        i for i in range(n) if i == 0 or is_gap[i] or in_spec[i] != in_spec[i - 1]
    ]
    runs = list(zip(starts, [*starts[1:], n], strict=True))

    nan = float("nan")
    for frame, side in ((w_in, True), (w_out, False)):
        mine = [(s, e) for s, e in runs if bool(in_spec[s]) == side]
        born = [s > 0 and not is_gap[s] for s, _ in mine]
        dies = [e < n and not is_gap[e] for _, e in mine]
        expected = {
            "n_reads": [e - s for s, e in mine],
            "t_birth_s": [t[s] for s, _ in mine],
            "t_last_s": [t[e - 1] for _, e in mine],
            "t_before_birth_s": [
                t[s - 1] if b else nan for (s, _), b in zip(mine, born)
            ],
            "t_death_s": [t[e] if d else t[e - 1] for (_, e), d in zip(mine, dies)],
            "birth_observed": born,
            "censored": [not d for d in dies],
            "extreme_margin": [float(np.min(margin_v[s:e])) for s, e in mine],
        }
        assert len(frame) == len(mine)
        for column, values in expected.items():
            np.testing.assert_array_equal(
                frame[column].to_numpy(), np.asarray(values), err_msg=column
            )

    in_runs = [(s, e) for s, e in runs if in_spec[s]]
    expected_index = np.full(n, -1)
    for index, (s, e) in enumerate(in_runs):
        expected_index[s:e] = index
    np.testing.assert_array_equal(
        result.reads["window_index"].fillna(-1).to_numpy(dtype=int), expected_index
    )

    gaps = result.diagnostics["gap_spans_s"]
    assert gaps == [(t[i - 1], t[i]) for i in np.flatnonzero(is_gap)]
    covered = math.fsum(
        [*w_in["duration_s"], *w_out["duration_s"], *(b - a for a, b in gaps)]
    )
    assert covered == pytest.approx(t[-1] - t[0], rel=1e-12, abs=1e-12)

    observed = ~is_gap[1:]
    for frame, on_side in ((w_in, in_spec), (w_out, ~in_spec)):
        charged = math.fsum(dt[observed & on_side[:-1]])
        assert math.fsum(frame["duration_s"]) == pytest.approx(charged, abs=1e-9)
