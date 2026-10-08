"""Rendering the same figure twice produces the same bytes.

Oracle: the bytes of the written files, inspected for the variable fields themselves.

WHY THIS IS WORTH A TEST. `output/` is append-only and holds hundreds of PDFs, so once the
bytes are a function of the drawing alone they become a regression signal for free: a figure
that changed without its inputs changing shows up as a differing hash. That only works if the
writer stops stamping the time and the renderer's own version into every file.

The three registered matplotlib targets are all covered, and the poster is the one that
matters most here: PNG carries none of PDF's metadata keys, so the obvious fix applied to all
three would have been a silent no-op on exactly that target.

WHY THIS ASSERTS THE FIELDS AND NOT TWO RENDERS AGREEING. A same-machine render-twice
comparison cannot see any of this. `/CreationDate` has one-second resolution, so two renders
inside the same second match whether or not it is suppressed; and the PNG's `Software` chunk
names the matplotlib version, which is constant on one machine. Measured: of the four
metadata mutations, such a test caught one reliably, one in five runs of six, and two never.
What it would have guarded uniquely is hash-ordering nondeterminism, and that lives inside
matplotlib where nothing in this repository can break it.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration

REPO_ROOT = Path(__file__).resolve().parents[1]


def _poster_font_or_skip() -> None:
    """Skip where the poster's required face is absent, the way the `r` tier skips.

    The face is an environment requirement, not a Python dependency: it cannot be pip
    installed and CI does not provision it. A test that FAILED without it would turn the
    default test run red on every machine but the one that has the font, which is the opposite
    of what a guard against silent substitution should cost.
    """
    from quebra.plots import theme

    try:
        theme.check_fonts("poster")
    except theme.MissingFontError as exc:
        pytest.skip(str(exc))


# Renders one target into a directory and prints nothing. Kept as a string so each render is
# a genuinely fresh interpreter rather than a subprocess of this one's import state.
_RENDER = """
import pathlib
import sys
import matplotlib
matplotlib.use("Agg")
from quebra.plots.targets import RENDER_TARGETS
from quebra.plots.base import BasePlot
import matplotlib.pyplot as plt


class _Tiny(BasePlot):
    def build_matplotlib(self, result, style="default"):
        fig, ax = plt.subplots(figsize=(3.0, 2.0))
        ax.plot([0, 1, 2], [2, 0, 1])
        ax.set_xlabel("Window age (h)")
        return fig


RENDER_TARGETS[sys.argv[1]](
    _Tiny(name="determinism_probe"), object(), pathlib.Path(sys.argv[2])
)
"""


def _render(target: str, out_dir: Path, seed: str = "0") -> bytes:
    if target == "poster":
        _poster_font_or_skip()
    out_dir.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        [sys.executable, "-c", _RENDER, target, str(out_dir)],
        cwd=REPO_ROOT,
        env={**os.environ, "PYTHONHASHSEED": seed},
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, (
        f"rendering {target} failed.\nstdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}"
    )
    written = sorted(p for p in out_dir.iterdir() if p.is_file())
    assert len(written) == 1, f"expected one file for {target}, got {written}"
    return written[0].read_bytes()


def test_no_written_file_names_the_renderer_or_the_time(tmp_path):
    """Oracle: the bytes of the written files, inspected for the fields themselves.

    The cross-process test above cannot see this. Both variable fields it guards against are
    CONSTANT within one machine and one matplotlib: the PNG's `Software` chunk names the
    matplotlib version, so two runs here agree whether or not it is suppressed, and the test
    passes either way. Measured: swapping the PNG's key for the PDF's leaves it green.

    What actually breaks is an upgrade, which no same-day test can stage. So this asserts the
    property directly rather than its consequence, and it is what makes the poster covered.
    """
    pdf_forbidden = (b"/CreationDate", b"Matplotlib")
    for target, forbidden in (
        ("static", pdf_forbidden),
        ("academic", pdf_forbidden),
        ("poster", (b"Software",)),
    ):
        written = _render(target, tmp_path / target)
        for token in forbidden:
            assert token not in written, (
                f"the `{target}` file carries {token!r}. That field names the machine's "
                f"clock or the renderer's version, so the bytes will move on an upgrade "
                f"even when the drawing does not."
            )
        if target != "poster":
            assert b"quebra" in written, (
                f"the `{target}` file should carry the pinned Creator/Producer instead"
            )
