"""Ratchet on hardcoded style literals in the render layer.

Colour and typography belong in plots/theme.py. Sweeping every existing call site in
one pass is not worth the churn, so this pins the current count instead: the number
may go DOWN, never up. A new panel that hardcodes its own hex or font sizes fails here.

Scope: a 6-digit hex literal, or a font size given as a LITERAL - numeric (`fontsize=9`) or
named (`fontsize="x-small"`) - in either the keyword or the dict spelling, and either quote
style. It does not catch 3- or 8-digit hex, named colours ("lightgray"), plt.cm.*, or
size=/labelsize=.

WHY THE VALUE MUST BE A LITERAL, not the key. `fontsize=theme.CAPTION["fontsize"]` is the
behaviour this test exists to encourage, so matching the bare key would count it as a
violation. Keying on the value instead means an identifier is never counted and a literal
always is, whichever way it is spelled. Measured on the tree as it stands: the two dict-form
sites both read from the theme and are correctly not counted, while `allan_plot.py`'s two
`fontsize="x-small"` are hardcoded and are.
"""

from __future__ import annotations

import re
from pathlib import Path
import pytest

pytestmark = pytest.mark.policy


REPO_ROOT = Path(__file__).resolve().parents[1]

# plots/theme.py is the one place these literals are supposed to live.
THEME_FILE = REPO_ROOT / "src" / "quebra" / "plots" / "theme.py"

PATTERNS = {
    "hex": re.compile(r"#[0-9A-Fa-f]{6}\b"),
    # Both spellings, both quote styles, and only when the VALUE is a literal.
    "fontsize": re.compile(r"""fontsize\s*=\s*["'\d]|["']fontsize["']\s*:\s*["'\d]"""),
}

# It ratchets DOWN only. 15 is the current floor. The pattern keys on any literal VALUE, so
# a named size such as `fontsize="x-small"` counts as the hardcoded size it is, and it reaches
# the dict spelling `"fontsize": 9` as well as the keyword; a size read from `theme` does not
# count.
BASELINE = 15


def _counts() -> dict[str, int]:
    counts: dict[str, int] = {}
    for directory in ("src/quebra/panels", "src/quebra/plots"):
        for path in sorted((REPO_ROOT / directory).glob("*.py")):
            if path == THEME_FILE:
                continue
            source = path.read_text(encoding="utf-8")
            n = sum(len(pattern.findall(source)) for pattern in PATTERNS.values())
            if n:
                counts[str(path.relative_to(REPO_ROOT))] = n
    return counts


def test_hardcoded_style_literals_do_not_grow() -> None:
    counts = _counts()
    total = sum(counts.values())
    breakdown = "\n".join(f"  {name}: {n}" for name, n in sorted(counts.items()))
    assert total <= BASELINE, (
        f"Hardcoded style literals rose to {total}, above the pinned baseline "
        f"{BASELINE}. This number ratchets DOWN, never up - move the new colour or "
        f"font size into plots/theme.py rather than raising the baseline.\n{breakdown}"
    )


def test_baseline_is_tight() -> None:
    # If the count has dropped, re-pin it. A baseline left slack above the real count
    # silently re-opens the room it was meant to close.
    total = sum(_counts().values())
    assert total == BASELINE, (
        f"Hardcoded style literals are down to {total}; re-pin BASELINE to {total} "
        f"in this file to keep the ratchet tight."
    )
