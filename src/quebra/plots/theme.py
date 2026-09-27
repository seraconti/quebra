"""Single source of colour and typography for the render layer.

Every data-derived quantity lives in the typed artifact; this module owns only the
functions-of-axes/theme half of that split. Constants and small helpers - no config
layer, no second theme module.

All colour accessors return a hex `str`. `mix_with_white` is the one exception and
predates the rule: it returns an RGB float tuple, so helpers built on it convert back.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from contextlib import contextmanager
from typing import Iterator

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors as mcolors
from matplotlib import font_manager


QUBIT_COLOR_MAP = {
    1: "#438ECB",
    2: "#E91B63",
    3: "#FDC113",
    4: "#1D988A",
    5: "#CDDC39",
    6: "#8F3F97",
}

# Ladder colours are sampled from one perceptually ordered map rather than taken from
# a fixed list, so a threshold ladder of ANY length gets distinct colours. The previous
# 8-entry list wrapped modulo its length and drew the 1 us and 9 us thresholds in the
# same blue on a 10-entry ladder.
THRESHOLD_CMAP = "cividis"

# Sampling is clipped away from both ends of the map: cividis bottoms out near-black
# and tops out in a yellow that is hard to see on white.
_THRESHOLD_SPAN = (0.05, 0.90)

# Hue for a dataset whose qubit id could not be resolved.
_UNKNOWN_QUBIT_COLOR = "tab:gray"

# Spec-state semantics. TWO CHANNELS, deliberately: HUE carries which side of the
# threshold a read fell (teal in, red out), LIGHTNESS carries whether the per-read
# error could resolve it (deep = resolved, pale = uncertain). Relative luminances,
# computed not asserted:
#     in_spec                #12776A  L=0.144    in_spec_uncertain      L=0.646
#     out_of_spec            #7E1A11  L=0.052    out_of_spec_uncertain  L=0.537
#     unresolved             #9B8F63  L=0.275
#     censored               #4D4D4D  L=0.074
#     unobserved             #FFFFFF  L=1.000  (readable only via its hatch)
# The within-hue gap is ~0.50, up from ~0.21 when the crisp states were mid-tone -
# at timeline bar widths the old pair was not tellable apart.
# TRADE-OFF, stated rather than hidden: the two DEEP states are close in greyscale
# (gap ~0.09), so a greyscale print separates certain from uncertain but leans on hue
# for in from out. Lightness cannot carry both without one of them becoming arbitrary.
STATE_COLORS = {
    "in_spec": "#12776A",
    "out_of_spec": "#7E1A11",
    "unresolved": "#9B8F63",
    "censored": "#4D4D4D",
    # Unobserved time (a read gap). Not a state of the system - a state of the record,
    # so it must not be readable as a pale version of any state. Drawn white with a
    # hatch (see UNOBSERVED_HATCH): the pale certainty tones sit at L~0.55-0.65, close
    # enough to a flat grey that colour alone was confusable.
    "unobserved": "#FFFFFF",
}

# A read gap is drawn as hatching, not a fill: no colour can be far enough from four
# state colours to be unmistakable, but a texture is categorically different.
UNOBSERVED_HATCH = "////"
UNOBSERVED_EDGE = "#9A9A9A"

# Confidence/tolerance bands: the fill carries the extent, the edge carries the method.
BAND_STYLE = {
    "fill_alpha": 0.18,
    "edge_alpha": 0.55,
    "edge_linewidth": 0.8,
}

# Censoring marks on a survival/duration axis.
CENSOR_TICK = {
    "marker": "|",
    "markersize": 5.0,
    "linewidth": 0.9,
}

# The region where the risk set has fallen below the point of being informative.
EPSILON_SHADE = {
    "color": "#B0B0B0",
    "alpha": 0.20,
}

# The poster's typeface, primary first. Roboto is not shipped with matplotlib and is not a
# Python dependency: it is an environment requirement, which is exactly why a missing one has
# to be reported rather than absorbed.
POSTER_FONT_STACK = ["Roboto", "DejaVu Sans"]

# Which targets REQUIRE their first-choice face, named explicitly. Keying the check on this
# rather than on whether `font.family` happens to be a list matters: the poster's previous
# spelling was a bare string, so a type test would let one reverted character disable the
# guard with no signal, which is the class of silence the guard exists to end.
REQUIRED_FIRST_FACE = {"poster": POSTER_FONT_STACK[0]}


# rcParams by RENDER TARGET, not by matplotlib style name. `static` keeps the sans
# default; `academic` goes serif at a smaller base for print.
RCPARAMS: dict[str, dict[str, object]] = {
    "static": {
        "font.family": "sans-serif",
        "font.size": 9.0,
        "axes.titlesize": 10.0,
        "axes.labelsize": 9.0,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "legend.fontsize": 8.0,
        "figure.titlesize": 12.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": "#CCCCCC",
        "grid.linewidth": 0.5,
        "grid.alpha": 0.4,
    },
    "academic": {
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "font.size": 8.0,
        "axes.titlesize": 9.0,
        "axes.labelsize": 8.0,
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.0,
        "legend.fontsize": 7.0,
        "figure.titlesize": 10.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": "#CCCCCC",
        "grid.linewidth": 0.4,
        "grid.alpha": 0.4,
    },
    "poster": {
        # A LIST, not a name, and the list is the contract. matplotlib silently falls back
        # to its default when a family is missing, so a poster rendered on a machine without
        # the first face is a different figure that looks like the same one. `check_fonts`
        # below refuses that rather than absorbing it.
        "font.family": POSTER_FONT_STACK,
        "font.size": 28.0,
        # Hierarchy
        "axes.titlesize": 30.0,
        "axes.titleweight": "bold",
        "axes.labelsize": 24.0,
        "axes.labelweight": "bold",
        # Ticks
        "xtick.labelsize": 20.0,
        "ytick.labelsize": 20.0,
        # Legend / figure title
        "legend.fontsize": 20.0,
        "figure.titlesize": 32.0,
        # Lines
        "axes.linewidth": 2.0,
        "lines.linewidth": 3.2,
        # Tick marks
        "xtick.major.width": 2.0,
        "ytick.major.width": 2.0,
        "xtick.major.size": 9.0,
        "ytick.major.size": 9.0,
        # Clean axes
        "axes.spines.top": False,
        "axes.spines.right": False,
        # Subtle grid: useful, but stays in the background
        "axes.grid": True,
        "grid.color": "#CCCCCC",
        "grid.linewidth": 0.8,
        "grid.alpha": 0.30,
    },
}


class MissingFontError(RuntimeError):
    """The configured typeface is not installed, so the figure would silently differ."""


def check_fonts(target: str) -> None:
    """Raise if this target's first-choice family is not installed.

    matplotlib reports a missing family by logging a warning and substituting its default,
    which means the figure still renders and still looks plausible. For a poster that is the
    wrong-but-plausible result AGENTS.md section 3 is about: the same script produces
    different typography on two machines and nothing says so.

    Only the FIRST family is required, and only for the targets named in
    `REQUIRED_FIRST_FACE`. The rest of the stack exists so a reader who cannot install it has
    a stated fallback rather than an arbitrary one.
    """
    primary = REQUIRED_FIRST_FACE.get(target)
    if primary is None:
        return
    try:
        font_manager.findfont(primary, fallback_to_default=False)
    except ValueError as exc:
        raise MissingFontError(
            f"the {target!r} target asks for {primary!r} and it is not installed, so "
            f"matplotlib would substitute its default and the figure would differ from the "
            f"one this project draws. Install it, or change POSTER_FONT_STACK in "
            f"plots/theme.py to a face this machine has."
        ) from exc


# Renderers are handed matplotlib style names by plots/targets.py; the theme is keyed
# by render target. One mapping, here, so no renderer has to know both vocabularies.
_TARGET_FOR_STYLE = {"default": "static", "paper": "academic", "poster": "poster"}


def extract_qubit_id(
    dataset_id: str | None = None, meta: dict[str, object] | None = None
) -> int | None:
    if isinstance(meta, dict) and meta.get("qubit") is not None:
        try:
            qubit = int(meta.get("qubit"))
            if qubit > 0:
                return qubit
        except (TypeError, ValueError):
            pass

    dataset = str(dataset_id) if dataset_id is not None else ""
    patterns = [r"(?:^|[_-])q(\d+)(?:[_-]|$)", r"qubit(\d+)"]
    for pattern in patterns:
        match = re.search(pattern, dataset, flags=re.IGNORECASE)
        if match:
            qubit = int(match.group(1))
            if qubit > 0:
                return qubit
    return None


def qubit_color(
    dataset_id: str | None = None, meta: dict[str, object] | None = None
) -> str:
    qubit = extract_qubit_id(dataset_id=dataset_id, meta=meta)
    if qubit in QUBIT_COLOR_MAP:
        return QUBIT_COLOR_MAP[qubit]
    return mcolors.to_hex(_UNKNOWN_QUBIT_COLOR)


def mix_with_white(color: str, amount: float = 0.4) -> tuple[float, float, float]:
    base = np.asarray(mcolors.to_rgb(color), dtype=float)
    amount = float(np.clip(amount, 0.0, 1.0))
    mixed = (1.0 - amount) * base + amount * np.ones(3, dtype=float)
    return tuple(mixed.tolist())


# Ramp stops short of white so the last dataset in a long series is still legible.
_RAMP_MAX_MIX = 0.55


def dataset_ramp(qubit: int | None, index: int, n: int) -> str:
    """Chronological lightness ramp within one qubit's hue.

    `index` is the dataset's position in time, `n` the number of datasets sharing the
    hue. Earliest is the saturated base colour, latest the most washed-out.
    """
    base = QUBIT_COLOR_MAP.get(qubit, _UNKNOWN_QUBIT_COLOR)
    if n <= 1:
        return mcolors.to_hex(base)
    position = float(np.clip(index / (n - 1), 0.0, 1.0))
    return mcolors.to_hex(mix_with_white(base, _RAMP_MAX_MIX * position))


def ordered_color(i: int, n: int) -> str:
    """Colour for item `i` of an ORDERED sequence of `n`, sampled across `THRESHOLD_CMAP`.

    For any quantity whose levels have a natural order and should read as a progression -
    a threshold ladder, an event-count sweep. A categorical palette would imply the levels
    are unrelated; this one makes "more" look like "further along the ramp".
    """
    lo, hi = _THRESHOLD_SPAN
    if n <= 1:
        position = lo
    else:
        position = lo + (hi - lo) * float(np.clip(i / (n - 1), 0.0, 1.0))
    return mcolors.to_hex(plt.get_cmap(THRESHOLD_CMAP)(position))


def threshold_color(i: int, n: int) -> str:
    """Colour for threshold `i` of a ladder of `n`. The ladder's name for `ordered_color`."""
    return ordered_color(i, n)


# Calibration method. Two ways of turning a statistic into a p-value, and the bench exists
# to show they are not interchangeable, so they must never be drawn in one colour. Warm =
# an approximation that degrades at small n; cool = exact under permutation.
CALIBRATION_COLORS = {
    "asymptotic": "#B4553C",
    "permutation": "#2C6E8F",
}


def calibration_color(calibration: str) -> str:
    if calibration not in CALIBRATION_COLORS:
        raise ValueError(
            f"Unknown calibration: {calibration!r}. Known: {sorted(CALIBRATION_COLORS)}"
        )
    return CALIBRATION_COLORS[calibration]


# Ledger verdicts. The palette encodes a THREE-WAY distinction, not a two-way one, because
# the whole point of the ledger is that "did not reject" is not "passed":
#   pass / fail          the check ran, was calibrated at this event count, and answered
#   underpowered         it answered, but the answer carries no information
#   not interpretable    ties or a missing dependency mean there is no answer at all
# The two "no information" verdicts are deliberately DESATURATED so they cannot be mistaken
# for a result at a glance; a reader skimming the table should see the grey and stop.
#   no row               `check_ledger.VERDICT_ABSENT`: no row for that cell at all. The
#                        ledger never emits it; a downstream output-builder does.
# The last one is separate from `not computed` on purpose. "The check ran and declined to
# answer" and "there is nothing here" look identical to a reader who is shown one tone for
# both, and SPEC 0008 R8.3 requires them to be distinguishable - a cell that was asked for and
# answered nothing must not read the same as one nobody asked for. It is the palest tone
# because it carries the least.
VERDICT_COLORS = {
    "pass": "#12776A",
    "fail": "#7E1A11",
    "underpowered": "#B9A87A",
    "not interpretable (ties)": "#9A9A9A",
    "not computed": "#D8D8D8",
    "no row": "#F0F0F0",
}


def verdict_color(verdict: str) -> str:
    if verdict not in VERDICT_COLORS:
        raise ValueError(
            f"Unknown verdict: {verdict!r}. Known: {sorted(VERDICT_COLORS)}"
        )
    return VERDICT_COLORS[verdict]


# A line the reader compares AGAINST rather than reads: nominal alpha, the 45-degree
# diagonal of a P-P plot, a y = x reference. Deliberately unsaturated - it is the ruler,
# not a series.
REFERENCE_LINE = {
    "color": "#5A5A5A",
    "linestyle": "--",
    "linewidth": 1.0,
    "alpha": 0.75,
}

# Where the REAL data sits on a swept axis. Shaded rather than outlined so it reads as a
# region of interest behind the curves instead of competing with them.
HIGHLIGHT_SHADE = {
    "color": "#C8A85A",
    "alpha": 0.22,
}


# How far an "uncertain" state is washed toward white. Enough to read as a distinct,
# weaker claim next to its crisp sibling without losing which sibling it is.
_UNCERTAIN_MIX = 0.72


def state_color(state: str) -> str:
    """Colour for a per-read spec state, including the two uncertain tones.

    `in_spec_uncertain` / `out_of_spec_uncertain` are their crisp siblings washed toward
    white, so a 4-state timeline still reads as green/red at a glance while marking
    which reads the error bars cannot separate from the threshold.

    `unresolved` is a different thing and keeps its own hue: it means the read cannot be
    classified at all (no usable sigma), not that it sits near the threshold.
    """
    if state in STATE_COLORS:
        return STATE_COLORS[state]
    if state.endswith("_uncertain"):
        base = state.removesuffix("_uncertain")
        if base in STATE_COLORS:
            return mcolors.to_hex(mix_with_white(STATE_COLORS[base], _UNCERTAIN_MIX))
    raise ValueError(
        f"Unknown spec state: {state!r}. Known: "
        f"{sorted(STATE_COLORS)} plus their _uncertain variants."
    )


def apply_rcparams(target: str) -> None:
    """Apply this render target's rcParams to the live matplotlib state.

    Call inside a style context - see `style_context`. Called outside one it still
    mutates global rcParams, but the effect is DISCARDED the moment a renderer enters
    its stylesheet, because entering one resets rcParams wholesale. That is why this
    is not called from plots/targets.py's render functions directly.
    """
    if target not in RCPARAMS:
        raise ValueError(
            f"Unknown render target: {target!r}. Known: {sorted(RCPARAMS)}"
        )
    plt.rcParams.update(RCPARAMS[target])


@contextmanager
def style_context(style: str) -> Iterator[None]:
    """The one stylesheet every renderer draws inside.

    Opens matplotlib's `default` stylesheet, then layers this project's rcParams on
    top. Ordering is load-bearing: the stylesheet must be entered FIRST, because
    entering one resets rcParams to that sheet's values and would otherwise discard
    the theme. `plt.style.context` restores the previous state on exit, so rendering
    several targets in one run cannot leak settings between them.
    """
    target = _TARGET_FOR_STYLE.get(style)
    if target is None:
        raise ValueError(
            f"Unknown style: {style!r}. Known: {sorted(_TARGET_FOR_STYLE)}"
        )
    check_fonts(target)
    with plt.style.context("default"):
        apply_rcparams(target)
        yield


# --- Report typography and neutrals ------------------------------------------------
# Added for the instrument-report figures, which carry more prose than a data panel does:
# a caption stating what is excluded, per-row justifications, and a methodological note
# beside the axis it qualifies. Those need consistent, named styles rather than a hex and
# a font size at each call site - which is exactly what tests/test_style_baseline.py
# ratchets against, and the ratchet only moves down.

# The small print under a figure: what was excluded, what this does not claim.
CAPTION = {"fontsize": 7, "color": "#555555"}

# A note attached to one axis rather than to the whole figure.
ANNOTATION = {"fontsize": 7, "color": "#444444"}

# Text drawn ON a filled cell whose fill is already carrying the meaning.
ON_FILL_TEXT = {"fontsize": 8, "color": "#333333"}

# Tick and legend text in the report figures, which run denser than the panels. Exposed as
# DICTS to be unpacked, not as bare sizes: the ratchet counts the string `fontsize=` at a
# call site regardless of what it points at, so `fontsize=theme.LABEL_SIZE` would still
# fail it. Unpacking keeps the literal here, which is the rule the ratchet encodes.
LABEL_TEXT = {"fontsize": 8}
LEGEND_TEXT = {"fontsize": 7}

# Neutral connectors: the line joining a matched pair, and the leader from an axis to a
# marker. They carry NO meaning of their own - the marker positions do - so they sit below
# every series in weight as well as in z-order.
PAIR_CONNECTOR = {"color": "#BBBBBB", "linewidth": 1.2}
LEADER_LINE = {"color": "#CCCCCC", "linewidth": 1.0}


# Poster ink. A poster figure is read as ONE object from a distance, not as a panel in a
# series, so it drops the per-qubit identity palette and carries a single ink per figure;
# series inside it separate by line style, which survives both distance and greyscale.
# Author-chosen, not derived - these two are the poster's own colours.
POSTER_INK_SURVIVAL = "#12293f"
POSTER_INK_DISTRIBUTION = "#8c2231"

# A label naming a marked location on an axis (a mean, a median). Bold and near-black:
# it names a value the figure is ABOUT, unlike ANNOTATION, which qualifies one.
GUIDE_LABEL = {"fontsize": 7, "color": "#1A1A1A", "fontweight": "bold"}


# The text styles above are ABSOLUTE point sizes, tuned for one panel among several on a
# page. The poster target draws a single figure at roughly three times that scale, where a
# 7-point caption is not small - it is invisible. Renderers therefore ask for a style by
# name and get it sized for the target they are drawing into, rather than each one picking
# its own poster size (which is what the ratchet in tests/test_style_baseline.py exists to
# prevent). rcParams-driven text - titles, axis labels, ticks, legend - scales through
# RCPARAMS and never comes through here.
_TEXT_SCALE = {"static": 1.0, "academic": 1.0, "poster": 2.2}


def scaled_text(base: Mapping[str, object], style: str) -> dict[str, object]:
    """One of the named text styles, sized for the render target behind `style`.

    `style` is the matplotlib style name a renderer is handed (`default`, `paper`,
    `poster`), the same vocabulary `style_context` takes, so a renderer never has to know
    the target name.
    """
    target = _TARGET_FOR_STYLE.get(style)
    if target is None:
        raise ValueError(
            f"Unknown style: {style!r}. Known: {sorted(_TARGET_FOR_STYLE)}"
        )
    out = dict(base)
    if "fontsize" in out:
        out["fontsize"] = float(out["fontsize"]) * _TEXT_SCALE[target]
    return out
