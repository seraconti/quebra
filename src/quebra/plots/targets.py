"""Render targets: how a figure sink writes a built figure to disk, registered by name.

Three are registered. `static` writes `{name}_static.pdf` in the default style, `academic`
writes `{name}_academic.pdf` in the paper style, and `poster` writes a 600 dpi PNG in the
poster style. A job names targets in `job.figure(..., targets=[...])`, and the runner refuses
an unregistered name before any step runs. A new target is one decorated function.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Protocol

import matplotlib.pyplot as plt

from quebra.plots import theme

# Written into every PDF instead of matplotlib's own version string, and the PNG's `Software`
# key is suppressed for the same reason: those fields name the matplotlib that rendered the
# file, so leaving them in makes the bytes move on every matplotlib release even when the
# drawing is identical. `CreationDate` is the other variable field and is dropped outright.
#
# The KEYS DIFFER BY FORMAT. PDF has CreationDate/Creator/Producer; PNG has none of them and
# seeds `Software` instead, so passing a PDF key to the PNG writer is a silent no-op.
_PDF_METADATA = {"CreationDate": None, "Creator": "quebra", "Producer": "quebra"}
_PNG_METADATA = {"Software": None}


class _RenderablePlot(Protocol):
    name: str

    def build_matplotlib(
        self, result: object, style: str = "default"
    ) -> plt.Figure: ...


RenderFn = Callable[[_RenderablePlot, object, Path], None]

RENDER_TARGETS: dict[str, RenderFn] = {}


def register_target(name: str) -> Callable[[RenderFn], RenderFn]:
    target_name = name.strip().lower()
    if not target_name:
        raise ValueError("Target name must be a non-empty string")

    def decorator(fn: RenderFn) -> RenderFn:
        existing = RENDER_TARGETS.get(target_name)
        if existing is not None and existing is not fn:
            raise ValueError(f"Target already registered: {target_name}")
        RENDER_TARGETS[target_name] = fn
        return fn

    return decorator


@register_target("static")
def render_static(plot: _RenderablePlot, result: object, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    # savefig is inside the style context too: a few rcParams (mathtext.fontset, the
    # font.serif fallback order) resolve at DRAW time, not when the artist is created,
    # so saving outside the context would silently drop them.
    with theme.style_context("default"):
        figure = plot.build_matplotlib(result, style="default")
        # Per-target filename: static and academic would otherwise write the same {name}.pdf,
        # so academic clobbered static while the prov record listed both targets.
        figure.savefig(
            out_dir / f"{plot.name}_static.pdf",
            dpi=300,
            bbox_inches="tight",
            metadata=_PDF_METADATA,
        )
    plt.close(figure)


@register_target("academic")
def render_academic(plot: _RenderablePlot, result: object, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with theme.style_context("paper"):
        figure = plot.build_matplotlib(result, style="paper")
        figure.savefig(
            out_dir / f"{plot.name}_academic.pdf",
            dpi=600,
            bbox_inches="tight",
            metadata=_PDF_METADATA,
        )
    plt.close(figure)


@register_target("poster")
def render_poster(plot: _RenderablePlot, result: object, out_dir: Path) -> None:
    """One figure, printed large: PNG at 600 dpi on the poster typography.

    PNG rather than PDF because this target exists for figures that go into a poster
    layout tool or a slide, which want a raster they can place. 600 dpi keeps a 10-inch
    figure legible at A0. `static` and `academic` stay PDF and are unaffected.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    with theme.style_context("poster"):
        figure = plot.build_matplotlib(result, style="poster")
        figure.savefig(
            out_dir / f"{plot.name}_poster.png",
            dpi=600,
            bbox_inches="tight",
            facecolor="white",
            metadata=_PNG_METADATA,
        )
    plt.close(figure)
