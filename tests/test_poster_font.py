"""The poster target refuses to render in a typeface it was not asked for.

Oracle: matplotlib's own font resolution, asked with `fallback_to_default=False`.

Roboto is not shipped with matplotlib, is not a Python dependency, and cannot be pip
installed. It is an environment requirement of the poster target, and matplotlib's response
to a missing family is to log a warning and substitute its default: the figure still renders
and still looks plausible, in different type. That is the silent fallback AGENTS.md section 3
exists to refuse.

Two traps, both load-bearing. `findfont` is memoised, so a test that asserted on the warning
would be order-dependent: an earlier resolution of the same family suppresses the second
warning and the assertion could not fail. These assert on the raise instead. And matplotlib
reports the miss through `logging`, not `warnings`, so `pytest.warns` would never fire.

Tests needing the real face SKIP where it is absent, the way the `r` tier does. Failing there
would turn the default test run red on every machine that does not happen to have it.
"""

from __future__ import annotations

import pytest

from quebra.plots import theme

pytestmark = pytest.mark.unit


def _poster_font_or_skip() -> None:
    """Skip where the poster's required face is absent. See the module docstring."""
    try:
        theme.check_fonts("poster")
    except theme.MissingFontError as exc:
        pytest.skip(str(exc))


def _clear_font_cache() -> None:
    """Clear matplotlib's memoised font lookup.

    The cache lives on the FontManager INSTANCE, not on the module. An earlier version of
    this helper looked for a module attribute, found nothing, and silently did nothing, while
    the docstring claimed the cache was cleared. No `getattr` default here: if matplotlib
    moves it again, this should raise rather than quietly stop working.
    """
    from matplotlib import font_manager

    font_manager.fontManager._findfont_cached.cache_clear()


def test_the_poster_stack_is_a_list_with_a_real_fallback():
    """A single-element stack would make the check unsatisfiable, not strict."""
    assert isinstance(theme.POSTER_FONT_STACK, list)
    assert len(theme.POSTER_FONT_STACK) >= 2
    assert theme.RCPARAMS["poster"]["font.family"] == theme.POSTER_FONT_STACK
    assert theme.REQUIRED_FIRST_FACE["poster"] == theme.POSTER_FONT_STACK[0]


def test_the_configured_poster_face_resolves_to_itself():
    """R9.2's acceptance: the face the target asks for is the face it gets.

    Structure tests cannot see this. Only asking matplotlib to resolve the name, refusing its
    default, says whether a poster rendered here is in the type this project chose.
    """
    _poster_font_or_skip()
    _clear_font_cache()
    from matplotlib import font_manager

    resolved = font_manager.findfont(
        theme.POSTER_FONT_STACK[0], fallback_to_default=False
    )
    assert theme.POSTER_FONT_STACK[0].replace(" ", "").lower() in (
        resolved.replace(" ", "").lower()
    ), f"{theme.POSTER_FONT_STACK[0]!r} resolved to {resolved!r}"
    _clear_font_cache()


def test_a_missing_poster_face_raises_instead_of_substituting(monkeypatch):
    """Point the target at a face no machine has, and it must refuse."""
    _clear_font_cache()
    monkeypatch.setitem(theme.REQUIRED_FIRST_FACE, "poster", "NoSuchFace ZZZ")
    with pytest.raises(theme.MissingFontError) as excinfo:
        theme.check_fonts("poster")
    message = str(excinfo.value)
    assert "NoSuchFace ZZZ" in message, "the error must name the face that is missing"
    assert "plots/theme.py" in message, "and where to change it"
    _clear_font_cache()


def test_entering_the_poster_style_context_performs_the_check(monkeypatch):
    """The wiring, not just the function.

    `check_fonts` is only a guard if something calls it. Deleting its call from
    `style_context` left the whole suite green, because every other test here drives the
    function directly.
    """
    _clear_font_cache()
    monkeypatch.setitem(theme.REQUIRED_FIRST_FACE, "poster", "NoSuchFace ZZZ")
    with pytest.raises(theme.MissingFontError):
        with theme.style_context("poster"):
            pass
    _clear_font_cache()


def test_a_target_with_no_required_face_is_not_checked(monkeypatch):
    """`static` and `academic` name no required face, so the guard must not fire for them."""
    _clear_font_cache()
    monkeypatch.setitem(theme.REQUIRED_FIRST_FACE, "poster", "NoSuchFace ZZZ")
    theme.check_fonts("static")
    theme.check_fonts("academic")
    with pytest.raises(theme.MissingFontError):
        theme.check_fonts("poster")
    _clear_font_cache()
