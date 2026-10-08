from pathlib import Path

"""Published data used to pin transcriptions against their sources."""

# Tier 4 fixtures, written by `rscripts/reference_values.R` and read by
# `tests/test_r_cross_implementation.py`. Committed so the suite never requires R.
# In jobs/reference/, NOT here: `tests/` is gitignored (see .gitignore), so a fixture
# living beside this file would not survive a clone and the suite would quietly start
# requiring R to regenerate it. rscripts/ is tracked, and keeping the CSVs next to the
# script that writes them also keeps their provenance obvious. `reference/` is the one
# tracked home for external validation data - published tables and R outputs alike -
# because a FIGURE now consumes them, not only the suite.
_REFERENCE_DIR = Path(__file__).resolve().parents[2] / "jobs" / "reference"
R_REFERENCE_INPUTS = _REFERENCE_DIR / "r_reference_inputs.csv"
R_REFERENCE_VALUES = _REFERENCE_DIR / "r_reference_values.csv"
# Written by `rscripts/survival_reference.R`, read by `tests/test_survival_r_reference.py`.
R_SURVIVAL_INPUTS = _REFERENCE_DIR / "r_survival_inputs.csv"
R_SURVIVAL_VALUES = _REFERENCE_DIR / "r_survival_values.csv"
R_TURNBULL_INPUTS = _REFERENCE_DIR / "r_turnbull_inputs.csv"
