"""The paper's dataset tables, as the figures and tables that flag their shaded rows quote them.

Capannelli et al., Supplementary notes, Supplementary Fig. 10 (6D2S, p. 10) and
Supplementary Fig. 11 (2x2, p. 11) list the released records and shade some of them. A leaf
that imports nothing, so any layer can import it. Jobs import it and pass the text to their
steps as a kwarg.
"""

from __future__ import annotations

# The text a record those tables shade carries, the same in every figure and table that
# flags it. Its rule: Supplementary notes, p. 9, the paragraph before Supplementary Fig. 10
# (more than 5% of a dataset's points filtered; unsuitable for the Fourier analysis of
# main-text Fig. 8).
FLAG_TEXT = (
    "shaded in the paper's dataset table (high filtered share; the paper excluded >5% "
    "filtered datasets from its Fourier analysis)"
)
