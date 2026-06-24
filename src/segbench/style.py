"""Shared figure style for all segmentation-benchmark plots.

Import and call ``apply_style()`` at the top of every figure script instead of
calling ``sns.set_theme()`` directly.  This ensures every figure uses the same
fonts, sizes, and background regardless of which script produced it.
"""

from __future__ import annotations

import matplotlib as mpl
import seaborn as sns


def apply_style(scatter: bool = False) -> None:
    """Set a consistent matplotlib/seaborn style for all figures.

    Parameters
    ----------
    scatter:
        Pass ``True`` for scatter/spatial figures (UMAP, tissue maps) where a
        grid background is distracting.  All other figures use ``whitegrid``.
    """
    style = "white" if scatter else "whitegrid"
    sns.set_theme(style=style, context="talk", font_scale=1.3)
    mpl.rcParams.update({
        "axes.titlesize":       20,
        "axes.labelsize":       18,
        "axes.titlepad":        12,
        "xtick.labelsize":      16,
        "ytick.labelsize":      16,
        "legend.fontsize":      16,
        "legend.title_fontsize": 16,
        "figure.titlesize":     22,
        "font.size":            16,
    })
