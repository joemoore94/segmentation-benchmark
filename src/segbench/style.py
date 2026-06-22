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
    sns.set_theme(style=style, context="poster", font_scale=1.1)
    mpl.rcParams.update({
        # axes
        "axes.titlesize":       17,
        "axes.labelsize":       14,
        "axes.titlepad":        10,
        # ticks
        "xtick.labelsize":      13,
        "ytick.labelsize":      13,
        # legend
        "legend.fontsize":      13,
        "legend.title_fontsize": 13,
        # figure title (suptitle)
        "figure.titlesize":     18,
        # general
        "font.size":            14,
    })
