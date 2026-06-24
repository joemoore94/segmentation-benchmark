"""Per-cell local Moran's I (LISA) for each pairwise disagreement table.

Reads the ``disagreement_table_10x_*.csv`` files written by
``run_comparison.py`` and appends two columns: ``local_morans_i`` (the LISA
statistic; positive = surrounded by like neighbours, negative = spatial
outlier) and ``lisa_cluster`` (HH / LL / HL / LH quadrant label). Writes
one ``local_morans_10x_*.csv`` per comparison to ``results/tables/``.

Usage::

    conda run -n segbench python scripts/run_local_morans.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from segbench.constants import METHOD_LABELS
from segbench.spatial import local_morans_i, local_morans_i_cluster

TABLES_DIR = Path("results/tables")

_METHODS = [
    "cellpose", "stardist", "mesmer",
    "voronoi", "voronoi_stardist", "voronoi_mesmer",
    "baysor", "baysor_prior", "baysor_prior_c08", "baysor_prior_c10", "baysor_stardist_prior_c10", "baysor_mesmer_prior_c10", "bidcell", "segger",
]
COMPARISONS = {
    f"10x native vs. {METHOD_LABELS[k]}": f"disagreement_table_10x_{k}.csv"
    for k in _METHODS
}


def main() -> None:
    for label, fname in COMPARISONS.items():
        if not (TABLES_DIR / fname).exists():
            print(f"\n=== {label}: skipped (file not found) ===")
            continue
        df = pd.read_csv(TABLES_DIR / fname)
        coords = df[["centroid_x", "centroid_y"]].to_numpy()
        values = df["disagree"].to_numpy(dtype=float)

        df["local_morans_i"] = local_morans_i(coords, values)
        df["lisa_cluster"] = local_morans_i_cluster(coords, values)

        out_name = fname.replace("disagreement_table_", "local_morans_")
        df.to_csv(TABLES_DIR / out_name, index=False)

        print(f"\n=== {label} ===")
        print(df["lisa_cluster"].value_counts().to_string())
        hh = (df["lisa_cluster"] == "HH").sum()
        ll = (df["lisa_cluster"] == "LL").sum()
        print(f"HH hotspots (disagreement clusters): {hh} ({hh / len(df):.1%})")
        print(f"LL coldspots (agreement clusters):   {ll} ({ll / len(df):.1%})")


if __name__ == "__main__":
    main()
