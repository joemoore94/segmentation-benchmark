"""Compute many-to-one cluster-matching accuracy for each comparison.

Hungarian matching enforces 1-to-1 cluster assignment. When one method
splits a cell type into more sub-clusters than the other, some of those
sub-clusters are forced to match unrelated 10x clusters, inflating the
disagreement rate. Many-to-one relaxes this: each comparison cluster is
independently assigned to its argmax 10x cluster (several comp clusters may
map to the same 10x cluster). The resulting accuracy is an upper bound that
isolates genuine cell-type misclassification from over-clustering artefacts.

  many-to-one accuracy = sum_j(max_i conf[i,j]) / sum(conf)

where rows i = 10x clusters, columns j = comparison clusters.

Reads:  results/tables/cell_type_confusion_10x_{method}.csv
Writes: results/tables/many_to_one_accuracy.csv
"""

from pathlib import Path

import pandas as pd

TABLES_DIR = Path("results/tables")

COMPARISONS = [
    ("cellpose",       "CellPose"),
    ("stardist",       "StarDist"),
    ("mesmer",         "Mesmer"),
    ("voronoi",        "Voronoi (CP)"),
    ("voronoi_mesmer", "Voronoi (M)"),
    ("baysor",         "Baysor"),
]

rows = []
for method, label in COMPARISONS:
    path = TABLES_DIR / f"cell_type_confusion_10x_{method}.csv"
    if not path.exists():
        print(f"  missing {path.name}, skipping")
        continue
    conf = pd.read_csv(path, index_col=0)
    total = conf.values.sum()
    m2o_correct = conf.max(axis=0).sum()
    m2o_acc = m2o_correct / total
    hungarian_acc = 1.0 - pd.read_csv(
        TABLES_DIR / f"disagreement_table_10x_{method}.csv"
    )["disagree"].mean()
    rows.append({
        "comparison": label,
        "total_pairs": int(total),
        "hungarian_accuracy": round(hungarian_acc, 4),
        "many_to_one_accuracy": round(m2o_acc, 4),
        "gap": round(m2o_acc - hungarian_acc, 4),
    })
    print(f"{label}: Hungarian {hungarian_acc:.3f}  M2O {m2o_acc:.3f}  gap {m2o_acc - hungarian_acc:+.3f}")

result = pd.DataFrame(rows)
result.to_csv(TABLES_DIR / "many_to_one_accuracy.csv", index=False)
print("\nSaved many_to_one_accuracy.csv")
print(result.to_string(index=False))
