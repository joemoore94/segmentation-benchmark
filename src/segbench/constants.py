"""Shared constants for the segmentation benchmark."""

# Manual annotation based on top Wilcoxon marker genes per cluster.
# Clusters 0, 1, 3, 8 are luminal epithelial subtypes; grouped under
# "Luminal epithelial" for cross-method comparisons. Clusters 2 and 7 are
# both myeloid/macrophage populations with slightly different marker profiles.
CLUSTER_ANNOTATIONS: dict[str, str] = {
    "0":  "Luminal epithelial",   # ESR1, FOXA1, GATA3, ANKRD30A
    "1":  "Luminal epithelial",   # MYBPC1, MUC1, CLIC6, SERPINA3, GATA3, PGR
    "2":  "Macrophages",          # FCER1G, LYZ, FCGR3A, HAVCR2
    "3":  "Luminal epithelial",   # SERPINA3, MUC1, CLIC6, S100A14, FASN, PGR
    "4":  "Myoepithelial",        # ACTA2, KRT14, MYH11, MYLK
    "5":  "T cells",              # CD3E, TRAC, CD96, GZMA, IL7R
    "6":  "B cells",              # MS4A1, BANK1, CD52, SELL
    "7":  "Macrophages",          # CD14, AIF1, LYZ, FCER1G, FGL2
    "8":  "Luminal epithelial",   # TACSTD2, KRT7, GATA3, MUC1, CCND1
    "9":  "CAFs",                 # LUM, SFRP4, FBLN1, THBS2, MMP2
    "10": "Smooth muscle",        # ACTA2, MYH11, MYLK, RGS5, CAV1
    "11": "Endothelial",          # PECAM1, VWF, AQP1, CD93, CLEC14A
    "12": "Plasma cells",         # MZB1, SLAMF7, ITM2C
    "13": "CAFs",                 # POSTN, PDGFRB, CTHRC1, LUM
    "14": "Adipocytes",           # ADIPOQ, PLIN1, PPARG, LPL
}

CELLTYPE_COLORS: dict[str, str] = {
    "Luminal epithelial": "#E377C2",
    "Myoepithelial":      "#8172B2",
    "T cells":            "#2CA02C",
    "B cells":            "#17BECF",
    "Macrophages":        "#D62728",
    "Plasma cells":       "#FF7F0E",
    "CAFs":               "#7F7F7F",
    "Smooth muscle":      "#BCBD22",
    "Endothelial":        "#1F77B4",
    "Adipocytes":         "#8C564B",
}
