"""Shared constants for the segmentation benchmark."""

# ---------------------------------------------------------------------------
# Cell-type annotations
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Segmentation method metadata
# ---------------------------------------------------------------------------

METHOD_COLORS: dict[str, str] = {
    "cellpose":         "#4C72B0",
    "stardist":         "#8172B2",
    "mesmer":           "#D62728",
    "10x_ranger":       "#2E8B57",
    "10x_native":       "#55A868",
    "voronoi":          "#17BECF",
    "voronoi_stardist": "#9467BD",
    "voronoi_mesmer":   "#BCBD22",
    "voronoi_10x_ranger": "#3CB371",
    "baysor":           "#DD8452",
    "baysor_prior":     "#937860",
    "baysor_prior_c08": "#A0522D",
    "baysor_prior_c10":          "#CD853F",
    "baysor_stardist_prior_c10": "#8B4513",
    "baysor_mesmer_prior_c10":      "#DAA520",
    "baysor_10x_ranger_prior_c10": "#006400",
    "bidcell":                      "#E377C2",
    "segger":           "#C49C94",
}

METHOD_LABELS: dict[str, str] = {
    "cellpose":         "CellPose",
    "stardist":         "StarDist",
    "mesmer":             "Mesmer",
    "10x_ranger":         "10x Ranger",
    "10x_native":         "10x native",
    "voronoi":            "Voronoi (CP)",
    "voronoi_stardist":   "Voronoi (SD)",
    "voronoi_mesmer":     "Voronoi (M)",
    "voronoi_10x_ranger": "Voronoi (10x)",
    "baysor":           "Baysor",
    "baysor_prior":     "Baysor (CP prior 0.2)",
    "baysor_prior_c08": "Baysor (CP prior 0.8)",
    "baysor_prior_c10":          "Baysor (CP prior 1.0)",
    "baysor_stardist_prior_c10": "Baysor (SD prior 1.0)",
    "baysor_mesmer_prior_c10":      "Baysor (M prior 1.0)",
    "baysor_10x_ranger_prior_c10": "Baysor (10x prior 1.0)",
    "bidcell":                      "BIDCell",
    "segger":           "Segger",
}

METHOD_FAMILIES: dict[str, str] = {
    "cellpose":         "Nuclear",
    "stardist":         "Nuclear",
    "mesmer":             "Nuclear",
    "10x_ranger":         "Nuclear",
    "10x_native":         "Reference",
    "voronoi":            "Voronoi",
    "voronoi_stardist":   "Voronoi",
    "voronoi_mesmer":     "Voronoi",
    "voronoi_10x_ranger": "Voronoi",
    "baysor":           "Transcript-density",
    "baysor_prior":     "Transcript-density",
    "baysor_prior_c08": "Transcript-density",
    "baysor_prior_c10":          "Transcript-density",
    "baysor_stardist_prior_c10": "Transcript-density",
    "baysor_mesmer_prior_c10":      "Transcript-density",
    "baysor_10x_ranger_prior_c10": "Transcript-density",
    "bidcell":                      "Multimodal",
    "segger":           "Multimodal",
}

NUCLEAR_ONLY: set[str] = {"cellpose", "stardist", "mesmer", "10x_ranger"}

# Methods included in main analysis figures (past the recovery section).
MAIN_METHODS: list[str] = [
    "10x_native",
    "voronoi", "voronoi_stardist", "voronoi_mesmer", "voronoi_10x_ranger",
    "baysor", "baysor_prior",
    "baysor_prior_c08", "baysor_prior_c10",
    "baysor_stardist_prior_c10", "baysor_mesmer_prior_c10",
    "baysor_10x_ranger_prior_c10",
]

# Pairwise comparisons shown in multi-panel figures (2×2+ grid order).
COMPARISON_ORDER: list[tuple[str, str]] = [
    ("voronoi",          "Voronoi (CP)"),
    ("voronoi_stardist", "Voronoi (SD)"),
    ("voronoi_mesmer",   "Voronoi (M)"),
    ("baysor",           "Baysor"),
    ("baysor_prior",     "Baysor (CP prior 0.2)"),
    ("baysor_prior_c08", "Baysor (CP prior 0.8)"),
    ("baysor_prior_c10",          "Baysor (CP prior 1.0)"),
    ("baysor_stardist_prior_c10", "Baysor (SD prior 1.0)"),
    ("baysor_mesmer_prior_c10",      "Baysor (M prior 1.0)"),
    ("voronoi_10x_ranger",           "Voronoi (10x)"),
    ("baysor_10x_ranger_prior_c10",  "Baysor (10x prior 1.0)"),
]

# ---------------------------------------------------------------------------
# Negative marker pairs for biology-grounded segmentation quality scoring
# ---------------------------------------------------------------------------

# Tier 1: high-confidence lineage exclusions - these cell types arise from
# different developmental lineages and never co-express in normal breast tissue.
NEGATIVE_PAIRS_TIER1: list[tuple[str, str, str, str]] = [
    ("CD3E",   "GATA3",  "T cell",       "Luminal epithelial"),
    ("CD3E",   "ESR1",   "T cell",       "Luminal epithelial"),
    ("TRAC",   "KRT14",  "T cell",       "Myoepithelial"),
    ("MS4A1",  "GATA3",  "B cell",       "Luminal epithelial"),
    ("MS4A1",  "KRT14",  "B cell",       "Myoepithelial"),
    ("PECAM1", "ESR1",   "Endothelial",  "Luminal epithelial"),
    ("VWF",    "CD3E",   "Endothelial",  "T cell"),
    ("CD14",   "ESR1",   "Macrophage",   "Luminal epithelial"),
    ("LYZ",    "GATA3",  "Macrophage",   "Luminal epithelial"),
    ("ADIPOQ", "CD3E",   "Adipocyte",    "T cell"),
    ("ADIPOQ", "GATA3",  "Adipocyte",    "Luminal epithelial"),
]

# Tier 2: moderate-confidence exclusions (generally exclusive but with rare
# biological exceptions; ACTA2+epithelial excluded due to EMT possibility).
NEGATIVE_PAIRS_TIER2: list[tuple[str, str, str, str]] = [
    ("LUM",  "CD3E",  "CAF",         "T cell"),
    ("MZB1", "GATA3", "Plasma cell", "Luminal epithelial"),
]
