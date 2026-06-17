# Segmentation Benchmarking on Xenium Spatial Transcriptomics Data

**Question:** Do segmentation methods developed for multiplexed imaging (CellPose, StarDist) transfer well to imaging-based spatial transcriptomics (10x Xenium), and does method choice meaningfully change downstream cell-type calls?

## Key findings

| Comparison | Matched pairs | Median corr | ARI | Disagreement rate | Moran's I |
|---|---|---|---|---|---|
| 10x native vs. CellPose | 18,966 | 0.822 | 0.547 | 30.8% | 0.178 |
| 10x native vs. StarDist | 21,429 | 0.826 | 0.545 | 33.5% | 0.215 |
| 10x native vs. Voronoi | 18,966 | 0.932 | 0.630 | 21.9% | 0.076 |
| 10x native vs. Baysor (prior) | 11,454 | 0.798 | 0.318 | 51.9% | 0.036 |
| 10x native vs. Baysor | 10,953 | 0.786 | 0.305 | 51.7% | 0.033 |

Three method families emerge. Nuclear methods (CellPose, StarDist): ARI ~0.55, ~31-34% disagreement, Moran's I 0.18-0.22 — spatially structured disagreement concentrated in tissue regions where nuclear detection is harder. Voronoi (CellPose nuclei, nearest-centroid transcript assignment): ARI 0.630, 21.9% disagreement, Moran's I 0.076, median expression correlation 0.932 — the best-performing non-reference method; its residual disagreement shows no phenotypic density effect (p=0.19), meaning remaining errors are purely geometric rather than cell-state-driven. Transcript-density methods (Baysor variants): ARI ~0.31, ~52% disagreement, Moran's I 0.033-0.036 — near-random spatial disagreement. Adding a CellPose-nucleus prior to Baysor marginally improves ARI (0.305 → 0.318) without changing the disagreement pattern, confirming that nucleus detection is not the bottleneck for Baysor. The Voronoi result also demonstrates that the nuclear-method gap (ARI 0.55 → 0.63) is explained almost entirely by the absence of cytoplasmic transcripts in nuclear-only segmentation.

This is Project 1 of a portfolio bridging imaging-based spatial biology into sequencing-based bioinformatics. Project 2 ([label-transfer-benchmark](https://github.com/joemoore94/label-transfer-benchmark)) uses this project's segmented cells to evaluate scRNA-seq label-transfer reliability.

## Dataset

**Xenium FFPE Human Breast (Custom Add-on Panel)**, Janesick et al. 2023, *Nature Communications* ([dataset page](https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-with-custom-add-on-panel-1-standard)). ~577,000 cells, ~78M transcripts, registered post-Xenium H&E. Matched scRNA-seq + Visium from the same tissue blocks: GEO [GSE243275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243275).

Segmentation runs on a ~2mm x 2mm ROI with a mix of tumor, stroma, and immune-infiltrated regions. See [`docs/dataset.md`](docs/dataset.md) for download and ROI details.

## Methods

| Method | Input | Notes |
|---|---|---|
| **10x native** | provided | Xenium Ranger's own segmentation, reshaped by `scripts/build_10x_adata.py` |
| **CellPose** | DAPI (2mm x 2mm ROI) | CellPose 3.x `nuclei` model, CPU |
| **StarDist** | DAPI (2mm x 2mm ROI) | `2D_versatile_fluo` model, separate `stardist` conda env |
| **Voronoi** | CellPose nuclear centroids | Nearest-centroid transcript assignment via scipy cKDTree; 100% capture, no additional model |
| **Mesmer** (DeepCell) | DAPI | not run; deepcell.org account system non-functional as of June 2026; env and wrapper ready |
| **Baysor** | transcripts (2mm x 2mm, 4 tiles) | transcript-density EM, Julia 1.10 |
| **Baysor (CellPose prior)** | transcripts + CellPose nuclei (4 tiles) | Baysor with `--prior-segmentation` from CellPose masks, `prior_segmentation_confidence=0.2` |

Per-cell transcript aggregation → AnnData → cell counts, transcript capture, expression correlation, Leiden clustering, and spatial structure of disagreement (Moran's I, Mellon density). Independent Leiden runs assign arbitrary cluster IDs, so cluster labels are aligned across methods using the Hungarian algorithm (linear sum assignment on the confusion matrix) before computing disagreement rate and ARI.

## Results

### Cell counts and transcript capture

| | CellPose | StarDist | Voronoi | Baysor | 10x native | Baysor (prior) |
|---|---|---|---|---|---|---|
| Cells | 20,166 | 24,745 | 20,166 | 18,321 | 23,629 | 19,061 |
| Median transcripts/cell | 49 | 45 | 168 | 53 | 124 | 53 |
| Transcript capture | 35.4% | 40.8% | 100.0% | 98.6% | 99.0% | 98.7% |

![Cell counts, transcripts/cell, and nucleus area by method](results/figures/cell_counts_and_sizes.png)

Nuclear-only methods (CellPose, StarDist) capture 35-41% of transcripts because cytoplasmic transcripts fall outside nucleus masks; whole-cell and transcript-based methods capture ~99%. Median nucleus area is nearly identical across CellPose, StarDist, and 10x native (27-30 µm²).

### Clustering structure

![PCA and UMAP embeddings colored by Leiden cluster, per method](results/figures/pca_umap_clusters.png)

All five methods produce well-separated UMAP clusters (12-24 Leiden clusters). Baysor variants produce more clusters (21-24) than nuclear methods (12-15), consistent with their higher per-cell transcript counts resolving finer expression differences.

### Pairwise comparisons (all vs. 10x native)

![Per-cell-pair expression correlation](results/figures/expression_correlation.png)
![Cell-type cluster correspondence](results/figures/cell_type_confusion.png)
![Disagreement mapped spatially](results/figures/disagreement_spatial_map.png)

All four comparisons use 10x native (Xenium Ranger's own segmentation) as the reference anchor.

**10x native vs. CellPose** (whole-cell vs. nuclear): 18,966 matched pairs, median expression correlation 0.822, ARI 0.547, 30.8% disagreement, Moran's I 0.178.

**10x native vs. StarDist** (whole-cell vs. nuclear): 21,429 matched pairs, correlation 0.826, ARI 0.545, 33.5% disagreement, Moran's I 0.215.

**10x native vs. Voronoi** (whole-cell vs. nearest-centroid expansion from CellPose nuclei): 18,966 matched pairs, correlation 0.932, ARI 0.630, 21.9% disagreement, Moran's I 0.076. Voronoi assigns all transcripts to the nearest CellPose nuclear centroid, capturing the full cytoplasmic signal with no additional model. The substantially higher ARI and correlation relative to nuclear CellPose (same 20,166 cells, same centroids) directly quantifies the contribution of cytoplasmic transcripts to cell-type identity.

**10x native vs. Baysor**: 10,953 matched pairs, correlation 0.786, ARI 0.305, 51.7% disagreement, Moran's I 0.033. More than half of matched cells land in different clusters, and the pattern is near-random spatially.

**10x native vs. Baysor (prior)**: 11,454 matched pairs, correlation 0.798, ARI 0.318, 51.9% disagreement, Moran's I 0.036. The CellPose-nucleus prior adds ~5% more matched pairs and marginally improves ARI (0.305 → 0.318) but leaves the fundamental disagreement pattern unchanged.

### Phenotypic density vs. disagreement (Mellon)

![10x native phenotypic density vs. disagreement](results/figures/density_vs_disagreement.png)

Each 10x-native cell gets a Mellon log-density estimate in PCA space; disagreeing vs. agreeing cells compared by Mann-Whitney U:

| Comparison | n agree / disagree | Median log-density (agree / disagree) | p |
|---|---|---|---|
| 10x native vs. CellPose | 13,121 / 5,845 | -21.31 / -20.78 | 2.9e-28 |
| 10x native vs. StarDist | 14,254 / 7,175 | -21.87 / -20.63 | 1.1e-90 |
| 10x native vs. Voronoi | 14,805 / 4,161 | -21.05 / -21.35 | 0.191 (n.s.) |
| 10x native vs. Baysor | 5,286 / 5,667 | -22.76 / -22.75 | 0.756 (n.s.) |
| 10x native vs. Baysor (prior) | 5,510 / 5,944 | -22.99 / -22.70 | 0.120 (n.s.) |

Nuclear methods disagree with 10x native on cells in *higher*-density phenotypic regions (p ≪ 0.001), explained by luminal breast epithelial cells (DE top genes: SERPINA3, MUC1, PGR, GATA3, FASN) whose cytoplasmic expression is captured by whole-cell segmentation but missed by nuclear-only methods. Voronoi's disagreement is density-neutral (p=0.19) — its remaining 21.9% error is geometric (Voronoi partition vs. true cell boundary) rather than cell-state-driven. Baysor variants show no density effect, and their disagreement concentrates on macrophages (CD14, MRC1, CD163) with T cells (TRAC, CD3E) robustly identified by all methods.

### Local Moran's I (LISA)

![LISA hotspot/coldspot maps per comparison](results/figures/local_morans_map.png)

HH clusters (local disagreement hotspots) and LL clusters (local agreement coldspots) per comparison:

| Comparison | HH hotspots | LL coldspots |
|---|---|---|
| 10x native vs. CellPose | 21.7% | 30.3% |
| 10x native vs. StarDist | 18.6% | 15.0% |
| 10x native vs. Baysor | 21.4% | 17.5% |
| 10x native vs. Baysor (prior) | 22.1% | 17.6% |

CellPose vs. 10x native has the most agreement coldspots (30.3% LL) — dense regions of tissue where both methods reliably agree — consistent with the global Moran's I finding that this comparison's disagreement is the most spatially concentrated.

### Differential expression: agree vs. disagree cells

![Volcano plots: DE genes in disagree vs. agree groups](results/figures/de_volcano.png)

## Repo layout

```
segmentation-benchmark/
├── environment.yml          # conda env (CellPose, Scanpy, Squidpy, SpatialData, ...)
├── data/
│   ├── raw/                 # downloaded Xenium bundle (gitignored)
│   └── processed/           # cropped ROI + derived files (gitignored)
├── notebooks/
├── src/segbench/
│   ├── io.py                # load Xenium bundle, ROI cropping
│   ├── segmentation/        # per-method wrappers
│   ├── quantify.py          # transcript aggregation -> per-cell AnnData
│   ├── compare.py           # cross-method comparison metrics
│   └── spatial.py           # spatial structure of disagreement
├── scripts/                 # CLI entry points
├── results/{figures,tables}/
└── tests/
```

## Environment setup

This project uses four toolchains: a main conda env for CellPose + Scanpy/Squidpy/SpatialData, a separate env for StarDist (TensorFlow-based), a separate env for Mesmer (incompatible TensorFlow pin), and Julia for Baysor.

### 1. Main env (CellPose + Scanpy stack)

```bash
conda env create -f environment.yml
conda activate segbench
```

### 2. StarDist

```bash
conda create -n stardist python=3.10
conda run -n stardist pip install stardist tensorflow-cpu
```

### 3. Mesmer (DeepCell)

```bash
conda create -n mesmer python=3.10
conda run -n mesmer pip install deepcell
```

Requires a `DEEPCELL_ACCESS_TOKEN` from [users.deepcell.org](https://users.deepcell.org); as of June 2026 that site's account system is non-functional.

### 4. Julia + Baysor

```bash
juliaup add 1.10
julia +1.10 -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git", rev="v0.7.1")); Pkg.build("Baysor")'
```

See [`scripts/run_baysor.sh`](scripts/run_baysor.sh) for the invocation.

## Status

- [x] Project scaffold + environments
- [x] Data acquisition (`scripts/download_data.sh`)
- [x] Segmentation: CellPose, Baysor, StarDist, 10x native
- [x] Quantification + cross-method comparison (10x native anchor)
- [x] Spatial disagreement analysis (global Moran's I)
- [x] PCA/UMAP per-method clustering
- [x] Baysor (CellPose prior) hybrid
- [x] Mellon phenotypic-density analysis (10x native anchor)
- [x] Local Moran's I (LISA) — disagreement hotspot/coldspot maps
- [x] DE: agree vs. disagree cells (Wilcoxon rank-sum)
- [ ] Mesmer: blocked on deepcell.org account system
