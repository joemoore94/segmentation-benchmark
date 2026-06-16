# Segmentation Benchmarking on Xenium Spatial Transcriptomics Data

**Question:** Do segmentation methods developed for multiplexed imaging (CellPose, StarDist) transfer well to imaging-based spatial transcriptomics (10x Xenium), and does method choice meaningfully change downstream cell-type calls?

## Key findings

| Comparison | ARI | Disagreement rate | Moran's I |
|---|---|---|---|
| CellPose vs. StarDist | 0.764 | 13.8% | 0.066 |
| CellPose vs. 10x native | 0.547 | 30.8% | 0.176 |
| CellPose vs. Baysor | 0.415 | 48.8% | 0.034 |
| CellPose vs. Baysor (CellPose prior) | 0.407 | 49.7% | 0.051 |

Algorithm choice within a modality (CellPose vs. StarDist, same DAPI image) matters far less than segmentation modality. Disagreement against the 10x reference is the most spatially structured (Moran's I 0.176); Baysor vs. CellPose disagreement is higher but nearly spatially uniform. Cross-modality disagreement preferentially hits phenotypically rare cells: CellPose-Baysor disagreeing pairs sit ~2 log-units lower in Mellon density (p = 2.3e-76); same-modality comparisons show no comparable effect.

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
| **Mesmer** (DeepCell) | DAPI | not run; deepcell.org account system non-functional as of June 2026; env and wrapper ready |
| **Baysor** | transcripts (2mm x 2mm, 4 tiles) | transcript-density EM, Julia 1.10 |
| **Baysor (CellPose prior)** | transcripts + CellPose nuclei (4 tiles) | Baysor with `--prior-segmentation` from CellPose masks, `prior_segmentation_confidence=0.2` |

Per-cell transcript aggregation → AnnData → cell counts, transcript capture, expression correlation, Leiden clustering, and spatial structure of disagreement (Moran's I, Mellon density).

## Results

### Cell counts and transcript capture

| | CellPose | Baysor | 10x native | StarDist | Baysor (prior) |
|---|---|---|---|---|---|
| Cells | 20,166 | 18,321 | 23,629 | 24,745 | 19,061 |
| Median transcripts/cell | 49 | 53 | 124 | 45 | 53 |
| Transcript capture | 35.4% | 98.6% | 99.0% | 40.8% | 98.7% |

![Cell counts, transcripts/cell, and nucleus area by method](results/figures/cell_counts_and_sizes.png)

Nuclear-only methods (CellPose, StarDist) capture 35-41% of transcripts because cytoplasmic transcripts fall outside nucleus masks; whole-cell and transcript-based methods capture ~99%. Median nucleus area is nearly identical across CellPose, StarDist, and 10x native (27-30 µm²).

### Clustering structure

![PCA and UMAP embeddings colored by Leiden cluster, per method](results/figures/pca_umap_clusters.png)

All five methods produce well-separated UMAP clusters (12-24 Leiden clusters). Baysor variants produce more clusters (21-24) than nuclear methods (12-15), consistent with their higher per-cell transcript counts resolving finer expression differences.

### Pairwise comparisons (all vs. CellPose)

![Per-cell-pair expression correlation](results/figures/expression_correlation.png)
![Cell-type cluster correspondence](results/figures/cell_type_confusion.png)
![Disagreement mapped spatially](results/figures/disagreement_spatial_map.png)

**CellPose vs. StarDist** (same DAPI input, different algorithm): 19,460 matched pairs, median expression correlation 0.96, ARI 0.764, 13.8% disagreement. The two algorithms trace the same nuclei closely.

**CellPose vs. 10x native** (nuclear vs. whole-cell): 18,966 matched pairs, correlation 0.82, ARI 0.547, 30.8% disagreement. Disagreement is the most spatially structured of any comparison (Moran's I 0.176), concentrated in specific tissue regions where nucleus segmentation is harder.

**CellPose vs. Baysor** (nuclear vs. transcript-density): 8,947 matched pairs, correlation 0.73, ARI 0.415, 48.8% disagreement, Moran's I 0.034. Cell counts are similar (20,166 vs. 18,321) but nearly half of matched cells land in different clusters, and that disagreement is spatially uniform.

**CellPose vs. Baysor (prior)**: 9,572 matched pairs, correlation 0.74, ARI 0.407, 49.7% disagreement, Moran's I 0.051. The CellPose-nucleus prior produces 7% more matched pairs than the non-prior run but leaves cell-type agreement essentially unchanged.

### Phenotypic density vs. disagreement (Mellon)

![CellPose phenotypic density vs. disagreement](results/figures/density_vs_disagreement.png)

Each CellPose cell gets a Mellon log-density estimate in PCA space; disagreeing vs. agreeing cells are compared via Mann-Whitney U:

| Comparison | n agree / disagree | Median log-density (agree / disagree) | p |
|---|---|---|---|
| CellPose vs. Baysor | 4,585 / 4,362 | -15.00 / -17.01 | 2.3e-76 |
| CellPose vs. Baysor (prior) | 4,812 / 4,760 | -14.77 / -17.03 | 4.4e-91 |
| CellPose vs. 10x native | 13,121 / 5,845 | -13.28 / -13.54 | 3.5e-11 |
| CellPose vs. StarDist | 16,780 / 2,680 | -13.39 / -13.27 | 6.1e-09 |

Cross-modality disagreement (both Baysor runs) concentrates on phenotypically rare cells (~2 log-unit gap); same-modality or same-input disagreement (10x native, StarDist) shows negligible density differences despite large sample sizes.

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
- [x] Quantification + cross-method comparison
- [x] Spatial disagreement analysis (Moran's I)
- [x] Mellon phenotypic-density analysis
- [x] PCA/UMAP per-method clustering
- [x] Baysor (CellPose prior) hybrid
- [ ] Mesmer: blocked on deepcell.org account system
- [ ] Future: local Moran's I (pockets of disagreement vs. global average); Kompot DE on agree/disagree cell groups
