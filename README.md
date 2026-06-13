# Segmentation Benchmarking on Xenium Spatial Transcriptomics Data

## Question

Do cell segmentation methods developed for multiplexed imaging (Mesmer, CellPose)
transfer well to imaging-based spatial transcriptomics (10x Xenium) data, and does the
choice of segmentation method meaningfully change downstream cell type calls? Where
methods disagree, is that disagreement spatially structured (e.g. concentrated at
tissue boundaries or specific niches) or essentially random noise?

This is Project 1 of a two-part series bridging imaging-based spatial biology into
sequencing-based spatial bioinformatics. Project 2 will use this project's output for
scRNA-seq label-transfer reliability analysis.

## Dataset

**Xenium FFPE Human Breast (Custom Add-on Panel)** — Janesick et al. 2023,
*"High resolution mapping of the breast cancer tumor microenvironment using integrated
single cell, spatial and in situ analysis of FFPE tissue"*, *Nature Communications*.

- 10x dataset page: [Xenium FFPE Human Breast with Custom Add-on Panel](https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-with-custom-add-on-panel-1-standard)
- ~577,000 cells, ~78M transcripts, registered post-Xenium H&E image
- Matched scRNA-seq + Visium from the same tissue blocks: GEO [GSE243275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243275)
  (this is what Project 2 will use as the label-transfer reference)

The full bundle (morphology image + transcripts) is tens of GB, so segmentation runs
operate on a cropped ROI (~2-4mm²) selected to contain a mix of tumor, stroma, and
immune-infiltrated regions. See [`docs/dataset.md`](docs/dataset.md) for download
details and ROI selection notes (filled in during data acquisition).

## Methods compared

| Method | Input | Notes |
|---|---|---|
| **10x native** | provided | Xenium's own nucleus/cell boundary segmentation, available as `cell_feature_matrix.h5ad` for the ROI — not yet incorporated into the cross-method comparison below |
| **CellPose** | DAPI (full 2mm x 2mm ROI) | CellPose 3.x classical `nuclei` U-Net model, CPU. CellPose 4.x dropped the lightweight U-Net models in favor of SAM-based foundation models, which are CPU-prohibitive |
| **Mesmer** (DeepCell) | DAPI (ROI crop) | not run — blocked on deepcell.org's account system (signup, login, and password reset all fail as of June 2026); wrapper and conda env are ready, see `docs/dataset.md` |
| **Baysor** | transcripts (centered 1mm x 1mm sub-region of the ROI) | transcript-density-based segmentation, run via Julia 1.10; the full 2mm x 2mm ROI is CPU-impractical for this method (see `docs/dataset.md`) |

For each method: per-cell transcript aggregation → AnnData (cells × genes) → compare
cell counts/density, size distributions, per-cell transcript counts, expression
correlation, and Leiden-based cell type calls (agreement via ARI / confusion matrices).
Finally, map where methods disagree and test whether that disagreement is spatially
structured.

## Results

Results cover **CellPose vs. Baysor** — Mesmer could not be run (see Status)
and this two-method comparison is the deliverable for Project 1.
CellPose ran on the full 2mm x 2mm ROI; Baysor ran on the centered 1mm x 1mm
sub-region, so raw cell counts aren't directly comparable —
[`cell_counts_and_sizes.png`](results/figures/cell_counts_and_sizes.png) reports
density (cells/mm²) instead.

| | CellPose | Baysor |
|---|---|---|
| Cells | 20,166 (4mm² ROI) | 4,514 (1mm² sub-region) |
| Density | ~5,040 cells/mm² | ~4,510 cells/mm² |
| Median size | 653 px² nucleus area (~29.5 µm²) | 50 transcripts/cell |
| Median transcripts/cell | 49 | 50 |

Despite very different segmentation strategies (nucleus-pixel masks vs.
transcript-density clustering), cell density and median transcripts-per-cell
are remarkably close.

**Matching**: 2,101 mutual-nearest-neighbor pairs (≤10 µm centroid distance)
out of 20,166 CellPose / 4,514 Baysor cells — limited by Baysor's smaller
footprint, so every matched-pair metric below is implicitly scoped to that
1mm² sub-region.

**Expression agreement**
([`expression_correlation.png`](results/figures/expression_correlation.png)):
median per-pair Pearson correlation = **0.74** across shared genes — fairly
strong per-cell expression agreement given the two methods use completely
different inputs.

**Cell-type agreement**: independent Leiden clustering on each method (13
CellPose clusters, 17 Baysor clusters) was Hungarian-aligned onto a shared
cluster vocabulary before direct comparison
([`cell_type_confusion.png`](results/figures/cell_type_confusion.png)).
**46.7%** of matched pairs were assigned to different cell-type clusters by
the two methods.

**Spatial structure of disagreement**
([`disagreement_spatial_map.png`](results/figures/disagreement_spatial_map.png)):
Moran's I = **0.0655** (permutation test, p = 0.001) — statistically
significant but weak spatial autocorrelation. Cross-method disagreement is
mostly scattered throughout the tissue rather than concentrated in specific
regions, with only a slight tendency to cluster.

## Repo layout

```
segmentation-benchmark/
├── environment.yml          # conda env (CellPose, Scanpy, Squidpy, SpatialData, ...)
├── data/
│   ├── raw/                 # downloaded Xenium bundle (gitignored)
│   └── processed/           # cropped ROI + derived files (gitignored)
├── notebooks/                # exploratory analysis
├── src/segbench/
│   ├── io.py                 # load Xenium bundle (spatialdata-io), ROI cropping
│   ├── segmentation/          # per-method wrappers (cellpose, mesmer/docker, baysor/julia)
│   ├── quantify.py            # transcript aggregation -> per-cell AnnData
│   ├── compare.py             # cross-method comparison metrics
│   └── spatial.py             # spatial structure of disagreement
├── scripts/                   # CLI entry points (download, run_mesmer.py, run_baysor.sh)
├── results/{figures,tables}/  # outputs (gitignored, regeneratable)
└── tests/
```

## Environment setup

This project uses three toolchains: a conda env for CellPose + the
Scanpy/Squidpy/SpatialData stack, a separate conda env for Mesmer (DeepCell's
TensorFlow pin is incompatible with the main env), and Julia for Baysor.

### 1. Python (conda, for CellPose + Scanpy/Squidpy/SpatialData stack)

```bash
conda env create -f environment.yml
conda activate segbench
```

### 2. Mesmer (DeepCell, separate conda env)

```bash
conda create -n mesmer python=3.10
conda run -n mesmer pip install deepcell
```

Requires a `DEEPCELL_ACCESS_TOKEN` (free account at
[users.deepcell.org](https://users.deepcell.org)) to download `Mesmer`'s
pretrained weights on first use. As of June 2026, deepcell.org's account
system is non-functional (signup returns a server error, login and password
reset both fail) — the env and wrapper (`segbench.segmentation.mesmer_run`,
`scripts/run_mesmer.py`) are ready to use once that's resolved.

### 3. Julia + Baysor (for transcript-based segmentation)

Baysor's `main`/`cpp` branch is a C++ rewrite with no `Project.toml`; install the
last Julia-package release (`v0.7.1`) instead:

```bash
juliaup add 1.10
julia +1.10 -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git", rev="v0.7.1")); Pkg.build("Baysor")'
```

See [`scripts/run_baysor.sh`](scripts/run_baysor.sh) for the invocation.

## Status

- [x] Project scaffold + environments
- [x] Data acquisition (`scripts/download_data.sh`, see `docs/dataset.md`)
- [x] ROI selection
- [x] Segmentation: CellPose, Baysor
- [x] Quantification + cross-method comparison (CellPose vs. Baysor)
- [x] Spatial analysis of disagreement
- [x] Figures + write-up
- [ ] Segmentation: Mesmer — blocked externally on deepcell.org's account
      system (signup/login/password-reset all fail); wrapper and conda env
      are ready, revisit if/when access is restored
- [ ] Stretch: incorporate 10x native segmentation as a third reference
