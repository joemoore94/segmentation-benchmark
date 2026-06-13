# Segmentation Benchmarking on Xenium Spatial Transcriptomics Data

## Question

Do cell segmentation methods developed for multiplexed imaging (Mesmer, CellPose)
transfer well to imaging-based spatial transcriptomics (10x Xenium) data, and does the
choice of segmentation method meaningfully change downstream cell type calls? Where
methods disagree, is that disagreement spatially structured (e.g. concentrated at
tissue boundaries or specific niches) or essentially random noise?

This is Project 1 of a portfolio bridging imaging-based spatial biology into
sequencing-based spatial bioinformatics and cancer genomics:

- **Project 2** — [label-transfer-benchmark](https://github.com/joemoore94/label-transfer-benchmark):
  uses this project's segmented/quantified cells to evaluate how reliably
  scRNA-seq reference cell-type labels transfer onto spatial data.
- **Project 3** (separate repo, not yet started): a multimodal cancer-progression
  risk model on bulk RNA-seq + DNA methylation (Barrett's esophagus → adenocarcinoma
  progression cohorts). Different domain (bulk cancer genomics rather than spatial
  imaging), but the same patient-level-split / cross-dataset-generalization rigor.

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

Both methods are compared on the **same centered 1mm x 1mm sub-region** of the
ROI — Baysor's full segmentation footprint, with CellPose's full-ROI result
subset to the identical bounds — so cell counts and per-cell transcript counts
below are directly, area-matched comparable (no density normalization needed).

| | CellPose (nuclear, DAPI) | Baysor (transcript density) |
|---|---|---|
| Cells (1mm x 1mm) | 4,459 | 4,514 |
| Median size | 619 px² nucleus area (~28.0 µm²) | 50 transcripts/cell |
| Median transcripts/cell | 48 | 50 |
| Mean transcripts/cell | 59.5 | 168.3 |
| Transcript capture rate | 34.4% | 98.6% |

[`cell_counts_and_sizes.png`](results/figures/cell_counts_and_sizes.png) shows
all three side by side: near-identical cell *counts*, a very similar *median*
transcripts/cell, but a starkly different *distribution* (compare the means)
and overall capture.

**Transcript capture rate** = fraction of all qv≥20 transcripts in the
sub-region (770,748 total) assigned to *any* cell. CellPose only sees the DAPI
image and only segments nuclei, so cytoplasmic/extranuclear transcripts — the
majority for most genes — are never assigned to a cell, capturing barely a
third of the total. Baysor, working directly on the transcript point cloud,
captures essentially all of them (98.6%). **This is the key QC takeaway**:
CellPose and Baysor agree closely on *how many cells* are present and even on
the *median* transcript count, but CellPose's per-cell expression profiles are
built from far fewer transcripts on average, because nuclear-only segmentation
structurally excludes most of the cytoplasm.

**Matching**: 2,101 mutual-nearest-neighbor pairs (≤10 µm centroid distance)
out of 4,459 CellPose / 4,514 Baysor cells in the 1mm x 1mm sub-region.

**Expression agreement**
([`expression_correlation.png`](results/figures/expression_correlation.png)):
median per-pair Pearson correlation = **0.74** across shared genes — fairly
strong per-cell expression agreement given the two methods use completely
different inputs and capture very different numbers of transcripts per cell.

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

**Bottom line**: segmentation method choice has a real, measurable effect on
both per-cell transcript capture (driven by nuclear vs. effectively-whole-cell
capture) and downstream cell-type calls (47% of matched cells land in
different clusters depending on method), but the *spatial pattern* of that
disagreement is close to noise in this tissue — it isn't concentrated at
tumor/stroma boundaries or any other obvious structure.

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
