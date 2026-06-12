# Segmentation Benchmarking on Xenium Spatial Transcriptomics Data

## Question

Do cell segmentation methods developed for multiplexed imaging (Mesmer, CellPose)
transfer well to imaging-based spatial transcriptomics (10x Xenium) data, and does the
choice of segmentation method meaningfully change downstream cell type calls? Where
methods disagree, is that disagreement spatially structured (e.g. concentrated at
tissue boundaries or specific niches) or essentially random noise?

This is Project 1 of a two-part series bridging imaging-based spatial biology
([Canonix](#)) into sequencing-based spatial bioinformatics. Project 2 will use this
project's output for scRNA-seq label-transfer reliability analysis.

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
| **10x native** | provided | Xenium's own nucleus/cell boundary segmentation — used as a reference, not "ground truth" |
| **CellPose** | DAPI + membrane channels (ROI crop) | run natively in the `segbench` conda env, MPS-accelerated where supported |
| **Mesmer** (DeepCell) | DAPI + membrane channels (ROI crop) | run via Docker (`vanvalenlab/deepcell-applications`) — no native Apple Silicon build exists |
| **Baysor** | transcripts within ROI (optionally DAPI-seeded) | transcript-density-based segmentation, run via Julia |

For each method: per-cell transcript aggregation → AnnData (cells × genes) → compare
cell counts, size distributions, per-cell transcript counts, spatial overlap (IoU) vs.
10x reference, expression correlation, and Leiden-based cell type calls (agreement via
ARI / confusion matrices). Finally, map where methods disagree and test whether that
disagreement is spatially structured.

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
├── scripts/                   # CLI entry points (download, run_mesmer.sh, run_baysor.sh)
├── results/{figures,tables}/  # outputs (gitignored, regeneratable)
└── tests/
```

## Environment setup

This project uses three toolchains, since Mesmer (TensorFlow/DeepCell) has no native
Apple Silicon build and Baysor is a Julia package.

### 1. Python (conda, for CellPose + Scanpy/Squidpy/SpatialData stack)

```bash
conda env create -f environment.yml
conda activate segbench
```

### 2. Docker (for Mesmer via `vanvalenlab/deepcell-applications`)

```bash
docker pull vanvalenlab/deepcell-applications:latest-cpu
```

See [`scripts/run_mesmer.sh`](scripts/run_mesmer.sh) for the invocation.

### 3. Julia + Baysor (for transcript-based segmentation)

```bash
juliaup add 1.10
julia +1.10 -e 'using Pkg; Pkg.add(url="https://github.com/kharchenkolab/Baysor.git")'
```

See [`scripts/run_baysor.sh`](scripts/run_baysor.sh) for the invocation.

## Status

- [x] Project scaffold + environments
- [ ] Data acquisition + ROI selection
- [ ] Segmentation runs (CellPose, Mesmer, Baysor, 10x reference)
- [ ] Quantification + cross-method comparison
- [ ] Spatial analysis of disagreement
- [ ] Write-up
