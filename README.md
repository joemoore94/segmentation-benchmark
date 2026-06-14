# Segmentation Benchmarking on Xenium Spatial Transcriptomics Data

## Question

Do cell segmentation methods developed for multiplexed imaging (Mesmer, CellPose,
StarDist) transfer well to imaging-based spatial transcriptomics (10x Xenium) data, and does the
choice of segmentation method meaningfully change downstream cell type calls? Where
methods disagree, is that disagreement spatially structured (e.g. concentrated at
tissue boundaries or specific niches) or essentially random noise?

This is Project 1 of a portfolio bridging imaging-based spatial biology into
sequencing-based spatial bioinformatics. Project 2
([label-transfer-benchmark](https://github.com/joemoore94/label-transfer-benchmark))
uses this project's segmented/quantified cells to evaluate how reliably
scRNA-seq reference cell-type labels transfer onto spatial data.

## Dataset

**Xenium FFPE Human Breast (Custom Add-on Panel)**, Janesick et al. 2023,
*"High resolution mapping of the breast cancer tumor microenvironment using integrated
single cell, spatial and in situ analysis of FFPE tissue"*, *Nature Communications*.

- 10x dataset page: [Xenium FFPE Human Breast with Custom Add-on Panel](https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-with-custom-add-on-panel-1-standard)
- ~577,000 cells, ~78M transcripts, registered post-Xenium H&E image
- Matched scRNA-seq + Visium from the same tissue blocks: GEO [GSE243275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243275)

Segmentation runs on a cropped ROI (~2-4mm²), selected for a mix of tumor, stroma,
and immune-infiltrated regions. See [`docs/dataset.md`](docs/dataset.md) for
download and ROI-selection details.

## Methods compared

| Method | Input | Notes |
|---|---|---|
| **10x native** | provided | Xenium Ranger's own nucleus/cell boundary segmentation (`cell_feature_matrix.h5ad` + `cells.parquet`), reshaped into `adata_10x.h5ad` by `scripts/build_10x_adata.py`; no new segmentation run needed, included below as the platform's reference |
| **CellPose** | DAPI (full 2mm x 2mm ROI) | CellPose 3.x classical `nuclei` U-Net model, CPU |
| **StarDist** | DAPI (full 2mm x 2mm ROI) | StarDist2D `2D_versatile_fluo` pretrained model (star-convex polygon regression), CPU, run via `scripts/run_stardist_roi.py` in a separate `stardist` conda env |
| **Mesmer** (DeepCell) | DAPI (ROI crop) | not run, blocked on deepcell.org's account system (signup, login, and password reset all fail as of June 2026); wrapper and conda env are ready, see `docs/dataset.md` |
| **Baysor** | transcripts (full 2mm x 2mm ROI, tiled into 4 quadrants) | transcript-density-based segmentation, run via Julia 1.10 |
| **Baysor (CellPose prior)** | transcripts + CellPose nucleus masks (full 2mm x 2mm ROI, tiled into 4 quadrants) | Baysor re-run with `--prior-segmentation` seeded from CellPose's nucleus masks (`configs/baysor_prior_config.toml`, `prior_segmentation_confidence=0.2`) |

For each method: per-cell transcript aggregation → AnnData (cells × genes) → compare
cell counts/density, size distributions, per-cell transcript counts, expression
correlation, and Leiden-based cell type calls (agreement via ARI / confusion matrices).
Finally, map where methods disagree and test whether that disagreement is spatially
structured.

## Results

Results cover **CellPose vs. Baysor vs. 10x native vs. StarDist vs. Baysor
(CellPose-prior hybrid)** (Mesmer could not be run; see Status). All five
methods now cover the same full 2mm x 2mm ROI (Baysor was extended from its
original centered 1mm x 1mm sub-region by tiling it into a 2x2 grid of
overlapping quadrants and merging the per-tile results, see
[`docs/dataset.md`](docs/dataset.md)), so cell counts, transcript capture, and
size are directly comparable across all five without density normalization.
Each method is then compared against CellPose in turn, followed by a summary
tying all four comparisons together.

### Cell counts, transcripts, and size (full 2mm x 2mm ROI)

| | CellPose (nuclear, DAPI) | Baysor (transcript density) | 10x native (whole-cell) | StarDist (nuclear, DAPI) | Baysor (CellPose prior) |
|---|---|---|---|---|---|
| Cells | 20,166 | 18,321 | 23,629 | 24,745 | 19,061 |
| Median transcripts/cell | 49 | 53 | 124 | 45 | 53 |
| Mean transcripts/cell | 59.6 | 182.6 | 142.1 | 56.0 | 175.6 |
| Median size | ~29.5 µm² nucleus | 53 transcripts/cell | ~28.1 µm² nucleus | ~27.2 µm² nucleus | 53 transcripts/cell |
| Transcript capture rate | 35.4% | 98.6% | 99.0% | 40.8% | 98.7% |

![Cell counts, transcripts/cell, and nucleus area by method](results/figures/cell_counts_and_sizes.png)

[`cell_counts_and_sizes.png`](results/figures/cell_counts_and_sizes.png) shows
all five side by side. **Transcript capture rate** = fraction of all qv≥20
transcripts in the ROI (3,392,051 total) assigned to *any* cell. CellPose and
StarDist, both nuclear-only DAPI segmentations, capture barely a third to
two-fifths of the total, since cytoplasmic transcripts (the majority for most
genes) are never assigned to a nucleus mask; Baysor, 10x native, and the
Baysor/CellPose-prior hybrid, all working from the transcript point cloud (with
or without a nuclear prior) or whole-cell boundaries, capture essentially all
of them (~99%). Median nucleus size is close across the two nuclear-mask
methods and 10x native's own nucleus segmentation (27.2-29.5 µm²): CellPose,
StarDist, and 10x native agree on *how big* a nucleus is, even though they
disagree on *how many* there are and, for 10x native, the cell boundary
extends well beyond the nucleus into the cytoplasm. **Key QC takeaway**: the
choice between nuclear-only and whole-cell/transcript-based segmentation
affects per-cell transcript counts far more than algorithm choice within a
modality (CellPose vs. StarDist differ by <6pp in capture rate; Baysor vs.
10x native vs. Baysor-prior all by <1pp of each other).

### Per-method clustering structure (PCA / UMAP)

![PCA and UMAP embeddings colored by Leiden cluster, per method](results/figures/pca_umap_clusters.png)

[`pca_umap_clusters.png`](results/figures/pca_umap_clusters.png) shows each
method's independent Leiden clustering in PCA (top row) and UMAP (bottom row)
space, the standard scRNA-seq QC view for how distinct each method's cell
populations are. All five produce well-separated UMAP clusters; cluster count
ranges from 12 (StarDist) to 24 (Baysor with the CellPose prior), with
non-prior Baysor (21), 10x native (15), and CellPose (13) in between. Baysor's
larger cluster counts -- and the CellPose-prior variant's even larger one --
likely reflect their much wider spread in transcripts/cell (mean 175.6-182.6
vs. 56-142 for the other three): with more signal per cell, finer expression
differences resolve into additional clusters.

### CellPose vs. Baysor (transcript-density segmentation)

Both methods now cover the full 2mm x 2mm ROI (20,166 CellPose / 18,321 Baysor
cells, see table above).

![Per-cell-pair expression correlation vs. CellPose, for each comparison](results/figures/expression_correlation.png)

**Matching and expression agreement**
([`expression_correlation.png`](results/figures/expression_correlation.png)):
8,947 mutual-nearest-neighbor pairs (≤10 µm centroid distance) out of 20,166
CellPose / 18,321 Baysor cells, with median per-pair Pearson correlation =
**0.73** across shared genes.

![Cell-type cluster correspondence (matched pairs)](results/figures/cell_type_confusion.png)

![Cell-type agreement (blue) vs. disagreement (red), mapped spatially](results/figures/disagreement_spatial_map.png)

**Cell-type agreement and spatial structure**
([`cell_type_confusion.png`](results/figures/cell_type_confusion.png) and
[`disagreement_spatial_map.png`](results/figures/disagreement_spatial_map.png),
left panels): independent Leiden clustering (13 CellPose, 21 Baysor clusters),
Hungarian-aligned onto a shared vocabulary, gives ARI = **0.415** with
**48.8%** of matched pairs landing in different clusters. Moran's I = **0.034**
(permutation test, p = 0.0001), barely above the noise floor: disagreement is
essentially unstructured spatially.

**Bottom line**: CellPose and Baysor agree closely on overall cell counts
(20,166 vs. 18,321, within 10%) and median transcripts/cell (49 vs. 53), but
nearly half of matched cells land in different cell-type clusters (ARI
0.415), and that disagreement is spread essentially uniformly across the ROI
rather than concentrated in any particular tissue region.

### CellPose vs. Baysor (CellPose-prior hybrid)

A natural follow-up: does seeding Baysor's transcript-density EM with
CellPose's own nucleus masks as a `--prior-segmentation`
(`configs/baysor_prior_config.toml`, `prior_segmentation_confidence=0.2`)
pull its segmentation closer to CellPose's? 35.4% of transcripts fall within a
CellPose nucleus and get a non-zero prior label; the rest are unconstrained,
as in the non-prior run. This hybrid covers the same full 2mm x 2mm ROI via
the same 4-tile scheme and produces **19,061 cells**, 98.7% transcript
capture -- both essentially unchanged from the non-prior run (18,321 cells,
98.6%, see table above).

**Matching and expression agreement**
([`expression_correlation.png`](results/figures/expression_correlation.png),
fourth panel): **9,572** mutual-nearest-neighbor pairs out of 20,166 CellPose
/ 19,061 Baysor (prior) cells, about 7% more than the non-prior run's 8,947
despite the hybrid having fewer total cells. Median Pearson correlation =
**0.74** (vs. 0.73 without the prior), essentially unchanged.

**Cell-type agreement and spatial structure**
([`cell_type_confusion.png`](results/figures/cell_type_confusion.png) and
[`disagreement_spatial_map.png`](results/figures/disagreement_spatial_map.png),
fourth panels): independent Leiden clustering gives 24 Baysor (prior) clusters
(vs. 21 without the prior), Hungarian-aligned as before, giving ARI = **0.407**
(vs. 0.415) with **49.7%** of matched pairs in different clusters (vs. 48.8%).
Moran's I = **0.051** (p = 0.0001), modestly higher than the non-prior run's
0.034 but still close to the noise floor.

**Bottom line**: a CellPose-nucleus prior pulls Baysor's segmentation into
slightly better spatial correspondence with CellPose (7% more matched pairs)
but leaves the underlying cell-type-call agreement essentially unchanged (ARI
0.407 vs. 0.415, ~49% disagreement either way) -- the two methods' roughly 50%
disagreement rate on matched cells reflects a deeper difference in how each
defines a "cell" (nuclear footprint vs. transcript-density neighborhood), not
something fixable by anchoring Baysor to CellPose's nuclei.

### CellPose vs. 10x native (platform reference)

10x's own (Xenium Ranger) segmentation is already part of the downloaded
dataset; `scripts/build_10x_adata.py` reshapes `cell_feature_matrix.h5ad` and
`cells.parquet` into the same AnnData schema used above, no new segmentation
run needed. Both cover the full 2mm x 2mm ROI (20,166 CellPose / 23,629 10x
native cells, see table above).

**Median nucleus area is nearly identical** (29.5 vs. 28.1 µm²)
([`cell_counts_and_sizes.png`](results/figures/cell_counts_and_sizes.png),
right panel): CellPose's nuclear segmentation produces nuclei essentially the
same size as 10x's own. The much higher transcripts/cell and capture rate for
10x native instead come from its *whole-cell* boundary, which extends into the
cytoplasm that CellPose's nuclear-only masks exclude.

**Matching and expression agreement**
([`expression_correlation.png`](results/figures/expression_correlation.png),
second panel): 18,966 mutual-nearest-neighbor pairs (≤10 µm centroid distance)
out of 20,166 CellPose / 23,629 10x native cells: **94%** of CellPose's cells
have a corresponding 10x cell. Median Pearson correlation = **0.82**, notably
higher than the 0.73 vs. Baysor, consistent with CellPose and 10x native both
being nuclear-pixel-mask segmentations of the *same* DAPI image.

**Cell-type agreement and spatial structure**
([`cell_type_confusion.png`](results/figures/cell_type_confusion.png) and
[`disagreement_spatial_map.png`](results/figures/disagreement_spatial_map.png),
second panels): ARI = **0.547** (vs. 0.415 for Baysor), with **30.8%** of
matched pairs landing in different clusters (vs. 48.8%): CellPose tracks the
platform's reference segmentation more closely than it tracks Baysor's.
Moran's I = **0.176** (p = 0.0001), the highest Moran's I of any comparison in
this project: disagreement here is the most spatially clustered, plausibly
reflecting tissue regions (e.g. dense tumor nests) where nuclear segmentation
is intrinsically harder and CellPose's calls drift from the platform
reference.

**Takeaway**: CellPose reproduces 10x's nucleus size and cell-type calls
reasonably well (ARI 0.55), but for analyses that depend on whole-cell
expression, nuclear-only segmentation is a bigger limitation than algorithm
disagreement.

### CellPose vs. StarDist (same input, different algorithm)

StarDist (`2D_versatile_fluo`, star-convex polygon regression) is run on the
*same* full 2mm x 2mm DAPI image as CellPose, the closest apples-to-apples
comparison in this project: two algorithms on identical input.

StarDist finds **more, slightly smaller nuclei** than CellPose (24,745 vs.
20,166 cells, 27.2 vs. 29.5 µm² median nucleus area): transcript capture
(40.8% vs. 35.4%) and median transcripts/cell (45 vs. 49) both sit close to
CellPose's, consistent with both methods segmenting the same nuclei with a
modest difference in how aggressively each splits touching ones.

**Matching and expression agreement**
([`expression_correlation.png`](results/figures/expression_correlation.png),
third panel): 19,460 mutual-nearest-neighbor pairs (≤10 µm centroid distance)
out of 20,166 CellPose / 24,745 StarDist cells: **96.5%** of CellPose's cells
have a corresponding StarDist cell, the highest overlap of any pair here.
Median Pearson correlation = **0.96**, by far the highest of any comparison in
this project (0.73 vs. Baysor, 0.82 vs. 10x native, 0.74 vs. Baysor prior),
reflecting how similarly the two algorithms trace nucleus boundaries on the
same DAPI image, not a difference in *what* each captures.

**Cell-type agreement and spatial structure**
([`cell_type_confusion.png`](results/figures/cell_type_confusion.png) and
[`disagreement_spatial_map.png`](results/figures/disagreement_spatial_map.png),
third panels): independent Leiden clustering gave 12 StarDist clusters (vs. 13
for CellPose), Hungarian-aligned as before, giving ARI = **0.764** (again the
highest of any comparison here) with only **13.8%** of matched pairs in
different clusters (vs. 30.8% for 10x native, 48.8% for Baysor, 49.7% for
Baysor prior). Moran's I = **0.066** (p = 0.0001), intermediate between
CellPose vs. Baysor (0.034) and CellPose vs. 10x native (0.176), despite
StarDist disagreeing far less often overall.

**Takeaway**: two algorithms on the same image agree closely (ARI 0.76,
expression correlation 0.96); algorithm choice is a minor source of
disagreement compared to differences in segmentation modality or input.

### Summary across method pairs

All four CellPose-anchored comparisons side by side:

| Comparison | Matched pairs | Median expression corr. | ARI | Disagreement rate | Moran's I (p) |
|---|---|---|---|---|---|
| CellPose vs. Baysor | 8,947 | 0.73 | 0.415 | 48.8% | 0.034 (0.0001) |
| CellPose vs. 10x native | 18,966 | 0.82 | 0.547 | 30.8% | 0.176 (0.0001) |
| CellPose vs. StarDist | 19,460 | 0.96 | 0.764 | 13.8% | 0.066 (0.0001) |
| CellPose vs. Baysor (prior) | 9,572 | 0.74 | 0.407 | 49.7% | 0.051 (0.0001) |

Expression correlation, ARI, and disagreement rate rank **StarDist (closest)
> 10x native > Baysor ≈ Baysor (prior) (furthest)**, tracking how much
CellPose and the other method share in input and segmentation modality (same
image + same algorithm → same image + different modality → different modality
entirely). The two Baysor variants are essentially tied with each other on all
three metrics (within 1pp of disagreement rate and 0.01 of ARI), confirming
that the CellPose-nucleus prior does not meaningfully change Baysor's
downstream cell-type agreement with CellPose.

Spatial structure of disagreement does not track that ranking. The two Baysor
variants have the highest disagreement rates (48.8% / 49.7%) but among the
lowest Moran's I values (0.034 / 0.051), barely above the permutation-test
noise floor; StarDist, with the lowest disagreement rate (13.8%), has an
intermediate Moran's I (0.066); and 10x native, intermediate on disagreement
rate (30.8%), has the *highest* spatial autocorrelation (0.176). Extending
Baysor to the full ROI roughly halved its Moran's I relative to the earlier
1mm² sub-region estimate (0.066 → 0.034), suggesting that estimate was partly
a small-sample artifact: with full-ROI coverage, Baysor's disagreement with
CellPose, though the most frequent of any comparison here, is also among the
most spatially uniform -- and the CellPose-prior hybrid, despite anchoring
directly to CellPose's own nuclei, doesn't change that. 10x native remains the
standout: disagreement with the platform's own reference segmentation is the
most spatially concentrated, plausibly reflecting tissue regions (e.g. dense
tumor nests) where nuclear segmentation is intrinsically hardest. Segmentation
method choice clearly affects cell-type calls, but whether that effect is
spatially "interesting" vs. just noise depends on *which* methods are
compared, not simply on how much they disagree.

### Phenotypic density vs. segmentation disagreement (Mellon)

A different question from *how much* the methods disagree is *where in
phenotypic space* disagreement happens. Using
[Mellon](https://github.com/settylab/mellon) (Otto et al. 2024, *Nature
Methods*), each CellPose cell gets a per-cell log-density estimate in PCA
space, a measure of how typical vs. phenotypically rare/transitional its
expression profile is. The hypothesis: cell-type disagreement should
concentrate on rare, ambiguous cells near cluster boundaries rather than on
well-defined, high-density cell types. Only CellPose cells with a matched
partner in the other method have a `disagree` label, so the unmatched
remainder of each CellPose cell are excluded from this analysis.

![CellPose phenotypic density (Mellon) vs. cell-type call disagreement](results/figures/density_vs_disagreement.png)

[`density_vs_disagreement.png`](results/figures/density_vs_disagreement.png)
compares the log-density distribution of agreeing vs. disagreeing cells for
each pair (Mann-Whitney U, two-sided):

| Comparison | n agree / disagree | Median log-density (agree / disagree) | p-value |
|---|---|---|---|
| CellPose vs. Baysor | 4,585 / 4,362 | -15.00 / -17.01 | 2.3e-76 |
| CellPose vs. 10x native | 13,121 / 5,845 | -13.28 / -13.54 | 3.5e-11 |
| CellPose vs. StarDist | 16,780 / 2,680 | -13.39 / -13.27 | 6.1e-09 |
| CellPose vs. Baysor (prior) | 4,812 / 4,760 | -14.77 / -17.03 | 4.4e-91 |

**Bottom line**: the hypothesis holds clearly for both Baysor comparisons
([`density_vs_disagreement.png`](results/figures/density_vs_disagreement.png),
first and fourth panels): disagreeing cells sit roughly two log-units lower in
density, by far the largest gap of any comparison here, with p = 2.3e-76
(non-prior) and p = 4.4e-91 (prior) -- the prior's larger matched set (9,572
vs. 8,947) only sharpens the effect. For CellPose vs. 10x native and CellPose
vs. StarDist the p-values are also significant (unsurprising with sample sizes
in the tens of thousands), but the median gaps are tiny (≤0.3 log-units), and
the StarDist comparison runs in the *opposite* direction (disagreeing cells
slightly *denser*, not sparser). Density-driven disagreement looks like a real
effect specific to cross-modality comparisons -- present in both the
transcript-density Baysor run and its CellPose-prior hybrid -- not a general
property of segmentation disagreement.

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
│   ├── segmentation/          # per-method wrappers (cellpose, stardist, mesmer, baysor)
│   ├── quantify.py            # transcript aggregation -> per-cell AnnData
│   ├── compare.py             # cross-method comparison metrics
│   └── spatial.py             # spatial structure of disagreement
├── scripts/                   # CLI entry points (download, run_stardist_roi.py, run_mesmer.py, run_baysor.sh)
├── results/{figures,tables}/  # outputs (gitignored, regeneratable)
└── tests/
```

## Environment setup

This project uses four toolchains: a conda env for CellPose + the
Scanpy/Squidpy/SpatialData stack, a separate conda env for StarDist
(TensorFlow-based), a separate conda env for Mesmer (DeepCell's TensorFlow
pin is incompatible with the main env), and Julia for Baysor.

### 1. Python (conda, for CellPose + Scanpy/Squidpy/SpatialData stack)

```bash
conda env create -f environment.yml
conda activate segbench
```

### 2. StarDist (separate conda env)

StarDist depends on TensorFlow, which conflicts with the PyTorch-based main
env, so it runs in its own env:

```bash
conda create -n stardist python=3.10
conda run -n stardist pip install stardist tensorflow-cpu
```

Run via
`conda run -n segbench python scripts/run_stardist_roi.py`, which calls
`segbench.segmentation.stardist_run.run_stardist` (a `conda run -n stardist
python scripts/run_stardist.py ...` subprocess wrapper, mirroring the Mesmer
pattern below).

### 3. Mesmer (DeepCell, separate conda env)

```bash
conda create -n mesmer python=3.10
conda run -n mesmer pip install deepcell
```

Requires a `DEEPCELL_ACCESS_TOKEN` (free account at
[users.deepcell.org](https://users.deepcell.org)) to download `Mesmer`'s
pretrained weights on first use. As of June 2026, deepcell.org's account
system is non-functional (signup returns a server error, login and password
reset both fail); the env and wrapper (`segbench.segmentation.mesmer_run`,
`scripts/run_mesmer.py`) are ready to use once that's resolved.

### 4. Julia + Baysor (for transcript-based segmentation)

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
- [x] Segmentation: CellPose, Baysor, StarDist
- [x] Quantification + cross-method comparison (CellPose vs. Baysor)
- [x] Spatial analysis of disagreement
- [x] Figures + write-up
- [ ] Segmentation: Mesmer, blocked externally on deepcell.org's account
      system (signup/login/password-reset all fail); wrapper and conda env
      are ready, revisit if/when access is restored
- [x] Stretch: incorporate 10x native segmentation as a third reference
      (CellPose vs. 10x native: 18,966 matched pairs, ARI = 0.547)
- [x] Stretch: incorporate StarDist as a fourth segmentation method
      (CellPose vs. StarDist: 19,460 matched pairs, ARI = 0.764, closest
      agreement of any pair in this project)
- [x] Stretch: extend Baysor from a centered 1mm x 1mm sub-region to the full
      2mm x 2mm ROI via 4-tile merge (18,321 cells, 98.6% transcript capture);
      CellPose vs. Baysor now has 8,947 matched pairs (up from 2,101)
- [x] Stretch: Mellon phenotypic-density analysis of disagreement (CellPose
      vs. Baysor: disagreeing cells ~2 log-units lower density, p = 2.3e-76;
      CellPose vs. 10x native / StarDist: significant but negligible effect)
- [x] Stretch: PCA/UMAP visualization of each method's independent Leiden
      clustering (`pca_umap_clusters.png`)
- [x] Stretch: Baysor re-run with a CellPose-nucleus `--prior-segmentation`
      hybrid (19,061 cells, 98.7% transcript capture); CellPose vs.
      Baysor (prior) has 9,572 matched pairs (ARI = 0.407, Moran's I = 0.051,
      p = 4.4e-91 for Mellon density vs. disagreement) -- essentially
      unchanged from non-prior Baysor, so the prior improves spatial
      correspondence with CellPose without changing cell-type agreement
- [ ] Future: Kompot differential-expression analysis between high- and
      low-density cell groups (installed, not yet run)
