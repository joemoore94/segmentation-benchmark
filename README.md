# Segmentation Benchmarking on Xenium Spatial Transcriptomics Data

**Question:** Do segmentation methods optimized for imaging — nuclear pixel-mask models (CellPose, StarDist, Mesmer), simple Voronoi nearest-centroid assignment, and transcript-density EM (Baysor) — transfer well to imaging-based spatial transcriptomics (10x Xenium), and does method choice meaningfully change downstream cell-type calls?

## Key findings

| Comparison | Matched pairs | Median corr | ARI | Disagreement rate | Moran's I |
|---|---|---|---|---|---|
| 10x native vs. CellPose | 18,966 | 0.822 | 0.547 | 30.8% | 0.178 |
| 10x native vs. StarDist | 21,429 | 0.826 | 0.545 | 33.5% | 0.215 |
| 10x native vs. Mesmer | 20,595 | 0.879 | 0.557 | 27.9% | 0.090 |
| 10x native vs. Voronoi (CellPose) | 18,966 | 0.959 | 0.630 | 21.9% | 0.076 |
| 10x native vs. Voronoi (Mesmer) | 20,595 | 0.964 | 0.686 | 18.8% | 0.161 |
| 10x native vs. Baysor | 10,953 | 0.786 | 0.305 | 51.7% | 0.033 |

*Matched pairs*: cells matched by nearest centroid across methods. *Median corr*: median per-pair Pearson correlation of log-normalised expression profiles. *ARI*: Adjusted Rand Index of cluster-label agreement after Hungarian alignment (0 = random, 1 = perfect). *Moran's I*: spatial autocorrelation of the disagree flag (0 = random, 1 = fully clustered).

Three method families emerge. Nuclear methods (CellPose, StarDist, Mesmer) score ARI ~0.55 with spatially structured disagreement (Moran's I 0.09–0.22) concentrated in the luminal epithelial population — which in this breast cancer tissue likely encompasses the malignant cells — where dense cell packing and cytoplasmic transcripts outside nuclear masks create consistent boundary ambiguity. Voronoi variants (nearest-centroid assignment from CellPose or Mesmer nuclei) reach ARI 0.63–0.69 with 100% transcript capture and diffuse, boundary-driven residual error. Baysor reaches ARI 0.31 with near-random spatial disagreement; relaxing the one-to-one cluster matching constraint (many-to-one assignment) recovers ~8% of its apparent disagreement as an over-clustering artefact, bringing effective accuracy from 48% to 56%.

**Voronoi (Mesmer) is the best non-reference method at ARI 0.686.** The gain over CellPose nuclear decomposes cleanly: adding cytoplasmic coverage (CellPose nuclear → Voronoi CP) contributes +0.083 ARI; improving nuclear centroid quality (Voronoi CP → Voronoi M) contributes an additional +0.056 ARI. Voronoi (Mesmer)'s residual disagreement retains a weak luminal-epithelial signal (density p=3.2e-12) versus nuclear Mesmer's strong one (p=3.8e-79), confirming that cytoplasmic coverage is the dominant driver.

<!-- Project 2 (label-transfer-benchmark): uses this project's segmented cells to evaluate scRNA-seq label-transfer reliability. Add link once repo is public. -->

## Dataset

**Xenium FFPE Human Breast (Custom Add-on Panel)**, Janesick et al. 2023, *Nature Communications* ([dataset page](https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-with-custom-add-on-panel-1-standard)). ~577,000 cells and ~78M transcripts across the full slide (invasive ductal carcinoma); matched scRNA-seq + Visium from the same tissue blocks: GEO [GSE243275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243275).

All segmentation and analysis runs on a 2mm x 2mm ROI (~23,600 cells, ~3.4M transcripts) with a mix of tumor, stroma, and immune-infiltrated regions. See [`docs/dataset.md`](docs/dataset.md) for download and ROI details.

## Methods

| Method | Input | Notes |
|---|---|---|
| **10x native** | provided | Xenium Ranger's own segmentation, reshaped by `scripts/build_10x_adata.py` |
| **CellPose** | DAPI (2mm x 2mm ROI) | CellPose 3.x `nuclei` model, CPU |
| **StarDist** | DAPI (2mm x 2mm ROI) | `2D_versatile_fluo` model, separate `stardist` conda env |
| **Voronoi (CellPose)** | CellPose nuclear centroids | Nearest-centroid transcript assignment via scipy cKDTree; 100% capture, no additional model |
| **Mesmer** (DeepCell) | DAPI | run via Docker (`vanvalenlab/deepcell-applications`); image bundles model weights, no access token needed |
| **Voronoi (Mesmer)** | Mesmer nuclear centroids | Same nearest-centroid assignment using Mesmer centroids; isolates nuclear detector quality from cytoplasmic coverage effect |
| **Baysor** | transcripts (2mm x 2mm, 4 tiles) | transcript-density EM, Julia 1.10 |

Per-cell transcript aggregation → AnnData → cell counts, transcript capture, expression correlation, Leiden clustering, and spatial structure of disagreement (Moran's I, Mellon density). Independent Leiden runs assign arbitrary cluster IDs with no shared meaning across methods; cluster labels are aligned using the Hungarian algorithm, which finds the one-to-one relabelling of comparison cluster IDs that maximises total overlap with 10x native, before computing disagreement rate and ARI.

## Results

### Cell counts and transcript capture

| | CellPose | StarDist | Mesmer | Voronoi (CP) | Voronoi (M) | Baysor | 10x native |
|---|---|---|---|---|---|---|---|
| Cells | 20,166 | 24,745 | 21,697 | 20,166 | 21,697 | 18,321 | 23,629 |
| Median transcripts/cell | 49 | 45 | 70 | 149 | 142 | 53 | 124 |
| Transcript capture | 35.4% | 40.8% | 51.8% | 100.0% | 100.0% | 98.6% | 99.0% |

![Cell counts, transcripts/cell, and nucleus area by method](results/figures/cell_counts_and_sizes.png)

Nuclear-only methods (CellPose, StarDist, Mesmer) capture 35-52% of transcripts; Mesmer's nuclear masks are calibrated larger than CellPose or StarDist, yielding higher capture without leaving nuclear-only mode. Voronoi variants capture 100% by construction; Baysor and 10x native capture 98-99%.

### Clustering structure

![UMAP embeddings colored by Leiden cluster, per method](results/figures/pca_umap_clusters.png)

Leiden clustering runs independently on each method's cells, grouping them into communities in PCA space. Because each method detects different cells with different transcript counts, the resulting cluster structures differ even when the underlying biology is the same — that divergence is what the pairwise comparisons measure. All methods produce well-separated UMAP embeddings. Cluster counts: CellPose 13, StarDist 12, Mesmer 15, Voronoi (CP) 14, Voronoi (M) 14, Baysor 21, 10x native 15. Baysor's higher count reflects its richer per-cell transcript profiles resolving finer expression differences; Mesmer and 10x native both converge on 15 clusters.

### Pairwise comparisons (all vs. 10x native)

All comparisons use 10x native (Xenium Ranger's own segmentation) as the reference anchor. Each comparison method's cells are matched to the nearest 10x-native cell by centroid; expression correlation is computed per matched pair; Leiden clusters are aligned by Hungarian algorithm before computing disagreement rate and ARI.

![Per-cell-pair expression correlation](results/figures/expression_correlation.png)
![Disagreement mapped spatially](results/figures/disagreement_spatial_map.png)
![Annotated cluster confusion matrices](results/figures/confusion_annotated.png)

**10x native vs. CellPose** (whole-cell vs. nuclear): 18,966 matched pairs, median expression correlation 0.822, ARI 0.547, 30.8% disagreement, Moran's I 0.178. Disagreement is spatially structured and concentrated in the luminal epithelial region; the confusion matrix shows reasonable one-to-one cluster alignment with several high-overlap pairs.

**10x native vs. StarDist** (whole-cell vs. nuclear): 21,429 matched pairs, correlation 0.826, ARI 0.545, 33.5% disagreement, Moran's I 0.215. Metrics and spatial disagreement pattern are nearly identical to CellPose; the marginally higher Moran's I reflects slightly more tightly clustered hotspots in the luminal epithelial region.

**10x native vs. Mesmer** (whole-cell vs. nuclear, DeepCell): 20,595 matched pairs, correlation 0.879, ARI 0.557, 27.9% disagreement, Moran's I 0.090. Mesmer outperforms CellPose and StarDist on every metric despite running in nuclear-only mode, largely because its larger nuclear masks capture ~52% of transcripts vs. 35-41% for the other nuclear methods. Its disagreement is spatially structured (Moran's I 0.090) with the same luminal-epithelial fingerprint (MYBPC1, SERPINA3, CLIC6, PGR, GATA3), and it has the highest fraction of agreement coldspots (32.5% LL) of any method.

**10x native vs. Voronoi (CellPose)** (whole-cell vs. nearest-centroid assignment from CellPose nuclei): 18,966 matched pairs, correlation 0.959, ARI 0.630, 21.9% disagreement, Moran's I 0.076. Voronoi assigns all transcripts to the nearest CellPose nuclear centroid, capturing the full cytoplasmic signal with no additional model. The substantially higher ARI and correlation relative to nuclear CellPose (same 20,166 cells, same centroids) directly quantifies the contribution of cytoplasmic transcripts to cell-type identity.

**10x native vs. Voronoi (Mesmer)** (whole-cell vs. nearest-centroid assignment from Mesmer nuclei): 20,595 matched pairs, correlation 0.964, ARI 0.686, 18.8% disagreement, Moran's I 0.161. Swapping Mesmer's higher-quality nuclear centroids into the same Voronoi framework improves ARI by 0.056 over Voronoi (CellPose), isolating the contribution of nuclear detector quality from cytoplasmic coverage. The Moran's I (0.161) is higher than nuclear Mesmer's (0.090) despite lower overall disagreement — Voronoi (Mesmer)'s remaining errors are spatially concentrated in luminal-epithelial patches (density p=3.2e-12, DE top genes: MUC1, SERPINA3, CLIC6, PGR, MYBPC1) rather than scattered randomly, but the effect is far weaker than nuclear Mesmer (p=3.8e-79), confirming that cytoplasmic coverage dominates.

**10x native vs. Baysor**: 10,953 matched pairs, correlation 0.786, ARI 0.305, 51.7% disagreement, Moran's I 0.033. More than half of matched cells land in different clusters, and the pattern is near-random spatially. Baysor resolves 21 clusters versus 10x native's 15; relaxing the one-to-one Hungarian constraint (many-to-one assignment, each comparison cluster independently assigned to its argmax 10x cluster) raises accuracy from 48.3% to 56.2%, confirming that ~8% of the apparent disagreement is a cluster-count artefact. The remaining ~44% is genuine: Baysor's transcript-density boundaries produce systematically different cell-type assignments concentrated on macrophages (CD14, MRC1, CD163).

### Which cell types disagree?

![Cell type vs. agreement with CellPose](results/figures/agreement_explainer.png)

Adipocytes and myoepithelial cells carry the highest per-cell disagreement rates (~50-68% and ~40-47% respectively across nuclear methods), but as relatively rare populations they contribute modestly to total disagreement events. Luminal epithelial cells are the dominant contributor by volume: at roughly 35% per-cell disagreement and ~8,500 cells (~37% of all 10x-native cells), they account for the majority of disagreements across every nuclear method. In this invasive ductal carcinoma tissue (Janesick et al. 2023), the luminal epithelial population at Leiden resolution 1.0 likely encompasses malignant cells alongside residual normal luminal epithelial — both share canonical markers (GATA3, PGR, ESR1, MUC1) and are not separable by nuclear morphology alone. These cells form dense, overlapping clusters where nuclear boundary ambiguity is highest and cytoplasmic transcripts invisible to nuclear-only masks carry the most discriminative expression. Comparing the cell type map (top left, pink region) with each method's spatial map (left column) shows the luminal epithelial territory overlapping directly with each method's dominant disagreement zone. T cells and B cells are robustly identified regardless of method, consistent with their isolated nuclei and highly distinctive transcriptional signatures (CD3E, TRAC, MS4A1).

### Phenotypic density vs. disagreement (Mellon)

![10x native phenotypic density vs. disagreement](results/figures/density_vs_disagreement.png)

Each 10x-native cell gets a Mellon log-density estimate in PCA space; disagreeing vs. agreeing cells compared by Mann-Whitney U:

| Comparison | n agree / disagree | Median log-density (agree / disagree) | p |
|---|---|---|---|
| 10x native vs. CellPose | 13,121 / 5,845 | -21.31 / -20.78 | 2.9e-28 |
| 10x native vs. StarDist | 14,254 / 7,175 | -21.87 / -20.63 | 1.1e-90 |
| 10x native vs. Mesmer | 14,850 / 5,745 | -21.73 / -20.14 | 3.8e-79 |
| 10x native vs. Voronoi (CellPose) | 14,805 / 4,161 | -21.05 / -21.35 | 0.191 (n.s.) |
| 10x native vs. Voronoi (Mesmer) | 16,720 / 3,875 | -21.38 / -20.68 | 3.2e-12 |
| 10x native vs. Baysor | 5,286 / 5,667 | -22.76 / -22.75 | 0.756 (n.s.) |

Nuclear methods (CellPose, StarDist, Mesmer) disagree with 10x native on cells in *higher*-density phenotypic regions (p ≪ 0.001), driven by luminal epithelial cells whose cytoplasmic expression is captured by whole-cell segmentation but missed by nuclear-only masks. Mesmer's effect size is the largest (median density gap 1.59 log units vs. 0.53 for CellPose), consistent with its larger nuclear masks creating more ambiguity at cell boundaries in dense epithelial regions. Voronoi (CellPose)'s disagreement is density-neutral (p=0.19) — its remaining error is geometric, not cell-state-driven. Voronoi (Mesmer) shows a weak but significant residual density effect (p=3.2e-12), still driven by luminal epithelial cells but much attenuated compared to nuclear Mesmer; the same luminal-epithelial marker signature (MUC1, SERPINA3, CLIC6) appears in its DE results. Baysor shows no density effect; its disagreement concentrates on macrophages (CD14, MRC1, CD163). T cells (TRAC, CD3E) are robustly identified by all methods.

### Local Moran's I (LISA)

![LISA hotspot/coldspot maps per comparison](results/figures/local_morans_map.png)

HH clusters (local disagreement hotspots) and LL clusters (local agreement coldspots) per comparison:

| Comparison | HH hotspots | LL coldspots |
|---|---|---|
| 10x native vs. CellPose | 21.7% | 30.3% |
| 10x native vs. StarDist | 18.6% | 15.0% |
| 10x native vs. Mesmer | 17.1% | 32.5% |
| 10x native vs. Voronoi (CellPose) | 11.1% | 27.2% |
| 10x native vs. Voronoi (Mesmer) | 9.5% | 20.4% |
| 10x native vs. Baysor | 21.4% | 17.5% |

Mesmer has the most agreement coldspots of any method (32.5% LL) — large contiguous tissue regions where Mesmer and 10x native call identical cell types — consistent with its overall lower disagreement rate and spatially coherent nuclear detection. Voronoi (Mesmer) has the fewest HH hotspots of any method (9.5%), reflecting that its remaining disagreement is diffuse boundary error rather than concentrated failure zones. Baysor's near-equal HH/LH split confirms the near-random spatial structure of transcript-density disagreement.

### Differential expression: agree vs. disagree cells

![Volcano plots: DE genes in disagree vs. agree groups](results/figures/de_volcano.png)

Wilcoxon rank-sum test comparing log-normalised expression in disagreeing vs. agreeing 10x-native cells for each comparison. For nuclear methods, top markers in disagreeing cells are consistently luminal epithelial genes (MYBPC1, SERPINA3, CLIC6, PGR, GATA3, MUC1) — cytoplasmic transcripts that nuclear-only masks undercount relative to 10x native's whole-cell segmentation. Baysor's disagreeing cells are enriched instead for macrophage markers (CD14, MRC1, CD163), consistent with the density analysis. Voronoi methods show fewer significant genes and smaller fold changes, confirming their residual disagreement is geometric rather than cell-state-driven.

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

This project uses three toolchains: a main conda env for CellPose + Scanpy/Squidpy/SpatialData, a separate env for StarDist (TensorFlow-based), and Julia for Baysor. Mesmer runs via Docker (no separate conda env required).

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
sudo apt-get install -y docker.io
sudo systemctl enable --now docker
sudo usermod -aG docker $USER
# log out and back in (or run: exec su -l $USER) for group membership to take effect
docker pull vanvalenlab/deepcell-applications:latest
```

The `vanvalenlab/deepcell-applications` image bundles pretrained model weights and does not require a `DEEPCELL_ACCESS_TOKEN`. See [`scripts/run_mesmer.sh`](scripts/run_mesmer.sh) for the invocation.

### 4. Julia + Baysor

```bash
juliaup add 1.10
julia +1.10 -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git", rev="v0.7.1")); Pkg.build("Baysor")'
```

See [`scripts/run_baysor.sh`](scripts/run_baysor.sh) for the invocation.

---

## Supplemental: exploratory methods

These methods were evaluated during development and informed the final method selection but are not part of the primary analysis.

### Baysor with CellPose nucleus prior

Baysor run with `--prior-segmentation` from CellPose nuclear masks at `prior_segmentation_confidence=0.2`. The nucleus prior adds ~5% more matched pairs (10,953 → 11,454) and marginally improves ARI (0.305 → 0.318) but leaves the fundamental disagreement pattern unchanged: 51.9% disagreement, Moran's I 0.036, near-random spatial structure. The nucleus prior is not the bottleneck for Baysor's performance gap; its disagreement is explained by macrophages and other low-transcript-count cells where EM-based transcript density segmentation differs systematically from morphology-based boundaries.

| Metric | Baysor | Baysor (CellPose prior) |
|---|---|---|
| Matched pairs | 10,953 | 11,454 |
| Median corr | 0.786 | 0.798 |
| ARI | 0.305 | 0.318 |
| Disagreement rate | 51.7% | 51.9% |
| Moran's I | 0.033 | 0.036 |

Scripts: `scripts/run_baysor.sh` (with `configs/baysor_prior_config.toml`).

### Nucleus expansion (10µm, 20µm)

CellPose nuclear masks expanded outward by 10µm (~47px) and 20µm (~94px) using `skimage.segmentation.expand_labels`, then re-quantified. This tests whether simple morphological dilation of nuclear masks closes the cytoplasmic-transcript gap.

| Method | ARI vs. 10x native |
|---|---|
| CellPose (nuclear only) | 0.547 |
| CellPose + 10µm expansion | 0.572 |
| CellPose + 20µm expansion | 0.592 |
| Voronoi (nearest-centroid) | 0.630 |

Expansion improves ARI monotonically but never reaches Voronoi, which assigns all transcripts without a fixed radius. Voronoi was selected for the primary analysis. Script: `scripts/build_expanded_adatas.py`.

### Baysor prior confidence sensitivity

Baysor rerun with `prior_segmentation_confidence` at 0.5 and 0.8 (default 0.2) to test whether stronger nuclear guidance improves agreement with 10x native.

| Config | ARI | Cells |
|---|---|---|
| Baysor (prior, c=0.2) | 0.318 | 19,061 |
| Baysor (prior, c=0.5) | 0.395 | — |
| Baysor (prior, c=0.8) | 0.488 | 29,771 |

c=0.8 inflates cell count to 29,771 (vs. 10x native's 23,629), a likely artifact of the prior overriding Baysor's own cell boundary inference at high confidence. Neither setting reaches Voronoi. Scripts: `scripts/run_baysor_conf_sensitivity.sh`, configs `configs/baysor_prior_config_c05.toml` / `configs/baysor_prior_config_c08.toml`.
