# Segmentation Benchmarking on Xenium Spatial Transcriptomics Data

**Question:** Do segmentation methods transfer well to Xenium spatial transcriptomics, does method choice meaningfully change downstream cell-type calls, and can we benchmark segmentation quality without relying on any single reference?

## Evaluation framework

Methods are evaluated along two orthogonal axes: a **factorial design** that decomposes segmentation into nuclear detection × expansion strategy, and a **three-anchor benchmark** that avoids privileging any single reference.

### N-detector × N-expansion factorial design

Every detector is paired with every expansion strategy, cleanly separating detection quality from expansion strategy. Four nuclear detectors (CellPose, StarDist, Mesmer, 10x Ranger) × expansion strategies (nuclear-only, geometric 10µm, geometric 20µm, Voronoi, watershed, Baysor PSC 0.2/0.5/0.8/1.0), plus whole-cell methods (Cellpose cyto3, Mesmer whole-cell) that produce masks directly.

### Three independent reference anchors

| Anchor | Information source | Strengths | Blind spots |
| --- | --- | --- | --- |
| **10x native** | Proprietary algorithm (DAPI-based) | Practical standard, most users' default | Black-box expansion, no ground truth |
| **H&E morphology** | DAPI + eosin whole-cell segmentation | Real cytoplasm boundaries from tissue morphology | Pre-IF registration artifacts, z-overlap |
| **scRNA-seq reference** | Dissociated single-cell data (GSE243275) | True single-cell expression profiles | Dissociation bias, loses spatial context |

Each method is evaluated against all three anchors. The anchors are also compared to each other: if all three agree that a method is good, that's robust. If they disagree, the pattern of disagreement reveals what each anchor is actually measuring — a method that scores high against 10x native but low against scRNA-seq is reproducing the platform's biases rather than recovering true cell states.

## Summary

The 4-detector × N-expansion factorial design cleanly separates nuclear detection quality from expansion strategy. Mesmer produces the best nuclear detections regardless of expansion method (Voronoi M: ARI 0.686, Baysor M prior: ARI 0.518), while 10x Ranger — despite being purpose-built for Xenium — is comparable to CellPose and StarDist. Voronoi consistently outperforms Baysor PSC=1.0 on ARI (0.58–0.69 vs 0.50–0.53) and cell-level displacement in scRNA-seq reference PCA space (0.92–1.41 vs 2.01–2.97 units). Voronoi methods also reproduce 10x native's cell type composition within 1–3 percentage points, while Baysor without a prior fundamentally reshapes it: luminal epithelial drops from 36.5% to 10.7% and CAFs rise from 28.3% to 37.2%.

The negative marker analysis — the one metric that requires no reference method — reveals a third ranking. Whole-cell NN methods, which score only modestly on Benchmark 1 (ARI 0.54–0.62), have by far the lowest cross-lineage violation rates: Cellpose cyto3 (DAPI) reaches 0.77% (0.14 per 1000 tx), a factor of 3–8 lower than any Voronoi or geometric method. Baysor prior variants at PSC ≥ 0.8 are next-lowest (1.79–2.34%), while Voronoi methods (4.49–6.45%) and 10x native itself (4.50%) cluster together. The pattern confirms that methods which reproduce 10x native's DAPI-expansion logic (high B1 ARI, high B2 negative delta) also produce more cross-lineage boundary contamination. The cross-lineage bleed analysis confirms these violations are boundary artifacts: T cells near CAFs show elevated CAF-marker expression under Voronoi but not Baysor, with bleed intensity correlated to proximity to the nearest CAF.

Baysor's prior_segmentation_confidence parameter exhibits a graduated threshold. Below PSC 0.2, the prior has almost no effect. Between PSC 0.2 and 0.5, violation rates begin dropping (4.71% → 3.28%) even though ARI remains near the no-prior baseline. At PSC ≥ 0.8, cell counts nearly double, ARI jumps from 0.35 to 0.49, and violations drop further to 2.17–2.23%. Among the four luminal epithelial subclusters, cluster 3 (NNMT/LUM/POSTN) was flagged as ambiguous — either genuine EMT or segmentation contamination. A cross-method comparison of stromal-to-luminal marker ratios (ranging from 2.90 for Baysor M prior to 3.69 for Voronoi CP) confirms that the stromal signal is at least partly a boundary artifact: methods with broader expansion pull proportionally more CAF transcripts into these cells.

Method ordering is stable across Leiden resolutions 0.3–2.0. In scRNA-seq reference PCA space (Benchmark 3), all nuclear-seeded expanded methods — Voronoi (3.48–3.54), geometric expansion (3.47–3.56), watershed (3.47–3.58) — land in essentially the same B3 band as 10x native (3.50). B3 primarily distinguishes three tiers: nuclear-seeded expansion methods are indistinguishable from 10x native; whole-cell NN methods form an intermediate tier (3.68–4.30, with Mesmer WC closest and Cellpose cyto3 DAPI-only furthest); and Baysor methods remain ~0.8 units further regardless of PSC or detector (Baysor PSC ≥ 0.8 at 4.22–4.38, PSC ≤ 0.5 at 4.61–4.65). B3 does not differentiate within the expanded-nuclear cluster.

The H&E morphology anchor (Benchmark 2) produces a fundamentally different ranking. Whole-cell NN methods — which scored modestly on Benchmark 1 (ARI 0.54–0.62) — top the H&E ranking: Cellpose cyto3 DAPI+density reaches ARI 0.807, Mesmer WC density 0.776. Critically, 10x native itself scores only 0.590 against the H&E anchor — lower than 22 of 31 tested methods — revealing that the platform's proprietary expansion does not follow cytoplasm boundaries. Methods with the largest positive Benchmark 2 − Benchmark 1 delta are capturing real morphology; methods with negative deltas (Voronoi CP, Voronoi M, 10x native) are reproducing the platform's DAPI-expansion logic rather than cell shape.

### Full results table (Benchmark 1: 10x native anchor)

| Comparison | Matched pairs | Median corr | ARI | Hungarian Disagree | Hungarian Moran's I | Argmax Disagree | Argmax Moran's I |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **Nuclear-only detectors** | | | | | | | |
| 10x native vs. CellPose | 18,966 | 0.822 | 0.547 | 33.5% | 0.169 | 32.8% | 0.189 |
| 10x native vs. StarDist | 21,429 | 0.826 | 0.545 | 31.2% | 0.094 | 30.6% | 0.111 |
| 10x native vs. Mesmer | 20,595 | 0.879 | 0.557 | 27.3% | 0.095 | 26.7% | 0.096 |
| 10x native vs. 10x Ranger | 23,155 | 0.822 | 0.504 | 41.3% | 0.209 | 40.5% | 0.293 |
| **Voronoi expansion** | | | | | | | |
| 10x native vs. Voronoi (CP) | 18,966 | 0.959 | 0.630 | 28.1% | 0.158 | 27.7% | 0.159 |
| 10x native vs. Voronoi (SD) | 21,428 | 0.959 | 0.584 | 32.7% | 0.203 | 27.2% | 0.212 |
| 10x native vs. Voronoi (M) | 20,595 | 0.964 | 0.686 | 27.7% | 0.220 | 27.1% | 0.212 |
| 10x native vs. Voronoi (10x) | 23,153 | 0.967 | 0.592 | 27.0% | 0.221 | 27.0% | 0.221 |
| **Geometric expansion 10µm** | | | | | | | |
| 10x native vs. Expansion 10µm (CP) | 19,343 | 0.981 | 0.604 | 27.0% | 0.198 | 21.7% | 0.134 |
| 10x native vs. Expansion 10µm (SD) | 21,818 | 0.980 | 0.631 | 27.6% | 0.226 | 27.2% | 0.224 |
| 10x native vs. Expansion 10µm (M) | 20,860 | 0.978 | 0.670 | 20.8% | 0.114 | 20.4% | 0.112 |
| 10x native vs. Expansion 10µm (10x) | 23,624 | 0.993 | 0.691 | 21.7% | 0.120 | 19.2% | 0.112 |
| **Geometric expansion 20µm** | | | | | | | |
| 10x native vs. Expansion 20µm (CP) | 19,343 | 0.982 | 0.603 | 27.6% | 0.113 | 22.6% | 0.133 |
| 10x native vs. Expansion 20µm (SD) | 21,825 | 0.981 | 0.644 | 26.2% | 0.223 | 25.8% | 0.221 |
| 10x native vs. Expansion 20µm (M) | 20,856 | 0.979 | 0.680 | 20.1% | 0.143 | 19.2% | 0.141 |
| 10x native vs. Expansion 20µm (10x) | 23,624 | 0.994 | 0.710 | 20.0% | 0.153 | 19.2% | 0.192 |
| **Watershed expansion** | | | | | | | |
| 10x native vs. Watershed (10x) | 23,142 | 0.979 | 0.664 | 25.1% | 0.221 | 24.1% | 0.220 |
| 10x native vs. Watershed (SD) | 21,393 | 0.966 | 0.656 | 21.6% | 0.095 | 20.7% | 0.098 |
| 10x native vs. Watershed (M) | 20,494 | 0.965 | 0.674 | 20.0% | 0.110 | 19.6% | 0.107 |
| **Baysor PSC sweep (CP prior)** | | | | | | | |
| 10x native vs. Baysor (no prior) | 10,953 | 0.786 | 0.305 | 49.7% | 0.049 | 44.9% | 0.106 |
| 10x native vs. Baysor (CP prior 0.2) | 11,454 | 0.798 | 0.318 | 55.4% | 0.049 | 42.9% | 0.103 |
| 10x native vs. Baysor (CP prior 0.5) | 15,316 | 0.798 | 0.349 | 52.0% | 0.068 | 44.2% | 0.100 |
| 10x native vs. Baysor (CP prior 0.8) | 20,108 | 0.884 | 0.488 | 41.9% | 0.134 | 40.0% | 0.164 |
| 10x native vs. Baysor (CP prior 1.0) | 20,308 | 0.902 | 0.501 | 37.8% | 0.146 | 36.2% | 0.162 |
| **Baysor other detectors at PSC=1.0** | | | | | | | |
| 10x native vs. Baysor (SD prior 1.0) | 21,814 | 0.905 | 0.498 | 39.7% | 0.116 | 35.3% | 0.147 |
| 10x native vs. Baysor (M prior 1.0) | 21,148 | 0.924 | 0.518 | 37.5% | 0.137 | 34.9% | 0.162 |
| 10x native vs. Baysor (10x prior 1.0) | 22,910 | 0.914 | 0.530 | 37.0% | 0.139 | 35.3% | 0.154 |
| **Whole-cell NN** | | | | | | | |
| 10x native vs. Cellpose cyto3 (DAPI) | 19,765 | 0.854 | 0.540 | 31.2% | 0.103 | 30.3% | 0.126 |
| 10x native vs. Cellpose cyto3 (DAPI+eosin) | 16,740 | 0.911 | 0.598 | 31.6% | 0.200 | 30.7% | 0.214 |
| 10x native vs. Cellpose cyto3 (DAPI+density) | 19,714 | 0.903 | 0.569 | 32.6% | 0.187 | 31.9% | 0.209 |
| 10x native vs. Mesmer WC (DAPI+eosin) | 20,579 | 0.927 | 0.590 | 31.2% | 0.199 | 30.7% | 0.203 |
| 10x native vs. Mesmer WC (DAPI+density) | 20,316 | 0.913 | 0.617 | 25.7% | 0.109 | 25.3% | 0.111 |

*Matched pairs*: nearest-centroid matching (max 10µm). *Median corr*: per-pair Pearson correlation of log-normalised expression. *ARI*: Adjusted Rand Index between each method's Leiden clustering and 10x native's Leiden clustering, computed on matched cells. *Disagreement*: fraction of matched pairs assigned to different clusters after label alignment. Nuclear-only methods capture 35–52% of transcripts and are excluded from downstream figures past the recovery section.

---

## Methods & Basic Statistics

### Dataset

**Xenium FFPE Human Breast (Custom Add-on Panel)**, Janesick et al. 2023, *Nature Communications* ([dataset page](https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-with-custom-add-on-panel-1-standard)). Invasive ductal carcinoma; matched scRNA-seq + Visium from the same tissue blocks: GEO [GSE243275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243275).

All analysis runs on a 2mm × 2mm ROI (~23,600 cells, ~3.4M transcripts, 380-gene panel) with a mix of tumor, stroma, and immune-infiltrated regions. See [`docs/dataset.md`](docs/dataset.md) for download and ROI details.

### Methods

| Method | Input | Notes |
| --- | --- | --- |
| **10x native** | provided | Xenium Ranger's full segmentation (nuclear detection + proprietary expansion); used as Benchmark 1 anchor. |
| **10x Ranger** | DAPI | Nuclear detection component of Xenium Ranger, extracted from `nucleus_boundaries.parquet`. Included as a fourth nuclear detector. |
| **CellPose** | DAPI | CellPose 3.x `nuclei` model, CPU |
| **StarDist** | DAPI | `2D_versatile_fluo` model, separate `stardist` env |
| **Mesmer** | DAPI | DeepCell via Docker; image bundles model weights |
| **Voronoi (CP / SD / M / 10x)** | nuclear centroids | Nearest-centroid transcript assignment; 100% capture by construction |
| **Geometric 10µm / 20µm (CP / SD / M / 10x)** | nuclear mask | `skimage.expand_labels` with 10 or 20 µm radius; deterministic, fills to midpoint between cells |
| **Watershed (10x / SD / M)** | nuclear mask + DAPI | DAPI gradient as landscape, nuclear mask as seeds; `compactness=0.01`; 100% capture |
| **Baysor (prior 0.2 / 0.5 / 0.8 / 1.0)** | transcripts + nuclear masks | Transcript-density EM with `prior_segmentation_confidence` sweep. Tested with all four nuclear detectors at PSC 1.0. |
| **Cellpose cyto3** | DAPI / DAPI+eosin / DAPI+density | CellPose whole-cell `cyto3` model with three input variants |
| **Mesmer WC** | DAPI+eosin / DAPI+density | DeepCell whole-cell application via Docker; eosin channel from registered H&E provides cytoplasm signal |

Nuclear-only methods (CellPose, StarDist, Mesmer, 10x Ranger) capture only 35–52% of transcripts and are included in the recovery section but excluded from downstream figures because their low transcript capture dominates any comparison.

### Prior strength sweep

Baysor's `prior_segmentation_confidence` (PSC) controls how strongly the nuclear prior constrains the density model. At PSC=0.0 (no prior), Baysor runs purely on transcript density; at PSC=1.0, nuclear transcripts are hard-locked and only cytoplasmic transcripts use density-adaptive expansion. All CellPose-prior variants were run at PSC 0.0 (no prior), 0.2, 0.5, 0.8, and 1.0.

| PSC | Cells | Matched pairs | Median tx/cell | Median corr | ARI | Disagree (H) | Moran's I (H) | Tier 1 violations |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.0 | 18,321 | 10,953 | 53 | 0.786 | 0.305 | 51.7% | 0.033 | 5.01% |
| 0.2 | 19,061 | 11,454 | 53 | 0.798 | 0.318 | 51.9% | 0.036 | 4.71% |
| 0.5 | 24,147 | — | 49 | — | — | — | — | 3.28% |
| 0.8 | 29,771 | 20,108 | 67 | 0.884 | 0.488 | 38.8% | 0.219 | 2.23% |
| 1.0 | 30,473 | 20,308 | 69 | 0.902 | 0.501 | 33.8% | 0.111 | 2.17% |

The sweep reveals a sharp transition between PSC 0.2 and 0.8, with PSC=0.5 sitting in the transition zone. At PSC ≤ 0.2, the prior barely constrains the density model: cell counts, transcript capture, and ARI are nearly identical to the no-prior baseline, and negative marker violations remain high (4.71–5.01%). Between PSC 0.2 and 0.5, violation rates drop notably from 4.71% to 3.28% — indicating the prior begins to limit cross-boundary transcript merging — but ARI and matched-pair counts remain near the no-prior regime. At PSC ≥ 0.8, the nuclear prior locks enough transcripts to prevent merging: cell counts jump to ~30,000, matched pairs nearly double, ARI increases from 0.32 to 0.49, and violations drop further to 2.17–2.23%. The Moran's I jump from 0.036 to 0.219 at PSC=0.8 indicates that the prior converts Baysor's diffuse, spatially random disagreement into structured disagreement concentrated at tissue boundaries. Full B1 comparison metrics for PSC=0.5 are pending.

### Cell and transcript recovery

**Nuclear detectors**

| Detector | Cells | Median tx/cell | Transcript capture |
| --- | ---: | ---: | ---: |
| CellPose | 20,166 | 49 | 35.4% |
| StarDist | 24,745 | 45 | 40.8% |
| Mesmer | 21,697 | 70 | 51.8% |
| 10x Ranger | 23,624 | 45 | 38.0% |

<p align="center"><img src="results/figures/nuclear_mask_sizes.png" alt="Nuclear mask size comparison across all four detectors" width="600"></p>

All four detectors operate on the same DAPI image but produce substantially different masks. Mesmer detects the largest nuclei (median ~45 µm², long tail to 200 µm²), capturing 51.8% of transcripts — nearly double CellPose's 35.4%. StarDist and 10x Ranger produce similar-sized masks but StarDist finds more nuclei (24,745 vs 23,624). 10x Ranger captures only 38% of transcripts despite detecting nearly as many cells as the 10x native whole-cell segmentation (23,624 vs 23,629), confirming that 10x native's 99% capture comes from its proprietary expansion, not from larger nuclei.

**Whole-cell and expansion methods**

| Method | Cells | Median tx/cell | Transcript capture |
| --- | ---: | ---: | ---: |
| 10x native (Ranger + proprietary expansion) | 23,629 | 124 | 99.0% |
| Voronoi (CP) | 20,166 | 149 | 100% |
| Voronoi (SD) | 24,745 | 122 | 100% |
| Voronoi (M) | 21,697 | 142 | 100% |
| Voronoi (10x) | 23,622 | 128 | 100% |
| Geometric 10µm (CP) | 20,166 | — | 95.9% |
| Geometric 20µm (CP) | 20,166 | — | 99.6% |
| Geometric 10µm (SD) | 24,745 | — | 95.9% |
| Geometric 20µm (SD) | 24,745 | — | 99.6% |
| Geometric 10µm (M) | 21,697 | — | 96.8% |
| Geometric 20µm (M) | 21,697 | — | 99.7% |
| Geometric 10µm (10x) | 23,624 | — | 96.2% |
| Geometric 20µm (10x) | 23,624 | — | 99.7% |
| Watershed (10x) | 23,142 | — | 100.0% |
| Watershed (SD) | 24,745 | — | 100.0% |
| Watershed (M) | 21,697 | — | 100.0% |
| Baysor (no prior) | 18,321 | 53 | 98.6% |
| Baysor (CP prior 0.2) | 19,061 | 53 | 98.7% |
| Baysor (CP prior 0.5) | 24,147 | 49 | 97.3% |
| Baysor (CP prior 0.8) | 29,771 | 67 | 99% |
| Baysor (CP prior 1.0) | 30,473 | 69 | 99% |
| Baysor (SD prior 1.0) | 34,230 | 63 | 98.9% |
| Baysor (M prior 1.0) | 31,764 | 74 | 98.9% |
| Baysor (10x prior 1.0) | 33,113 | 65 | 98.9% |
| Cellpose cyto3 (DAPI) | 21,641 | 55 | 40.3% |
| Cellpose cyto3 (DAPI+eosin) | 17,450 | 77 | 45.0% |
| Cellpose cyto3 (DAPI+density) | 20,782 | 76 | 53.8% |
| Mesmer WC (DAPI+eosin) | 21,841 | 91 | 68.3% |
| Mesmer WC (DAPI+density) | 21,493 | 86 | 63.7% |

![Transcripts/cell distribution across expansion methods](results/figures/transcripts_per_cell.png)

Voronoi expansion captures 100% of transcripts by construction. Geometric expansion at 10µm captures ~96% and at 20µm ~99.6–99.7% regardless of detector; watershed fills the entire ROI (100%). Whole-cell NN methods are more conservative: Cellpose cyto3 captures only 41–55% (real cytoplasm boundaries leave extracellular gaps), while Mesmer WC captures 64–68% — substantially higher than Cellpose, suggesting Mesmer's training data better matches the tissue architecture. Baysor without a prior captures 98.6%, collapsing to ~11,000 matched cells due to cell merging.

---

## Benchmark 1: 10x Native Anchor

The 10x native segmentation is the platform's default output, produced by Xenium Ranger's proprietary DAPI-based algorithm. All methods are compared against it via centroid matching. Disagreement and spatial structure analyses throughout this section anchor on 10x native labels. This is the most direct measure of compatibility with the platform standard, but it treats 10x native as ground truth — any systematic biases in the platform's own algorithm are invisible to this benchmark.

### Per-cell expression correlation

<p align="center"><img src="results/figures/expression_correlation.png" alt="Per-cell expression correlation vs. 10x native" width="600"></p>

Per-cell expression correlation is high for all methods (median 0.79–0.97). Voronoi methods lead (median r 0.96–0.97), with Voronoi (10x) highest at 0.967 because it shares nuclear seeds with 10x native and differs only in expansion algorithm. Baysor PSC=1.0 variants cluster at 0.90–0.92, and nuclear-only methods at 0.82–0.88, reflecting missing cytoplasmic transcripts rather than misassignment.

### Clustering comparison

Leiden clustering runs independently on each method's cells (normalize → PCA → neighbors → Leiden at resolution 1.0). Cluster labels are aligned across methods before computing confusion matrices and disagreement, using two algorithms: Hungarian (one-to-one) and argmax (many-to-one).

#### Resolution stability

<p align="center"><img src="results/figures/resolution_ari.png" alt="ARI across Leiden resolutions (alignment-agnostic)" width="600"></p>

ARI is partition-based and does not depend on cluster alignment, so it is the same under Hungarian and argmax. The method ordering is stable across Leiden resolutions 0.3–2.0. Voronoi (Mesmer) leads at most resolutions (0.3, 0.6, 0.8–1.2); at resolutions 0.5 and 0.7, Voronoi (StarDist) briefly takes the lead, and at 1.5+ StarDist's higher cell count gives it a durable advantage as finer clustering demands more cells per cluster. Baysor without a prior is consistently lowest.

![Disagreement and Moran's I across Leiden resolutions - Hungarian alignment](results/figures/resolution_disagree_morans_hungarian.png)

![Disagreement and Moran's I across Leiden resolutions - argmax alignment](results/figures/resolution_disagree_morans_argmax.png)

Disagreement and Moran's I do depend on alignment. The Hungarian alignment forces unmatched clusters into poor pairings when cluster counts differ, inflating disagreement for methods that produce more clusters. The argmax alignment lets multiple clusters map to the same reference cluster, reducing this artifact. The Moran's I panel confirms that the spatial-structure gap is resolution-invariant under both algorithms: Voronoi and Baysor prior methods maintain spatially structured disagreement while Baysor without a prior stays near zero regardless of cluster granularity.

UMAP embeddings colored by aligned cluster labels illustrate how the alignment algorithm reshapes cluster identity. Baysor without a prior shows the starkest contrast: Hungarian forces 6 of its 21 clusters into empty pairings, leaving large regions unmatched (gray), while argmax lets multiple Baysor clusters map to the same reference cluster, producing coherent coloring across the manifold.

![Baysor UMAP: Hungarian vs argmax alignment](results/figures/umap_baysor_alignment_comparison.png)

**All method UMAPs - Hungarian alignment**

![10x native](results/figures/umap/umap_10x_native.png)
![Voronoi (CellPose)](results/figures/umap/umap_voronoi_hungarian.png)
![Voronoi (StarDist)](results/figures/umap/umap_voronoi_stardist_hungarian.png)
![Voronoi (Mesmer)](results/figures/umap/umap_voronoi_mesmer_hungarian.png)
![Voronoi (10x Ranger)](results/figures/umap/umap_voronoi_10x_ranger_hungarian.png)
![Baysor (no prior)](results/figures/umap/umap_baysor_hungarian.png)
![Baysor (CellPose prior)](results/figures/umap/umap_baysor_prior_hungarian.png)
![Baysor (CellPose prior, PSC=0.8)](results/figures/umap/umap_baysor_prior_c08_hungarian.png)
![Baysor (CellPose prior, PSC=1.0)](results/figures/umap/umap_baysor_prior_c10_hungarian.png)
![Baysor (StarDist prior, PSC=1.0)](results/figures/umap/umap_baysor_stardist_prior_c10_hungarian.png)
![Baysor (Mesmer prior, PSC=1.0)](results/figures/umap/umap_baysor_mesmer_prior_c10_hungarian.png)
![Baysor (10x Ranger prior, PSC=1.0)](results/figures/umap/umap_baysor_10x_ranger_prior_c10_hungarian.png)

**All method UMAPs - argmax alignment**

![10x native](results/figures/umap/umap_10x_native.png)
![Voronoi (CellPose)](results/figures/umap/umap_voronoi_argmax.png)
![Voronoi (StarDist)](results/figures/umap/umap_voronoi_stardist_argmax.png)
![Voronoi (Mesmer)](results/figures/umap/umap_voronoi_mesmer_argmax.png)
![Voronoi (10x Ranger)](results/figures/umap/umap_voronoi_10x_ranger_argmax.png)
![Baysor (no prior)](results/figures/umap/umap_baysor_argmax.png)
![Baysor (CellPose prior)](results/figures/umap/umap_baysor_prior_argmax.png)
![Baysor (CellPose prior, PSC=0.8)](results/figures/umap/umap_baysor_prior_c08_argmax.png)
![Baysor (CellPose prior, PSC=1.0)](results/figures/umap/umap_baysor_prior_c10_argmax.png)
![Baysor (StarDist prior, PSC=1.0)](results/figures/umap/umap_baysor_stardist_prior_c10_argmax.png)
![Baysor (Mesmer prior, PSC=1.0)](results/figures/umap/umap_baysor_mesmer_prior_c10_argmax.png)
![Baysor (10x Ranger prior, PSC=1.0)](results/figures/umap/umap_baysor_10x_ranger_prior_c10_argmax.png)

#### Cluster alignment

**Confusion matrices (Hungarian and argmax)**

![Confusion matrices with Hungarian and argmax alignment](results/figures/confusion_clusters.png)

Each row is one 10x native cluster; columns are the comparison method's clusters. Red cells mark Hungarian (one-to-one) matched pairs, blue cells mark argmax (many-to-one) matches, and purple cells mark pairs selected by both algorithms. Voronoi methods produce clean matches under both algorithms. Baysor's 15×21 matrix shows the key difference: under Hungarian, 6 clusters are forced into empty pairings; under argmax, every column maps to the highest-overlap reference cluster with no wasted assignments.

**Clustering agreement vs. 10x native**

![Clustering agreement vs. 10x native](results/figures/clustering_agreement.png)

Voronoi methods achieve the highest ARI (0.584–0.686), with Voronoi (M) leading. Nuclear-only methods cluster at ARI 0.504–0.557, and Baysor without a prior is lowest at 0.305. Argmax reduces Baysor's disagreement by ~8pp (51.7% to 43.8%) by eliminating forced mismatches from unmatched clusters. Voronoi methods with matched cluster counts are barely affected. The Moran's I increase for Baysor under argmax (0.033 to 0.079) shows that removing alignment noise reveals spatially structured disagreement that was previously masked.

**Per-cluster pseudobulk**

<p align="center"><img src="results/figures/pseudobulk_by_cluster.png" alt="Per-cluster pseudobulk correlation vs. 10x native" width="600"></p>

To test whether cluster-level expression profiles agree, matched cells are grouped by 10x native's 15 Leiden clusters and pseudobulked per method. Nuclear methods drop to r = 0.86–0.87 on luminal epithelial clusters (0, 1, 3, 8) — the same populations driving single-cell disagreement — while Voronoi variants stay above 0.99 across all clusters. Baysor shows a comparable luminal dip plus reduced correlation on macrophage clusters (2, 7), consistent with transcript-density boundaries partitioning those populations differently.

### Cell type composition vs. 10x native

To directly test whether segmentation method changes the apparent cell type composition of the tissue, each method's cells are annotated independently using `sc.tl.score_genes` on canonical marker panels. This bypasses Leiden clustering entirely — each cell is assigned to the cell type whose marker panel it scores highest on.

| Method | Lum. epi. | CAFs | Macro. | Endo. | T cells | Myoepi. | Sm. muscle | B cells | Plasma | Adipo. |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10x native | 36.5 | 28.3 | 10.4 | 6.3 | 6.1 | 3.8 | 4.4 | 1.9 | 1.4 | 0.9 |
| Voronoi (CP) | 39.6 | 27.2 | 10.2 | 6.5 | 4.4 | 3.9 | 4.4 | 1.5 | 1.4 | 0.9 |
| Voronoi (SD) | 36.1 | 29.3 | 9.9 | 6.5 | 5.6 | 3.9 | 4.4 | 1.9 | 1.4 | 0.9 |
| Voronoi (M) | 37.1 | 28.9 | 10.3 | 6.4 | 5.3 | 3.7 | 4.2 | 1.9 | 1.3 | 0.9 |
| Voronoi (10x) | 36.8 | 29.1 | 9.7 | 6.5 | 5.7 | 3.9 | 4.1 | 1.9 | 1.4 | 0.9 |
| Baysor (no prior) | 10.7 | 37.2 | 20.0 | 5.7 | 9.1 | 5.8 | 3.7 | 3.0 | 2.3 | 2.4 |
| Baysor (CP prior 1.0) | 28.7 | 28.4 | 14.8 | 5.1 | 8.0 | 4.8 | 4.5 | 2.6 | 1.6 | 1.4 |
| Baysor (M prior 1.0) | 28.4 | 28.9 | 14.3 | 5.3 | 8.2 | 5.0 | 4.3 | 2.6 | 1.7 | 1.5 |

Voronoi methods reproduce 10x native's cell type composition within 1–3 percentage points across all cell types. The largest Voronoi deviation is Voronoi (CP) enriching luminal epithelial cells to 39.6% (vs 36.5%), consistent with fewer, larger cells concentrating more epithelial transcripts per cell.

Baysor without a prior fundamentally reshapes the tissue's apparent composition: luminal epithelial cells drop from 36.5% to 10.7% while CAFs rise to 37.2% and macrophages double to 20.0%. The density model's smaller cells (median 53 tx/cell) accumulate too few luminal markers to outscore the ambient stromal signature, so cells at tissue boundaries are reclassified from epithelial to stromal. Baysor with nuclear priors partially recovers (luminal 28.4–28.7%), but even at PSC=1.0 the composition remains shifted.

**Spatial cell type maps**

![Cell type annotation across segmentation methods](results/figures/celltype_spatial_by_method.png)

The spatial maps confirm the proportion shifts are not diffuse — they have clear tissue structure. Under 10x native and Voronoi methods, the tumor nests are dominated by luminal epithelial cells (pink), surrounded by CAF-rich stroma (gray) with immune infiltrate at the margins. Baysor without a prior fragments these epithelial regions: the density model's smaller cells lose enough luminal signal that many are reclassified as CAFs or macrophages, breaking up the coherent tumor structure.

### Spatial structure of disagreement

![Spatial autocorrelation of disagreement - all methods](results/figures/spatial_morans_dotplot.png)

**Spatial disagreement maps (Hungarian alignment)**

![Disagreement mapped spatially - Hungarian](results/figures/disagreement_spatial_map.png)

**Spatial disagreement maps (Argmax alignment)**

![Disagreement mapped spatially - argmax](results/figures/disagreement_spatial_map_argmax.png)

**LISA hotspot/coldspot maps (Hungarian)**

![LISA hotspot/coldspot maps - Hungarian](results/figures/local_morans_map.png)

**LISA hotspot/coldspot maps (argmax)**

![LISA hotspot/coldspot maps - argmax](results/figures/local_morans_map_argmax.png)

| Comparison | Global Moran's I | HH hotspots | LL coldspots |
| --- | --- | --- | --- |
| 10x native vs. CellPose | 0.178 | 21.7% | 30.3% |
| 10x native vs. StarDist | 0.215 | 18.6% | 15.0% |
| 10x native vs. Mesmer | 0.090 | 17.1% | 32.5% |
| 10x native vs. Voronoi (CP) | 0.076 | 11.1% | 27.2% |
| 10x native vs. Voronoi (SD) | 0.194 | 23.2% | 30.8% |
| 10x native vs. Voronoi (M) | 0.161 | 9.5% | 20.4% |
| 10x native vs. Baysor | 0.033 | 21.4% | 17.5% |

Nuclear and Voronoi disagreements are spatially structured (Moran's I 0.076–0.215 under Hungarian), concentrated in luminal epithelial territory. Mesmer has the most agreement coldspots (32.5% LL); Voronoi (Mesmer) has the fewest disagreement hotspots (9.5% HH), consistent with residual errors being diffuse boundary noise. Under Hungarian alignment Baysor's near-zero Moran's I (0.033) reflects noise from forced cluster mismatches; under argmax alignment Moran's I increases to 0.079, revealing that Baysor's genuine disagreements are spatially structured — just less so than morphological methods.

*LISA labels*: Local Moran's I (Anselin 1995) decomposes the global statistic to a per-cell level. HH = disagreeing cell surrounded by disagreeing neighbors (hotspot), LL = agreeing cell surrounded by agreeing neighbors (coldspot), HL/LH = spatial outliers.

### Cell-type sensitivity

**Cell type vs. agreement (Hungarian)**

![Cell type vs. agreement - Hungarian](results/figures/agreement_explainer.png)

**Cell type vs. agreement (Argmax)**

![Cell type vs. agreement - argmax](results/figures/agreement_explainer_argmax.png)

Adipocytes and myoepithelial cells have the highest per-cell disagreement (~50–68% and ~40–47%) but are rare. Luminal epithelial cells dominate by volume: ~35% disagreement across ~8,500 cells drives the majority of total disagreement events. T cells and B cells are robustly identified regardless of method or alignment algorithm.

### Disagreement drivers: cell state vs. geometry

**Phenotypic density vs. disagreement (Hungarian)**

![Phenotypic density vs. disagreement - Hungarian](results/figures/density_vs_disagreement.png)

**Phenotypic density vs. disagreement (Argmax)**

![Phenotypic density vs. disagreement - argmax](results/figures/density_vs_disagreement_argmax.png)

**DE: disagree vs. agree cells (Hungarian)**

![DE: disagree vs. agree cells - Hungarian](results/figures/de_volcano.png)

| Comparison | n agree / disagree | Median log-density (agree / disagree) | p |
| --- | --- | --- | --- |
| 10x native vs. CellPose | 13,121 / 5,845 | -21.31 / -20.78 | 2.9e-28 |
| 10x native vs. StarDist | 14,254 / 7,175 | -21.87 / -20.63 | 1.1e-90 |
| 10x native vs. Mesmer | 14,850 / 5,745 | -21.73 / -20.14 | 3.8e-79 |
| 10x native vs. Voronoi (CP) | 14,805 / 4,161 | -21.05 / -21.35 | 0.191 n.s. |
| 10x native vs. Voronoi (SD) | 14,597 / 6,831 | -21.74 / -20.56 | 5.5e-51 |
| 10x native vs. Voronoi (M) | 16,720 / 3,875 | -21.38 / -20.68 | 3.2e-12 |
| 10x native vs. Baysor | 5,286 / 5,667 | -22.76 / -22.75 | 0.756 n.s. |

Nuclear methods disagree on cells in higher-density phenotypic regions (Mann-Whitney p ≪ 0.001). The DE volcano confirms this: disagreeing cells are enriched for luminal epithelial markers (MYBPC1, SERPINA3, CLIC6, PGR, GATA3, MUC1), cytoplasmic transcripts underrepresented in nuclear-only masks. Voronoi (CellPose) disagreement is density-neutral (p = 0.19) with few DE genes, indicating residual errors are geometric. Baysor disagreement is also density-neutral but enriched for macrophage markers (CD14, MRC1, CD163), consistent with transcript-density boundaries partitioning macrophage-rich regions differently.

### Cell size and disagreement

**Cell size vs. disagreement (Hungarian)**

![Cell size vs. disagreement - Hungarian](results/figures/cell_size_disagreement.png)

**Cell size vs. disagreement (Argmax)**

![Cell size vs. disagreement - argmax](results/figures/cell_size_disagreement_argmax.png)

| Comparison | Median area (agree) | Median area (disagree) | p |
| --- | --- | --- | --- |
| 10x native vs. CellPose | 123.9 µm² | 121.4 µm² | 7.2e-07 |
| 10x native vs. StarDist | 121.2 µm² | 116.9 µm² | 2.4e-12 |
| 10x native vs. Mesmer | 126.2 µm² | 119.1 µm² | 3.2e-07 |
| 10x native vs. Voronoi (CP) | 126.8 µm² | 111.5 µm² | 4.9e-32 |
| 10x native vs. Voronoi (SD) | 123.8 µm² | 112.4 µm² | 8.8e-29 |
| 10x native vs. Voronoi (M) | 125.9 µm² | 117.7 µm² | 3.7e-19 |
| 10x native vs. Baysor | 167.0 µm² | 173.4 µm² | 0.28 n.s. |

Smaller 10x-native cells are significantly more likely to disagree with every morphological method (p ≪ 0.001). Smaller cells likely correspond to densely packed regions where any method's cluster assignment is noisier. Baysor shows no size dependence (p = 0.28); its boundaries are insensitive to morphologically defined cell area.

### Pairwise method agreement

<p align="center"><img src="results/figures/pairwise_consensus.png" alt="Pairwise ARI between all segmentation methods" width="600"></p>

Voronoi expansion does not improve cross-method reproducibility. CellPose and StarDist agree with each other at ARI 0.764 — higher than the Voronoi (CP) vs Voronoi (SD) pair at 0.661 — because both are nuclear-morphology methods operating on the same DAPI image. What Voronoi does raise is agreement with the 10x-native whole-cell reference (0.63–0.69): compatibility with the platform's own segmentation, not cross-method reproducibility. Baysor remains isolated from all morphological methods (ARI 0.30–0.46 regardless of partner).

---

## Benchmark 2: H&E Morphology Anchor

The H&E anchor is grounded in real cytoplasm and membrane signal from the registered H&E image. Eosin stains cytoplasm, stroma, and ECM; hematoxylin stains nuclei. After color deconvolution (`skimage.color.rgb2hed()`), the eosin channel provides a cytoplasm signal unavailable in the DAPI-only Xenium image. **Mesmer WC DAPI+eosin** is used as the H&E reference segmentation: it is the highest-performing eosin-based method (68.3% transcript capture vs. Cellpose cyto3's 45.6%) and uses the only available real cytoplasm channel.

Methods that score highly against this anchor capture real morphological boundaries; methods that score high against 10x native but low against this anchor are reproducing the platform's DAPI-expansion logic rather than tissue morphology.

### Registration quality

The H&E image was acquired pre-IF from registered tissue sections. Color deconvolution and pixel-level correlation were computed by `scripts/build_eosin_channel.py`. Hematoxylin vs. DAPI Pearson r = **0.896** (n=100,000 tissue pixels), confirming the registration is sufficient to use eosin as a second input channel. Registration is pre-IF (potential alignment artifacts) and does not resolve cells that overlap in z.

<p align="center"><img src="results/figures/registration_overlay.png" alt="DAPI + hematoxylin registration overlay" width="600"></p>

Nuclei appear white where DAPI (green) and hematoxylin (magenta) overlap, with minimal color fringing in both tumor-nest and stroma crops. The H&E anchor is trustworthy at the resolution of the 380-gene panel analysis; sub-cellular alignment errors are not expected to materially affect cell-level transcript assignment.

### Results: all methods vs. H&E morphology anchor

| Method | n_matched | Median corr | ARI (H&E anchor) | ARI (10x anchor) | Δ (H&E − 10x) |
| --- | ---: | ---: | ---: | ---: | ---: |
| **Whole-cell NN** | | | | | |
| Cellpose cyto3 (DAPI+density) | 19,625 | 0.962 | 0.807 | 0.569 | +0.238 |
| Mesmer WC (DAPI+density) | 20,747 | 0.977 | 0.776 | 0.617 | +0.159 |
| Cellpose cyto3 (DAPI+eosin) | 16,467 | 0.958 | 0.709 | 0.598 | +0.111 |
| Cellpose cyto3 (DAPI only) | 19,146 | 0.905 | 0.675 | 0.540 | +0.135 |
| **Nuclear-only** | | | | | |
| Mesmer | 20,951 | 0.949 | 0.727 | 0.557 | +0.170 |
| StarDist | 20,859 | 0.888 | 0.669 | 0.545 | +0.124 |
| CellPose | 18,624 | 0.873 | 0.701 | 0.547 | +0.154 |
| 10x Ranger | 20,641 | 0.870 | 0.689 | 0.504 | +0.185 |
| **Geometric expansion** | | | | | |
| Expansion 10µm (SD) | 20,895 | 0.930 | 0.728 | 0.631 | +0.097 |
| Expansion 20µm (SD) | 20,634 | 0.928 | 0.714 | 0.644 | +0.070 |
| Expansion 20µm (M) | 20,740 | 0.944 | 0.636 | 0.680 | −0.044 |
| Expansion 10µm (M) | 20,954 | 0.945 | 0.616 | 0.670 | −0.054 |
| Expansion 10µm (10x) | 20,687 | 0.928 | 0.621 | 0.691 | −0.070 |
| Expansion 20µm (CP) | 18,508 | 0.927 | 0.626 | 0.603 | +0.023 |
| Expansion 10µm (CP) | 18,709 | 0.930 | 0.552 | 0.604 | −0.052 |
| Expansion 20µm (10x) | 20,469 | 0.926 | 0.590 | 0.710 | −0.120 |
| **Watershed** | | | | | |
| Watershed (10x) | 19,960 | 0.933 | 0.711 | 0.664 | +0.047 |
| Watershed (M) | 20,172 | 0.949 | 0.658 | 0.674 | −0.016 |
| Watershed (SD) | 20,139 | 0.935 | 0.598 | 0.656 | −0.058 |
| **Voronoi** | | | | | |
| Voronoi (SD) | 20,859 | 0.914 | 0.664 | 0.584 | +0.080 |
| Voronoi (M) | 20,951 | 0.927 | 0.605 | 0.686 | −0.081 |
| Voronoi (10x) | 20,641 | 0.915 | 0.632 | 0.592 | +0.040 |
| Voronoi (CP) | 18,624 | 0.916 | 0.542 | 0.630 | −0.088 |
| **10x native** | 20,579 | 0.927 | **0.590** | 1.000 | −0.410 |
| **Baysor** | | | | | |
| Baysor (M prior 1.0) | 20,937 | 0.945 | 0.698 | 0.518 | +0.180 |
| Baysor (10x prior 1.0) | 20,645 | 0.911 | 0.681 | 0.530 | +0.151 |
| Baysor (SD prior 1.0) | 20,905 | 0.917 | 0.646 | 0.498 | +0.148 |
| Baysor (CP prior 1.0) | 19,465 | 0.912 | 0.654 | 0.501 | +0.153 |
| Baysor (CP prior 0.8) | 19,255 | 0.897 | 0.590 | 0.488 | +0.102 |
| Baysor (CP prior 0.5) | 14,644 | 0.843 | 0.473 | 0.349 | +0.124 |
| Baysor (CP prior 0.2) | 10,567 | 0.847 | 0.415 | 0.318 | +0.097 |
| Baysor (no prior) | 9,998 | 0.839 | 0.471 | 0.305 | +0.166 |

The H&E anchor reveals a strikingly different method ranking than the 10x native anchor. Three patterns stand out.

**Whole-cell NN methods ranked at the top here but not on B1.** Cellpose cyto3 DAPI+density scores ARI 0.807 — highest of all — compared to 0.569 against 10x native. Mesmer WC density jumps from 0.617 to 0.776. These methods were trained to segment whole cells against cytoplasm signal; against an anchor that is itself eosin-derived, they perform far better than their B1 scores suggested. Their apparently weak B1 performance reflects that 10x native's proprietary expansion does not follow cytoplasm boundaries.

**10x native scores only 0.590 against the H&E anchor.** This is the lowest among all expansion methods and below every Voronoi variant. The platform's own segmentation agrees moderately at best with what real tissue morphology looks like — its DAPI-based expansion assigns transcripts in a pattern that diverges substantially from cytoplasm extent. Methods designed to agree with 10x native (all B1 analysis) are therefore not necessarily capturing real cell morphology.

**Nuclear-only detectors score higher against H&E (0.67–0.73) than against 10x native (0.50–0.56).** Mesmer nuclear alone achieves ARI 0.727 against the H&E anchor. This makes sense: Mesmer's nucleus-only masks are closer in extent to Mesmer WC eosin masks than 10x native's arbitrarily expanded cells are. Accurate nuclear detection partially recovers morphological boundaries even without cytoplasm segmentation.

**Voronoi (CP) and Voronoi (M) rank lower on H&E than on 10x native** (CP: 0.542 vs. 0.630; M: 0.605 vs. 0.686). Voronoi tessellation assigns transcripts to the nearest nucleus centroid regardless of cytoplasm extent, so it diverges from eosin-defined boundaries whenever cell shape is asymmetric or cells are touching.

The Δ column (H&E − 10x ARI) isolates methods that produce genuinely morphology-grounded segmentation from those that reproduce DAPI-expansion logic: positive Δ means the method agrees better with tissue morphology than with the platform standard; negative Δ (Voronoi CP, Voronoi M, 10x native itself) means the method is closer to the platform's arbitrary expansion than to real cell shape.

---

## Benchmark 3: scRNA-seq Reference Anchor

The scRNA-seq anchor is grounded in what single cells "should" look like transcriptionally, independent of any spatial segmentation. PCA is fit on companion scRNA-seq from the same tissue blocks (GSE243275, 3' chemistry, 7,329 cells after QC) subsetted to the 374 Xenium panel genes present in both datasets. Each segmentation method's cells — including **10x native**, which is scored as a method here, not as an anchor — are projected into this reference space (30 PCs, 62% variance explained), clustered via Leiden, and compared.

Methods that score poorly here are capturing transcripts in patterns that don't match what isolated single cells look like, regardless of how well they agree with 10x native.

### Projection and density

![Xenium cells projected into scRNA-seq reference UMAP](results/figures/ref_projection_umap.png)

![Density distortion vs. scRNA-seq reference](results/figures/ref_projection_density.png)

All methods project into the same regions of the scRNA-seq reference landscape, but Baysor without a prior collapses into a narrow band — consistent with its lower transcript-per-cell counts compressing the phenotypic range. Voronoi methods and 10x native show nearly identical density profiles. Baysor PSC=1.0 variants are intermediate.

### Clustering convergence in reference space

![Reference-space Leiden clustering](results/figures/ref_projection_ref_clustering.png)

![Reference-space clusters mapped to tissue](results/figures/ref_projection_spatial.png)

Leiden clustering in the shared reference PCA space (resolution 1.0) produces 12–14 clusters for Voronoi methods and 10x native, but 22–25 for Baysor variants. The spatial maps confirm that reference-space clusters are biologically coherent: tissue structures (ducts, stroma, immune infiltrate) are visible across methods despite the clustering being performed in an independent coordinate system.

| Method | Ref clusters | ARI (own-PCA vs. 10x) | ARI (ref-PCA vs. 10x) | ΔARI |
| --- | ---: | ---: | ---: | ---: |
| **10x native** | 14 | 1.000 | 0.494 | −0.506 |
| **Voronoi** | | | | |
| Voronoi (CP) | 13 | 0.568 | 0.466 | −0.102 |
| Voronoi (SD) | 14 | 0.555 | 0.452 | −0.103 |
| Voronoi (M) | 13 | 0.614 | 0.477 | −0.137 |
| Voronoi (10x) | 11 | 0.624 | 0.489 | −0.135 |
| **Geometric 10µm** | | | | |
| Expansion 10µm (CP) | 14 | 0.591 | 0.450 | −0.141 |
| Expansion 10µm (SD) | 13 | 0.595 | 0.510 | −0.085 |
| Expansion 10µm (M) | 12 | 0.653 | 0.511 | −0.142 |
| Expansion 10µm (10x) | 14 | 0.691 | 0.518 | −0.172 |
| **Geometric 20µm** | | | | |
| Expansion 20µm (CP) | 13 | 0.593 | 0.453 | −0.141 |
| Expansion 20µm (SD) | 14 | 0.606 | 0.453 | −0.154 |
| Expansion 20µm (M) | 13 | 0.663 | 0.496 | −0.167 |
| Expansion 20µm (10x) | 15 | 0.710 | 0.510 | −0.199 |
| **Watershed** | | | | |
| Watershed (10x) | 12 | 0.658 | 0.468 | −0.190 |
| Watershed (SD) | 12 | 0.606 | 0.488 | −0.118 |
| Watershed (M) | 13 | 0.651 | 0.477 | −0.174 |
| **Baysor PSC sweep** | | | | |
| Baysor (no prior) | 24 | 0.172 | 0.179 | +0.007 |
| Baysor (CP prior 0.2) | 24 | 0.140 | 0.192 | +0.052 |
| Baysor (CP prior 0.5) | 20 | 0.206 | 0.229 | +0.023 |
| Baysor (CP prior 0.8) | 23 | 0.297 | 0.283 | −0.015 |
| Baysor (CP prior 1.0) | 21 | 0.327 | 0.290 | −0.036 |
| **Baysor other detectors** | | | | |
| Baysor (SD prior 1.0) | 25 | 0.306 | 0.275 | −0.031 |
| Baysor (M prior 1.0) | 26 | 0.324 | 0.296 | −0.028 |
| Baysor (10x prior 1.0) | 22 | 0.343 | 0.298 | −0.045 |
| **Whole-cell NN** | | | | |
| Cellpose cyto3 (DAPI) | 15 | 0.512 | 0.443 | −0.069 |
| Cellpose cyto3 (DAPI+eosin) | 11 | 0.583 | 0.488 | −0.094 |
| Cellpose cyto3 (DAPI+density) | 12 | 0.553 | 0.458 | −0.095 |
| Mesmer WC (DAPI+eosin) | 12 | 0.566 | 0.471 | −0.095 |
| Mesmer WC (DAPI+density) | 12 | 0.595 | 0.478 | −0.117 |

![Method convergence: internal vs. reference PCA pairwise ARI](results/figures/convergence_delta_ari.png)

Even 10x native's own reference-space clustering only achieves ARI 0.522 against its internal clustering, and every method except prior-free Baysor moves away from the original 10x native cell types after projection. The pairwise convergence observed between methods in the reference space (mean off-diagonal ARI 0.570 → 0.631) therefore reflects a shared loss of resolution rather than convergence toward ground truth.

### Cell-level displacement

![Matched-cell displacement by method](results/figures/ref_projection_displacement.png)

![Matched-cell displacement by cell type](results/figures/ref_projection_displacement_by_ct.png)

| Method | Matched cells | Median displacement |
| --- | ---: | ---: |
| Expansion 20µm (10x) | 23,628 | 0.32 |
| Expansion 10µm (10x) | 23,628 | 0.37 |
| Watershed (10x) | 23,581 | 0.73 |
| Expansion 20µm (SD) | 23,596 | 0.75 |
| Expansion 10µm (SD) | 23,562 | 0.79 |
| Expansion 10µm (M) | 23,606 | 0.87 |
| Expansion 20µm (M) | 23,617 | 0.89 |
| Voronoi (10x) | 23,626 | 0.92 |
| Expansion 10µm (CP) | 23,413 | 1.00 |
| Watershed (SD) | 23,535 | 1.02 |
| Expansion 20µm (CP) | 23,516 | 1.06 |
| Voronoi (SD) | 23,515 | 1.08 |
| Watershed (M) | 23,563 | 1.11 |
| Voronoi (M) | 23,580 | 1.18 |
| Voronoi (CP) | 23,267 | 1.41 |
| Mesmer WC (DAPI+eosin) | 23,603 | 1.71 |
| Mesmer WC (DAPI+density) | 23,574 | 1.89 |
| Baysor (M prior 1.0) | 23,628 | 2.01 |
| Cellpose cyto3 (DAPI+density) | 23,412 | 2.15 |
| Cellpose cyto3 (DAPI+eosin) | 22,632 | 2.24 |
| Baysor (CP prior 1.0) | 23,626 | 2.27 |
| Baysor (10x prior 1.0) | 23,629 | 2.29 |
| Baysor (CP prior 0.8) | 23,626 | 2.41 |
| Baysor (SD prior 1.0) | 23,628 | 2.41 |
| Cellpose cyto3 (DAPI) | 23,340 | 2.92 |
| Baysor (CP prior 0.2) | 22,487 | 2.93 |
| Baysor (no prior) | 22,344 | 2.98 |
| Baysor (CP prior 0.5) | 23,413 | 3.14 |

For each matched cell pair (nearest centroid, <15 µm), displacement measures how far the same cell moves in the 30-PC reference space depending on which segmentation method assigned its transcripts. Geometric expansion and watershed methods with 10x Ranger detection displace cells as little as 0.32–0.73 units — using the same nuclear seeds as 10x native, only the expansion algorithm differs. Voronoi methods sit at 0.92–1.41. Whole-cell NN methods (Mesmer WC) achieve 1.71–1.89, while Cellpose cyto3 splits sharply by input: DAPI+eosin/density at 2.15–2.24 and DAPI-only at 2.92 (as high as Baysor no-prior). Baysor PSC ≥ 0.8 displaces 2.01–2.41, and Baysor PSC ≤ 0.5 sits at 2.93–3.14.

### Population centroid distances

Two complementary measures quantify how far each method's cell populations sit from the scRNA-seq reference in the shared 30-PC space: cluster centroid distance (Leiden-resolution-dependent) and cell type centroid distance (resolution-independent, using `sc.tl.score_genes` on canonical marker panels).

<p align="center"><img src="results/figures/ref_projection_ari.png" alt="Pairwise ARI in reference PCA space" width="600"></p>

**Cell type centroid distance to scRNA-seq reference**

| Method | Mean | Median | Max |
| --- | ---: | ---: | ---: |
| **Nuclear-seeded expansion (Voronoi / Geometric / Watershed)** | | | |
| Watershed (M) | 3.47 | 3.45 | 4.90 |
| Expansion 10µm (M) | 3.47 | 3.45 | 4.91 |
| Expansion 20µm (M) | 3.47 | 3.48 | 4.90 |
| Expansion 10µm (CP) | 3.47 | 3.38 | 4.88 |
| Voronoi (M) | 3.48 | 3.48 | 4.87 |
| Expansion 20µm (CP) | 3.49 | 3.46 | 4.87 |
| Expansion 20µm (10x) | 3.50 | 3.35 | 4.97 |
| Voronoi (10x) | 3.50 | 3.42 | 4.93 |
| Expansion 10µm (10x) | 3.50 | 3.38 | 4.98 |
| **10x native** | **3.50** | **3.37** | **4.97** |
| Watershed (10x) | 3.51 | 3.38 | 4.98 |
| Voronoi (CP) | 3.53 | 3.53 | 4.86 |
| Voronoi (SD) | 3.54 | 3.42 | 4.97 |
| Expansion 10µm (SD) | 3.56 | 3.43 | 5.04 |
| Expansion 20µm (SD) | 3.55 | 3.45 | 5.02 |
| Watershed (SD) | 3.58 | 3.45 | 5.03 |
| **Whole-cell NN** | | | |
| Mesmer WC (DAPI+eosin) | 3.68 | 3.46 | 5.21 |
| Mesmer WC (DAPI+density) | 3.72 | 3.49 | 5.27 |
| Cellpose cyto3 (DAPI+eosin) | 3.84 | 3.58 | 5.37 |
| Cellpose cyto3 (DAPI+density) | 3.87 | 3.63 | 5.39 |
| Cellpose cyto3 (DAPI) | 4.30 | 4.02 | 5.90 |
| **Baysor (PSC ≥ 0.8)** | | | |
| Baysor (M prior 1.0) | 4.22 | 3.95 | 5.58 |
| Baysor (CP prior 1.0) | 4.28 | 4.01 | 5.65 |
| Baysor (CP prior 0.8) | 4.29 | 4.03 | 5.69 |
| Baysor (10x prior 1.0) | 4.33 | 4.06 | 5.74 |
| Baysor (SD prior 1.0) | 4.38 | 4.10 | 5.81 |
| **Baysor (PSC ≤ 0.5) / no prior** | | | |
| Baysor (CP prior 0.5) | 4.64 | 4.28 | 6.35 |
| Baysor (no prior) | 4.65 | 4.33 | 6.33 |
| Baysor (CP prior 0.2) | 4.61 | 4.31 | 6.28 |

All nuclear-seeded expanded methods — Voronoi, geometric, and watershed — land in a tight band (3.47–3.58), nearly indistinguishable from 10x native (3.50). The expansion strategy itself does not move methods relative to the scRNA-seq reference once transcripts are fully captured. Whole-cell NN methods show a two-level split: Mesmer WC is closest to the expanded-nuclear range (3.68–3.72), while Cellpose cyto3 eosin/density sits slightly further (3.84–3.87), and Cellpose cyto3 DAPI-only falls into the Baysor-prior range (4.30) — consistent with DAPI-only whole-cell boundaries failing to capture cytoplasm signal reliably. Baysor PSC ≥ 0.8 sits at 4.22–4.38 and Baysor PSC ≤ 0.5 at 4.61–4.65.

**Cluster centroid distance to scRNA-seq** (swept across Leiden resolutions):

| Method | 0.3 | 0.5 | 0.7 | 1.0 | 1.5 | 2.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 10x native | 4.09 | 4.11 | 4.09 | 4.06 | 4.00 | 4.01 |
| Voronoi (CP) | 4.05 | 4.05 | 4.00 | 3.91 | 3.81 | 3.78 |
| Voronoi (M) | 4.08 | 4.04 | 3.98 | 3.88 | 3.80 | 3.81 |
| Voronoi (10x) | 4.13 | 4.08 | 4.07 | 4.06 | 4.00 | 3.98 |
| Baysor | 4.57 | 4.58 | 4.58 | 4.57 | 4.57 | 4.58 |
| Baysor (CP prior 1.0) | 4.25 | 4.22 | 4.22 | 4.23 | 4.23 | 4.23 |
| Baysor (M prior 1.0) | 4.24 | 4.22 | 4.22 | 4.23 | 4.25 | 4.24 |

Voronoi methods improve modestly with finer resolution (−0.27 from 0.3→2.0), while Baysor methods are completely flat (Δ ≈ 0.00). Splitting Baysor from 11 to 37 clusters does not move centroids closer to the reference — they are fixed in the wrong locations in PCA space.

**Pairwise centroid distances** (Hungarian one-to-one matching, resolution 1.0):

|  | 10x native | V (CP) | V (SD) | V (M) | V (10x) | Baysor | B (CP) | B (SD) | B (M) | B (10x) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10x native | - | 0.96 | 0.40 | 0.92 | 0.46 | 2.77 | 2.39 | 2.17 | 1.90 | 2.20 |
| V (CP) | 0.96 | - | 0.88 | 0.66 | 0.93 | 3.35 | 2.98 | 2.77 | 2.47 | 2.80 |
| V (SD) | 0.40 | 0.88 | - | 0.99 | 0.46 | 2.61 | 2.25 | 2.03 | 1.74 | 2.09 |
| V (M) | 0.92 | 0.66 | 0.99 | - | 0.95 | 3.28 | 2.92 | 2.70 | 2.43 | 2.74 |
| V (10x) | 0.46 | 0.93 | 0.46 | 0.95 | - | 2.84 | 2.46 | 2.24 | 1.94 | 2.30 |
| Baysor | 2.77 | 3.35 | 2.61 | 3.28 | 2.84 | - | 0.63 | 0.57 | 0.72 | 0.60 |
| B (CP) | 2.39 | 2.98 | 2.25 | 2.92 | 2.46 | 0.63 | - | 0.41 | 0.46 | 0.43 |
| B (SD) | 2.17 | 2.77 | 2.03 | 2.70 | 2.24 | 0.57 | 0.41 | - | 0.52 | 0.46 |
| B (M) | 1.90 | 2.47 | 1.74 | 2.43 | 1.94 | 0.72 | 0.46 | 0.52 | - | 0.63 |
| B (10x) | 2.20 | 2.80 | 2.09 | 2.74 | 2.30 | 0.60 | 0.43 | 0.46 | 0.63 | - |

The pairwise matrix shows clear block structure: within-Voronoi distances 0.40–0.99, within-Baysor 0.41–0.72, cross-family 1.74–3.35. Voronoi and Baysor occupy distinct regions of the reference PCA space.

---

## Cross-Benchmark Synthesis

The three-anchor framework is designed to surface exactly the disagreements that a single-anchor analysis hides. The table below shows each method's score across all three benchmarks (B1: ARI vs. 10x native; B2: ARI vs. H&E morphology; B3: cell-type centroid distance to scRNA-seq reference, lower = better). Methods are sorted by B2 ARI.

| Method | B1 ARI (vs. 10x native) | B2 ARI (vs. H&E) | B3 centroid dist. | B2 − B1 |
| --- | ---: | ---: | ---: | ---: |
| **Whole-cell NN** | | | | |
| Cellpose cyto3 (DAPI+density) | 0.569 | 0.807 | 3.87 | +0.238 |
| Mesmer WC (DAPI+density) | 0.617 | 0.776 | 3.72 | +0.159 |
| Cellpose cyto3 (DAPI+eosin) | 0.598 | 0.709 | 3.84 | +0.111 |
| Cellpose cyto3 (DAPI only) | 0.540 | 0.675 | 4.30 | +0.135 |
| Mesmer WC (DAPI+eosin) | 0.590 | — | 3.68 | — |
| **Geometric 10µm** | | | | |
| Expansion 10µm (SD) | 0.631 | 0.728 | 3.56 | +0.097 |
| Expansion 10µm (M) | 0.670 | 0.616 | 3.47 | −0.054 |
| Expansion 10µm (10x) | 0.691 | 0.621 | 3.50 | −0.070 |
| Expansion 10µm (CP) | 0.604 | 0.552 | 3.47 | −0.052 |
| **Geometric 20µm** | | | | |
| Expansion 20µm (SD) | 0.644 | 0.714 | 3.55 | +0.070 |
| Expansion 20µm (M) | 0.680 | 0.636 | 3.47 | −0.044 |
| Expansion 20µm (10x) | 0.710 | 0.590 | 3.50 | −0.120 |
| Expansion 20µm (CP) | 0.603 | 0.626 | 3.49 | +0.023 |
| **Watershed** | | | | |
| Watershed (10x) | 0.664 | 0.711 | 3.51 | +0.047 |
| Watershed (M) | 0.674 | 0.658 | 3.47 | −0.016 |
| Watershed (SD) | 0.656 | 0.598 | 3.58 | −0.058 |
| **Voronoi** | | | | |
| Voronoi (M) | 0.686 | 0.605 | 3.48 | −0.081 |
| Voronoi (10x) | 0.592 | 0.632 | 3.50 | +0.040 |
| Voronoi (SD) | 0.584 | 0.664 | 3.54 | +0.080 |
| Voronoi (CP) | 0.630 | 0.542 | 3.53 | −0.088 |
| **10x native** | 1.000 | 0.590 | 3.50 | −0.410 |
| **Baysor PSC sweep** | | | | |
| Baysor (CP prior 1.0) | 0.501 | 0.654 | 4.28 | +0.153 |
| Baysor (CP prior 0.8) | 0.488 | 0.590 | 4.29 | +0.102 |
| Baysor (M prior 1.0) | 0.518 | 0.698 | 4.22 | +0.180 |
| Baysor (SD prior 1.0) | 0.498 | 0.646 | 4.38 | +0.148 |
| Baysor (10x prior 1.0) | 0.530 | 0.681 | 4.33 | +0.151 |
| Baysor (CP prior 0.5) | 0.349 | 0.473 | 4.64 | +0.124 |
| Baysor (CP prior 0.2) | 0.318 | 0.415 | 4.61 | +0.097 |
| Baysor (no prior) | 0.305 | 0.471 | 4.65 | +0.166 |

Four patterns emerge from the complete three-anchor view.

**Methods that score high on all three anchors don't exist.** Voronoi (M) has the highest B1 ARI (0.686) but a negative B2 delta (−0.081 vs. H&E). Cellpose cyto3 density tops B2 (0.807) but is lower on B1 (0.569) and B3 (3.87). No single method dominates across all three axes — which is the point. Each anchor is measuring something real and different.

**B2 − B1 sign divides the method space cleanly.** Methods with positive B2 − B1 (whole-cell NN, Baysor prior variants, StarDist-based geometric/watershed) agree more with H&E morphology than with 10x native. Methods with negative B2 − B1 (Voronoi CP, Voronoi M, geometric expansion with 10x Ranger, 10x native itself) agree more with the platform's DAPI-expansion logic. This is not a performance ranking — it is a characterization of what each method captures.

**10x native's B2 score (0.590) is lower than 26 of 33 tested methods.** The platform standard diverges substantially from tissue morphology. Its proprietary expansion does not follow cytoplasm boundaries; benchmarks anchored solely on 10x native measure agreement with this expansion logic, not with cell biology.

**B3 resolves the space into three tiers, not a continuous ranking.** All nuclear-seeded expanded methods — Voronoi (3.48–3.54), geometric (3.47–3.56), watershed (3.47–3.58) — land in essentially the same B3 band as 10x native (3.50). The expansion strategy does not move methods relative to the scRNA-seq reference once transcripts are fully captured; it is the presence of cytoplasm signal that matters. Whole-cell NN methods show a split: Mesmer WC (3.68–3.72) falls just above the expanded-nuclear cluster, while Cellpose cyto3 eosin/density (3.84–3.87) sits slightly further, and Cellpose cyto3 DAPI-only (4.30) falls into the Baysor-prior range — DAPI-only whole-cell boundaries fail to add cytoplasm signal. Baysor PSC ≥ 0.8 sits at 4.22–4.38 and Baysor PSC ≤ 0.5 at 4.61–4.65. Critically, Baysor's B3 gap relative to Voronoi (~0.8 units) persists regardless of which detector provides the prior.

**The negative marker analysis adds a fourth axis.** Whole-cell NN methods, which score mid-table on B1 and top B2, have by far the lowest cross-lineage violation rates (0.77–2.53%). This is a biologically independent signal: they produce fewer cells with impossible co-expression patterns, suggesting their learned boundaries cause less transcript leakage between adjacent lineages. High B1 ARI is weakly negatively correlated with low violation rate — methods that match 10x native most closely tend to have slightly higher violation rates, consistent with both reproducing the same boundary placement errors.

### Three-Way Anchor Correlation

To characterize how similar the three ranking systems are as evaluation instruments, Spearman rank correlations were computed across all 28 non-anchor, non-nuclear-only methods. B3 centroid distance is inverted so that lower distance = higher rank, consistent with B1 and B2 where higher = better.

| Anchor pair | Spearman ρ | p-value |
| --- | ---: | ---: |
| B1 vs. B2 | 0.139 | 0.49 |
| B1 vs. B3 | 0.832 | <0.0001 |
| B2 vs. B3 | 0.006 | 0.98 |

B1 and B3 are strongly correlated: methods that best reproduce 10x native clustering also project closest to the scRNA-seq reference cell-type centroids. Both anchors reward high transcript capture per cell and accurate boundary placement — the signals they encode are largely redundant. B2 is orthogonal to both. H&E morphology ranks methods by how well their boundaries follow visible cytoplasm and membrane signal, a property that is independent of whether those cells cluster like dissociated single cells. A method that draws accurate boundaries around complex morphology (whole-cell NN, StarDist-based expansion) can have mediocre B1 and mid-range B3; conversely, a method that reproduces 10x native clustering well can have a B2 far below the morphologically-accurate methods.

The top-5 methods per anchor share almost no overlap:

| Rank | B1 top 5 | B2 top 5 | B3 top 5 |
| ---: | --- | --- | --- |
| 1 | 10x native — 1.000 | Cellpose cyto3 density — 0.807 | Expansion 10µm (CP) — 3.47 |
| 2 | Expansion 20µm (10x) — 0.710 | Mesmer WC density — 0.776 | Watershed (M) — 3.47 |
| 3 | Expansion 10µm (10x) — 0.691 | Expansion 10µm (SD) — 0.728 | Expansion 10µm (M) — 3.47 |
| 4 | Voronoi (M) — 0.686 | Expansion 20µm (SD) — 0.714 | Expansion 20µm (M) — 3.47 |
| 5 | Expansion 20µm (M) — 0.680 | Watershed (10x) — 0.711 | Voronoi (M) — 3.48 |

Across the full top-10 lists, only one method appears in all three: Watershed (10x) (B1=0.664, B2=0.711, B3=3.51). It scores in the upper quarter on each anchor without topping any of them — a generalist that avoids the failure modes each individual anchor exposes. The three anchors are genuinely measuring different things, and their near-zero pairwise correlation with B2 confirms that H&E morphology provides orthogonal signal that neither the platform comparison nor the scRNA-seq reference captures.

### Full-Waterfall Promotion

Stages 7–16 (resolution sweep, cell type annotation, composition, spatial disagreement, phenotypic distortion, negative marker analysis) are computationally expensive and analytically deep. Three criteria determine which methods are promoted beyond the triage benchmarks:

1. Not nuclear-only (≥ 60% transcript capture required; nuclear-only methods recover 35–52%).
2. Pass the negative marker screen (violation rate ≤ 4%; Voronoi CP at 6.45% is borderline but included for comparison).
3. Represent a distinct region of the B1×B2×B3 space — not just a near-duplicate of an already-promoted method.

| Method | B1 | B2 | Neg. rate | Status |
| --- | ---: | ---: | ---: | --- |
| 10x native | 1.000 | 0.590 | 4.50% | Full waterfall (baseline) |
| Voronoi (CP/SD/M/10x) | 0.58–0.69 | 0.54–0.66 | 4.49–6.45% | Full waterfall (B1 leaders) |
| Baysor PSC sweep (CP/SD/M/10x prior) | 0.30–0.53 | 0.47–0.70 | 1.79–5.01% | Full waterfall (density-adaptive) |
| **Watershed (10x)** | **0.664** | **0.711** | **3.30%** | **Promoted — only method top-10 on all three anchors** |
| **Cellpose cyto3 (DAPI+density)** | **0.569** | **0.807** | **1.84%** | **Promoted — B2 leader, lowest neg. rate among WC NN** |
| **Mesmer WC (DAPI+density)** | **0.617** | **0.776** | **2.17%** | **Promoted — second B2, adds Mesmer WC boundary style** |
| Geometric ×8, Watershed SD/M | 0.60–0.71 | 0.60–0.73 | 2.4–4.9% | Triage only — near-duplicates of promoted methods on B3 |
| Cellpose cyto3 (DAPI/DAPI+eosin) | 0.54–0.60 | 0.67–0.71 | 0.77–1.64% | Triage only — dominated by DAPI+density on all anchors |

The three newly promoted methods cover the distinct regions of the B2×B3 plane that the original waterfall methods do not: Watershed (10x) is the only method strong on both the platform-like B1 anchor and the morphology B2 anchor simultaneously; Cellpose cyto3 DAPI+density and Mesmer WC DAPI+density sample the high-B2, mid-B1 space where whole-cell NN methods with strong morphological priors sit. The geometric expansion methods (×8) cluster tightly with the already-promoted Voronoi methods on B3 and offer no new information about the reference PCA landscape; they are retained in the synthesis tables but not run through stages 7–16.

---

## Reference-Free & Intrinsic Analyses

These analyses do not depend on the choice of reference anchor. Cell type annotation is derived from within each method's own clustering; the negative marker analysis requires no reference at all. Results here are comparable across anchors and serve as a check on findings from Benchmarks 1–3.

### Cell type annotation

Cell types are annotated on the 10x native segmentation. Leiden clustering (resolution 1.0) partitions the 10x native cells into 15 clusters, then differential expression identifies each cluster's distinguishing genes, matched to canonical breast tissue markers. The raw Xenium output carries no cell type labels — only coordinates and transcript counts.

| Cluster | Cells | Annotation |
| ---: | ---: | --- |
| 0 | 2,799 | Luminal epithelial |
| 1 | 2,926 | Luminal epithelial |
| 2 | 335 | Macrophages |
| 3 | 902 | Luminal epithelial |
| 4 | 1,025 | Myoepithelial |
| 5 | 2,861 | T cells |
| 6 | 758 | B cells |
| 7 | 2,612 | Macrophages |
| 8 | 1,921 | Luminal epithelial |
| 9 | 3,314 | CAFs |
| 10 | 870 | Smooth muscle |
| 11 | 1,508 | Endothelial |
| 12 | 266 | Plasma cells |
| 13 | 1,333 | CAFs |
| 14 | 199 | Adipocytes |

The dotplot below shows which canonical markers are expressed in each cluster, confirming the assignments. Clusters 0, 1, 3, and 8 all annotate as luminal epithelial but are distinguished by different marker profiles: cluster 0 is ESR1/FOXA1-dominant (ER+ hormone-responsive), cluster 1 is PGR/MUC1-dominant, cluster 3 expresses stromal-adjacent markers (NNMT, LUM), and cluster 8 is TACSTD2/KRT7/STC2-dominant (proliferative/stress-response). Clusters 2 and 7 are both macrophage populations: cluster 2 (335 cells) expresses FCGR3A and HAVCR2 (non-classical/M2-like), while cluster 7 (2,612 cells) expresses CD14 and AIF1 (classical monocyte-derived). Clusters 9 and 13 are both CAFs: cluster 9 expresses SFRP4 and FBLN1 (matrix-producing), while cluster 13 expresses POSTN and CTHRC1 (myofibroblastic).

**Annotated UMAP**

![10x native UMAP — cell type annotation and Leiden cluster IDs](results/figures/umap_annotated.png)

**Canonical marker expression by Leiden cluster**

![Canonical marker expression by Leiden cluster](results/figures/annotation_dotplot_flipped.png)

### Luminal epithelial subtypes

The tissue sample is an invasive ductal carcinoma (IDC) of the breast, with the ROI spanning several tumor nests (DCIS and invasive) surrounded by stroma and immune infiltrate. At Leiden resolution 1.0, four of the 15 clusters (0, 1, 3, 8) annotate as luminal epithelial, together accounting for 8,548 cells — 36% of the ROI. Wilcoxon DE within this subpopulation reveals that the four clusters are not redundant: each has a distinct marker profile corresponding to known luminal biology.

| Subcluster | Cells | Top markers | Interpretation |
| ---: | ---: | --- | --- |
| 0 | 2,799 | FLNB, DHRS2, MLPH, ESR1, PDZK1, AGR3 | ER+ hormone-responsive (luminal A-like) |
| 1 | 2,926 | MYBPC1, PGR, CLIC6, FASN, MUC1, ELOVL5, SCD | Secretory luminal with lipid metabolism |
| 3 | 902 | NNMT, ENAH, LUM, POSTN, AQP3 | Stromal-adjacent; EMT-like or boundary contamination |
| 8 | 1,921 | EGR1, STC2, KRT7, CCND1, TACSTD2, TCIM | Stress-response / proliferative |

![Top markers per luminal epithelial subcluster](results/figures/lum_ep_subtypes/lum_ep_subcluster_dotplot.png)

![Luminal epithelial subcluster heatmap](results/figures/lum_ep_subtypes/lum_ep_subcluster_heatmap.png)

Clusters 0 and 1 map to recognized luminal subtypes in breast cancer: cluster 0 expresses the canonical ER+ program (ESR1, MLPH, PDZK1, AGR3) while cluster 1 emphasizes progesterone signaling and lipid biosynthesis (PGR, FASN, ELOVL5, SCD). Cluster 8 is dominated by immediate-early genes (EGR1 logFC +3.4, STC2 logFC +4.2) and CCND1, consistent with an actively cycling or stress-responsive state. Cluster 3 is the most ambiguous — its top markers include NNMT (metabolic stress), plus LUM and POSTN, which are canonical CAF markers. This could reflect genuine EMT in a subset of tumor cells, or it could indicate that 10x native's expansion algorithm is pulling stromal transcripts into epithelial cell boundaries in densely packed regions.

**Resolution sweep**

![Leiden resolution sweep - luminal epithelial splitting](results/figures/lum_ep_subtypes/lum_ep_resolution_sweep.png)

A Leiden resolution sweep from 0.1 to 2.0 reveals that the four subclusters do not emerge simultaneously. Cluster 8 (EGR1/STC2/KRT7) splits off first at resolution 0.3. Clusters 0 and 1 separate from each other at resolution 0.6. Cluster 3 (NNMT/LUM/POSTN) is the last to emerge, not appearing as a distinct cluster until resolution 0.9. The late splitting and stromal marker expression together make cluster 3 the weakest of the four subtypes.

![Luminal epithelial subclusters in UMAP context](results/figures/lum_ep_subtypes/lum_ep_umap_context.png)

### Cluster 3: segmentation contamination test

Cluster 3's top markers (LUM, POSTN, NNMT) are canonical CAF genes that also define clusters 9 and 13. To test whether this stromal signal arises from segmentation boundaries pulling CAF transcripts into epithelial cells, the 1,277 cluster 3 cells are spatially matched across all segmentation methods and their mean expression compared. If the signal is contamination, methods with broader boundaries should show a stronger stromal signature; if it is genuine EMT, the signal should be stable.

| Method | Matched | LUM | POSTN | NNMT | GATA3 | ESR1 | PGR | Stromal/luminal ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10x native | 1,277 | 0.733 | 0.258 | 0.244 | 0.309 | 0.036 | 0.043 | 3.18 |
| Voronoi (CP) | 1,276 | 1.081 | 0.441 | 0.384 | 0.389 | 0.057 | 0.070 | 3.69 |
| Voronoi (SD) | 1,277 | 0.790 | 0.301 | 0.263 | 0.323 | 0.044 | 0.049 | 3.26 |
| Voronoi (M) | 1,277 | 0.925 | 0.357 | 0.319 | 0.362 | 0.047 | 0.053 | 3.47 |
| Voronoi (10x) | 1,277 | 0.803 | 0.295 | 0.268 | 0.323 | 0.038 | 0.052 | 3.31 |
| Baysor (CP prior 1.0) | 1,277 | 0.469 | 0.183 | 0.155 | 0.211 | 0.024 | 0.023 | 3.13 |
| Baysor (M prior 1.0) | 1,277 | 0.481 | 0.184 | 0.162 | 0.231 | 0.026 | 0.028 | 2.90 |
| Baysor (no prior) | 1,275 | 0.544 | 0.211 | 0.217 | 0.147 | 0.019 | 0.019 | 5.25 |

Both stromal and luminal markers vary in absolute terms across methods — this reflects differences in total transcript capture per cell, not marker-specific effects. The stromal-to-luminal ratio isolates the relative contribution. Voronoi methods have higher ratios (3.26–3.69) than 10x native (3.18), consistent with broader nearest-centroid boundaries pulling in proportionally more stromal transcripts. Baysor prior methods have lower ratios (2.90–3.13), consistent with density-adaptive expansion capturing less cross-boundary contamination. The pattern matches the negative marker analysis: methods with lower violation rates also show less stromal enrichment in cluster 3. The three promoted methods are not yet in this table pending a dedicated analysis pass. Based on their Tier 1 violation rates, Watershed (10x) at 4.88% is expected to land near the Voronoi (10x) range (ratio ~3.3), while Cellpose cyto3 density (1.84%) and Mesmer WC density (2.17%) — below all Voronoi and geometric methods — are expected to fall between the Baysor prior and Voronoi bands.

### Phenotypic landscape distortion

**All methods in shared PCA/UMAP space**

![All methods in shared PCA/UMAP space](results/figures/manifold_shared_umap.png)

**Density distortion vs. 10x native**

![Phenotypic landscape distortion vs. 10x native](results/figures/manifold_distortion.png)

All methods are projected into a shared PCA space fit on 10x native (30 PCs, 55% variance explained) and embedded in a joint UMAP. Density ratio maps (log₂ method/10x) show which phenotypic regions each method enriches or depletes. Nuclear methods show depleted regions in high-density luminal epithelial areas, consistent with missed cytoplasmic transcripts pulling cells toward lower-expression PCA states. Voronoi methods track 10x native closely. Baysor shows enrichment in a distinct region corresponding to its finer resolution of macrophage and stromal subtypes.

### Marker gene recovery

<p align="center"><img src="results/figures/marker_recovery.png" alt="Marker gene recovery relative to 10x native" width="600"></p>

Using 10x-native cell-type annotations as ground truth, nuclear methods recover 75–92% of cytoplasmic marker expression relative to 10x native, with the largest deficits for extranuclear markers like MUC1, SERPINA3, and LYZ. Voronoi methods recover near-100% across all cell types. Baysor recovers macrophage markers (LYZ, CD14) at or above 10x-native levels while showing slightly reduced T cell marker (CD3E) recovery.

### Population-level convergence (pseudobulk)

<p align="center"><img src="results/figures/pseudobulk_correlation.png" alt="Pseudobulk correlation vs. 10x native" width="600"></p>

| Method | Per-cell-type pseudobulk r (range) | Aggregate r | Single-cell ARI |
| --- | --- | --- | --- |
| CellPose | 0.87-0.98 | 0.970 | 0.547 |
| StarDist | 0.88-0.99 | 0.975 | 0.545 |
| Mesmer | 0.92-0.99 | 0.983 | 0.557 |
| Voronoi (CP) | 0.98-1.00 | 0.9999 | 0.630 |
| Voronoi (SD) | 0.98-1.00 | 0.9999 | 0.584 |
| Voronoi (M) | 0.98-1.00 | 0.9999 | 0.686 |
| Watershed (10x) | 0.97-1.00 | 0.9999 | 0.664 |
| Baysor | 0.94-1.00 | 0.999 | 0.305 |
| Cellpose cyto3 (DAPI+density) | 0.91-0.98 | 0.983 | 0.569 |
| Mesmer WC (DAPI+density) | 0.93-0.99 | 0.989 | 0.617 |

Pseudobulk is computed within each of 10 annotated cell types, testing whether each method's cell-type compartments recover the same expression programs as 10x native. Despite its low single-cell ARI of 0.305, Baysor is competitive with nuclear methods at the cell-type level — its aggregate r of 0.999 sits above CellPose (0.970) and StarDist (0.975). Nuclear methods show reduced pseudobulk r (0.97–0.98) because missing cytoplasmic transcripts suppress marker signal systematically across all cells of a type. Voronoi methods achieve both high single-cell ARI and near-perfect pseudobulk agreement. Watershed (10x) matches the Voronoi tier exactly (aggregate r 0.9999, range 0.97–1.00), expected given its shared 10x Ranger nuclear seeds. Cellpose cyto3 density and Mesmer WC density fall at 0.983 and 0.989 respectively — above the nuclear-only methods but below the Voronoi tier — with the per-cell-type spread widest at plasma cells and adipocytes, the two rarest and smallest populations. The aggregate r compresses this spread; at the cell-type level, rare populations with fewer cells per pseudobulk sample drive the per-method range.

### Negative marker analysis

Every comparison above uses an anchor segmentation as reference. The negative marker analysis provides a reference-free alternative: instead of asking "how well does method X agree with anchor Y?", it asks "how often does method X produce cells with biologically impossible co-expression?" A cell co-expressing CD3E (T cell) and GATA3 (luminal epithelial) almost certainly contains transcripts from two adjacent cells that were merged by the segmentation algorithm. The rate of these violations reflects boundary quality without privileging any single method as ground truth.

Eleven Tier 1 pairs are defined from the 380-gene panel, each pairing a lineage-specific marker from one cell type with a lineage-specific marker from a developmentally unrelated cell type. Two additional Tier 2 pairs are included with lower confidence. A cell is flagged as a violation if it expresses both markers at ≥2 raw counts.

| Pair | Lineage A | Lineage B |
| --- | --- | --- |
| CD3E + GATA3 | T cell | Luminal epithelial |
| CD3E + ESR1 | T cell | Luminal epithelial |
| TRAC + KRT14 | T cell | Myoepithelial |
| MS4A1 + GATA3 | B cell | Luminal epithelial |
| MS4A1 + KRT14 | B cell | Myoepithelial |
| PECAM1 + ESR1 | Endothelial | Luminal epithelial |
| VWF + CD3E | Endothelial | T cell |
| CD14 + ESR1 | Macrophage | Luminal epithelial |
| LYZ + GATA3 | Macrophage | Luminal epithelial |
| ADIPOQ + CD3E | Adipocyte | T cell |
| ADIPOQ + GATA3 | Adipocyte | Luminal epithelial |

**Violation rates across methods**

![Negative marker analysis](results/figures/negative_marker.png)

| Method | Cells | Tier 1 violations | Violation rate | Violations per 1000 tx |
| --- | ---: | ---: | ---: | ---: |
| **10x native** | | | | |
| 10x native | 23,629 | 1,064 | 4.50% | 0.36 |
| **Voronoi expansion** | | | | |
| Voronoi (CP) | 20,166 | 1,300 | 6.45% | 0.43 |
| Voronoi (SD) | 24,743 | 1,111 | 4.49% | 0.37 |
| Voronoi (M) | 21,697 | 1,189 | 5.48% | 0.39 |
| Voronoi (10x) | 23,622 | 1,153 | 4.88% | 0.38 |
| **Geometric expansion 10µm** | | | | |
| Expansion 10µm (CP) | 20,166 | 1,147 | 5.69% | 0.40 |
| Expansion 10µm (SD) | 24,745 | 978 | 3.95% | 0.34 |
| Expansion 10µm (M) | 21,697 | 1,101 | 5.07% | 0.38 |
| Expansion 10µm (10x) | 23,624 | 1,008 | 4.27% | 0.35 |
| **Geometric expansion 20µm** | | | | |
| Expansion 20µm (CP) | 20,166 | 1,249 | 6.19% | 0.42 |
| Expansion 20µm (SD) | 24,745 | 1,052 | 4.25% | 0.36 |
| Expansion 20µm (M) | 21,697 | 1,153 | 5.31% | 0.39 |
| Expansion 20µm (10x) | 23,624 | 1,069 | 4.53% | 0.36 |
| **Watershed expansion** | | | | |
| Watershed (10x) | 23,624 | 1,154 | 4.88% | 0.40 |
| Watershed (SD) | 24,745 | 1,126 | 4.55% | 0.40 |
| Watershed (M) | 21,697 | 1,207 | 5.56% | 0.41 |
| **Baysor PSC sweep (CP prior)** | | | | |
| Baysor (no prior) | 18,321 | 918 | 5.01% | 0.95 |
| Baysor (CP prior 0.2) | 19,061 | 898 | 4.71% | 0.89 |
| Baysor (CP prior 0.5) | 24,147 | 792 | 3.28% | 0.67 |
| Baysor (CP prior 0.8) | 29,771 | 664 | 2.23% | 0.33 |
| Baysor (CP prior 1.0) | 30,473 | 661 | 2.17% | 0.31 |
| **Baysor other detectors at PSC=1.0** | | | | |
| Baysor (SD prior 1.0) | 34,230 | 613 | 1.79% | 0.28 |
| Baysor (M prior 1.0) | 31,764 | 742 | 2.34% | 0.32 |
| Baysor (10x prior 1.0) | 33,113 | 631 | 1.91% | 0.29 |
| **Whole-cell NN** | | | | |
| Cellpose cyto3 (DAPI) | 21,641 | 167 | 0.77% | 0.14 |
| Cellpose cyto3 (DAPI+eosin) | 17,450 | 287 | 1.64% | 0.21 |
| Cellpose cyto3 (DAPI+density) | 20,782 | 382 | 1.84% | 0.24 |
| Mesmer WC (DAPI+eosin) | 21,841 | 552 | 2.53% | 0.28 |
| Mesmer WC (DAPI+density) | 21,493 | 467 | 2.17% | 0.25 |

Nuclear-only methods (CellPose, StarDist, Mesmer, 10x Ranger nuclei) are excluded because their low transcript capture (35–52%) makes violations trivially rare. The expanded analysis across all 29 methods reveals a clear family-level pattern. Whole-cell NN methods have by far the lowest violation rates — Cellpose cyto3 (DAPI) reaches just 0.77% (0.14 per 1000 tx), lower than any geometric or Voronoi method by a factor of 3–8. This is a qualitatively different result from Benchmark 1 where whole-cell NN methods scored only modestly on ARI (0.54–0.62). When boundaries are evaluated against biology rather than against 10x native, methods that learned real cell morphology outperform centroid-based expansion.

Among geometric and watershed methods, Voronoi (CP) has the highest raw violation rate (6.45%, 0.43 per 1000 tx) and geometric expansion at 20µm (CP) is nearly as high (6.19%). Watershed methods (4.5–5.6%) are similar to Voronoi/geometric despite using the DAPI gradient — the gradient provides limited discrimination between adjacent cells of different lineages. Watershed (10x), the method promoted to the full waterfall as the only entry in the top-10 on all three anchors, sits at 4.88%, confirming that DAPI-gradient expansion does not reduce cross-lineage boundary leakage compared to simpler centroid-based approaches. Geometric expansion using StarDist nuclei at 10µm is the best among centroid-based methods (3.95%), likely because StarDist finds the most cells and its nuclei are well-calibrated in size, leaving less inter-nuclear space to be incorrectly expanded into.

The Baysor PSC sweep shows a monotonic improvement with increasing prior strength. The transition between PSC=0.2 (4.71%) and PSC=0.5 (3.28%) marks where the prior begins to limit cross-lineage merging, with PSC=0.8 and 1.0 (2.17–2.23%) approaching the whole-cell NN range. Baysor (SD prior 1.0) achieves the lowest rate among all Baysor variants (1.79%), probably because StarDist detects more cells and its prior provides more precise nuclear anchoring. Baysor without a prior has a low raw rate (5.01%) but the highest transcript-normalized rate (0.95 per 1000 tx) because its tendency to merge cells produces very large segments that accumulate diverse transcripts.

The LYZ+GATA3 pair (macrophage + luminal epithelial) dominates violations across all methods, consistent with macrophages infiltrating tumor epithelium where segmentation boundaries are most ambiguous.

### Cross-lineage bleed: spatial structure

The negative marker analysis counts violations at the cell level. The cross-lineage bleed analysis maps these violations spatially to test whether they concentrate at lineage boundaries — the expected pattern if violations arise from segmentation mis-assignment rather than genuine co-expression. For each method, T cells are scored by their mean expression of CAF markers (LUM, SFRP4, FBLN1, CCDC80, THBS2, MMP2), and those with above-median bleed are plotted alongside CAF positions.

![Cross-lineage bleed across methods](results/figures/cross_lineage_bleed.png)

Under 10x native, Voronoi methods, Watershed (10x), and Cellpose cyto3 density, approximately 50% of T cells register above the high-bleed threshold and median T-cell-to-nearest-CAF distances span 15–33µm — consistent with boundary-mediated leakage in which T cells adjacent to CAF-rich stroma accumulate CAF marker signal regardless of expansion algorithm. Baysor (CP prior 1.0) reduces this modestly: 37.8% of T cells are high-bleed with a median nearest-CAF distance of 33.2µm, reflecting density-adaptive expansion that partially constrains cross-lineage capture without eliminating it.

Mesmer WC density is qualitatively different from every other method. Only 9.6% of its T cells score above the high-bleed threshold — a factor of four to five lower than any other method in the benchmark — and the median distance from all T cells to the nearest CAF is 131µm, versus 15–33µm for all other methods. At 131µm, the high-bleed T cells are not adjacent to CAF boundaries at all; the small fraction of violations that do occur are not a proximity artifact. The contrast with Cellpose cyto3 density (49.9% high-bleed, 33µm) is striking: both are whole-cell NN methods at similar transcript capture rates (54% vs 64%), but Mesmer's training on nuclear-plus-density signal produces boundaries that exclude the narrow zones between CAFs and immune cells that every centroid-based and density-gradient method treats as fillable space. This is the largest single-method spatial signal in the benchmark — not a marginal improvement but a qualitative step-change that only becomes visible at the per-cell spatial level.

![Cross-lineage bleed detail — Voronoi (CP)](results/figures/cross_lineage_bleed_detail.png)

The detail panel for Voronoi (CP) confirms the proximity dependence: T cells closer to CAFs have higher CAF-marker expression (center), and the distance distribution of high-bleed cells is shifted toward shorter distances compared to low-bleed cells (right). This spatial gradient is the expected signature of transcript leakage across segmentation boundaries and would not arise from genuine T cell biology.

---

## Repo layout

```text
segmentation-benchmark/
├── environment.yml          # conda env (CellPose, Scanpy, Squidpy, SpatialData, ...)
├── data/
│   ├── raw/                 # downloaded Xenium bundle (gitignored)
│   └── processed/           # cropped ROI + derived files (gitignored)
├── notebooks/
├── src/segbench/
│   ├── constants.py         # method metadata, cell-type annotations, negative marker pairs
│   ├── io.py                # load Xenium bundle, ROI cropping
│   ├── segmentation/        # per-method wrappers (CellPose, StarDist, Mesmer, Baysor)
│   ├── quantify.py          # transcript aggregation -> per-cell AnnData
│   ├── compare.py           # cross-method comparison metrics
│   ├── spatial.py           # spatial structure of disagreement
│   └── style.py             # shared matplotlib theme
├── scripts/                 # CLI entry points
├── results/{figures,tables}/
└── tests/
```

## Environment setup

This project uses three toolchains: a main conda env for CellPose + Scanpy/Squidpy/SpatialData, a separate env for StarDist (TensorFlow-based), and Julia for Baysor. Mesmer runs via Docker.

### 1. Main env

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
docker pull vanvalenlab/deepcell-applications:latest
```

The image bundles pretrained model weights and does not require a `DEEPCELL_ACCESS_TOKEN`. See [`scripts/run_mesmer.sh`](scripts/run_mesmer.sh).

### 4. Julia + Baysor

```bash
juliaup add 1.10
julia +1.10 -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git", rev="v0.7.1")); Pkg.build("Baysor")'
```

See [`scripts/run_baysor.sh`](scripts/run_baysor.sh).
