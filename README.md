# Segmentation Benchmarking on Xenium Spatial Transcriptomics Data

**Question:** Do nuclear-mask (CellPose, StarDist, Mesmer), Voronoi, and transcript-density (Baysor) segmentation methods transfer well to Xenium spatial transcriptomics, and does method choice meaningfully change downstream cell-type calls?

## Summary

Cluster labels are aligned via two algorithms before computing disagreement and Moran's I. Hungarian finds the optimal one-to-one assignment; when cluster counts differ, unmatched clusters are forced into poor pairings. Argmax maps each method's clusters to the 10x native cluster with plurality overlap, allowing many-to-one mapping. Matched pairs, median correlation, and ARI do not depend on cluster alignment.

| Comparison | Matched pairs | Median corr | ARI | Hungarian Disagree | Hungarian Moran's I | Argmax Disagree | Argmax Moran's I |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **Nuclear-only** | | | | | | | |
| 10x native vs. CellPose | 18,966 | 0.822 | 0.547 | 30.8% | 0.178 | 30.5% | 0.189 |
| 10x native vs. StarDist | 21,429 | 0.826 | 0.545 | 33.5% | 0.215 | 33.4% | 0.221 |
| 10x native vs. Mesmer | 20,595 | 0.879 | 0.557 | 27.9% | 0.090 | 27.4% | 0.098 |
| 10x native vs. 10x Ranger | 23,155 | 0.822 | 0.504 | 35.0% | 0.191 | 34.5% | 0.203 |
| **Voronoi** | | | | | | | |
| 10x native vs. Voronoi (CP) | 18,966 | 0.959 | 0.630 | 21.9% | 0.076 | 21.9% | 0.076 |
| 10x native vs. Voronoi (SD) | 21,428 | 0.959 | 0.584 | 31.9% | 0.194 | 27.7% | 0.229 |
| 10x native vs. Voronoi (M) | 20,595 | 0.964 | 0.686 | 18.8% | 0.161 | 18.8% | 0.161 |
| 10x native vs. Voronoi (10x) | 23,153 | - | 0.592 | 28.3% | 0.172 | 25.8% | 0.168 |
| **Baysor** | | | | | | | |
| 10x native vs. Baysor | 10,953 | 0.786 | 0.305 | 51.7% | 0.033 | 43.8% | 0.079 |
| 10x native vs. Baysor (CP prior 0.2) | 11,454 | 0.798 | 0.318 | 51.9% | 0.036 | 39.2% | 0.086 |
| 10x native vs. Baysor (CP prior 1.0) | 20,308 | 0.902 | 0.501 | 33.8% | 0.111 | 32.1% | 0.122 |
| 10x native vs. Baysor (SD prior 1.0) | 21,814 | 0.905 | 0.498 | 37.7% | 0.136 | 32.9% | 0.170 |
| 10x native vs. Baysor (M prior 1.0) | 21,148 | 0.924 | 0.518 | 32.3% | 0.115 | 30.7% | 0.119 |
| 10x native vs. Baysor (10x prior 1.0) | 22,910 | 0.914 | 0.530 | 34.7% | 0.208 | 33.1% | 0.204 |

*Matched pairs*: nearest-centroid matching. *Median corr*: per-pair Pearson correlation of log-normalised expression. *ARI*: Adjusted Rand Index (partition-based, independent of cluster alignment). *Disagreement*: fraction of matched cell pairs assigned to different clusters after alignment. *Moran's I*: spatial autocorrelation of the disagree flag. Voronoi (10x) correlation will be populated on next pipeline run (gene-name encoding mismatch has been fixed).

The results cleanly separate the contributions of nuclear detection quality and expansion strategy. Among Voronoi methods, Mesmer centroids produce the highest ARI (0.686) and lowest disagreement (18.8%), while 10x Ranger centroids - despite being purpose-built for Xenium - score lower (ARI 0.592). Among Baysor PSC=1.0 variants, the same detector ordering holds: Mesmer prior leads at ARI 0.518, followed by 10x Ranger (0.530), CellPose (0.501), and StarDist (0.498). Baysor (10x prior 1.0) achieves the highest ARI of any Baysor variant (0.530) and the most matched pairs (22,910), benefiting from 10x Ranger detecting almost as many nuclei as the 10x native reference.

Across both expansion strategies, Voronoi consistently outperforms Baysor PSC=1.0 on ARI (0.58-0.69 vs 0.50-0.53) - geometric nearest-centroid assignment agrees more with 10x native than density-adaptive expansion. However, Baysor prior variants have lower negative marker violation rates (0.31 per 1000 tx vs 0.37-0.43 for Voronoi), indicating that density-adaptive boundaries produce fewer cross-lineage contamination artifacts even when they disagree with the 10x reference.

Nuclear-only methods capture too few transcripts for meaningful downstream comparison and are excluded from figures past the recovery section.

<!-- Project 2 (label-transfer-benchmark): uses this project's segmented cells to evaluate scRNA-seq label-transfer reliability. Add link once repo is public. -->

## Dataset

**Xenium FFPE Human Breast (Custom Add-on Panel)**, Janesick et al. 2023, *Nature Communications* ([dataset page](https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-with-custom-add-on-panel-1-standard)). Invasive ductal carcinoma; matched scRNA-seq + Visium from the same tissue blocks: GEO [GSE243275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243275).

All analysis runs on a 2mm × 2mm ROI (~23,600 cells, ~3.4M transcripts, 380-gene panel) with a mix of tumor, stroma, and immune-infiltrated regions. See [`docs/dataset.md`](docs/dataset.md) for download and ROI details.

## Methods

| Method | Input | Notes |
| --- | --- | --- |
| **10x native** | provided | Xenium Ranger's full segmentation (nuclear detection + proprietary expansion); used as reference anchor. The expansion algorithm is closed-source. |
| **10x Ranger** | DAPI | The nuclear detection component of Xenium Ranger, extracted from `nucleus_boundaries.parquet` and rasterized into a label mask. Included alongside CellPose/StarDist/Mesmer as a fourth nuclear detector to test whether 10x's purpose-built detector outperforms general-purpose models. |
| **CellPose** | DAPI | CellPose 3.x `nuclei` model, CPU |
| **StarDist** | DAPI | `2D_versatile_fluo` model, separate `stardist` env |
| **Mesmer** | DAPI | DeepCell via Docker; image bundles model weights |
| **Voronoi (CP / SD / M / 10x)** | nuclear centroids | Nearest-centroid transcript assignment using CellPose, StarDist, Mesmer, or 10x Ranger centroids; 100% transcript capture by construction |
| **Baysor** | transcripts | Transcript-density EM (no prior), Julia 1.10, 4 tiles |
| **Baysor (prior 0.2 / 0.8 / 1.0)** | transcripts + nuclear masks | Baysor with `prior_segmentation_confidence` controlling the blend between density model and nuclear prior. At PSC 1.0, nuclear transcripts are hard-locked and only cytoplasmic transcripts use density-adaptive expansion. Tested with all four nuclear detectors at PSC 1.0 to isolate detector quality from expansion strategy. |

Nuclear-only methods (CellPose, StarDist, Mesmer, 10x Ranger) capture only 35-52% of transcripts and are included in the cell/transcript recovery section but excluded from downstream figures because their low transcript capture dominates any comparison. Cells are matched by nearest centroid across methods. Leiden clustering runs independently on each method's cells; cluster labels are aligned via Hungarian algorithm before computing ARI and disagreement rate.

---

## Cell and transcript recovery

### Nuclear detectors

| Detector | Cells | Median tx/cell | Transcript capture |
| --- | ---: | ---: | ---: |
| CellPose | 20,166 | 49 | 35.4% |
| StarDist | 24,745 | 45 | 40.8% |
| Mesmer | 21,697 | 70 | 51.8% |
| 10x Ranger | 23,624 | 45 | 38.0% |

![Nuclear mask size comparison across all four detectors](results/figures/nuclear_mask_sizes.png)

All four detectors operate on the same DAPI image but produce substantially different masks. Mesmer detects the largest nuclei (median ~45 µm², long tail to 200 µm²), capturing 51.8% of transcripts - nearly double CellPose's 35.4%. StarDist and 10x Ranger produce similar-sized masks but StarDist finds more nuclei (24,745 vs 23,624). 10x Ranger captures only 38% of transcripts despite detecting nearly as many cells as the 10x native whole-cell segmentation (23,624 vs 23,629), confirming that 10x native's 99% capture comes from its proprietary expansion, not from larger nuclei.

### Expansion methods

| Method | Cells | Median tx/cell | Transcript capture |
| --- | ---: | ---: | ---: |
| 10x native (Ranger + proprietary expansion) | 23,629 | 124 | 99.0% |
| Voronoi (CP) | 20,166 | 149 | 100% |
| Voronoi (SD) | 24,745 | 122 | 100% |
| Voronoi (M) | 21,697 | 142 | 100% |
| Voronoi (10x) | 23,622 | 128 | 100% |
| Baysor (no prior) | 18,321 | 53 | 98.6% |
| Baysor (CP prior 0.2) | 19,061 | 53 | 98.7% |
| Baysor (CP prior 0.8) | 29,771 | 67 | 99% |
| Baysor (CP prior 1.0) | 30,473 | 69 | 99% |
| Baysor (SD prior 1.0) | 34,230 | 63 | 98.9% |
| Baysor (M prior 1.0) | 31,764 | 74 | 98.9% |
| Baysor (10x prior 1.0) | 33,113 | 65 | 98.9% |

![Transcripts/cell distribution across expansion methods](results/figures/transcripts_per_cell.png)

Voronoi expansion captures 100% of transcripts by construction regardless of detector. Median transcripts per cell under Voronoi varies with detector quality: Voronoi (CP) leads at 149 tx/cell because CellPose detects fewer, larger nuclei, concentrating more transcripts per cell. Voronoi (10x) at 128 tx/cell is closest to 10x native (124), consistent with using the same nuclear seeds - the gap between them reflects 10x native's proprietary expansion assigning some transcripts differently than nearest-centroid.

Baysor without a prior captures 98.6% but detects fewer cells (18,321) - the density model merges adjacent cells freely. At PSC 0.2, the prior barely changes behavior (19,061 cells, 53 tx/cell). At PSC 0.8-1.0, the hard-locked nuclear seeds prevent merging and cell count jumps to ~30,000+. The four PSC=1.0 variants reveal how detector choice propagates through density-adaptive expansion: Baysor (SD prior) produces the most cells (34,230) because StarDist detects the most nuclei, while Baysor (M prior) has the highest median tx/cell (74) because Mesmer's larger nuclei anchor more cytoplasmic transcripts per cell. Baysor (10x prior) at 33,113 cells and 65 tx/cell falls between the two, and Baysor (CP prior) at 30,473 cells and 69 tx/cell reflects CellPose's more conservative nuclear detection. All four achieve ~99% transcript capture - the expansion strategy saturates regardless of which detector seeds it.

---

## Clustering comparison

Leiden clustering runs independently on each method's cells (normalize → PCA → neighbors → Leiden at resolution 1.0). Cluster labels are aligned across methods before computing confusion matrices and disagreement, using two algorithms: Hungarian (one-to-one) and argmax (many-to-one).

### Resolution stability

![ARI, disagreement, and Moran's I across Leiden resolutions - Hungarian alignment](results/figures/resolution_sensitivity_hungarian.png)

![ARI, disagreement, and Moran's I across Leiden resolutions - argmax alignment](results/figures/resolution_sensitivity_argmax.png)

The method ordering is stable across Leiden resolutions 0.3-2.0 under both alignment algorithms. Voronoi (Mesmer) leads at most resolutions (0.3, 0.6, 0.8-1.2); at resolutions 0.5 and 0.7, Voronoi (StarDist) briefly takes the lead, and at 1.5+ StarDist's higher cell count gives it a durable advantage as finer clustering demands more cells per cluster. Baysor without a prior is consistently lowest. The Hungarian alignment forces unmatched clusters into poor pairings when cluster counts differ, inflating disagreement for methods that produce more clusters. The argmax alignment lets multiple clusters map to the same reference cluster, reducing this artifact. The Moran's I panel confirms that the spatial-structure gap is resolution-invariant under both algorithms: Voronoi and Baysor prior methods maintain spatially structured disagreement while Baysor without a prior stays near zero regardless of cluster granularity.

![Leiden clustering comparison across methods](results/figures/cluster_comparison.png)

10x native and Voronoi methods converge on 14-16 clusters with median sizes above 1,000 cells. Baysor without a prior and at PSC 0.2 produce 21-24 smaller clusters, consistent with over-segmentation. At PSC 0.8-1.0, Baysor prior variants produce 20-23 clusters with higher cell counts (29,000-34,000) because the hard-locked nuclear seeds prevent merging; their median cluster sizes approach the Voronoi range.

UMAP embeddings colored by aligned cluster labels illustrate how the alignment algorithm reshapes cluster identity. Baysor without a prior shows the starkest contrast: Hungarian forces 6 of its 21 clusters into empty pairings, leaving large regions unmatched (gray), while argmax lets multiple Baysor clusters map to the same reference cluster, producing coherent coloring across the manifold.

| Baysor (Hungarian) | Baysor (argmax) |
|:---:|:---:|
| ![Baysor Hungarian](results/figures/umap/umap_baysor_hungarian.png) | ![Baysor argmax](results/figures/umap/umap_baysor_argmax.png) |

<details>
<summary>All individual UMAP plots</summary>

**Reference**

![10x native](results/figures/umap/umap_10x_native.png)

**Voronoi**

| Method | Hungarian | Argmax |
|---|---|---|
| Voronoi (CP) | ![](results/figures/umap/umap_voronoi_hungarian.png) | ![](results/figures/umap/umap_voronoi_argmax.png) |
| Voronoi (SD) | ![](results/figures/umap/umap_voronoi_stardist_hungarian.png) | ![](results/figures/umap/umap_voronoi_stardist_argmax.png) |
| Voronoi (M) | ![](results/figures/umap/umap_voronoi_mesmer_hungarian.png) | ![](results/figures/umap/umap_voronoi_mesmer_argmax.png) |
| Voronoi (10x) | ![](results/figures/umap/umap_voronoi_10x_ranger_hungarian.png) | ![](results/figures/umap/umap_voronoi_10x_ranger_argmax.png) |

**Baysor**

| Method | Hungarian | Argmax |
|---|---|---|
| Baysor | ![](results/figures/umap/umap_baysor_hungarian.png) | ![](results/figures/umap/umap_baysor_argmax.png) |
| Baysor (CP prior 0.2) | ![](results/figures/umap/umap_baysor_prior_hungarian.png) | ![](results/figures/umap/umap_baysor_prior_argmax.png) |
| Baysor (CP prior 0.8) | ![](results/figures/umap/umap_baysor_prior_c08_hungarian.png) | ![](results/figures/umap/umap_baysor_prior_c08_argmax.png) |
| Baysor (CP prior 1.0) | ![](results/figures/umap/umap_baysor_prior_c10_hungarian.png) | ![](results/figures/umap/umap_baysor_prior_c10_argmax.png) |
| Baysor (SD prior 1.0) | ![](results/figures/umap/umap_baysor_stardist_prior_c10_hungarian.png) | ![](results/figures/umap/umap_baysor_stardist_prior_c10_argmax.png) |
| Baysor (M prior 1.0) | ![](results/figures/umap/umap_baysor_mesmer_prior_c10_hungarian.png) | ![](results/figures/umap/umap_baysor_mesmer_prior_c10_argmax.png) |
| Baysor (10x prior 1.0) | ![](results/figures/umap/umap_baysor_10x_ranger_prior_c10_hungarian.png) | ![](results/figures/umap/umap_baysor_10x_ranger_prior_c10_argmax.png) |

</details>

### Cluster alignment

![Confusion matrices with Hungarian and argmax alignment](results/figures/confusion_clusters.png)

Each row is one 10x native cluster; columns are the comparison method's clusters. Red cells mark Hungarian (one-to-one) matched pairs, blue cells mark argmax (many-to-one) matches, and purple cells mark pairs selected by both algorithms. Voronoi methods produce clean matches under both algorithms. Baysor's 15×21 matrix shows the key difference: under Hungarian, 6 clusters are forced into empty pairings; under argmax, every column maps to the highest-overlap reference cluster with no wasted assignments.

#### Clustering agreement vs. 10x native

| Method | ARI | Hungarian Disagree | Hungarian Moran's I | Argmax Disagree | Argmax Moran's I |
| --- | ---: | ---: | ---: | ---: | ---: |
| **Nuclear-only** | | | | | |
| CellPose | 0.547 | 30.8% | 0.178 | 30.5% | 0.189 |
| StarDist | 0.545 | 33.5% | 0.215 | 33.4% | 0.221 |
| Mesmer | 0.557 | 27.9% | 0.090 | 27.4% | 0.098 |
| 10x Ranger | 0.504 | 35.0% | 0.191 | 34.5% | 0.203 |
| **Voronoi** | | | | | |
| Voronoi (CP) | 0.630 | 21.9% | 0.076 | 21.9% | 0.076 |
| Voronoi (SD) | 0.584 | 31.9% | 0.194 | 27.7% | 0.229 |
| Voronoi (M) | 0.686 | 18.8% | 0.161 | 18.8% | 0.161 |
| Voronoi (10x) | 0.592 | 28.3% | 0.172 | 25.8% | 0.168 |
| **Baysor** | | | | | |
| Baysor | 0.305 | 51.7% | 0.033 | 43.8% | 0.079 |
| Baysor (CP prior 0.2) | 0.318 | 51.9% | 0.036 | 39.2% | 0.086 |
| Baysor (CP prior 0.8) | 0.488 | 38.8% | 0.219 | 37.4% | 0.234 |
| Baysor (CP prior 1.0) | 0.501 | 33.8% | 0.112 | 32.1% | 0.122 |
| Baysor (SD prior 1.0) | 0.498 | 37.7% | 0.136 | 32.9% | 0.170 |
| Baysor (M prior 1.0) | 0.518 | 32.3% | 0.115 | 30.7% | 0.119 |
| Baysor (10x prior 1.0) | 0.530 | 34.7% | 0.208 | 33.1% | 0.204 |

Voronoi methods achieve the highest ARI (0.584-0.686), with Voronoi (M) leading. Nuclear-only methods cluster at ARI 0.504-0.557, and Baysor without a prior is lowest at 0.305. Argmax reduces Baysor's disagreement by ~8pp (51.7% to 43.8%) by eliminating forced mismatches from unmatched clusters. Voronoi methods with matched cluster counts are barely affected. The Moran's I increase for Baysor under argmax (0.033 to 0.079) shows that removing alignment noise reveals spatially structured disagreement that was previously masked.

![Per-cell-pair expression correlation](results/figures/expression_correlation.png)

Per-cell expression correlation is high for all methods (median 0.79-0.96), but cluster-label agreement tells a different story. Voronoi methods disagree with 10x native on 19-32% of matched cells. Baysor without a prior disagrees on 52% (Hungarian) / 44% (argmax); at PSC 0.2, disagreement is unchanged (52% / 39%), but PSC 1.0 variants drop to 32-38% (Hungarian) / 31-33% (argmax), approaching the Voronoi range.

### Per-cluster pseudobulk

![Per-cluster pseudobulk correlation vs. 10x native](results/figures/pseudobulk_by_cluster.png)

To test whether cluster-level expression profiles agree, matched cells are grouped by 10x native's 15 Leiden clusters and pseudobulked per method. Nuclear methods drop to r = 0.86-0.87 on luminal epithelial clusters (0, 1, 3, 8) - the same populations driving single-cell disagreement - while Voronoi variants stay above 0.99 across all clusters. Baysor shows a comparable luminal dip plus reduced correlation on macrophage clusters (2, 7), consistent with transcript-density boundaries partitioning those populations differently.

---

## Spatial structure of disagreement

### Hungarian alignment

![Disagreement mapped spatially - Hungarian](results/figures/disagreement_spatial_map.png)

![LISA hotspot/coldspot maps - Hungarian](results/figures/local_morans_map.png)

| Comparison | Global Moran's I | HH hotspots | LL coldspots |
| --- | --- | --- | --- |
| 10x native vs. CellPose | 0.178 | 21.7% | 30.3% |
| 10x native vs. StarDist | 0.215 | 18.6% | 15.0% |
| 10x native vs. Mesmer | 0.090 | 17.1% | 32.5% |
| 10x native vs. Voronoi (CP) | 0.076 | 11.1% | 27.2% |
| 10x native vs. Voronoi (SD) | 0.194 | 23.2% | 30.8% |
| 10x native vs. Voronoi (M) | 0.161 | 9.5% | 20.4% |
| 10x native vs. Baysor | 0.033 | 21.4% | 17.5% |

### Argmax alignment

![Disagreement mapped spatially - argmax](results/figures/disagreement_spatial_map_argmax.png)

![LISA hotspot/coldspot maps - argmax](results/figures/local_morans_map_argmax.png)

| Comparison | Global Moran's I | HH hotspots | LL coldspots |
| --- | --- | --- | --- |
| 10x native vs. CellPose | 0.189 | 21.5% | 29.3% |
| 10x native vs. StarDist | 0.221 | 18.6% | 14.7% |
| 10x native vs. Mesmer | 0.098 | 16.6% | 31.4% |
| 10x native vs. Voronoi (CP) | 0.076 | 11.1% | 27.2% |
| 10x native vs. Voronoi (SD) | 0.229 | 18.7% | 25.6% |
| 10x native vs. Voronoi (M) | 0.161 | 9.5% | 20.4% |
| 10x native vs. Baysor | 0.079 | 26.1% | 26.2% |

Nuclear and Voronoi disagreements are spatially structured (Moran's I 0.076-0.215 under Hungarian), concentrated in luminal epithelial territory. Mesmer has the most agreement coldspots (32.5% LL); Voronoi (Mesmer) has the fewest disagreement hotspots (9.5% HH), consistent with residual errors being diffuse boundary noise. Under Hungarian alignment Baysor's near-zero Moran's I (0.033) reflects noise from forced cluster mismatches; under argmax alignment Moran's I increases to 0.079, revealing that Baysor's genuine disagreements are spatially structured - just less so than morphological methods.

*LISA labels*: HH = disagreement hotspot; LL = agreement coldspot. Global Moran's I summarises spatial autocorrelation of the disagree flag (0 = random, 1 = fully clustered).

---

## Cell-type sensitivity

![Cell type vs. agreement - Hungarian](results/figures/agreement_explainer.png)

![Cell type vs. agreement - argmax](results/figures/agreement_explainer_argmax.png)

Adipocytes and myoepithelial cells have the highest per-cell disagreement (~50-68% and ~40-47%) but are rare. Luminal epithelial cells dominate by volume: ~35% disagreement across ~8,500 cells drives the majority of total disagreement events. These clusters likely encompass malignant and normal epithelial cells; both share canonical markers (GATA3, PGR, ESR1, MUC1) and are inseparable by nuclear morphology alone. T cells and B cells are robustly identified regardless of method or alignment algorithm.

---

## Disagreement drivers: cell state vs. geometry

![Phenotypic density vs. disagreement - Hungarian](results/figures/density_vs_disagreement.png)

![Phenotypic density vs. disagreement - argmax](results/figures/density_vs_disagreement_argmax.png)

![DE: disagree vs. agree cells - Hungarian](results/figures/de_volcano.png)

![DE: disagree vs. agree cells - argmax](results/figures/de_volcano_argmax.png)

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

---

## Phenotypic landscape distortion

![All methods in shared PCA/UMAP space](results/figures/manifold_shared_umap.png)

![Phenotypic landscape distortion vs. 10x native](results/figures/manifold_distortion.png)

All methods are projected into a shared PCA space fit on 10x native (30 PCs, 55% variance explained) and embedded in a joint UMAP. Density ratio maps (log₂ method/10x) show which phenotypic regions each method enriches or depletes. Nuclear methods show depleted regions in high-density luminal epithelial areas, consistent with missed cytoplasmic transcripts pulling cells toward lower-expression PCA states. Voronoi methods track 10x native closely. Baysor shows enrichment in a distinct region corresponding to its finer resolution of macrophage and stromal subtypes.

---

## Pairwise method agreement

![Pairwise ARI between all segmentation methods](results/figures/pairwise_consensus.png)

| | 10x native | CellPose | StarDist | Mesmer | Voronoi (CP) | Voronoi (SD) | Voronoi (M) | Baysor |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **10x native** | 1.0 | 0.547 | 0.545 | 0.557 | 0.630 | 0.584 | 0.686 | 0.305 |
| **CellPose** | | 1.0 | 0.764 | 0.606 | 0.533 | 0.557 | 0.556 | 0.415 |
| **StarDist** | | | 1.0 | 0.639 | 0.525 | 0.568 | 0.540 | 0.411 |
| **Mesmer** | | | | 1.0 | 0.480 | 0.535 | 0.528 | 0.455 |
| **Voronoi (CP)** | | | | | 1.0 | 0.581 | 0.661 | 0.336 |
| **Voronoi (SD)** | | | | | | 1.0 | 0.633 | 0.361 |
| **Voronoi (M)** | | | | | | | 1.0 | 0.377 |

No, at least not within the Voronoi family. CellPose and StarDist agree with each other at ARI 0.764 (higher than the Voronoi pair at 0.661) because both are nuclear-morphology methods on the same DAPI image. Switching to Voronoi assignment lowers within-paradigm agreement because the two Voronoi variants use different centroids, shifting boundaries even where centroids are close. What Voronoi does raise is agreement with the 10x-native whole-cell reference (0.63-0.69): compatibility with the platform's own segmentation, not cross-method reproducibility. Baysor remains isolated from all morphological methods (ARI 0.30-0.46 regardless of partner).

---

## Cell size and disagreement

![Cell size vs. disagreement - Hungarian](results/figures/cell_size_disagreement.png)

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

Smaller 10x-native cells are significantly more likely to disagree with every morphological method (p ≪ 0.001). The direction is counter-intuitive: larger cells have more cytoplasm, yet it is smaller cells that disagree more. The pattern holds for Voronoi methods too, ruling out transcript capture as the cause. Smaller cells likely correspond to densely packed regions where any method's cluster assignment is noisier. Baysor shows no size dependence (p = 0.28); its boundaries are insensitive to morphologically defined cell area.

---

## Marker gene recovery

![Marker gene recovery relative to 10x native](results/figures/marker_recovery.png)

Using 10x-native cell-type annotations as ground truth, nuclear methods recover 75-92% of cytoplasmic marker expression relative to 10x native, with the largest deficits for extranuclear markers like MUC1, SERPINA3, and LYZ. Voronoi methods recover near-100% across all cell types. Baysor recovers macrophage markers (LYZ, CD14) at or above 10x-native levels while showing slightly reduced T cell marker (CD3E) recovery.

---

## Population-level convergence

![Pseudobulk correlation vs. 10x native](results/figures/pseudobulk_correlation.png)

| Method | Per-cell-type pseudobulk r (range) | Aggregate r | Single-cell ARI |
| --- | --- | --- | --- |
| CellPose | 0.87-0.98 | 0.970 | 0.547 |
| StarDist | 0.88-0.99 | 0.975 | 0.545 |
| Mesmer | 0.92-0.99 | 0.983 | 0.557 |
| Voronoi (CP) | 0.98-1.00 | 0.9999 | 0.630 |
| Voronoi (SD) | 0.98-1.00 | 0.9999 | 0.584 |
| Voronoi (M) | 0.98-1.00 | 0.9999 | 0.686 |
| Baysor | 0.94-1.00 | 0.999 | 0.305 |

Pseudobulk is computed within each of 10 annotated cell types (not as a whole-ROI sum), so the correlation tests whether each method's cell-type compartments recover the same expression programs as 10x native. Baysor's per-cell-type correlations range from 0.94 (plasma cells) to 0.997 (CAFs), degrading predictably on rare populations with fewer cells. Despite its low single-cell ARI of 0.305, Baysor is competitive with nuclear methods at the cell-type level - its aggregate r of 0.999 sits above CellPose (0.970) and StarDist (0.975). Nuclear methods show reduced pseudobulk r (0.97-0.98) because missing cytoplasmic transcripts suppress marker signal systematically across all cells of a type. Voronoi methods achieve both high single-cell ARI and near-perfect pseudobulk agreement.

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
