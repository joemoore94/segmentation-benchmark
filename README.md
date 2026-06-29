# Segmentation Benchmarking on Xenium Spatial Transcriptomics Data

**Question:** Do nuclear-mask (CellPose, StarDist, Mesmer), Voronoi, and transcript-density (Baysor) segmentation methods transfer well to Xenium spatial transcriptomics, and does method choice meaningfully change downstream cell-type calls?

## Summary

Cluster labels are aligned via two algorithms before computing disagreement and Moran's I. Hungarian finds the optimal one-to-one assignment; when cluster counts differ, unmatched clusters are forced into poor pairings. Argmax maps each method's clusters to the 10x native cluster with plurality overlap, allowing many-to-one mapping. Matched pairs, median correlation, and ARI do not depend on cluster alignment.

<details>
<summary><b>Full results table</b> — click to expand</summary>

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

</details>

The results cleanly separate the contributions of nuclear detection quality and expansion strategy. Among Voronoi methods, Mesmer centroids produce the highest ARI (0.686) and lowest disagreement (18.8%), while 10x Ranger centroids - despite being purpose-built for Xenium - score lower (ARI 0.592). Among Baysor PSC=1.0 variants, the same detector ordering holds: Mesmer prior leads at ARI 0.518, followed by 10x Ranger (0.530), CellPose (0.501), and StarDist (0.498). Baysor (10x prior 1.0) achieves the highest ARI of any Baysor variant (0.530) and the most matched pairs (22,910), benefiting from 10x Ranger detecting almost as many nuclei as the 10x native reference.

Across both expansion strategies, Voronoi consistently outperforms Baysor PSC=1.0 on ARI (0.58-0.69 vs 0.50-0.53) - geometric nearest-centroid assignment agrees more with 10x native than density-adaptive expansion. However, Baysor prior variants have lower negative marker violation rates (0.31 per 1000 tx vs 0.37-0.43 for Voronoi), indicating that density-adaptive boundaries produce fewer cross-lineage contamination artifacts even when they disagree with the 10x reference.

Nuclear-only methods capture too few transcripts for meaningful downstream comparison and are excluded from figures past the recovery section.

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

<p align="center"><img src="results/figures/nuclear_mask_sizes.png" alt="Nuclear mask size comparison across all four detectors" width="600"></p>

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

### Per-cell expression correlation

<p align="center"><img src="results/figures/expression_correlation.png" alt="Per-cell expression correlation vs. 10x native" width="600"></p>

Per-cell expression correlation is high for all methods (median 0.79-0.96). Voronoi methods lead (median r 0.96) because nearest-centroid expansion assigns largely the same transcripts as 10x native. Baysor PSC=1.0 variants cluster at 0.90-0.92, and nuclear-only methods at 0.82-0.88, reflecting missing cytoplasmic transcripts rather than misassignment.

---

## Clustering comparison

Leiden clustering runs independently on each method's cells (normalize → PCA → neighbors → Leiden at resolution 1.0). Cluster labels are aligned across methods before computing confusion matrices and disagreement, using two algorithms: Hungarian (one-to-one) and argmax (many-to-one).

### Resolution stability

<p align="center"><img src="results/figures/resolution_ari.png" alt="ARI across Leiden resolutions (alignment-agnostic)" width="600"></p>

ARI is partition-based and does not depend on cluster alignment, so it is the same under Hungarian and argmax. The method ordering is stable across Leiden resolutions 0.3-2.0. Voronoi (Mesmer) leads at most resolutions (0.3, 0.6, 0.8-1.2); at resolutions 0.5 and 0.7, Voronoi (StarDist) briefly takes the lead, and at 1.5+ StarDist's higher cell count gives it a durable advantage as finer clustering demands more cells per cluster. Baysor without a prior is consistently lowest.

![Disagreement and Moran's I across Leiden resolutions - Hungarian alignment](results/figures/resolution_disagree_morans_hungarian.png)

![Disagreement and Moran's I across Leiden resolutions - argmax alignment](results/figures/resolution_disagree_morans_argmax.png)

Disagreement and Moran's I do depend on alignment. The Hungarian alignment forces unmatched clusters into poor pairings when cluster counts differ, inflating disagreement for methods that produce more clusters. The argmax alignment lets multiple clusters map to the same reference cluster, reducing this artifact. The Moran's I panel confirms that the spatial-structure gap is resolution-invariant under both algorithms: Voronoi and Baysor prior methods maintain spatially structured disagreement while Baysor without a prior stays near zero regardless of cluster granularity.

UMAP embeddings colored by aligned cluster labels illustrate how the alignment algorithm reshapes cluster identity. Baysor without a prior shows the starkest contrast: Hungarian forces 6 of its 21 clusters into empty pairings, leaving large regions unmatched (gray), while argmax lets multiple Baysor clusters map to the same reference cluster, producing coherent coloring across the manifold.

![Baysor UMAP: Hungarian vs argmax alignment](results/figures/umap_baysor_alignment_comparison.png)

<details>
<summary><b>All method UMAPs — Hungarian alignment</b> — click to expand</summary>

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

</details>

<details>
<summary><b>All method UMAPs — argmax alignment</b> — click to expand</summary>

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

</details>

<details>
<summary><b>Per-cluster pseudobulk</b> — click to expand</summary>

<p align="center"><img src="results/figures/pseudobulk_by_cluster.png" alt="Per-cluster pseudobulk correlation vs. 10x native" width="600"></p>

To test whether cluster-level expression profiles agree, matched cells are grouped by 10x native's 15 Leiden clusters and pseudobulked per method. Nuclear methods drop to r = 0.86-0.87 on luminal epithelial clusters (0, 1, 3, 8) - the same populations driving single-cell disagreement - while Voronoi variants stay above 0.99 across all clusters. Baysor shows a comparable luminal dip plus reduced correlation on macrophage clusters (2, 7), consistent with transcript-density boundaries partitioning those populations differently.

</details>

---

## Cluster alignment

<details>
<summary><b>Confusion matrices (Hungarian and argmax)</b> — click to expand</summary>

![Confusion matrices with Hungarian and argmax alignment](results/figures/confusion_clusters.png)

Each row is one 10x native cluster; columns are the comparison method's clusters. Red cells mark Hungarian (one-to-one) matched pairs, blue cells mark argmax (many-to-one) matches, and purple cells mark pairs selected by both algorithms. Voronoi methods produce clean matches under both algorithms. Baysor's 15x21 matrix shows the key difference: under Hungarian, 6 clusters are forced into empty pairings; under argmax, every column maps to the highest-overlap reference cluster with no wasted assignments.

</details>

Hungarian produces min(n_10x, n_method) matched pairs; when cluster counts differ, the excess clusters are forced into poor pairings. Argmax maps every method cluster to its best-overlap 10x cluster (many-to-one), so all method clusters are assigned but some 10x clusters may receive multiple mappings while others receive none. At resolution 1.0, 10x native produces 15 clusters; Voronoi methods produce 14-16 (near-perfect matching under both algorithms) while Baysor without a prior produces 21 (6 unmatched under Hungarian, all 21 mapped under argmax).

### Clustering agreement vs. 10x native

![Clustering agreement vs. 10x native](results/figures/clustering_agreement.png)

Voronoi methods achieve the highest ARI (0.584-0.686), with Voronoi (M) leading. Nuclear-only methods cluster at ARI 0.504-0.557, and Baysor without a prior is lowest at 0.305. Argmax reduces Baysor's disagreement by ~8pp (51.7% to 43.8%) by eliminating forced mismatches from unmatched clusters. Voronoi methods with matched cluster counts are barely affected. The Moran's I increase for Baysor under argmax (0.033 to 0.079) shows that removing alignment noise reveals spatially structured disagreement that was previously masked.

---

## Cell type annotation

Cell types are annotated on the 10x native segmentation only. Leiden clustering (resolution 1.0) partitions the 10x native cells into 15 clusters, then differential expression (DE) identifies each cluster's distinguishing genes. DE compares the expression of every gene in the cells of one cluster against all other cells using a Wilcoxon rank-sum test, ranking genes by how strongly and specifically they are upregulated in that cluster. The top DE genes are matched to canonical breast tissue markers to assign a cell type label. These 10x native annotations serve as the reference for all downstream cross-method comparisons. The raw Xenium output carries no cell type labels - only coordinates and transcript counts.

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

The dotplot below shows which canonical markers are expressed in each cluster, confirming the assignments. Dot size encodes the percentage of cells in the cluster expressing the marker; color encodes mean expression level. Clusters 0, 1, 3, and 8 all annotate as luminal epithelial but are distinguished by different marker profiles: cluster 0 is ESR1/FOXA1-dominant (ER+ hormone-responsive), cluster 1 is PGR/MUC1-dominant, cluster 3 expresses stromal-adjacent markers (NNMT, LUM), and cluster 8 is TACSTD2/KRT7/STC2-dominant (proliferative/stress-response). Clusters 2 and 7 are both macrophage populations: cluster 2 (335 cells) expresses FCGR3A and HAVCR2 (non-classical/M2-like), while cluster 7 (2,612 cells) expresses CD14 and AIF1 (classical monocyte-derived). Clusters 9 and 13 are both CAFs: cluster 9 expresses SFRP4 and FBLN1 (matrix-producing), while cluster 13 expresses POSTN and CTHRC1 (myofibroblastic).

<details>
<summary><b>Canonical marker expression by Leiden cluster</b> — click to expand</summary>

![Canonical marker expression by Leiden cluster](results/figures/annotation_dotplot_flipped.png)

</details>

## Luminal epithelial subtypes

The tissue sample is an invasive ductal carcinoma (IDC) of the breast, with the ROI spanning several tumor nests (DCIS and invasive) surrounded by stroma and immune infiltrate. At Leiden resolution 1.0, four of the 15 clusters (0, 1, 3, 8) annotate as luminal epithelial, together accounting for 8,548 cells — 36% of the ROI. Wilcoxon DE within this subpopulation reveals that the four clusters are not redundant: each has a distinct marker profile corresponding to known luminal biology, consistent with the heterogeneity expected in IDC where hormone-receptor-positive and proliferative tumor populations coexist within the same lesion.

| Subcluster | Cells | Top markers | Interpretation |
| ---: | ---: | --- | --- |
| 0 | 2,799 | FLNB, DHRS2, MLPH, ESR1, PDZK1, AGR3 | ER+ hormone-responsive (luminal A-like) |
| 1 | 2,926 | MYBPC1, PGR, CLIC6, FASN, MUC1, ELOVL5, SCD | Secretory luminal with lipid metabolism |
| 3 | 902 | NNMT, ENAH, LUM, POSTN, AQP3 | Stromal-adjacent; EMT-like or boundary contamination |
| 8 | 1,921 | EGR1, STC2, KRT7, CCND1, TACSTD2, TCIM | Stress-response / proliferative |

![Top markers per luminal epithelial subcluster](results/figures/lum_ep_subtypes/lum_ep_subcluster_dotplot.png)

![Luminal epithelial subcluster heatmap](results/figures/lum_ep_subtypes/lum_ep_subcluster_heatmap.png)

Clusters 0 and 1 map to recognized luminal subtypes in breast cancer: cluster 0 expresses the canonical ER+ program (ESR1, MLPH, PDZK1, AGR3) while cluster 1 emphasizes progesterone signaling and lipid biosynthesis (PGR, FASN, ELOVL5, SCD). Cluster 8 is dominated by immediate-early genes (EGR1 logFC +3.4, STC2 logFC +4.2) and CCND1, consistent with an actively cycling or stress-responsive state. Cluster 3 is the most ambiguous — its top markers include NNMT (metabolic stress), plus LUM and POSTN, which are canonical CAF markers and in fact define clusters 9 and 13 in the full annotation. This could reflect genuine epithelial-mesenchymal transition in a subset of tumor cells, or it could indicate that 10x native's expansion algorithm is pulling stromal transcripts into epithelial cell boundaries in densely packed regions. Either way, this subcluster sits at the lineage boundary where segmentation errors are most consequential.

### Resolution sweep

![Leiden resolution sweep — luminal epithelial splitting](results/figures/lum_ep_subtypes/lum_ep_resolution_sweep.png)

A Leiden resolution sweep from 0.1 to 2.0 reveals that the four subclusters do not emerge simultaneously. Cluster 8 (EGR1/STC2/KRT7) splits off first at resolution 0.3, consistent with its strong immediate-early gene signature making it the most transcriptionally distinct subpopulation. Clusters 0 and 1 separate from each other at resolution 0.6 as the ESR1-high and PGR/FASN-high programs become resolvable. Cluster 3 (NNMT/LUM/POSTN) is the last to emerge, not appearing as a distinct cluster until resolution 0.9. The late splitting and stromal marker expression together make cluster 3 the weakest of the four subtypes — the one most likely to represent a segmentation boundary effect rather than a discrete biological state.

![Luminal epithelial subclusters in UMAP context](results/figures/lum_ep_subtypes/lum_ep_umap_context.png)

The UMAP confirms this hierarchy: cluster 8 occupies a clearly separated island within the luminal epithelial territory, while clusters 0 and 1 are adjacent but distinguishable, and cluster 3 cells are scattered along the boundary between luminal epithelial and CAF/stromal regions. The four-subcluster structure within luminal epithelial — and particularly the graded transcriptional boundaries between them — helps explain why this cell type drives the plurality of cross-method disagreement: a handful of misassigned boundary transcripts can shift a cell between subtypes without crossing a sharp lineage boundary.

---

## Spatial structure of disagreement

![Spatial autocorrelation of disagreement - all methods](results/figures/spatial_morans_dotplot.png)

<details>
<summary><b>Spatial disagreement maps (Hungarian alignment)</b> - click to expand</summary>

![Disagreement mapped spatially - Hungarian](results/figures/disagreement_spatial_map.png)

</details>

<details>
<summary><b>Spatial disagreement maps (Argmax alignment)</b> - click to expand</summary>

![Disagreement mapped spatially - argmax](results/figures/disagreement_spatial_map_argmax.png)

</details>

<details>
<summary><b>LISA hotspot/coldspot maps (Hungarian)</b> — click to expand</summary>

![LISA hotspot/coldspot maps - Hungarian](results/figures/local_morans_map.png)

</details>

<details>
<summary><b>LISA hotspot/coldspot maps (argmax)</b> — click to expand</summary>

![LISA hotspot/coldspot maps - argmax](results/figures/local_morans_map_argmax.png)

</details>

| Comparison | Global Moran's I | HH hotspots | LL coldspots |
| --- | --- | --- | --- |
| 10x native vs. CellPose | 0.178 | 21.7% | 30.3% |
| 10x native vs. StarDist | 0.215 | 18.6% | 15.0% |
| 10x native vs. Mesmer | 0.090 | 17.1% | 32.5% |
| 10x native vs. Voronoi (CP) | 0.076 | 11.1% | 27.2% |
| 10x native vs. Voronoi (SD) | 0.194 | 23.2% | 30.8% |
| 10x native vs. Voronoi (M) | 0.161 | 9.5% | 20.4% |
| 10x native vs. Baysor | 0.033 | 21.4% | 17.5% |

Nuclear and Voronoi disagreements are spatially structured (Moran's I 0.076-0.215 under Hungarian), concentrated in luminal epithelial territory. Mesmer has the most agreement coldspots (32.5% LL); Voronoi (Mesmer) has the fewest disagreement hotspots (9.5% HH), consistent with residual errors being diffuse boundary noise. Under Hungarian alignment Baysor's near-zero Moran's I (0.033) reflects noise from forced cluster mismatches; under argmax alignment Moran's I increases to 0.079, revealing that Baysor's genuine disagreements are spatially structured - just less so than morphological methods.

*LISA labels*: Local Moran's I (Anselin 1995) decomposes the global statistic to a per-cell level. Each cell's local I measures whether it and its k-nearest spatial neighbors share similar disagree/agree status. Cells are classified by the sign of their local I and their own value: HH = disagreeing cell surrounded by disagreeing neighbors (hotspot), LL = agreeing cell surrounded by agreeing neighbors (coldspot), HL/LH = spatial outliers. Global Moran's I summarises spatial autocorrelation of the disagree flag across the whole ROI (0 = random, 1 = fully clustered).

---

## Cell-type sensitivity

<details>
<summary><b>Cell type vs. agreement (Hungarian)</b> — click to expand</summary>

![Cell type vs. agreement - Hungarian](results/figures/agreement_explainer.png)

</details>

<details>
<summary><b>Cell type vs. agreement (Argmax)</b> — click to expand</summary>

![Cell type vs. agreement - argmax](results/figures/agreement_explainer_argmax.png)

</details>

Adipocytes and myoepithelial cells have the highest per-cell disagreement (~50-68% and ~40-47%) but are rare. Luminal epithelial cells dominate by volume: ~35% disagreement across ~8,500 cells drives the majority of total disagreement events. These clusters likely encompass malignant and normal epithelial cells; both share canonical markers (GATA3, PGR, ESR1, MUC1) and are inseparable by nuclear morphology alone. T cells and B cells are robustly identified regardless of method or alignment algorithm.

---

## Disagreement drivers: cell state vs. geometry

<details>
<summary><b>Phenotypic density vs. disagreement (Hungarian)</b> — click to expand</summary>

![Phenotypic density vs. disagreement - Hungarian](results/figures/density_vs_disagreement.png)

</details>

<details>
<summary><b>Phenotypic density vs. disagreement (Argmax)</b> — click to expand</summary>

![Phenotypic density vs. disagreement - argmax](results/figures/density_vs_disagreement_argmax.png)

</details>

<details>
<summary><b>DE: disagree vs. agree cells (Hungarian)</b> — click to expand</summary>

![DE: disagree vs. agree cells - Hungarian](results/figures/de_volcano.png)

</details>

<details>
<summary><b>DE: disagree vs. agree cells (Argmax)</b> — click to expand</summary>

![DE: disagree vs. agree cells - argmax](results/figures/de_volcano_argmax.png)

</details>

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

<details>
<summary><b>All methods in shared PCA/UMAP space</b> — click to expand</summary>

![All methods in shared PCA/UMAP space](results/figures/manifold_shared_umap.png)

</details>

<details>
<summary><b>Density distortion vs. 10x native</b> — click to expand</summary>

![Phenotypic landscape distortion vs. 10x native](results/figures/manifold_distortion.png)

</details>

All methods are projected into a shared PCA space fit on 10x native (30 PCs, 55% variance explained) and embedded in a joint UMAP. Density ratio maps (log₂ method/10x) show which phenotypic regions each method enriches or depletes. Nuclear methods show depleted regions in high-density luminal epithelial areas, consistent with missed cytoplasmic transcripts pulling cells toward lower-expression PCA states. Voronoi methods track 10x native closely. Baysor shows enrichment in a distinct region corresponding to its finer resolution of macrophage and stromal subtypes.

---

## Comparison to scRNA-seq reference

The phenotypic landscape distortion analysis above uses a PCA space fit on 10x native, so the coordinate system is itself a segmentation output. To test whether method choice shifts cells in a segmentation-independent reference, PCA is fit on companion scRNA-seq from the same tissue blocks (GSE243275, 3' chemistry, 7,329 cells after QC) subsetted to the 374 Xenium panel genes present in both datasets. Each segmentation method's cells are then projected into this reference PCA space (30 PCs, 62% variance explained), clustered via Leiden in the shared reference embedding, and compared.

<details>
<summary><b>Xenium cells in scRNA-seq reference UMAP</b> — click to expand</summary>

![Xenium cells projected into scRNA-seq reference UMAP](results/figures/ref_projection_umap.png)

</details>

<details>
<summary><b>Density distortion vs. scRNA-seq reference</b> — click to expand</summary>

![Density distortion vs. scRNA-seq reference](results/figures/ref_projection_density.png)

</details>

All methods project into the same regions of the scRNA-seq reference landscape, but Baysor without a prior collapses into a narrow band — consistent with its lower transcript-per-cell counts compressing the phenotypic range. Voronoi methods and 10x native show nearly identical density profiles, with enrichment in the luminal epithelial region and depletion in populations that are better represented in the dissociated scRNA-seq (e.g. immune subtypes lost during tissue sectioning). Baysor PSC=1.0 variants are intermediate.

### Clustering in reference space

<details>
<summary><b>Reference-space Leiden clustering UMAPs</b> — click to expand</summary>

![Reference-space Leiden clustering](results/figures/ref_projection_ref_clustering.png)

</details>

<details>
<summary><b>Reference-space clusters mapped to tissue</b> — click to expand</summary>

![Reference-space clusters mapped to tissue](results/figures/ref_projection_spatial.png)

</details>

Leiden clustering in the shared reference PCA space (resolution 1.0) produces 12–14 clusters for Voronoi methods and 10x native, but 22–25 for Baysor variants. The spatial maps confirm that reference-space clusters are biologically coherent: tissue structures (ducts, stroma, immune infiltrate) are visible across methods despite the clustering being performed in an independent coordinate system.

<p align="center"><img src="results/figures/ref_projection_ari.png" alt="Pairwise ARI in reference PCA space" width="600"></p>

To test whether the reference projection genuinely improves cell typing or merely compresses method differences ("graying"), ARI is computed against a fixed anchor: 10x native's own-PCA clustering. The "before" column clusters each method in its own PCA space; the "after" column clusters in the reference PCA space. Both are compared to the same 10x native own-space labels.

| Method | Ref clusters | ARI (before) | ARI (after) | ΔARI |
| --- | ---: | ---: | ---: | ---: |
| 10x native | 14 | 1.000 | 0.522 | −0.478 |
| Voronoi (CP) | 12 | 0.568 | 0.454 | −0.114 |
| Voronoi (SD) | 13 | 0.555 | 0.462 | −0.093 |
| Voronoi (M) | 13 | 0.614 | 0.480 | −0.135 |
| Voronoi (10x) | 13 | 0.624 | 0.500 | −0.124 |
| Baysor | 23 | 0.172 | 0.181 | +0.009 |
| Baysor (CP prior 1.0) | 22 | 0.327 | 0.266 | −0.061 |
| Baysor (SD prior 1.0) | 25 | 0.306 | 0.290 | −0.017 |
| Baysor (M prior 1.0) | 25 | 0.324 | 0.288 | −0.036 |
| Baysor (10x prior 1.0) | 24 | 0.343 | 0.291 | −0.052 |

Even 10x native's own reference-space clustering only achieves ARI 0.522 against its internal clustering, and every method except prior-free Baysor moves away from the original 10x native cell types after projection. The pairwise convergence observed between methods in the reference space (mean off-diagonal ARI 0.570 → 0.631) therefore reflects a shared loss of resolution rather than convergence toward ground truth. The 374-gene panel captures 62% of variance in the reference but emphasizes different axes of variation than each method's own PCA, compressing cell-type boundaries in the process. The density and displacement analyses above remain valid since they measure geometry rather than clustering, but the reference PCA is not a better basis for cell typing than each method's own space.

### Cell-level displacement

<details>
<summary><b>Matched-cell displacement by method</b> — click to expand</summary>

![Matched-cell displacement by method](results/figures/ref_projection_displacement.png)

</details>

<details>
<summary><b>Matched-cell displacement by cell type</b> — click to expand</summary>

![Matched-cell displacement by cell type](results/figures/ref_projection_displacement_by_ct.png)

</details>

| Method | Matched cells | Median displacement |
| --- | ---: | ---: |
| Voronoi (CP) | 23,267 | 1.41 |
| Voronoi (SD) | 23,515 | 1.08 |
| Voronoi (M) | 23,580 | 1.18 |
| Voronoi (10x) | 23,626 | 0.92 |
| Baysor | 22,344 | 2.97 |
| Baysor (CP prior 1.0) | 23,626 | 2.27 |
| Baysor (SD prior 1.0) | 23,628 | 2.41 |
| Baysor (M prior 1.0) | 23,628 | 2.01 |
| Baysor (10x prior 1.0) | 23,629 | 2.28 |

For each matched cell pair (nearest centroid, <15 µm), displacement measures how far the same cell moves in the 30-PC reference space depending on which segmentation method assigned its transcripts. Voronoi methods displace cells 0.92–1.41 units from their 10x native position; Baysor displaces them 2.01–2.97 — roughly 2× further. Voronoi (10x) has the smallest displacement (0.92) because it shares nuclear seeds with 10x native, differing only in the expansion algorithm. The cell-type breakdown shows that Baysor's displacement is elevated across all populations but peaks on plasma cells (3.8), T cells (3.4), and B cells (3.1) — rare or small-cell populations where the density model has few transcripts per cell and boundary assignments are most uncertain. Luminal epithelial displacement is moderate for Baysor (2.2–2.4) but drives the largest absolute number of displaced cells given its 36% share of the ROI. Baysor (M prior 1.0) has the lowest displacement among Baysor variants (2.01), consistent with Mesmer's larger nuclear masks anchoring more transcripts and reducing cytoplasmic boundary noise.

### Centroid distance to scRNA-seq reference

To measure how closely each segmentation method's clusters match scRNA-seq cell states, all methods are projected into the shared reference PCA space, clustered at each Leiden resolution, and the Euclidean distance from each method cluster centroid to the nearest scRNA-seq cluster centroid is computed. Both the reference and method resolution are swept independently (0.3–2.0), giving a 10×10 grid of resolution combinations per method.

The two key slices through this grid reveal qualitatively different behavior between method families. Fixing the reference at resolution 1.0 (17 clusters) and sweeping the method resolution asks: does splitting method clusters finer bring their centroids closer to the reference?

| Method | 0.3 | 0.5 | 0.7 | 1.0 | 1.5 | 2.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 10x native | 4.09 | 4.11 | 4.09 | 4.06 | 4.00 | 4.01 |
| Voronoi (CP) | 4.05 | 4.05 | 4.00 | 3.91 | 3.81 | 3.78 |
| Voronoi (SD) | 4.24 | 4.16 | 4.14 | 4.10 | 4.09 | 4.01 |
| Voronoi (M) | 4.08 | 4.04 | 3.98 | 3.88 | 3.80 | 3.81 |
| Voronoi (10x) | 4.13 | 4.08 | 4.07 | 4.06 | 4.00 | 3.98 |
| Baysor | 4.57 | 4.58 | 4.58 | 4.57 | 4.57 | 4.58 |
| Baysor (CP prior 1.0) | 4.25 | 4.22 | 4.22 | 4.23 | 4.23 | 4.23 |
| Baysor (SD prior 1.0) | 4.32 | 4.31 | 4.30 | 4.29 | 4.29 | 4.29 |
| Baysor (M prior 1.0) | 4.24 | 4.22 | 4.22 | 4.23 | 4.25 | 4.24 |
| Baysor (10x prior 1.0) | 4.28 | 4.28 | 4.27 | 4.27 | 4.26 | 4.27 |

Voronoi methods improve modestly with finer resolution (−0.27 from 0.3→2.0), while Baysor methods are completely flat (Δ ≈ 0.00). Splitting Baysor from 11 to 37 clusters does not move centroids closer to the reference — they are fixed in the wrong locations in PCA space.

Fixing the method at resolution 1.0 and sweeping the reference resolution asks: does a denser reference grid help methods find nearby centroids?

| Method | 0.3 | 0.5 | 0.7 | 1.0 | 1.5 | 2.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 10x native | 3.96 | 3.90 | 4.13 | 4.06 | 4.05 | 3.98 |
| Voronoi (CP) | 3.91 | 3.90 | 3.94 | 3.91 | 3.90 | 3.82 |
| Voronoi (SD) | 3.96 | 3.91 | 4.12 | 4.10 | 4.10 | 3.98 |
| Voronoi (M) | 3.90 | 3.87 | 3.94 | 3.88 | 3.87 | 3.81 |
| Voronoi (10x) | 3.98 | 3.92 | 4.12 | 4.06 | 4.05 | 3.98 |
| Baysor | 4.29 | 4.43 | 4.38 | 4.57 | 4.58 | 3.51 |
| Baysor (CP prior 1.0) | 3.94 | 4.03 | 4.13 | 4.23 | 4.23 | 3.58 |
| Baysor (SD prior 1.0) | 3.97 | 4.08 | 4.14 | 4.29 | 4.30 | 3.48 |
| Baysor (M prior 1.0) | 3.97 | 4.03 | 4.13 | 4.23 | 4.24 | 3.64 |
| Baysor (10x prior 1.0) | 3.97 | 4.08 | 4.15 | 4.27 | 4.28 | 3.52 |

Voronoi methods are nearly flat across reference resolution (−0.09 from 0.3→2.0) — their centroids already sit near reference centroids regardless of granularity. Baysor methods drop sharply at ref_res 2.0 (−0.78 for Baysor, −0.33 for Baysor M prior), meaning their centroids sit in gaps between reference cell states and only a dense enough reference grid (26 clusters) has a centroid nearby. This asymmetry — Baysor benefits from finer reference but not from finer self-resolution — indicates that Baysor's density-adaptive boundaries shift enough transcripts to displace population centroids away from where scRNA-seq places them, and this displacement cannot be corrected by clustering at a different resolution.

### Cell type centroid distance

As a complementary measure, each cell in both the scRNA-seq reference and every segmentation method is annotated with a cell type using `sc.tl.score_genes` on the same canonical marker panels, then per-cell-type centroids are computed in the shared 30-PC space. This bypasses Leiden resolution entirely and directly asks: for a given cell type, how close is the method's centroid to the reference's?

| Method | Mean | Median | Max |
| --- | ---: | ---: | ---: |
| Voronoi (M) | 3.49 | 3.44 | 4.86 |
| 10x native | 3.51 | 3.38 | 4.96 |
| Voronoi (10x) | 3.52 | 3.39 | 4.92 |
| Voronoi (CP) | 3.55 | 3.51 | 4.85 |
| Voronoi (SD) | 3.55 | 3.44 | 4.96 |
| Baysor (M prior 1.0) | 4.23 | 3.94 | 5.58 |
| Baysor (CP prior 1.0) | 4.29 | 4.00 | 5.64 |
| Baysor (10x prior 1.0) | 4.35 | 4.05 | 5.73 |
| Baysor (SD prior 1.0) | 4.39 | 4.09 | 5.80 |
| Baysor | 4.66 | 4.32 | 6.32 |

<details>
<summary><b>Per-cell-type breakdown</b> — click to expand</summary>

| Cell type | 10x native | Voronoi (CP) | Voronoi (SD) | Voronoi (M) | Voronoi (10x) | Baysor | Baysor (CP) | Baysor (SD) | Baysor (M) | Baysor (10x) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Smooth muscle | 2.34 | 2.12 | 2.33 | 2.18 | 2.27 | 4.05 | 3.61 | 3.76 | 3.51 | 3.72 |
| Adipocytes | 3.20 | 3.74 | 3.27 | 3.63 | 3.38 | 2.99 | 2.96 | 2.92 | 2.91 | 2.90 |
| B cells | 2.96 | 3.40 | 3.01 | 3.07 | 3.03 | 4.03 | 3.66 | 3.72 | 3.58 | 3.77 |
| CAFs | 3.41 | 2.91 | 3.45 | 3.01 | 3.31 | 5.45 | 5.04 | 5.21 | 4.94 | 5.13 |
| Endothelial | 3.62 | 3.37 | 3.69 | 3.39 | 3.61 | 5.57 | 4.97 | 5.15 | 4.85 | 5.10 |
| Luminal epithelial | 4.96 | 4.85 | 4.96 | 4.86 | 4.92 | 6.32 | 5.64 | 5.80 | 5.57 | 5.73 |
| Macrophages | 3.28 | 3.37 | 3.34 | 3.34 | 3.28 | 4.31 | 3.99 | 4.14 | 3.94 | 4.07 |
| Myoepithelial | 3.50 | 3.72 | 3.53 | 3.67 | 3.55 | 4.34 | 4.01 | 4.03 | 3.94 | 4.03 |
| Plasma cells | 3.35 | 3.62 | 3.42 | 3.50 | 3.39 | 3.72 | 3.45 | 3.57 | 3.49 | 3.56 |
| T cells | 4.49 | 4.38 | 4.52 | 4.29 | 4.43 | 5.79 | 5.54 | 5.60 | 5.58 | 5.44 |

</details>

Voronoi methods and 10x native are tightly grouped (mean 3.49–3.55), with Voronoi (M) marginally closest to the reference. Baysor methods are consistently ~0.8 units further (mean 4.23–4.66), with Baysor (M prior 1.0) closest among them. The gap is driven by CAFs, endothelial, and luminal epithelial — populations where Baysor's density-adaptive boundaries assign transcripts differently enough to shift the population centroid in PCA space. Adipocytes are the one cell type where Baysor methods sit closer to the reference than Voronoi (2.90 vs 3.20–3.74), possibly because adipocytes' diffuse morphology is better captured by transcript-density expansion than by Voronoi tessellation around compact nuclei.

### Hungarian-matched cluster centroids

The nearest-centroid tables above allow many-to-one matching — multiple method clusters can snap to the same reference centroid. Hungarian (one-to-one) matching enforces a strict bijection, giving a more principled measure of structural alignment. At each Leiden resolution (same for both sides), cluster centroids are matched one-to-one between each method and the scRNA-seq reference, between each method and 10x native, and pairwise between methods.

**Centroid distance to 10x native** (Hungarian matching, Euclidean in 30-PC reference space):

| Method | 0.3 | 0.5 | 0.7 | 1.0 | 1.5 | 2.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Voronoi (SD) | 0.34 | 0.22 | 0.75 | 0.40 | 0.73 | 0.69 |
| Voronoi (10x) | 0.39 | 0.68 | 0.80 | 0.46 | 0.49 | 0.55 |
| Voronoi (M) | 0.73 | 1.02 | 1.00 | 0.92 | 0.78 | 0.82 |
| Voronoi (CP) | 0.92 | 0.87 | 1.19 | 0.96 | 0.88 | 1.10 |
| Baysor (M prior 1.0) | 1.87 | 1.85 | 2.13 | 1.90 | 2.08 | 1.97 |
| Baysor (CP prior 1.0) | 2.06 | 2.17 | 2.27 | 2.39 | 2.30 | 2.19 |
| Baysor (10x prior 1.0) | 2.20 | 2.21 | 2.45 | 2.20 | 2.42 | 2.27 |
| Baysor (SD prior 1.0) | 2.29 | 2.25 | 2.49 | 2.17 | 2.43 | 2.32 |
| Baysor | 2.44 | 2.51 | 2.66 | 2.77 | 2.63 | 2.54 |

Voronoi cluster centroids are very close to 10x native's (0.22–1.19) across all resolutions, with Voronoi (SD) closest. Baysor centroids are 2–3× further away (1.78–2.77), confirming the population-level displacement seen in the cell-type analysis. The gap is stable across resolutions — it reflects a structural difference in where methods place cells in PCA space, not a resolution artifact.

**Centroid distance to 10x native** (argmax matching — each method cluster mapped to the 10x native cluster containing the plurality of its spatially matched cells):

| Method | 0.3 | 0.5 | 0.7 | 1.0 | 1.5 | 2.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Voronoi (SD) | 0.34 | 0.22 | 0.81 | 0.41 | 0.69 | 0.55 |
| Voronoi (10x) | 0.39 | 0.68 | 0.89 | 0.52 | 0.52 | 0.47 |
| Voronoi (M) | 0.73 | 1.07 | 0.95 | 0.92 | 0.74 | 0.81 |
| Voronoi (CP) | 0.92 | 0.87 | 1.26 | 1.04 | 0.95 | 1.02 |
| Baysor (M prior 1.0) | 2.59 | 2.84 | 3.03 | 2.92 | 3.07 | 3.21 |
| Baysor (CP prior 1.0) | 2.73 | 2.73 | 2.98 | 3.09 | 3.16 | 3.14 |
| Baysor (10x prior 1.0) | 2.91 | 2.62 | 2.94 | 3.04 | 3.04 | 3.20 |
| Baysor (SD prior 1.0) | 2.90 | 2.99 | 3.16 | 3.13 | 3.14 | 3.23 |
| Baysor | 3.04 | 3.19 | 3.32 | 3.33 | 3.27 | 3.44 |

Argmax matching uses cell overlap to determine which clusters correspond, allowing many-to-one mapping (multiple method clusters can map to the same 10x native cluster). The Voronoi distances are nearly identical to Hungarian (0.22–1.26), confirming that their clusters have clean one-to-one correspondence with 10x native. Baysor distances are ~30% larger under argmax than Hungarian (2.59–3.44 vs 1.78–2.77), indicating that when clusters are matched by actual cell content rather than geometric proximity, Baysor's centroids are even further from where 10x native places them. Baysor's many fine-grained clusters (22–37 at res 1.0–2.0) map many-to-one onto 10x native's coarser clusters, and the surplus clusters sit far from their assigned target centroids.

**ARI vs 10x native** (cell-level spatial matching, same reference PCA space):

| Method | 0.3 | 0.5 | 0.7 | 1.0 | 1.5 | 2.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Voronoi (10x) | 0.866 | 0.761 | 0.713 | 0.728 | 0.576 | 0.543 |
| Voronoi (M) | 0.820 | 0.698 | 0.705 | 0.670 | 0.572 | 0.529 |
| Voronoi (CP) | 0.771 | 0.768 | 0.747 | 0.674 | 0.519 | 0.472 |
| Voronoi (SD) | 0.763 | 0.745 | 0.694 | 0.689 | 0.520 | 0.518 |
| Baysor (10x prior 1.0) | 0.758 | 0.635 | 0.631 | 0.578 | 0.449 | 0.327 |
| Baysor (M prior 1.0) | 0.753 | 0.650 | 0.629 | 0.600 | 0.440 | 0.402 |
| Baysor (CP prior 1.0) | 0.750 | 0.627 | 0.623 | 0.531 | 0.423 | 0.343 |
| Baysor (SD prior 1.0) | 0.748 | 0.637 | 0.631 | 0.607 | 0.428 | 0.394 |
| Baysor | 0.365 | 0.397 | 0.386 | 0.365 | 0.277 | 0.215 |

ARI tracks centroid distance inversely: methods with closer centroids have higher cell-level agreement. All methods converge toward high ARI at low resolution (0.3, where 5–6 clusters make agreement easier) and diverge at high resolution where fine-grained cluster boundaries amplify method-specific boundary differences. Baysor PSC=1.0 variants start comparable to Voronoi at resolution 0.3 (ARI ~0.75) but decay faster, reaching 0.33–0.40 at resolution 2.0 vs 0.47–0.54 for Voronoi.

<details>
<summary><b>Pairwise centroid distances (resolution 1.0)</b> — click to expand</summary>

|  | 10x native | V (CP) | V (SD) | V (M) | V (10x) | Baysor | B (CP) | B (SD) | B (M) | B (10x) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10x native | — | 0.96 | 0.40 | 0.92 | 0.46 | 2.77 | 2.39 | 2.17 | 1.90 | 2.20 |
| V (CP) | 0.96 | — | 0.88 | 0.66 | 0.93 | 3.35 | 2.98 | 2.77 | 2.47 | 2.80 |
| V (SD) | 0.40 | 0.88 | — | 0.99 | 0.46 | 2.61 | 2.25 | 2.03 | 1.74 | 2.09 |
| V (M) | 0.92 | 0.66 | 0.99 | — | 0.95 | 3.28 | 2.92 | 2.70 | 2.43 | 2.74 |
| V (10x) | 0.46 | 0.93 | 0.46 | 0.95 | — | 2.84 | 2.46 | 2.24 | 1.94 | 2.30 |
| Baysor | 2.77 | 3.35 | 2.61 | 3.28 | 2.84 | — | 0.63 | 0.57 | 0.72 | 0.60 |
| B (CP) | 2.39 | 2.98 | 2.25 | 2.92 | 2.46 | 0.63 | — | 0.41 | 0.46 | 0.43 |
| B (SD) | 2.17 | 2.77 | 2.03 | 2.70 | 2.24 | 0.57 | 0.41 | — | 0.52 | 0.46 |
| B (M) | 1.90 | 2.47 | 1.74 | 2.43 | 1.94 | 0.72 | 0.46 | 0.52 | — | 0.63 |
| B (10x) | 2.20 | 2.80 | 2.09 | 2.74 | 2.30 | 0.60 | 0.43 | 0.46 | 0.63 | — |

</details>

The pairwise matrix shows clear block structure: within-Voronoi distances are 0.40–0.99, within-Baysor 0.41–0.72, and cross-family 1.74–3.35. Voronoi and Baysor occupy distinct regions of the reference PCA space. Baysor (M prior 1.0) is the Baysor variant closest to the Voronoi family (1.74–2.47), consistent with Mesmer's larger nuclear masks anchoring more transcripts and pulling the density model's centroids toward morphology-based positions.

---

## Pairwise method agreement

<p align="center"><img src="results/figures/pairwise_consensus.png" alt="Pairwise ARI between all segmentation methods" width="600"></p>

No, at least not within the Voronoi family. CellPose and StarDist agree with each other at ARI 0.764 (higher than the Voronoi pair at 0.661) because both are nuclear-morphology methods on the same DAPI image. Switching to Voronoi assignment lowers within-paradigm agreement because the two Voronoi variants use different centroids, shifting boundaries even where centroids are close. What Voronoi does raise is agreement with the 10x-native whole-cell reference (0.63-0.69): compatibility with the platform's own segmentation, not cross-method reproducibility. Baysor remains isolated from all morphological methods (ARI 0.30-0.46 regardless of partner).

---

## Cell size and disagreement

<details>
<summary><b>Cell size vs. disagreement (Hungarian)</b> — click to expand</summary>

![Cell size vs. disagreement - Hungarian](results/figures/cell_size_disagreement.png)

</details>

<details>
<summary><b>Cell size vs. disagreement (Argmax)</b> — click to expand</summary>

![Cell size vs. disagreement - argmax](results/figures/cell_size_disagreement_argmax.png)

</details>

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

<p align="center"><img src="results/figures/marker_recovery.png" alt="Marker gene recovery relative to 10x native" width="600"></p>

Using 10x-native cell-type annotations as ground truth, nuclear methods recover 75-92% of cytoplasmic marker expression relative to 10x native, with the largest deficits for extranuclear markers like MUC1, SERPINA3, and LYZ. Voronoi methods recover near-100% across all cell types. Baysor recovers macrophage markers (LYZ, CD14) at or above 10x-native levels while showing slightly reduced T cell marker (CD3E) recovery.

---

## Population-level convergence

<p align="center"><img src="results/figures/pseudobulk_correlation.png" alt="Pseudobulk correlation vs. 10x native" width="600"></p>

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
