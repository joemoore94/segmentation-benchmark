# Dataset notes

## Source

**Xenium FFPE Human Breast (Custom Add-on Panel)**, Janesick et al. 2023
(*Nat Commun*), 10x Genomics dataset page:
<https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-with-custom-add-on-panel-1-standard>

Internal 10x sample name: `Xenium_V1_FFPE_Human_Breast_IDC_With_Addon` (Xenium
ranger 1.3.0), hosted at
`https://cf.10xgenomics.com/samples/xenium/1.3.0/Xenium_V1_FFPE_Human_Breast_IDC_With_Addon/`.
The dataset landing page is JS-rendered, but the per-file URLs follow the
`<base>_<filename>` pattern below (resolved by fetching the page's raw HTML and
grepping for `cf.10xgenomics.com` links).

[`scripts/download_data.sh`](../scripts/download_data.sh) downloads the files below
into `data/raw/` (~8.4 GB total; the full `_outs.zip` bundle is ~27 GB and is not used):

| File                         | Size   | Use                                                                  |
| ---------------------------- | ------ | -------------------------------------------------------------------- |
| `morphology_focus.ome.tif`   | 1.5 GB | Multi-channel IF image (DAPI + boundary markers) for CellPose/Mesmer |
| `transcripts.parquet`        | 1.9 GB | Per-transcript detections for Baysor                                 |
| `cell_feature_matrix.h5`     | 38 MB  | 10x quantified cell x gene matrix                                    |
| `cells.parquet`              | 12 MB  | 10x cell centroids/metadata                                          |
| `cell_boundaries.parquet`    | 31 MB  | 10x native cell segmentation (reference)                             |
| `nucleus_boundaries.parquet` | 30 MB  | 10x native nucleus segmentation (reference)                          |
| `experiment.xenium`          | 1.6 KB | Run metadata                                                         |
| `metrics_summary.csv`        | <1 KB  | Run QC metrics                                                       |
| `he_image.ome.tif`           | 5.0 GB | Paired H&E (unregistered)                                            |

Note: a registered/aligned H&E + alignment matrix were not found as standalone
downloads for this sample (only for the separate `Rep1` pre-designed-panel sample).
If H&E/IF registration is needed, it may have to be done manually (e.g. with VALIS),
or `he_image.ome.tif` can be used for qualitative tissue-region context only.

## Matched scRNA-seq / Visium (for Project 2)

Same study, GEO accession [GSE243275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243275).
Companion analysis code: <https://github.com/10XGenomics/janesick_nature_comms_2023_companion>

## ROI selection (TODO)

Once the full bundle is downloaded, pick a ~2-4mm² region containing a mix of tumor,
stroma, and immune infiltrate (e.g. by overlaying the H&E + 10x cell type annotations
from the companion repo's label-transfer results). Record the chosen ROI bounding box
(in the Xenium coordinate system) here once selected.
