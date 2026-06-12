# Dataset notes

## Source

**Xenium FFPE Human Breast (Custom Add-on Panel)**, Janesick et al. 2023
(*Nat Commun*), 10x Genomics dataset page:
https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-with-custom-add-on-panel-1-standard

The dataset landing page is JS-rendered, so direct file links aren't easily scraped.
To get the actual download URLs:

1. Open the dataset page in a browser.
2. Find the "Output and supplemental files" section and copy the links for:
   - Xenium Output Bundle (`*_outs.zip` or the individual files below, if listed
     separately)
   - `morphology_focus` OME-TIFF(s) (multi-channel: DAPI + boundary/membrane markers)
   - `transcripts.parquet`
   - `cell_feature_matrix.h5`
   - `cells.parquet` (centroids/metadata)
   - `cell_boundaries.parquet` / `nucleus_boundaries.parquet` (10x native segmentation)
   - registered post-Xenium H&E image + alignment file
   - `experiment.xenium` (run metadata)
3. Record the URLs in [`scripts/download_data.sh`](../scripts/download_data.sh) and
   run it to populate `data/raw/`.

## Matched scRNA-seq / Visium (for Project 2)

Same study, GEO accession [GSE243275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243275).
Companion analysis code: https://github.com/10XGenomics/janesick_nature_comms_2023_companion

## ROI selection (TODO)

Once the full bundle is downloaded, pick a ~2-4mm² region containing a mix of tumor,
stroma, and immune infiltrate (e.g. by overlaying the H&E + 10x cell type annotations
from the companion repo's label-transfer results). Record the chosen ROI bounding box
(in the Xenium coordinate system) here once selected.
