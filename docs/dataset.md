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
| `morphology_focus.ome.tif`   | 1.5 GB | Single-channel DAPI image for CellPose/Mesmer (see gotchas below)    |
| `transcripts.parquet`        | 1.9 GB | Per-transcript detections for Baysor                                 |
| `cell_feature_matrix.h5`     | 38 MB  | 10x quantified cell x gene matrix                                    |
| `cells.parquet`              | 12 MB  | 10x cell centroids/metadata                                          |
| `cell_boundaries.parquet`    | 31 MB  | 10x native cell segmentation (reference)                             |
| `nucleus_boundaries.parquet` | 30 MB  | 10x native nucleus segmentation (reference)                          |
| `experiment.xenium`          | 1.6 KB | Run metadata                                                         |
| `metrics_summary.csv`        | <1 KB  | Run QC metrics                                                       |
| `he_image.ome.tif`           | 5.0 GB | Paired H&E, deconvolved (Hematoxylin/Eosin/Residual)                 |

Note: no standalone alignment-matrix file was found for this sample (only for the
separate `Rep1` pre-designed-panel sample). However, `he_image.ome.tif` has the same
`SizeX`/`SizeY` (48441 x 53833) and `PhysicalSizeX/Y` (0.2125 um) as
`morphology_focus.ome.tif`, so the two appear to already share the same pixel grid —
no separate registration step was needed to crop matching ROIs from both.

## Image format gotchas

These cost real debugging time and matter for any code that touches the OME-TIFFs:

1. **`morphology_focus.ome.tif` is single-channel DAPI only** (shape `53833 x 48441`,
   `uint16`, axes `YX`, OME channel name `"DAPI"`) — *not* multi-channel
   DAPI+membrane/boundary as the original plan assumed. CellPose/Mesmer segmentation
   on this dataset is therefore nuclear/DAPI-based only; Mesmer's `--membrane-image`
   flag is omitted (nuclear-only or `--compartment nuclear`).
2. **`morphology_focus.ome.tif` is JPEG2000-compressed** (TIFF `Compression=34712`,
   tiled 1024x1024). `pyvips`/libvips' libtiff build can't decode this — it silently
   returns an all-zero image (`dapi.avg() == 0`, no error) instead of failing loudly.
   `tifffile` (via `imagecodecs`) decodes it correctly. [`io.py`](../src/segbench/io.py)
   reads it through `tifffile.imread(..., aszarr=True)` + `zarr`, which also allows
   cropping an ROI without decoding the full ~5.2 GB array.
3. **`he_image.ome.tif` stores its 3 channels (Hematoxylin/Eosin/Residual) as 3
   separate TIFF *pages*, not as 3 bands of one page** (compression=8, deflate, which
   pyvips reads fine). Loading with `pyvips.Image.new_from_file(path, n=3)` produces a
   "toilet roll" image (`width x height*3`, 1 band) rather than an RGB image.
   [`io.py`](../src/segbench/io.py) loads each page separately
   (`page=0,1,2`) and `bandjoin`s them into an `(H, W, 3)` array.
4. **`spatialdata_io.xenium()` doesn't work with this file set** — it requires the
   Xenium Explorer `*.zarr.zip` bundles (`cells.zarr.zip`, `transcripts.zarr.zip`,
   `cell_feature_matrix.zarr.zip`, `analysis.zarr.zip`), which aren't part of the
   per-file download used here (we have the `.parquet`/`.h5`/`.ome.tif` equivalents
   instead). [`io.py`](../src/segbench/io.py) therefore reads the raw files directly
   rather than going through `spatialdata`/`spatialdata_io`.
5. **Julia's `Parquet.jl` can't read pandas/pyarrow-written parquet files** — Baysor
   (via `Parquet.jl`) throws inside its Thrift metadata reader on
   `transcripts_baysor.parquet`. CSV is Baysor's most robust input format, so the ROI
   transcript table is also written as `transcripts_baysor.csv` and that's what's
   passed to `run_baysor.sh`.

## Matched scRNA-seq / Visium (for Project 2)

Same study, GEO accession [GSE243275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243275).
Companion analysis code: <https://github.com/10XGenomics/janesick_nature_comms_2023_companion>

## ROI selection

**Chosen ROI**: `x = [5600, 7600] um`, `y = [8000, 10000] um` (2mm x 2mm, "global"
Xenium coordinate system; pixel bounds `x_px = [26353, 35765]`, `y_px = [37647, 47059]`
at `pixel_size = 0.2125 um/px`).

Selection method (ad hoc analysis script, not checked in):

1. Binned all 576,963 cells into a 200um x 200um grid using `cells.parquet` +
   `cell_feature_matrix.h5`, summing per-bin counts for three marker panels:
   - epithelial/tumor: `EPCAM`, `KRT8`, `KRT14`, `KRT5`, `ERBB2`, `ESR1`
   - immune: `PTPRC`, `CD3E`, `CD8A`, `CD4`, `CD68`, `ITGAM`, `MS4A1`, `CD19`, `CD14`
   - stroma: `PDGFRB`, `ACTA2`, `PECAM1`, `VWF`
2. Computed a tissue-coverage mask from the smallest `morphology_focus.ome.tif`
   pyramid level (level 7, 420x378px, ~27.2 um/px) to exclude slide background.
3. Slid a 10x10-bin (2mm x 2mm) window over the grid, keeping only windows where
   every bin has >=95% tissue coverage and the window has >=5000 cells total, then
   ranked by `balance = min(epi_frac, imm_frac, str_frac)` (the per-window fraction
   of total epi+imm+stroma marker counts in each category).
4. Top candidate: `n_cells = 23629`, `balance = 0.291`
   (`epi_frac = 0.40`, `imm_frac = 0.29`, `str_frac = 0.30`), 100% tissue coverage.

A first candidate (`x=[200,2200], y=[7400,9400]`, `balance=0.312`) had a higher balance
score but ~40% of the window fell on empty slide background at the tissue edge
(visible as a black block in the DAPI crop / blank region in H&E) — the tissue-coverage
constraint above was added specifically to rule out windows like this. The chosen ROI's
H&E shows several rounded tumor (DCIS/invasive) nests surrounded by stroma and immune
infiltrate — good heterogeneity for segmentation comparison.

Extracted via `segbench.io.extract_roi()` into `data/processed/roi/`:
`dapi.tif`, `he.tif`, `transcripts.parquet`, `cells.parquet`, `cell_boundaries.parquet`,
`nucleus_boundaries.parquet`, `cell_feature_matrix.h5ad` (all coordinates shifted so the
ROI's top-left corner is the origin).

## Segmentation runs

CellPose and Mesmer segment the DAPI image (`dapi.tif`, 9412x9412px at
`pixel_size=0.2125`); Baysor segments transcripts directly
(`transcripts_baysor.csv`/`.parquet`, the qv>=20 non-control-probe subset,
3,392,051 of 4,360,872 rows).

### CellPose

- Nuclear segmentation on `dapi.tif`, MPS-accelerated (`scripts/run_cellpose_roi.py`).
- Runtime: ~94 min (5625.6s).
- Result: 23,660 cells -> `masks_cellpose.tif`.

### Mesmer

- Nuclear-only segmentation on `dapi.tif` via the `vanvalenlab/deepcell-applications`
  Docker image (no membrane channel available, see image-format gotcha #1).
- Runtime: ~58 min.
- Result: 21,693 cells -> `mesmer_out/mask.tif`.

### Baysor

- Transcript-based segmentation, config in `configs/baysor_config.toml`
  (`scale=12.5`, `scale_std="25%"`, `n_clusters=4`).
- **Smoke test** on a 500x500um corner (175,411 transcripts, 5.17% of the ROI's
  transcripts): completed in ~7m48s, 1,155 cells -> `baysor_test_out/`.
- **Full-ROI run** (3,392,051 transcripts) was killed after ~75 min, at iteration
  104/500 of the main EM step (21% progress, ETA still >4 hours). Baysor's main-EM
  runtime scales worse than linearly with molecule count (roughly N^1.8, based on
  the smoke-test vs full-ROI iteration timing), so the full 2mm x 2mm ROI is
  impractical for this method on this hardware.
- **Plan**: rerun Baysor on a smaller, centered 1mm x 1mm sub-region of the same ROI
  (~850K transcripts, est. ~35-40 min by the same scaling), and restrict any
  cross-method comparison involving Baysor to that sub-region.

### Resource notes

Running CellPose (MPS), Mesmer (Docker/x86_64 emulation via Rosetta), and Baysor
(Julia) concurrently on a 16GB M1 caused severe swap thrashing (up to 14.3/15GB
swap used), stalling CellPose and Baysor at near-0% CPU for extended periods.
Quitting Docker Desktop after Mesmer finished freed enough memory for both to
resume normal progress.
