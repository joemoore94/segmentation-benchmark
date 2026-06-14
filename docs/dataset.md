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
`morphology_focus.ome.tif`, so the two appear to already share the same pixel grid;
no separate registration step was needed to crop matching ROIs from both.

## Image format gotchas

These cost real debugging time and matter for any code that touches the OME-TIFFs:

1. **`morphology_focus.ome.tif` is single-channel DAPI only** (shape `53833 x 48441`,
   `uint16`, axes `YX`, OME channel name `"DAPI"`), *not* multi-channel
   DAPI+membrane/boundary as the original plan assumed. CellPose/Mesmer segmentation
   on this dataset is therefore nuclear/DAPI-based only; Mesmer's `--membrane-image`
   flag is omitted (nuclear-only or `--compartment nuclear`).
2. **`morphology_focus.ome.tif` is JPEG2000-compressed** (TIFF `Compression=34712`,
   tiled 1024x1024). `pyvips`/libvips' libtiff build can't decode this: it silently
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
4. **`spatialdata_io.xenium()` doesn't work with this file set**: it requires the
   Xenium Explorer `*.zarr.zip` bundles (`cells.zarr.zip`, `transcripts.zarr.zip`,
   `cell_feature_matrix.zarr.zip`, `analysis.zarr.zip`), which aren't part of the
   per-file download used here (we have the `.parquet`/`.h5`/`.ome.tif` equivalents
   instead). [`io.py`](../src/segbench/io.py) therefore reads the raw files directly
   rather than going through `spatialdata`/`spatialdata_io`.
5. **Julia's `Parquet.jl` can't read pandas/pyarrow-written parquet files**: Baysor
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
(visible as a black block in the DAPI crop / blank region in H&E); the tissue-coverage
constraint above was added specifically to rule out windows like this. The chosen ROI's
H&E shows several rounded tumor (DCIS/invasive) nests surrounded by stroma and immune
infiltrate, good heterogeneity for segmentation comparison.

Extracted via `segbench.io.extract_roi()` into `data/processed/roi/`:
`dapi.tif`, `he.tif`, `transcripts.parquet`, `cells.parquet`, `cell_boundaries.parquet`,
`nucleus_boundaries.parquet`, `cell_feature_matrix.h5ad` (all coordinates shifted so the
ROI's top-left corner is the origin).

## Segmentation runs

All segmentation runs on this project are on a Linux x86_64 workstation.
CellPose and Mesmer segment the DAPI image (`dapi.tif`, 9412x9412px at
`pixel_size=0.2125`); Baysor segments transcripts directly
(`transcripts_baysor.csv`/`.parquet`, the qv>=20 non-control-probe subset,
3,392,051 of 4,360,872 rows).

### CellPose

- CellPose 3.1.1.3. CellPose 4.x dropped the lightweight classical U-Net
  models (`nuclei`, `cyto2`, `cyto3`) in favor of SAM-based foundation models,
  which are CPU-prohibitive; 3.x's `nuclei` model runs comfortably on CPU.
- Nuclear segmentation on the full `dapi.tif` (9412x9412px, 2mm x 2mm ROI), CPU
  only (`scripts/run_cellpose_roi.py`, `model_type="nuclei"`).
- Runtime: 633.4s (~10.5 min).
- Result: **20,166 cells** -> `masks_cellpose.tif`, median nucleus area
  653 px² (~29.5 µm²).

### Baysor

- Julia 1.10 via juliaup (`julia +1.10`, see Environment setup); Baysor
  v0.7.1 isn't compatible with the default Julia channel.
- Baysor's main-EM runtime scales worse than linearly with transcript count
  (roughly N^1.8), so the full ROI (3,392,051 transcripts) is run as **4
  overlapping ~1mm x 1mm tiles** rather than one job.
  `scripts/tile_baysor_transcripts.py` splits the ROI into a 2x2 grid of
  quadrants, each padded by 50um on its interior edges so Baysor sees full
  local context near tile boundaries. Config in `configs/baysor_config.toml`
  (`scale=12.5`, `scale_std="25%"`, `n_clusters=4`) is unchanged per tile --
  the scale parameters are physical and area-independent.
- Runtime: ~12-16 min main run per tile, ~60 min total for all 4 (after a
  one-time ~8 min Julia precompilation on the first run).
- `scripts/build_baysor_adata.py` merges the 4 tiles' `segmentation.csv`
  outputs into `adata_baysor.h5ad`: for each tile, keeps only cells whose
  centroid (mean molecule position) falls in that tile's non-padded "core"
  region. The 4 cores exactly partition the ROI, so every cell is counted
  once while still having its full set of assigned molecules (including any
  in the tile's padding). Transcripts assigned to a kept cell by two tiles'
  overlapping runs are deduplicated by `transcript_id` (4,711 of ~3.73M
  tile-molecule rows).
- Result: **18,321 cells**, 3,344,675 / 3,392,051 transcripts assigned
  (98.6% capture rate) -> `baysor_tiles/<tile>/segmentation*.{csv,loom,json}`,
  merged into `adata_baysor.h5ad`.
- Gene-name gotcha: `transcripts_baysor*.csv`'s `gene` column is written from
  a `bytes`-typed `feature_name` column, so values come out as literal
  `"b'GENENAME'"` strings. Cleaned with a regex (`r"^b'(.*)'$"` -> `r"\1"`)
  in `build_baysor_adata.py`; without this, CellPose and Baysor shared 0/379
  genes.

### Baysor (CellPose-prior hybrid)

- Same transcript-density EM as above, but seeded with CellPose's nucleus
  masks as a `--prior-segmentation` (Baysor's `:column_name` syntax).
  `scripts/add_cellpose_prior.py` looks up each transcript's
  `(x_location, y_location)` in `masks_cellpose.tif` (converting um -> pixel
  via `PIXEL_SIZE=0.2125`) and writes a `cellpose_prior` column (0 =
  background, matching Baysor's `unassigned_prior_label="0"` default);
  1,202,133 / 3,392,051 transcripts (35.4%) fall within a CellPose nucleus.
- `configs/baysor_prior_config.toml` adds `prior_segmentation_confidence =
  0.2` (Baysor's default, range [0,1]) on top of the same `scale=12.5,
  scale_std="25%", n_clusters=4` config as the non-prior run. Because `scale`
  is set explicitly, Baysor does not auto-estimate scale from the prior's
  nucleus sizes (`estimate_scale_from_centers` only fires when `scale<=0`),
  so the prior nudges the transcript-density EM toward CellPose's nuclei
  without rescaling the whole model toward nuclear dimensions.
- Same 4-tile (2x2, 50um padding) tiling/merge scheme as the non-prior run
  (`scripts/tile_baysor_prior_transcripts.py`,
  `scripts/build_baysor_prior_adata.py`), with `:cellpose_prior` passed as
  the 4th arg to `scripts/run_baysor.sh`.
- Runtime: comparable to the non-prior run, ~13-18 min main run per tile.
- Result: **19,061 cells**, 3,346,292 / 3,392,051 transcripts assigned
  (98.7% capture rate) -> merged into `adata_baysor_prior.h5ad`.

### Mesmer

- Not run, blocked externally on deepcell.org, not on anything in this repo.
  A native `mesmer` conda env (Python 3.10, DeepCell 0.12.10, TensorFlow
  2.8.4, no Docker needed on this x86_64 machine) is set up, and
  `segbench.segmentation.mesmer_run.run_mesmer` / `scripts/run_mesmer.py`
  call `deepcell.applications.Mesmer()` directly via `conda run -n mesmer`.
- `Mesmer()` requires a `DEEPCELL_ACCESS_TOKEN` to fetch its pretrained
  weights (`models/MultiplexSegmentation-9.tar.gz`, md5
  `a1dfbce2594f927b9112f23a0a1739e0`) from `https://users.deepcell.org/api/getData/`.
  As of June 2026, `users.deepcell.org`'s account system is broken end to
  end: signup returns a server error (HTTP 500, tried with two different
  emails), login fails ("username and password didn't match") even with the
  documented `email-local-part` username convention, and "forgot password"
  never sends a reset email. No public mirror of the weights archive (which
  would let `fetch_data`'s local md5-cache check succeed without a token) was
  found. The DeepCell.org "PREDICT" cloud page is a token-free alternative
  but requires browser upload, which was ruled out in favor of staying CLI-only.
- Plan if/when access is restored: nuclear-only segmentation on `dapi.tif`
  (no membrane channel, see image-format gotcha #1), on the same 1mm x 1mm
  sub-region used for Baysor to keep its compute footprint small and its
  results directly comparable.

### Cross-method comparison scope

All five methods (CellPose, Baysor, 10x native, StarDist, and the
CellPose-prior Baysor hybrid) now cover the same full 2mm x 2mm ROI, so cell
counts, size distributions, and transcript capture rates are directly
comparable without density normalization. `match_cells_by_centroid` (max
10 µm) pairs 8,947 CellPose/Baysor cells, 18,966 CellPose/10x native cells,
19,460 CellPose/StarDist cells, and 9,572 CellPose/Baysor(prior) cells; all
matched-pair metrics (expression correlation, cell-type agreement, spatial
disagreement) are computed over these full-ROI matched sets. See
[`../README.md`](../README.md#results) for results.
