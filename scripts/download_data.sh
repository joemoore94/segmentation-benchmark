#!/usr/bin/env bash
# Download the Xenium FFPE Human Breast (Custom Add-on Panel) dataset
# (Janesick et al. 2023; "IDC With Addon", Xenium ranger 1.3.0) into data/raw/.
#
# Source: https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-with-custom-add-on-panel-1-standard
# See docs/dataset.md for details. Total download is ~8.4 GB.

set -euo pipefail

RAW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data/raw"
mkdir -p "${RAW_DIR}"

BASE="https://cf.10xgenomics.com/samples/xenium/1.3.0/Xenium_V1_FFPE_Human_Breast_IDC_With_Addon/Xenium_V1_FFPE_Human_Breast_IDC_With_Addon"

FILES="
morphology_focus.ome.tif
transcripts.parquet
cell_feature_matrix.h5
cells.parquet
cell_boundaries.parquet
nucleus_boundaries.parquet
experiment.xenium
metrics_summary.csv
he_image.ome.tif
"

for fname in ${FILES}; do
  dest="${RAW_DIR}/${fname}"
  if [[ -f "${dest}" ]]; then
    echo "Skipping ${fname} (already exists)"
    continue
  fi
  echo "Downloading ${fname} ..."
  curl -L --fail -o "${dest}.partial" "${BASE}_${fname}"
  mv "${dest}.partial" "${dest}"
done

echo "Downloaded raw Xenium bundle to ${RAW_DIR}"
