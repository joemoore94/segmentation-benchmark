#!/usr/bin/env bash
# Download the Xenium FFPE Human Breast (Custom Add-on Panel) dataset
# (Janesick et al. 2023) into data/raw/.
#
# 10x's dataset page is JS-rendered, so the URLs below must be filled in
# manually from the "Output and supplemental files" section of:
# https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-with-custom-add-on-panel-1-standard
# See docs/dataset.md for details.

set -euo pipefail

RAW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data/raw"
mkdir -p "${RAW_DIR}"

# TODO: replace with real URLs from the dataset page.
MORPHOLOGY_URL="REPLACE_ME"
TRANSCRIPTS_URL="REPLACE_ME"
CELL_FEATURE_MATRIX_URL="REPLACE_ME"
CELLS_URL="REPLACE_ME"
CELL_BOUNDARIES_URL="REPLACE_ME"
NUCLEUS_BOUNDARIES_URL="REPLACE_ME"
HE_IMAGE_URL="REPLACE_ME"
EXPERIMENT_URL="REPLACE_ME"

for var in MORPHOLOGY_URL TRANSCRIPTS_URL CELL_FEATURE_MATRIX_URL CELLS_URL \
           CELL_BOUNDARIES_URL NUCLEUS_BOUNDARIES_URL HE_IMAGE_URL EXPERIMENT_URL; do
  if [[ "${!var}" == "REPLACE_ME" ]]; then
    echo "Error: ${var} is not set. See docs/dataset.md for how to find these URLs." >&2
    exit 1
  fi
done

curl -L -o "${RAW_DIR}/morphology_focus.ome.tif" "${MORPHOLOGY_URL}"
curl -L -o "${RAW_DIR}/transcripts.parquet" "${TRANSCRIPTS_URL}"
curl -L -o "${RAW_DIR}/cell_feature_matrix.h5" "${CELL_FEATURE_MATRIX_URL}"
curl -L -o "${RAW_DIR}/cells.parquet" "${CELLS_URL}"
curl -L -o "${RAW_DIR}/cell_boundaries.parquet" "${CELL_BOUNDARIES_URL}"
curl -L -o "${RAW_DIR}/nucleus_boundaries.parquet" "${NUCLEUS_BOUNDARIES_URL}"
curl -L -o "${RAW_DIR}/he_image.ome.tif" "${HE_IMAGE_URL}"
curl -L -o "${RAW_DIR}/experiment.xenium" "${EXPERIMENT_URL}"

echo "Downloaded raw Xenium bundle to ${RAW_DIR}"
