#!/usr/bin/env bash
# Run Mesmer (DeepCell) segmentation via the vanvalenlab/deepcell-applications
# Docker image. The image bundles pretrained model weights and does not require
# a DEEPCELL_ACCESS_TOKEN, bypassing the auth requirement added in the pip
# package.
#
# Usage:
#   bash scripts/run_mesmer.sh <data_dir> <nuclear_file> <output_dir> [compartment] [image_mpp]
#
#   data_dir      absolute path to data directory; mounted as /data in the container
#   nuclear_file  nuclear-channel TIFF, relative to data_dir
#   output_dir    output directory, relative to data_dir; mask.tif written here
#   compartment   nuclear (default) | whole-cell | both
#   image_mpp     microns per pixel (default: 0.2125)
#
# Example:
#   bash scripts/run_mesmer.sh \
#       data/processed/roi dapi.tif mesmer_out nuclear 0.2125

set -euo pipefail

DATA_DIR=$(realpath "${1}")
NUCLEAR_FILE="${2}"
OUTPUT_DIR="${3}"
COMPARTMENT="${4:-nuclear}"
IMAGE_MPP="${5:-0.2125}"

mkdir -p "${DATA_DIR}/${OUTPUT_DIR}"

docker run --rm \
    -v "${DATA_DIR}:/data" \
    vanvalenlab/deepcell-applications:latest \
    python /run_app.py mesmer \
    --nuclear "/data/${NUCLEAR_FILE}" \
    --output-dir "/data/${OUTPUT_DIR}" \
    --compartment "${COMPARTMENT}" \
    --image-mpp "${IMAGE_MPP}"

echo "Mesmer done -> ${DATA_DIR}/${OUTPUT_DIR}/mask.tif"
