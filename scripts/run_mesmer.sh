#!/usr/bin/env bash
# Run Mesmer (DeepCell) segmentation via the vanvalenlab/deepcell-applications
# Docker image. No native arm64 build exists (TF2.8-based image is x86_64
# only), so it runs under emulation via --platform linux/amd64.
#
# Usage: run_mesmer.sh <data_dir> <nuclear_file> <output_dir> <compartment> <image_mpp> [membrane_file]
#
# <nuclear_file>, <membrane_file>, and <output_dir> are paths relative to
# <data_dir>, which is bind-mounted into the container at /data.
# <compartment> is one of: whole-cell, nuclear, both.
# [membrane_file] is optional; if omitted, the membrane channel input is left
# blank (our morphology_focus image is DAPI-only, so this is the normal case).

set -euo pipefail

DATA_DIR=$1
NUCLEAR_FILE=$2
OUTPUT_DIR=$3
COMPARTMENT=$4
IMAGE_MPP=$5
MEMBRANE_FILE=${6:-}

MOUNT_DIR=/data

mkdir -p "${DATA_DIR}/${OUTPUT_DIR}"

ARGS=(
  mesmer
  --nuclear-image "${MOUNT_DIR}/${NUCLEAR_FILE}"
  --output-directory "${MOUNT_DIR}/${OUTPUT_DIR}"
  --compartment "${COMPARTMENT}"
  --image-mpp "${IMAGE_MPP}"
  --squeeze
)
if [[ -n "${MEMBRANE_FILE}" ]]; then
  ARGS+=(--membrane-image "${MOUNT_DIR}/${MEMBRANE_FILE}")
fi

docker run --rm --platform linux/amd64 \
  -v "${DATA_DIR}:${MOUNT_DIR}" \
  vanvalenlab/deepcell-applications:latest \
  "${ARGS[@]}"
