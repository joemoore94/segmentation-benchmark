#!/usr/bin/env bash
# Run Mesmer (DeepCell) segmentation via the vanvalenlab/deepcell-applications
# Docker image. No native arm64 build exists (TF2.8-based image is x86_64
# only), so it runs under emulation via --platform linux/amd64.
#
# Usage: run_mesmer.sh <data_dir> <nuclear_file> <membrane_file> <output_dir> [compartment]
#
# <nuclear_file>, <membrane_file>, and <output_dir> are paths relative to
# <data_dir>, which is bind-mounted into the container at /data.
# [compartment] is one of: whole-cell (default), nuclear.

set -euo pipefail

DATA_DIR=$1
NUCLEAR_FILE=$2
MEMBRANE_FILE=$3
OUTPUT_DIR=$4
COMPARTMENT=${5:-whole-cell}

MOUNT_DIR=/data

mkdir -p "${DATA_DIR}/${OUTPUT_DIR}"

docker run --rm --platform linux/amd64 \
  -v "${DATA_DIR}:${MOUNT_DIR}" \
  vanvalenlab/deepcell-applications:latest \
  mesmer \
  --nuclear-image "${MOUNT_DIR}/${NUCLEAR_FILE}" \
  --membrane-image "${MOUNT_DIR}/${MEMBRANE_FILE}" \
  --output-directory "${MOUNT_DIR}/${OUTPUT_DIR}" \
  --compartment "${COMPARTMENT}"
