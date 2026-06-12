#!/usr/bin/env bash
# Run Mesmer (DeepCell) segmentation via the vanvalenlab/deepcell-applications
# Docker image (CPU build - no native Apple Silicon deepcell-tf install exists).
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

docker run --rm \
  -v "${DATA_DIR}:${MOUNT_DIR}" \
  vanvalenlab/deepcell-applications:latest-cpu \
  mesmer \
  --nuclear-image "${MOUNT_DIR}/${NUCLEAR_FILE}" \
  --membrane-image "${MOUNT_DIR}/${MEMBRANE_FILE}" \
  --output-directory "${MOUNT_DIR}/${OUTPUT_DIR}" \
  --compartment "${COMPARTMENT}"
