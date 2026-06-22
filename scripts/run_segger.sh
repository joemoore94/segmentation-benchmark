#!/usr/bin/env bash
# Run Segger (GNN multimodal segmentation) on the Xenium dataset.
#
# Segger uses a Graph Neural Network to combine nuclear masks (from the 10x
# native segmentation bundled in the Xenium output) with transcript point
# clouds, learning cell boundaries from both modalities.
#
# Requirements:
#   - CUDA GPU with sufficient VRAM (~8 GB)
#   - segger installed: pip install "git+https://github.com/dpeerlab/segger.git"
#   - torch-geometric, cupy, rmm (RAPIDS) installed in the segbench env
#
# Usage:
#   conda run -n segbench bash scripts/run_segger.sh
#
# After running, build the AnnData for the benchmark pipeline:
#   conda run -n segbench python scripts/build_segger_adata.py

set -euo pipefail

DATA_DIR="data/raw"
OUT_DIR="data/processed/segger"

mkdir -p "$OUT_DIR"

echo "=== Running Segger on Xenium data ==="
echo "Input:  $DATA_DIR"
echo "Output: $OUT_DIR"

conda run -n segbench python -m segger segment \
    -i "$DATA_DIR" \
    -o "$OUT_DIR" \
    --prediction-mode nucleus \
    --save-anndata \
    --n-epochs 20

echo "=== Segger complete ==="
echo "Check $OUT_DIR for output AnnData."
