#!/usr/bin/env bash
# Run Baysor with prior_segmentation_confidence=1.0 on the same 4 tiles used
# for the c05/c08 runs, then merge into an AnnData.
#
# Usage: bash scripts/run_baysor_c10.sh
#        (run from repo root)

set -euo pipefail

ROI_DIR="data/processed/roi"
SRC_TILES="${ROI_DIR}/baysor_prior_tiles"
CONFIG="configs/baysor_prior_config_c10.toml"
OUT_TILES="${ROI_DIR}/baysor_prior_c10_tiles"

mkdir -p "${OUT_TILES}"
cp "${SRC_TILES}/manifest.json" "${OUT_TILES}/manifest.json"

echo "=== confidence=1.0: ${CONFIG} ==="
for TILE in x0_y0 x1_y0 x0_y1 x1_y1; do
    TRANSCRIPT_FILE="${SRC_TILES}/transcripts_${TILE}.csv"
    OUT_DIR="${OUT_TILES}/${TILE}"
    echo "  tile ${TILE} -> ${OUT_DIR}"
    bash scripts/run_baysor.sh "${TRANSCRIPT_FILE}" "${CONFIG}" "${OUT_DIR}" cellpose_prior
done

echo "  merging tiles..."
/home/joe/miniforge3/bin/conda run -n segbench \
    python scripts/build_baysor_conf_adata.py 1.0
echo "  done: c10"
