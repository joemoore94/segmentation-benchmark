#!/usr/bin/env bash
# Run Baysor with prior_segmentation_confidence = 0.5 and 0.8 on the same
# 4 tiles and transcript files used for the confidence=0.2 baseline.
# Creates baysor_prior_c05_tiles/ and baysor_prior_c08_tiles/, each with
# the manifest.json from baysor_prior_tiles/ and 4 Baysor output subdirs.
#
# Usage: bash scripts/run_baysor_conf_sensitivity.sh
#        (run from repo root, no arguments)

set -euo pipefail

ROI_DIR="data/processed/roi"
SRC_TILES="${ROI_DIR}/baysor_prior_tiles"
TRANSCRIPTS_SRC="${ROI_DIR}/transcripts_baysor_prior.csv"

for CONF_LABEL in c05 c08; do
    if [[ "${CONF_LABEL}" == "c05" ]]; then
        CONFIG="configs/baysor_prior_config_c05.toml"
        CONF_FLOAT="0.5"
    else
        CONFIG="configs/baysor_prior_config_c08.toml"
        CONF_FLOAT="0.8"
    fi

    OUT_TILES="${ROI_DIR}/baysor_prior_${CONF_LABEL}_tiles"
    mkdir -p "${OUT_TILES}"
    cp "${SRC_TILES}/manifest.json" "${OUT_TILES}/manifest.json"

    echo "=== confidence=${CONF_LABEL}: ${CONFIG} ==="
    for TILE in x0_y0 x1_y0 x0_y1 x1_y1; do
        TRANSCRIPT_FILE="${SRC_TILES}/transcripts_${TILE}.csv"
        OUT_DIR="${OUT_TILES}/${TILE}"
        echo "  tile ${TILE} -> ${OUT_DIR}"
        bash scripts/run_baysor.sh "${TRANSCRIPT_FILE}" "${CONFIG}" "${OUT_DIR}" cellpose_prior
    done

    echo "  merging tiles for ${CONF_LABEL}..."
    /home/joe/miniforge3/bin/conda run -n segbench \
        python scripts/build_baysor_conf_adata.py "${CONF_FLOAT}"
    echo "  done: ${CONF_LABEL}"
done

echo "both confidence runs complete"
