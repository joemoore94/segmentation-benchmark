#!/usr/bin/env bash
# Run Baysor PSC=1.0 with StarDist and Mesmer nuclear priors.
#
# Compares density-adaptive expansion (Baysor) with different nuclear
# detectors as seeds, all at PSC=1.0 (hard-locked nuclear transcripts).
# The CellPose prior already exists from run_baysor_c10.sh.
#
# Usage: bash scripts/run_baysor_nuclear_priors.sh
#        (run from repo root)

set -euo pipefail

PYTHON="/home/joe/miniforge3/envs/segbench/bin/python"
CONFIG="configs/baysor_prior_config_c10.toml"
ROI_DIR="data/processed/roi"

for METHOD in stardist mesmer; do
    PRIOR_COL="${METHOD}_prior"
    TILES_DIR="${ROI_DIR}/baysor_${METHOD}_prior_c10_tiles"

    echo ""
    echo "========================================"
    echo "  Baysor PSC=1.0 with ${METHOD} prior"
    echo "========================================"

    # Step 1: Add nuclear prior column to transcripts
    echo "  Adding ${METHOD} prior column..."
    ${PYTHON} scripts/add_nuclear_prior.py "${METHOD}"

    # Step 2: Tile transcripts
    echo "  Tiling transcripts..."
    ${PYTHON} scripts/tile_nuclear_prior_transcripts.py "${METHOD}"

    # Step 3: Run Baysor on each tile
    for TILE in x0_y0 x1_y0 x0_y1 x1_y1; do
        TRANSCRIPT_FILE="${TILES_DIR}/transcripts_${TILE}.csv"
        OUT_DIR="${TILES_DIR}/${TILE}"
        echo "  tile ${TILE} -> ${OUT_DIR}"
        bash scripts/run_baysor.sh "${TRANSCRIPT_FILE}" "${CONFIG}" "${OUT_DIR}" "${PRIOR_COL}"
    done

    # Step 4: Merge tiles into AnnData
    echo "  Merging tiles..."
    ${PYTHON} -c "
from scripts.build_baysor_conf_adata import conf_label
import json, pandas as pd
from pathlib import Path
from segbench.quantify import quantify_baysor

ROI_DIR = Path('${ROI_DIR}')
TILES_DIR = ROI_DIR / 'baysor_${METHOD}_prior_c10_tiles'
GENE_RE = r\"^b'(.*)'$\"

with open(TILES_DIR / 'manifest.json') as f:
    manifest = json.load(f)

tables = []
for name, info in manifest.items():
    seg = pd.read_csv(TILES_DIR / name / 'segmentation.csv')
    seg['gene'] = seg['gene'].str.replace(GENE_RE, r'\1', regex=True)
    assigned = seg.dropna(subset=['cell'])
    centroids = assigned.groupby('cell')[['x', 'y']].mean()
    core_x, core_y = info['core_x'], info['core_y']
    core_cells = centroids.index[
        centroids['x'].between(*core_x) & centroids['y'].between(*core_y)
    ]
    tile_table = assigned[assigned['cell'].isin(core_cells)].copy()
    tile_table['cell'] = name + '_' + tile_table['cell'].astype(str)
    tables.append(tile_table)
    print(f'{name}: {len(core_cells)} core cells, {len(tile_table)} assigned molecules')

merged = pd.concat(tables, ignore_index=True)
before = len(merged)
merged = merged.drop_duplicates(subset='transcript_id', keep='first')
print(f'dropped {before - len(merged)} duplicate transcripts from overlap bands')

adata = quantify_baysor(merged)
out_path = ROI_DIR / 'adata_baysor_${METHOD}_prior_c10.h5ad'
adata.write_h5ad(out_path)
print(f'wrote {out_path}: {adata.n_obs} cells x {adata.n_vars} genes, {int(adata.X.sum())} transcripts')
"
    echo "  done: ${METHOD}"
done

echo ""
echo "All nuclear prior runs complete."
