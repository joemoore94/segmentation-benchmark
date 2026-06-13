#!/usr/bin/env bash
# Run Baysor transcript-based segmentation.
#
# No standalone `baysor` executable is built (we installed the Baysor.jl
# package via Pkg.add, not PackageCompiler), so the CLI is invoked through
# `julia -e 'using Baysor; Baysor.command_main()' -- <args>`.
#
# Usage: run_baysor.sh <transcripts_file> <config_toml> <output_dir> [prior_segmentation_column]
#
# [prior_segmentation_column] is the name of a column in <transcripts_file>
# holding a prior cell assignment (e.g. from 10x or a DAPI-based segmentation)
# used to seed Baysor.

set -euo pipefail

TRANSCRIPTS=$1
CONFIG=$2
OUTPUT_DIR=$3
PRIOR_COLUMN=${4:-}

mkdir -p "${OUTPUT_DIR}"

if [[ -n "${PRIOR_COLUMN}" ]]; then
  julia -e 'using Baysor; Baysor.command_main()' -- run -c "${CONFIG}" -o "${OUTPUT_DIR}" "${TRANSCRIPTS}" ":${PRIOR_COLUMN}"
else
  julia -e 'using Baysor; Baysor.command_main()' -- run -c "${CONFIG}" -o "${OUTPUT_DIR}" "${TRANSCRIPTS}"
fi
