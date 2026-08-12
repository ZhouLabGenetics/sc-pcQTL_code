#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="${SC_PCQTL_GIMAP_SPECIFICITY_ROOT:?Set SC_PCQTL_GIMAP_SPECIFICITY_ROOT.}"
RSCRIPT="${COQTL_RSCRIPT:-${R_SCRIPT:-Rscript}}"

mkdir -p "${ANALYSIS_DIR}/logs"

"${RSCRIPT}" "${SCRIPT_DIR}/run_gimap_specificity_analysis.R" \
  2>&1 | tee "${ANALYSIS_DIR}/logs/run_gimap_specificity_analysis.log"
