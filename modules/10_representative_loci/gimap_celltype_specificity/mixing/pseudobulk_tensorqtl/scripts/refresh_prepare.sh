#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="${SC_PCQTL_GIMAP_MIXING_ROOT:?Set SC_PCQTL_GIMAP_MIXING_ROOT.}"
RSCRIPT="${COQTL_RSCRIPT:-${R_SCRIPT:-Rscript}}"
mkdir -p "${BASE}/logs"
"${RSCRIPT}" "${SCRIPT_DIR}/01_prepare_pseudobulk_tensorqtl_inputs.R" 2>&1 | tee "${BASE}/logs/01_prepare_pseudobulk_tensorqtl_inputs.log"
bash "${SCRIPT_DIR}/02_extract_gimap_genotype_region.sh" 2>&1 | tee "${BASE}/logs/02_extract_gimap_genotype_region.log"
