#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export COQTL_SIM_CODE_ROOT="${COQTL_SIM_CODE_ROOT:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
export SC_PCQTL_HURDLE_SCRIPT_DIR="${SC_PCQTL_HURDLE_SCRIPT_DIR:-${COQTL_SIM_CODE_ROOT}/../01_pairwise_hurdle_clustering}"
: "${COQTL_SIM_ROOT:?Set COQTL_SIM_ROOT to the external simulation-output root}"
export COQTL_SIM_DATA_ROOT="${COQTL_SIM_DATA_ROOT:-${COQTL_SIM_ROOT}/data}"
export COQTL_RAW_COUNTS_FILE="${COQTL_RAW_COUNTS_FILE:-}"
export COQTL_GENE_INFO_FILE="${COQTL_GENE_INFO_FILE:-}"

if [[ -z "${COQTL_RSCRIPT:-}" ]]; then
  export COQTL_RSCRIPT="${R_SCRIPT:-$(command -v Rscript || true)}"
fi
[[ -n "${COQTL_RSCRIPT}" && -x "${COQTL_RSCRIPT}" ]] || {
  echo "COQTL_RSCRIPT is not an executable Rscript path" >&2
  return 1 2>/dev/null || exit 1
}
