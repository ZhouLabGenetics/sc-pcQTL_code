#!/bin/bash
#SBATCH --job-name=sc_shuffle_full
#SBATCH --output=logs/sc_hurdle_%j.out
#SBATCH --error=logs/sc_hurdle_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

set -euo pipefail

CONFIG_LABEL_ARG="${1:-}"
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${CODE_DIR}/../shared/module_env.sh"
: "${COQTL_SIM_NULL_ROOT:?Set COQTL_SIM_NULL_ROOT}"
WORK_DIR="${COQTL_SIM_NULL_ROOT}"
source "${CODE_DIR}/config/paper_defaults.env"
CONFIG_LABEL="${CONFIG_LABEL_ARG:-${CONFIG_LABEL}}"
SC_SCRIPT="${CODE_DIR}/scripts/test_sc_hurdle_null.R"
SC_FILE="${COQTL_SIM_DATA_ROOT}/shuffle_full/${CONFIG_LABEL}/sc_counts_shuffle_null.rds"
OUT_DIR="${WORK_DIR}/results/${CONFIG_LABEL}/sc_hurdle"

mkdir -p "${WORK_DIR}/logs" "${OUT_DIR}"

"${COQTL_RSCRIPT}" "${SC_SCRIPT}" \
  --sc_null_file "${SC_FILE}" \
  --output_dir "${OUT_DIR}" \
  --nonzero_cutoff "${NONZERO_CUTOFF}" \
  --p_thresholds "${P_THRESHOLDS}"
