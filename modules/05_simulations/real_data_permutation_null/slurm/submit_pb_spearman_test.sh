#!/bin/bash
#SBATCH --job-name=pb_shuffle_full
#SBATCH --output=logs/pb_spearman_%j.out
#SBATCH --error=logs/pb_spearman_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -euo pipefail

CONFIG_LABEL_ARG="${1:-}"
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${CODE_DIR}/../shared/module_env.sh"
: "${COQTL_SIM_NULL_ROOT:?Set COQTL_SIM_NULL_ROOT}"
WORK_DIR="${COQTL_SIM_NULL_ROOT}"
source "${CODE_DIR}/config/paper_defaults.env"
CONFIG_LABEL="${CONFIG_LABEL_ARG:-${CONFIG_LABEL}}"
PB_SCRIPT="${CODE_DIR}/scripts/test_pb_spearman_null.R"
PB_FILE="${COQTL_SIM_DATA_ROOT}/shuffle_full/${CONFIG_LABEL}/pb_counts_shuffle_null.rds"
OUT_DIR="${WORK_DIR}/results/${CONFIG_LABEL}/pb_spearman"

mkdir -p "${WORK_DIR}/logs" "${OUT_DIR}"

"${COQTL_RSCRIPT}" "${PB_SCRIPT}" \
  --pb_null_file "${PB_FILE}" \
  --output_dir "${OUT_DIR}" \
  --p_thresholds "${P_THRESHOLDS}"
