#!/usr/bin/env bash
#SBATCH --job-name=null_compare
#SBATCH --output=logs/null_compare_%j.out
#SBATCH --error=logs/null_compare_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1

set -euo pipefail

CONFIG_LABEL_ARG="${1:-}"
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${CODE_DIR}/../shared/module_env.sh"
: "${COQTL_SIM_NULL_ROOT:?Set COQTL_SIM_NULL_ROOT}"
WORK_DIR="${COQTL_SIM_NULL_ROOT}"
source "${CODE_DIR}/config/paper_defaults.env"
CONFIG_LABEL="${CONFIG_LABEL_ARG:-${CONFIG_LABEL}}"

RESULTS_DIR="${WORK_DIR}/results/${CONFIG_LABEL}"
mkdir -p "${WORK_DIR}/logs" "${RESULTS_DIR}"

"${COQTL_RSCRIPT}" "${CODE_DIR}/scripts/compare_fpr_results.R" \
  --results_dir "${RESULTS_DIR}" \
  --sc_dir "${RESULTS_DIR}/sc_hurdle" \
  --pb_dir "${RESULTS_DIR}/pb_spearman" \
  --setting_label "${CONFIG_LABEL}" \
  --p_thresholds "${P_THRESHOLDS}" \
  --output_dir "${RESULTS_DIR}"
