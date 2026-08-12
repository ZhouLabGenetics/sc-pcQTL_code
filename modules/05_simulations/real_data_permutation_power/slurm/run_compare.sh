#!/usr/bin/env bash
#SBATCH --job-name=power_compare
#SBATCH --output=logs/power_compare_%j.out
#SBATCH --error=logs/power_compare_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1

set -euo pipefail

CONFIG_LABEL_ARG="${1:-}"
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${CODE_DIR}/../shared/module_env.sh"
: "${COQTL_SIM_POWER_ROOT:?Set COQTL_SIM_POWER_ROOT}"
WORK_DIR="${COQTL_SIM_POWER_ROOT}"
source "${CODE_DIR}/config/paper_defaults.env"
CONFIG_LABEL="${CONFIG_LABEL_ARG:-${CONFIG_LABEL}}"

DATA_DIR="${COQTL_SIM_DATA_ROOT}/power_full/${CONFIG_LABEL}"
RESULTS_DIR="${WORK_DIR}/results/${CONFIG_LABEL}"
COMPARE_DIR="${RESULTS_DIR}/compare"
mkdir -p "${WORK_DIR}/logs" "${COMPARE_DIR}"

for strength_dir in "${DATA_DIR}/signal"/strength_*; do
  [[ -d "${strength_dir}" ]] || continue
  label=$(basename "${strength_dir}")
  out_dir="${COMPARE_DIR}/${label}"
  [[ ! -e "${out_dir}" ]] || {
    echo "Refusing to overwrite existing comparison output: ${out_dir}" >&2
    exit 1
  }
  "${COQTL_RSCRIPT}" "${CODE_DIR}/scripts/compare_power.R" \
    --truth_file "${DATA_DIR}/truth_pairs.tsv" \
    --sc_null_dir "${RESULTS_DIR}/sc/null" \
    --sc_signal_dir "${RESULTS_DIR}/sc/${label}" \
    --pb_null_dir "${RESULTS_DIR}/pb/null" \
    --pb_signal_dir "${RESULTS_DIR}/pb/${label}" \
    --p_thresholds "${P_THRESHOLDS}" \
    --output_dir "${out_dir}" \
    --label "${label}"
done
