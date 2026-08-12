#!/bin/bash
#SBATCH --job-name=pb_power_full_single
#SBATCH --output=logs/pb_single_%x_%j.out
#SBATCH --error=logs/pb_single_%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1

set -euo pipefail

CONFIG_LABEL_ARG="${1:-}"
TARGET_LABEL="${2:-null}"

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${CODE_DIR}/../shared/module_env.sh"
: "${COQTL_SIM_POWER_ROOT:?Set COQTL_SIM_POWER_ROOT}"
WORK_DIR="${COQTL_SIM_POWER_ROOT}"
source "${CODE_DIR}/config/paper_defaults.env"
CONFIG_LABEL="${CONFIG_LABEL_ARG:-${CONFIG_LABEL}}"
DATA_DIR="${COQTL_SIM_DATA_ROOT}/power_full/${CONFIG_LABEL}"

if [[ ! -d "${DATA_DIR}" ]]; then
  echo "Missing data dir: ${DATA_DIR}" >&2
  exit 1
fi

mkdir -p "${WORK_DIR}/logs"

if [[ "${TARGET_LABEL}" == "null" ]]; then
  PB_FILE="${DATA_DIR}/null/pb_counts_null.rds"
  OUT_DIR="${WORK_DIR}/results/${CONFIG_LABEL}/pb/null"
else
  PB_FILE="${DATA_DIR}/signal/${TARGET_LABEL}/pb_counts_signal.rds"
  OUT_DIR="${WORK_DIR}/results/${CONFIG_LABEL}/pb/${TARGET_LABEL}"
fi

if [[ ! -f "${PB_FILE}" ]]; then
  echo "Missing PB file: ${PB_FILE}" >&2
  exit 1
fi

echo "================================================================"
echo "PB Spearman single job"
echo "Config : ${CONFIG_LABEL}"
echo "Target : ${TARGET_LABEL}"
echo "Start  : $(date)"
echo "Node   : $(hostname)"
echo "================================================================"

"${COQTL_RSCRIPT}" "${CODE_DIR}/scripts/pb_run_spearman.R" \
  --pb_file "${PB_FILE}" \
  --output_dir "${OUT_DIR}"

echo "================================================================"
echo "End time: $(date)"
echo "================================================================"
