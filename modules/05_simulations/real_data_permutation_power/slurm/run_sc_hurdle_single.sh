#!/bin/bash
#SBATCH --job-name=sc_power_full_single
#SBATCH --output=logs/sc_single_%x_%j.out
#SBATCH --error=logs/sc_single_%x_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=32G
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
  SC_FILE="${DATA_DIR}/null/sc_counts_null.rds"
  OUT_DIR="${WORK_DIR}/results/${CONFIG_LABEL}/sc/null"
else
  SC_FILE="${DATA_DIR}/signal/${TARGET_LABEL}/sc_counts_signal.rds"
  OUT_DIR="${WORK_DIR}/results/${CONFIG_LABEL}/sc/${TARGET_LABEL}"
fi

if [[ ! -f "${SC_FILE}" ]]; then
  echo "Missing SC file: ${SC_FILE}" >&2
  exit 1
fi

echo "================================================================"
echo "SC fasthurdle single job"
echo "Config : ${CONFIG_LABEL}"
echo "Target : ${TARGET_LABEL}"
echo "Start  : $(date)"
echo "Node   : $(hostname)"
echo "================================================================"

"${COQTL_RSCRIPT}" "${CODE_DIR}/scripts/sc_run_hurdle.R" \
  --sc_file "${SC_FILE}" \
  --output_dir "${OUT_DIR}" \
  --nonzero_cutoff "${NONZERO_CUTOFF}" \
  --p_thresholds "${P_THRESHOLDS}"

echo "================================================================"
echo "End time: $(date)"
echo "================================================================"
