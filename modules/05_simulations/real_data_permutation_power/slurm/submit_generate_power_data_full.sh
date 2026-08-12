#!/bin/bash
#SBATCH --job-name=gen_power_full
#SBATCH --output=logs/generate_power_full_%j.out
#SBATCH --error=logs/generate_power_full_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=1

set -euo pipefail

CONFIG_LABEL_ARG="${1:-}"
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${CODE_DIR}/../shared/module_env.sh"
: "${COQTL_SIM_POWER_ROOT:?Set COQTL_SIM_POWER_ROOT}"
WORK_DIR="${COQTL_SIM_POWER_ROOT}"
source "${CODE_DIR}/config/paper_defaults.env"
CONFIG_LABEL="${CONFIG_LABEL_ARG:-${CONFIG_LABEL}}"
SCRIPT="${CODE_DIR}/scripts/generate_power_data_full.R"

mkdir -p "${WORK_DIR}/logs"
mkdir -p "${COQTL_SIM_DATA_ROOT}/power_full"

echo "================================================================================"
echo "Generate power_full datasets"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Host:   $(hostname)"
echo "Start:  $(date)"
echo "================================================================================"

"${COQTL_RSCRIPT}" "${SCRIPT}" \
  --raw_counts "${COQTL_RAW_COUNTS_FILE:?Set COQTL_RAW_COUNTS_FILE}" \
  --gene_info "${COQTL_GENE_INFO_FILE:?Set COQTL_GENE_INFO_FILE}" \
  --output "${COQTL_SIM_DATA_ROOT}/power_full" \
  --config_label "${CONFIG_LABEL}" \
  --n_genes "${N_GENES}" \
  --n_modules "${N_MODULES}" \
  --module_size "${MODULE_SIZE}" \
  --effect_strengths "${EFFECT_STRENGTHS}"

echo "================================================================================"
echo "End time: $(date)"
echo "================================================================================"
