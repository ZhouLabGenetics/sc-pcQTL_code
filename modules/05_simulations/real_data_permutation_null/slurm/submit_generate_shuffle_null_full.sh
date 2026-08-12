#!/bin/bash
#SBATCH --job-name=gen_shuffle_full
#SBATCH --output=logs/generate_shuffle_full_%j.out
#SBATCH --error=logs/generate_shuffle_full_%j.err
#SBATCH --time=08:00:00
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
SCRIPT="${CODE_DIR}/scripts/generate_shuffle_null_full.R"

mkdir -p "${WORK_DIR}/logs"
mkdir -p "${COQTL_SIM_DATA_ROOT}/shuffle_full"

"${COQTL_RSCRIPT}" "${SCRIPT}" \
  --raw_counts "${COQTL_RAW_COUNTS_FILE:?Set COQTL_RAW_COUNTS_FILE to the OneK1K count matrix}" \
  --gene_info "${COQTL_GENE_INFO_FILE:?Set COQTL_GENE_INFO_FILE}" \
  --output "${COQTL_SIM_DATA_ROOT}/shuffle_full" \
  --label "${CONFIG_LABEL}" \
  --n_genes "${N_GENES}"
