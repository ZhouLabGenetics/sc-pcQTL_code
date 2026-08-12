#!/usr/bin/env bash
#SBATCH --job-name=scpcqtl_sim_ref
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm_ref_%j.out
#SBATCH --error=logs/slurm_ref_%j.err

set -euo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CODE_DIR}/../shared/module_env.sh"
: "${SC_PCQTL_MODEL_SIM_WORK_ROOT:?Set SC_PCQTL_MODEL_SIM_WORK_ROOT}"
WORK_DIR="${SC_PCQTL_MODEL_SIM_WORK_ROOT}"

: "${COQTL_RAW_COUNTS_FILE:?Set COQTL_RAW_COUNTS_FILE to the reference count matrix}"
mkdir -p "${WORK_DIR}/logs" "${WORK_DIR}/reference_large"

"${COQTL_RSCRIPT}" "${CODE_DIR}/01_fit_reference_from_real_data.R" \
  --count_file "${COQTL_RAW_COUNTS_FILE}" \
  --n_reference_cells "${N_REFERENCE_CELLS:-20000}" \
  --n_reference_genes "${N_REFERENCE_GENES:-3000}" \
  --out_dir "${WORK_DIR}/reference_large"
