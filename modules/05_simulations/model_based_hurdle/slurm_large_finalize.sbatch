#!/usr/bin/env bash
#SBATCH --job-name=scpcqtl_sim_final
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/slurm_finalize_%j.out
#SBATCH --error=logs/slurm_finalize_%j.err

set -euo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CODE_DIR}/../shared/module_env.sh"
: "${SC_PCQTL_MODEL_SIM_WORK_ROOT:?Set SC_PCQTL_MODEL_SIM_WORK_ROOT}"
WORK_DIR="${SC_PCQTL_MODEL_SIM_WORK_ROOT}"
mkdir -p "${WORK_DIR}/logs" "${WORK_DIR}/plots_target"

"${COQTL_RSCRIPT}" "${CODE_DIR}/10_make_target_figures.R" \
  --sim_root "${WORK_DIR}/simulations_large" \
  --analysis_root "${WORK_DIR}/analysis_large_compare" \
  --out_dir "${WORK_DIR}/plots_target"
