#!/usr/bin/env bash
#SBATCH --cpus-per-task=1

set -euo pipefail

MODULE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RSCRIPT="${SC_PCQTL_RSCRIPT:-Rscript}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

"${RSCRIPT}" --vanilla "${MODULE_DIR}/03_run_joint_score_chunk.R" "${SLURM_ARRAY_TASK_ID}"
