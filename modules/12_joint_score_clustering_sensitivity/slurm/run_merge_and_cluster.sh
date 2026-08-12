#!/usr/bin/env bash
#SBATCH --cpus-per-task=1

set -euo pipefail

MODULE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RSCRIPT="${SC_PCQTL_RSCRIPT:-Rscript}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

"${RSCRIPT}" --vanilla "${MODULE_DIR}/04_merge_joint_score_chunks.R"
for chromosome in $(seq 1 22); do
  "${RSCRIPT}" --vanilla "${MODULE_DIR}/05_call_clusters_chr.R" "${chromosome}"
done
"${RSCRIPT}" --vanilla "${MODULE_DIR}/06_merge_clusters.R"
