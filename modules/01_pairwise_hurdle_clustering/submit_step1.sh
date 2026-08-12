#!/bin/bash
# Run sparsity filtering step (single SLURM job)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${SC_PCQTL_HURDLE_WORK_DIR:?Set SC_PCQTL_HURDLE_WORK_DIR before submission.}"
LOG_DIR="${WORK_DIR}/logs"
mkdir -p "${LOG_DIR}"
CELL_TYPE="${SC_PCQTL_CELL_TYPE:?Set SC_PCQTL_CELL_TYPE before submission.}"
RUN_R="${COQTL_RSCRIPT:-${R_SCRIPT:-Rscript}}"
python3 "${SCRIPT_DIR}/../../bin/celltype_eligibility.py" --require-primary "${CELL_TYPE}" >/dev/null

sbatch --job-name=scfilt_${CELL_TYPE} \
       --output="${LOG_DIR}/step1_${CELL_TYPE}_%j.out" \
       --error="${LOG_DIR}/step1_${CELL_TYPE}_%j.err" \
       --time=04:00:00 \
       --mem=64G \
       --wrap="cd ${SCRIPT_DIR} && ${RUN_R} step1_filter_sparse_genes.R"
