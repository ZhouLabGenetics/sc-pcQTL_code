#!/bin/bash
# Submit Step 2: Cluster PCA — prints job ID to stdout
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
CELL_TYPE="${SC_PCQTL_CELL_TYPE:?Set SC_PCQTL_CELL_TYPE before submission.}"
python3 "${SCRIPT_DIR}/../../bin/celltype_eligibility.py" --require-primary "${CELL_TYPE}" >/dev/null
RUN_R="${COQTL_RSCRIPT:-${R_SCRIPT:-Rscript}}"
PCQTL_DIR="${SC_PCQTL_PCQTL_DIR:?Set SC_PCQTL_PCQTL_DIR before submission.}"
LOG_DIR="${PCQTL_DIR}/logs"
mkdir -p "${LOG_DIR}"

sbatch --parsable \
    --job-name=pcqtl_pca_${CELL_TYPE} \
    --output="${LOG_DIR}/step2_pca_${CELL_TYPE}_%j.out" \
    --error="${LOG_DIR}/step2_pca_${CELL_TYPE}_%j.err" \
    --time=04:00:00 \
    --mem=64G \
    --cpus-per-task=1 \
    --wrap="cd ${SCRIPT_DIR} && ${RUN_R} step2_cluster_pca.R"
