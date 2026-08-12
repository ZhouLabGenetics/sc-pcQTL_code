#!/bin/bash
# Submit Step 3 prep (depends on Step 2 job) — prints job ID to stdout
# Usage: bash submit_step3_prepare.sh <step2_job_id>
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
CELL_TYPE="${SC_PCQTL_CELL_TYPE:?Set SC_PCQTL_CELL_TYPE before submission.}"
python3 "${SCRIPT_DIR}/../../bin/celltype_eligibility.py" --require-primary "${CELL_TYPE}" >/dev/null
RUN_R="${COQTL_RSCRIPT:-${R_SCRIPT:-Rscript}}"
PCQTL_DIR="${SC_PCQTL_PCQTL_DIR:?Set SC_PCQTL_PCQTL_DIR before submission.}"
LOG_DIR="${PCQTL_DIR}/logs"
STEP2_JOB=${1:?Usage: $0 <step2_job_id>}

sbatch --parsable \
    --dependency=afterok:${STEP2_JOB} \
    --job-name=pcqtl_prep_${CELL_TYPE} \
    --output="${LOG_DIR}/step3_prep_${CELL_TYPE}_%j.out" \
    --error="${LOG_DIR}/step3_prep_${CELL_TYPE}_%j.err" \
    --time=01:00:00 \
    --mem=16G \
    --cpus-per-task=1 \
    --wrap="cd ${SCRIPT_DIR} && ${RUN_R} step3_prepare_saige_inputs.R"
