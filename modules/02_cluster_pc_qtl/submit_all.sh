#!/bin/bash
# Full pcQTL pipeline: Step2 PCA → Step3 Prepare → Step3 SAIGE-QTL
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
CELL_TYPE="${SC_PCQTL_CELL_TYPE:?Set SC_PCQTL_CELL_TYPE before submission.}"
python3 "${SCRIPT_DIR}/../../bin/celltype_eligibility.py" --require-primary "${CELL_TYPE}" >/dev/null

echo "=== pcQTL: ${CELL_TYPE} ==="

STEP2_JOB=$(bash submit_step2.sh)
echo "  Step2 PCA       : job ${STEP2_JOB}"

PREP_JOB=$(bash submit_step3_prepare.sh ${STEP2_JOB})
echo "  Step3 Prepare   : job ${PREP_JOB}"

SAIGE_JOB=$(bash submit_step3_saige.sh ${PREP_JOB})
echo "  Step3 SAIGE-QTL : job ${SAIGE_JOB}"
