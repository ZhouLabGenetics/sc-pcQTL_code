#!/bin/bash
# Submit SC hurdle association jobs for chromosomes 1-22

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${SC_PCQTL_HURDLE_WORK_DIR:?Set SC_PCQTL_HURDLE_WORK_DIR before submission.}"
LOG_DIR="${WORK_DIR}/logs"
RESULT_DIR="${WORK_DIR}/results/method2_sc_hurdle"
TMP_DIR="${RESULT_DIR}/tmp"
mkdir -p "${LOG_DIR}"
mkdir -p "${TMP_DIR}"
CELL_TYPE="${SC_PCQTL_CELL_TYPE:?Set SC_PCQTL_CELL_TYPE before submission.}"
RUN_R="${COQTL_RSCRIPT:-${R_SCRIPT:-Rscript}}"
python3 "${SCRIPT_DIR}/../../bin/celltype_eligibility.py" --require-primary "${CELL_TYPE}" >/dev/null

CHR_LIST="${1:-1-22}"

if [[ "${CHR_LIST}" == "1-22" ]]; then
  CHRS=$(seq 1 22)
else
  CHRS=${CHR_LIST}
fi

for chr in ${CHRS}; do
  sbatch --job-name=sc_${CELL_TYPE}_chr${chr} \
         --output="${LOG_DIR}/step2_chr${chr}_${CELL_TYPE}_%j.out" \
         --error="${LOG_DIR}/step2_chr${chr}_${CELL_TYPE}_%j.err" \
         --time=24:00:00 \
         --mem=48G \
         --cpus-per-task=1 \
         --export=ALL,TMPDIR="${TMP_DIR}" \
         --wrap="cd ${SCRIPT_DIR} && ${RUN_R} step2_calculate_sc_associations_single_chr.R --chr ${chr}"
done
