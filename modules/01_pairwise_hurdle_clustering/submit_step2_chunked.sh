#!/bin/bash
# Submit SC hurdle association jobs using chunked approach (optimized for <6h per job)

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

echo "=== Method 2 Step 2 (${CELL_TYPE}): Submit Chunked SC Hurdle Jobs ==="
echo ""

# Check if Step 1 is complete
FILTER_FILE="${RESULT_DIR}/filtered_genes.tsv"
if [[ ! -f "${FILTER_FILE}" ]]; then
  echo "❌ Error: filtered_genes.tsv not found"
  echo "Please run Step 1 first: bash submit_step1.sh"
  exit 1
fi

# Calculate chunk plan
echo "Step 1: Calculating optimal chunks..."
cd "${SCRIPT_DIR}"
"${RUN_R}" calculate_chunks.R

CHUNK_PLAN="${WORK_DIR}/chunk_plan.tsv"
if [[ ! -f "${CHUNK_PLAN}" ]]; then
  echo "❌ Error: chunk_plan.tsv not created"
  exit 1
fi

# Count total chunks
TOTAL_CHUNKS=$(tail -n +2 "${CHUNK_PLAN}" | wc -l)
echo "Total chunks to submit: ${TOTAL_CHUNKS}"
echo ""

# Submit jobs
echo "Step 2: Submitting jobs..."
SUBMITTED=0

while IFS=$'\t' read -r chr chunk_id total_chunks gene_start gene_end n_genes; do
  # Skip header
  if [[ "$chr" == "chromosome" ]]; then continue; fi

  sbatch --job-name=sc_${CELL_TYPE}_ch${chr}_ck${chunk_id} \
         --output="${LOG_DIR}/step2_chr${chr}_chunk${chunk_id}_${CELL_TYPE}_%j.out" \
         --error="${LOG_DIR}/step2_chr${chr}_chunk${chunk_id}_${CELL_TYPE}_%j.err" \
         --time=06:00:00 \
         --mem=32G \
         --cpus-per-task=1 \
         --export=ALL,TMPDIR="${TMP_DIR}" \
         --wrap="cd ${SCRIPT_DIR} && ${RUN_R} step2_calculate_sc_associations_chunked.R \
                --chr ${chr} \
                --gene_start ${gene_start} \
                --gene_end ${gene_end} \
                --chunk_id ${chunk_id}"

  SUBMITTED=$((SUBMITTED + 1))
  if [[ $((SUBMITTED % 10)) -eq 0 ]]; then
    echo "  Submitted ${SUBMITTED}/${TOTAL_CHUNKS} jobs..."
  fi
done < "${CHUNK_PLAN}"

echo ""
echo "=== All ${SUBMITTED} jobs submitted ==="
echo ""
echo "Check status: squeue -u \$USER | grep sc_${CELL_TYPE}"
echo "Check logs: ls -lt ${LOG_DIR}/step2_*_${CELL_TYPE}_*.out | head"
echo ""
echo "After all jobs complete, run: bash submit_step2b_merge.sh"
