#!/bin/bash
# Submit the entire SC hurdle pipeline (chunked) using SLURM dependencies only

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

WORK_DIR="${SC_PCQTL_HURDLE_WORK_DIR:?Set SC_PCQTL_HURDLE_WORK_DIR before submission.}"
LOG_DIR="${WORK_DIR}/logs"
RESULT_DIR="${WORK_DIR}/results/method2_sc_hurdle"
TMP_DIR="${RESULT_DIR}/tmp"
RUN_R="${COQTL_RSCRIPT:-${R_SCRIPT:-Rscript}}"
CELL_TYPE="${SC_PCQTL_CELL_TYPE:?Set SC_PCQTL_CELL_TYPE before submission.}"
ELIGIBILITY_MANIFEST="${SC_PCQTL_CELLTYPE_MANIFEST:-${SCRIPT_DIR}/../../config/celltype_eligibility.tsv}"

eligibility=$(awk -F '\t' -v ct="${CELL_TYPE}" 'NR > 1 && $1 == ct {print $4 "\t" $5 "\t" $6}' "${ELIGIBILITY_MANIFEST}")
if [[ -z "${eligibility}" ]]; then
  echo "Cell type ${CELL_TYPE} is absent from ${ELIGIBILITY_MANIFEST}" >&2
  exit 1
fi
IFS=$'\t' read -r N_CELLS INCLUDE_PRIMARY EXCLUSION_REASON <<< "${eligibility}"
if [[ "${INCLUDE_PRIMARY}" != "TRUE" ]]; then
  echo "Skipping ${CELL_TYPE}: ${N_CELLS} cells; ${EXCLUSION_REASON}."
  exit 0
fi

CHUNK_PLAN="${WORK_DIR}/chunk_plan.tsv"

mkdir -p "${LOG_DIR}"
mkdir -p "${TMP_DIR}"

echo "=== Method 2 SC Hurdle (${CELL_TYPE}): Full Chunked Pipeline ==="
echo ""

# ---------------------------------------------------------------------------
# Step 1: sparsity filtering
echo "[Step 1] Filtering sparse genes..."
STEP1_JOB=$(sbatch --parsable \
                   --job-name=scfilt_${CELL_TYPE} \
                   --output=${LOG_DIR}/step1_${CELL_TYPE}_%j.out \
                   --error=${LOG_DIR}/step1_${CELL_TYPE}_%j.err \
                   --time=04:00:00 \
                   --mem=64G \
                   --cpus-per-task=1 \
                   --export=ALL,TMPDIR="${TMP_DIR}" \
                   --wrap="mkdir -p ${TMP_DIR} && export TMPDIR=${TMP_DIR} && ${RUN_R} step1_filter_sparse_genes.R")
echo "  Submitted job: ${STEP1_JOB}"

# ---------------------------------------------------------------------------
# Stage 2+: orchestrate remaining steps from a dependency-controlled job
STAGE2_JOB=$(sbatch --parsable \
                    --dependency=afterok:${STEP1_JOB} \
                    --job-name=scstage2_${CELL_TYPE} \
                    --output=${LOG_DIR}/stage2_${CELL_TYPE}_%j.out \
                    --error=${LOG_DIR}/stage2_${CELL_TYPE}_%j.err \
                    --time=03:00:00 \
                    --mem=16G \
                    --cpus-per-task=1 \
                    --export=ALL,SCRIPT_DIR="${SCRIPT_DIR}",LOG_DIR="${LOG_DIR}",RUN_R="${RUN_R}",CELL_TYPE="${CELL_TYPE}",CHUNK_PLAN="${CHUNK_PLAN}",TMP_DIR="${TMP_DIR}",TMPDIR="${TMP_DIR}" <<'EOF'
#!/bin/bash
set -euo pipefail

cd "${SCRIPT_DIR}"
mkdir -p "${TMP_DIR}"
export TMPDIR="${TMP_DIR}"

echo "[Stage 2] Calculating chunk plan..."
${RUN_R} calculate_chunks.R
if [[ ! -f "${CHUNK_PLAN}" ]]; then
  echo "❌ Error: chunk plan ${CHUNK_PLAN} not found." >&2
  exit 1
fi

echo "[Stage 2] Precomputing library size cache..."
CACHE_JOB=$(sbatch --parsable \
                   --dependency=afterok:${SLURM_JOB_ID} \
                   --job-name=sccache_${CELL_TYPE} \
                   --output=${LOG_DIR}/cache_${CELL_TYPE}_%j.out \
                   --error=${LOG_DIR}/cache_${CELL_TYPE}_%j.err \
                   --time=08:00:00 \
                   --mem=64G \
                   --cpus-per-task=1 \
                   --export=ALL,TMPDIR="${TMP_DIR}" \
                   --wrap="mkdir -p ${TMP_DIR} && export TMPDIR=${TMP_DIR} && ${RUN_R} compute_library_size_cache.R")
echo "    Cache job: ${CACHE_JOB}"

echo "[Stage 2] Submitting chunk jobs..."
chunk_ids=()
while IFS=$'\t' read -r chromosome chunk_id total_chunks gene_start gene_end n_genes_chunk; do
  if [[ "${chromosome}" == "chromosome" ]]; then
    continue
  fi
  jobid=$(sbatch --parsable \
                 --dependency=afterok:${CACHE_JOB} \
                 --job-name=sc_${CELL_TYPE}_ch${chromosome}_ck${chunk_id} \
                 --output=${LOG_DIR}/step2_chr${chromosome}_chunk${chunk_id}_${CELL_TYPE}_%j.out \
                 --error=${LOG_DIR}/step2_chr${chromosome}_chunk${chunk_id}_${CELL_TYPE}_%j.err \
                 --time=06:00:00 \
                 --mem=64G \
                 --cpus-per-task=1 \
                 --export=ALL,TMPDIR="${TMP_DIR}" \
                 --wrap="${RUN_R} step2_calculate_sc_associations_chunked.R \
                         --chr ${chromosome} \
                         --gene_start ${gene_start} \
                         --gene_end ${gene_end} \
                         --chunk_id ${chunk_id}")
  chunk_ids+=("${jobid}")
  printf "    Chr%-2s chunk %03d -> %s\n" "${chromosome}" "${chunk_id}" "${jobid}"
done < "${CHUNK_PLAN}"

if [[ ${#chunk_ids[@]} -eq 0 ]]; then
  echo "❌ No chunk jobs submitted." >&2
  exit 1
fi

chunk_dep=$(IFS=:; echo "${chunk_ids[*]}")

echo "[Stage 2c] Merging chunk outputs..."
MERGE_JOB=$(sbatch --parsable \
                   --dependency=afterok:${chunk_dep} \
                   --job-name=scmerge_${CELL_TYPE} \
                   --output=${LOG_DIR}/step2c_${CELL_TYPE}_%j.out \
                   --error=${LOG_DIR}/step2c_${CELL_TYPE}_%j.err \
                   --time=02:00:00 \
                   --mem=16G \
                   --cpus-per-task=1 \
                   --export=ALL,TMPDIR="${TMP_DIR}" \
                   --wrap="mkdir -p ${TMP_DIR} && export TMPDIR=${TMP_DIR} && ${RUN_R} step2b_merge_chunks.R --chr 0")
echo "    Merge job: ${MERGE_JOB}"

echo "[Stage 3] Submitting sliding-window cluster jobs..."
step3_ids=()
for chr in $(seq 1 22); do
  jobid=$(sbatch --parsable \
                 --dependency=afterok:${MERGE_JOB} \
                 --job-name=sccl_${CELL_TYPE}_chr${chr} \
                 --output=${LOG_DIR}/step3_chr${chr}_${CELL_TYPE}_%j.out \
                 --error=${LOG_DIR}/step3_chr${chr}_${CELL_TYPE}_%j.err \
                 --time=04:00:00 \
                 --mem=16G \
                 --cpus-per-task=1 \
                 --export=ALL,TMPDIR="${TMP_DIR}" \
                 --wrap="mkdir -p ${TMP_DIR} && export TMPDIR=${TMP_DIR} && ${RUN_R} step3_identify_clusters_single_chr.R ${chr}")
  step3_ids+=("${jobid}")
  printf "    Chr%-2d -> %s\n" "${chr}" "${jobid}"
done

if [[ ${#step3_ids[@]} -eq 0 ]]; then
  echo "❌ No Step 3 jobs submitted." >&2
  exit 1
fi

step3_dep=$(IFS=:; echo "${step3_ids[*]}")

echo "[Stage 4] Submitting final merge..."
FINAL_JOB=$(sbatch --parsable \
                   --dependency=afterok:${step3_dep} \
                   --job-name=scfinal_${CELL_TYPE} \
                   --output=${LOG_DIR}/step4_${CELL_TYPE}_%j.out \
                   --error=${LOG_DIR}/step4_${CELL_TYPE}_%j.err \
                   --time=01:00:00 \
                   --mem=8G \
                   --cpus-per-task=1 \
                   --export=ALL,TMPDIR="${TMP_DIR}" \
                   --wrap="mkdir -p ${TMP_DIR} && export TMPDIR=${TMP_DIR} && ${RUN_R} step4_merge_clusters.R")

echo "  Final merge job: ${FINAL_JOB}"
EOF
)

echo "[Info] Stage 2 orchestration job: ${STAGE2_JOB}"
echo ""
echo "All jobs submitted with SLURM dependencies. Monitor via:"
echo "  squeue -u $USER | grep sc_${CELL_TYPE}"
