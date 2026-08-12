#!/bin/bash
# Submit Step 3 SAIGE array job (depends on prep job) — prints job ID to stdout
# Usage: bash submit_step3_saige.sh <prep_job_id>
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
CELL_TYPE="${SC_PCQTL_CELL_TYPE:?Set SC_PCQTL_CELL_TYPE before submission.}"
python3 "${SCRIPT_DIR}/../../bin/celltype_eligibility.py" --require-primary "${CELL_TYPE}" >/dev/null
PCQTL_DIR="${SC_PCQTL_PCQTL_DIR:?Set SC_PCQTL_PCQTL_DIR before submission.}"
LOG_DIR="${PCQTL_DIR}/logs"
PREP_JOB=${1:?Usage: $0 <prep_job_id>}

# array size = number of clusters in cluster_summary (some may be skipped at runtime)
UPSTREAM_CELLTYPES_ROOT="${SC_PCQTL_UPSTREAM_CELLTYPES_ROOT:-${SC_PCQTL_UPSTREAM_ROOT:?Set SC_PCQTL_UPSTREAM_ROOT.}/celltypes}"
CLUSTER_FILE="${SC_PCQTL_CLUSTER_FILE:-${UPSTREAM_CELLTYPES_ROOT}/${CELL_TYPE}/cluster_identification/results/method2_sc_hurdle/merged_clusters/cluster_summary.tsv}"
N_CLUSTERS=$(( $(wc -l < "${CLUSTER_FILE}") - 1 ))

sbatch --parsable \
    --dependency=afterok:${PREP_JOB} \
    --job-name=saige_${CELL_TYPE} \
    --output="${LOG_DIR}/step3_saige_${CELL_TYPE}_%a.out" \
    --error="${LOG_DIR}/step3_saige_${CELL_TYPE}_%a.err" \
    --time=12:00:00 \
    --mem=32G \
    --cpus-per-task=1 \
    --array=1-${N_CLUSTERS} \
    --wrap="cd ${SCRIPT_DIR} && bash step3_run_saige_cluster.sh \${SLURM_ARRAY_TASK_ID}"
