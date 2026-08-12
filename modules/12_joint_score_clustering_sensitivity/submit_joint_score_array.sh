#!/usr/bin/env bash

set -euo pipefail

MODULE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${SC_PCQTL_JOINT_SCORE_WORK_DIR:?Set SC_PCQTL_JOINT_SCORE_WORK_DIR}"
: "${SC_PCQTL_CELL_TYPE:?Set SC_PCQTL_CELL_TYPE}"

stage_dir="${SC_PCQTL_JOINT_SCORE_WORK_DIR}/results/method2_sc_hurdle/joint_score_stage"
manifest="${stage_dir}/joint_score_task_manifest.tsv"
[[ -s "${stage_dir}/STAGE_COMPLETE" ]] || {
  echo "Run 01_stage_inputs.R and 02_plan_chunks.R before submission." >&2
  exit 1
}
[[ -s "${manifest}" ]] || { echo "Missing task manifest: ${manifest}" >&2; exit 1; }

n_tasks="$(($(wc -l < "${manifest}") - 1))"
[[ "${n_tasks}" -gt 0 ]] || { echo "Task manifest is empty." >&2; exit 1; }

partition="${SC_PCQTL_SLURM_PARTITION:-normal}"
max_concurrent="${SC_PCQTL_SLURM_MAX_CONCURRENT:-20}"
score_time="${SC_PCQTL_SLURM_SCORE_TIME:-04:00:00}"
score_mem="${SC_PCQTL_SLURM_SCORE_MEM:-8G}"
merge_time="${SC_PCQTL_SLURM_MERGE_TIME:-02:00:00}"
merge_mem="${SC_PCQTL_SLURM_MERGE_MEM:-16G}"
log_dir="${SC_PCQTL_JOINT_SCORE_WORK_DIR}/logs/joint_score"
mkdir -p "${log_dir}"

score_job="$(sbatch --parsable \
  --partition="${partition}" \
  --time="${score_time}" \
  --mem="${score_mem}" \
  --array="1-${n_tasks}%${max_concurrent}" \
  --job-name="js_${SC_PCQTL_CELL_TYPE}" \
  --output="${log_dir}/score_%A_%a.out" \
  --error="${log_dir}/score_%A_%a.err" \
  "${MODULE_DIR}/slurm/run_joint_score_chunk.sh")"

merge_job="$(sbatch --parsable \
  --partition="${partition}" \
  --time="${merge_time}" \
  --mem="${merge_mem}" \
  --dependency="afterok:${score_job}" \
  --job-name="jsm_${SC_PCQTL_CELL_TYPE}" \
  --output="${log_dir}/merge_%j.out" \
  --error="${log_dir}/merge_%j.err" \
  "${MODULE_DIR}/slurm/run_merge_and_cluster.sh")"

printf 'stage\tjob_id\tdependency\tn_tasks\n' > "${stage_dir}/submission_chain.tsv"
printf 'joint_score\t%s\tnone\t%s\n' "${score_job}" "${n_tasks}" >> "${stage_dir}/submission_chain.tsv"
printf 'merge_cluster\t%s\tafterok:%s\t1\n' "${merge_job}" "${score_job}" >> "${stage_dir}/submission_chain.tsv"
printf 'score_job=%s\nmerge_cluster_job=%s\n' "${score_job}" "${merge_job}"
