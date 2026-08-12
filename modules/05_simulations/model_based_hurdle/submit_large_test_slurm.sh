#!/usr/bin/env bash
set -euo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CODE_DIR}/../shared/module_env.sh"
: "${SC_PCQTL_MODEL_SIM_WORK_ROOT:?Set SC_PCQTL_MODEL_SIM_WORK_ROOT}"

WORK_DIR="${SC_PCQTL_MODEL_SIM_WORK_ROOT}"
MANIFEST="${WORK_DIR}/manifests/large_test_manifest.tsv"
ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-20}"
N_REPS="${N_REPS:-50}"
SCENARIOS=(
  null zero_only_low count_only_low both_low
  zero_only_high count_only_high both_high
)

mkdir -p "${WORK_DIR}/logs" "${WORK_DIR}/manifests"
{
  printf "scenario\treplicate\n"
  for scenario in "${SCENARIOS[@]}"; do
    for rep in $(seq 1 "${N_REPS}"); do
      printf "%s\t%s\n" "${scenario}" "${rep}"
    done
  done
} > "${MANIFEST}"

task_count=$(( ${#SCENARIOS[@]} * N_REPS ))
cd "${WORK_DIR}"
ref_job=$(sbatch --parsable "${CODE_DIR}/slurm_large_reference.sbatch")
array_job=$(sbatch --parsable --dependency="afterok:${ref_job}" \
  --array="1-${task_count}%${ARRAY_CONCURRENCY}" \
  "${CODE_DIR}/slurm_large_task.sbatch" "${MANIFEST}")
fin_job=$(sbatch --parsable --dependency="afterok:${array_job}" \
  "${CODE_DIR}/slurm_large_finalize.sbatch")

printf 'reference_job=%s\narray_job=%s\nfinalize_job=%s\nmanifest=%s\n' \
  "${ref_job}" "${array_job}" "${fin_job}" "${MANIFEST}"
