#!/usr/bin/env bash
set -euo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CODE_DIR}/../shared/module_env.sh"
source "${CODE_DIR}/config/paper_defaults.env"

RUN_ROOT="${COQTL_WITHIN_DONOR_ROOT:-${COQTL_SIM_ROOT}/03_sensitivity_within_donor_permutation}"
if [[ -e "${RUN_ROOT}/data" || -e "${RUN_ROOT}/results" ]]; then
  echo "Refusing to mix with an existing run beneath ${RUN_ROOT}" >&2
  exit 1
fi
mkdir -p "${RUN_ROOT}/logs"

job_id=$(sbatch --parsable \
  --output="${RUN_ROOT}/logs/pipeline_%j.out" \
  --error="${RUN_ROOT}/logs/pipeline_%j.err" \
  "${CODE_DIR}/slurm/run_pipeline.sbatch")
printf 'within_donor_permutation_job=%s\n' "${job_id}"
