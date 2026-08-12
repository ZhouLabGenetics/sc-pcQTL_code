#!/usr/bin/env bash
set -euo pipefail

CONFIG_LABEL_ARG="${1:-}"
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${CODE_DIR}/../shared/module_env.sh"
: "${COQTL_SIM_NULL_ROOT:?Set COQTL_SIM_NULL_ROOT}"
WORK_DIR="${COQTL_SIM_NULL_ROOT}"
source "${CODE_DIR}/config/paper_defaults.env"
CONFIG_LABEL="${CONFIG_LABEL_ARG:-${CONFIG_LABEL}}"
DATA_DIR="${COQTL_SIM_DATA_ROOT}/shuffle_full/${CONFIG_LABEL}"
RESULTS_DIR="${WORK_DIR}/results/${CONFIG_LABEL}"
if [[ -e "${DATA_DIR}" || -e "${RESULTS_DIR}" ]]; then
  echo "Refusing to mix with an existing run: ${DATA_DIR} or ${RESULTS_DIR}" >&2
  echo "Use a new config label or archive/remove the old run explicitly." >&2
  exit 1
fi
mkdir -p "${WORK_DIR}/logs"
cd "${WORK_DIR}"

gen_job=$(sbatch --parsable "${CODE_DIR}/slurm/submit_generate_shuffle_null_full.sh" "${CONFIG_LABEL}")
sc_job=$(sbatch --parsable --dependency="afterok:${gen_job}" \
  "${CODE_DIR}/slurm/submit_sc_hurdle_test.sh" "${CONFIG_LABEL}")
pb_job=$(sbatch --parsable --dependency="afterok:${gen_job}" \
  "${CODE_DIR}/slurm/submit_pb_spearman_test.sh" "${CONFIG_LABEL}")
compare_job=$(sbatch --parsable --dependency="afterok:${sc_job}:${pb_job}" \
  "${CODE_DIR}/slurm/run_compare.sh" "${CONFIG_LABEL}")

printf 'generation_job=%s\nsc_job=%s\npb_job=%s\ncompare_job=%s\n' \
  "${gen_job}" "${sc_job}" "${pb_job}" "${compare_job}"
