#!/usr/bin/env bash
set -euo pipefail

CONFIG_LABEL_ARG="${1:-}"
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${CODE_DIR}/../shared/module_env.sh"
: "${COQTL_SIM_POWER_ROOT:?Set COQTL_SIM_POWER_ROOT}"
WORK_DIR="${COQTL_SIM_POWER_ROOT}"
source "${CODE_DIR}/config/paper_defaults.env"
CONFIG_LABEL="${CONFIG_LABEL_ARG:-${CONFIG_LABEL}}"

DATA_DIR="${COQTL_SIM_DATA_ROOT}/power_full/${CONFIG_LABEL}"
RESULTS_DIR="${WORK_DIR}/results/${CONFIG_LABEL}"
if [[ -e "${DATA_DIR}" || -e "${RESULTS_DIR}" ]]; then
  echo "Refusing to mix with an existing run: ${DATA_DIR} or ${RESULTS_DIR}" >&2
  echo "Use a new config label or archive/remove the old run explicitly." >&2
  exit 1
fi

build_targets() {
  local IFS=,
  local strength pct
  printf 'null\n'
  for strength in ${EFFECT_STRENGTHS}; do
    pct=$(awk -v x="${strength}" 'BEGIN{printf "%03d", int(x * 100 + 0.5)}')
    printf 'strength_s%s\n' "${pct}"
  done
}
mapfile -t TARGETS < <(build_targets)

mkdir -p "${WORK_DIR}/logs"
cd "${WORK_DIR}"
gen_job=$(sbatch --parsable "${CODE_DIR}/slurm/submit_generate_power_data_full.sh" "${CONFIG_LABEL}")

declare -a analysis_jobs=()
for label in "${TARGETS[@]}"; do
  sc_job=$(sbatch --parsable --dependency="afterok:${gen_job}" \
    "${CODE_DIR}/slurm/run_sc_hurdle_single.sh" "${CONFIG_LABEL}" "${label}")
  pb_job=$(sbatch --parsable --dependency="afterok:${gen_job}" \
    "${CODE_DIR}/slurm/run_pb_spearman_single.sh" "${CONFIG_LABEL}" "${label}")
  analysis_jobs+=("${sc_job}" "${pb_job}")
done

dependency=$(IFS=:; echo "${analysis_jobs[*]}")
compare_job=$(sbatch --parsable --dependency="afterok:${dependency}" \
  "${CODE_DIR}/slurm/run_compare.sh" "${CONFIG_LABEL}")
printf 'generation_job=%s\nanalysis_jobs=%s\ncompare_job=%s\n' \
  "${gen_job}" "${dependency}" "${compare_job}"
