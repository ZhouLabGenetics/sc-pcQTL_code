#!/usr/bin/env bash
#SBATCH --job-name=within_donor_perm
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

set -euo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${CODE_DIR}/../shared/module_env.sh"
source "${CODE_DIR}/config/paper_defaults.env"

RUN_ROOT="${COQTL_WITHIN_DONOR_ROOT:-${COQTL_SIM_ROOT}/03_sensitivity_within_donor_permutation}"
DATA_DIR="${RUN_ROOT}/data"
RESULTS_DIR="${RUN_ROOT}/results"
GLOBAL_TEMPLATE="${COQTL_SIM_DATA_ROOT}/shuffle_full/${CONFIG_LABEL}/sc_counts_shuffle_null.rds"
GLOBAL_RESULTS="${COQTL_SIM_NULL_ROOT}/results/${CONFIG_LABEL}"

: "${COQTL_RAW_COUNTS_FILE:?Set COQTL_RAW_COUNTS_FILE}"
[[ -f "${GLOBAL_TEMPLATE}" ]] || { echo "Missing global template: ${GLOBAL_TEMPLATE}" >&2; exit 1; }
[[ -d "${GLOBAL_RESULTS}" ]] || { echo "Missing global results: ${GLOBAL_RESULTS}" >&2; exit 1; }
if [[ -e "${DATA_DIR}" || -e "${RESULTS_DIR}" ]]; then
  echo "Refusing to overwrite an existing within-donor run beneath ${RUN_ROOT}" >&2
  exit 1
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

"${COQTL_RSCRIPT}" "${CODE_DIR}/scripts/generate_within_donor_permutation.R" \
  --raw_counts "${COQTL_RAW_COUNTS_FILE}" \
  --global_template "${GLOBAL_TEMPLATE}" \
  --output_dir "${DATA_DIR}" \
  --seed "${SEED}"

"${COQTL_RSCRIPT}" "${CODE_DIR}/scripts/run_within_donor_tests.R" \
  --input_file "${DATA_DIR}/sc_counts_within_donor_permuted.rds" \
  --output_dir "${RESULTS_DIR}" \
  --workers "${N_WORKERS}" \
  --nonzero_cutoff "${NONZERO_CUTOFF}" \
  --fasthurdle_version "${FASTHURDLE_VERSION}" \
  --fastglm_version "${FASTGLM_VERSION}" \
  --p_thresholds "${P_THRESHOLDS}"

cp "${DATA_DIR}/permutation_validation.tsv" "${RESULTS_DIR}/permutation_validation.tsv"

"${COQTL_RSCRIPT}" "${CODE_DIR}/scripts/plot_global_vs_within_qq.R" \
  --global_dir "${GLOBAL_RESULTS}" \
  --within_dir "${RESULTS_DIR}" \
  --output_dir "${RESULTS_DIR}" \
  --p_thresholds "${P_THRESHOLDS}"

printf 'Completed within-donor permutation diagnostic: %s\n' "${RESULTS_DIR}"
