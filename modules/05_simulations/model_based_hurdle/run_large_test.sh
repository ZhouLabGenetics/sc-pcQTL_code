#!/usr/bin/env bash
set -euo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CODE_DIR}/../shared/module_env.sh"
: "${SC_PCQTL_MODEL_SIM_WORK_ROOT:?Set SC_PCQTL_MODEL_SIM_WORK_ROOT}"
: "${COQTL_RAW_COUNTS_FILE:?Set COQTL_RAW_COUNTS_FILE}"

WORK_DIR="${SC_PCQTL_MODEL_SIM_WORK_ROOT}"
REFERENCE_DIR="${WORK_DIR}/reference_large"
SIM_ROOT="${WORK_DIR}/simulations_large"
ANALYSIS_ROOT="${WORK_DIR}/analysis_large_compare"
SOURCE_ROOT="${WORK_DIR}/plots_target"

N_REFERENCE_CELLS="${N_REFERENCE_CELLS:-20000}"
N_REFERENCE_GENES="${N_REFERENCE_GENES:-3000}"
N_DONORS="${N_DONORS:-300}"
N_GENES="${N_GENES:-100}"
SIGNAL_PAIR_FRACTION="${SIGNAL_PAIR_FRACTION:-0.05}"
N_REPS="${N_REPS:-50}"
SCENARIOS=(
  null zero_only_low count_only_low both_low
  zero_only_high count_only_high both_high
)

mkdir -p "${REFERENCE_DIR}" "${SIM_ROOT}" "${ANALYSIS_ROOT}" "${SOURCE_ROOT}"

if [[ ! -s "${REFERENCE_DIR}/reference_stats.rds" ]]; then
  "${COQTL_RSCRIPT}" "${CODE_DIR}/01_fit_reference_from_real_data.R" \
    --count_file "${COQTL_RAW_COUNTS_FILE}" \
    --n_reference_cells "${N_REFERENCE_CELLS}" \
    --n_reference_genes "${N_REFERENCE_GENES}" \
    --out_dir "${REFERENCE_DIR}"
fi

for scenario in "${SCENARIOS[@]}"; do
  for rep in $(seq 1 "${N_REPS}"); do
    rep_dir=$(printf "replicate_%03d" "${rep}")
    sim_dir="${SIM_ROOT}/${scenario}/${rep_dir}"
    hurdle_dir="${ANALYSIS_ROOT}/hurdle/${scenario}/${rep_dir}"
    donor_dir="${ANALYSIS_ROOT}/donor_agg/${scenario}/${rep_dir}"

    if [[ ! -s "${sim_dir}/simulation_metadata.tsv" ]]; then
      "${COQTL_RSCRIPT}" "${CODE_DIR}/02_simulate_hurdle_counts.R" \
        --reference_rds "${REFERENCE_DIR}/reference_stats.rds" \
        --scenario "${scenario}" \
        --replicate "${rep}" \
        --n_donors "${N_DONORS}" \
        --n_genes "${N_GENES}" \
        --signal_pair_fraction "${SIGNAL_PAIR_FRACTION}" \
        --out_dir "${SIM_ROOT}"
    fi

    if [[ ! -s "${hurdle_dir}/results/method2_sc_hurdle/gene_associations_chunked/chr1_chunk001/summary.tsv" ]]; then
      "${COQTL_RSCRIPT}" "${CODE_DIR}/03_run_hurdle_association.R" \
        --sim_dir "${sim_dir}" --out_dir "${hurdle_dir}"
    fi

    if [[ ! -s "${donor_dir}/results/donor_agg/summary.tsv" ]]; then
      "${COQTL_RSCRIPT}" "${CODE_DIR}/06_run_donor_aggregate_cor.R" \
        --sim_dir "${sim_dir}" --out_dir "${donor_dir}"
    fi
  done
done

"${COQTL_RSCRIPT}" "${CODE_DIR}/10_make_target_figures.R" \
  --sim_root "${SIM_ROOT}" \
  --analysis_root "${ANALYSIS_ROOT}" \
  --out_dir "${SOURCE_ROOT}"
