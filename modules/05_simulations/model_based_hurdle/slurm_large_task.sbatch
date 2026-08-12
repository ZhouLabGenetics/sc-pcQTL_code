#!/usr/bin/env bash
#SBATCH --job-name=scpcqtl_sim_task
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=10:00:00
#SBATCH --output=logs/slurm_task_%A_%a.out
#SBATCH --error=logs/slurm_task_%A_%a.err

set -euo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CODE_DIR}/../shared/module_env.sh"
: "${SC_PCQTL_MODEL_SIM_WORK_ROOT:?Set SC_PCQTL_MODEL_SIM_WORK_ROOT}"
WORK_DIR="${SC_PCQTL_MODEL_SIM_WORK_ROOT}"

MANIFEST="${1:?Pass the simulation manifest path}"
[[ -s "${MANIFEST}" ]] || { echo "Manifest not found: ${MANIFEST}" >&2; exit 1; }

line=$(awk -F'\t' -v idx="${SLURM_ARRAY_TASK_ID}" 'NR == idx + 1 { print $0 }' "${MANIFEST}")
[[ -n "${line}" ]] || { echo "No manifest entry for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IFS=$'\t' read -r scenario rep <<< "${line}"

rep_dir=$(printf "replicate_%03d" "${rep}")
sim_dir="${WORK_DIR}/simulations_large/${scenario}/${rep_dir}"
hurdle_dir="${WORK_DIR}/analysis_large_compare/hurdle/${scenario}/${rep_dir}"
donor_dir="${WORK_DIR}/analysis_large_compare/donor_agg/${scenario}/${rep_dir}"
mkdir -p "${WORK_DIR}/logs" "${WORK_DIR}/simulations_large" "${WORK_DIR}/analysis_large_compare"

if [[ ! -s "${sim_dir}/simulation_metadata.tsv" ]]; then
  "${COQTL_RSCRIPT}" "${CODE_DIR}/02_simulate_hurdle_counts.R" \
    --reference_rds "${WORK_DIR}/reference_large/reference_stats.rds" \
    --scenario "${scenario}" \
    --replicate "${rep}" \
    --n_donors 300 \
    --n_genes 100 \
    --signal_pair_fraction 0.05 \
    --out_dir "${WORK_DIR}/simulations_large"
fi

if [[ ! -s "${hurdle_dir}/results/method2_sc_hurdle/gene_associations_chunked/chr1_chunk001/summary.tsv" ]]; then
  "${COQTL_RSCRIPT}" "${CODE_DIR}/03_run_hurdle_association.R" --sim_dir "${sim_dir}" --out_dir "${hurdle_dir}"
fi

if [[ ! -s "${donor_dir}/results/donor_agg/summary.tsv" ]]; then
  "${COQTL_RSCRIPT}" "${CODE_DIR}/06_run_donor_aggregate_cor.R" --sim_dir "${sim_dir}" --out_dir "${donor_dir}"
fi
