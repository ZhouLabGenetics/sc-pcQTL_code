#!/usr/bin/env bash
#SBATCH --job-name=gimap_coloc_susie
#SBATCH --partition=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00

set -euo pipefail

CODE_ROOT="${SC_PCQTL_CODE_ROOT:?Set SC_PCQTL_CODE_ROOT.}"
MIX_ROOT="${SC_PCQTL_GIMAP_MIXING_ROOT:?Set SC_PCQTL_GIMAP_MIXING_ROOT.}"
RSCRIPT="${COQTL_RSCRIPT:-${R_SCRIPT:-Rscript}}"
SCRIPT="${CODE_ROOT}/modules/10_representative_loci/gimap_celltype_specificity/mixing/pseudobulk_tensorqtl/scripts/06_coloc_mixture_vs_lymphocyte_susie.R"

mkdir -p "${MIX_ROOT}/coloc_susie" "${MIX_ROOT}/logs"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

"${RSCRIPT}" "${SCRIPT}" \
  --sumstats "${MIX_ROOT}/coloc_susie/mix_all_sumstats.tsv" \
  --out "${MIX_ROOT}/coloc_susie/mix_coloc_susie_summary.tsv"
