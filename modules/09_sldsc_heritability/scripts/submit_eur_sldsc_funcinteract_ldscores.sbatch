#!/usr/bin/env bash
#SBATCH --job-name=scpcqtl_function_l2
#SBATCH --partition=normal
#SBATCH --array=1-22%12
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=logs/function_l2_%A_%a.out
#SBATCH --error=logs/function_l2_%A_%a.err

set -euo pipefail
CODE_DIR="${SC_PCQTL_SLDSC_CODE_DIR:-}"
[[ -d "${CODE_DIR}/scripts" ]] || {
  echo "Set SC_PCQTL_SLDSC_CODE_DIR to the absolute module 09 code directory" >&2
  exit 1
}
MODULE_DIR="${SC_PCQTL_SLDSC_ROOT:-${CODE_DIR}}"
export SC_PCQTL_SLDSC_ROOT="${MODULE_DIR}"
cd "${MODULE_DIR}"
: "${EUR_SLDSC_RESOURCE_DIR:=${MODULE_DIR}/resources/eur_sldsc_ref}"
[[ -s "${EUR_SLDSC_RESOURCE_DIR}/env.sh" ]] && source "${EUR_SLDSC_RESOURCE_DIR}/env.sh"
: "${LDSC_PYTHON:=${MODULE_DIR}/tools/ldsc_py3_venv/bin/python}"
: "${LDSC_PY:=${MODULE_DIR}/tools/ldsc/ldsc.py}"
: "${EUR_PLINK_DIR:=${EUR_SLDSC_RESOURCE_DIR}/1000G_EUR_Phase3_plink}"
FI="${MODULE_DIR}/resources/eur_sldsc_custom_annotations/single_maf05_funcinteract"
chr="${SLURM_ARRAY_TASK_ID}"
mkdir -p logs qc/eur_sldsc_annotations
"${LDSC_PYTHON}" "${CODE_DIR}/scripts/65_build_functional_interaction_annotations.py" "${chr}"

print_snps_arg=()
snps="${EUR_SLDSC_RESOURCE_DIR}/1000G_Phase3_baselineLD_v2.2/baselineLD.${chr}.snps"
[[ -s "${snps}" ]] && print_snps_arg=(--print-snps "${snps}")
build_one() {
  local ann="$1" stem="$2" dir="${FI}/$1" prefix="${FI}/$1/$2.${chr}"
  "${LDSC_PYTHON}" "${LDSC_PY}" --l2 \
    --bfile "${EUR_PLINK_DIR}/1000G.EUR.QC.${chr}" --ld-wind-cm 1 \
    --annot "${dir}/${stem}.${chr}.annot.gz" --out "${prefix}" "${print_snps_arg[@]}"
  "${LDSC_PYTHON}" "${CODE_DIR}/scripts/39_align_single_annotation_ldscores_to_baseline.py" \
    --annotation-name "${ann}" --custom-prefix "${dir}/${stem}." --chrom "${chr}"
}
build_one pcQTL_prom fi_pcQTL_prom
build_one pcQTL_enh fi_pcQTL_enh
build_one eQTL_prom fi_eQTL_prom
build_one eQTL_enh fi_eQTL_enh
