#!/usr/bin/env bash
#SBATCH --job-name=scpcqtl_qtlsig_l2
#SBATCH --partition=normal
#SBATCH --array=1-22%12
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=logs/qtlsig_l2_%A_%a.out
#SBATCH --error=logs/qtlsig_l2_%A_%a.err

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
: "${SINGLE_DIR:=${MODULE_DIR}/resources/eur_sldsc_custom_annotations/single_maf05}"
: "${QTLSIG_PCQTL_SOURCE:=${MODULE_DIR}/annotations/sig_cisqtl/sig_pcQTL_snps_maf05.tsv}"
: "${QTLSIG_EQTL_SOURCE:=${MODULE_DIR}/annotations/sig_cisqtl/sig_eQTL_snps_maf05.tsv}"
export EUR_PLINK_DIR EUR_QTLSIG_SINGLE_DIR="${SINGLE_DIR}" QTLSIG_PCQTL_SOURCE QTLSIG_EQTL_SOURCE
chr="${SLURM_ARRAY_TASK_ID}"
mkdir -p logs qc/eur_sldsc_annotations
"${LDSC_PYTHON}" "${CODE_DIR}/scripts/54_build_eur_sig_annotations.py" "${chr}"

print_snps_arg=()
snps="${EUR_SLDSC_RESOURCE_DIR}/1000G_Phase3_baselineLD_v2.2/baselineLD.${chr}.snps"
[[ -s "${snps}" ]] && print_snps_arg=(--print-snps "${snps}")

build_one() {
  local ann="$1" stem="$2" dir="${SINGLE_DIR}/$1" prefix="${SINGLE_DIR}/$1/$2.${chr}"
  "${LDSC_PYTHON}" "${LDSC_PY}" --l2 \
    --bfile "${EUR_PLINK_DIR}/1000G.EUR.QC.${chr}" --ld-wind-cm 1 \
    --annot "${dir}/${stem}.${chr}.annot.gz" --out "${prefix}" "${print_snps_arg[@]}"
  "${LDSC_PYTHON}" "${CODE_DIR}/scripts/39_align_single_annotation_ldscores_to_baseline.py" \
    --annotation-name "${ann}" --custom-prefix "${dir}/${stem}." --chrom "${chr}"
}
build_one QTLsig_pcQTL_gt0 qtlsig_pcQTL
build_one QTLsig_eQTL_gt0 qtlsig_eQTL
