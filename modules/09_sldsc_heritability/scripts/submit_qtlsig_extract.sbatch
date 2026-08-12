#!/usr/bin/env bash
#SBATCH --job-name=scpcqtl_qtlsig_extract
#SBATCH --partition=normal
#SBATCH --array=1-10%10
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/qtlsig_extract_%x_%A_%a.out
#SBATCH --error=logs/qtlsig_extract_%x_%A_%a.err

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
CODE_DIR="${SC_PCQTL_SLDSC_CODE_DIR:?Set SC_PCQTL_SLDSC_CODE_DIR}"
MODULE_DIR="${SC_PCQTL_SLDSC_ROOT:?Set SC_PCQTL_SLDSC_ROOT}"
QTL_TYPE="${QTLSIG_EXTRACT_QTL:?Set QTLSIG_EXTRACT_QTL to pcQTL or eQTL}"
ELIGIBILITY="${CODE_DIR}/../../config/celltype_eligibility.tsv"
mapfile -t PC_CELLTYPES < <(awk -F '\t' 'NR>1 && $5=="TRUE"{print $1}' "${ELIGIBILITY}")
mapfile -t EQ_CELLTYPES < <(awk -F '\t' 'NR>1 && $5=="TRUE"{print $2}' "${ELIGIBILITY}")
index=$((SLURM_ARRAY_TASK_ID - 1))
[[ ${#PC_CELLTYPES[@]} -eq 10 && ${#EQ_CELLTYPES[@]} -eq 10 ]] || {
  echo "Eligibility manifest must contain exactly 10 primary cell types" >&2
  exit 1
}
if [[ "${QTL_TYPE}" == "pcQTL" ]]; then
  celltype="${PC_CELLTYPES[${index}]}"
elif [[ "${QTL_TYPE}" == "eQTL" ]]; then
  celltype="${EQ_CELLTYPES[${index}]}"
else
  echo "Unknown QTL type: ${QTL_TYPE}" >&2
  exit 2
fi

export SC_PCQTL_SLDSC_ROOT="${MODULE_DIR}"
cd "${MODULE_DIR}"
mkdir -p logs annotations/sig_cisqtl/partials
PYTHON="${SC_PCQTL_PYTHON:-${MODULE_DIR}/tools/ldsc_py3_venv/bin/python}"
"${PYTHON}" "${CODE_DIR}/scripts/53_build_sig_cisqtl_snpsets.py" \
  --qtl "${QTL_TYPE}" --celltype "${celltype}"
