#!/usr/bin/env bash
#SBATCH --job-name=scpcqtl_qtlsig_prepare
#SBATCH --partition=normal
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=logs/qtlsig_prepare_%j.out
#SBATCH --error=logs/qtlsig_prepare_%j.err

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
CODE_DIR="${SC_PCQTL_SLDSC_CODE_DIR:-}"
[[ -d "${CODE_DIR}/scripts" ]] || {
  echo "Set SC_PCQTL_SLDSC_CODE_DIR to the absolute module 09 code directory" >&2
  exit 1
}
MODULE_DIR="${SC_PCQTL_SLDSC_ROOT:-${CODE_DIR}}"
export SC_PCQTL_SLDSC_ROOT="${MODULE_DIR}"
cd "${MODULE_DIR}"
mkdir -p logs
PYTHON="${SC_PCQTL_PYTHON:-${MODULE_DIR}/tools/ldsc_py3_venv/bin/python}"
"${PYTHON}" "${CODE_DIR}/scripts/53_build_sig_cisqtl_snpsets.py" --merge
test -s annotations/sig_cisqtl/sig_pcQTL_snps_maf05.tsv
test -s annotations/sig_cisqtl/sig_eQTL_snps_maf05.tsv
if [[ ! -s resources/eur_sldsc_ref/baselineLD_noQTL/baselineLD_noQTL.1.l2.ldscore.gz ]]; then
  "${PYTHON}" "${CODE_DIR}/scripts/30_build_baselineLD_noQTL_resource.py"
fi
