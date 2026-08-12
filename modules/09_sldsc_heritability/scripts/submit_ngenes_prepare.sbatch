#!/usr/bin/env bash
#SBATCH --job-name=scpcqtl_ngenes_prepare
#SBATCH --partition=normal
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=logs/ngenes_prepare_%j.out
#SBATCH --error=logs/ngenes_prepare_%j.err

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
"${SC_PCQTL_PYTHON:-${MODULE_DIR}/tools/ldsc_py3_venv/bin/python}" \
  "${CODE_DIR}/scripts/66_build_pcqtl_ngenes_stratified_snpsets.py"
