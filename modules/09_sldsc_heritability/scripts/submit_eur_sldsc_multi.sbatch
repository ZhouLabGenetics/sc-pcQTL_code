#!/usr/bin/env bash
#SBATCH --job-name=scpcqtl_sldsc
#SBATCH --partition=normal
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/eur_sldsc_multi_%A_%a.out
#SBATCH --error=logs/eur_sldsc_multi_%A_%a.err

set -euo pipefail
CODE_DIR="${SC_PCQTL_SLDSC_CODE_DIR:-}"
[[ -d "${CODE_DIR}/scripts" ]] || {
  echo "Set SC_PCQTL_SLDSC_CODE_DIR to the absolute module 09 code directory" >&2
  exit 1
}
MODULE_DIR="${SC_PCQTL_SLDSC_ROOT:-${CODE_DIR}}"
export SC_PCQTL_SLDSC_ROOT="${MODULE_DIR}"
cd "${MODULE_DIR}"
: "${QTLSIG_MULTI_MANIFEST:?Set QTLSIG_MULTI_MANIFEST to a task manifest}"
: "${LDSC_PYTHON:=${MODULE_DIR}/tools/ldsc_py3_venv/bin/python}"
mkdir -p logs
"${LDSC_PYTHON}" "${CODE_DIR}/scripts/62_run_eur_sldsc_multi_task.py" "${SLURM_ARRAY_TASK_ID}" --manifest "${QTLSIG_MULTI_MANIFEST}"
