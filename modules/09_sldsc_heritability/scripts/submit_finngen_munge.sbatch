#!/usr/bin/env bash
#SBATCH --job-name=fg247_munge
#SBATCH --partition=normal
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/finngen_munge_%A_%a.out
#SBATCH --error=logs/finngen_munge_%A_%a.err

set -euo pipefail
: "${SC_PCQTL_SLDSC_CODE_DIR:?Set SC_PCQTL_SLDSC_CODE_DIR}"
: "${SC_PCQTL_SLDSC_ROOT:?Set SC_PCQTL_SLDSC_ROOT}"
: "${FINNGEN_MUNGE_MANIFEST:?Set FINNGEN_MUNGE_MANIFEST}"
cd "${SC_PCQTL_SLDSC_ROOT}"
"${LDSC_PYTHON:-${SC_PCQTL_SLDSC_ROOT}/tools/ldsc_py3_venv/bin/python}" \
  "${SC_PCQTL_SLDSC_CODE_DIR}/scripts/76_munge_finngen_r12_for_ldsc.py" \
  "${SLURM_ARRAY_TASK_ID}" --manifest "${FINNGEN_MUNGE_MANIFEST}"
