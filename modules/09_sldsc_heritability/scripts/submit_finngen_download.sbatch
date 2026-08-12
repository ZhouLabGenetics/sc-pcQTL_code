#!/usr/bin/env bash
#SBATCH --job-name=fg247_download
#SBATCH --partition=filemove
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/finngen_download_%A_%a.out
#SBATCH --error=logs/finngen_download_%A_%a.err

set -euo pipefail
: "${SC_PCQTL_SLDSC_CODE_DIR:?Set SC_PCQTL_SLDSC_CODE_DIR}"
: "${FINNGEN_DOWNLOAD_MANIFEST:?Set FINNGEN_DOWNLOAD_MANIFEST}"
python3 "${SC_PCQTL_SLDSC_CODE_DIR}/scripts/75_download_finngen_sumstats.py" \
  "${SLURM_ARRAY_TASK_ID}" --manifest "${FINNGEN_DOWNLOAD_MANIFEST}"
