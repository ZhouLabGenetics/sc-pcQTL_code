#!/usr/bin/env bash
set -euo pipefail

CODE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ROOT_DIR="${SC_PCQTL_FORMAL_COLOC_ROOT:?Set SC_PCQTL_FORMAL_COLOC_ROOT to the external formal-colocalization work directory.}"
LOG_DIR="${ROOT_DIR}/05_post_coloc_susie_official_finngen_all_finemapped/logs"
mkdir -p "${LOG_DIR}"

job_id=$(sbatch --parsable \
  --partition=normal \
  --job-name=post_susie_fg \
  --output="${LOG_DIR}/post_susie_fg_%j.out" \
  --error="${LOG_DIR}/post_susie_fg_%j.err" \
  --time=06:00:00 \
  --mem=32G \
  --cpus-per-task=4 \
  --wrap="bash ${CODE_ROOT}/05_post_coloc_susie_official_finngen_all_finemapped/scripts/run_post_coloc_susie_official_finngen.sh")

echo "Post-coloc job: ${job_id}"
