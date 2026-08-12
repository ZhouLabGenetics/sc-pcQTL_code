#!/bin/bash
# Merge chunks for all chromosomes

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CELL_TYPE="${SC_PCQTL_CELL_TYPE:?Set SC_PCQTL_CELL_TYPE before submission.}"
RUN_R="${COQTL_RSCRIPT:-${R_SCRIPT:-Rscript}}"
WORK_DIR="${SC_PCQTL_HURDLE_WORK_DIR:?Set SC_PCQTL_HURDLE_WORK_DIR before submission.}"
python3 "${SCRIPT_DIR}/../../bin/celltype_eligibility.py" --require-primary "${CELL_TYPE}" >/dev/null

echo "=== Method 2 Step 2b (${CELL_TYPE}): Merge Chunks ==="
echo ""

# Check if chunks exist
CHUNK_DIR="${WORK_DIR}/results/method2_sc_hurdle/gene_associations_chunked"
ASSOC_DIR="${WORK_DIR}/results/method2_sc_hurdle/gene_associations"
if [[ ! -d "${CHUNK_DIR}" ]]; then
  echo "❌ Error: No chunked results found"
  echo "Please run submit_step2_chunked.sh first"
  exit 1
fi

NCHUNKS=$(find "${CHUNK_DIR}" -mindepth 1 -maxdepth 1 -type d | wc -l)
echo "Found ${NCHUNKS} chunks to merge"

if [[ ${NCHUNKS} -eq 0 ]]; then
  echo "❌ No chunks to merge"
  exit 1
fi

echo ""
echo "Merging all chromosomes..."
cd "${SCRIPT_DIR}"
"${RUN_R}" step2b_merge_chunks.R --chr 0

echo ""
echo "=== Complete ===\n"
echo "Merged results available in:"
echo "  ${ASSOC_DIR}/"
echo ""
echo "Next step: bash submit_step3.sh"
