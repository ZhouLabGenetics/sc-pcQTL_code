#!/bin/bash
# ============================================================================
# run_all.sh
# Master script to run complete pcQTL comparison analysis
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT="${COQTL_RSCRIPT:-${R_SCRIPT:-Rscript}}"

echo "=========================================="
echo "  pcQTL phenotype summary"
echo "=========================================="
echo ""
echo "Working directory: ${SCRIPT_DIR}"
echo "Start time: $(date)"
echo ""

# Step 1: Collect results
echo "=== Step 1: Collecting results from primary-analysis cell types ==="
"${RSCRIPT}" "${SCRIPT_DIR}/01_collect_results.R"
echo ""

# Step 2: Identify significant QTLs
echo "=== Step 2: FDR correction and QTL identification ==="
"${RSCRIPT}" "${SCRIPT_DIR}/02_identify_qtls.R"
echo ""

if [[ "${SC_PCQTL_RUN_EGENE_OVERLAP:-1}" == "1" ]]; then
  echo "=== Step 3: Extracting matched primary-cell-type eGenes ==="
  "${RSCRIPT}" "${SCRIPT_DIR}/03_extract_primary_egenes.R"
  echo ""

  echo "=== Step 4: Summarizing pcQTL/eGene overlap ==="
  "${RSCRIPT}" "${SCRIPT_DIR}/04_summarize_pcqtl_egene_overlap.R"
  echo ""
fi

echo "=========================================="
echo "  Analysis Complete!"
echo "=========================================="
echo ""
echo "End time: $(date)"
echo ""
echo "Results:"
echo "  Data: ${SC_PCQTL_PCQTL_SUMMARY_ROOT}/data/"
echo ""
echo "Key outputs:"
echo "  - all_pcqtl_results.tsv: Results from all primary-analysis cell types"
echo "  - sig_qtls.tsv: Significant QTLs (FDR < 0.05)"
echo "  - ${SC_PCQTL_PCQTL_EGENE_ROOT:-<summary-parent>/pcqtl_compare_saigeqtl}/data/summary_by_celltype.tsv: pcQTL/eGene overlap summary"
echo ""
