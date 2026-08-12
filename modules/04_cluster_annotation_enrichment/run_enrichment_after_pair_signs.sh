#!/usr/bin/env bash
set -euo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${SC_PCQTL_CLUSTER_ENRICHMENT_ROOT:-${CODE_DIR}}"
export SC_PCQTL_CLUSTER_ENRICHMENT_ROOT="${ROOT}"
cd "${CODE_DIR}"
RSCRIPT="${COQTL_RSCRIPT:-${R_SCRIPT:-Rscript}}"

"${RSCRIPT}" 04_build_null_sets.R
bash 05_download_annotations.sh
"${RSCRIPT}" 06_prepare_annotations.R
"${RSCRIPT}" 07_annotate_cluster_sets.R
"${RSCRIPT}" 08_run_enrichment.R
"${RSCRIPT}" 12_plot_publication_all_strata.R
