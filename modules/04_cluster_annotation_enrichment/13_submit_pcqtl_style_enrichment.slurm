#!/usr/bin/env bash
#SBATCH --job-name=cluster_enrich_pub
#SBATCH --partition=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=cluster_enrich_pub_%j.out
#SBATCH --error=cluster_enrich_pub_%j.err

set -euo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${SC_PCQTL_CLUSTER_ENRICHMENT_ROOT:-${CODE_DIR}}"
export SC_PCQTL_CLUSTER_ENRICHMENT_ROOT="${ROOT}"
FIGURE_OUT="${SC_PCQTL_CLUSTER_ENRICHMENT_FIGURE_OUTPUT_DIR:-${ROOT}/figures/publication}"
cd "${CODE_DIR}"

echo "[INFO] Start: $(date)"
echo "[INFO] Host: $(hostname)"
echo "[INFO] Working directory: ${ROOT}"

mkdir -p "${ROOT}/logs"
RSCRIPT="${COQTL_RSCRIPT:-${R_SCRIPT:-Rscript}}"
"${RSCRIPT}" 06_prepare_annotations.R
"${RSCRIPT}" 07_annotate_cluster_sets.R
"${RSCRIPT}" 08_run_enrichment.R
"${RSCRIPT}" 12_plot_publication_all_strata.R

echo "[INFO] Outputs:"
ls -lh \
  "${ROOT}/results/annotations/cluster_annotation_matrix.tsv" \
  "${ROOT}/results/enrichment/enrichment_by_method.tsv" \
  "${FIGURE_OUT}/cluster_annotation_enrichment_add_cov.pdf"

echo "[INFO] Done: $(date)"
