#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RELEASE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
ELIGIBILITY="${SC_PCQTL_CELLTYPE_MANIFEST:-${RELEASE_ROOT}/config/celltype_eligibility.tsv}"
ELIGIBILITY_HELPER="${RELEASE_ROOT}/bin/celltype_eligibility.py"
ROOT="${SC_PCQTL_CLUSTER_ENRICHMENT_ROOT:-${SCRIPT_DIR}}"
SRC_ADD="${SC_PCQTL_UPSTREAM_ROOT:?Set SC_PCQTL_UPSTREAM_ROOT to the add-covariate upstream pcQTL workflow root.}"
CACHE="${ROOT}/cache"
DST_ADD="${CACHE}/01_upstream_main_pipeline_add_cov"

mkdir -p "${DST_ADD}/celltypes" "${CACHE}"

mapfile -t PRIMARY_CELLTYPES < <(python3 "${ELIGIBILITY_HELPER}" --manifest "${ELIGIBILITY}" --list-primary)
if [[ ${#PRIMARY_CELLTYPES[@]} -eq 0 ]]; then
  echo "[stage] No primary cell types in ${ELIGIBILITY}" >&2
  exit 1
fi
declare -A PRIMARY_SET=()
for ct in "${PRIMARY_CELLTYPES[@]}"; do PRIMARY_SET["${ct}"]=1; done

# Remove only generated cache entries that are no longer eligible. Source data
# and eligible staged data are never deleted.
for dst_ct in "${DST_ADD}/celltypes/"*; do
  [[ -d "${dst_ct}" ]] || continue
  ct="$(basename "${dst_ct}")"
  if [[ -z "${PRIMARY_SET[${ct}]:-}" ]]; then
    echo "[stage] Removing excluded stale cache: ${ct}"
    rm -rf "${dst_ct}"
  fi
done

echo "[stage] Copying add-cov cluster inputs and PCA files for ${#PRIMARY_CELLTYPES[@]} primary cell types"
for ct_dir in "${SRC_ADD}/celltypes/"*; do
  [[ -d "${ct_dir}" ]] || continue
  ct="$(basename "${ct_dir}")"
  [[ -n "${PRIMARY_SET[${ct}]:-}" ]] || continue
  dst_ct="${DST_ADD}/celltypes/${ct}"
  mkdir -p "${dst_ct}/cluster_identification/results/method2_sc_hurdle" "${dst_ct}/pcQTL/step2_pca"
  for f in \
    "${ct_dir}/cluster_identification/results/method2_sc_hurdle/filtered_genes.tsv" \
    "${ct_dir}/cluster_identification/results/method2_sc_hurdle/merged_clusters/cluster_summary.tsv" \
    "${ct_dir}/cluster_identification/results/method2_sc_hurdle/merged_clusters/cluster_gene_assignments.tsv"; do
    if [[ -s "${f}" ]]; then
      rel="${f#${ct_dir}/}"
      mkdir -p "${dst_ct}/$(dirname "${rel}")"
      cp -f "${f}" "${dst_ct}/${rel}"
    fi
  done
  if [[ -d "${ct_dir}/cluster_identification/results/method2_sc_hurdle/gene_associations_chunked" ]]; then
    mkdir -p "${dst_ct}/cluster_identification/results/method2_sc_hurdle/gene_associations_chunked"
    find "${ct_dir}/cluster_identification/results/method2_sc_hurdle/gene_associations_chunked" -name significant_pairs.tsv -print0 |
      while IFS= read -r -d '' f; do
        rel="${f#${ct_dir}/}"
        mkdir -p "${dst_ct}/$(dirname "${rel}")"
        cp -f "${f}" "${dst_ct}/${rel}"
      done
  fi
  if [[ -d "${ct_dir}/pcQTL/step2_pca" ]]; then
    find "${ct_dir}/pcQTL/step2_pca" -mindepth 2 -maxdepth 2 \( -name pca_results.rda -o -name gene_loadings.tsv \) -print0 |
      while IFS= read -r -d '' f; do
        rel="${f#${ct_dir}/}"
        mkdir -p "${dst_ct}/$(dirname "${rel}")"
        cp -f "${f}" "${dst_ct}/${rel}"
      done
  fi
done

echo "[stage] Done"
