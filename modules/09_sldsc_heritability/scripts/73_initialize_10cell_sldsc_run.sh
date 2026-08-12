#!/usr/bin/env bash
set -euo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASE_ROOT="${SC_PCQTL_SLDSC_BASE_ROOT:?Set SC_PCQTL_SLDSC_BASE_ROOT to the existing resource tree}"
RUN_ROOT="${SC_PCQTL_SLDSC_ROOT:?Set SC_PCQTL_SLDSC_ROOT to a new versioned run directory}"
ENDPOINT_MANIFEST="${FINNGEN_SLDSC_ENDPOINT_MANIFEST:?Set FINNGEN_SLDSC_ENDPOINT_MANIFEST to the prespecified 247-trait manifest}"

[[ "$(realpath -m "${BASE_ROOT}")" != "$(realpath -m "${RUN_ROOT}")" ]] || {
  echo "RUN_ROOT must differ from BASE_ROOT" >&2
  exit 1
}

ANNOTATION_CACHE_ROOT="${SC_PCQTL_SLDSC_ANNOTATION_CACHE_ROOT:-}"
mkdir -p "${RUN_ROOT}"/{figures,logs,manifests,qc,results/heritability_enrichment,resources}

link_shared() {
  local source="$1" target="$2"
  [[ -e "${source}" ]] || { echo "Missing shared resource: ${source}" >&2; exit 1; }
  if [[ -L "${target}" ]]; then
    [[ "$(readlink -f "${target}")" == "$(readlink -f "${source}")" ]] || {
      echo "Existing link points elsewhere: ${target}" >&2
      exit 1
    }
  elif [[ -e "${target}" ]]; then
    echo "Refusing to replace existing path: ${target}" >&2
    exit 1
  else
    ln -s "${source}" "${target}"
  fi
}

link_shared "${BASE_ROOT}/resources/eur_sldsc_ref" "${RUN_ROOT}/resources/eur_sldsc_ref"
link_shared "${BASE_ROOT}/resources/finngen_r12_sldsc_sumstats" "${RUN_ROOT}/resources/finngen_r12_sldsc_sumstats"
link_shared "${BASE_ROOT}/tools" "${RUN_ROOT}/tools"
if [[ -n "${ANNOTATION_CACHE_ROOT}" ]]; then
  link_shared "${ANNOTATION_CACHE_ROOT}/annotations" "${RUN_ROOT}/annotations"
  link_shared "${ANNOTATION_CACHE_ROOT}/resources/eur_sldsc_custom_annotations" \
    "${RUN_ROOT}/resources/eur_sldsc_custom_annotations"
else
  mkdir -p "${RUN_ROOT}/annotations"
fi
endpoint_count="$(awk 'END{print NR-1}' "${ENDPOINT_MANIFEST}")"
[[ "${endpoint_count}" == 247 ]] || {
  echo "Expected the prespecified 247-trait universe; found ${endpoint_count} rows" >&2
  exit 1
}
cp "${ENDPOINT_MANIFEST}" "${RUN_ROOT}/manifests/finngen_r12_sldsc_prespecified247_endpoints.tsv"

eligibility="${CODE_DIR}/../../config/celltype_eligibility.tsv"
cp "${eligibility}" "${RUN_ROOT}/manifests/celltype_eligibility.tsv"
{
  printf 'field\tvalue\n'
  printf 'run_root\t%s\n' "$(realpath -m "${RUN_ROOT}")"
  printf 'base_resource_root\t%s\n' "$(realpath "${BASE_ROOT}")"
  printf 'code_dir\t%s\n' "$(realpath "${CODE_DIR}")"
  printf 'endpoint_manifest_sha256\t%s\n' "$(sha256sum "${ENDPOINT_MANIFEST}" | awk '{print $1}')"
  printf 'endpoint_universe\tprespecified_core_disease_traits\n'
  printf 'endpoint_universe_count\t%s\n' "${endpoint_count}"
  printf 'gwas_finemapping_used_for_endpoint_selection\tFALSE\n'
  printf 'annotation_cache_root\t%s\n' "${ANNOTATION_CACHE_ROOT:-NONE}"
  printf 'eligibility_manifest_sha256\t%s\n' "$(sha256sum "${eligibility}" | awk '{print $1}')"
  printf 'celltype_minimum_cells\t10000\n'
  printf 'primary_celltype_count\t%s\n' "$(awk -F '\t' 'NR>1 && $5=="TRUE"{n++} END{print n+0}' "${eligibility}")"
} > "${RUN_ROOT}/run_provenance.tsv"

echo "Initialized ${RUN_ROOT}"
