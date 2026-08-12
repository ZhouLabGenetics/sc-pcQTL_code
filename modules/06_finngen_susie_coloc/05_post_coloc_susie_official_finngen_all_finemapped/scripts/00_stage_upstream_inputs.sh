#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RELEASE_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
ELIGIBILITY="${SC_PCQTL_CELLTYPE_MANIFEST:-${RELEASE_ROOT}/config/celltype_eligibility.tsv}"
ELIGIBILITY_HELPER="${RELEASE_ROOT}/bin/celltype_eligibility.py"
mapfile -t PRIMARY_CELLTYPES < <(python3 "${ELIGIBILITY_HELPER}" --manifest "${ELIGIBILITY}" --list-primary)
mapfile -t EXCLUDED_CELLTYPES < <(python3 "${ELIGIBILITY_HELPER}" --manifest "${ELIGIBILITY}" --list-excluded)
if [[ ${#PRIMARY_CELLTYPES[@]} -eq 0 ]]; then
  echo "No primary cell types in ${ELIGIBILITY}" >&2
  exit 1
fi

ROOT_DIR="${SC_PCQTL_FORMAL_COLOC_ROOT:?Set SC_PCQTL_FORMAL_COLOC_ROOT to the external formal-colocalization work directory.}"
PCQTL_ROOT="${SC_PCQTL_UPSTREAM_ROOT:?Set SC_PCQTL_UPSTREAM_ROOT to the add-covariate upstream pcQTL workflow root.}"
EQTL_ROOT="${ONEK1K_EQTL_ROOT:?Set ONEK1K_EQTL_ROOT to the OneK1K single-gene eQTL result directory.}"

CACHE_DIR="${ROOT_DIR}/05_post_coloc_susie_official_finngen_all_finemapped/cache"
ASSIGN_DIR="${CACHE_DIR}/cluster_gene_assignments"
LOAD_DIR="${CACHE_DIR}/pc_loadings"
NOM_DIR="${CACHE_DIR}/nominal_eqtl"
mkdir -p "${ASSIGN_DIR}" "${LOAD_DIR}" "${NOM_DIR}"

declare -A EQTL_TAR=()
declare -A PRIMARY_SET=()
while IFS=$'\t' read -r ct eqtl_ct; do
  [[ -n "${ct}" && -n "${eqtl_ct}" ]] || continue
  PRIMARY_SET["${ct}"]=1
  EQTL_TAR["${ct}"]="cis_${eqtl_ct}.tar.gz"
done < <(python3 "${ELIGIBILITY_HELPER}" --manifest "${ELIGIBILITY}" --list-primary-mapping)

# Remove excluded entries from the generated staging cache so rerunning in an
# Recreate staged inputs so excluded-cell artifacts cannot persist.
for ct in "${EXCLUDED_CELLTYPES[@]}"; do
  rm -f "${ASSIGN_DIR}/${ct}.cluster_gene_assignments.tsv" \
        "${CACHE_DIR}/${ct}.nominal_members.tsv" \
        "${CACHE_DIR}/${ct}.nominal_members.list"
  rm -rf "${LOAD_DIR:?}/${ct}" "${NOM_DIR:?}/${ct}" "${CACHE_DIR}/tmp_extract_${ct}"
done

echo "[stage] Copying cluster gene assignments"
for ct_dir in "${PCQTL_ROOT}/celltypes/"*; do
  [[ -d "${ct_dir}" ]] || continue
  ct="$(basename "${ct_dir}")"
  [[ -n "${PRIMARY_SET[${ct}]:-}" ]] || continue
  assign="${ct_dir}/cluster_identification/results/method2_sc_hurdle/merged_clusters/cluster_gene_assignments.tsv"
  if [[ -s "${assign}" ]]; then
    cp -f "${assign}" "${ASSIGN_DIR}/${ct}.cluster_gene_assignments.tsv"
  fi
done

echo "[stage] Building required cluster and gene lists from non-empty QTL credible sets"
required="${CACHE_DIR}/required_nominal_genes.tsv"
required_clusters="${CACHE_DIR}/required_clusters.tsv"
primary_list="${CACHE_DIR}/primary_celltypes.list"
: > "${required}"
printf "%s\n" "${PRIMARY_CELLTYPES[@]}" > "${primary_list}"
find "${ROOT_DIR}/results/fine_mapping/qtl" -name "*.credible_sets.tsv" -type f -size +200c -print |
  awk -v prefix="${ROOT_DIR}/results/fine_mapping/qtl/" -v allowed="${primary_list}" \
    'BEGIN{FS="/"; OFS="\t"; while ((getline x < allowed) > 0) keep[x]=1} {sub(prefix, "", $0); split($0, a, "/"); if (a[1] in keep) print a[1], a[2]}' |
  sort -u > "${required_clusters}"

echo "[stage] Copying PC loadings for required clusters"
while IFS=$'\t' read -r ct cluster; do
  src="${PCQTL_ROOT}/celltypes/${ct}/pcQTL/step2_pca/${cluster}/gene_loadings.tsv"
  if [[ -s "${src}" ]]; then
    mkdir -p "${LOAD_DIR}/${ct}/${cluster}"
    cp -f "${src}" "${LOAD_DIR}/${ct}/${cluster}/gene_loadings.tsv"
  fi
done < "${required_clusters}"

while IFS=$'\t' read -r ct cluster; do
  [[ -n "${ct}" && -n "${cluster}" ]] || continue
  assign="${ASSIGN_DIR}/${ct}.cluster_gene_assignments.tsv"
  [[ -s "${assign}" ]] || continue
  awk -v ct="${ct}" -v cl="${cluster}" 'BEGIN{FS=OFS="\t"} NR>1 && $1==cl {print ct, $1, $4}' "${assign}"
done < "${required_clusters}" | sort -u > "${required}"

echo "[stage] Extracting nominal eQTL files needed for gene-level effects"
extract_ct() {
  local ct="$1"
  local tar_name="${EQTL_TAR[${ct}]:-}"
  [[ -n "${tar_name}" ]] || return 0
  local tar_path="${EQTL_ROOT}/${tar_name}"
  [[ -s "${tar_path}" ]] || return 0
  local inner="${tar_name%.tar.gz}"
  local suffix="${inner#cis_}"
  local out_dir="${NOM_DIR}/${ct}"
  mkdir -p "${out_dir}"

  local member_map="${CACHE_DIR}/${ct}.nominal_members.tsv"
  awk -v ct="${ct}" -v inner="${inner}" -v suffix="${suffix}" 'BEGIN{FS=OFS="\t"} $1==ct {gene=$3; member=inner "/" gene "_" suffix "_count_saigeqtl_cis_window_1000000.singleVar.txt"; print gene, member}' "${required}" > "${member_map}"
  [[ -s "${member_map}" ]] || return 0

  local member_list="${CACHE_DIR}/${ct}.nominal_members.list"
  cut -f2 "${member_map}" > "${member_list}"
  local tmp_dir="${CACHE_DIR}/tmp_extract_${ct}"
  rm -rf "${tmp_dir}"
  mkdir -p "${tmp_dir}"

  tar --ignore-failed-read -xzf "${tar_path}" -C "${tmp_dir}" -T "${member_list}" 2>/dev/null || true
  while IFS=$'\t' read -r gene member; do
    local src="${tmp_dir}/${member}"
    local out_file="${out_dir}/${gene}.singleVar.txt"
    [[ -s "${src}" ]] && cp -f "${src}" "${out_file}"
  done < "${member_map}"
  rm -rf "${tmp_dir}"
}

max_jobs=4
for ct in "${PRIMARY_CELLTYPES[@]}"; do
  extract_ct "${ct}" &
  while (( $(jobs -rp | wc -l) >= max_jobs )); do
    wait -n || true
  done
done
wait || true

echo "[stage] Done"
echo "[stage] Required gene rows: $(wc -l < "${required}")"
echo "[stage] Cached nominal files: $(find "${NOM_DIR}" -name '*.singleVar.txt' | wc -l)"
