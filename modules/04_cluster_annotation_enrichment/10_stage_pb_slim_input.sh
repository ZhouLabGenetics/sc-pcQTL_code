#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RELEASE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
ROOT="${SC_PCQTL_CLUSTER_ENRICHMENT_ROOT:-${SCRIPT_DIR}}"
UP_DATA="${SC_PCQTL_UPSTREAM_DATA_ROOT:-${SC_PCQTL_UPSTREAM_ROOT:?Set SC_PCQTL_UPSTREAM_ROOT to the add-covariate upstream pcQTL workflow root.}/data}"
OUT_DIR="${ROOT}/cache/pb_inputs"
GENE_UNIVERSE="${ROOT}/results/cluster_sets/expressed_gene_universe.tsv"
mkdir -p "${OUT_DIR}"

ELIGIBILITY="${SC_PCQTL_CELLTYPE_MANIFEST:-${RELEASE_ROOT}/config/celltype_eligibility.tsv}"
ELIGIBILITY_HELPER="${RELEASE_ROOT}/bin/celltype_eligibility.py"
mapfile -t CELLTYPES < <(python3 "${ELIGIBILITY_HELPER}" --manifest "${ELIGIBILITY}" --list-primary)
idx="${SLURM_ARRAY_TASK_ID:-}"
ct="${1:-}"
if [[ -z "${ct}" ]]; then
  if [[ -z "${idx}" ]]; then
    echo "Usage: $0 <celltype> or set SLURM_ARRAY_TASK_ID" >&2
    exit 2
  fi
  ct="${CELLTYPES[$((idx-1))]}"
fi
python3 "${ELIGIBILITY_HELPER}" --manifest "${ELIGIBILITY}" --require-primary "${ct}" >/dev/null

src="${UP_DATA}/${ct}_readcounts.tsv"
src_gz="${UP_DATA}/${ct}_readcounts.tsv.gz"
if [[ -e "${src}" ]]; then
  reader=(cat "${src}")
elif [[ -s "${src_gz}" ]]; then
  reader=(zcat "${src_gz}")
else
  echo "Missing source readcount file for ${ct}" >&2
  exit 1
fi

gene_list="${OUT_DIR}/${ct}.keep_genes.txt"
awk -v ct="${ct}" 'BEGIN{FS=OFS="\t"} NR>1 && $1==ct {print $2}' "${GENE_UNIVERSE}" | sort -u > "${gene_list}"
out="${OUT_DIR}/${ct}_pb_input.tsv.gz"
tmp="${out}.tmp"

echo "[stage_pb] ${ct}: extracting $(wc -l < "${gene_list}") genes -> ${out}"
"${reader[@]}" | awk -v genes="${gene_list}" 'BEGIN{
  FS=OFS="\t";
  while ((getline g < genes) > 0) keep_gene[g]=1;
}
NR==1{
  for (i=1; i<=NF; i++) {
    if ($i=="barcode" || $i=="CellID" || $i=="individual" || $i=="IndividualID" || keep_gene[$i]) {
      keep[++n]=i;
    }
  }
  for (j=1; j<=n; j++) printf "%s%s", $keep[j], (j<n?OFS:ORS);
  next;
}
{
  for (j=1; j<=n; j++) printf "%s%s", $keep[j], (j<n?OFS:ORS);
}' | gzip -c > "${tmp}"
mv "${tmp}" "${out}"
echo "[stage_pb] ${ct}: done"
