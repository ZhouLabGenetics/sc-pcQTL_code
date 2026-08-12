#!/usr/bin/env bash
set -euo pipefail

SELECTOR="${1:?Usage: $0 <job_id|all|scenario_name> [cis|cis_nominal|both]}"
MODE="${2:-both}"
BASE="${SC_PCQTL_GIMAP_MIXING_ROOT:?Set SC_PCQTL_GIMAP_MIXING_ROOT.}"
MANIFEST="${BASE}/data/tensorqtl_job_manifest.tsv"
[[ "${MODE}" =~ ^(cis|cis_nominal|both)$ ]] || { echo "ERROR: mode must be cis, cis_nominal, or both." >&2; exit 2; }
[[ -s "${MANIFEST}" ]] || { echo "ERROR: missing manifest ${MANIFEST}" >&2; exit 1; }
GENO_PREFIX="$(cat "${BASE}/genotypes/gimap_tensorqtl_plink_prefix.txt")"
PY="${TENSORQTL_PYTHON:-${BASE}/env/tensorqtl_venv/bin/python}"
mkdir -p "${BASE}/logs/tensorqtl" "${BASE}/tensorqtl/cis" "${BASE}/tensorqtl/cis_nominal"

if [[ ! -x "${PY}" ]]; then
  echo "ERROR: TensorQTL python not found at ${PY}. Run scripts/00_setup_tensorqtl_env.sh first or set TENSORQTL_PYTHON." >&2
  exit 1
fi
"${PY}" - <<'PY'
import tensorqtl
PY

run_one() {
  local jid="$1"
  local line scenario pset bed covars
  line="$(awk -F'\t' -v id="${jid}" 'NR>1 && $1==id {print; exit}' "${MANIFEST}")"
  if [[ -z "${line}" ]]; then
    echo "ERROR: job_id ${jid} not found" >&2
    return 1
  fi
  IFS=$'\t' read -r _ scenario pset bed covars <<< "${line}"
  echo "============================================================"
  echo "TensorQTL job ${jid}: ${scenario} / ${pset} / mode=${MODE}"
  echo "phenotypes: ${bed}"
  echo "covariates: ${covars}"
  echo "genotype: ${GENO_PREFIX}"
  echo "start: $(date)"
  echo "============================================================"

  if [[ "${MODE}" == "cis" || "${MODE}" == "both" ]]; then
    local cis_dir="${BASE}/tensorqtl/cis/${scenario}"
    mkdir -p "${cis_dir}"
    local prefix="${cis_dir}/${pset}"
    if ls "${prefix}"*.cis_qtl.txt.gz >/dev/null 2>&1 || ls "${prefix}"*.cis_qtl.txt >/dev/null 2>&1; then
      echo "cis output exists; skipping ${prefix}"
    else
      "${PY}" -m tensorqtl "${GENO_PREFIX}" "${bed}" "${prefix}" \
        --covariates "${covars}" \
        --mode cis \
        --maf_threshold 0.05
    fi
  fi

  if [[ "${MODE}" == "cis_nominal" || "${MODE}" == "both" ]]; then
    local nom_dir="${BASE}/tensorqtl/cis_nominal/${scenario}"
    mkdir -p "${nom_dir}"
    local prefix="${nom_dir}/${pset}"
    if ls "${prefix}"*.cis_qtl_pairs.* >/dev/null 2>&1; then
      echo "cis_nominal output exists; skipping ${prefix}"
    else
      "${PY}" -m tensorqtl "${GENO_PREFIX}" "${bed}" "${prefix}" \
        --covariates "${covars}" \
        --mode cis_nominal \
        --maf_threshold 0.05
    fi
  fi
  echo "finish: $(date)"
}

select_jobs() {
  if [[ "${SELECTOR}" == "all" ]]; then
    awk -F'\t' 'NR>1 {print $1}' "${MANIFEST}"
  elif [[ "${SELECTOR}" =~ ^[0-9]+$ ]]; then
    echo "${SELECTOR}"
  else
    awk -F'\t' -v sc="${SELECTOR}" 'NR>1 && $2==sc {print $1}' "${MANIFEST}"
  fi
}

select_jobs | while read -r jid; do
  [[ -z "${jid}" ]] && continue
  run_one "${jid}" 2>&1 | tee "${BASE}/logs/tensorqtl/job_${jid}_${MODE}.log"
done
