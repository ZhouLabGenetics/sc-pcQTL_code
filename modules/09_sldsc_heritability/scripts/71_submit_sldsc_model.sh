#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

MODE="${1:?Use main, functional, or ngenes}"
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODULE_DIR="${SC_PCQTL_SLDSC_ROOT:-${CODE_DIR}}"
export SC_PCQTL_SLDSC_ROOT="${MODULE_DIR}"
export SC_PCQTL_SLDSC_CODE_DIR="${CODE_DIR}"
PYTHON="${SC_PCQTL_PYTHON:-python3}"
ENDPOINT_MANIFEST="${FINNGEN_SLDSC_ENDPOINT_MANIFEST:-${MODULE_DIR}/manifests/finngen_r12_sldsc_prespecified247_endpoints.tsv}"
ANALYSIS_PREFIX="${SLDSC_ANALYSIS_PREFIX:-prespecified247_h2qc}"
ARRAY_CONCURRENCY="${SLDSC_ARRAY_CONCURRENCY:-24}"
PARTITION="${SLDSC_PARTITION:-normal}"
RG_MATRIX="${FINNGEN_R12_RG_MATRIX:-}"
TRAIT_DIR="${MODULE_DIR}/manifests/trait_independence"
UPSTREAM_DEPENDENCY="${SLDSC_UPSTREAM_DEPENDENCY:-}"
SUMSTATS_DEPENDENCY="${SLDSC_SUMSTATS_DEPENDENCY:-}"
REUSE_LDSCORES="${SLDSC_REUSE_LDSCORES:-0}"
cd "${MODULE_DIR}"
mkdir -p logs manifests

[[ -s "${ENDPOINT_MANIFEST}" ]] || { echo "Missing endpoint manifest: ${ENDPOINT_MANIFEST}" >&2; exit 1; }
if [[ "${PREPARE_TRAIT_SELECTION:-1}" == 1 ]]; then
  [[ -n "${RG_MATRIX}" && -s "${RG_MATRIX}" ]] || {
    echo "Set FINNGEN_R12_RG_MATRIX to a complete official FinnGen R12 FIN rg matrix" >&2
    exit 1
  }
  "${PYTHON}" "${CODE_DIR}/scripts/72_select_independent_finngen_traits.py" \
    --endpoint-manifest "${ENDPOINT_MANIFEST}" --rg-matrix "${RG_MATRIX}" --out-dir "${TRAIT_DIR}" \
    --h2-z-min 4 --rg-thresholds 0.6,0.7,0.8 --primary-rg-threshold 0.7 --rg-fdr 0.05
fi
MODEL_ENDPOINT_MANIFEST="${TRAIT_DIR}/h2_eligible_endpoints.tsv"
[[ -s "${MODEL_ENDPOINT_MANIFEST}" ]] || { echo "Missing h2-QC endpoint manifest: ${MODEL_ENDPOINT_MANIFEST}" >&2; exit 1; }
"${PYTHON}" "${CODE_DIR}/scripts/69_build_sldsc_model_manifests.py" "${MODE}" \
  --endpoint-manifest "${MODEL_ENDPOINT_MANIFEST}" --analysis-prefix "${ANALYSIS_PREFIX}"

sb_retry() {
  local jid
  for _attempt in $(seq 1 "${SBATCH_RETRIES:-12}"); do
    jid=$(sbatch --parsable "$@" 2>/dev/null || true)
    if [[ "${jid}" =~ ^[0-9]+$ ]]; then printf '%s\n' "${jid}"; return 0; fi
    sleep "${SBATCH_RETRY_SECONDS:-15}"
  done
  echo "Slurm submission failed: $*" >&2
  return 1
}

exclude_args=()
[[ -n "${SLDSC_EXCLUDE_NODES:-}" ]] && exclude_args=(--exclude="${SLDSC_EXCLUDE_NODES}")

submit_array() {
  local manifest="$1" dependency="$2"
  local n_tasks
  n_tasks="$(awk 'END{print NR-1}' "${manifest}")"
  [[ "${n_tasks}" -gt 0 ]] || { echo "No tasks in ${manifest}" >&2; return 1; }
  sb_retry --partition="${PARTITION}" --array="1-${n_tasks}%${ARRAY_CONCURRENCY}" --dependency="afterok:${dependency}" \
    "${exclude_args[@]}" --export="ALL,QTLSIG_MULTI_MANIFEST=${manifest}" \
    "${CODE_DIR}/scripts/submit_eur_sldsc_multi.sbatch"
}

join_dependencies() {
  local joined="" value
  for value in "$@"; do
    [[ -n "${value}" ]] || continue
    joined="${joined:+${joined}:}${value}"
  done
  [[ -n "${joined}" ]] || { echo "At least one upstream dependency is required" >&2; return 1; }
  printf '%s\n' "${joined}"
}

case "${MODE}" in
  main)
    prep_dep=()
    if [[ "${PREPARE_INPUTS:-0}" == 1 ]]; then
      shopt -s nullglob
      old_partials=(annotations/sig_cisqtl/partials/*.tsv)
      shopt -u nullglob
      [[ ${#old_partials[@]} -eq 0 ]] || {
        echo "Refusing to mix a fresh extraction with existing partials under ${MODULE_DIR}/annotations/sig_cisqtl/partials" >&2
        exit 1
      }
      pc_extract=$(sb_retry --partition="${PARTITION}" --array="1-10%${QTLSIG_EXTRACT_CONCURRENCY:-10}" \
        --job-name=scpcqtl_extract_pc "${exclude_args[@]}" \
        --export="ALL,QTLSIG_EXTRACT_QTL=pcQTL" "${CODE_DIR}/scripts/submit_qtlsig_extract.sbatch")
      eq_extract=$(sb_retry --partition="${PARTITION}" --array="1-10%${QTLSIG_EXTRACT_CONCURRENCY:-10}" \
        --job-name=scpcqtl_extract_eq "${exclude_args[@]}" \
        --export="ALL,QTLSIG_EXTRACT_QTL=eQTL" "${CODE_DIR}/scripts/submit_qtlsig_extract.sbatch")
      prep=$(sb_retry --partition="${PARTITION}" --dependency="afterok:${pc_extract}:${eq_extract}" \
        "${exclude_args[@]}" "${CODE_DIR}/scripts/submit_qtlsig_prepare.sbatch")
      prep_dep=(--dependency="afterok:${prep}")
    elif [[ "${REUSE_LDSCORES}" != 1 ]]; then
      for file in annotations/sig_cisqtl/sig_pcQTL_snps_maf05.tsv annotations/sig_cisqtl/sig_eQTL_snps_maf05.tsv resources/eur_sldsc_ref/baselineLD_noQTL/baselineLD_noQTL.1.l2.ldscore.gz; do
        [[ -s "${file}" ]] || { echo "Missing ${file}; rerun with PREPARE_INPUTS=1" >&2; exit 1; }
      done
    fi
    if [[ "${REUSE_LDSCORES}" == 1 ]]; then
      [[ -n "${UPSTREAM_DEPENDENCY}" ]] || { echo "Set SLDSC_UPSTREAM_DEPENDENCY when reusing LD scores" >&2; exit 1; }
      l2_qtl="${UPSTREAM_DEPENDENCY}"
    else
      l2_qtl=$(sb_retry --partition="${PARTITION}" "${prep_dep[@]}" "${exclude_args[@]}" "${CODE_DIR}/scripts/submit_eur_sldsc_qtlsig_ldscores.sbatch")
    fi
    dep="$(join_dependencies "${l2_qtl}" "${SUMSTATS_DEPENDENCY}")"
    jobs=(
      "$(submit_array manifests/${ANALYSIS_PREFIX}_joint_tasks.tsv "${dep}")"
      "$(submit_array manifests/${ANALYSIS_PREFIX}_marg_pc_tasks.tsv "${dep}")"
      "$(submit_array manifests/${ANALYSIS_PREFIX}_marg_eq_tasks.tsv "${dep}")"
    )
    ;;
  functional)
    functional_dep=()
    if [[ "${REUSE_LDSCORES}" == 1 ]]; then
      [[ -n "${UPSTREAM_DEPENDENCY}" ]] || { echo "Set SLDSC_UPSTREAM_DEPENDENCY when reusing LD scores" >&2; exit 1; }
      l2="${UPSTREAM_DEPENDENCY}"
    elif [[ -n "${UPSTREAM_DEPENDENCY}" ]]; then
      functional_dep=(--dependency="afterok:${UPSTREAM_DEPENDENCY}")
      l2=$(sb_retry --partition="${PARTITION}" "${functional_dep[@]}" "${exclude_args[@]}" "${CODE_DIR}/scripts/submit_eur_sldsc_funcinteract_ldscores.sbatch")
    else
      [[ -s resources/eur_sldsc_custom_annotations/single_maf05/QTLsig_pcQTL_gt0/qtlsig_pcQTL.1.annot.gz ]] || { echo "Run the main annotation chain first or set SLDSC_UPSTREAM_DEPENDENCY" >&2; exit 1; }
      l2=$(sb_retry --partition="${PARTITION}" "${exclude_args[@]}" "${CODE_DIR}/scripts/submit_eur_sldsc_funcinteract_ldscores.sbatch")
    fi
    dep="$(join_dependencies "${l2}" "${SUMSTATS_DEPENDENCY}")"
    jobs=("$(submit_array manifests/${ANALYSIS_PREFIX}_funcinteract_tasks.tsv "${dep}")")
    ;;
  ngenes)
    if [[ "${REUSE_LDSCORES}" == 1 ]]; then
      [[ -n "${UPSTREAM_DEPENDENCY}" ]] || { echo "Set SLDSC_UPSTREAM_DEPENDENCY when reusing LD scores" >&2; exit 1; }
      l2="${UPSTREAM_DEPENDENCY}"
    elif [[ "${PREPARE_INPUTS:-1}" == 1 ]]; then
      prep=$(sb_retry --partition="${PARTITION}" "${exclude_args[@]}" "${CODE_DIR}/scripts/submit_ngenes_prepare.sbatch")
      l2=$(sb_retry --partition="${PARTITION}" --dependency="afterok:${prep}" "${exclude_args[@]}" "${CODE_DIR}/scripts/submit_eur_sldsc_ngenes_ldscores.sbatch")
    else
      [[ -s annotations/sig_cisqtl/pcQTL_ngenes_small_maf05.tsv ]] || { echo "Missing n-genes SNP sets" >&2; exit 1; }
      l2=$(sb_retry --partition="${PARTITION}" "${exclude_args[@]}" "${CODE_DIR}/scripts/submit_eur_sldsc_ngenes_ldscores.sbatch")
    fi
    dep="$(join_dependencies "${l2}" "${SUMSTATS_DEPENDENCY}")"
    jobs=("$(submit_array manifests/${ANALYSIS_PREFIX}_ngenes_tasks.tsv "${dep}")")
    ;;
  *) echo "Unknown mode: ${MODE}" >&2; exit 2 ;;
esac

job_dep=$(IFS=:; echo "${jobs[*]}")
if [[ "${SUBMIT_COLLECT:-1}" == 1 ]]; then
  for slug in rg0p6 rg0p7 rg0p8; do
    [[ -s "${TRAIT_DIR}/selected_endpoints_${slug}.tsv" ]] || {
      echo "Missing trait selection for collection; rerun with PREPARE_TRAIT_SELECTION=1" >&2
      exit 1
    }
  done
  collect=$(sb_retry --dependency="afterok:${job_dep}" "${exclude_args[@]}" \
    --job-name="scpcqtl_sldsc_${MODE}_collect" --partition="${PARTITION}" --time=00:30:00 --mem=8G \
    --output="logs/${MODE}_collect_%j.out" --error="logs/${MODE}_collect_%j.err" \
    --export="ALL,SLDSC_ANALYSIS_PREFIX=${ANALYSIS_PREFIX}" \
    --wrap="bash ${CODE_DIR}/scripts/70_collect_sldsc_models.sh ${MODE}")
else
  collect="not_submitted"
fi
printf 'mode=%s\nanalysis_jobs=%s\ncollect_job=%s\n' "${MODE}" "${job_dep}" "${collect}"
