#!/usr/bin/env bash
set -euo pipefail

MODE="${1:?Use main, functional, or ngenes}"
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODULE_DIR="${SC_PCQTL_SLDSC_ROOT:-${CODE_DIR}}"
export SC_PCQTL_SLDSC_ROOT="${MODULE_DIR}"
PYTHON="${LDSC_PYTHON:-${MODULE_DIR}/tools/ldsc_py3_venv/bin/python}"
RESULTS="${MODULE_DIR}/results/heritability_enrichment"
TRAIT_DIR="${MODULE_DIR}/manifests/trait_independence"
ANALYSIS_PREFIX="${SLDSC_ANALYSIS_PREFIX:-prespecified247_h2qc}"
cd "${MODULE_DIR}"

collect() {
  "${PYTHON}" "${CODE_DIR}/scripts/63_collect_eur_sldsc_multi_results.py" "$@"
}
meta() {
  "${PYTHON}" "${CODE_DIR}/scripts/64_meta_enrichment_taustar_paperstyle.py" "$@"
}

require_trait_manifests() {
  local slug
  for slug in rg0p6 rg0p7 rg0p8; do
    [[ -s "${TRAIT_DIR}/selected_endpoints_${slug}.tsv" ]] || {
      echo "Missing independent-trait manifest: ${TRAIT_DIR}/selected_endpoints_${slug}.tsv" >&2
      exit 1
    }
  done
}

meta_sets() {
  local collector="$1" out="$2" prefix="$3" label_map="${4:-}"
  local label_args=()
  [[ -n "${label_map}" ]] && label_args=(--label-map "${label_map}")
  meta --collector "${collector}" --out-dir "${out}" --prefix "${prefix}" \
    --trait-manifest "${TRAIT_DIR}/selected_endpoints_rg0p7.tsv" \
    --analysis-label independent_traits_rg0p7 "${label_args[@]}"
  meta --collector "${collector}" --out-dir "${out}" --prefix "${prefix}" \
    --trait-manifest "${TRAIT_DIR}/selected_endpoints_rg0p6.tsv" \
    --analysis-label independent_traits_rg0p6 --output-suffix rg0p6_sensitivity "${label_args[@]}"
  meta --collector "${collector}" --out-dir "${out}" --prefix "${prefix}" \
    --trait-manifest "${TRAIT_DIR}/selected_endpoints_rg0p8.tsv" \
    --analysis-label independent_traits_rg0p8 --output-suffix rg0p8_sensitivity "${label_args[@]}"
}

delta_meta_sets() {
  local out="$1" prefix="$2"
  "${PYTHON}" "${CODE_DIR}/scripts/59_meta_pcqtl_eqtl_tau_star.py" \
    --summary-dir "${out}" --prefix "${prefix}" \
    --trait-manifest "${TRAIT_DIR}/selected_endpoints_rg0p7.tsv" \
    --analysis-label independent_traits_rg0p7
  "${PYTHON}" "${CODE_DIR}/scripts/59_meta_pcqtl_eqtl_tau_star.py" \
    --summary-dir "${out}" --prefix "${prefix}" \
    --trait-manifest "${TRAIT_DIR}/selected_endpoints_rg0p6.tsv" \
    --analysis-label independent_traits_rg0p6 --output-suffix rg0p6_sensitivity
  "${PYTHON}" "${CODE_DIR}/scripts/59_meta_pcqtl_eqtl_tau_star.py" \
    --summary-dir "${out}" --prefix "${prefix}" \
    --trait-manifest "${TRAIT_DIR}/selected_endpoints_rg0p8.tsv" \
    --analysis-label independent_traits_rg0p8 --output-suffix rg0p8_sensitivity
}

require_trait_manifests

case "${MODE}" in
  main)
    joint="${RESULTS}/eur_sldsc_qtlsig_${ANALYSIS_PREFIX}_joint_summary"
    marg_pc="${RESULTS}/eur_sldsc_qtlsig_${ANALYSIS_PREFIX}_marg_pc_summary"
    marg_eq="${RESULTS}/eur_sldsc_qtlsig_${ANALYSIS_PREFIX}_marg_eq_summary"
    collect --manifest "manifests/${ANALYSIS_PREFIX}_joint_tasks.tsv" --names 'QTLsig_pcQTL_gt0;QTLsig_eQTL_gt0' --delta-pair 'QTLsig_pcQTL_gt0,QTLsig_eQTL_gt0' --out-dir "${joint}" --prefix "${ANALYSIS_PREFIX}_joint"
    collect --manifest "manifests/${ANALYSIS_PREFIX}_marg_pc_tasks.tsv" --names 'QTLsig_pcQTL_gt0' --out-dir "${marg_pc}" --prefix "${ANALYSIS_PREFIX}_marg_pc"
    collect --manifest "manifests/${ANALYSIS_PREFIX}_marg_eq_tasks.tsv" --names 'QTLsig_eQTL_gt0' --out-dir "${marg_eq}" --prefix "${ANALYSIS_PREFIX}_marg_eq"
    delta_meta_sets "${joint}" "${ANALYSIS_PREFIX}_joint"
    meta_sets "${joint}/${ANALYSIS_PREFIX}_joint.tsv" "${joint}" "${ANALYSIS_PREFIX}_joint"
    meta_sets "${marg_pc}/${ANALYSIS_PREFIX}_marg_pc.tsv" "${marg_pc}" "${ANALYSIS_PREFIX}_marg_pc"
    meta_sets "${marg_eq}/${ANALYSIS_PREFIX}_marg_eq.tsv" "${marg_eq}" "${ANALYSIS_PREFIX}_marg_eq"
    ;;
  functional)
    out="${RESULTS}/eur_sldsc_qtlsig_${ANALYSIS_PREFIX}_funcinteract_summary"
    names='pcQTL_prom;pcQTL_enh;eQTL_prom;eQTL_enh'
    collect --manifest "manifests/${ANALYSIS_PREFIX}_funcinteract_tasks.tsv" --names "${names}" --delta-pair 'pcQTL_prom,eQTL_prom' --out-dir "${out}" --prefix fi_prom
    collect --manifest "manifests/${ANALYSIS_PREFIX}_funcinteract_tasks.tsv" --names "${names}" --delta-pair 'pcQTL_enh,eQTL_enh' --out-dir "${out}" --prefix fi_enh
    delta_meta_sets "${out}" fi_prom
    delta_meta_sets "${out}" fi_enh
    meta_sets "${out}/fi_prom.tsv" "${out}" fi_funcinteract
    ;;
  ngenes)
    out="${RESULTS}/eur_sldsc_qtlsig_${ANALYSIS_PREFIX}_ngenes_summary"
    collect --manifest "manifests/${ANALYSIS_PREFIX}_ngenes_tasks.tsv" --names 'ng_small;ng_med;ng_large' --delta-pair 'ng_large,ng_small' --out-dir "${out}" --prefix ng
    delta_meta_sets "${out}" ng
    meta_sets "${out}/ng.tsv" "${out}" ng 'ng_small=2genes;ng_med=3-4genes;ng_large=5+genes'
    ;;
  *) echo "Unknown mode: ${MODE}" >&2; exit 2 ;;
esac
