#!/usr/bin/env python3
"""Create lightweight supplementary TSV tables from fixed workflow outputs."""

from __future__ import annotations

import csv
import gzip
import math
import os
import sys
from argparse import ArgumentParser
from itertools import zip_longest
from pathlib import Path

RELEASE_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(RELEASE_ROOT / "bin"))
sys.dont_write_bytecode = True
from celltype_eligibility import MIN_CELLTYPE_CELLS, load_eligibility, primary_celltypes

PRIMARY_CELLTYPES = set(primary_celltypes())


WORKFLOW_ROOT = Path(os.environ.get("COQTL_WORKFLOW_ROOT", os.environ.get("SC_PCQTL_WORKFLOW_ROOT", "work/coQTL_workflow")))
RESULT_ROOT = Path(os.environ.get("SC_PCQTL_PRIMARY_RESULT_ROOT", WORKFLOW_ROOT))
UPSTREAM_CLUSTER_ROOT = Path(os.environ.get(
    "COQTL_UPSTREAM_CELLTYPES_DIR",
    WORKFLOW_ROOT / "03_analysis_celltypes/01_upstream_main_pipeline_add_cov/celltypes",
))
SLDSC_RUN_ROOT = Path(os.environ.get(
    "SC_PCQTL_SLDSC_RESULT_ROOT",
    RESULT_ROOT / (
        "04_formal_colocalization/08_heritability_enrichment_maxpip/runs/"
        "primary_10cell_prespecified247_rgpruned"
    ),
))
SLDSC_ROOT = SLDSC_RUN_ROOT / "results/heritability_enrichment"
CROSSMAP_ROOT = Path(os.environ.get(
    "SC_PCQTL_CROSSMAP_ROOT",
    RESULT_ROOT / "05_review_resistance/wo05_cross_mappability",
))
JOINT_SCORE_COMPARISON_ROOT = Path(os.environ.get(
    "SC_PCQTL_JOINT_SCORE_COMPARISON_ROOT",
    RESULT_ROOT / "joint_p_test/comparison/tables",
))
SLDSC_ANNOTATION_ROOT = (
    SLDSC_RUN_ROOT / "resources/eur_sldsc_custom_annotations/single_maf05"
)
SLDSC_FRQ_ROOT = SLDSC_RUN_ROOT / "resources/eur_sldsc_ref/1000G_Phase3_frq"

SOURCES = {
    "pcqtl_summary": RESULT_ROOT
    / "03_analysis_celltypes/02_downstream_analysis_modules_add_cov_fdr/pcqtl_compare_saigeqtl/data/summary_by_celltype.tsv",
    "cluster_counts": RESULT_ROOT
    / "03_analysis_celltypes/04_cluster_annotation_enrichment_add_cov/results/cluster_sets/add_cov_cluster_count_summary.tsv",
    "qtl_counts": RESULT_ROOT
    / "03_analysis_celltypes/02_downstream_analysis_modules_add_cov_fdr/pcqtl_compare/data/qtl_counts_per_celltype.tsv",
    "pcqtl_all_results": RESULT_ROOT
    / "03_analysis_celltypes/02_downstream_analysis_modules_add_cov_fdr/pcqtl_compare/data/all_pcqtl_results.tsv",
    "enrichment": RESULT_ROOT
    / "03_analysis_celltypes/04_cluster_annotation_enrichment_add_cov/results/enrichment/enrichment_by_method.tsv",
    "strict_susie_edges": RESULT_ROOT
    / "04_formal_colocalization/06_strict_signal_grouping/results/strict_graph/strict_graph_edges.tsv",
    "strict_susie_groups": RESULT_ROOT
    / "04_formal_colocalization/06_strict_signal_grouping/results/strict_graph/strict_signal_groups.tsv",
    "strict_pcqtl_mechanism": RESULT_ROOT
    / (
        "04_formal_colocalization/09_mechanistic_celltype_analysis/results/regulatory_annotation/"
        "strict_pcqtl_specific_mechanism_table_with_regulatory_annotation.tsv"
    ),
    "strict_susie_sensitivity": RESULT_ROOT
    / "04_formal_colocalization/06_strict_signal_grouping/results/strict_graph/strict_additional_hit_summary.tsv",
    "crossmap_summary": CROSSMAP_ROOT / "results/tables/crossmap_sensitivity_summary.tsv",
    "crossmap_gimap_pairs": CROSSMAP_ROOT / "results/gimap_case/gimap_crossmap_pair_audit.tsv",
    "joint_score_cluster_sensitivity": JOINT_SCORE_COMPARISON_ROOT
    / "table_s2_joint_score_cluster_sensitivity.tsv",
    "sldsc_marg_pc_meta": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_marg_pc_summary/prespecified247_h2qc_marg_pc_annotation_meta.tsv",
    "sldsc_marg_eq_meta": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_marg_eq_summary/prespecified247_h2qc_marg_eq_annotation_meta.tsv",
    "sldsc_joint_per_endpoint": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_joint_summary/prespecified247_h2qc_joint.tsv",
    "sldsc_joint_meta": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_joint_summary/prespecified247_h2qc_joint_annotation_meta.tsv",
    "sldsc_joint_delta": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_joint_summary/prespecified247_h2qc_joint_delta.tsv",
    "sldsc_joint_delta_meta": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_joint_summary/prespecified247_h2qc_joint_delta_meta.tsv",
    "sldsc_joint_meta_rg0p6": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_joint_summary/prespecified247_h2qc_joint_annotation_meta_rg0p6_sensitivity.tsv",
    "sldsc_joint_meta_rg0p8": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_joint_summary/prespecified247_h2qc_joint_annotation_meta_rg0p8_sensitivity.tsv",
    "sldsc_joint_delta_meta_rg0p6": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_joint_summary/prespecified247_h2qc_joint_delta_meta_rg0p6_sensitivity.tsv",
    "sldsc_joint_delta_meta_rg0p8": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_joint_summary/prespecified247_h2qc_joint_delta_meta_rg0p8_sensitivity.tsv",
    "sldsc_functional_meta": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_funcinteract_summary/fi_funcinteract_annotation_meta.tsv",
    "sldsc_promoter_delta": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_funcinteract_summary/fi_prom_delta.tsv",
    "sldsc_promoter_delta_meta": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_funcinteract_summary/fi_prom_delta_meta.tsv",
    "sldsc_enhancer_delta": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_funcinteract_summary/fi_enh_delta.tsv",
    "sldsc_enhancer_delta_meta": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_funcinteract_summary/fi_enh_delta_meta.tsv",
    "sldsc_ngenes_per_endpoint": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_ngenes_summary/ng.tsv",
    "sldsc_ngenes_meta": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_ngenes_summary/ng_annotation_meta.tsv",
    "sldsc_ngenes_delta": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_ngenes_summary/ng_delta.tsv",
    "sldsc_ngenes_delta_meta": SLDSC_ROOT
    / "eur_sldsc_qtlsig_prespecified247_h2qc_ngenes_summary/ng_delta_meta.tsv",
    "sldsc_primary_traits": SLDSC_RUN_ROOT
    / "manifests/trait_independence/selected_endpoints_rg0p7.tsv",
    "sldsc_trait_selection_rg0p7": SLDSC_RUN_ROOT
    / "manifests/trait_independence/trait_selection_rg0p7.tsv",
    "sldsc_trait_selection_summary": SLDSC_RUN_ROOT
    / "manifests/trait_independence/trait_selection_summary.tsv",
    "sldsc_manifest_provenance": SLDSC_RUN_ROOT
    / "manifests/prespecified247_manifest_provenance.tsv",
    "sldsc_snpset_qc": SLDSC_RUN_ROOT
    / "annotations/sig_cisqtl/sig_cisqtl_snpset_qc.tsv",
    "sldsc_pcqtl_snps": SLDSC_RUN_ROOT
    / "annotations/sig_cisqtl/sig_pcQTL_snps_maf05.tsv",
    "sldsc_eqtl_snps": SLDSC_RUN_ROOT
    / "annotations/sig_cisqtl/sig_eQTL_snps_maf05.tsv",
}
FORMAL_TABLES = {
    "s1": "table_s1_celltype_eligibility_qc.tsv",
    "s2": "table_s2_joint_score_cluster_sensitivity.tsv",
    "s3": "table_s3_gimap_cross_mappability_pair_audit.tsv",
    "s4": "table_s4_acat_sensitivity.tsv",
    "s5": "table_s5_sldsc_primary_summary.tsv",
}

DATA_FILES = {
    "celltype_cluster_pcqtl_summary": "celltype_cluster_pcqtl_summary.tsv",
    "cluster_annotation_enrichment": "cluster_annotation_enrichment.tsv",
    "susie_strict_graph_edges": "susie_strict_graph_edges.tsv",
    "susie_strict_signal_groups": "susie_strict_signal_groups.tsv",
    "strict_pcqtl_specific_loading_effects": "strict_pcqtl_specific_loading_effects.tsv",
    "susie_strict_threshold_sensitivity": "susie_strict_threshold_sensitivity.tsv",
    "cross_mappability_sensitivity": "cross_mappability_sensitivity.tsv",
}

SLDSC_TABLES = {
    "sldsc_marg_pc_meta": "prespecified247_marg_pc_annotation_meta.tsv",
    "sldsc_marg_eq_meta": "prespecified247_marg_eq_annotation_meta.tsv",
    "sldsc_joint_per_endpoint": "prespecified247_joint_per_endpoint.tsv",
    "sldsc_joint_meta": "prespecified247_joint_annotation_meta.tsv",
    "sldsc_joint_delta": "prespecified247_joint_delta.tsv",
    "sldsc_joint_delta_meta": "prespecified247_joint_delta_meta.tsv",
    "sldsc_joint_meta_rg0p6": "prespecified247_joint_annotation_meta_rg0p6_sensitivity.tsv",
    "sldsc_joint_meta_rg0p8": "prespecified247_joint_annotation_meta_rg0p8_sensitivity.tsv",
    "sldsc_joint_delta_meta_rg0p6": "prespecified247_joint_delta_meta_rg0p6_sensitivity.tsv",
    "sldsc_joint_delta_meta_rg0p8": "prespecified247_joint_delta_meta_rg0p8_sensitivity.tsv",
    "sldsc_functional_meta": "prespecified247_promoter_enhancer_annotation_meta.tsv",
    "sldsc_promoter_delta": "prespecified247_promoter_delta_per_endpoint.tsv",
    "sldsc_promoter_delta_meta": "prespecified247_promoter_delta_meta.tsv",
    "sldsc_enhancer_delta": "prespecified247_enhancer_delta_per_endpoint.tsv",
    "sldsc_enhancer_delta_meta": "prespecified247_enhancer_delta_meta.tsv",
    "sldsc_ngenes_per_endpoint": "prespecified247_cluster_size_per_endpoint.tsv",
    "sldsc_ngenes_meta": "prespecified247_cluster_size_annotation_meta.tsv",
    "sldsc_ngenes_delta": "prespecified247_cluster_size_delta_per_endpoint.tsv",
    "sldsc_ngenes_delta_meta": "prespecified247_cluster_size_delta_meta.tsv",
    "sldsc_trait_selection_rg0p7": "prespecified247_trait_selection_rg0p7.tsv",
    "sldsc_trait_selection_summary": "prespecified247_trait_selection_summary.tsv",
    "sldsc_manifest_provenance": "prespecified247_manifest_provenance.tsv",
}

SLDSC_DERIVED_TABLES = {
    "prespecified247_annotation_overlap_summary.tsv",
}

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def require_sources() -> None:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    cluster_files = list(UPSTREAM_CLUSTER_ROOT.glob("*/pcQTL/step2_pca/all_clusters_summary.tsv"))
    if not cluster_files:
        missing.append(str(UPSTREAM_CLUSTER_ROOT / "*/pcQTL/step2_pca/all_clusters_summary.tsv"))
    if missing:
        raise FileNotFoundError("Missing workflow source table(s):\n" + "\n".join(missing))
    validate_final_sources()


def require_sldsc_sources() -> None:
    keys = set(SLDSC_TABLES) | {
        "sldsc_joint_delta",
        "sldsc_trait_selection_rg0p7",
        "sldsc_trait_selection_summary",
        "sldsc_snpset_qc",
        "sldsc_pcqtl_snps",
        "sldsc_eqtl_snps",
    }
    required_reference_files = []
    for chromosome in range(1, 23):
        required_reference_files.extend([
            SLDSC_ANNOTATION_ROOT / "QTLsig_pcQTL_gt0" / f"qtlsig_pcQTL.{chromosome}.annot.gz",
            SLDSC_ANNOTATION_ROOT / "QTLsig_pcQTL_gt0" / f"qtlsig_pcQTL.{chromosome}.l2.M_5_50",
            SLDSC_ANNOTATION_ROOT / "QTLsig_eQTL_gt0" / f"qtlsig_eQTL.{chromosome}.annot.gz",
            SLDSC_ANNOTATION_ROOT / "QTLsig_eQTL_gt0" / f"qtlsig_eQTL.{chromosome}.l2.M_5_50",
            SLDSC_FRQ_ROOT / f"1000G.EUR.QC.{chromosome}.frq",
        ])
    missing = [str(SOURCES[key]) for key in keys if not SOURCES[key].is_file()]
    missing.extend(str(path) for path in required_reference_files if not path.is_file())
    if missing:
        raise FileNotFoundError("Missing S-LDSC source table(s):\n" + "\n".join(missing))
    validate_sldsc_sources()


def assert_primary_celltypes_only(source_key: str, column: str = "celltype") -> None:
    rows = read_tsv(SOURCES[source_key])
    observed: set[str] = set()
    for row in rows:
        observed.update(value for value in row.get(column, "").split(",") if value)
    excluded = observed - PRIMARY_CELLTYPES
    if excluded:
        raise ValueError(
            f"{SOURCES[source_key]} retains ineligible cell types {sorted(excluded)}. "
            "Rerun the producing workflow module with config/celltype_eligibility.tsv; "
            "post hoc row filtering is unsafe because strict signal-group IDs are rebuilt."
        )


def require_int(row: dict[str, str], field: str, expected: int, label: str) -> None:
    observed = int(float(row[field]))
    if observed != expected:
        raise ValueError(f"Stale {label}: observed {field}={observed}, expected {expected}")


def validate_final_sources() -> None:
    """Refuse to publish tables from pre-eligibility aggregate outputs."""
    for source_key in (
        "pcqtl_summary",
        "qtl_counts",
        "pcqtl_all_results",
        "strict_susie_edges",
        "strict_pcqtl_mechanism",
    ):
        assert_primary_celltypes_only(source_key)
    assert_primary_celltypes_only("strict_susie_groups", column="celltypes")

    enrichment_rows = read_tsv(SOURCES["enrichment"])
    enrichment_key = next(
        (
            row for row in enrichment_rows
            if row.get("method") == "add_cov_sc_hurdle"
            and row.get("stratum") == "all_correlated"
            and row.get("annotation") == "has_shared_go_bp"
        ),
        None,
    )
    if enrichment_key is None:
        raise ValueError("Final add-covariate sc-pcQTL enrichment row is missing")
    require_int(enrichment_key, "n_correlated", 2464, "cluster-enrichment source")
    require_int(enrichment_key, "n_null", 87714, "cluster-enrichment source")

    sensitivity = {row["threshold"]: row for row in read_tsv(SOURCES["strict_susie_sensitivity"])}
    main = sensitivity.get("0.75")
    if main is None:
        raise ValueError("Strict SuSiE threshold 0.75 summary is missing")
    require_int(main, "n_gwas_components", 394, "strict SuSiE summary")
    require_int(main, "n_pcQTL_specific_events", 46, "strict SuSiE summary")

    crossmap_rows = read_tsv(SOURCES["crossmap_summary"])
    crossmap = {(row["section"], row["metric"]): row["value"] for row in crossmap_rows}
    expected_crossmap = {
        ("resource", "final_clusters"): "2485",
        ("pcqtl_counts", "all_clusters::significant_pc_phenotypes"): "2040",
        ("pcqtl_counts", "exclude_known_flagged::significant_pc_phenotypes"): "1778",
        ("strict_susie_h4_0.75", "all_groups::signal_groups"): "394",
        ("strict_susie_h4_0.75", "all_groups::pcqtl_specific"): "46",
    }
    for key, expected in expected_crossmap.items():
        if crossmap.get(key) != expected:
            raise ValueError(
                f"Stale cross-mappability summary for {key}: "
                f"observed={crossmap.get(key)!r}, expected={expected!r}"
            )
    crossmap_p = float(crossmap.get(
        ("independent_t_test", "two_sided_p_value"),
        "nan",
    ))
    if not math.isclose(crossmap_p, 0.158515104850428, rel_tol=0, abs_tol=1e-12):
        raise ValueError(
            "Stale cross-mappability credible-set independent-test p-value: "
            f"observed={crossmap_p!r}"
        )
    gimap_pairs = read_tsv(SOURCES["crossmap_gimap_pairs"])
    if len(gimap_pairs) != 21 or any(row.get("crossmap_score") != "0" for row in gimap_pairs):
        raise ValueError("Final GIMAP cross-mappability table must contain 21 zero-score pairs")

    joint_rows = read_tsv(SOURCES["joint_score_cluster_sensitivity"])
    primary_rows = [row for row in joint_rows if row.get("celltype") != "ALL"]
    if {row.get("celltype") for row in primary_rows} != PRIMARY_CELLTYPES:
        raise ValueError("Joint-score sensitivity table must contain the 10 primary cell types")
    aggregate = next((row for row in joint_rows if row.get("celltype") == "ALL"), None)
    if aggregate is None:
        raise ValueError("Joint-score sensitivity table is missing its aggregate row")
    require_int(aggregate, "primary_pairs", 64295, "joint-score sensitivity source")
    require_int(aggregate, "joint_score_pairs", 121821, "joint-score sensitivity source")
    require_int(aggregate, "shared_pairs", 58851, "joint-score sensitivity source")
    require_int(aggregate, "primary_clusters", 2485, "joint-score sensitivity source")
    require_int(aggregate, "joint_score_clusters", 3036, "joint-score sensitivity source")
    require_int(aggregate, "exact_clusters", 1956, "joint-score sensitivity source")

    validate_sldsc_sources()


def validate_sldsc_sources() -> None:
    """Reject pre-prespecified-trait S-LDSC summaries."""

    marg_pc = read_tsv(SOURCES["sldsc_marg_pc_meta"])
    if len(marg_pc) != 1:
        raise ValueError("Expected one marginal pcQTL S-LDSC meta-analysis row")
    require_int(marg_pc[0], "k_traits", 35, "S-LDSC source")
    observed_enrichment = float(marg_pc[0]["enrichment_re"])
    expected_enrichment = 1.8040765240091787
    if not math.isclose(observed_enrichment, expected_enrichment, rel_tol=0, abs_tol=1e-12):
        raise ValueError(
            "Stale S-LDSC source: pcQTL marginal enrichment="
            f"{observed_enrichment}, expected {expected_enrichment}"
        )

    selection_rows = read_tsv(SOURCES["sldsc_trait_selection_summary"])
    primary = next((row for row in selection_rows if row.get("analysis") == "primary"), None)
    if primary is None:
        raise ValueError("Primary rg=0.7 FinnGen trait-selection summary is missing")
    require_int(primary, "n_manifest_endpoints", 247, "S-LDSC trait selection")
    require_int(primary, "n_h2_eligible", 91, "S-LDSC trait selection")
    require_int(primary, "n_selected_for_meta", 35, "S-LDSC trait selection")


def make_celltype_cluster_pcqtl_data(outdir: Path) -> None:
    pcqtl_rows = {
        row["celltype"]: row
        for row in read_tsv(SOURCES["pcqtl_summary"])
        if row["celltype"] in PRIMARY_CELLTYPES
    }
    # Total cluster-PC phenotypes with valid SAIGE-QTL output per cell type.
    # An included cell type absent from this file had no successful test.
    qtl_test_counts = {
        row["celltype"]: row["n_tests"]
        for row in read_tsv(SOURCES["qtl_counts"])
        if row["celltype"] in PRIMARY_CELLTYPES
    }

    all_cluster_totals: dict[str, int] = {}
    all_cluster_size_counts: dict[str, dict[int, int]] = {}
    selected_pc_counts: dict[str, int] = {}
    for path in sorted(UPSTREAM_CLUSTER_ROOT.glob("*/pcQTL/step2_pca/all_clusters_summary.tsv")):
        celltype = path.parts[-4]
        if celltype not in PRIMARY_CELLTYPES:
            continue
        all_cluster_size_counts[celltype] = {}
        for row in read_tsv(path):
            if row.get("status") and row.get("status") != "SUCCESS":
                continue
            n_genes = int(float(row["n_genes"]))
            all_cluster_totals[celltype] = all_cluster_totals.get(celltype, 0) + 1
            all_cluster_size_counts[celltype][n_genes] = all_cluster_size_counts[celltype].get(n_genes, 0) + 1
            selected_pc_counts[celltype] = selected_pc_counts.get(celltype, 0) + int(float(row["n_pcs_95pct"]))

    eligible_totals: dict[str, int] = {}
    eligible_size_parts: dict[str, list[str]] = {}
    for row in read_tsv(SOURCES["cluster_counts"]):
        if row.get("method") != "add_cov_sc_hurdle":
            continue
        celltype = row["celltype"]
        if celltype not in PRIMARY_CELLTYPES:
            continue
        n = int(float(row["N"]))
        num_genes = row["num_genes"]
        eligible_totals[celltype] = eligible_totals.get(celltype, 0) + n
        eligible_size_parts.setdefault(celltype, []).append(f"{num_genes}:{n}")

    rows: list[dict[str, object]] = []
    for celltype in sorted(all_cluster_totals):
        pc = pcqtl_rows.get(celltype, {})
        n_selected = selected_pc_counts.get(celltype, 0)
        n_successfully_tested = int(qtl_test_counts.get(celltype, "0"))
        if n_successfully_tested > n_selected:
            raise ValueError(f"Successfully tested count exceeds selected count for {celltype}")
        all_size_parts = [
            f"{num_genes}:{count}"
            for num_genes, count in sorted(all_cluster_size_counts.get(celltype, {}).items())
        ]
        rows.append(
            {
                "celltype": celltype,
                "n_sc_pcqtl_clusters_total": all_cluster_totals[celltype],
                "cluster_size_count_by_num_genes_all": ";".join(all_size_parts),
                "n_enrichment_eligible_clusters_size_2_to_5": eligible_totals.get(celltype, 0),
                "enrichment_eligible_size_count_by_num_genes": ";".join(eligible_size_parts.get(celltype, [])),
                "n_selected_cluster_pc_phenotypes": n_selected,
                "n_successfully_tested_cluster_pc_phenotypes": n_successfully_tested,
                "n_failed_convergence_cluster_pc_phenotypes": n_selected - n_successfully_tested,
                "n_within_phenotype_fdr_hit_cluster_pc_phenotypes": pc.get("total_clusters", "0"),
                "n_egene_negative_pcqtl_phenotypes": pc.get("novel", "0"),
                "n_partial_egene_overlap_pcqtl_phenotypes": pc.get("partially_novel", "0"),
                "n_mostly_known_pcqtl_phenotypes": pc.get("mostly_known", "0"),
                "n_all_known_pcqtl_phenotypes": pc.get("all_known", "0"),
                "pct_egene_negative_pcqtl_phenotypes": pc.get("pct_novel", "0"),
                "pct_partial_egene_overlap_pcqtl_phenotypes": pc.get("pct_partially_novel", "0"),
            }
        )

    write_tsv(
        outdir / DATA_FILES["celltype_cluster_pcqtl_summary"],
        rows,
        [
            "celltype",
            "n_sc_pcqtl_clusters_total",
            "cluster_size_count_by_num_genes_all",
            "n_enrichment_eligible_clusters_size_2_to_5",
            "enrichment_eligible_size_count_by_num_genes",
            "n_selected_cluster_pc_phenotypes",
            "n_successfully_tested_cluster_pc_phenotypes",
            "n_failed_convergence_cluster_pc_phenotypes",
            "n_within_phenotype_fdr_hit_cluster_pc_phenotypes",
            "n_egene_negative_pcqtl_phenotypes",
            "n_partial_egene_overlap_pcqtl_phenotypes",
            "n_mostly_known_pcqtl_phenotypes",
            "n_all_known_pcqtl_phenotypes",
            "pct_egene_negative_pcqtl_phenotypes",
            "pct_partial_egene_overlap_pcqtl_phenotypes",
        ],
    )


def make_table_s1(outdir: Path) -> None:
    rows = []
    for row in load_eligibility():
        include = row["include_primary"].upper() == "TRUE"
        rows.append(
            {
                "celltype": row["celltype"],
                "onek1k_label": row["eqtl_celltype"],
                "display_label": row["display_label"],
                "n_cells": int(row["n_cells"]),
                "minimum_cells_for_primary_analysis": MIN_CELLTYPE_CELLS,
                "primary_analysis_status": "included" if include else "excluded",
                "exclusion_reason": "NA" if include else row.get("exclusion_reason", ""),
            }
        )
    write_tsv(
        outdir / FORMAL_TABLES["s1"],
        rows,
        [
            "celltype",
            "onek1k_label",
            "display_label",
            "n_cells",
            "minimum_cells_for_primary_analysis",
            "primary_analysis_status",
            "exclusion_reason",
        ],
    )


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    """Return Benjamini-Hochberg adjusted values in the original row order."""
    n = len(pvalues)
    order = sorted(range(n), key=pvalues.__getitem__)
    adjusted = [1.0] * n
    running_min = 1.0
    for rank_index in range(n - 1, -1, -1):
        row_index = order[rank_index]
        rank = rank_index + 1
        running_min = min(running_min, pvalues[row_index] * n / rank)
        adjusted[row_index] = min(1.0, running_min)
    return adjusted


def make_table_s4(outdir: Path) -> None:
    """Compare the primary within-phenotype variant FDR rule with ACAT-BH."""
    source_rows = [
        row
        for row in read_tsv(SOURCES["pcqtl_all_results"])
        if row["celltype"] in PRIMARY_CELLTYPES
    ]
    if len(source_rows) != 4353:
        raise ValueError(
            "Final ACAT sensitivity source must contain 4,353 successfully "
            f"tested phenotypes, found {len(source_rows)}"
        )

    acat_pvalues: list[float] = []
    for row in source_rows:
        try:
            acat_pvalue = float(row["ACAT_p"])
            min_snp_fdr = float(row["min_snp_fdr"])
        except (KeyError, ValueError) as exc:
            raise ValueError(
                "ACAT sensitivity source requires finite ACAT_p and min_snp_fdr values"
            ) from exc
        if not (math.isfinite(acat_pvalue) and math.isfinite(min_snp_fdr)):
            raise ValueError("ACAT sensitivity source contains non-finite test values")
        acat_pvalues.append(acat_pvalue)

    acat_qvalues = benjamini_hochberg(acat_pvalues)
    summaries = {
        celltype: {
            "celltype": celltype,
            "n_successfully_tested": 0,
            "n_primary_within_phenotype_fdr": 0,
            "n_acat_bh": 0,
            "n_both": 0,
            "n_primary_only": 0,
            "n_acat_only": 0,
        }
        for celltype in sorted(PRIMARY_CELLTYPES)
    }

    for row, acat_qvalue in zip(source_rows, acat_qvalues):
        summary = summaries[row["celltype"]]
        primary_hit = float(row["min_snp_fdr"]) < 0.05
        acat_hit = acat_qvalue < 0.05
        summary["n_successfully_tested"] += 1
        summary["n_primary_within_phenotype_fdr"] += int(primary_hit)
        summary["n_acat_bh"] += int(acat_hit)
        summary["n_both"] += int(primary_hit and acat_hit)
        summary["n_primary_only"] += int(primary_hit and not acat_hit)
        summary["n_acat_only"] += int(acat_hit and not primary_hit)

    rows = list(summaries.values())
    total = {"celltype": "All"}
    count_fields = [
        "n_successfully_tested",
        "n_primary_within_phenotype_fdr",
        "n_acat_bh",
        "n_both",
        "n_primary_only",
        "n_acat_only",
    ]
    for field in count_fields:
        total[field] = sum(int(row[field]) for row in rows)
    rows.append(total)

    expected_total = {
        "n_successfully_tested": 4353,
        "n_primary_within_phenotype_fdr": 2040,
        "n_acat_bh": 2017,
        "n_both": 1992,
        "n_primary_only": 48,
        "n_acat_only": 25,
    }
    for field, expected in expected_total.items():
        if total[field] != expected:
            raise ValueError(
                f"Stale ACAT sensitivity result: {field}={total[field]}, expected {expected}"
            )

    write_tsv(
        outdir / FORMAL_TABLES["s4"],
        rows,
        ["celltype", *count_fields],
    )


def make_table_s5(outdir: Path) -> None:
    """Create the displayed primary S-LDSC summary from final meta-analysis outputs."""
    marginal = {
        "pcQTL": read_tsv(SOURCES["sldsc_marg_pc_meta"])[0],
        "eQTL": read_tsv(SOURCES["sldsc_marg_eq_meta"])[0],
    }
    joint = {
        row["annotation"]: row
        for row in read_tsv(SOURCES["sldsc_joint_meta"])
    }
    contrast = next(
        (
            row
            for row in read_tsv(SOURCES["sldsc_joint_delta_meta"])
            if row.get("model") == "random_effect_DL"
        ),
        None,
    )
    if set(joint) != {"pcQTL", "eQTL"} or contrast is None:
        raise ValueError("Incomplete primary S-LDSC summary sources")

    rows: list[dict[str, object]] = []
    for annotation in ("pcQTL", "eQTL"):
        row = marginal[annotation]
        rows.append(
            {
                "model": "marginal",
                "annotation_or_contrast": annotation,
                "metric": "enrichment",
                "estimate": row["enrichment_re"],
                "ci95_lo": row["enrichment_ci_lo"],
                "ci95_hi": row["enrichment_ci_hi"],
                "pvalue": "NA",
            }
        )
    for annotation in ("pcQTL", "eQTL"):
        row = joint[annotation]
        rows.append(
            {
                "model": "joint",
                "annotation_or_contrast": annotation,
                "metric": "tau_star",
                "estimate": row["tau_star_re"],
                "ci95_lo": row["tau_star_ci_lo"],
                "ci95_hi": row["tau_star_ci_hi"],
                "pvalue": row["tau_star_p"],
            }
        )
    rows.append(
        {
            "model": "joint_contrast",
            "annotation_or_contrast": "pcQTL_minus_eQTL",
            "metric": "delta_tau_star",
            "estimate": contrast["delta_tau_star"],
            "ci95_lo": contrast["ci95_lo"],
            "ci95_hi": contrast["ci95_hi"],
            "pvalue": contrast["p"],
        }
    )
    write_tsv(
        outdir / FORMAL_TABLES["s5"],
        rows,
        [
            "model",
            "annotation_or_contrast",
            "metric",
            "estimate",
            "ci95_lo",
            "ci95_hi",
            "pvalue",
        ],
    )


def copy_table(source: Path, target: Path) -> None:
    if source.resolve() == target.resolve():
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    text = source.read_text(encoding="utf-8", errors="replace")
    # Preserve trailing tabs because they encode empty terminal TSV fields.
    lines = [line.rstrip(" \r") for line in text.splitlines()]
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_table_with_explicit_missing(source: Path, target: Path) -> None:
    """Copy a TSV while representing empty fields explicitly as NA."""
    with source.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in {source}")
        rows = [
            {field: value if value not in ("", None) else "NA" for field, value in row.items()}
            for row in reader
        ]
    write_tsv(target, rows, list(reader.fieldnames))


def copy_sldsc_table_for_submission(source_key: str, source: Path, target: Path) -> None:
    """Remove execution-environment paths from public S-LDSC tables."""
    with source.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in {source}")
        fieldnames = list(reader.fieldnames)
        rows = list(reader)

    if source_key == "sldsc_trait_selection_rg0p7":
        fieldnames.remove("summary_stats_path")
        for row in rows:
            row.pop("summary_stats_path", None)
    elif source_key == "sldsc_manifest_provenance":
        fieldnames = [
            "source_manifest_file" if field == "source_manifest" else field
            for field in fieldnames
        ]
        for row in rows:
            source_manifest = row.pop("source_manifest", "")
            row["source_manifest_file"] = Path(source_manifest).name

    write_tsv(target, rows, fieldnames)


def copy_crossmap_table_for_submission(source: Path, target: Path) -> None:
    """Replace internal cross-mappability status labels with explicit terms."""
    rows = read_tsv(source)
    replacements = (
        ("exclude_known_flagged", "exclude_cross_mappable"),
        ("complete_clean", "no_cross_mappable_pair"),
        ("flagged_clusters_gt100", "cross_mappable_clusters_gt100"),
        ("flagged_minus_clean", "cross_mappable_minus_no_cross_mappable_pair"),
        ("fraction_flagged", "fraction_cross_mappable"),
        ("fraction_clean", "fraction_no_cross_mappable_pair"),
        ("flagged", "cross_mappable"),
    )
    for row in rows:
        metric = row["metric"]
        for old, new in replacements:
            metric = metric.replace(old, new)
        row["metric"] = metric
    write_tsv(target, rows, ["section", "metric", "value"])


def make_sldsc_annotation_overlap_table(target: Path) -> int:
    """Summarize source and EUR-reference pcQTL/eQTL annotation overlap."""
    def variant_key(row: dict[str, str]) -> tuple[str, str, str, str]:
        left, right = sorted((row["a1"], row["a2"]))
        return row["chr"], row["pos"], left, right

    with SOURCES["sldsc_pcqtl_snps"].open(newline="", encoding="utf-8") as handle:
        pc_keys = {variant_key(row) for row in csv.DictReader(handle, delimiter="\t")}
    eq_total = 0
    shared = 0
    with SOURCES["sldsc_eqtl_snps"].open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            eq_total += 1
            shared += variant_key(row) in pc_keys

    qc = {row["qtl"]: row for row in read_tsv(SOURCES["sldsc_snpset_qc"])}
    expected_pc = int(qc["pcQTL"]["n_sig_snps_union_maf_ge_0.05"])
    expected_eq = int(qc["eQTL"]["n_sig_snps_union_maf_ge_0.05"])
    if len(pc_keys) != expected_pc or eq_total != expected_eq:
        raise ValueError(
            "S-LDSC annotation SNP counts disagree with sig_cisqtl_snpset_qc.tsv: "
            f"pcQTL={len(pc_keys)}/{expected_pc}, eQTL={eq_total}/{expected_eq}"
        )

    reference_counts = {
        "pcQTL": 0,
        "eQTL": 0,
        "shared": 0,
        "pcQTL_common": 0,
        "eQTL_common": 0,
        "shared_common": 0,
    }
    for chromosome in range(1, 23):
        pc_path = (
            SLDSC_ANNOTATION_ROOT / "QTLsig_pcQTL_gt0"
            / f"qtlsig_pcQTL.{chromosome}.annot.gz"
        )
        eq_path = (
            SLDSC_ANNOTATION_ROOT / "QTLsig_eQTL_gt0"
            / f"qtlsig_eQTL.{chromosome}.annot.gz"
        )
        frq_path = SLDSC_FRQ_ROOT / f"1000G.EUR.QC.{chromosome}.frq"

        with gzip.open(pc_path, "rt", newline="") as pc_handle, gzip.open(
            eq_path, "rt", newline=""
        ) as eq_handle, frq_path.open(encoding="utf-8") as frq_handle:
            pc_header = pc_handle.readline().rstrip("\r\n").split("\t")
            eq_header = eq_handle.readline().rstrip("\r\n").split("\t")
            pc_snp_index = pc_header.index("SNP")
            eq_snp_index = eq_header.index("SNP")
            pc_value_index = pc_header.index("QTLsig_pcQTL_gt0")
            eq_value_index = eq_header.index("QTLsig_eQTL_gt0")
            frq_header = frq_handle.readline().split()
            snp_index = frq_header.index("SNP")
            maf_index = frq_header.index("MAF")
            for pc_line, eq_line, frq_line in zip_longest(pc_handle, eq_handle, frq_handle):
                if pc_line is None or eq_line is None or frq_line is None:
                    raise ValueError(f"Reference annotation row-count mismatch on chromosome {chromosome}")
                pc_fields = pc_line.rstrip("\r\n").split("\t")
                eq_fields = eq_line.rstrip("\r\n").split("\t")
                frq_fields = frq_line.split()
                pc_snp = pc_fields[pc_snp_index]
                eq_snp = eq_fields[eq_snp_index]
                snp = frq_fields[snp_index]
                if not (pc_snp == eq_snp == snp):
                    raise ValueError(
                        f"Reference SNP-order mismatch on chromosome {chromosome}: "
                        f"{pc_snp}, {eq_snp}, {snp}"
                    )
                maf = float(frq_fields[maf_index])
                pc_value = int(float(pc_fields[pc_value_index]))
                eq_value = int(float(eq_fields[eq_value_index]))
                reference_counts["pcQTL"] += pc_value
                reference_counts["eQTL"] += eq_value
                reference_counts["shared"] += pc_value & eq_value
                if 0.05 <= maf <= 0.5:
                    reference_counts["pcQTL_common"] += pc_value
                    reference_counts["eQTL_common"] += eq_value
                    reference_counts["shared_common"] += pc_value & eq_value

    for annotation, stem in (("pcQTL", "qtlsig_pcQTL"), ("eQTL", "qtlsig_eQTL")):
        m_5_50 = 0.0
        annotation_dir = SLDSC_ANNOTATION_ROOT / f"QTLsig_{annotation}_gt0"
        for chromosome in range(1, 23):
            values = (
                annotation_dir / f"{stem}.{chromosome}.l2.M_5_50"
            ).read_text(encoding="utf-8").split()
            m_5_50 += float(values[0])
        expected_common = reference_counts[f"{annotation}_common"]
        if not math.isclose(m_5_50, expected_common, rel_tol=0, abs_tol=1e-6):
            raise ValueError(
                f"{annotation} reference-common count disagrees with LDSC M_5_50: "
                f"{expected_common} versus {m_5_50}"
            )

    rows = [
        {
            "annotation": "pcQTL",
            "n_source_maf_ge_0.05_snps": len(pc_keys),
            "n_eur_reference_matched_snps": reference_counts["pcQTL"],
            "n_eur_reference_common_snps": reference_counts["pcQTL_common"],
            "n_shared_source_snps": shared,
            "n_shared_eur_reference_snps": reference_counts["shared"],
            "n_shared_eur_reference_common_snps": reference_counts["shared_common"],
            "pct_eur_reference_common_annotation_shared": (
                100.0 * reference_counts["shared_common"] / reference_counts["pcQTL_common"]
            ),
        },
        {
            "annotation": "eQTL",
            "n_source_maf_ge_0.05_snps": eq_total,
            "n_eur_reference_matched_snps": reference_counts["eQTL"],
            "n_eur_reference_common_snps": reference_counts["eQTL_common"],
            "n_shared_source_snps": shared,
            "n_shared_eur_reference_snps": reference_counts["shared"],
            "n_shared_eur_reference_common_snps": reference_counts["shared_common"],
            "pct_eur_reference_common_annotation_shared": (
                100.0 * reference_counts["shared_common"] / reference_counts["eQTL_common"]
            ),
        },
    ]
    write_tsv(
        target,
        rows,
        [
            "annotation", "n_source_maf_ge_0.05_snps",
            "n_eur_reference_matched_snps", "n_eur_reference_common_snps",
            "n_shared_source_snps", "n_shared_eur_reference_snps",
            "n_shared_eur_reference_common_snps",
            "pct_eur_reference_common_annotation_shared",
        ],
    )
    return len(rows)


def make_sldsc_data(outdir: Path) -> dict[str, int]:
    """Collect the machine-readable S-LDSC outputs used by the manuscript."""
    sldsc_out = outdir / "sldsc_heritability"
    sldsc_out.mkdir(parents=True, exist_ok=True)
    overlap_name = "prespecified247_annotation_overlap_summary.tsv"
    expected_names = set(SLDSC_TABLES.values()) | SLDSC_DERIVED_TABLES
    for existing in sldsc_out.glob("*.tsv"):
        if existing.name not in expected_names:
            existing.unlink()
    counts: dict[str, int] = {}
    for source_key, file_name in SLDSC_TABLES.items():
        target = sldsc_out / file_name
        copy_sldsc_table_for_submission(source_key, SOURCES[source_key], target)
        counts[file_name] = max(0, sum(1 for _ in target.open(encoding="utf-8")) - 1)

    overlap_file = sldsc_out / overlap_name
    counts[overlap_name] = make_sldsc_annotation_overlap_table(overlap_file)

    readme = sldsc_out / "README.md"
    readme.write_text(
        "# S-LDSC Heritability Results\n\n"
        "These machine-readable files support Supplementary Figures S8--S9 and\n"
        "Supplementary Table S5. Annotations are dense sets of all per-feature\n"
        "BH-significant cis-QTL variants, with both QTL types restricted to MAF >= 0.05,\n"
        "and conditioned on baseline-LD v2.2 after removal of its four built-in\n"
        "molecular-QTL MaxCPP columns. Trait selection starts from 247 prespecified\n"
        "FinnGen R12 core-disease traits; 91 pass h2 Z > 4 and 35\n"
        "pairwise-rg-clumped traits enter the primary meta-analysis.\n\n"
        "Included analyses are marginal pcQTL/eQTL models, the joint conditional "
        "model, per-trait and meta-analyzed pcQTL-minus-eQTL delta tau*, "
        "promoter/enhancer interaction models, pcQTL cluster-size strata, and "
        "genetic-correlation-threshold sensitivity summaries.\n\n"
        "## File manifest\n\n"
        "| File | Rows | Content |\n"
        "|---|---:|---|\n"
        f"| `prespecified247_marg_pc_annotation_meta.tsv` | {counts['prespecified247_marg_pc_annotation_meta.tsv']} | Marginal pcQTL enrichment and tau* meta-analysis |\n"
        f"| `prespecified247_marg_eq_annotation_meta.tsv` | {counts['prespecified247_marg_eq_annotation_meta.tsv']} | Marginal single-gene eQTL enrichment and tau* meta-analysis |\n"
        f"| `prespecified247_joint_per_endpoint.tsv` | {counts['prespecified247_joint_per_endpoint.tsv']} | Joint pcQTL/eQTL conditional estimates for each trait |\n"
        f"| `prespecified247_joint_annotation_meta.tsv` | {counts['prespecified247_joint_annotation_meta.tsv']} | Joint conditional annotation meta-analysis |\n"
        f"| `prespecified247_joint_delta.tsv` | {counts['prespecified247_joint_delta.tsv']} | Per-trait pcQTL-minus-eQTL delta tau* |\n"
        f"| `prespecified247_joint_delta_meta.tsv` | {counts['prespecified247_joint_delta_meta.tsv']} | Genome-wide delta tau* meta-analysis and sensitivity estimators |\n"
        f"| `prespecified247_annotation_overlap_summary.tsv` | {counts['prespecified247_annotation_overlap_summary.tsv']} | Source, EUR-reference-matched, and EUR-reference-common pcQTL/eQTL SNP-set sizes and overlap |\n"
        f"| `prespecified247_joint_annotation_meta_rg0p6_sensitivity.tsv` | {counts['prespecified247_joint_annotation_meta_rg0p6_sensitivity.tsv']} | Joint-model sensitivity at absolute rg 0.6 |\n"
        f"| `prespecified247_joint_annotation_meta_rg0p8_sensitivity.tsv` | {counts['prespecified247_joint_annotation_meta_rg0p8_sensitivity.tsv']} | Joint-model sensitivity at absolute rg 0.8 |\n"
        f"| `prespecified247_joint_delta_meta_rg0p6_sensitivity.tsv` | {counts['prespecified247_joint_delta_meta_rg0p6_sensitivity.tsv']} | pcQTL-minus-eQTL delta sensitivity at absolute rg 0.6 |\n"
        f"| `prespecified247_joint_delta_meta_rg0p8_sensitivity.tsv` | {counts['prespecified247_joint_delta_meta_rg0p8_sensitivity.tsv']} | pcQTL-minus-eQTL delta sensitivity at absolute rg 0.8 |\n"
        f"| `prespecified247_promoter_enhancer_annotation_meta.tsv` | {counts['prespecified247_promoter_enhancer_annotation_meta.tsv']} | Promoter/enhancer interaction annotation meta-analysis |\n"
        f"| `prespecified247_promoter_delta_per_endpoint.tsv` | {counts['prespecified247_promoter_delta_per_endpoint.tsv']} | Per-trait promoter delta tau* |\n"
        f"| `prespecified247_promoter_delta_meta.tsv` | {counts['prespecified247_promoter_delta_meta.tsv']} | Promoter delta tau* meta-analysis |\n"
        f"| `prespecified247_enhancer_delta_per_endpoint.tsv` | {counts['prespecified247_enhancer_delta_per_endpoint.tsv']} | Per-trait enhancer delta tau* |\n"
        f"| `prespecified247_enhancer_delta_meta.tsv` | {counts['prespecified247_enhancer_delta_meta.tsv']} | Enhancer delta tau* meta-analysis |\n"
        f"| `prespecified247_cluster_size_per_endpoint.tsv` | {counts['prespecified247_cluster_size_per_endpoint.tsv']} | Per-trait pcQTL estimates stratified by cluster size |\n"
        f"| `prespecified247_cluster_size_annotation_meta.tsv` | {counts['prespecified247_cluster_size_annotation_meta.tsv']} | Cluster-size-stratified annotation meta-analysis |\n"
        f"| `prespecified247_cluster_size_delta_per_endpoint.tsv` | {counts['prespecified247_cluster_size_delta_per_endpoint.tsv']} | Per-trait cluster-size delta tau* |\n"
        f"| `prespecified247_cluster_size_delta_meta.tsv` | {counts['prespecified247_cluster_size_delta_meta.tsv']} | Cluster-size delta tau* meta-analysis |\n"
        f"| `prespecified247_trait_selection_rg0p7.tsv` | {counts['prespecified247_trait_selection_rg0p7.tsv']} | Primary trait selection |\n"
        f"| `prespecified247_trait_selection_summary.tsv` | {counts['prespecified247_trait_selection_summary.tsv']} | Trait counts across rg thresholds |\n"
        f"| `prespecified247_manifest_provenance.tsv` | {counts['prespecified247_manifest_provenance.tsv']} | Prespecified trait-manifest provenance |\n\n"
        "## Key fields\n\n"
        "- `phenocode` and `phenotype`: FinnGen trait identifiers.\n"
        "- `annotation`: pcQTL/eQTL annotation or functional/cluster-size stratum.\n"
        "- `enrichment` and `tau_star`: S-LDSC enrichment and standardized per-SNP effect.\n"
        "- `delta_tau_star`: paired pcQTL-minus-eQTL standardized effect.\n"
        "- `se`, `ci95_lo`, `ci95_hi`, and `p`: uncertainty and test summaries.\n"
        "- `model` or `se_method`: fixed-effect, DerSimonian-Laird, Knapp-Hartung, or sign-test estimator where applicable.\n"
        "- `source_url`: public FinnGen R12 summary-statistics URL; local analysis paths are not included.\n"
        "- `source_manifest_file` and `source_manifest_sha256`: source-manifest filename and checksum.\n\n"
        "## Provenance\n\n"
        "Source tables are prespecified-trait outputs from workflow module 09. The\n"
        "model uses baseline-LD v2.2 after removal of its four built-in molecular-QTL\n"
        "MaxCPP annotations, 1000 Genomes Phase 3 European LD resources, HapMap3\n"
        "regression SNPs, and MHC exclusion.\n",
        encoding="utf-8",
    )
    return counts


def copy_table_with_renamed_columns(source: Path, target: Path, rename_map: dict[str, str]) -> None:
    with source.open(newline="", encoding="utf-8") as inp:
        reader = csv.DictReader(inp, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in {source}")
        fieldnames = [rename_map.get(field, field) for field in reader.fieldnames]
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", newline="", encoding="utf-8") as out:
            writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            for row in reader:
                writer.writerow({rename_map.get(field, field): value for field, value in row.items()})


def copy_enrichment_table_for_submission(source: Path, target: Path) -> None:
    """Keep source workflow keys internal while using submission-facing method names."""
    fieldnames = [
        "method",
        "stratum",
        "annotation",
        "status",
        "min_expected",
        "n_sets",
        "n_correlated",
        "n_null",
        "annotated_correlated",
        "annotated_null",
        "beta",
        "se",
        "OR",
        "CI_low",
        "CI_high",
        "pvalue",
        "fdr",
    ]
    with source.open(newline="", encoding="utf-8") as inp:
        reader = csv.DictReader(inp, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in {source}")
        missing = set(fieldnames) - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing enrichment columns in {source}: {sorted(missing)}")
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", newline="", encoding="utf-8") as out:
            writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            for row in reader:
                if row.get("method") == "add_cov_sc_hurdle":
                    row = dict(row)
                    row["method"] = "add_cov_sc_pcqtl"
                writer.writerow({field: row[field] if row[field] else "NA" for field in fieldnames})


def filter_table(source: Path, target: Path, predicate) -> int:
    with source.open(newline="", encoding="utf-8") as inp:
        reader = csv.DictReader(inp, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in {source}")
        target.parent.mkdir(parents=True, exist_ok=True)
        count = 0
        with target.open("w", newline="", encoding="utf-8") as out:
            writer = csv.DictWriter(out, delimiter="\t", fieldnames=reader.fieldnames, lineterminator="\n")
            writer.writeheader()
            for row in reader:
                if predicate(row):
                    writer.writerow(row)
                    count += 1
    return count


def row_count(path: Path) -> int:
    return max(0, sum(1 for _ in path.open(encoding="utf-8")) - 1)


def write_table_readme(outdir: Path) -> None:
    lines = [
        "# Displayed Supplementary Tables",
        "",
        "Only tables rendered in `supplementary_material.pdf` receive Supplementary",
        "Table numbers. They are numbered consecutively in order of appearance.",
        "Machine-readable result files that are not displayed in the document are stored",
        "under `../data/` and are not Supplementary Tables.",
        "",
        "| Table | Machine-readable file | Description |",
        "|---|---|---|",
        f"| S1 | `{FORMAL_TABLES['s1']}` | Cell-type eligibility and exclusion QC |",
        f"| S2 | `{FORMAL_TABLES['s2']}` | Joint-score sensitivity of pairwise associations and cluster calls |",
        f"| S3 | `{FORMAL_TABLES['s3']}` | Cross-mappability results for the representative GIMAP cluster |",
        f"| S4 | `{FORMAL_TABLES['s4']}` | Concordance of the primary and ACAT-BH pcQTL discovery criteria |",
        f"| S5 | `{FORMAL_TABLES['s5']}` | Primary S-LDSC meta-analysis summary |",
    ]
    (outdir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_data_readme(outdir: Path) -> None:
    lines = [
        "# Machine-Readable Supplementary Data",
        "",
        "These files provide complete result tables supporting the manuscript figures",
        "and displayed Supplementary Tables. They are archived as Supplementary Data and",
        "do not carry Supplementary Table numbers.",
        "All files contain aggregate or public summary statistics. No individual-level",
        "genotypes or expression measurements, donor or cell identifiers, restricted",
        "participant data, credentials, or local analysis paths are included.",
        "",
        "| File | Content |",
        "|---|---|",
        f"| `{DATA_FILES['celltype_cluster_pcqtl_summary']}` | Per-cell-type cluster and pcQTL output summary |",
        f"| `{DATA_FILES['cluster_annotation_enrichment']}` | Cluster annotation-enrichment model results |",
        f"| `{DATA_FILES['susie_strict_graph_edges']}` | Strict SuSiE graph edges |",
        f"| `{DATA_FILES['susie_strict_signal_groups']}` | Strict SuSiE connected-component signal groups |",
        f"| `{DATA_FILES['strict_pcqtl_specific_loading_effects']}` | pcQTL-specific loading, nominal-effect, and annotation summaries |",
        f"| `{DATA_FILES['susie_strict_threshold_sensitivity']}` | Strict SuSiE threshold-sensitivity summary |",
        f"| `{DATA_FILES['cross_mappability_sensitivity']}` | Cross-mappability credible-set yield and exclusion sensitivity results |",
        "| `sldsc_heritability/` | Complete S-LDSC marginal, joint, trait-level, and secondary analyses |",
    ]
    (outdir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_outputs(root: Path, tables_dir: Path, data_dir: Path) -> dict[str, int]:
    require_sources()
    tables_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    make_table_s1(tables_dir)
    copy_table(
        SOURCES["joint_score_cluster_sensitivity"],
        tables_dir / FORMAL_TABLES["s2"],
    )
    copy_table(
        SOURCES["crossmap_gimap_pairs"],
        tables_dir / FORMAL_TABLES["s3"],
    )
    make_table_s4(tables_dir)
    make_table_s5(tables_dir)

    make_celltype_cluster_pcqtl_data(data_dir)
    copy_enrichment_table_for_submission(
        SOURCES["enrichment"],
        data_dir / DATA_FILES["cluster_annotation_enrichment"],
    )
    copy_table(
        SOURCES["strict_susie_edges"],
        data_dir / DATA_FILES["susie_strict_graph_edges"],
    )
    copy_table(
        SOURCES["strict_susie_groups"],
        data_dir / DATA_FILES["susie_strict_signal_groups"],
    )
    copy_table_with_explicit_missing(
        SOURCES["strict_pcqtl_mechanism"],
        data_dir / DATA_FILES["strict_pcqtl_specific_loading_effects"],
    )
    copy_table(
        SOURCES["strict_susie_sensitivity"],
        data_dir / DATA_FILES["susie_strict_threshold_sensitivity"],
    )
    copy_crossmap_table_for_submission(
        SOURCES["crossmap_summary"],
        data_dir / DATA_FILES["cross_mappability_sensitivity"],
    )
    sldsc_counts = make_sldsc_data(data_dir)

    write_table_readme(tables_dir)
    write_data_readme(data_dir)
    counts = {
        key: row_count(tables_dir / file_name)
        for key, file_name in FORMAL_TABLES.items()
    }
    counts.update({
        key: row_count(data_dir / file_name)
        for key, file_name in DATA_FILES.items()
    })
    counts["sldsc_files"] = len(sldsc_counts)
    return counts


def main() -> int:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("repo_root", nargs="?", default=".")
    parser.add_argument(
        "--sldsc-only",
        action="store_true",
        help="regenerate only the machine-readable S-LDSC Supplementary Data",
    )
    args = parser.parse_args()

    root = Path(args.repo_root).resolve()
    tables_dir = root / "supplementary/tables"
    data_dir = root / "supplementary/data"
    if args.sldsc_only:
        require_sldsc_sources()
        counts = make_sldsc_data(data_dir)
        print(f"data.sldsc={data_dir / 'sldsc_heritability'}\tfiles={len(counts)}")
        return 0

    counts = generate_outputs(root, tables_dir, data_dir)
    for key, file_name in FORMAL_TABLES.items():
        print(f"{key}={tables_dir / file_name}\trows={counts[key]}")
    for key, file_name in DATA_FILES.items():
        print(f"data.{key}={data_dir / file_name}\trows={counts[key]}")
    print(f"data.sldsc={data_dir / 'sldsc_heritability'}\tfiles={counts['sldsc_files']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
