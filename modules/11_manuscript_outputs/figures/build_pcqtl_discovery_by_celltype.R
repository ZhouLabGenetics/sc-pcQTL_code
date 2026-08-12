#!/usr/bin/env Rscript

# Build the cell-type pcQTL summary consumed by Figure 3.

script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
release_root <- normalizePath(file.path(dirname(script_arg), "..", "..", ".."), mustWork = TRUE)
eligibility <- read.delim(
  file.path(release_root, "config", "celltype_eligibility.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
eligibility <- eligibility[eligibility$include_primary, , drop = FALSE]

workflow_root <- Sys.getenv(
  "COQTL_WORKFLOW_ROOT",
  unset = Sys.getenv("SC_PCQTL_WORKFLOW_ROOT", unset = "")
)
if (!nzchar(workflow_root)) {
  stop("Set COQTL_WORKFLOW_ROOT or SC_PCQTL_WORKFLOW_ROOT.")
}
result_root <- Sys.getenv("SC_PCQTL_PRIMARY_RESULT_ROOT", unset = workflow_root)

manuscript_root <- Sys.getenv("SC_PCQTL_MANUSCRIPT_ROOT", unset = "")
if (!nzchar(manuscript_root) && file.exists(file.path(getwd(), "main.tex"))) {
  manuscript_root <- getwd()
}
if (!nzchar(manuscript_root) || !file.exists(file.path(manuscript_root, "main.tex"))) {
  stop("Set SC_PCQTL_MANUSCRIPT_ROOT to the manuscript checkout.")
}

input_file <- file.path(
  result_root,
  "03_analysis_celltypes",
  "02_downstream_analysis_modules_add_cov_fdr",
  "pcqtl_compare",
  "data",
  "qtl_counts_per_celltype.tsv"
)
if (!file.exists(input_file)) stop("Missing pcQTL count summary: ", input_file)
upstream_celltype_root <- Sys.getenv(
  "COQTL_UPSTREAM_CELLTYPES_DIR",
  unset = file.path(
    workflow_root,
    "03_analysis_celltypes",
    "01_upstream_main_pipeline_add_cov",
    "celltypes"
  )
)

canonical <- data.frame(
  celltype = eligibility$celltype,
  display_label = eligibility$eqtl_celltype,
  stringsAsFactors = FALSE
)

counts <- read.delim(input_file, stringsAsFactors = FALSE, check.names = FALSE)
counts <- counts[counts$celltype %in% canonical$celltype, , drop = FALSE]
required <- c("celltype", "n_tests", "n_sig")
if (!all(required %in% names(counts))) {
  stop("Input is missing required columns: ", paste(setdiff(required, names(counts)), collapse = ", "))
}
if (anyDuplicated(counts$celltype)) stop("Duplicate cell types in ", input_file)
if (any(!counts$celltype %in% canonical$celltype)) {
  stop("Unexpected cell types in input: ", paste(setdiff(counts$celltype, canonical$celltype), collapse = ", "))
}

selected_counts <- do.call(rbind, lapply(canonical$celltype, function(celltype) {
  pc_file <- file.path(
    upstream_celltype_root,
    celltype,
    "pcQTL",
    "step2_pca",
    "pc_list_for_saige.tsv"
  )
  if (!file.exists(pc_file)) stop("Missing selected cluster-PC list: ", pc_file)
  pc_dt <- read.delim(pc_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("cluster_id", "PC") %in% names(pc_dt))) {
    stop("Selected cluster-PC list is missing cluster_id/PC: ", pc_file)
  }
  if (anyDuplicated(paste(pc_dt$cluster_id, pc_dt$PC, sep = "\t"))) {
    stop("Duplicate cluster-PC phenotypes in ", pc_file)
  }
  data.frame(celltype = celltype, n_selected = nrow(pc_dt), stringsAsFactors = FALSE)
}))

plot_dt <- merge(canonical, selected_counts, by = "celltype", all.x = TRUE)
plot_dt <- merge(plot_dt, counts[, c("celltype", "n_tests", "n_sig")], by = "celltype", all.x = TRUE)
plot_dt$n_tests[is.na(plot_dt$n_tests)] <- 0L
plot_dt$n_sig[is.na(plot_dt$n_sig)] <- 0L
plot_dt$n_selected <- as.integer(plot_dt$n_selected)
plot_dt$n_tests <- as.integer(plot_dt$n_tests)
plot_dt$n_sig <- as.integer(plot_dt$n_sig)
if (any(plot_dt$n_tests > plot_dt$n_selected)) stop("Successfully tested count exceeds selected count.")
if (any(plot_dt$n_sig > plot_dt$n_tests)) stop("FDR-significant count exceeds tested count.")
plot_dt$n_failed_convergence <- plot_dt$n_selected - plot_dt$n_tests
if (
  sum(plot_dt$n_selected) != 4430L ||
  sum(plot_dt$n_tests) != 4353L ||
  sum(plot_dt$n_failed_convergence) != 77L ||
  sum(plot_dt$n_sig) != 2040L
) {
  stop(
    "Primary 10-cell-type totals differ from the analysis contract: selected=", sum(plot_dt$n_selected),
    ", successfully tested=", sum(plot_dt$n_tests),
    ", failed convergence=", sum(plot_dt$n_failed_convergence),
    ", within-phenotype FDR hits=", sum(plot_dt$n_sig)
  )
}

plot_dt$percent_fdr_significant <- ifelse(
  plot_dt$n_tests > 0L,
  100 * plot_dt$n_sig / plot_dt$n_tests,
  NA_real_
)
plot_dt <- plot_dt[order(plot_dt$n_sig, plot_dt$n_tests, plot_dt$display_label), ]
plot_dt$display_label <- factor(plot_dt$display_label, levels = plot_dt$display_label)

data_out <- file.path(manuscript_root, "data", "supp_pcqtl_discovery_by_celltype.tsv")
dir.create(dirname(data_out), recursive = TRUE, showWarnings = FALSE)

source_dt <- data.frame(
  celltype = plot_dt$celltype,
  display_label = as.character(plot_dt$display_label),
  n_selected_cluster_pc_phenotypes = plot_dt$n_selected,
  n_successfully_tested_cluster_pc_phenotypes = plot_dt$n_tests,
  n_failed_convergence_cluster_pc_phenotypes = plot_dt$n_failed_convergence,
  n_within_phenotype_fdr_hit_cluster_pc_phenotypes = plot_dt$n_sig,
  percent_within_phenotype_fdr_hit_among_successfully_tested = plot_dt$percent_fdr_significant
)
write.table(source_dt, data_out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
message("Wrote source data: ", data_out)
