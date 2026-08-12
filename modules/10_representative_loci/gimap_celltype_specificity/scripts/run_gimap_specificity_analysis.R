#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
release_root <- normalizePath(file.path(dirname(script_arg), "..", "..", "..", ".."), mustWork = TRUE)
source(file.path(release_root, "config", "celltype_eligibility.R"))

require_dir <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value) || !dir.exists(value)) stop("Set ", name, " to an existing directory.")
  normalizePath(value, mustWork = TRUE)
}

upstream_root <- require_dir("SC_PCQTL_UPSTREAM_CELLTYPES_ROOT")
pcqtl_root <- require_dir("SC_PCQTL_PCQTL_SUMMARY_ROOT")
output_root <- Sys.getenv("SC_PCQTL_GIMAP_SPECIFICITY_ROOT", unset = "")
if (!nzchar(output_root)) stop("Set SC_PCQTL_GIMAP_SPECIFICITY_ROOT.")
output_root <- normalizePath(output_root, mustWork = FALSE)
dir.create(file.path(output_root, "data"), recursive = TRUE, showWarnings = FALSE)

pcqtl_file <- file.path(pcqtl_root, "data/all_pcqtl_results.tsv")
if (!file.exists(pcqtl_file)) stop("Missing final pcQTL table: ", pcqtl_file)

gimap_genes <- c("GIMAP8", "GIMAP7", "GIMAP4", "GIMAP6", "GIMAP2", "GIMAP1", "GIMAP5")
celltypes <- primary_celltypes(file.path(release_root, "config", "celltype_eligibility.tsv"))
split_genes <- function(value) {
  value <- gsub("[; ]+", ",", as.character(value))
  unique(Filter(nzchar, unlist(strsplit(value, ",", fixed = TRUE))))
}

pcqtl <- fread(pcqtl_file)
required_pcqtl <- c("celltype", "cluster_id", "genes", "PC", "min_snp_p", "min_snp_fdr", "ACAT_p")
missing_pcqtl <- setdiff(required_pcqtl, names(pcqtl))
if (length(missing_pcqtl)) stop("Missing pcQTL columns: ", paste(missing_pcqtl, collapse = ", "))

pcqtl[, gimap_gene_list := vapply(
  genes,
  function(value) paste(intersect(split_genes(value), gimap_genes), collapse = ","),
  character(1)
)]
pcqtl[, n_gimap_genes := vapply(
  gimap_gene_list,
  function(value) if (nzchar(value)) length(strsplit(value, ",", fixed = TRUE)[[1]]) else 0L,
  integer(1)
)]
pcqtl[, `:=`(
  min_snp_p_num = suppressWarnings(as.numeric(min_snp_p)),
  min_snp_fdr_num = suppressWarnings(as.numeric(min_snp_fdr)),
  ACAT_p_num = suppressWarnings(as.numeric(ACAT_p))
)]

gimap <- pcqtl[n_gimap_genes >= 2L]
if (!nrow(gimap)) stop("No multi-gene GIMAP cluster-PC phenotypes were found.")

best_by_cluster <- gimap[, {
  ranking <- fifelse(is.finite(min_snp_p_num), min_snp_p_num, Inf)
  .SD[order(ranking, PC)[1L]]
}, by = .(celltype, cluster_id)]

best_by_celltype <- best_by_cluster[, {
  ranking <- fifelse(is.finite(min_snp_p_num), min_snp_p_num, Inf)
  .SD[order(ranking, -n_gimap_genes, cluster_id)[1L]]
}, by = celltype]
best_output <- best_by_celltype[, .(
  celltype,
  cluster_id,
  n_gimap_genes,
  gimap_genes = gimap_gene_list,
  best_PC = PC,
  best_pc_min_snp_p = min_snp_p_num,
  best_pc_min_snp_fdr = min_snp_fdr_num,
  best_pc_ACAT_p = ACAT_p_num,
  best_pc_has_sig_snp = is.finite(min_snp_fdr_num) & min_snp_fdr_num < 0.05
)]
setorder(best_output, celltype)
fwrite(
  best_output,
  file.path(output_root, "data/sc_gimap_best_pcqtl_by_celltype.tsv"),
  sep = "\t", na = ""
)

module_rows <- lapply(celltypes, function(celltype_id) {
  assignment_file <- file.path(
    upstream_root,
    celltype_id,
    "cluster_identification/results/method2_sc_hurdle/merged_clusters/cluster_gene_assignments.tsv"
  )
  if (!file.exists(assignment_file)) {
    stop("Missing final cluster assignment table for ", celltype_id, ": ", assignment_file)
  }
  assignments <- fread(assignment_file)
  missing_assignment <- setdiff(c("cluster_id", "gene_name"), names(assignments))
  if (length(missing_assignment)) {
    stop(
      "Missing cluster-assignment columns for ", celltype_id, ": ",
      paste(missing_assignment, collapse = ", ")
    )
  }
  subset <- assignments[, {
    present <- gimap_genes[gimap_genes %in% unique(gene_name)]
    .(n_gimap_genes = length(present), gimap_gene_list = paste(present, collapse = ","))
  }, by = cluster_id][n_gimap_genes > 0L]
  if (!nrow(subset)) {
    return(data.table(
      celltype = celltype_id,
      max_gimap_genes_in_any_cluster = 0L,
      has_full_7_gene_cluster = FALSE,
      n_tied_largest_clusters = 0L,
      largest_gimap_cluster_id = NA_character_,
      largest_gimap_cluster_genes = NA_character_
    ))
  }
  max_count <- max(subset$n_gimap_genes)
  largest <- subset[n_gimap_genes == max_count][order(cluster_id)]
  data.table(
    celltype = celltype_id,
    max_gimap_genes_in_any_cluster = max_count,
    has_full_7_gene_cluster = max_count == length(gimap_genes),
    n_tied_largest_clusters = nrow(largest),
    largest_gimap_cluster_id = paste(largest$cluster_id, collapse = ";"),
    largest_gimap_cluster_genes = paste(largest$gimap_gene_list, collapse = ";")
  )
})
module_completeness <- rbindlist(module_rows, fill = TRUE)
setorder(module_completeness, celltype)
fwrite(
  module_completeness,
  file.path(output_root, "data/sc_gimap_module_completeness_by_celltype.tsv"),
  sep = "\t", na = ""
)

message("Wrote final GIMAP specificity tables to ", file.path(output_root, "data"))
