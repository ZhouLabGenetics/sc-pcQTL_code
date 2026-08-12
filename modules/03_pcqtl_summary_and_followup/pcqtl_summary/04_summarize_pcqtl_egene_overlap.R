#!/usr/bin/env Rscript

.local_libs <- c(
  Sys.getenv("SC_PCQTL_R_LIBS", unset = ""),
  Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = "")
)
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
module_dir <- normalizePath(dirname(args_file), mustWork = TRUE)
release_root <- normalizePath(file.path(module_dir, "..", "..", ".."), mustWork = TRUE)
source(file.path(release_root, "config", "celltype_eligibility.R"))
primary_ct <- primary_celltypes(file.path(release_root, "config", "celltype_eligibility.tsv"))

summary_root <- Sys.getenv("SC_PCQTL_PCQTL_SUMMARY_ROOT", unset = "")
if (!nzchar(summary_root)) stop("Set SC_PCQTL_PCQTL_SUMMARY_ROOT.")
overlap_root <- Sys.getenv(
  "SC_PCQTL_PCQTL_EGENE_ROOT",
  unset = file.path(dirname(normalizePath(summary_root, mustWork = FALSE)), "pcqtl_compare_saigeqtl")
)
data_dir <- file.path(overlap_root, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

sig_file <- file.path(summary_root, "data", "sig_qtls.tsv")
egene_file <- file.path(data_dir, "all_egenes_combined.tsv")
if (!file.exists(sig_file)) stop("Missing final pcQTL table: ", sig_file)
if (!file.exists(egene_file)) stop("Missing primary eGene table: ", egene_file)

sig_qtls <- fread(sig_file)
all_egenes <- fread(egene_file)
excluded_pc <- setdiff(unique(sig_qtls$celltype), primary_ct)
excluded_eq <- setdiff(unique(all_egenes$celltype), primary_ct)
if (length(excluded_pc) || length(excluded_eq)) {
  stop(
    "Input tables retain ineligible cell types. Rebuild module 03 first: ",
    paste(sort(unique(c(excluded_pc, excluded_eq))), collapse = ", ")
  )
}
if (!nrow(sig_qtls)) stop("No significant pcQTL phenotypes were found")

egene_lookup <- unique(all_egenes[is_egene == TRUE, .(celltype, gene)])
overlap <- rbindlist(lapply(seq_len(nrow(sig_qtls)), function(i) {
  row <- sig_qtls[i]
  genes <- trimws(strsplit(as.character(row$genes), ",", fixed = TRUE)[[1]])
  genes <- genes[nzchar(genes)]
  is_egene <- genes %chin% egene_lookup[celltype == row$celltype, gene]
  n_genes <- length(genes)
  n_egenes <- sum(is_egene)
  pct_egenes <- if (n_genes) 100 * n_egenes / n_genes else NA_real_
  category <- if (!n_genes) {
    "empty"
  } else if (n_egenes == 0L) {
    "novel"
  } else if (pct_egenes <= 50) {
    "partially_novel"
  } else if (pct_egenes < 100) {
    "mostly_known"
  } else {
    "all_known"
  }
  data.table(
    celltype = row$celltype,
    cluster_id = row$cluster_id,
    PC = row$PC,
    chr = row$chr,
    n_genes = n_genes,
    n_egenes = n_egenes,
    pct_egenes = pct_egenes,
    category = category,
    genes = paste(genes, collapse = ","),
    egene_genes = paste(genes[is_egene], collapse = ","),
    non_egene_genes = paste(genes[!is_egene], collapse = ","),
    min_snp_fdr = as.numeric(row$fdr)
  )
}), fill = TRUE)

category_dist <- overlap[, .N, by = category][order(-N)]
category_dist[, pct := round(100 * N / nrow(overlap), 2)]
summary_by_celltype <- overlap[, .(
  total_clusters = .N,
  novel = sum(category == "novel"),
  partially_novel = sum(category == "partially_novel"),
  mostly_known = sum(category == "mostly_known"),
  all_known = sum(category == "all_known"),
  pct_novel = round(100 * mean(category == "novel"), 2),
  pct_partially_novel = round(100 * mean(category == "partially_novel"), 2)
), by = celltype][order(-novel)]

fwrite(overlap, file.path(data_dir, "cluster_egene_overlap.tsv"), sep = "\t")
fwrite(summary_by_celltype, file.path(data_dir, "summary_by_celltype.tsv"), sep = "\t")
fwrite(category_dist, file.path(data_dir, "category_distribution.tsv"), sep = "\t")
message("Wrote primary-analysis pcQTL/eGene overlap summaries to ", data_dir)
