#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))
source("config/config.R")

out_dir <- file.path(ROOT_DIR, "results", "cluster_sets")
dir_create(out_dir)

all_genes <- list()
all_clusters <- list()
all_assign <- list()

for (ct in CELLTYPES) {
  base <- file.path(ADD_COV_ROOT, "celltypes", ct, "cluster_identification", "results", "method2_sc_hurdle")
  filt_file <- file.path(base, "filtered_genes.tsv")
  cl_sum_file <- file.path(base, "merged_clusters", "cluster_summary.tsv")
  cl_assign_file <- file.path(base, "merged_clusters", "cluster_gene_assignments.tsv")
  if (!file.exists(filt_file)) {
    warning("Missing filtered genes for ", ct)
    next
  }
  genes <- fread(filt_file)
  genes <- genes[keep == TRUE & !is.na(chr_numeric) & chr_numeric %in% 1:22 & !is.na(start) & !is.na(end)]
  genes[, `:=`(
    celltype = ct,
    chr = std_chr(chr_numeric),
    start = as.integer(start),
    end = as.integer(end)
  )]
  all_genes[[ct]] <- unique(genes[, .(celltype, gene_name, chr, start, end)])

  if (file.exists(cl_sum_file)) {
    cl <- fread(cl_sum_file)
    cl <- cl[cluster_size %in% MAIN_CLUSTER_SIZES]
    cl[, `:=`(
      method = "add_cov_sc_hurdle",
      celltype = ct,
      chr = std_chr(chromosome),
      num_genes = as.integer(cluster_size),
      cluster_start = as.integer(start_position),
      cluster_end = as.integer(end_position),
      cluster_length = pmax(1L, as.integer(cluster_span_bp)),
      is_correlated_cluster = TRUE
    )]
    all_clusters[[ct]] <- cl[, .(
      method, celltype, cluster_id, chr, num_genes, cluster_start, cluster_end,
      cluster_length, genes, is_correlated_cluster
    )]
  }

  if (file.exists(cl_assign_file)) {
    asn <- fread(cl_assign_file)
    asn <- asn[cluster_size %in% MAIN_CLUSTER_SIZES]
    asn[, `:=`(method = "add_cov_sc_hurdle", celltype = ct)]
    all_assign[[ct]] <- asn[, .(method, celltype, cluster_id, gene_name)]
  }
}

gene_universe <- rbindlist(all_genes, fill = TRUE)
setorder(gene_universe, celltype, chr, start, end)
fwrite(gene_universe, file.path(out_dir, "expressed_gene_universe.tsv"), sep = "\t")

clusters <- rbindlist(all_clusters, fill = TRUE)
assignments <- rbindlist(all_assign, fill = TRUE)
fwrite(clusters, file.path(out_dir, "add_cov_correlated_clusters.tsv"), sep = "\t")
fwrite(assignments, file.path(out_dir, "add_cov_cluster_gene_assignments.tsv"), sep = "\t")

summary_dt <- clusters[, .N, by = .(method, celltype, num_genes)]
fwrite(summary_dt, file.path(out_dir, "add_cov_cluster_count_summary.tsv"), sep = "\t")
print(summary_dt)
