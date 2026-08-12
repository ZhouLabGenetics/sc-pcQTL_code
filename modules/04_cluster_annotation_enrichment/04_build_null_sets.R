#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))
source("config/config.R")

cluster_dir <- file.path(ROOT_DIR, "results", "cluster_sets")
out_dir <- file.path(ROOT_DIR, "results", "null_sets")
dir_create(out_dir)

genes <- fread(file.path(cluster_dir, "expressed_gene_universe.tsv"))

is_window_correlated <- function(window_genes, sig_set) {
  if (length(window_genes) < 2) return(FALSE)
  cmb <- combn(window_genes, 2)
  keys <- apply(cmb, 2, function(x) paste(sort(x), collapse = "\t"))
  mean(keys %in% sig_set) >= CLUSTER_THRESHOLD
}

load_sig_set_add_cov <- function(ct) {
  files <- list.files(
    file.path(ADD_COV_ROOT, "celltypes", ct, "cluster_identification", "results", "method2_sc_hurdle", "gene_associations_chunked"),
    pattern = "^significant_pairs\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  x <- rbindlist(lapply(files, safe_fread), fill = TRUE)
  if (!nrow(x)) return(character())
  unique(apply(x[, .(Gene1, Gene2)], 1, function(z) paste(sort(z), collapse = "\t")))
}

method <- "add_cov_sc_hurdle"
clusters <- safe_fread(file.path(cluster_dir, "add_cov_correlated_clusters.tsv"))
real_keys <- if (nrow(clusters)) unique(clusters[, paste(celltype, chr, genes, sep = "__")]) else character()
nulls <- list()
ni <- 0L
for (ct in unique(genes$celltype)) {
  sig_set <- load_sig_set_add_cov(ct)
  for (chr_i in sort(unique(genes[celltype == ct, chr]))) {
    chr_genes <- genes[celltype == ct & chr == chr_i][order(start), gene_name]
    chr_info <- genes[celltype == ct & chr == chr_i]
    if (length(chr_genes) < 2) next
    for (size in MAIN_CLUSTER_SIZES) {
      if (length(chr_genes) < size) next
      for (start_idx in seq_len(length(chr_genes) - size + 1L)) {
        window_genes <- chr_genes[start_idx:(start_idx + size - 1L)]
        gene_string <- paste(window_genes, collapse = ",")
        if (paste(ct, chr_i, gene_string, sep = "__") %in% real_keys) next
        if (is_window_correlated(window_genes, sig_set)) next
        gp <- chr_info[gene_name %in% window_genes]
        ni <- ni + 1L
        nulls[[ni]] <- data.table(
          method = method,
          celltype = ct,
          cluster_id = sprintf("NULL_%s_%s_chr%s_%05d", method, ct, chr_i, ni),
          chr = chr_i,
          num_genes = size,
          cluster_start = min(gp$start),
          cluster_end = max(gp$end),
          cluster_length = pmax(1L, max(gp$end) - min(gp$start)),
          genes = gene_string,
          is_correlated_cluster = FALSE
        )
      }
    }
  }
}
null_dt <- rbindlist(nulls, fill = TRUE)
fwrite(null_dt, file.path(out_dir, "add_cov_sc_hurdle_null_neighbor_clusters.tsv"), sep = "\t")
message(method, ": null sets=", nrow(null_dt))
