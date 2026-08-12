#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))
source("config/config.R")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}
ct <- get_arg("--celltype", Sys.getenv("CELLTYPE", NA_character_))
if (is.na(ct) || !nzchar(ct)) {
  idx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA_character_))
  if (is.na(idx)) stop("Provide --celltype or SLURM_ARRAY_TASK_ID")
  ct <- CELLTYPES[idx]
}
if (!ct %in% CELLTYPES) stop("Unknown celltype: ", ct)

out_dir <- file.path(ROOT_DIR, "results", "pb_spearman", ct)
dir_create(out_dir)

gene_universe_file <- file.path(ROOT_DIR, "results", "cluster_sets", "expressed_gene_universe.tsv")
if (!file.exists(gene_universe_file)) stop("Run 01_collect_add_cov_clusters.R first")
genes <- fread(gene_universe_file)[celltype == ct]
genes <- genes[chr %in% as.character(1:22)]
setorder(genes, chr, start, end)
gene_names <- unique(genes$gene_name)

count_file <- pb_input_file(ct)
if (!file.exists(count_file)) stop("Missing slim PB input for ", ct, ": ", count_file)

header <- names(fread_maybe_gz(count_file, nrows = 0, sep = "\t", showProgress = FALSE))
select_cols <- intersect(c("barcode", "CellID", "individual", "IndividualID", gene_names), header)
missing_genes <- setdiff(gene_names, header)
if (length(missing_genes)) {
  warning(sprintf("%s: %d expressed genes absent from count file", ct, length(missing_genes)))
}

message(sprintf("[%s] Reading slim PB input with %d selected genes", ct, length(intersect(gene_names, header))))
dt <- fread_maybe_gz(count_file, select = select_cols, sep = "\t", showProgress = TRUE)
if (!nrow(dt)) stop("No cells for celltype: ", ct)

use_genes <- intersect(gene_names, names(dt))
donor_col <- intersect(c("individual", "IndividualID"), names(dt))[1]
if (is.na(donor_col)) stop("No donor column found in ", count_file)
setnames(dt, donor_col, "individual")
pb <- dt[, lapply(.SD, mean, na.rm = TRUE), by = individual, .SDcols = use_genes]
pb_mat <- as.matrix(pb[, ..use_genes])
storage.mode(pb_mat) <- "numeric"
rownames(pb_mat) <- pb$individual
donor_n <- nrow(pb_mat)
if (donor_n < 4) stop("Too few donors for stable Spearman correlations: ", donor_n)

pair_rows <- list()
pair_i <- 0L
max_size <- max(MAIN_CLUSTER_SIZES)

for (chr_i in sort(unique(genes$chr))) {
  chr_genes <- genes[chr == chr_i][order(start), gene_name]
  chr_genes <- chr_genes[chr_genes %in% use_genes]
  n <- length(chr_genes)
  if (n < 2) next
  for (i in seq_len(n - 1L)) {
    jmax <- min(n, i + max_size - 1L)
    if (jmax <= i) next
    for (j in (i + 1L):jmax) {
      g1 <- chr_genes[i]
      g2 <- chr_genes[j]
      rho <- suppressWarnings(cor(pb_mat[, g1], pb_mat[, g2], method = "spearman", use = "pairwise.complete.obs"))
      pair_i <- pair_i + 1L
      pair_rows[[pair_i]] <- data.table(Gene1 = g1, Gene2 = g2, chr = chr_i, rho = rho)
    }
  }
}

pairs <- rbindlist(pair_rows, fill = TRUE)
fwrite(pairs, file.path(out_dir, "pb_spearman_local_pairs.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(data.table(
  celltype = ct,
  n_cells = nrow(dt),
  n_donors = donor_n,
  n_genes = length(use_genes),
  local_pairs_with_sign = nrow(pairs)
), file.path(out_dir, "pb_spearman_summary.tsv"), sep = "\t")

message(sprintf("[%s] Done: donor-level Spearman signs=%d", ct, nrow(pairs)))
