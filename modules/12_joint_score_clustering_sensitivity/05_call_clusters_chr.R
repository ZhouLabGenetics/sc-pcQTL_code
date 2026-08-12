#!/usr/bin/env Rscript

.local_lib <- Sys.getenv("SC_PCQTL_R_LIBS", unset = "")
if (nzchar(.local_lib)) .libPaths(c(.local_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args_all[grep("--file=", args_all)[1]]))
module_dir <- dirname(script_path)
source(file.path(module_dir, "load_config.R"))
cfg <- load_joint_score_config(module_dir)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 05_call_clusters_chr.R <chromosome>")
chr <- as.integer(args[1])
if (!chr %in% 1:22) stop("Chromosome must be an integer from 1 to 22")

max_window_size <- 50L
min_window_size <- 2L
cluster_threshold <- 0.70
dir.create(cfg$clusters_dir, recursive = TRUE, showWarnings = FALSE)

gene_info <- fread(cfg$gene_info_file)
filtered <- fread(file.path(cfg$method_dir, "filtered_genes.tsv"))[keep == TRUE, .(gene_name)]
gene_info <- unique(merge(gene_info, filtered, by = "gene_name"), by = "gene_name")
chr_genes <- gene_info[chr_numeric == chr][order(start, end, gene_name), gene_name]
n_genes <- length(chr_genes)

empty_clusters <- data.table(
  cluster_id = character(), chromosome = integer(), cluster_size = integer(),
  start_position = numeric(), end_position = numeric(), cluster_span_bp = numeric(),
  genes = character()
)
if (!n_genes) {
  fwrite(empty_clusters, file.path(cfg$clusters_dir, sprintf("chr%d_clusters.tsv", chr)), sep = "\t")
  fwrite(data.table(
    Chromosome = chr, NumClusters = 0L, NumGenes = 0L, AssignedGenes = 0L
  ), file.path(cfg$clusters_dir, sprintf("chr%d_cluster_summary.tsv", chr)), sep = "\t")
  message(sprintf("[%s] Chromosome %d: no filtered genes", cfg$celltype, chr))
  quit(status = 0L)
}

sig_file <- file.path(
  cfg$assoc_dir, sprintf("chr%d", chr), sprintf("chr%d_significant_pairs.tsv", chr)
)
if (!file.exists(sig_file)) stop("Missing chromosome association file: ", sig_file)
sig_data <- fread(sig_file)

sig_matrix <- matrix(FALSE, nrow = n_genes, ncol = n_genes)
rownames(sig_matrix) <- chr_genes
colnames(sig_matrix) <- chr_genes
diag(sig_matrix) <- TRUE
if (nrow(sig_data)) {
  keep <- sig_data$Gene1 %in% chr_genes & sig_data$Gene2 %in% chr_genes
  sig_data <- sig_data[keep]
  if (nrow(sig_data)) {
    left <- match(sig_data$Gene1, chr_genes)
    right <- match(sig_data$Gene2, chr_genes)
    sig_matrix[cbind(left, right)] <- TRUE
    sig_matrix[cbind(right, left)] <- TRUE
  }
}

check_cluster <- function(window_genes) {
  n <- length(window_genes)
  if (n < 2L) return(FALSE)
  sub_mat <- sig_matrix[window_genes, window_genes, drop = FALSE]
  mean(sub_mat[upper.tri(sub_mat)], na.rm = TRUE) >= cluster_threshold
}

assigned <- setNames(rep(FALSE, n_genes), chr_genes)
clusters <- list()
cluster_index <- 1L
for (window_size in max_window_size:min_window_size) {
  if (window_size > n_genes) next
  for (start_idx in seq_len(n_genes - window_size + 1L)) {
    end_idx <- start_idx + window_size - 1L
    window_genes <- chr_genes[start_idx:end_idx]
    if (any(assigned[window_genes])) next
    if (!check_cluster(window_genes)) next
    positions <- gene_info[gene_name %in% window_genes, .(gene_name, start, end)]
    clusters[[cluster_index]] <- data.table(
      cluster_id = sprintf("SC_chr%d_cluster_%03d", chr, cluster_index),
      chromosome = chr,
      cluster_size = window_size,
      start_position = min(positions$start),
      end_position = max(positions$end),
      cluster_span_bp = max(positions$end) - min(positions$start),
      genes = paste(window_genes, collapse = ",")
    )
    assigned[window_genes] <- TRUE
    cluster_index <- cluster_index + 1L
  }
}

cluster_dt <- if (length(clusters)) rbindlist(clusters) else empty_clusters
fwrite(cluster_dt, file.path(cfg$clusters_dir, sprintf("chr%d_clusters.tsv", chr)), sep = "\t")
fwrite(data.table(
  Chromosome = chr,
  NumClusters = nrow(cluster_dt),
  NumGenes = n_genes,
  AssignedGenes = sum(assigned)
), file.path(cfg$clusters_dir, sprintf("chr%d_cluster_summary.tsv", chr)), sep = "\t")
message(sprintf("[%s] Chromosome %d: %d clusters", cfg$celltype, chr, nrow(cluster_dt)))
