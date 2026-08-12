#!/usr/bin/env Rscript

.local_lib <- Sys.getenv("SC_PCQTL_R_LIBS", unset = "")
if (nzchar(.local_lib)) .libPaths(c(.local_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args_all[grep("--file=", args_all)[1]]))
module_dir <- dirname(script_path)
release_root <- normalizePath(file.path(module_dir, "..", ".."), mustWork = TRUE)
source(file.path(module_dir, "R", "joint_score_utils.R"))
setDTthreads(1L)

args <- commandArgs(trailingOnly = TRUE)
primary_root <- if (length(args) >= 1L) args[[1]] else Sys.getenv("SC_PCQTL_PRIMARY_CELLTYPE_ROOT", "")
joint_root <- if (length(args) >= 2L) args[[2]] else Sys.getenv("SC_PCQTL_JOINT_SCORE_CELLTYPE_ROOT", "")
output_dir <- if (length(args) >= 3L) args[[3]] else Sys.getenv("SC_PCQTL_JOINT_SCORE_COMPARISON_DIR", "")
if (!nzchar(primary_root) || !nzchar(joint_root) || !nzchar(output_dir)) {
  stop(
    "Usage: Rscript 07_compare_primary_joint_score.R ",
    "<primary_celltype_root> <joint_celltype_root> <output_dir>"
  )
}

manifest_path <- Sys.getenv(
  "SC_PCQTL_CELLTYPE_MANIFEST",
  unset = file.path(release_root, "config", "celltype_eligibility.tsv")
)
manifest <- fread(manifest_path)
celltypes <- manifest[include_primary == TRUE, celltype]
if (length(celltypes) != 10L) stop("Expected exactly 10 primary cell types")

method_dir <- function(root, celltype) {
  file.path(root, celltype, "cluster_identification", "results", "method2_sc_hurdle")
}

read_pair_set <- function(root, celltype) {
  base <- method_dir(root, celltype)
  files <- file.path(base, "gene_associations", sprintf("chr%d", 1:22), sprintf("chr%d_significant_pairs.tsv", 1:22))
  if (any(!file.exists(files))) {
    stop("Missing significant-pair files for ", celltype, ": ", paste(basename(files[!file.exists(files)]), collapse = ", "))
  }
  values <- rbindlist(lapply(files, function(path) {
    tab <- fread(path, select = c("Gene1", "Gene2"))
    if (!nrow(tab)) return(data.table(pair_key = character()))
    tab[, .(pair_key = canonical_pair_key(Gene1, Gene2))]
  }))
  unique(values$pair_key)
}

read_clusters <- function(root, celltype) {
  path <- file.path(method_dir(root, celltype), "merged_clusters", "cluster_summary.tsv")
  if (!file.exists(path)) stop("Missing cluster summary: ", path)
  clusters <- fread(path)
  required <- c("cluster_id", "chromosome", "genes")
  if (!all(required %in% names(clusters))) stop("Cluster summary lacks required columns: ", path)
  clusters[, gene_key := vapply(genes, canonical_gene_set, character(1))]
  clusters[, gene_list := lapply(gene_key, function(value) strsplit(value, ",", fixed = TRUE)[[1]])]
  clusters[, exact_key := paste(chromosome, gene_key, sep = "::")]
  clusters
}

rows <- list()
primary_best_all <- numeric()
joint_best_all <- numeric()
cluster_cache <- list()

for (celltype_id in celltypes) {
  primary_pairs <- read_pair_set(primary_root, celltype_id)
  joint_pairs <- read_pair_set(joint_root, celltype_id)
  shared_pairs <- intersect(primary_pairs, joint_pairs)

  primary_clusters <- read_clusters(primary_root, celltype_id)
  joint_clusters <- read_clusters(joint_root, celltype_id)
  primary_best <- best_jaccard_matches(
    primary_clusters$gene_list, primary_clusters$chromosome,
    joint_clusters$gene_list, joint_clusters$chromosome
  )
  joint_best <- best_jaccard_matches(
    joint_clusters$gene_list, joint_clusters$chromosome,
    primary_clusters$gene_list, primary_clusters$chromosome
  )
  exact_clusters <- length(intersect(primary_clusters$exact_key, joint_clusters$exact_key))

  primary_best_all <- c(primary_best_all, primary_best)
  joint_best_all <- c(joint_best_all, joint_best)
  cluster_cache[[celltype_id]] <- list(primary = primary_clusters, joint = joint_clusters)

  rows[[celltype_id]] <- data.table(
    celltype = celltype_id,
    primary_pairs = length(primary_pairs),
    joint_score_pairs = length(joint_pairs),
    shared_pairs = length(shared_pairs),
    primary_pair_recovered_fraction = length(shared_pairs) / length(primary_pairs),
    joint_pair_shared_fraction = length(shared_pairs) / length(joint_pairs),
    pair_jaccard = length(shared_pairs) / length(union(primary_pairs, joint_pairs)),
    primary_clusters = nrow(primary_clusters),
    joint_score_clusters = nrow(joint_clusters),
    exact_clusters = exact_clusters,
    primary_exact_cluster_fraction = exact_clusters / nrow(primary_clusters),
    joint_exact_cluster_fraction = exact_clusters / nrow(joint_clusters),
    primary_best_jaccard_median = median(primary_best),
    primary_best_jaccard_ge_0.5_fraction = mean(primary_best >= 0.5),
    joint_best_jaccard_median = median(joint_best),
    joint_best_jaccard_ge_0.5_fraction = mean(joint_best >= 0.5)
  )
}

summary <- rbindlist(rows)
aggregate_row <- summary[, .(
  celltype = "ALL",
  primary_pairs = sum(primary_pairs),
  joint_score_pairs = sum(joint_score_pairs),
  shared_pairs = sum(shared_pairs),
  primary_pair_recovered_fraction = sum(shared_pairs) / sum(primary_pairs),
  joint_pair_shared_fraction = sum(shared_pairs) / sum(joint_score_pairs),
  pair_jaccard = sum(shared_pairs) / (sum(primary_pairs) + sum(joint_score_pairs) - sum(shared_pairs)),
  primary_clusters = sum(primary_clusters),
  joint_score_clusters = sum(joint_score_clusters),
  exact_clusters = sum(exact_clusters),
  primary_exact_cluster_fraction = sum(exact_clusters) / sum(primary_clusters),
  joint_exact_cluster_fraction = sum(exact_clusters) / sum(joint_score_clusters),
  primary_best_jaccard_median = median(primary_best_all),
  primary_best_jaccard_ge_0.5_fraction = mean(primary_best_all >= 0.5),
  joint_best_jaccard_median = median(joint_best_all),
  joint_best_jaccard_ge_0.5_fraction = mean(joint_best_all >= 0.5)
)]
summary <- rbind(summary, aggregate_row, fill = TRUE)

fraction_cols <- grep("fraction$|jaccard$|jaccard_median$", names(summary), value = TRUE)
for (column in fraction_cols) {
  set(summary, j = column, value = sprintf("%.6f", summary[[column]]))
}

representative_path <- file.path(module_dir, "config", "representative_clusters.tsv")
representatives <- fread(representative_path)
representative_retention <- rbindlist(lapply(seq_len(nrow(representatives)), function(i) {
  item <- representatives[i]
  cache <- cluster_cache[[item$celltype]]
  target_key <- canonical_gene_set(item$genes)
  primary_hit <- cache$primary[gene_key == target_key]
  joint_hit <- cache$joint[gene_key == target_key]
  max_joint_jaccard <- if (nrow(cache$joint)) {
    max(vapply(cache$joint$gene_list, function(value) {
      jaccard_gene_sets(strsplit(target_key, ",", fixed = TRUE)[[1]], value)
    }, numeric(1)))
  } else {
    0
  }
  data.table(
    label = item$label,
    celltype = item$celltype,
    genes = target_key,
    present_in_primary = nrow(primary_hit) > 0,
    primary_cluster_id = paste(primary_hit$cluster_id, collapse = ";"),
    exact_in_joint_score = nrow(joint_hit) > 0,
    joint_score_cluster_id = paste(joint_hit$cluster_id, collapse = ";"),
    best_joint_score_jaccard = round(max_joint_jaccard, 6)
  )
}))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
summary_file <- file.path(output_dir, "table_s2_joint_score_cluster_sensitivity.tsv")
retention_file <- file.path(output_dir, "representative_cluster_retention.tsv")
fwrite(summary, summary_file, sep = "\t")
fwrite(representative_retention, retention_file, sep = "\t")

message("Wrote ", summary_file)
message("Wrote ", retention_file)
