#!/usr/bin/env Rscript

.local_lib <- Sys.getenv("SC_PCQTL_R_LIBS", unset = "")
if (nzchar(.local_lib)) .libPaths(c(.local_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args_all[grep("--file=", args_all)[1]]))
module_dir <- dirname(script_path)
source(file.path(module_dir, "load_config.R"))
cfg <- load_joint_score_config(module_dir)

cluster_files <- file.path(cfg$clusters_dir, sprintf("chr%d_clusters.tsv", 1:22))
if (any(!file.exists(cluster_files))) {
  stop("Missing chromosome cluster files: ", paste(basename(cluster_files[!file.exists(cluster_files)]), collapse = ", "))
}
all_clusters <- rbindlist(lapply(cluster_files, fread), fill = TRUE)
setorder(all_clusters, chromosome, start_position, end_position, cluster_id)

if (nrow(all_clusters)) {
  gene_assignments <- all_clusters[, .(
    gene_name = unlist(strsplit(genes, ",", fixed = TRUE))
  ), by = .(cluster_id, chromosome, cluster_size)]
} else {
  gene_assignments <- data.table(
    cluster_id = character(), chromosome = integer(), cluster_size = integer(),
    gene_name = character()
  )
}

dir.create(cfg$merged_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(all_clusters, file.path(cfg$merged_dir, "cluster_summary.tsv"), sep = "\t")
fwrite(gene_assignments, file.path(cfg$merged_dir, "cluster_gene_assignments.tsv"), sep = "\t")
writeLines("OK", file.path(cfg$merged_dir, "MERGE_COMPLETE"))
message(sprintf("[%s] Merged %d joint-score clusters", cfg$celltype, nrow(all_clusters)))
