#!/usr/bin/env Rscript

# =============================================================================
# Method 2 - Step 4: Merge SC hurdle clusters across chromosomes
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

get_this_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_idx <- grep("--file=", cmd_args)
  if (length(file_idx) > 0) {
    return(normalizePath(sub("--file=", "", cmd_args[file_idx])))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
  normalizePath(".")
}

script_dir <- dirname(get_this_path())
source(file.path(script_dir, "load_config.R"))
cfg <- load_config(script_dir)

cat(sprintf("\n=== Method 2 Step 4: Merge SC Clusters (%s) ===\n", cfg$cell_type))
cat(sprintf("Start time: %s\n\n", Sys.time()))

CLUSTER_DIR <- cfg$clusters_dir
OUTPUT_DIR <- file.path(cfg$method2_results, "merged_clusters")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cluster_files <- list.files(CLUSTER_DIR, pattern = "^chr\\d+_clusters\\.tsv$", full.names = TRUE)
cat(sprintf("Found %d chromosome cluster files\n", length(cluster_files)))

if (!length(cluster_files)) stop("No cluster files detected")

all_clusters <- rbindlist(lapply(cluster_files, fread), fill = TRUE)
cat(sprintf("Total clusters: %d\n", nrow(all_clusters)))

gene_assignments <- all_clusters[, .(
  gene_name = unlist(strsplit(genes, ","))
), by = .(cluster_id, chromosome, cluster_size)]

fwrite(all_clusters,
       file.path(OUTPUT_DIR, "cluster_summary.tsv"),
       sep = "\t")
fwrite(gene_assignments,
       file.path(OUTPUT_DIR, "cluster_gene_assignments.tsv"),
       sep = "\t")

cat("\nCluster counts per chromosome:\n")
print(all_clusters[, .N, by = chromosome][order(chromosome)])

cat("\nCluster size distribution:\n")
print(all_clusters[, .N, by = cluster_size][order(cluster_size)])

cat(sprintf("\nOutputs written to: %s\n", OUTPUT_DIR))
cat(sprintf("End time: %s\n", Sys.time()))
