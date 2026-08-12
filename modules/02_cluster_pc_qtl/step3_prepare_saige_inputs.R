#!/usr/bin/env Rscript
# ============================================================================
# Step 3 preparation: region files + cluster_pc_map for SAIGE-QTL
# ============================================================================

suppressPackageStartupMessages({ library(data.table) })
cat(sprintf("\n=== Prepare SAIGE Inputs ===\nStart: %s\n\n", Sys.time()))

source("config.R")

STEP2_DIR  <- file.path(PCQTL_DIR, "step2_pca")
SAIGE_DIR  <- file.path(PCQTL_DIR, "step3_saige")
REGION_DIR <- file.path(SAIGE_DIR, "regions")
dir.create(REGION_DIR, showWarnings = FALSE, recursive = TRUE)

CIS_WINDOW <- 500000   # 500 kb each side

clusters <- fread(CLUSTER_FILE)
pc_list  <- fread(file.path(STEP2_DIR, "pc_list_for_saige.tsv"))

cat(sprintf("Clusters in summary : %d\n", nrow(clusters)))
cat(sprintf("Cluster-PCs (95%%)  : %d\n\n", nrow(pc_list)))

# keep only clusters that have PCA output
clusters <- clusters[cluster_id %in% pc_list$cluster_id]
cat(sprintf("Clusters with PCA output : %d\n\n", nrow(clusters)))

# ── region files (one per cluster, no header) ───────────────────────────────
cat("Writing region files ...\n")
for (i in seq_len(nrow(clusters))) {
  cid   <- clusters$cluster_id[i]
  chr   <- clusters$chromosome[i]
  start <- max(0L, as.integer(clusters$start_position[i]) - CIS_WINDOW)
  end   <- as.integer(clusters$end_position[i])   + CIS_WINDOW
  writeLines(paste(chr, start, end, sep = "\t"),
             file.path(REGION_DIR, paste0(cid, "_region.txt")))
}

# ── cluster_pc_map.tsv ──────────────────────────────────────────────────────
# columns: cluster_id  PC  pheno_file  chromosome  region_file
pc_map <- merge(pc_list, clusters[, .(cluster_id, chromosome)], by = "cluster_id")
pc_map[, region_file := file.path(REGION_DIR, paste0(cluster_id, "_region.txt"))]
fwrite(pc_map, file.path(SAIGE_DIR, "cluster_pc_map.tsv"), sep = "\t", quote = FALSE)

cat(sprintf("Region files written : %d\n", nrow(clusters)))
cat(sprintf("cluster_pc_map rows  : %d\n", nrow(pc_map)))
cat(sprintf("Output dir           : %s\n", SAIGE_DIR))
cat(sprintf("End                  : %s\n", Sys.time()))
