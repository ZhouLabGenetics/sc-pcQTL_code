#!/usr/bin/env Rscript

# =============================================================================
# Method 2 - Step 3: Identify clusters via sliding window (SC hurdle)
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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript step3_identify_clusters_single_chr.R <CHR>")

CHR <- as.integer(args[1])

cat(sprintf("\n=== Method 2 Step 3: Identify Clusters - Chr%d (%s) ===\n",
            CHR, cfg$cell_type))
cat(sprintf("Start time: %s\n\n", Sys.time()))

GENE_INFO_FILE <- cfg$gene_info_file
ASSOC_DIR <- cfg$assoc_dir
OUTPUT_DIR <- cfg$clusters_dir
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

MAX_WINDOW_SIZE <- 50
MIN_WINDOW_SIZE <- 2
CLUSTER_THRESHOLD <- 0.70

cat(sprintf("Parameters: window %d-%d, threshold %.0f%%\n",
            MIN_WINDOW_SIZE, MAX_WINDOW_SIZE, CLUSTER_THRESHOLD * 100))

sig_file <- file.path(ASSOC_DIR, sprintf("chr%d", CHR), sprintf("chr%d_significant_pairs.tsv", CHR))
if (!file.exists(sig_file)) {
  cat("  No significant pair file; writing empty summary.\n")
  fwrite(data.table(Chromosome = CHR, NumClusters = 0),
         file.path(OUTPUT_DIR, sprintf("chr%d_cluster_summary.tsv", CHR)),
         sep = "\t")
  quit(status = 0)
}

gene_info <- fread(GENE_INFO_FILE)
chr_genes <- gene_info[chr_numeric == CHR][order(start), gene_name]
chr_genes <- unique(chr_genes)
n_genes <- length(chr_genes)

cat(sprintf("  Chromosome %d genes: %d\n", CHR, n_genes))

sig_data <- fread(sig_file)
if (nrow(sig_data) == 0) {
  cat("  No significant pairs; nothing to cluster.\n")
  fwrite(data.table(Chromosome = CHR, NumClusters = 0, NumGenes = n_genes),
         file.path(OUTPUT_DIR, sprintf("chr%d_cluster_summary.tsv", CHR)), sep = "\t")
  quit(status = 0)
}

sig_matrix <- matrix(FALSE, nrow = n_genes, ncol = n_genes)
rownames(sig_matrix) <- chr_genes
colnames(sig_matrix) <- chr_genes
diag(sig_matrix) <- TRUE

for (i in seq_len(nrow(sig_data))) {
  g1 <- sig_data$Gene1[i]
  g2 <- sig_data$Gene2[i]
  if (g1 %in% chr_genes && g2 %in% chr_genes) {
    sig_matrix[g1, g2] <- TRUE
    sig_matrix[g2, g1] <- TRUE
  }
}

check_cluster <- function(window_genes, sig_matrix) {
  n <- length(window_genes)
  if (n < 2) return(FALSE)
  sub_mat <- sig_matrix[window_genes, window_genes]
  pct <- mean(sub_mat[upper.tri(sub_mat)], na.rm = TRUE)
  pct >= CLUSTER_THRESHOLD
}

assigned <- rep(FALSE, n_genes)
names(assigned) <- chr_genes
clusters <- list()
cluster_id <- 1

for (window_size in MAX_WINDOW_SIZE:MIN_WINDOW_SIZE) {
  if (window_size > n_genes) next
  n_windows <- n_genes - window_size + 1
  for (start_idx in 1:n_windows) {
    end_idx <- start_idx + window_size - 1
    window_genes <- chr_genes[start_idx:end_idx]
    if (any(assigned[window_genes])) next
    if (check_cluster(window_genes, sig_matrix)) {
      gene_positions <- gene_info[gene_name %in% window_genes, .(gene_name, start, end)]
      clusters[[cluster_id]] <- data.table(
        cluster_id = sprintf("SC_chr%d_cluster_%03d", CHR, cluster_id),
        chromosome = CHR,
        cluster_size = window_size,
        start_position = min(gene_positions$start),
        end_position = max(gene_positions$end),
        cluster_span_bp = max(gene_positions$end) - min(gene_positions$start),
        genes = paste(window_genes, collapse = ",")
      )
      assigned[window_genes] <- TRUE
      cluster_id <- cluster_id + 1
    }
  }
}

n_clusters <- length(clusters)
cat(sprintf("Chromosome %d: %d clusters\n", CHR, n_clusters))

if (n_clusters > 0) {
  cluster_dt <- rbindlist(clusters)
  fwrite(cluster_dt,
         file.path(OUTPUT_DIR, sprintf("chr%d_clusters.tsv", CHR)),
         sep = "\t")
}

summary_dt <- data.table(
  Chromosome = CHR,
  NumClusters = n_clusters,
  NumGenes = n_genes,
  AssignedGenes = sum(assigned)
)
fwrite(summary_dt,
       file.path(OUTPUT_DIR, sprintf("chr%d_cluster_summary.tsv", CHR)),
       sep = "\t")

cat("\n=== Complete ===\n")
cat(sprintf("End time: %s\n", Sys.time()))
