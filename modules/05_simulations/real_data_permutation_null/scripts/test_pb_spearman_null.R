#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
sim_code_root <- Sys.getenv(
  "COQTL_SIM_CODE_ROOT",
  unset = normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
)
source(file.path(sim_code_root, "shared", "sim_paths.R"))

option_list <- list(
  make_option("--pb_null_file",
              default = file.path(get_data_root(), "shuffle_full", "full_realcounts", "pb_counts_shuffle_null.rds")),
  make_option("--output_dir",
              default = file.path(get_sim_root(), "01_main_text_results", "01_null_calibration_final_sc", "results", "full_realcounts", "pb_spearman")),
  make_option("--p_thresholds", default = "0.05,0.01,0.005,0.001")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
P_THRESHOLDS <- as.numeric(strsplit(opt$p_thresholds, ",")[[1]])

cat("\n=== Testing pseudobulk Spearman on the per-gene permutation null ===\n")
cat(sprintf("Start time: %s\n\n", Sys.time()))

pb_null <- readRDS(opt$pb_null_file)
count_matrix <- pb_null$counts
genes <- pb_null$genes
donors <- pb_null$donors

n_genes <- length(genes)
n_donors <- length(donors)
n_pairs <- n_genes * (n_genes - 1) / 2

cat(sprintf("  Donors: %d\n  Genes: %d\n  Pairs: %d\n\n", n_donors, n_genes, n_pairs))

corr_start <- Sys.time()
corr_matrix <- cor(count_matrix, method = "spearman", use = "pairwise.complete.obs")
corr_end <- Sys.time()
t_stat <- corr_matrix * sqrt((n_donors - 2) / (1 - corr_matrix^2 + 1e-10))
pval_mat <- 2 * pt(-abs(t_stat), df = n_donors - 2)
diag(pval_mat) <- 0
upper_vals <- pval_mat[upper.tri(pval_mat)]

fpr_results <- data.table()
for (pth in P_THRESHOLDS) {
  sig_count <- sum(upper_vals < pth, na.rm = TRUE)
  fpr <- sig_count / n_pairs
  fpr_results <- rbind(fpr_results, data.table(
    p_threshold = pth,
    n_pairs_tested = n_pairs,
    n_significant = sig_count,
    fpr = fpr
  ))
  cat(sprintf("  p < %.3f -> FPR %.4f (expected %.3f)\n", pth, fpr, pth))
}

pval_dist <- data.table(
  quantile_level = seq(0, 1, 0.1),
  observed = as.numeric(quantile(upper_vals, probs = seq(0, 1, 0.1), na.rm = TRUE)),
  expected_uniform = seq(0, 1, 0.1)
)
pval_dist[, difference := observed - expected_uniform]

fwrite(as.data.table(corr_matrix, keep.rownames = "Gene"),
       file.path(opt$output_dir, "correlation_matrix.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(as.data.table(pval_mat, keep.rownames = "Gene"),
       file.path(opt$output_dir, "pvalue_matrix.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(fpr_results, file.path(opt$output_dir, "fpr_summary.tsv"), sep = "\t")
fwrite(pval_dist, file.path(opt$output_dir, "pvalue_distribution.tsv"), sep = "\t")

sink(file.path(opt$output_dir, "summary_report.txt"))
cat("=== PB Spearman per-gene permutation-null test ===\n\n")
cat(sprintf("Data: %s\n", opt$pb_null_file))
cat(sprintf("Donors: %d  Genes: %d  Unordered pairs: %d\n", n_donors, n_genes, n_pairs))
cat(sprintf("Runtime: %.2f seconds\n\n",
            as.numeric(difftime(corr_end, corr_start, units = "secs"))))
cat("False positive rates:\n")
print(fpr_results)
cat("\nP-value distribution (observed vs expected):\n")
print(pval_dist)
sink()

cat("\n=== Analysis Complete ===\n")
cat(sprintf("End time: %s\n", Sys.time()))
cat(sprintf("Results saved to: %s\n", opt$output_dir))
