#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(fasthurdle)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
sim_code_root <- Sys.getenv(
  "COQTL_SIM_CODE_ROOT",
  unset = normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
)
source(file.path(sim_code_root, "shared", "sim_paths.R"))

option_list <- list(
  make_option("--sc_null_file",
              default = file.path(get_data_root(), "shuffle_full", "full_realcounts", "sc_counts_shuffle_null.rds")),
  make_option("--output_dir",
              default = file.path(get_sim_root(), "01_main_text_results", "01_null_calibration_final_sc", "results", "full_realcounts", "sc_hurdle")),
  make_option("--p_thresholds", default = "0.05,0.01,0.005,0.001"),
  make_option("--nonzero_cutoff", type = "double", default = 0.01,
              help = "Minimum nonzero-cell proportion to treat pair as dense [default %default]"),
  make_option("--n_pairs", type = "integer", default = 0,
              help = "Number of gene pairs to test (0 = all pairs)"),
  make_option("--seed_pairs", type = "integer", default = 20240113)
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
P_THRESHOLDS <- as.numeric(strsplit(opt$p_thresholds, ",")[[1]])

message("=== Testing unadjusted SC hurdle models on the per-gene permutation null ===")
message(sprintf("fasthurdle version: %s", as.character(utils::packageVersion("fasthurdle"))))
sc_null <- readRDS(opt$sc_null_file)
count_matrix <- sc_null$counts
genes <- sc_null$genes
n_genes <- length(genes)
n_cells <- nrow(count_matrix)

nonzero_prop <- colMeans(count_matrix > 0)
message(sprintf("  Cells: %d, Genes: %d", n_cells, n_genes))

pair_indices <- combn(n_genes, 2, simplify = FALSE)
if (opt$n_pairs > 0 && opt$n_pairs < length(pair_indices)) {
  set.seed(opt$seed_pairs)
  pair_indices <- sample(pair_indices, opt$n_pairs)
}
n_pairs <- length(pair_indices)
message(sprintf("  Pairs to test: %d", n_pairs))

pval_count_mat <- matrix(NA_real_, n_genes, n_genes, dimnames = list(genes, genes))
pval_zero_mat <- matrix(NA_real_, n_genes, n_genes, dimnames = list(genes, genes))
pair_filter <- matrix("filtered", n_genes, n_genes, dimnames = list(genes, genes))
directional_count_dt <- data.table(Gene1 = character(), Gene2 = character(), Pvalue = numeric(), Direction = character())
directional_zero_dt <- data.table(Gene1 = character(), Gene2 = character(), Pvalue = numeric(), Direction = character())

hurdle_pair <- function(i, j) {
  df <- data.frame(
    count_i = count_matrix[, i],
    count_j = count_matrix[, j]
  )
  if (min(df$count_i) > 0 || min(df$count_j) > 0) {
    return(list(
      p_count = NA_real_, p_zero = NA_real_,
      p_count_ij = NA_real_, p_count_ji = NA_real_,
      p_zero_ij = NA_real_, p_zero_ji = NA_real_
    ))
  }
  model_ij <- tryCatch(
    fasthurdle(count_i ~ count_j, data = df, dist = "poisson", zero.dist = "binomial"),
    error = function(e) NULL
  )
  model_ji <- tryCatch(
    fasthurdle(count_j ~ count_i, data = df, dist = "poisson", zero.dist = "binomial"),
    error = function(e) NULL
  )
  if (is.null(model_ij) || is.null(model_ji)) {
    return(list(
      p_count = NA_real_, p_zero = NA_real_,
      p_count_ij = NA_real_, p_count_ji = NA_real_,
      p_zero_ij = NA_real_, p_zero_ji = NA_real_
    ))
  }
  sum_ij <- summary(model_ij)$coefficients
  sum_ji <- summary(model_ji)$coefficients
  get_p <- function(smry, term) {
    if (!is.null(smry) && term %in% rownames(smry)) smry[term, "Pr(>|z|)"] else NA_real_
  }
  p_count_ij <- get_p(sum_ij$count, "count_j")
  p_count_ji <- get_p(sum_ji$count, "count_i")
  p_zero_ij <- get_p(sum_ij$zero, "count_j")
  p_zero_ji <- get_p(sum_ji$zero, "count_i")
  list(
    p_count = mean(c(p_count_ij, p_count_ji), na.rm = TRUE),
    p_zero = mean(c(p_zero_ij, p_zero_ji), na.rm = TRUE),
    p_count_ij = p_count_ij,
    p_count_ji = p_count_ji,
    p_zero_ij = p_zero_ij,
    p_zero_ji = p_zero_ji
  )
}

eligible_pairs <- 0L
filtered_pairs <- 0L
start_time <- Sys.time()
for (idx in seq_along(pair_indices)) {
  pair <- pair_indices[[idx]]
  i <- pair[1]
  j <- pair[2]
  if (min(nonzero_prop[i], nonzero_prop[j]) < opt$nonzero_cutoff) {
    filtered_pairs <- filtered_pairs + 1L
    next
  }
  eligible_pairs <- eligible_pairs + 1L
  pair_filter[i, j] <- pair_filter[j, i] <- "tested"

  res <- hurdle_pair(i, j)

  pval_count_mat[i, j] <- pval_count_mat[j, i] <- res$p_count
  pval_zero_mat[i, j] <- pval_zero_mat[j, i] <- res$p_zero

  directional_count_dt <- rbind(
    directional_count_dt,
    data.table(
      Gene1 = c(genes[i], genes[j]),
      Gene2 = c(genes[j], genes[i]),
      Pvalue = c(res$p_count_ij, res$p_count_ji),
      Direction = c("i_to_j", "j_to_i")
    ),
    use.names = TRUE
  )
  directional_zero_dt <- rbind(
    directional_zero_dt,
    data.table(
      Gene1 = c(genes[i], genes[j]),
      Gene2 = c(genes[j], genes[i]),
      Pvalue = c(res$p_zero_ij, res$p_zero_ji),
      Direction = c("i_to_j", "j_to_i")
    ),
    use.names = TRUE
  )

  if (idx %% 250 == 0 || idx == n_pairs) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    rate <- idx / max(elapsed, 1e-6)
    eta <- (n_pairs - idx) / max(rate, 1e-6)
    message(sprintf("  %d/%d pairs (%.1f%%)  %.2f pairs/sec, ETA %.1f min",
                    idx, n_pairs, 100 * idx / n_pairs, rate, eta / 60))
  }
}

calc_fpr <- function(vals, thresh) mean(vals < thresh, na.rm = TRUE)
results_dt <- data.table()
for (pth in P_THRESHOLDS) {
  count_flags <- directional_count_dt$Pvalue < pth
  zero_flags <- directional_zero_dt$Pvalue < pth
  combined_flags <- count_flags | zero_flags
  count_flags[is.na(count_flags)] <- FALSE
  zero_flags[is.na(zero_flags)] <- FALSE
  combined_flags[is.na(combined_flags)] <- FALSE
  results_dt <- rbind(results_dt, data.table(
    p_threshold = pth,
    count_fpr = mean(count_flags),
    zero_fpr = mean(zero_flags),
    combined_fpr = mean(combined_flags),
    eligible_pairs = eligible_pairs,
    filtered_pairs = filtered_pairs
  ))
}

fwrite(as.data.table(pval_count_mat, keep.rownames = "Gene"),
       file.path(opt$output_dir, "pvalues_count_mean.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(as.data.table(pval_zero_mat, keep.rownames = "Gene"),
       file.path(opt$output_dir, "pvalues_zero_mean.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(as.data.table(pair_filter, keep.rownames = "Gene"),
       file.path(opt$output_dir, "pair_filter.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(data.table(Gene = genes, nonzero_prop = nonzero_prop),
       file.path(opt$output_dir, "gene_nonzero_prop.tsv"), sep = "\t")
fwrite(directional_count_dt,
       file.path(opt$output_dir, "sc_hurdle_count_pvalues.tsv"), sep = "\t")
fwrite(directional_zero_dt,
       file.path(opt$output_dir, "sc_hurdle_zero_pvalues.tsv"), sep = "\t")
fwrite(results_dt, file.path(opt$output_dir, "fpr_summary.tsv"), sep = "\t")
fwrite(data.table(
  metric = c("eligible_pairs", "filtered_pairs", "nonzero_cutoff"),
  value = c(eligible_pairs, filtered_pairs, opt$nonzero_cutoff)
), file.path(opt$output_dir, "summary_metrics.tsv"), sep = "\t")

sink(file.path(opt$output_dir, "summary_report.txt"))
cat("=== SC hurdle per-gene permutation-null summary ===\n")
cat(sprintf("Data file: %s\n", opt$sc_null_file))
cat(sprintf("Pairs tested: %d\n", n_pairs))
cat(sprintf("Eligible pairs (>=%.2f%% nonzero): %d\n", 100 * opt$nonzero_cutoff, eligible_pairs))
cat(sprintf("Filtered pairs: %d\n\n", filtered_pairs))
cat("Covariates and library-size terms: none\n")
cat("Directional rule: both ordered fits retained as separate tests\n\n")
cat("FPR results:\n")
print(results_dt)
sink()

message("Done. Results saved to ", opt$output_dir)
