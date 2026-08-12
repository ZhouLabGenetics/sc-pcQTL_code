#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--sim_dir", type = "character"),
  make_option("--out_dir", type = "character"),
  make_option("--template_dir", type = "character",
              default = Sys.getenv(
                "SC_PCQTL_HURDLE_SCRIPT_DIR",
                unset = normalizePath(file.path(getwd(), "..", "..", "01_pairwise_hurdle_clustering"), mustWork = FALSE)
              )),
  make_option("--cell_type", type = "character", default = "simulation"),
  make_option("--p_threshold", type = "double", default = 0.05),
  make_option("--nonzero_cutoff", type = "double", default = 0.01)
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$sim_dir) || is.null(opt$out_dir)) {
  stop("--sim_dir and --out_dir are required.")
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$out_dir, "results"), recursive = TRUE, showWarnings = FALSE)

count_file <- file.path(opt$sim_dir, "sim_counts.tsv.gz")
gene_info_file <- file.path(opt$sim_dir, "gene_info.tsv")
if (!file.exists(count_file) || !file.exists(gene_info_file)) {
  stop("Missing simulation inputs in ", opt$sim_dir)
}

config_path <- file.path(opt$out_dir, "config.R")
config_lines <- c(
  sprintf('CELL_TYPE <- "%s"', opt$cell_type),
  sprintf('DATA_ROOT <- "%s"', normalizePath(opt$sim_dir)),
  sprintf('WORK_ROOT <- "%s"', normalizePath(opt$out_dir)),
  sprintf('WORK_DIR <- "%s"', normalizePath(opt$out_dir)),
  sprintf('COUNT_FILE <- "%s"', normalizePath(count_file)),
  sprintf('GENE_INFO_FILE <- "%s"', normalizePath(gene_info_file)),
  sprintf('RESULTS_ROOT <- "%s"', normalizePath(file.path(opt$out_dir, "results"))),
  sprintf('METHOD2_RESULTS <- "%s"', normalizePath(file.path(opt$out_dir, "results", "method2_sc_hurdle"), mustWork = FALSE)),
  sprintf('LOG_DIR <- "%s"', normalizePath(file.path(opt$out_dir, "logs"), mustWork = FALSE))
)
writeLines(config_lines, config_path)

template_files <- c(
  "load_config.R",
  "step1_filter_sparse_genes.R",
  "step2_calculate_sc_associations_chunked.R"
)
for (f in template_files) {
  ok <- file.copy(file.path(opt$template_dir, f), file.path(opt$out_dir, f), overwrite = TRUE)
  if (!ok) {
    stop("Failed to copy template file: ", f)
  }
}

step1_script <- file.path(opt$out_dir, "step1_filter_sparse_genes.R")
step2_script <- file.path(opt$out_dir, "step2_calculate_sc_associations_chunked.R")

# Simulation calibration keeps both directions as separate tests.
step2_lines <- readLines(step2_script)
step2_lines <- sub("global_pairs <- global_counts\\[, sum\\(N \\* \\(N - 1\\) / 2\\)\\]",
                   "global_pairs <- global_counts[, sum(N * (N - 1))]",
                   step2_lines)
step2_lines <- sub("n_pairs <- length\\(pair_idx\\)",
                   "n_pair_jobs <- length(pair_idx)\nn_pairs <- n_genes * (n_genes - 1)",
                   step2_lines)
step2_lines <- sub("idx == n_pairs", "idx == n_pair_jobs", step2_lines, fixed = TRUE)
step2_lines <- sub("n_pairs - idx", "n_pair_jobs - idx", step2_lines, fixed = TRUE)
step2_lines <- sub("100 * idx / n_pairs", "100 * idx / n_pair_jobs", step2_lines, fixed = TRUE)
step2_lines <- sub("  p_count <- min\\(p_count_vals, na.rm = TRUE\\)",
                   "  p_count <- fit_ij$p_count",
                   step2_lines)
step2_lines <- sub("  if \\(!is.finite\\(p_count\\)\\) p_count <- NA_real_",
                   "  p_count_ji <- fit_ji$p_count\n  if (!is.finite(p_count)) p_count <- NA_real_\n  if (!is.finite(p_count_ji)) p_count_ji <- NA_real_",
                   step2_lines)
step2_lines <- sub("  if \\(!is.finite\\(p_zero\\)\\) p_zero <- NA_real_",
                   "  p_zero_ji <- fit_ji$p_zero\n  if (!is.finite(p_zero)) p_zero <- NA_real_\n  if (!is.finite(p_zero_ji)) p_zero_ji <- NA_real_",
                   step2_lines)
step2_lines <- sub("  p_zero <- min\\(p_zero_vals, na.rm = TRUE\\)",
                   "  p_zero <- fit_ij$p_zero",
                   step2_lines)
step2_lines <- sub("  list\\(p_count = p_count, p_zero = p_zero\\)",
                   "  list(p_count = p_count, p_zero = p_zero, p_count_ji = p_count_ji, p_zero_ji = p_zero_ji)",
                   step2_lines)
step2_lines <- sub("  pval_count\\[i, j\\] <- pval_count\\[j, i\\] <- res\\$p_count",
                   "  pval_count[i, j] <- res$p_count\n  pval_count[j, i] <- res$p_count_ji",
                   step2_lines)
step2_lines <- sub("  pval_zero\\[i, j\\] <- pval_zero\\[j, i\\] <- res\\$p_zero",
                   "  pval_zero[i, j] <- res$p_zero\n  pval_zero[j, i] <- res$p_zero_ji",
                   step2_lines)

sig_start <- grep("^sig_pairs <- data.table\\(\\)$", step2_lines)
save_start <- grep("^# Save results$", step2_lines)
if (length(sig_start) != 1 || length(save_start) != 1 || sig_start >= save_start) {
  stop("Upstream hurdle worker layout changed; cannot construct directional simulation worker")
}
step2_lines <- c(
  step2_lines[seq_len(sig_start - 1)],
  "sig_pairs <- data.table()",
  "for (i in seq_len(n_genes)) {",
  "  for (j in seq_len(n_genes)) {",
  "    if (i == j) next",
  "    p_count <- pval_count[i, j]",
  "    p_zero <- pval_zero[i, j]",
  "    if ((!is.na(p_count) && p_count < bonferroni_threshold) ||",
  "        (!is.na(p_zero) && p_zero < bonferroni_threshold)) {",
  "      sig_pairs <- rbind(sig_pairs, data.table(",
  "        Gene1 = chr_genes[i],",
  "        Gene2 = chr_genes[j],",
  "        Pvalue_count = p_count,",
  "        Pvalue_zero = p_zero",
  "      ))",
  "    }",
  "  }",
  "}",
  step2_lines[save_start:length(step2_lines)]
)

required_directional_lines <- c(
  "global_pairs <- global_counts[, sum(N * (N - 1))]",
  "n_pair_jobs <- length(pair_idx)",
  "n_pairs <- n_genes * (n_genes - 1)",
  "p_count <- fit_ij$p_count",
  "p_count_ji <- fit_ji$p_count",
  "p_zero <- fit_ij$p_zero",
  "p_zero_ji <- fit_ji$p_zero",
  "pval_count[j, i] <- res$p_count_ji",
  "pval_zero[j, i] <- res$p_zero_ji"
)
missing_directional_lines <- required_directional_lines[!vapply(
  required_directional_lines,
  function(line) any(grepl(line, step2_lines, fixed = TRUE)),
  logical(1)
)]
if (length(missing_directional_lines)) {
  stop(
    "Failed to adapt upstream hurdle worker for directional simulation tests: ",
    paste(missing_directional_lines, collapse = "; ")
  )
}
writeLines(step2_lines, step2_script)

rscript_bin <- Sys.getenv("COQTL_RSCRIPT", unset = Sys.getenv("R_SCRIPT", unset = "Rscript"))
run_rscript <- function(args, log_file) {
  output <- system2(
    rscript_bin, args, stdout = TRUE, stderr = TRUE,
    env = "R_PROFILE_USER="
  )
  writeLines(output, log_file)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Child R script failed with status ", status, "; see ", log_file)
  }
}

run_rscript(step1_script, file.path(opt$out_dir, "logs", "step1.log"))

gene_info <- fread(gene_info_file)
n_genes_chr1 <- nrow(gene_info[chr_numeric == 1])
if (n_genes_chr1 < 2) {
  stop("Need at least 2 genes on chr1.")
}

cmd2 <- c(
  step2_script,
  "--chr", "1",
  "--gene_start", "1",
  "--gene_end", as.character(n_genes_chr1),
  "--chunk_id", "1",
  "--p_threshold", as.character(opt$p_threshold),
  "--nonzero_cutoff", as.character(opt$nonzero_cutoff)
)
run_rscript(cmd2, file.path(opt$out_dir, "logs", "step2.log"))

message("Finished hurdle analysis in: ", opt$out_dir)
