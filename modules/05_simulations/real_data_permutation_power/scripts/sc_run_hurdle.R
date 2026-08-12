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
  make_option("--sc_file", help = "RDS with list(counts, genes, ...)", metavar = "FILE"),
  make_option("--output_dir", default = "results/power_full/sc",
              help = "Directory for SC hurdle outputs"),
  make_option("--nonzero_cutoff", type = "double", default = 0.01,
              help = "Minimum nonzero-cell proportion to run fasthurdle [default %default]"),
  make_option("--p_thresholds", default = "0.001,0.01,0.05"),
  make_option("--max_pairs", type = "integer", default = 0,
              help = "Optional random subset of gene pairs (0 = all)"),
  make_option("--seed_pairs", type = "integer", default = 20240120)
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$sc_file) || !file.exists(opt$sc_file)) {
  stop("Must provide --sc_file")
}
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
p_thresholds <- as.numeric(strsplit(opt$p_thresholds, ",")[[1]])

message("=== Running SC fasthurdle ===")
message(sprintf("fasthurdle version: %s", as.character(utils::packageVersion("fasthurdle"))))
message(sprintf("Input: %s", opt$sc_file))
message(sprintf("Output: %s", opt$output_dir))

sc_obj <- readRDS(opt$sc_file)
counts <- sc_obj$counts
genes <- sc_obj$genes
if (is.null(counts) || is.null(genes)) stop("counts/genes missing in sc_file")

n_genes <- length(genes)
n_cells <- nrow(counts)
nonzero_prop <- colMeans(counts > 0)

pair_idx <- combn(n_genes, 2, simplify = FALSE)
if (opt$max_pairs > 0 && opt$max_pairs < length(pair_idx)) {
  set.seed(opt$seed_pairs)
  pair_idx <- sample(pair_idx, opt$max_pairs)
}
n_pairs <- length(pair_idx)
message(sprintf("Genes: %d  Cells: %d  Pairs queued: %d",
                n_genes, n_cells, n_pairs))

pval_count <- matrix(NA_real_, n_genes, n_genes, dimnames = list(genes, genes))
pval_zero <- matrix(NA_real_, n_genes, n_genes, dimnames = list(genes, genes))
pair_tag <- matrix("filtered", n_genes, n_genes, dimnames = list(genes, genes))
directional_count_dt <- data.table(Gene1 = character(), Gene2 = character(), Pvalue = numeric(), Direction = character())
directional_zero_dt <- data.table(Gene1 = character(), Gene2 = character(), Pvalue = numeric(), Direction = character())

run_pair <- function(i, j) {
  df <- data.frame(
    count_i = counts[, i],
    count_j = counts[, j]
  )
  if (all(df$count_i == 0) || all(df$count_j == 0)) {
    return(list(p_count = NA_real_, p_zero = NA_real_))
  }
  fit_ij <- tryCatch({
    fasthurdle(count_i ~ count_j, data = df,
               dist = "poisson", zero.dist = "binomial")
  }, error = function(e) NULL)
  fit_ji <- tryCatch({
    fasthurdle(count_j ~ count_i, data = df,
               dist = "poisson", zero.dist = "binomial")
  }, error = function(e) NULL)
  if (is.null(fit_ij) || is.null(fit_ji)) {
    return(list(p_count = NA_real_, p_zero = NA_real_))
  }
  sum_ij <- summary(fit_ij)$coefficients
  sum_ji <- summary(fit_ji)$coefficients
  grab <- function(smry, component, term) {
    if (!is.null(smry[[component]]) && term %in% rownames(smry[[component]])) {
      smry[[component]][term, "Pr(>|z|)"]
    } else {
      NA_real_
    }
  }
  p_count_ij <- grab(sum_ij, "count", "count_j")
  p_count_ji <- grab(sum_ji, "count", "count_i")
  p_zero_ij <- grab(sum_ij, "zero", "count_j")
  p_zero_ji <- grab(sum_ji, "zero", "count_i")
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
for (idx in seq_along(pair_idx)) {
  pair <- pair_idx[[idx]]
  i <- pair[1]; j <- pair[2]
  nz_min <- min(nonzero_prop[i], nonzero_prop[j])
  if (nz_min < opt$nonzero_cutoff) {
    filtered_pairs <- filtered_pairs + 1L
    next
  }
  eligible_pairs <- eligible_pairs + 1L
  pair_tag[i, j] <- pair_tag[j, i] <- "tested"
  res <- run_pair(i, j)
  pval_count[i, j] <- pval_count[j, i] <- res$p_count
  pval_zero[i, j] <- pval_zero[j, i] <- res$p_zero
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
  if (idx %% 1000 == 0 || idx == n_pairs) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    rate <- idx / max(elapsed, 1e-6)
    eta <- (n_pairs - idx) / max(rate, 1e-6)
    message(sprintf("  %d/%d pairs processed (%.1f%%)  %.2f pairs/s  ETA %.1f min",
                    idx, n_pairs, 100 * idx / n_pairs, rate, eta / 60))
  }
}

fwrite(as.data.table(pval_count, keep.rownames = "Gene"),
       file.path(opt$output_dir, "pvalues_count_mean.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(as.data.table(pval_zero, keep.rownames = "Gene"),
       file.path(opt$output_dir, "pvalues_zero_mean.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(as.data.table(pair_tag, keep.rownames = "Gene"),
       file.path(opt$output_dir, "pair_filter.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(directional_count_dt,
       file.path(opt$output_dir, "sc_hurdle_count_pvalues.tsv"),
       sep = "\t")
fwrite(directional_zero_dt,
       file.path(opt$output_dir, "sc_hurdle_zero_pvalues.tsv"),
       sep = "\t")

summary_dt <- data.table(
  metric = c("eligible_pairs", "filtered_pairs", "nonzero_cutoff"),
  value = c(eligible_pairs, filtered_pairs, opt$nonzero_cutoff)
)
fwrite(summary_dt, file.path(opt$output_dir, "summary_metrics.tsv"), sep = "\t")

sink(file.path(opt$output_dir, "summary_report.txt"))
cat("=== SC fasthurdle summary ===\n")
cat(sprintf("Input file: %s\n", opt$sc_file))
cat(sprintf("Genes: %d  Cells: %d\n", n_genes, n_cells))
cat(sprintf("Pairs processed: %d (eligible %d, filtered %d)\n",
            n_pairs, eligible_pairs, filtered_pairs))
cat(sprintf("Nonzero cutoff: %.3f\n\n", opt$nonzero_cutoff))
sink()

message("SC hurdle outputs written to ", opt$output_dir)
