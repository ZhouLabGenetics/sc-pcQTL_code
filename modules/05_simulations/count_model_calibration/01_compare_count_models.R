#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fasthurdle)
})

required_fasthurdle <- base::package_version("1.2.0")
if (utils::packageVersion("fasthurdle") != required_fasthurdle) {
  stop(
    "This diagnostic requires fasthurdle 1.2.0; found ",
    as.character(utils::packageVersion("fasthurdle")),
    call. = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 4L) {
  stop(
    paste(
      "Usage: Rscript 01_compare_count_models.R",
      "NULL_RDS OUTPUT_DIR [NONZERO_CUTOFF] [SETTINGS_TSV]"
    ),
    call. = FALSE
  )
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
null_rds <- args[[1L]]
output_dir <- args[[2L]]
nonzero_cutoff <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 0.01
settings_file <- if (length(args) == 4L) args[[4L]] else file.path(script_dir, "settings.tsv")

if (!file.exists(null_rds)) stop("Missing permutation-null RDS: ", null_rds)
if (!file.exists(settings_file)) stop("Missing settings table: ", settings_file)
if (!is.finite(nonzero_cutoff) || nonzero_cutoff < 0 || nonzero_cutoff > 1) {
  stop("NONZERO_CUTOFF must be between 0 and 1")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

settings <- fread(settings_file)
required_columns <- c("setting_label", "count_dist", "zero_dist", "score_test")
if (!all(required_columns %in% names(settings))) {
  stop("SETTINGS_TSV must contain: ", paste(required_columns, collapse = ", "))
}

null_data <- readRDS(null_rds)
if (is.null(null_data$counts)) stop("NULL_RDS must contain a cells-by-genes 'counts' matrix")
counts <- as.matrix(null_data$counts)
storage.mode(counts) <- "numeric"
genes <- null_data$genes
if (is.null(genes)) genes <- colnames(counts)
if (is.null(genes) || length(genes) != ncol(counts)) {
  stop("NULL_RDS must contain one gene label per count-matrix column")
}

nonzero_fraction <- colMeans(counts > 0)
eligible_genes <- which(nonzero_fraction >= nonzero_cutoff)
if (length(eligible_genes) < 2L) stop("Fewer than two genes pass NONZERO_CUTOFF")
eligible_pairs <- utils::combn(eligible_genes, 2L, simplify = FALSE)
p_thresholds <- c(0.05, 0.01, 0.005, 0.001)

as_flag <- function(value) tolower(as.character(value)) %in% c("true", "1", "yes")

extract_wald_p <- function(coefficients, component, predictor) {
  table <- coefficients[[component]]
  if (is.null(table) || !predictor %in% rownames(table)) return(NA_real_)
  as.numeric(table[predictor, "Pr(>|z|)"])
}

fit_direction <- function(data, response, predictor, count_dist, zero_dist, score_enabled) {
  fit <- tryCatch(
    fasthurdle(
      stats::reformulate(predictor, response = response),
      data = data,
      dist = count_dist,
      zero.dist = zero_dist,
      score_test = if (score_enabled) predictor else NULL
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(c(count = NA_real_, zero = NA_real_))

  coefficients <- summary(fit)$coefficients
  count_p <- if (
    score_enabled && !is.null(fit$score_test) && !is.null(fit$score_test$pvalue)
  ) {
    as.numeric(fit$score_test$pvalue)
  } else {
    extract_wald_p(coefficients, "count", predictor)
  }
  c(
    count = count_p,
    zero = extract_wald_p(coefficients, "zero", predictor)
  )
}

run_setting <- function(setting) {
  setting_label <- setting$setting_label[[1L]]
  count_dist <- setting$count_dist[[1L]]
  zero_dist <- setting$zero_dist[[1L]]
  score_enabled <- as_flag(setting$score_test[[1L]])
  setting_dir <- file.path(output_dir, setting_label)
  dir.create(setting_dir, recursive = TRUE, showWarnings = FALSE)

  n_records <- 2L * length(eligible_pairs)
  gene1 <- character(n_records)
  gene2 <- character(n_records)
  direction <- character(n_records)
  count_p <- rep(NA_real_, n_records)
  zero_p <- rep(NA_real_, n_records)

  record <- 0L
  for (pair_index in seq_along(eligible_pairs)) {
    pair <- eligible_pairs[[pair_index]]
    i <- pair[[1L]]
    j <- pair[[2L]]
    pair_data <- data.frame(count_i = counts[, i], count_j = counts[, j])
    result_ij <- fit_direction(
      pair_data, "count_i", "count_j", count_dist, zero_dist, score_enabled
    )
    result_ji <- fit_direction(
      pair_data, "count_j", "count_i", count_dist, zero_dist, score_enabled
    )

    rows <- record + seq_len(2L)
    gene1[rows] <- c(genes[[i]], genes[[j]])
    gene2[rows] <- c(genes[[j]], genes[[i]])
    direction[rows] <- c("i_to_j", "j_to_i")
    count_p[rows] <- c(result_ij[["count"]], result_ji[["count"]])
    zero_p[rows] <- c(result_ij[["zero"]], result_ji[["zero"]])
    record <- record + 2L

    if (pair_index %% 250L == 0L || pair_index == length(eligible_pairs)) {
      message(sprintf("%s: %d/%d pairs", setting_label, pair_index, length(eligible_pairs)))
    }
  }

  count_results <- data.table(
    Gene1 = gene1, Gene2 = gene2, Pvalue = count_p, Direction = direction
  )
  zero_results <- data.table(
    Gene1 = gene1, Gene2 = gene2, Pvalue = zero_p, Direction = direction
  )
  fwrite(count_results, file.path(setting_dir, "sc_hurdle_count_pvalues.tsv"), sep = "\t")
  fwrite(zero_results, file.path(setting_dir, "sc_hurdle_zero_pvalues.tsv"), sep = "\t")

  hurdle_p <- pmin(count_p, zero_p)
  fpr <- rbindlist(lapply(
    list(Count = count_p, Zero = zero_p, Hurdle = hurdle_p),
    function(values) {
      data.table(
        p_threshold = p_thresholds,
        fpr = vapply(
          p_thresholds,
          function(threshold) mean(values < threshold, na.rm = TRUE),
          numeric(1)
        )
      )
    },
    idcol = "method"
  ))
  fwrite(fpr, file.path(setting_dir, "fpr_summary.tsv"), sep = "\t")
  fwrite(
    data.table(
      metric = c(
        "fasthurdle_version", "count_dist", "zero_dist", "score_test",
        "nonzero_cutoff", "n_cells", "n_input_genes", "n_eligible_genes",
        "n_unordered_pairs"
      ),
      value = c(
        as.character(utils::packageVersion("fasthurdle")), count_dist, zero_dist,
        score_enabled, nonzero_cutoff, nrow(counts), ncol(counts),
        length(eligible_genes), length(eligible_pairs)
      )
    ),
    file.path(setting_dir, "summary_metrics.tsv"),
    sep = "\t"
  )
  fpr[, setting_label := setting_label]
  fpr
}

fpr_summary <- rbindlist(lapply(seq_len(nrow(settings)), function(i) {
  run_setting(settings[i])
}))
fwrite(fpr_summary, file.path(output_dir, "fpr_all_settings.tsv"), sep = "\t")
fwrite(settings, file.path(output_dir, "settings_used.tsv"), sep = "\t")

message("Count-model calibration results written to ", normalizePath(output_dir))
