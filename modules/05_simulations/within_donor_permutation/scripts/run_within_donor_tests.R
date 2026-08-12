#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(script_file), "common.R"))
prepend_project_libraries()
suppressPackageStartupMessages({
  library(fasthurdle)
})

args <- parse_cli()
opt <- list(
  input_file = require_arg(args, "input_file"),
  output_dir = require_arg(args, "output_dir"),
  workers = integer_arg(args, "workers", 12L),
  nonzero_cutoff = numeric_arg(args, "nonzero_cutoff", 0.01),
  fasthurdle_version = if (is.null(args$fasthurdle_version)) "1.1.1" else args$fasthurdle_version,
  fastglm_version = if (is.null(args$fastglm_version)) "0.0.4" else args$fastglm_version,
  p_thresholds = if (is.null(args$p_thresholds)) {
    "0.05,0.01,0.005,0.001"
  } else {
    args$p_thresholds
  }
)
if (opt$workers < 1L) stop("--workers must be positive")
loaded_fasthurdle <- as.character(packageVersion("fasthurdle"))
loaded_fastglm <- as.character(packageVersion("fastglm"))
if (!identical(loaded_fasthurdle, opt$fasthurdle_version) ||
    !identical(loaded_fastglm, opt$fastglm_version)) {
  stop(
    "Published diagnostic requires fasthurdle ", opt$fasthurdle_version,
    " and fastglm ", opt$fastglm_version, "; loaded ", loaded_fasthurdle,
    " and ", loaded_fastglm
  )
}
if (dir.exists(opt$output_dir)) stop("Refusing to overwrite result directory: ", opt$output_dir)
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

permuted_data <- readRDS(opt$input_file)
count_matrix <- permuted_data$counts
genes <- as.character(permuted_data$genes)
donors <- as.character(permuted_data$donors)
if (nrow(count_matrix) != length(donors) || ncol(count_matrix) != length(genes)) {
  stop("Input dimensions do not match donor and gene metadata")
}

pair_matrix <- t(combn(seq_along(genes), 2L))
nonzero_prop <- colMeans(count_matrix > 0)
workers <- min(opt$workers, nrow(pair_matrix))
worker_indices <- split(
  seq_len(nrow(pair_matrix)),
  ((seq_len(nrow(pair_matrix)) - 1L) %% workers) + 1L
)

fit_direction <- function(response, predictor) {
  model_data <- data.frame(response = response, predictor = predictor)
  fit <- tryCatch(
    fasthurdle(
      response ~ predictor,
      data = model_data,
      dist = "poisson",
      zero.dist = "binomial"
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(ok = FALSE, p = c(count = NA_real_, detection = NA_real_)))
  }
  coefs <- summary(fit)$coefficients
  list(
    ok = TRUE,
    p = c(
      count = extract_term_p(coefs$count, "predictor"),
      detection = extract_term_p(coefs$zero, "predictor")
    )
  )
}

run_worker <- function(indices) {
  records <- vector("list", length(indices))
  for (k in seq_along(indices)) {
    pair_index <- indices[[k]]
    i <- pair_matrix[pair_index, 1L]
    j <- pair_matrix[pair_index, 2L]
    eligible <- min(nonzero_prop[c(i, j)]) >= opt$nonzero_cutoff
    p_ij <- c(count = NA_real_, detection = NA_real_)
    p_ji <- c(count = NA_real_, detection = NA_real_)
    status <- "filtered"
    if (eligible) {
      if (min(count_matrix[, i]) > 0 || min(count_matrix[, j]) > 0) {
        status <- "no_zero"
      } else {
        fit_ij <- fit_direction(count_matrix[, i], count_matrix[, j])
        fit_ji <- fit_direction(count_matrix[, j], count_matrix[, i])
        if (fit_ij$ok && fit_ji$ok) {
          p_ij <- fit_ij$p
          p_ji <- fit_ji$p
          status <- "tested"
        } else {
          status <- "fit_failed"
        }
      }
    }
    records[[k]] <- data.frame(
      pair_index = pair_index,
      gene_i = genes[[i]],
      gene_j = genes[[j]],
      nonzero_i = nonzero_prop[[i]],
      nonzero_j = nonzero_prop[[j]],
      eligible = eligible,
      status = status,
      p_count_ij = p_ij[["count"]],
      p_count_ji = p_ji[["count"]],
      p_detection_ij = p_ij[["detection"]],
      p_detection_ji = p_ji[["detection"]],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, records)
}

if (.Platform$OS.type == "windows" && workers > 1L) {
  stop("Parallel workers require a Unix-like platform; set --workers 1")
}
parts <- if (workers == 1L) {
  lapply(worker_indices, run_worker)
} else {
  parallel::mclapply(worker_indices, run_worker, mc.cores = workers, mc.preschedule = TRUE)
}
pair_results <- do.call(rbind, parts)
pair_results <- pair_results[order(pair_results$pair_index), , drop = FALSE]
if (nrow(pair_results) != nrow(pair_matrix) || anyDuplicated(pair_results$pair_index)) {
  stop("Incomplete or duplicated pair results")
}

eligible <- pair_results[pair_results$eligible, , drop = FALSE]
directional_count <- rbind(
  data.frame(
    Gene1 = eligible$gene_i, Gene2 = eligible$gene_j,
    Pvalue = eligible$p_count_ij, Direction = "i_to_j",
    stringsAsFactors = FALSE
  ),
  data.frame(
    Gene1 = eligible$gene_j, Gene2 = eligible$gene_i,
    Pvalue = eligible$p_count_ji, Direction = "j_to_i",
    stringsAsFactors = FALSE
  )
)
directional_detection <- rbind(
  data.frame(
    Gene1 = eligible$gene_i, Gene2 = eligible$gene_j,
    Pvalue = eligible$p_detection_ij, Direction = "i_to_j",
    stringsAsFactors = FALSE
  ),
  data.frame(
    Gene1 = eligible$gene_j, Gene2 = eligible$gene_i,
    Pvalue = eligible$p_detection_ji, Direction = "j_to_i",
    stringsAsFactors = FALSE
  )
)
directional_count <- directional_count[order(directional_count$Gene1, directional_count$Gene2), ]
directional_detection <- directional_detection[
  order(directional_detection$Gene1, directional_detection$Gene2),
]
keys <- c("Gene1", "Gene2", "Direction")
if (!identical(directional_count[, keys], directional_detection[, keys])) {
  stop("Count and detection directional records are misaligned")
}

thresholds <- as.numeric(strsplit(opt$p_thresholds, ",", fixed = TRUE)[[1L]])
hurdle_summary <- do.call(rbind, lapply(thresholds, function(threshold) {
  count_flag <- directional_count$Pvalue < threshold
  detection_flag <- directional_detection$Pvalue < threshold
  data.frame(
    p_threshold = threshold,
    count_fraction = safe_fraction(count_flag),
    detection_fraction = safe_fraction(detection_flag),
    component_union_fraction = safe_fraction(count_flag | detection_flag),
    eligible_pairs = nrow(eligible),
    fit_failed_pairs = sum(eligible$status != "tested"),
    stringsAsFactors = FALSE
  )
}))

donor_levels <- unique(donors)
donor_factor <- factor(donors, levels = donor_levels)
donor_sums <- rowsum(count_matrix, donor_factor, reorder = FALSE)
donor_n <- as.numeric(table(donor_factor))
donor_means <- sweep(donor_sums, 1L, donor_n, "/")
rho <- suppressWarnings(cor(donor_means, method = "spearman", use = "pairwise.complete.obs"))
t_stat <- rho * sqrt((length(donor_levels) - 2) / (1 - rho^2 + 1e-10))
p_matrix <- 2 * pt(-abs(t_stat), df = length(donor_levels) - 2)
diag(p_matrix) <- NA_real_
eligible_pair_matrix <- pair_matrix[
  apply(pair_matrix, 1L, function(pair) min(nonzero_prop[pair]) >= opt$nonzero_cutoff),
  , drop = FALSE
]
pb_p <- p_matrix[eligible_pair_matrix]
pb_directional <- rbind(
  data.frame(
    Gene1 = genes[eligible_pair_matrix[, 1L]],
    Gene2 = genes[eligible_pair_matrix[, 2L]],
    Pvalue = pb_p,
    Direction = "i_to_j",
    stringsAsFactors = FALSE
  ),
  data.frame(
    Gene1 = genes[eligible_pair_matrix[, 2L]],
    Gene2 = genes[eligible_pair_matrix[, 1L]],
    Pvalue = pb_p,
    Direction = "j_to_i",
    stringsAsFactors = FALSE
  )
)
pb_directional <- pb_directional[order(pb_directional$Gene1, pb_directional$Gene2), ]
pb_summary <- do.call(rbind, lapply(thresholds, function(threshold) {
  data.frame(
    p_threshold = threshold,
    n_pairs_tested = length(pb_p),
    positive_fraction = safe_fraction(pb_p < threshold),
    stringsAsFactors = FALSE
  )
}))

sc_dir <- file.path(opt$output_dir, "sc_hurdle")
pb_dir <- file.path(opt$output_dir, "pb_spearman")
atomic_write_table(pair_results, file.path(sc_dir, "pair_results.tsv.gz"))
atomic_write_table(directional_count, file.path(sc_dir, "sc_hurdle_count_pvalues.tsv"))
atomic_write_table(directional_detection, file.path(sc_dir, "sc_hurdle_zero_pvalues.tsv"))
atomic_write_table(hurdle_summary, file.path(sc_dir, "fraction_summary.tsv"))
atomic_write_table(pb_directional, file.path(pb_dir, "pb_spearman_pvalues.tsv"))
atomic_write_table(pb_summary, file.path(pb_dir, "fraction_summary.tsv"))
atomic_write_table(data.frame(
  package = c("R", "fasthurdle", "fastglm"),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    loaded_fasthurdle,
    loaded_fastglm
  )
), file.path(opt$output_dir, "analysis_versions.tsv"))

message(sprintf(
  "Completed %d pairs (%d eligible) with %d workers",
  nrow(pair_matrix), nrow(eligible), workers
))
