#!/usr/bin/env Rscript

.local_lib <- Sys.getenv("SC_PCQTL_R_LIBS", unset = "")
if (nzchar(.local_lib)) .libPaths(c(.local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(fasthurdle)
})

if (utils::packageVersion("fasthurdle") != package_version("1.2.0")) {
  stop("This sensitivity analysis requires fasthurdle 1.2.0")
}

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args_all[grep("--file=", args_all)[1]]))
module_dir <- dirname(script_path)
source(file.path(module_dir, "load_config.R"))
source(file.path(module_dir, "R", "joint_score_utils.R"))
cfg <- load_joint_score_config(module_dir)
setDTthreads(1L)

args <- commandArgs(trailingOnly = TRUE)
task_id_value <- if (length(args)) as.integer(args[1]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", ""))
if (!is.finite(task_id_value)) stop("Provide task_id or SLURM_ARRAY_TASK_ID")

manifest_file <- file.path(cfg$stage_dir, "joint_score_task_manifest.tsv")
if (!file.exists(manifest_file)) stop("Run 02_plan_chunks.R first")
manifest <- fread(manifest_file)
task <- manifest[task_id == task_id_value]
if (nrow(task) != 1L) stop("Task id absent or duplicated: ", task_id_value)

out_dir <- file.path(cfg$chunk_dir, task$chunk_id)
if (file.exists(file.path(out_dir, "COMPLETE")) && file.exists(file.path(out_dir, "chunk_summary.tsv"))) {
  message(sprintf("[%s] %s already complete; skipping", cfg$celltype, task$chunk_id))
  quit(status = 0L)
}

chr <- task$chromosome
genes_all <- fread(file.path(cfg$stage_dir, "filtered_autosomal_genes.tsv"))
chr_info <- genes_all[chr_numeric == chr][order(chr_rank)]
gene_names <- chr_info$gene_name
counts <- readRDS(file.path(cfg$stage_dir, sprintf("chr%d_counts.rds", chr)))
covars <- readRDS(file.path(cfg$stage_dir, "covariates.rds"))
if (!identical(colnames(counts), gene_names)) stop("Staged count columns do not match filtered-gene order")
if (nrow(counts) != nrow(covars)) stop("Count/covariate row mismatch")

count_formula <- ~ age + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pf1 + pf2
zero_formula <- ~ log_total_read_counts + age + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pf1 + pf2
X_null <- model.matrix(count_formula, data = covars)
Z_null <- model.matrix(zero_formula, data = covars)
offsetx <- covars$log_total_read_counts
weights <- rep.int(1, nrow(covars))

stage_summary <- fread(file.path(cfg$stage_dir, "stage_summary.tsv"))
threshold <- stage_summary$bonferroni_joint_threshold[[1]]
log_threshold <- log(threshold)
batch_count_cpp <- getFromNamespace("score_test_count_batch_cpp", "fasthurdle")

run_response <- function(response_idx) {
  predictor_idx <- setdiff(
    seq.int(max(1L, response_idx - 49L), min(length(gene_names), response_idx + 49L)),
    response_idx
  )
  y <- as.numeric(counts[, response_idx])
  predictor_mat <- counts[, predictor_idx, drop = FALSE]
  base <- data.table(
    celltype = cfg$celltype,
    chromosome = chr,
    response = gene_names[response_idx],
    predictor = gene_names[predictor_idx],
    response_rank = response_idx,
    predictor_rank = predictor_idx
  )

  if (!any(y > 0) || length(unique(y)) < 2L) {
    base[, `:=`(
      beta_count = NA_real_, se_count = NA_real_, stat_count = NA_real_,
      beta_detection = NA_real_, se_detection = NA_real_, stat_detection = NA_real_,
      stat_joint = NA_real_, log_p_joint = NA_real_, neglog10_p_joint = NA_real_,
      p_joint = NA_real_, significant_direction = FALSE,
      count_null_converged = FALSE, detection_null_converged = FALSE,
      count_score_cache_valid = FALSE, detection_score_cache_valid = FALSE
    )]
    return(base)
  }

  count_null <- tryCatch(
    fit_null_count(X_null, y, offsetx = offsetx, weights = weights, dist = "poisson"),
    error = function(e) NULL
  )
  detection_null <- tryCatch(
    fit_null_zero(Z_null, y, weights = weights),
    error = function(e) NULL
  )
  count_ok <- !is.null(count_null) && isTRUE(count_null$convergence == 0)
  detection_ok <- !is.null(detection_null) && isTRUE(detection_null$convergence == 0)

  count_res <- NULL
  count_cache_ok <- FALSE
  if (count_ok) {
    count_null <- tryCatch(
      prepare_score_cache_count(count_null, y, X_null, offsetx = offsetx, weights = weights),
      error = function(e) NULL
    )
    sc <- if (!is.null(count_null)) count_null$score_cache else NULL
    if (!is.null(sc) && isTRUE(sc$valid)) {
      count_cache_ok <- TRUE
      count_res <- tryCatch(
        batch_count_cpp(
          predictor_mat[sc$Y1 + 1L, , drop = FALSE],
          sc$Y1, sc$grad_weights, sc$v_ee, sc$Y_pos,
          sc$I_nn_inv, sc$I_nn_beta_inv, sc$beta_inv_ok,
          sc$Xnull_vee_t, sc$X_null_pos, sc$w_pos, sc$theta,
          sc$beta_null, sc$eta_null_pos, sc$mu_pos, sc$p0_pos,
          sc$log_p1_pos, sc$kx_null, FALSE, FALSE, 1e30, NULL
        ),
        error = function(e) NULL
      )
    }
  }

  detection_beta <- detection_se <- detection_stat <- rep(NA_real_, length(predictor_idx))
  detection_cache_ok <- FALSE
  if (detection_ok) {
    detection_null <- tryCatch(
      prepare_score_cache_zero(detection_null, y, Z_null, weights = weights),
      error = function(e) NULL
    )
    cache <- if (!is.null(detection_null)) detection_null$score_cache else NULL
    if (!is.null(cache) && isTRUE(cache$valid)) {
      detection_cache_ok <- TRUE
      for (j in seq_along(predictor_idx)) {
        result <- tryCatch(
          score_test_zero(
            Z_null, predictor_mat[, j], y,
            weights = weights, null_fit_zero = detection_null, spa_cutoff = NULL
          ),
          error = function(e) NULL
        )
        if (!is.null(result)) {
          detection_beta[j] <- result$beta[[1]]
          detection_se[j] <- result$se[[1]]
          detection_stat[j] <- result$statistic[[1]]
        }
      }
    }
  }

  if (is.null(count_res)) {
    count_beta <- count_se <- count_stat <- rep(NA_real_, length(predictor_idx))
  } else {
    count_beta <- as.numeric(count_res$beta)
    count_se <- as.numeric(count_res$se)
    count_stat <- as.numeric(count_res$statistic)
  }

  stat_joint <- count_stat + detection_stat
  log_p_joint <- joint_score_p(count_stat, detection_stat, log.p = TRUE)
  base[, `:=`(
    beta_count = count_beta,
    se_count = count_se,
    stat_count = count_stat,
    beta_detection = detection_beta,
    se_detection = detection_se,
    stat_detection = detection_stat,
    stat_joint = stat_joint,
    log_p_joint = log_p_joint,
    neglog10_p_joint = -log_p_joint / log(10),
    p_joint = exp(log_p_joint),
    significant_direction = is.finite(log_p_joint) & log_p_joint < log_threshold,
    count_null_converged = count_ok,
    detection_null_converged = detection_ok,
    count_score_cache_valid = count_cache_ok,
    detection_score_cache_valid = detection_cache_ok
  )]
  base
}

start_time <- Sys.time()
result <- rbindlist(lapply(seq.int(task$response_start, task$response_end), run_response), fill = TRUE)
result[, rank_low := pmin(response_rank, predictor_rank)]
result[, rank_high := pmax(response_rank, predictor_rank)]
result[, Gene1 := gene_names[rank_low]]
result[, Gene2 := gene_names[rank_high]]
result[, direction := paste0(response, "<-", predictor)]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(result, file.path(out_dir, "directional_joint_tests.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(result[significant_direction == TRUE], file.path(out_dir, "significant_directions.tsv"), sep = "\t")
summary <- data.table(
  task_id = task_id_value,
  chunk_id = task$chunk_id,
  chromosome = chr,
  response_start = task$response_start,
  response_end = task$response_end,
  n_directional_tests = nrow(result),
  n_complete_joint_tests = sum(is.finite(result$stat_joint)),
  n_significant_directions = sum(result$significant_direction),
  count_null_failures = uniqueN(result[count_null_converged == FALSE, response]),
  detection_null_failures = uniqueN(result[detection_null_converged == FALSE, response]),
  count_score_cache_failures = uniqueN(result[count_score_cache_valid == FALSE, response]),
  detection_score_cache_failures = uniqueN(result[detection_score_cache_valid == FALSE, response]),
  elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
  threshold = threshold
)
fwrite(summary, file.path(out_dir, "chunk_summary.tsv"), sep = "\t")
writeLines("OK", file.path(out_dir, "COMPLETE"))
message(sprintf(
  "[%s] %s complete: %d directional tests, %d significant",
  cfg$celltype, task$chunk_id, nrow(result), sum(result$significant_direction)
))
