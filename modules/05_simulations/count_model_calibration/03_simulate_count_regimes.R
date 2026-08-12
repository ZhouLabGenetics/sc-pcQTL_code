#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fasthurdle)
  library(parallel)
})

required_fasthurdle <- base::package_version("1.2.0")
if (utils::packageVersion("fasthurdle") != required_fasthurdle) {
  stop(
    "This simulation requires fasthurdle 1.2.0; found ",
    as.character(utils::packageVersion("fasthurdle")),
    call. = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 5L) {
  stop(
    paste(
      "Usage: Rscript 03_simulate_count_regimes.R OUTPUT_DIR",
      "[N_REPLICATES] [N_CELLS] [N_CORES] [SEED]"
    ),
    call. = FALSE
  )
}

output_dir <- args[[1L]]
n_replicates <- if (length(args) >= 2L) as.integer(args[[2L]]) else 2000L
n_cells <- if (length(args) >= 3L) as.integer(args[[3L]]) else 20000L
n_cores <- if (length(args) >= 4L) as.integer(args[[4L]]) else 1L
base_seed <- if (length(args) >= 5L) as.integer(args[[5L]]) else 20260812L

if (any(!is.finite(c(n_replicates, n_cells, n_cores, base_seed)))) {
  stop("Simulation arguments must be finite integers", call. = FALSE)
}
if (n_replicates < 1L || n_cells < 100L || n_cores < 1L) {
  stop("N_REPLICATES and N_CORES must be positive; N_CELLS must be at least 100")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

detection_probability <- 0.12
positive_mean <- 2
nb_size <- 2
shifted_poisson_rate_x <- 1.0
shifted_poisson_rate_y <- 0.6
p_thresholds <- c(0.05, 0.01, 0.005, 0.001)

solve_ztpois_rate <- function(target_mean) {
  stats::uniroot(
    function(rate) rate / (1 - exp(-rate)) - target_mean,
    interval = c(1e-8, target_mean)
  )$root
}

solve_ztnb_mean <- function(target_mean, size) {
  stats::uniroot(
    function(mu) {
      p_zero <- stats::dnbinom(0, size = size, mu = mu)
      mu / (1 - p_zero) - target_mean
    },
    interval = c(1e-8, target_mean)
  )$root
}

ztpois_rate <- solve_ztpois_rate(positive_mean)
ztnb_mean <- solve_ztnb_mean(positive_mean, nb_size)

sample_zero_truncated <- function(n, generator) {
  values <- generator(n)
  zero <- which(values == 0L)
  while (length(zero) > 0L) {
    values[zero] <- generator(length(zero))
    zero <- which(values == 0L)
  }
  values
}

generate_gene <- function(regime, gene) {
  detected <- stats::rbinom(n_cells, 1L, detection_probability)
  values <- integer(n_cells)
  n_positive <- sum(detected)
  if (n_positive == 0L) return(values)

  if (regime == "poisson_generated") {
    values[detected == 1L] <- sample_zero_truncated(
      n_positive,
      function(n) stats::rpois(n, lambda = ztpois_rate)
    )
  } else if (regime == "shifted_poisson_generated") {
    rate <- if (gene == "x") shifted_poisson_rate_x else shifted_poisson_rate_y
    values[detected == 1L] <- 1L + stats::rpois(n_positive, lambda = rate)
  } else if (regime == "negative_binomial_generated") {
    values[detected == 1L] <- sample_zero_truncated(
      n_positive,
      function(n) stats::rnbinom(n, size = nb_size, mu = ztnb_mean)
    )
  } else {
    stop("Unknown regime: ", regime)
  }
  values
}

extract_count_p <- function(fit, predictor = "x") {
  coefficients <- summary(fit)$coefficients$count
  if (is.null(coefficients) || !predictor %in% rownames(coefficients)) return(NA_real_)
  as.numeric(coefficients[predictor, "Pr(>|z|)"])
}

fit_model <- function(data, count_dist) {
  warning_messages <- character()
  fit <- tryCatch(
    withCallingHandlers(
      fasthurdle(
        y ~ x,
        data = data,
        dist = count_dist,
        zero.dist = "binomial"
      ),
      warning = function(warning) {
        warning_messages <<- c(warning_messages, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) NULL
  )

  if (is.null(fit)) {
    return(list(p = NA_real_, theta = NA_real_, converged = FALSE, warned = length(warning_messages) > 0L))
  }

  list(
    p = extract_count_p(fit),
    theta = if (count_dist == "negbin") unname(fit$theta[["count"]]) else NA_real_,
    converged = isTRUE(fit$converged),
    warned = length(warning_messages) > 0L
  )
}

run_replicate <- function(index, regime) {
  seed_offset <- c(
    poisson_generated = 0L,
    shifted_poisson_generated = 10000000L,
    negative_binomial_generated = 20000000L
  )[[regime]]
  set.seed(base_seed + index + seed_offset)
  data <- data.frame(
    x = generate_gene(regime, "x"),
    y = generate_gene(regime, "y")
  )
  poisson <- fit_model(data, "poisson")
  negbin <- fit_model(data, "negbin")

  data.table(
    regime = regime,
    replicate = index,
    model = c("Poisson", "Negative-Binomial"),
    p_value = c(poisson$p, negbin$p),
    theta = c(poisson$theta, negbin$theta),
    converged = c(poisson$converged, negbin$converged),
    warned = c(poisson$warned, negbin$warned),
    detection_x = mean(data$x > 0),
    detection_y = mean(data$y > 0),
    positive_mean_x = mean(data$x[data$x > 0]),
    positive_mean_y = mean(data$y[data$y > 0])
  )
}

run_regime <- function(regime) {
  message("Running ", regime, " regime")
  results <- mclapply(
    seq_len(n_replicates),
    run_replicate,
    regime = regime,
    mc.cores = n_cores,
    mc.preschedule = TRUE
  )
  rbindlist(results, use.names = TRUE)
}

results <- rbindlist(lapply(
  c(
    "poisson_generated",
    "shifted_poisson_generated",
    "negative_binomial_generated"
  ),
  run_regime
))
fwrite(results, file.path(output_dir, "controlled_count_regime_replicates.tsv"), sep = "\t")

calibration <- results[
  is.finite(p_value),
  .(
    empirical_rejection_rate = vapply(
      p_thresholds,
      function(threshold) mean(p_value < threshold),
      numeric(1)
    ),
    nominal_alpha = p_thresholds
  ),
  by = .(regime, model)
]
calibration[, rejection_ratio := empirical_rejection_rate / nominal_alpha]
fwrite(calibration, file.path(output_dir, "controlled_count_regime_calibration.tsv"), sep = "\t")

fit_summary <- results[, .(
  n_replicates = .N,
  n_finite_p = sum(is.finite(p_value)),
  convergence_fraction = mean(converged),
  warning_fraction = mean(warned),
  median_theta = if (all(is.na(theta))) NA_real_ else median(theta, na.rm = TRUE),
  theta_q95 = if (all(is.na(theta))) NA_real_ else unname(stats::quantile(theta, 0.95, na.rm = TRUE)),
  theta_ge_1000_fraction = if (all(is.na(theta))) NA_real_ else mean(theta >= 1000, na.rm = TRUE)
), by = .(regime, model)]
fwrite(fit_summary, file.path(output_dir, "controlled_count_regime_fit_summary.tsv"), sep = "\t")

fwrite(
  data.table(
    parameter = c(
      "fasthurdle_version", "n_replicates", "n_cells", "n_cores", "base_seed",
      "detection_probability", "positive_mean", "nb_size", "ztpois_rate", "ztnb_mean",
      "shifted_poisson_rate_x", "shifted_poisson_rate_y"
    ),
    value = c(
      as.character(utils::packageVersion("fasthurdle")), n_replicates, n_cells,
      n_cores, base_seed, detection_probability, positive_mean, nb_size,
      ztpois_rate, ztnb_mean, shifted_poisson_rate_x, shifted_poisson_rate_y
    )
  ),
  file.path(output_dir, "controlled_count_regime_parameters.tsv"),
  sep = "\t"
)

message("Controlled count-regime results written to ", normalizePath(output_dir))
