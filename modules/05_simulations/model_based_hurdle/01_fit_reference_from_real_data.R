#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--count_file", type = "character",
              default = Sys.getenv("COQTL_RAW_COUNTS_FILE", unset = "")),
  make_option("--out_dir", type = "character", default = "reference"),
  make_option("--n_reference_cells", type = "integer", default = 20000L),
  make_option("--n_reference_genes", type = "integer", default = 3000L),
  make_option("--seed", type = "integer", default = 1L)
)

opt <- parse_args(OptionParser(option_list = option_list))
if (!nzchar(opt$count_file)) {
  stop("--count_file is required, or set COQTL_RAW_COUNTS_FILE.")
}
set.seed(opt$seed)

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading header")
header <- names(fread(opt$count_file, nrows = 0))

meta_cols <- c(
  "CellID", "barcode", "IndividualID", "individual", "CellType",
  "sex", paste0("pc", 1:6), "age", "pf1", "pf2",
  "total_read_counts", "log_total_read_counts"
)
gene_cols <- setdiff(header, meta_cols)
if (!"individual" %in% header) {
  stop("Expected column 'individual' in count_file.")
}

message("Reading individual and library-size columns")
select_cols <- intersect(c("individual", "total_read_counts", "log_total_read_counts"), header)
meta_dt <- fread(opt$count_file, select = select_cols)

if (!"total_read_counts" %in% names(meta_dt)) {
  stop("Expected column 'total_read_counts' in count_file.")
}
if (!"log_total_read_counts" %in% names(meta_dt)) {
  meta_dt[, log_total_read_counts := log(pmax(total_read_counts, 1))]
}

message("Summarizing donor cell counts")
donor_counts <- meta_dt[, .N, by = individual][order(individual)]$N

n_cells_use <- min(opt$n_reference_cells, nrow(meta_dt))
cell_idx <- sort(sample.int(nrow(meta_dt), n_cells_use))
lib_sizes <- meta_dt$total_read_counts[cell_idx]
log_lib_sizes <- meta_dt$log_total_read_counts[cell_idx]

n_genes_use <- min(opt$n_reference_genes, length(gene_cols))
gene_sample <- sort(sample(gene_cols, n_genes_use))

message("Reading sampled gene matrix")
ref_dt <- fread(
  opt$count_file,
  select = c("individual", "total_read_counts", "log_total_read_counts", gene_sample)
)
ref_dt <- ref_dt[cell_idx]

count_mat <- as.matrix(ref_dt[, ..gene_sample])
nonzero_rate <- colMeans(count_mat > 0)
positive_mean <- vapply(seq_along(gene_sample), function(i) {
  vals <- count_mat[, i]
  nz <- vals[vals > 0]
  if (!length(nz)) {
    return(NA_real_)
  }
  mean(nz)
}, numeric(1))
overall_mean <- colMeans(count_mat)
lib_center <- mean(ref_dt$log_total_read_counts, na.rm = TRUE)
log_lib_centered <- ref_dt$log_total_read_counts - lib_center

estimate_ztnb_params <- function(y_pos, exposure) {
  keep <- is.finite(y_pos) & is.finite(exposure) & y_pos > 0 & exposure > 0
  y_pos <- y_pos[keep]
  exposure <- exposure[keep]
  if (!length(y_pos)) {
    return(list(log_rate = NA_real_, theta = NA_real_))
  }

  rate_init <- max(sum(y_pos) / sum(exposure), 1e-8)
  var_y <- stats::var(y_pos)
  mean_y <- mean(y_pos)
  theta_init <- if (is.finite(var_y) && var_y > mean_y) {
    max(mean_y^2 / (var_y - mean_y), 0.5)
  } else {
    10
  }

  nll <- function(par) {
    log_rate <- par[1]
    log_theta <- par[2]
    theta <- exp(log_theta)
    mu <- pmax(exp(log_rate) * exposure, 1e-8)
    p0 <- (theta / (theta + mu))^theta
    loglik <- stats::dnbinom(y_pos, size = theta, mu = mu, log = TRUE) -
      log1p(-p0)
    if (any(!is.finite(loglik))) {
      return(Inf)
    }
    -sum(loglik)
  }

  fit <- tryCatch(
    stats::optim(
      par = c(log(rate_init), log(theta_init)),
      fn = nll,
      method = "BFGS",
      control = list(maxit = 200)
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || !is.finite(fit$value)) {
    return(list(log_rate = log(rate_init), theta = theta_init))
  }

  list(
    log_rate = fit$par[1],
    theta = exp(fit$par[2])
  )
}

message("Fitting gene-wise hurdle reference parameters")
alpha_intercept <- numeric(length(gene_sample))
alpha_lib <- numeric(length(gene_sample))
log_rate <- numeric(length(gene_sample))
theta <- numeric(length(gene_sample))

for (i in seq_along(gene_sample)) {
  y <- count_mat[, i]
  detected <- as.integer(y > 0)

  if (length(unique(detected)) < 2) {
    p_hat <- mean(detected)
    alpha_intercept[i] <- qlogis(pmin(pmax(p_hat, 1e-4), 1 - 1e-4))
    alpha_lib[i] <- 0
  } else {
    fit_zero <- suppressWarnings(
      glm.fit(
        x = cbind("(Intercept)" = 1, log_lib_centered = log_lib_centered),
        y = detected,
        family = stats::binomial()
      )
    )
    alpha_intercept[i] <- fit_zero$coefficients[1]
    alpha_lib[i] <- fit_zero$coefficients[2]
  }

  y_pos <- y[y > 0]
  exp_pos <- ref_dt$total_read_counts[y > 0]
  fit_count <- estimate_ztnb_params(y_pos, exp_pos)
  log_rate[i] <- fit_count$log_rate
  theta[i] <- fit_count$theta
}

alpha_lib[!is.finite(alpha_lib)] <- 0
alpha_intercept[!is.finite(alpha_intercept)] <- qlogis(pmin(pmax(nonzero_rate[!is.finite(alpha_intercept)], 1e-4), 1 - 1e-4))
log_rate[!is.finite(log_rate)] <- log(pmax(overall_mean[!is.finite(log_rate)] / pmax(mean(ref_dt$total_read_counts), 1), 1e-10))
theta[!is.finite(theta) | theta <= 0] <- 10

gene_ref <- data.table(
  gene_name = gene_sample,
  nonzero_rate = nonzero_rate,
  positive_mean = positive_mean,
  overall_mean = overall_mean,
  alpha_intercept = alpha_intercept,
  alpha_lib = alpha_lib,
  log_rate = log_rate,
  theta = theta
)
gene_ref <- gene_ref[is.finite(alpha_intercept) & is.finite(alpha_lib) & is.finite(log_rate) & is.finite(theta)]

reference <- list(
  source_count_file = opt$count_file,
  n_total_cells = nrow(meta_dt),
  n_total_genes = length(gene_cols),
  donor_cell_counts = donor_counts,
  library_sizes = lib_sizes,
  log_library_sizes = log_lib_sizes,
  lib_center = lib_center,
  gene_reference = gene_ref,
  defaults = list(
    sigma_zero = 0.00,
    sigma_count = 0.00,
    nonzero_cutoff = 0.01,
    stable_nonzero_min = 0.02,
    stable_nonzero_max = 0.80,
    stable_log_rate_min = -12,
    stable_log_rate_max = -8,
    stable_theta_min = 1,
    stable_theta_max = 100,
    stable_alpha_lib_min = -2,
    stable_alpha_lib_max = 8
  )
)

summary_dt <- data.table(
  metric = c(
    "n_total_cells", "n_total_genes", "n_donors",
    "donor_cells_min", "donor_cells_median", "donor_cells_mean", "donor_cells_max",
    "lib_size_min", "lib_size_median", "lib_size_mean", "lib_size_max",
    "gene_nonzero_min", "gene_nonzero_median", "gene_nonzero_mean", "gene_nonzero_max"
  ),
  value = c(
    reference$n_total_cells,
    reference$n_total_genes,
    length(reference$donor_cell_counts),
    min(reference$donor_cell_counts),
    median(reference$donor_cell_counts),
    mean(reference$donor_cell_counts),
    max(reference$donor_cell_counts),
    min(reference$library_sizes),
    median(reference$library_sizes),
    mean(reference$library_sizes),
    max(reference$library_sizes),
    min(reference$gene_reference$nonzero_rate),
    median(reference$gene_reference$nonzero_rate),
    mean(reference$gene_reference$nonzero_rate),
    max(reference$gene_reference$nonzero_rate)
  )
)

saveRDS(reference, file.path(opt$out_dir, "reference_stats.rds"))
fwrite(summary_dt, file.path(opt$out_dir, "reference_summary.tsv"), sep = "\t")
fwrite(reference$gene_reference, file.path(opt$out_dir, "reference_gene_parameters.tsv"), sep = "\t")

message("Saved reference to: ", file.path(opt$out_dir, "reference_stats.rds"))
