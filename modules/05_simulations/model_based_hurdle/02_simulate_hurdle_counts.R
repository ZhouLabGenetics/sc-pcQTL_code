#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--reference_rds", type = "character", default = "reference/reference_stats.rds"),
  make_option("--scenario", type = "character", default = "null"),
  make_option("--replicate", type = "integer", default = 1L),
  make_option("--out_dir", type = "character", default = "simulations"),
  make_option("--n_donors", type = "integer", default = 500L),
  make_option("--n_genes", type = "integer", default = 200L),
  make_option("--signal_pair_fraction", type = "double", default = 0.05),
  make_option("--sigma_zero", type = "double", default = NA_real_),
  make_option("--sigma_count", type = "double", default = NA_real_),
  make_option("--shared_zero_lib_coef", type = "double", default = 0.15),
  make_option("--seed", type = "integer", default = 1L)
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed + opt$replicate)

reference <- readRDS(opt$reference_rds)

if (is.na(opt$sigma_zero)) opt$sigma_zero <- reference$defaults$sigma_zero
if (is.na(opt$sigma_count)) opt$sigma_count <- reference$defaults$sigma_count

scenario_map <- list(
  null = list(gamma_zero = 0, gamma_count = 0),
  zero_only_low = list(gamma_zero = 0.35, gamma_count = 0),
  zero_only_high = list(gamma_zero = 0.60, gamma_count = 0),
  count_only_low = list(gamma_zero = 0, gamma_count = 0.20),
  count_only_high = list(gamma_zero = 0, gamma_count = 0.35),
  both_low = list(gamma_zero = 0.35, gamma_count = 0.20),
  both_high = list(gamma_zero = 0.60, gamma_count = 0.35)
)
if (!opt$scenario %in% names(scenario_map)) {
  stop("Unknown scenario: ", opt$scenario)
}
gammas <- scenario_map[[opt$scenario]]

sample_ztpois <- function(lambda) {
  lambda <- pmax(lambda, 1e-8)
  p0 <- exp(-lambda)
  u <- stats::runif(length(lambda), min = p0, max = 1)
  stats::qpois(u, lambda = lambda)
}

rep_dir <- file.path(opt$out_dir, opt$scenario, sprintf("replicate_%03d", opt$replicate))
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

donor_ids <- sprintf("donor_%03d", seq_len(opt$n_donors))
donor_cell_counts <- sample(reference$donor_cell_counts, size = opt$n_donors, replace = TRUE)
cell_donor <- rep(donor_ids, donor_cell_counts)
n_cells <- length(cell_donor)

library_pool <- reference$library_sizes[is.finite(reference$library_sizes) & reference$library_sizes > 0]
lib_q <- stats::quantile(library_pool, probs = c(0.10, 0.90), na.rm = TRUE)
library_pool <- library_pool[library_pool >= lib_q[1] & library_pool <= lib_q[2]]
library_size <- sample(library_pool, size = n_cells, replace = TRUE)
log_library_size <- log(pmax(library_size, 1))
centered_log_library_size <- log_library_size - mean(log(library_pool))

gene_ref <- reference$gene_reference
gene_ref <- gene_ref[
  nonzero_rate >= max(reference$defaults$nonzero_cutoff, reference$defaults$stable_nonzero_min) &
  nonzero_rate <= reference$defaults$stable_nonzero_max &
  log_rate >= reference$defaults$stable_log_rate_min &
  log_rate <= reference$defaults$stable_log_rate_max &
  theta >= reference$defaults$stable_theta_min &
  theta <= reference$defaults$stable_theta_max &
  alpha_lib >= reference$defaults$stable_alpha_lib_min &
  alpha_lib <= reference$defaults$stable_alpha_lib_max
]

central_keep <- rep(TRUE, nrow(gene_ref))
for (nm in c("alpha_intercept", "log_rate", "theta", "nonzero_rate")) {
  qs <- stats::quantile(gene_ref[[nm]], probs = c(0.10, 0.90), na.rm = TRUE)
  central_keep <- central_keep & gene_ref[[nm]] >= qs[1] & gene_ref[[nm]] <= qs[2]
}
gene_ref <- gene_ref[central_keep]
if (nrow(gene_ref) < opt$n_genes) {
  stop("Not enough genes in reference to sample ", opt$n_genes, " genes.")
}
gene_idx <- sample(seq_len(nrow(gene_ref)), size = opt$n_genes, replace = FALSE)
gene_par <- copy(gene_ref[gene_idx])
gene_par[, gene_name := sprintf("gene_%03d", seq_len(.N))]
gene_par[, alpha_intercept := pmin(pmax(alpha_intercept, -3), -1)]
gene_par[, log_rate := pmin(pmax(log_rate, -10.2), -9.0)]

count_mat <- matrix(0L, nrow = n_cells, ncol = opt$n_genes)
colnames(count_mat) <- gene_par$gene_name

u_zero <- matrix(0, nrow = opt$n_donors, ncol = opt$n_genes)
u_count <- matrix(0, nrow = opt$n_donors, ncol = opt$n_genes)
donor_index <- match(cell_donor, donor_ids)

all_pairs <- as.data.table(t(combn(gene_par$gene_name, 2)))
setnames(all_pairs, c("gene1", "gene2"))
all_pairs[, pair_id := .I]

signal_gene_set <- gene_par[
  nonzero_rate >= stats::median(nonzero_rate) &
  log_rate >= stats::median(log_rate)
]$gene_name
signal_pairs_universe <- all_pairs[gene1 %in% signal_gene_set & gene2 %in% signal_gene_set]
if (!nrow(signal_pairs_universe)) {
  signal_pairs_universe <- all_pairs
}

n_signal_pairs <- if (opt$scenario == "null") 0L else max(1L, floor(nrow(all_pairs) * opt$signal_pair_fraction))
signal_pairs <- if (n_signal_pairs > 0) {
  copy(signal_pairs_universe[sample(.N, min(.N, n_signal_pairs))])
} else {
  copy(all_pairs[0])
}
signal_pairs[, `:=`(
  scenario = opt$scenario,
  gamma_zero = gammas$gamma_zero,
  gamma_count = gammas$gamma_count
)]

signal_zero <- matrix(0, nrow = n_cells, ncol = opt$n_genes)
signal_count <- matrix(0, nrow = n_cells, ncol = opt$n_genes)

if (nrow(signal_pairs) > 0) {
  gene_to_idx <- setNames(seq_len(opt$n_genes), gene_par$gene_name)
  for (k in seq_len(nrow(signal_pairs))) {
    h <- rnorm(n_cells)
    i <- gene_to_idx[[signal_pairs$gene1[k]]]
    j <- gene_to_idx[[signal_pairs$gene2[k]]]
    signal_zero[, i] <- signal_zero[, i] + gammas$gamma_zero * h
    signal_zero[, j] <- signal_zero[, j] + gammas$gamma_zero * h
    signal_count[, i] <- signal_count[, i] + gammas$gamma_count * h
    signal_count[, j] <- signal_count[, j] + gammas$gamma_count * h
  }
}

for (g in seq_len(opt$n_genes)) {
  eta_zero <- gene_par$alpha_intercept[g] +
    opt$shared_zero_lib_coef * centered_log_library_size +
    u_zero[donor_index, g] +
    signal_zero[, g]
  pi_g <- plogis(eta_zero)
  detected <- rbinom(n_cells, size = 1, prob = pi_g)

  eta_count <- gene_par$log_rate[g] +
    log_library_size +
    u_count[donor_index, g] +
    signal_count[, g]
  lambda_g <- exp(eta_count)

  y_pos <- sample_ztpois(lambda_g)
  count_mat[, g] <- ifelse(detected == 1L, y_pos, 0L)
}

sex_by_donor <- sample(c("M", "F"), size = opt$n_donors, replace = TRUE)
age_by_donor <- rep(45, opt$n_donors)
pcs <- matrix(0, nrow = opt$n_donors, ncol = 6)
pfs <- matrix(0, nrow = opt$n_donors, ncol = 2)
observed_total_counts <- pmax(round(library_size), 1L)
observed_log_total_counts <- log(observed_total_counts)

cell_dt <- data.table(
  CellID = sprintf("cell_%06d", seq_len(n_cells)),
  individual = cell_donor,
  sex = sex_by_donor[donor_index],
  age = age_by_donor[donor_index],
  latent_library_size = library_size,
  latent_log_library_size = log_library_size,
  total_read_counts = observed_total_counts,
  log_total_read_counts = observed_log_total_counts
)
for (i in 1:6) {
  cell_dt[[sprintf("pc%d", i)]] <- pcs[donor_index, i]
}
for (i in 1:2) {
  cell_dt[[sprintf("pf%d", i)]] <- pfs[donor_index, i]
}
cell_dt <- cbind(cell_dt, as.data.table(count_mat))

gene_info <- data.table(
  gene_name = gene_par$gene_name,
  chr_numeric = 1L,
  start = seq(1L, by = 1000L, length.out = opt$n_genes),
  end = seq(500L, by = 1000L, length.out = opt$n_genes)
)

meta_summary <- data.table(
  scenario = opt$scenario,
  replicate = opt$replicate,
  n_donors = opt$n_donors,
  n_cells = n_cells,
  n_genes = opt$n_genes,
  n_signal_pairs = nrow(signal_pairs),
  mean_latent_library_size = mean(library_size),
  mean_observed_total_counts = mean(observed_total_counts),
  sigma_zero = opt$sigma_zero,
  sigma_count = opt$sigma_count,
  shared_zero_lib_coef = opt$shared_zero_lib_coef,
  gamma_zero = gammas$gamma_zero,
  gamma_count = gammas$gamma_count
)

fwrite(cell_dt, file.path(rep_dir, "sim_counts.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(gene_info, file.path(rep_dir, "gene_info.tsv"), sep = "\t")
fwrite(signal_pairs, file.path(rep_dir, "truth_pairs.tsv"), sep = "\t")
fwrite(meta_summary, file.path(rep_dir, "simulation_metadata.tsv"), sep = "\t")

message("Saved replicate to: ", rep_dir)
