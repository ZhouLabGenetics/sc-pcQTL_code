#!/usr/bin/env Rscript

# Build the model-based calibration and power source tables consumed by the
# final Figure 2 compositor. Figure styling belongs to module 11.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--sim_root", type = "character", default = "simulations_large"),
  make_option("--analysis_root", type = "character", default = "analysis_large_compare"),
  make_option("--out_dir", type = "character", default = "plots_target"),
  make_option("--chunk_name", type = "character", default = "chr1_chunk001"),
  make_option("--alpha_grid", type = "character", default = "0.05,0.01,0.005,0.001")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (!dir.exists(opt$sim_root)) stop("Missing simulation root: ", opt$sim_root)
if (!dir.exists(opt$analysis_root)) stop("Missing analysis root: ", opt$analysis_root)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

scenario_levels <- c(
  "null",
  "zero_only_low", "count_only_low", "both_low",
  "zero_only_high", "count_only_high", "both_high"
)
method_levels <- c(
  "Hurdle Count", "Hurdle Zero", "Hurdle Union",
  "Donor Pearson", "Donor Spearman"
)
threshold_levels <- c("0.05", "0.01", "0.005", "0.001")

alpha_grid <- as.numeric(strsplit(opt$alpha_grid, ",", fixed = TRUE)[[1]])
alpha_grid <- unique(alpha_grid[is.finite(alpha_grid) & alpha_grid > 0 & alpha_grid <= 1])
alpha_grid <- alpha_grid[order(-alpha_grid)]
if (!length(alpha_grid)) stop("--alpha_grid contains no valid thresholds")

label_threshold <- function(x) {
  out <- rep(NA_character_, length(x))
  out[abs(x - 0.05) < 1e-12] <- "0.05"
  out[abs(x - 0.01) < 1e-12] <- "0.01"
  out[abs(x - 0.005) < 1e-12] <- "0.005"
  out[abs(x - 0.001) < 1e-12] <- "0.001"
  out
}

method_component_label <- function(method, component) {
  if (method == "hurdle") {
    paste("Hurdle", tools::toTitleCase(component))
  } else {
    paste("Donor", tools::toTitleCase(component))
  }
}

extract_pair_table <- function(mat_dt) {
  genes <- mat_dt[[1]]
  mat <- as.matrix(mat_dt[, -1, with = FALSE])
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    return(data.table(pair_key = character(), pvalue = numeric()))
  }
  idx <- which(upper.tri(mat, diag = FALSE), arr.ind = TRUE)
  pvals <- mat[idx]
  keep <- is.finite(pvals) & pvals > 0 & pvals <= 1
  idx <- idx[keep, , drop = FALSE]
  data.table(
    pair_key = paste(
      pmin(genes[idx[, 1]], genes[idx[, 2]]),
      pmax(genes[idx[, 1]], genes[idx[, 2]]),
      sep = "||"
    ),
    pvalue = pvals[keep]
  )
}

extract_pair_table_directional <- function(mat_dt) {
  genes <- mat_dt[[1]]
  mat <- as.matrix(mat_dt[, -1, with = FALSE])
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    return(data.table(pair_key = character(), pvalue = numeric()))
  }
  idx <- which(row(mat) != col(mat), arr.ind = TRUE)
  pvals <- mat[idx]
  keep <- is.finite(pvals) & pvals > 0 & pvals <= 1
  idx <- idx[keep, , drop = FALSE]
  data.table(
    pair_key = paste(genes[idx[, 1]], genes[idx[, 2]], sep = "||"),
    pvalue = pvals[keep]
  )
}

get_tested_genes_from_matrix <- function(path) {
  if (!file.exists(path)) return(character())
  setdiff(names(fread(path, nrows = 0)), "Gene")
}

load_method_tables <- function(scenario, replicate) {
  sim_dir <- file.path(opt$sim_root, scenario, replicate)
  truth_file <- file.path(sim_dir, "truth_pairs.tsv")
  truth_dt <- if (file.exists(truth_file)) fread(truth_file) else data.table()
  if (!all(c("gene1", "gene2") %in% names(truth_dt))) {
    stop("Missing gene1/gene2 truth columns: ", truth_file)
  }

  out <- list()
  hurdle_dir <- file.path(
    opt$analysis_root, "hurdle", scenario, replicate,
    "results", "method2_sc_hurdle", "gene_associations_chunked", opt$chunk_name
  )
  count_file <- file.path(hurdle_dir, "pvalues_count.tsv.gz")
  zero_file <- file.path(hurdle_dir, "pvalues_zero.tsv.gz")
  if (file.exists(count_file) && file.exists(zero_file)) {
    hurdle_genes <- intersect(
      get_tested_genes_from_matrix(count_file),
      get_tested_genes_from_matrix(zero_file)
    )
    truth_sub <- truth_dt[gene1 %in% hurdle_genes & gene2 %in% hurdle_genes]
    truth_keys <- unique(c(
      paste(truth_sub$gene1, truth_sub$gene2, sep = "||"),
      paste(truth_sub$gene2, truth_sub$gene1, sep = "||")
    ))
    merged <- merge(
      extract_pair_table_directional(fread(count_file)),
      extract_pair_table_directional(fread(zero_file)),
      by = "pair_key", all = TRUE, suffixes = c("_count", "_zero")
    )
    merged[, truth := pair_key %in% truth_keys]
    merged[, pvalue_union := pmin(pvalue_count, pvalue_zero, na.rm = TRUE)]
    merged[!is.finite(pvalue_union), pvalue_union := NA_real_]
    out[[1]] <- list(component = "count", dt = merged[, .(pvalue = pvalue_count, truth)])
    out[[2]] <- list(component = "zero", dt = merged[, .(pvalue = pvalue_zero, truth)])
    out[[3]] <- list(component = "union", dt = merged[, .(pvalue = pvalue_union, truth)])
  }

  donor_dir <- file.path(
    opt$analysis_root, "donor_agg", scenario, replicate, "results", "donor_agg"
  )
  pearson_file <- file.path(donor_dir, "pvalues_pearson.tsv.gz")
  spearman_file <- file.path(donor_dir, "pvalues_spearman.tsv.gz")
  if (file.exists(pearson_file) && file.exists(spearman_file)) {
    donor_genes <- intersect(
      get_tested_genes_from_matrix(pearson_file),
      get_tested_genes_from_matrix(spearman_file)
    )
    truth_sub <- truth_dt[gene1 %in% donor_genes & gene2 %in% donor_genes]
    truth_keys <- unique(paste(
      pmin(truth_sub$gene1, truth_sub$gene2),
      pmax(truth_sub$gene1, truth_sub$gene2),
      sep = "||"
    ))
    pearson <- extract_pair_table(fread(pearson_file))
    spearman <- extract_pair_table(fread(spearman_file))
    pearson[, truth := pair_key %in% truth_keys]
    spearman[, truth := pair_key %in% truth_keys]
    out[[4]] <- list(component = "pearson", dt = pearson[, .(pvalue, truth)])
    out[[5]] <- list(component = "spearman", dt = spearman[, .(pvalue, truth)])
  }
  out
}

sim_dirs <- list.dirs(opt$sim_root, recursive = TRUE, full.names = TRUE)
sim_dirs <- sim_dirs[grepl("replicate_[0-9]+$", sim_dirs)]
records <- list()

for (sim_dir in sim_dirs) {
  rel_path <- sub(
    paste0("^", normalizePath(opt$sim_root), "/?"), "", normalizePath(sim_dir)
  )
  parts <- strsplit(rel_path, "/", fixed = TRUE)[[1]]
  if (length(parts) < 2 || !parts[1] %in% scenario_levels) next
  scenario <- parts[1]
  replicate <- parts[2]

  tables <- load_method_tables(scenario, replicate)
  for (entry in tables) {
    dt <- entry$dt
    if (!nrow(dt)) next
    method <- if (entry$component %in% c("count", "zero", "union")) "hurdle" else "donor"
    label <- method_component_label(method, entry$component)
    n_signal <- sum(dt$truth, na.rm = TRUE)
    n_null <- sum(!dt$truth, na.rm = TRUE)
    for (alpha in alpha_grid) {
      detected <- is.finite(dt$pvalue) & dt$pvalue < alpha
      records[[length(records) + 1L]] <- data.table(
        scenario = scenario,
        replicate = replicate,
        method_component = label,
        p_threshold = alpha,
        power = if (n_signal > 0) sum(detected & dt$truth) / n_signal else NA_real_,
        fpr = if (n_null > 0) sum(detected & !dt$truth) / n_null else NA_real_
      )
    }
  }
}

curve_dt <- rbindlist(records, fill = TRUE)
if (!nrow(curve_dt)) stop("No model-simulation source data were generated")
curve_dt[, scenario := factor(scenario, levels = scenario_levels)]
curve_dt[, method_component := factor(method_component, levels = method_levels)]
curve_dt[, threshold_label := factor(label_threshold(p_threshold), levels = threshold_levels)]
curve_dt <- curve_dt[!is.na(threshold_label)]

curve_summary <- curve_dt[, .(
  mean_power = mean(power, na.rm = TRUE),
  mean_fpr = mean(fpr, na.rm = TRUE)
), by = .(scenario, method_component, p_threshold, threshold_label)]

fwrite(curve_dt, file.path(opt$out_dir, "alpha_curve_replicate_metrics.tsv"), sep = "\t")
fwrite(curve_summary, file.path(opt$out_dir, "alpha_curve_summary.tsv"), sep = "\t")
message("Saved model-simulation source tables to: ", opt$out_dir)
