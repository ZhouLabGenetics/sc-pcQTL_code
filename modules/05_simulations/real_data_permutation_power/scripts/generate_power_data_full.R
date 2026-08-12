#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
sim_code_root <- Sys.getenv(
  "COQTL_SIM_CODE_ROOT",
  unset = normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
)
source(file.path(sim_code_root, "shared", "sim_paths.R"))

option_list <- list(
  make_option("--raw_counts",
              default = get_raw_counts_file()),
  make_option("--gene_info",
              default = get_gene_info_file()),
  make_option("--output", default = file.path(get_data_root(), "power_full"),
              help = "Base output directory [default %default]"),
  make_option("--config_label", default = "full_realcounts",
              help = "Subdirectory name under output"),
  make_option("--n_genes", type = "integer", default = 200),
  make_option("--n_modules", type = "integer", default = 4),
  make_option("--module_size", type = "integer", default = 20),
  make_option("--effect_strengths", default = "0.02,0.04,0.06"),
  make_option("--seed", type = "integer", default = 20240121)
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

effect_strengths <- as.numeric(strsplit(opt$effect_strengths, ",")[[1]])
effect_strengths <- effect_strengths[!is.na(effect_strengths)]
if (!length(effect_strengths)) stop("At least one effect strength needed")

base_dir <- file.path(opt$output, opt$config_label)
null_dir <- file.path(base_dir, "null")
signal_dir <- file.path(base_dir, "signal")
if (dir.exists(base_dir) && length(list.files(base_dir, all.files = TRUE, no.. = TRUE))) {
  stop(
    "Refusing to overwrite an existing simulation dataset: ", base_dir,
    ". Use a new --config_label or archive/remove the old run explicitly."
  )
}
dir.create(null_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(signal_dir, recursive = TRUE, showWarnings = FALSE)

message("=== power_full data generator ===")
message(sprintf("Output: %s", base_dir))

message("Loading gene info...")
gene_info <- fread(opt$gene_info, select = "gene_name")
if (nrow(gene_info) < opt$n_genes) stop("n_genes exceeds available genes")
genes <- sample(gene_info$gene_name, opt$n_genes)

message("Loading entire raw count matrix for selected genes...")
counts_dt <- fread(opt$raw_counts, select = c("CellID", genes))
counts_dt[, DonorID := tstrsplit(CellID, "_", keep = 1)]
setcolorder(counts_dt, c("CellID", "DonorID", genes))
cell_meta <- counts_dt[, .(CellID, DonorID)]
sc_mat <- as.matrix(counts_dt[, ..genes])
rownames(sc_mat) <- cell_meta$CellID
n_cells <- nrow(sc_mat)
n_donors <- length(unique(cell_meta$DonorID))
message(sprintf("Retained donors: %d   cells: %d", n_donors, n_cells))

permute_matrix <- function(mat) {
  res <- mat
  for (j in seq_len(ncol(mat))) res[, j] <- sample(res[, j])
  res
}

message("Building strict null via per-gene permutation...")
sc_null <- permute_matrix(sc_mat)

pb_null_dt <- data.table(DonorID = cell_meta$DonorID, sc_null)
pb_null <- pb_null_dt[, lapply(.SD, mean), by = DonorID, .SDcols = genes]
pb_null_mat <- as.matrix(pb_null[, ..genes])
rownames(pb_null_mat) <- pb_null$DonorID

saveRDS(list(
  counts = sc_null,
  donors = cell_meta$DonorID,
  cells = cell_meta$CellID,
  genes = genes
), file.path(null_dir, "sc_counts_null.rds"))
saveRDS(list(
  counts = pb_null_mat,
  donors = pb_null$DonorID,
  genes = genes
), file.path(null_dir, "pb_counts_null.rds"))
fwrite(cell_meta, file.path(base_dir, "cell_metadata.tsv"), sep = "\t")

message("Sampling signal modules...")
module_genes <- vector("list", opt$n_modules)
available <- genes
for (m in seq_len(opt$n_modules)) {
  if (length(available) < opt$module_size) available <- genes
  module_genes[[m]] <- sample(available, opt$module_size)
  available <- setdiff(available, module_genes[[m]])
}
module_assign <- rbindlist(lapply(seq_along(module_genes), function(idx) {
  data.table(Gene = module_genes[[idx]], Module = paste0("M", idx))
}))
fwrite(module_assign, file.path(base_dir, "module_assignments.tsv"), sep = "\t")

pair_idx <- combn(genes, 2)
pair_dt <- data.table(
  Gene1 = pair_idx[1, ],
  Gene2 = pair_idx[2, ],
  Label = "null"
)
module_lookup <- module_assign[, .(Module), keyby = Gene]
pair_dt[, Module1 := module_lookup[Gene1, Module]]
pair_dt[, Module2 := module_lookup[Gene2, Module]]
pair_dt[!is.na(Module1) & Module1 == Module2, Label := "signal"]
pair_dt[, c("Module1", "Module2") := NULL]
fwrite(pair_dt, file.path(base_dir, "truth_pairs.tsv"), sep = "\t")

simulate_signal <- function(base_counts, strength) {
  latent <- rnorm(n_cells)
  mat <- base_counts
  for (genes_in_mod in module_genes) {
    for (g in genes_in_mod) {
      idx <- match(g, genes)
      base_vec <- base_counts[, idx]
      log_mu <- log1p(base_vec) + strength * latent
      mu <- pmax(exp(log_mu) - 1, 0)
      mat[, idx] <- rpois(n_cells, lambda = mu)
    }
    latent <- rnorm(n_cells)
  }
  mat
}

message("Generating signal datasets...")
effect_labels <- sprintf("strength_s%03d", round(effect_strengths * 100))
for (k in seq_along(effect_strengths)) {
  strength <- effect_strengths[k]
  label <- effect_labels[k]
  dir_k <- file.path(signal_dir, label)
  dir.create(dir_k, recursive = TRUE, showWarnings = FALSE)
  message(sprintf("  Strength %.3f -> %s", strength, dir_k))
  sc_sig <- simulate_signal(sc_null, strength)
  pb_sig_dt <- data.table(DonorID = cell_meta$DonorID, sc_sig)
  pb_sig <- pb_sig_dt[, lapply(.SD, mean), by = DonorID, .SDcols = genes]
  pb_sig_mat <- as.matrix(pb_sig[, ..genes])
  rownames(pb_sig_mat) <- pb_sig$DonorID

  saveRDS(list(
    counts = sc_sig,
    donors = cell_meta$DonorID,
    cells = cell_meta$CellID,
    genes = genes,
    effect_strength = strength
  ), file.path(dir_k, "sc_counts_signal.rds"))
  saveRDS(list(
    counts = pb_sig_mat,
    donors = pb_sig$DonorID,
    genes = genes,
    effect_strength = strength
  ), file.path(dir_k, "pb_counts_signal.rds"))
}

meta_dt <- data.table(
  config_label = opt$config_label,
  n_genes = opt$n_genes,
  n_donors = n_donors,
  n_cells = n_cells,
  n_modules = opt$n_modules,
  module_size = opt$module_size,
  effect_strength = effect_strengths,
  effect_label = effect_labels
)
fwrite(meta_dt, file.path(base_dir, "config_summary.tsv"), sep = "\t")

message("All datasets written to ", base_dir)
