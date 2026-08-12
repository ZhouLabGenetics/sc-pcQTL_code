#!/usr/bin/env Rscript

.local_lib <- Sys.getenv("SC_PCQTL_R_LIBS", unset = "")
if (nzchar(.local_lib)) .libPaths(c(.local_lib, .libPaths()))

suppressPackageStartupMessages(library(data.table))

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args_all[grep("--file=", args_all)[1]]))
module_dir <- dirname(script_path)
source(file.path(module_dir, "load_config.R"))
cfg <- load_joint_score_config(module_dir)
setDTthreads(1L)

dir.create(cfg$method_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$stage_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cfg$primary_filtered_genes_file)) {
  stop("Missing primary filtered-gene file: ", cfg$primary_filtered_genes_file)
}
filter_target <- file.path(cfg$method_dir, "filtered_genes.tsv")
source_md5 <- unname(tools::md5sum(cfg$primary_filtered_genes_file))
if (file.exists(filter_target)) {
  target_md5 <- unname(tools::md5sum(filter_target))
  if (!identical(source_md5, target_md5)) {
    stop("Existing staged filtered_genes.tsv differs from the primary input")
  }
} else if (!file.copy(cfg$primary_filtered_genes_file, filter_target)) {
  stop("Failed to stage primary filtered-gene file")
}

gene_info <- fread(
  cfg$gene_info_file,
  select = c("gene_name", "chr_numeric", "start", "end")
)
filtered <- fread(filter_target)
if (!all(c("gene_name", "keep") %in% names(filtered))) {
  stop("filtered_genes.tsv must contain gene_name and keep columns")
}
filtered <- filtered[keep == TRUE]
if (!"nonzero_prop" %in% names(filtered)) filtered[, nonzero_prop := NA_real_]

genes <- merge(
  gene_info[chr_numeric %in% 1:22],
  filtered[, .(gene_name, nonzero_prop)],
  by = "gene_name"
)
genes <- unique(genes, by = "gene_name")
setorder(genes, chr_numeric, start, end, gene_name)
genes[, chr_rank := seq_len(.N), by = chr_numeric]
if (!nrow(genes)) stop("No autosomal genes passed the primary nonzero filter")

header <- names(fread(cfg$count_file, nrows = 0L))
id_col <- intersect(c("barcode", "CellID"), header)[1]
donor_col <- intersect(c("individual", "IndividualID"), header)[1]
if (is.na(id_col) || is.na(donor_col)) stop("Missing cell or donor identifier")

base_covars <- c("sex", "age", sprintf("pc%d", 1:6), "pf1", "pf2")
missing_covars <- setdiff(base_covars, header)
if (length(missing_covars)) stop("Missing covariates: ", paste(missing_covars, collapse = ", "))
missing_genes <- setdiff(genes$gene_name, header)
if (length(missing_genes)) {
  stop("Filtered genes absent from count matrix: ", paste(head(missing_genes, 20), collapse = ", "))
}

embedded_library <- all(c("total_read_counts", "log_total_read_counts") %in% header)
cache_file <- cfg$library_size_cache
if (!nzchar(cache_file)) cache_file <- file.path(cfg$method_dir, "library_size_cache.tsv.gz")
cached_library <- file.exists(cache_file)
metadata_cols <- unique(c(
  "CellID", "barcode", "IndividualID", "individual", "CellType",
  base_covars, "total_read_counts", "log_total_read_counts"
))
all_gene_cols <- setdiff(header, metadata_cols)
read_genes <- if (embedded_library || cached_library) genes$gene_name else all_gene_cols
select_cols <- unique(c(
  id_col, donor_col, base_covars,
  if (embedded_library) c("total_read_counts", "log_total_read_counts"),
  read_genes
))

message(sprintf("[%s] Reading staged matrix from %s", cfg$celltype, cfg$count_file))
dt <- fread(cfg$count_file, select = select_cols, showProgress = TRUE)
if (nrow(dt) < cfg$min_celltype_cells) {
  stop("Cell count below the central primary-analysis threshold")
}

if (embedded_library) {
  log_library <- as.numeric(dt$log_total_read_counts)
  library_source <- "embedded_input_columns"
} else if (cached_library) {
  lib <- fread(cache_file)
  lib_id <- intersect(c(id_col, "cell_id", "barcode", "CellID"), names(lib))[1]
  if (is.na(lib_id)) stop("Library-size cache has no cell identifier")
  idx <- match(dt[[id_col]], lib[[lib_id]])
  if (anyNA(idx)) stop("Library-size cache does not cover every staged cell")
  log_library <- as.numeric(lib$log_total_read_counts[idx])
  library_source <- "existing_cache"
} else {
  message(sprintf("[%s] Computing library size from %d complete gene columns", cfg$celltype, length(all_gene_cols)))
  total_library <- rowSums(as.matrix(dt[, ..all_gene_cols]))
  log_library <- log(pmax(total_library, 1))
  library_source <- "computed_from_complete_input_matrix"
  fwrite(
    data.table(
      cell_id = dt[[id_col]],
      total_read_counts = total_library,
      log_total_read_counts = log_library
    ),
    file.path(cfg$method_dir, "library_size_cache.tsv.gz"),
    sep = "\t",
    compress = "gzip"
  )
  rm(total_library)
}

covars <- dt[, c(id_col, donor_col, base_covars), with = FALSE]
setnames(covars, c(id_col, donor_col), c("cell_id", "individual"))
covars[, log_total_read_counts := log_library]
covars[, age := as.numeric(age)]
for (nm in c(sprintf("pc%d", 1:6), "pf1", "pf2")) {
  set(covars, j = nm, value = as.numeric(covars[[nm]]))
}
covars[, sex := as.factor(sex)]
model_covars <- c("sex", "age", sprintf("pc%d", 1:6), "pf1", "pf2", "log_total_read_counts")
if (any(!complete.cases(covars[, ..model_covars]))) stop("Missing values in model covariates")

saveRDS(covars, file.path(cfg$stage_dir, "covariates.rds"), compress = "gzip", version = 3)
fwrite(genes, file.path(cfg$stage_dir, "filtered_autosomal_genes.tsv"), sep = "\t")

for (chr in 1:22) {
  chr_genes <- genes[chr_numeric == chr, gene_name]
  if (!length(chr_genes)) next
  mat <- as.matrix(dt[, ..chr_genes])
  if (!all(is.finite(mat))) stop("Non-finite expression value on chromosome ", chr)
  storage.mode(mat) <- "integer"
  saveRDS(mat, file.path(cfg$stage_dir, sprintf("chr%d_counts.rds", chr)), compress = FALSE, version = 3)
  rm(mat)
  gc(verbose = FALSE)
}

chr_counts <- genes[, .(n_genes = .N), by = chr_numeric]
m_unordered <- chr_counts[, sum(n_genes * (n_genes - 1) / 2)]
local_unordered <- sum(
  genes[, sum(pmax(0L, pmin(49L, .N - seq_len(.N)))), by = chr_numeric]$V1
)
summary <- data.table(
  celltype = cfg$celltype,
  n_cells = nrow(covars),
  n_filtered_autosomal_genes = nrow(genes),
  m_unordered_all_chr_pairs = m_unordered,
  n_directional_tests_denominator = 2 * m_unordered,
  n_local_unordered_pairs_computed = local_unordered,
  n_local_directional_tests_computed = 2 * local_unordered,
  bonferroni_joint_threshold = 0.05 / (2 * m_unordered),
  library_size_source = library_source,
  primary_filtered_genes_md5 = source_md5
)
fwrite(summary, file.path(cfg$stage_dir, "stage_summary.tsv"), sep = "\t")
writeLines("OK", file.path(cfg$stage_dir, "STAGE_COMPLETE"))
message(sprintf(
  "[%s] Staged %d cells; M=%s; local pairs=%s",
  cfg$celltype, nrow(covars), format(m_unordered, scientific = FALSE),
  format(local_unordered, scientific = FALSE)
))
