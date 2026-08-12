#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

get_this_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_idx <- grep("--file=", cmd_args)
  if (length(file_idx) > 0) {
    return(normalizePath(sub("--file=", "", cmd_args[file_idx])))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
  normalizePath(".")
}

script_dir <- dirname(get_this_path())
source(file.path(script_dir, "load_config.R"))
cfg <- load_config(script_dir)

COUNT_FILE <- cfg$count_file
OUT_FILE <- file.path(cfg$method2_results, "library_size_cache.tsv.gz")
GENE_BLOCK_SIZE <- 500

cat(sprintf("\n=== Compute library size cache (%s) ===\n", cfg$cell_type))
cat(sprintf("Start time: %s\n\n", Sys.time()))

if (file.exists(OUT_FILE)) {
  cat(sprintf("Cache already exists: %s\n", OUT_FILE))
  quit(status = 0)
}

header_cols <- names(fread(COUNT_FILE, nrows = 0))
id_col <- intersect(c("CellID", "barcode"), header_cols)
if (!length(id_col)) stop("Count file must contain CellID or barcode column")
id_col <- id_col[1]

metadata_cols <- c(
  id_col, "CellID", "barcode", "IndividualID", "individual", "CellType",
  "sex", sprintf("pc%d", 1:6), "age", "pf1", "pf2",
  "total_read_counts", "log_total_read_counts"
)

if (all(c("total_read_counts", "log_total_read_counts") %in% header_cols)) {
  cat("Reading existing library size columns from count file...\n")
  cache_dt <- fread(COUNT_FILE, select = c(id_col, "total_read_counts", "log_total_read_counts"))
  fwrite(cache_dt, OUT_FILE, sep = "\t", compress = "gzip")
  cat(sprintf("Saved: %s\n", OUT_FILE))
  quit(status = 0)
}

gene_cols <- setdiff(header_cols, metadata_cols)
gene_blocks <- split(gene_cols, ceiling(seq_along(gene_cols) / GENE_BLOCK_SIZE))

cat(sprintf("Gene columns: %d in %d blocks\n", length(gene_cols), length(gene_blocks)))
block_sum <- function(dt) {
  block_mat <- as.matrix(dt)
  storage.mode(block_mat) <- "numeric"
  rowSums(block_mat, na.rm = TRUE)
}
probe_dt <- fread(COUNT_FILE, select = c(id_col, gene_blocks[[1]]))
total_counts <- block_sum(probe_dt[, ..gene_blocks[[1]]])
cache_dt <- probe_dt[, ..id_col]
rm(probe_dt)
gc(verbose = FALSE)

if (length(gene_blocks) > 1) {
  for (block_idx in 2:length(gene_blocks)) {
    block_dt <- fread(COUNT_FILE, select = gene_blocks[[block_idx]])
    total_counts <- total_counts + block_sum(block_dt)
    rm(block_dt)
    gc(verbose = FALSE)

    if (block_idx %% 10 == 0 || block_idx == length(gene_blocks)) {
      cat(sprintf("  Processed %d / %d blocks\n", block_idx, length(gene_blocks)))
    }
  }
}

cache_dt[, total_read_counts := total_counts]
cache_dt[, log_total_read_counts := log(pmax(total_read_counts, 1))]
fwrite(cache_dt, OUT_FILE, sep = "\t", compress = "gzip")

cat(sprintf("Saved: %s\n", OUT_FILE))
cat(sprintf("End time: %s\n", Sys.time()))
