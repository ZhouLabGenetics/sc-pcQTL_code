#!/usr/bin/env Rscript

# =============================================================================
# Method 2 - Step 1: Filter sparse genes using 1% nonzero cutoff
# =============================================================================

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

cat(sprintf("\n=== Method 2 Step 1: Filter sparse genes (%s) ===\n", cfg$cell_type))
cat(sprintf("Start time: %s\n\n", Sys.time()))

COUNT_FILE <- cfg$count_file
GENE_INFO_FILE <- cfg$gene_info_file
OUTPUT_DIR <- cfg$method2_results
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

NONZERO_CUTOFF <- 0.01
GENE_BLOCK_SIZE <- 500

cat("Inspecting count matrix header...\n")
header_cols <- names(fread(COUNT_FILE, nrows = 0))
id_col <- intersect(c("CellID", "barcode"), header_cols)
if (!length(id_col)) stop("Count file must contain CellID or barcode column")
id_col <- id_col[1]
metadata_cols <- c(
  id_col, "CellID", "barcode", "IndividualID", "individual", "CellType",
  "sex", sprintf("pc%d", 1:6), "age", "pf1", "pf2",
  "total_read_counts", "log_total_read_counts"
)
genes <- setdiff(header_cols, metadata_cols)

cat("Loading one block to determine cell count...\n")
probe_cols <- unique(c(id_col, genes[seq_len(min(length(genes), GENE_BLOCK_SIZE))]))
probe_dt <- fread(COUNT_FILE, select = probe_cols)
n_cells <- nrow(probe_dt)
cat(sprintf("  Cells: %d  Genes: %d\n", n_cells, length(genes)))
if (n_cells < cfg$min_celltype_cells) {
  stop(sprintf(
    "Cell type %s has %d cells, below the primary-analysis threshold of %d; no hurdle tasks were generated.",
    cfg$cell_type, n_cells, cfg$min_celltype_cells
  ))
}

cat(sprintf("Calculating nonzero proportions in gene blocks of %d...\n", GENE_BLOCK_SIZE))
gene_blocks <- split(genes, ceiling(seq_along(genes) / GENE_BLOCK_SIZE))
nonzero_list <- vector("list", length(gene_blocks))

for (block_idx in seq_along(gene_blocks)) {
  gene_block <- gene_blocks[[block_idx]]
  if (block_idx == 1L) {
    block_dt <- probe_dt[, ..gene_block]
  } else {
    block_dt <- fread(COUNT_FILE, select = gene_block)
  }
  block_mat <- as.matrix(block_dt)
  nonzero_list[[block_idx]] <- colMeans(block_mat > 0)
  rm(block_dt, block_mat)
  gc(verbose = FALSE)

  if (block_idx %% 10 == 0 || block_idx == length(gene_blocks)) {
    cat(sprintf("  Processed %d / %d blocks\n", block_idx, length(gene_blocks)))
  }
}

rm(probe_dt)
gc(verbose = FALSE)

nonzero_prop <- unlist(nonzero_list, use.names = TRUE)
nonzero_prop <- nonzero_prop[genes]

nonzero_dt <- data.table(
  gene_name = genes,
  nonzero_prop = nonzero_prop
)

gene_info <- fread(GENE_INFO_FILE, select = c("gene_name", "chr_numeric", "start", "end"))
result_dt <- merge(nonzero_dt, gene_info, by = "gene_name", all.x = TRUE)
result_dt[, keep := nonzero_prop >= NONZERO_CUTOFF]

cat(sprintf("  Genes passing cutoff: %d / %d (%.1f%%)\n",
            sum(result_dt$keep), nrow(result_dt),
            100 * mean(result_dt$keep)))

out_file <- file.path(OUTPUT_DIR, "filtered_genes.tsv")
fwrite(result_dt, out_file, sep = "\t")

cat("\n=== Complete ===\n")
cat(sprintf("Saved: %s\n", out_file))
cat(sprintf("End time: %s\n", Sys.time()))
