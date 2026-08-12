#!/usr/bin/env Rscript

# Calculate optimal chunking for SC hurdle to keep jobs under 6 hours

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

GENE_INFO_FILE <- cfg$gene_info_file
FILTER_FILE <- file.path(cfg$method2_results, "filtered_genes.tsv")

MAX_HOURS <- 5.5  # Target max time per job
SECS_PER_PAIR <- 0.5  # Conservative estimate (0.2-1.0 range)
CHUNK_SIZE_GENES <- 50  # Process ~50 genes per chunk

cat("=== Calculating optimal chunks for SC hurdle ===\n\n")

# Load gene info
gene_info <- fread(GENE_INFO_FILE)

# Load filtered genes if available
if (file.exists(FILTER_FILE)) {
  filt <- fread(FILTER_FILE)
  gene_info <- merge(gene_info, filt[keep == TRUE, .(gene_name)], by = "gene_name")
  cat("Using filtered genes (sparsity cutoff applied)\n")
} else {
  cat("Filtered genes not available, using 60% estimate\n")
}

# Calculate chunks per chromosome
chunk_plan <- data.table()

for (chr in 1:22) {
  chr_genes <- gene_info[chr_numeric == chr, gene_name]
  n_genes <- length(chr_genes)

  if (n_genes == 0) next

  # Calculate number of chunks needed
  n_chunks <- ceiling(n_genes / CHUNK_SIZE_GENES)
  genes_per_chunk <- ceiling(n_genes / n_chunks)

  # Estimate time for one chunk
  pairs_per_chunk <- genes_per_chunk * (genes_per_chunk - 1) / 2
  est_hours <- (pairs_per_chunk * SECS_PER_PAIR) / 3600

  cat(sprintf("Chr%2d: %4d genes → %2d chunks (%2d genes/chunk, ~%.1f hours/chunk)\n",
              chr, n_genes, n_chunks, genes_per_chunk, est_hours))

  # Generate chunk boundaries
  for (i in 1:n_chunks) {
    start_idx <- (i - 1) * genes_per_chunk + 1
    end_idx <- min(i * genes_per_chunk, n_genes)

    chunk_plan <- rbind(chunk_plan, data.table(
      chromosome = chr,
      chunk_id = i,
      total_chunks = n_chunks,
      gene_start = start_idx,
      gene_end = end_idx,
      n_genes_chunk = end_idx - start_idx + 1
    ))
  }
}

# Save chunk plan
out_file <- file.path(cfg$work_dir, "chunk_plan.tsv")
fwrite(chunk_plan, out_file, sep = "\t")

cat(sprintf("\nTotal chunks to submit: %d\n", nrow(chunk_plan)))
cat(sprintf("Chunk plan saved to: %s\n", out_file))
