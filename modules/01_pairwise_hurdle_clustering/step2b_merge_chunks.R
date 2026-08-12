#!/usr/bin/env Rscript

# =============================================================================
# Method 2 - Step 2b: Merge chunks for each chromosome
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
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

option_list <- list(
  make_option("--chr", type = "integer", help = "Chromosome to merge (or 0 for all)")
)

opt <- parse_args(OptionParser(option_list = option_list))
CHR_TO_MERGE <- if (is.null(opt$chr)) 0 else opt$chr

cat(sprintf("\n=== Method 2 Step 2b: Merge Chunks (%s) ===\n\n", cfg$cell_type))

CHUNK_DIR <- cfg$chunk_dir
OUTPUT_DIR <- cfg$assoc_dir
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Determine which chromosomes to merge
if (CHR_TO_MERGE == 0) {
  chrs_to_process <- 1:22
  cat("Merging all chromosomes 1-22\n")
} else {
  chrs_to_process <- CHR_TO_MERGE
  cat(sprintf("Merging chromosome %d\n", CHR_TO_MERGE))
}

for (chr in chrs_to_process) {
  cat(sprintf("\nChromosome %d:\n", chr))

  # Find all chunks for this chromosome
  chunk_dirs <- list.dirs(CHUNK_DIR, recursive = FALSE)
  chr_chunks <- grep(sprintf("chr%d_chunk", chr), basename(chunk_dirs), value = TRUE)

  if (length(chr_chunks) == 0) {
    cat("  No chunks found, skipping\n")
    next
  }

  cat(sprintf("  Found %d chunks\n", length(chr_chunks)))

  # Load all significant pairs
  all_sig_pairs <- list()

  for (chunk_name in chr_chunks) {
    chunk_path <- file.path(CHUNK_DIR, chunk_name)
    sig_file <- file.path(chunk_path, "significant_pairs.tsv")

    if (file.exists(sig_file)) {
      sig_data <- fread(sig_file)
      all_sig_pairs[[chunk_name]] <- sig_data
    }
  }

  if (length(all_sig_pairs) == 0) {
    cat("  No significant pairs found\n")
    next
  }

  # Combine all significant pairs
  combined_sig <- rbindlist(all_sig_pairs, fill = TRUE)

  # Remove duplicates (same pair might appear in multiple chunks)
  combined_sig[, pair_id := paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_")]
  combined_sig <- combined_sig[!duplicated(pair_id)]
  combined_sig[, pair_id := NULL]

  cat(sprintf("  Total unique significant pairs: %d\n", nrow(combined_sig)))

  # Save merged results
  out_base <- file.path(OUTPUT_DIR, sprintf("chr%d", chr))
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

  fwrite(combined_sig,
         file.path(out_base, sprintf("chr%d_significant_pairs.tsv", chr)),
         sep = "\t")

  # Save summary
  summary_dt <- data.table(
    Chromosome = chr,
    NumChunks = length(chr_chunks),
    SignificantPairs = nrow(combined_sig)
  )
  fwrite(summary_dt,
         file.path(out_base, sprintf("chr%d_summary.tsv", chr)),
         sep = "\t")

  cat(sprintf("  Saved to: %s\n", out_base))
}

cat("\n=== Complete ===\n")
cat(sprintf("Merged results in: %s\n", OUTPUT_DIR))
cat("\nNext step: bash submit_step3.sh\n")
