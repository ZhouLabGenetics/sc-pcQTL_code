#!/usr/bin/env Rscript
# ============================================================================
# 02_identify_qtls.R
# Per-phenotype FDR correction and identify significant QTLs
# ============================================================================

.local_libs <- c(
  Sys.getenv("SC_PCQTL_R_LIBS", unset = ""),
  Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = "")
)
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# Setup
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  normalizePath(getwd())
}
MODULE_DIR <- get_script_dir()
SUMMARY_ROOT <- Sys.getenv("SC_PCQTL_PCQTL_SUMMARY_ROOT", unset = "")
if (!nzchar(SUMMARY_ROOT)) stop("Set SC_PCQTL_PCQTL_SUMMARY_ROOT.")
DATA_DIR <- file.path(SUMMARY_ROOT, "data")

FDR_THRESHOLD <- 0.05

cat("=== Identifying significant QTLs ===\n")
cat("FDR threshold:", FDR_THRESHOLD, "\n")
cat("Method: per phenotype (cluster × PC) BH using SNP-level tests\n\n")

# Read results
in_file <- file.path(DATA_DIR, "all_pcqtl_results.tsv")
if (!file.exists(in_file)) {
  stop("ERROR: all_pcqtl_results.tsv not found. Run 01_collect_results.R first.")
}

all_data <- fread(in_file)
cat("Loaded", nrow(all_data), "phenotype-level entries\n\n")

required_cols <- c("n_snps", "min_snp_fdr", "has_sig_snp")
missing_cols <- setdiff(required_cols, names(all_data))
if (length(missing_cols)) {
  stop("Missing columns in all_pcqtl_results.tsv: ", paste(missing_cols, collapse = ", "),
       ". Re-run 01_collect_results.R to include SNP summaries.")
}

# Phenotype-level significance (already BH within each cluster×PC)
all_data[, fdr := min_snp_fdr]
all_data[, is_sig := !is.na(fdr) & fdr < FDR_THRESHOLD]

# Summary
cat("\nResults per celltype:\n")
summary_ct <- all_data[, .(
  n_tests = .N,
  n_sig = sum(is_sig, na.rm = TRUE),
  pct_sig = round(100 * sum(is_sig, na.rm = TRUE) / .N, 2),
  min_p = suppressWarnings(min(min_snp_p, na.rm = TRUE)),
  median_p = median(min_snp_p, na.rm = TRUE)
), by = celltype][order(-n_sig)]

print(summary_ct)

# Identify significant QTLs
sig_qtls <- all_data[is_sig == TRUE]
cat("\n=== Significant QTLs ===\n")
cat("Total significant phenotypes:", nrow(sig_qtls), "\n")
cat("Across", uniqueN(sig_qtls$celltype), "celltypes\n")
cat("Unique clusters with QTLs:", uniqueN(sig_qtls$cluster_id), "\n")

# QTL details
cat("\nTop 10 strongest QTLs (by min SNP FDR):\n")
top_qtls <- sig_qtls[order(fdr)][1:10, .(
  celltype, cluster_id, PC, chr, genes,
  lead_snp_marker, lead_snp_chr, lead_snp_pos,
  min_snp_p, fdr, n_sig_snps
)]
print(top_qtls)

# Save results
cat("\n=== Saving results ===\n")

# All results with phenotype-level FDR
out_file1 <- file.path(DATA_DIR, "all_results_with_fdr.tsv")
fwrite(all_data, out_file1, sep = "\t")
cat("✓ Saved all results:", out_file1, "\n")

# Significant QTLs only
out_file2 <- file.path(DATA_DIR, "sig_qtls.tsv")
fwrite(sig_qtls, out_file2, sep = "\t")
cat("✓ Saved significant QTLs:", out_file2, "\n")

# Per-celltype summary
out_file3 <- file.path(DATA_DIR, "qtl_counts_per_celltype.tsv")
fwrite(summary_ct, out_file3, sep = "\t")
cat("✓ Saved summary:", out_file3, "\n")

cat("\n✓ Step 2 complete!\n")
