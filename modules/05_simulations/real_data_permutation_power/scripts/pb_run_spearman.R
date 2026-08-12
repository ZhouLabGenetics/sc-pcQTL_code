#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--pb_file", help = "RDS with list(counts, genes)"),
  make_option("--output_dir", default = "results/power_full/pb"),
  make_option("--method", default = "spearman")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$pb_file) || !file.exists(opt$pb_file)) {
  stop("Must provide --pb_file")
}
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

message("=== Running PB Spearman correlations ===")
pb_obj <- readRDS(opt$pb_file)
counts <- pb_obj$counts
genes <- pb_obj$genes
if (is.null(counts) || is.null(genes)) stop("counts/genes missing in pb_file")

n_genes <- length(genes)
pval_mat <- matrix(NA_real_, n_genes, n_genes, dimnames = list(genes, genes))
cor_mat <- matrix(NA_real_, n_genes, n_genes, dimnames = list(genes, genes))

for (i in seq_len(n_genes - 1)) {
  for (j in (i + 1):n_genes) {
    v1 <- counts[, i]
    v2 <- counts[, j]
    test <- suppressWarnings(cor.test(v1, v2, method = opt$method, exact = FALSE))
    pval_mat[i, j] <- pval_mat[j, i] <- test$p.value
    cor_mat[i, j] <- cor_mat[j, i] <- test$estimate
  }
}

diag(pval_mat) <- NA_real_
diag(cor_mat) <- 1

fwrite(as.data.table(pval_mat, keep.rownames = "Gene"),
       file.path(opt$output_dir, "pvalue_matrix.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(as.data.table(cor_mat, keep.rownames = "Gene"),
       file.path(opt$output_dir, "correlation_matrix.tsv.gz"),
       sep = "\t", compress = "gzip")

sink(file.path(opt$output_dir, "summary_report.txt"))
cat("=== PB Spearman summary ===\n")
cat(sprintf("Input file: %s\n", opt$pb_file))
cat(sprintf("Genes: %d  Samples: %d\n", n_genes, nrow(counts)))
sink()

message("PB outputs written to ", opt$output_dir)
