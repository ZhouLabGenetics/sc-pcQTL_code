#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--sim_dir", type = "character"),
  make_option("--out_dir", type = "character"),
  make_option("--p_threshold", type = "double", default = 0.05),
  make_option("--nonzero_cutoff", type = "double", default = 0.01)
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$sim_dir) || is.null(opt$out_dir)) {
  stop("--sim_dir and --out_dir are required.")
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
result_dir <- file.path(opt$out_dir, "results", "donor_agg")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

count_file <- file.path(opt$sim_dir, "sim_counts.tsv.gz")
if (!file.exists(count_file)) {
  stop("Missing simulation count file: ", count_file)
}

count_dt <- fread(count_file)
meta_cols <- c(
  "CellID", "barcode", "IndividualID", "individual", "CellType",
  "sex", paste0("pc", 1:6), "age", "pf1", "pf2",
  "total_read_counts", "log_total_read_counts"
)
gene_cols <- setdiff(names(count_dt), meta_cols)
if (!"individual" %in% names(count_dt)) {
  stop("Expected 'individual' column in simulation counts.")
}

agg_dt <- count_dt[, lapply(.SD, mean), by = individual, .SDcols = gene_cols]
agg_mat <- as.matrix(agg_dt[, ..gene_cols])
nonzero_prop <- colMeans(agg_mat > 0)
keep_genes <- gene_cols[nonzero_prop >= opt$nonzero_cutoff]
if (length(keep_genes) < 2) {
  stop("Fewer than 2 genes remain after donor-level nonzero filter.")
}

agg_mat <- as.matrix(agg_dt[, ..keep_genes])
n_genes <- ncol(agg_mat)
n_pairs <- choose(n_genes, 2)
bonf <- opt$p_threshold / n_pairs

p_pearson <- matrix(NA_real_, n_genes, n_genes, dimnames = list(keep_genes, keep_genes))
p_spearman <- matrix(NA_real_, n_genes, n_genes, dimnames = list(keep_genes, keep_genes))

sig_pearson <- data.table()
sig_spearman <- data.table()

for (i in seq_len(n_genes - 1)) {
  xi <- agg_mat[, i]
  for (j in (i + 1):n_genes) {
    xj <- agg_mat[, j]

    pearson_p <- tryCatch(
      stats::cor.test(xi, xj, method = "pearson")$p.value,
      error = function(e) NA_real_
    )
    spearman_p <- tryCatch(
      stats::cor.test(xi, xj, method = "spearman", exact = FALSE)$p.value,
      error = function(e) NA_real_
    )

    p_pearson[i, j] <- p_pearson[j, i] <- pearson_p
    p_spearman[i, j] <- p_spearman[j, i] <- spearman_p

    if (!is.na(pearson_p) && pearson_p < bonf) {
      sig_pearson <- rbind(sig_pearson, data.table(
        Gene1 = keep_genes[i],
        Gene2 = keep_genes[j],
        Pvalue = pearson_p
      ))
    }
    if (!is.na(spearman_p) && spearman_p < bonf) {
      sig_spearman <- rbind(sig_spearman, data.table(
        Gene1 = keep_genes[i],
        Gene2 = keep_genes[j],
        Pvalue = spearman_p
      ))
    }
  }
}

fwrite(as.data.table(p_pearson, keep.rownames = "Gene"),
       file.path(result_dir, "pvalues_pearson.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(as.data.table(p_spearman, keep.rownames = "Gene"),
       file.path(result_dir, "pvalues_spearman.tsv.gz"),
       sep = "\t", compress = "gzip")

if (nrow(sig_pearson)) {
  fwrite(sig_pearson, file.path(result_dir, "significant_pairs_pearson.tsv"), sep = "\t")
}
if (nrow(sig_spearman)) {
  fwrite(sig_spearman, file.path(result_dir, "significant_pairs_spearman.tsv"), sep = "\t")
}

summary_dt <- data.table(
  NumDonors = nrow(agg_dt),
  NumGenesInput = length(gene_cols),
  NumGenesTested = n_genes,
  NumPairs = n_pairs,
  NonzeroCutoff = opt$nonzero_cutoff,
  BonferroniThreshold = bonf,
  SignificantPairsPearson = nrow(sig_pearson),
  SignificantPairsSpearman = nrow(sig_spearman)
)
fwrite(summary_dt, file.path(result_dir, "summary.tsv"), sep = "\t")

message("Finished donor-level Pearson/Spearman analysis in: ", opt$out_dir)
