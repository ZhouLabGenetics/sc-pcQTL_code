#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(script_file), "common.R"))
prepend_project_libraries()
args <- parse_cli()
opt <- list(
  raw_counts = require_arg(args, "raw_counts"),
  global_template = require_arg(args, "global_template"),
  output_dir = require_arg(args, "output_dir"),
  seed = integer_arg(args, "seed", 20240121L)
)

sc_out <- file.path(opt$output_dir, "sc_counts_within_donor_permuted.rds")
if (dir.exists(opt$output_dir) || file.exists(sc_out)) {
  stop("Refusing to overwrite within-donor runtime data: ", opt$output_dir)
}
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

template <- readRDS(opt$global_template)
required_fields <- c("counts", "donors", "cells", "genes")
if (!all(required_fields %in% names(template))) {
  stop("Global template lacks: ", paste(setdiff(required_fields, names(template)), collapse = ", "))
}
genes <- as.character(template$genes)
cells <- as.character(template$cells)
donors <- as.character(template$donors)
if (length(cells) != length(donors)) stop("Template cell and donor lengths differ")

header_con <- gzfile(opt$raw_counts, open = "rt")
header <- strsplit(readLines(header_con, n = 1L), "\t", fixed = TRUE)[[1L]]
close(header_con)
needed <- c("CellID", genes)
if (!all(needed %in% header)) {
  stop("Raw matrix lacks columns: ", paste(setdiff(needed, header), collapse = ", "))
}
col_classes <- rep("NULL", length(header))
col_classes[match("CellID", header)] <- "character"
col_classes[match(genes, header)] <- "numeric"
raw <- read.delim(
  gzfile(opt$raw_counts, open = "rt"), header = TRUE, sep = "\t", quote = "",
  comment.char = "", check.names = FALSE, colClasses = col_classes
)
if (!identical(as.character(raw$CellID), cells)) {
  stop("Raw matrix rows do not match the primary global-permutation template")
}
parsed_donors <- sub("_.*$", "", cells)
if (!identical(parsed_donors, donors)) stop("Parsed donor IDs do not match the template")

original <- as.matrix(raw[, genes, drop = FALSE])
storage.mode(original) <- "double"
rownames(original) <- cells
colnames(original) <- genes
rm(raw)

template_counts <- template$counts[, genes, drop = FALSE]
marginal_match <- vapply(seq_along(genes), function(j) {
  isTRUE(all.equal(
    unname(sort(original[, j])),
    unname(sort(template_counts[, j])),
    tolerance = 0
  ))
}, logical(1L))
if (!all(marginal_match)) {
  stop("Raw and global-template gene marginals differ for: ",
       paste(genes[!marginal_match], collapse = ", "))
}
template_nonzero <- colSums(template_counts > 0)
original_nonzero <- colSums(original > 0)
rm(template_counts)

donor_levels <- unique(donors)
donor_factor <- factor(donors, levels = donor_levels)
donor_indices <- split(seq_along(donors), donor_factor, drop = TRUE)
set.seed(opt$seed)
within <- original
for (j in seq_along(genes)) {
  for (idx in donor_indices) {
    if (length(idx) > 1L) within[idx, j] <- sample(original[idx, j], replace = FALSE)
  }
}

original_sums <- rowsum(original, donor_factor, reorder = FALSE)
within_sums <- rowsum(within, donor_factor, reorder = FALSE)
if (!isTRUE(all.equal(original_sums, within_sums, tolerance = 0))) {
  stop("Within-donor permutation changed a donor-gene sum")
}
donor_n <- as.numeric(table(donor_factor))
original_means <- sweep(original_sums, 1L, donor_n, "/")
within_means <- sweep(within_sums, 1L, donor_n, "/")
if (!isTRUE(all.equal(original_means, within_means, tolerance = 0))) {
  stop("Within-donor permutation changed a donor-gene mean")
}

atomic_save_rds(list(
  counts = within,
  donors = donors,
  cells = cells,
  genes = genes,
  permutation = "within_donor",
  seed = opt$seed,
  global_template = normalizePath(opt$global_template)
), sc_out)
atomic_write_table(data.frame(
  metric = c(
    "n_cells", "n_donors", "n_genes", "seed",
    "max_abs_donor_gene_sum_difference",
    "max_abs_donor_gene_mean_difference",
    "fraction_cell_gene_entries_changed",
    "max_abs_nonzero_count_difference_vs_global_template",
    "genes_with_exact_marginal_match_to_global_template"
  ),
  value = c(
    nrow(within), length(donor_levels), ncol(within), opt$seed,
    max(abs(original_sums - within_sums)),
    max(abs(original_means - within_means)),
    mean(within != original),
    max(abs(original_nonzero - template_nonzero)),
    sum(marginal_match)
  )
), file.path(opt$output_dir, "permutation_validation.tsv"))

message("Within-donor diagnostic permutation written to ", sc_out)
