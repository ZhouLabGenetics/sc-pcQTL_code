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
  make_option("--output", default = file.path(get_data_root(), "shuffle_full"),
              help = "Base output directory [default %default]"),
  make_option("--label", default = "full_realcounts",
              help = "Subfolder name under output [default %default]"),
  make_option("--n_genes", type = "integer", default = 200),
  make_option("--seed", type = "integer", default = 20240121)
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

base_dir <- file.path(opt$output, opt$label)
if (dir.exists(base_dir) && length(list.files(base_dir, all.files = TRUE, no.. = TRUE))) {
  stop(
    "Refusing to overwrite an existing null dataset: ", base_dir,
    ". Use a new --label or archive/remove the old run explicitly."
  )
}
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

message("=== real-data permutation-based null generator ===")
message(sprintf("Input file: %s", opt$raw_counts))
message(sprintf("Output dir: %s", base_dir))

gene_info <- fread(opt$gene_info, select = "gene_name")
if (nrow(gene_info) < opt$n_genes) stop("n_genes exceeds available genes")
genes <- sample(gene_info$gene_name, opt$n_genes)

counts_dt <- fread(opt$raw_counts, select = c("CellID", genes))
counts_dt[, DonorID := tstrsplit(CellID, "_", keep = 1)]
setcolorder(counts_dt, c("CellID", "DonorID", genes))
cell_meta <- counts_dt[, .(CellID, DonorID)]
sc_mat <- as.matrix(counts_dt[, ..genes])
rownames(sc_mat) <- cell_meta$CellID

permute_columns <- function(mat) {
  out <- mat
  for (j in seq_len(ncol(mat))) out[, j] <- sample(out[, j])
  out
}

message("Permuting each gene independently across cells...")
sc_perm <- permute_columns(sc_mat)

pb_dt <- data.table(DonorID = cell_meta$DonorID, sc_perm)
pb_perm <- pb_dt[, lapply(.SD, mean), by = DonorID, .SDcols = genes]
pb_mat <- as.matrix(pb_perm[, ..genes])
rownames(pb_mat) <- pb_perm$DonorID

saveRDS(list(
  counts = sc_perm,
  donors = cell_meta$DonorID,
  cells = cell_meta$CellID,
  genes = genes
), file.path(base_dir, "sc_counts_shuffle_null.rds"))

saveRDS(list(
  counts = pb_mat,
  donors = pb_perm$DonorID,
  genes = genes
), file.path(base_dir, "pb_counts_shuffle_null.rds"))

fwrite(cell_meta, file.path(base_dir, "cell_metadata.tsv"), sep = "\t")
fwrite(data.table(gene = genes), file.path(base_dir, "selected_genes.tsv"), sep = "\t")

message(sprintf("Done. Donors: %d  cells: %d  genes: %d",
                length(unique(cell_meta$DonorID)), nrow(cell_meta), length(genes)))
message("Artifacts written to ", base_dir)
