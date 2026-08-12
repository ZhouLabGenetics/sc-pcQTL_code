#!/usr/bin/env Rscript

.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- normalizePath(dirname(args_file), mustWork = TRUE)
release_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
source(file.path(release_root, "config", "celltype_eligibility.R"))
eligibility <- load_celltype_eligibility(file.path(release_root, "config", "celltype_eligibility.tsv"))
upstream <- Sys.getenv("SC_PCQTL_UPSTREAM_ROOT", unset = "")
if (!nzchar(upstream)) stop("Set SC_PCQTL_UPSTREAM_ROOT to the final add-covariate upstream workflow.")
celltype_root <- file.path(normalizePath(upstream, mustWork = TRUE), "celltypes")
manuscript_root <- Sys.getenv("SC_PCQTL_MANUSCRIPT_ROOT", unset = "")
out_default <- if (nzchar(manuscript_root)) file.path(manuscript_root, "data", "cluster_landscape") else ""
out_configured <- Sys.getenv("SC_PCQTL_CLUSTER_LANDSCAPE_DATA", unset = out_default)
if (!nzchar(out_configured)) {
  stop("Set SC_PCQTL_CLUSTER_LANDSCAPE_DATA or SC_PCQTL_MANUSCRIPT_ROOT.")
}
out_dir <- normalizePath(
  out_configured,
  mustWork = FALSE
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

celltypes <- eligibility[include_primary == TRUE, celltype]
labels <- canonical_celltype_map(file.path(release_root, "config", "celltype_eligibility.tsv"))

cluster_rows <- rbindlist(lapply(celltypes, function(celltype) {
  path <- file.path(celltype_root, celltype, "pcQTL", "step2_pca", "all_clusters_summary.tsv")
  if (!file.exists(path)) stop("Missing cluster PCA summary: ", path)
  x <- fread(path)
  required <- c("status", "n_genes", "n_samples")
  if (length(setdiff(required, names(x)))) stop("Missing required columns in ", path)
  x <- x[status == "SUCCESS"]
  if (!nrow(x)) return(NULL)
  x[, `:=`(cell_type = celltype, cell_type_label = unname(labels[celltype]))]
  x
}), fill = TRUE)

summary_dt <- cluster_rows[, .(
  cells = as.integer(n_samples[which.max(is.finite(n_samples))]),
  clusters = .N,
  assigned_genes = as.integer(sum(as.integer(n_genes))),
  median_size = as.numeric(median(as.integer(n_genes))),
  max_size = as.integer(max(as.integer(n_genes))),
  mean_size = as.numeric(mean(as.integer(n_genes)))
), by = .(cell_type, cell_type_label)]
setorder(summary_dt, -cells)

cluster_rows[, size_group := fifelse(as.integer(n_genes) >= 6L, "6+", as.character(as.integer(n_genes)))]
size_dt <- cluster_rows[, .(n_clusters = .N), by = size_group]
size_dt[, size_order := suppressWarnings(as.integer(size_group))]
size_dt[is.na(size_order) & size_group == "6+", size_order := 6L]
setorder(size_dt, size_order)
size_dt[, size_order := NULL]

fwrite(summary_dt, file.path(out_dir, "sc_cluster_summary_by_celltype.tsv"), sep = "\t")
fwrite(size_dt, file.path(out_dir, "sc_cluster_size_distribution.tsv"), sep = "\t")
message("Wrote cluster-landscape source tables to ", out_dir)
