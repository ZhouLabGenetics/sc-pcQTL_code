#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))

script <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1])
root <- normalizePath(file.path(dirname(script), ".."))
source(file.path(root, "config", "config.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

out_dir <- get_arg("--out-dir", file.path(ROOT_DIR, "results/fine_mapping/gwas_all_finemapped"))
chunk_dir <- file.path(out_dir, "official_import_chunks")
status_files <- list.files(chunk_dir, pattern = "^finngen_official_cs_import_status.*[.]tsv$", full.names = TRUE)
cs_files <- list.files(chunk_dir, pattern = "^finngen_official_cs_by_cluster.*[.]tsv$", full.names = TRUE)
if (!length(status_files)) stop("No official all-finemapped import status chunks found in ", chunk_dir)

read_chunk <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(data.table())
  tryCatch(fread(path), error = function(e) data.table())
}

nonempty <- function(dts) {
  dts[vapply(dts, function(x) ncol(x) > 0L, logical(1))]
}

empty_status_dt <- function() {
  data.table(
    phenocode = character(),
    phenotype = character(),
    endpoint_source = character(),
    trait_type = character(),
    category = character(),
    status = character(),
    n_overlapping_cluster_windows = integer(),
    n_imported_cluster_windows = integer(),
    message = character()
  )
}

empty_cs_dt <- function() {
  data.table(
    phenocode = character(),
    phenotype = character(),
    endpoint_source = character(),
    trait_type = character(),
    category = character(),
    celltype = character(),
    cluster_id = character(),
    gwas_region = character(),
    gwas_cs = integer(),
    cs_log10bf = numeric(),
    cs_avg_r2 = numeric(),
    cs_min_r2 = numeric(),
    low_purity = logical(),
    cs_size = integer(),
    gwas_susie_rds = character(),
    gwas_finemap_source = character(),
    ld_source = character()
  )
}

status_list <- nonempty(lapply(status_files, read_chunk))
cs_list <- nonempty(lapply(cs_files, read_chunk))
status_dt <- if (length(status_list)) rbindlist(status_list, fill = TRUE) else empty_status_dt()
cs_dt <- if (length(cs_list)) rbindlist(cs_list, fill = TRUE) else empty_cs_dt()
if (nrow(cs_dt)) cs_dt <- unique(cs_dt)

fwrite(status_dt, file.path(out_dir, "finngen_official_cs_import_status.tsv"), sep = "\t", quote = FALSE)
fwrite(cs_dt, file.path(out_dir, "finngen_official_cs_by_cluster.tsv"), sep = "\t", quote = FALSE)

qc <- rbindlist(list(
  status_dt[, .(N = .N), by = .(value = status)][, metric := "phenotype_import_status"],
  status_dt[, .(N = sum(n_imported_cluster_windows, na.rm = TRUE)), by = .(value = endpoint_source)][, metric := "imported_cluster_windows_by_source"],
  cs_dt[, .(N = .N), by = .(value = endpoint_source)][, metric := "cs_rows_by_source"],
  cs_dt[, .(N = uniqueN(phenocode)), by = .(value = endpoint_source)][, metric := "imported_phenotypes_by_source"]
), fill = TRUE)
setcolorder(qc, c("metric", "value", "N"))
fwrite(qc, file.path(out_dir, "finngen_official_cs_import_qc.tsv"), sep = "\t", quote = FALSE)

cat("Merged", length(status_files), "status chunks and", length(cs_files), "CS chunks\n")
cat("Imported cluster windows:", sum(status_dt$n_imported_cluster_windows, na.rm = TRUE), "\n")
cat("CS rows:", nrow(cs_dt), "\n")
print(qc)
