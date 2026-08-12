#!/usr/bin/env Rscript

.local_libs <- c(
  Sys.getenv("SC_PCQTL_R_LIBS", unset = ""),
  Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = "")
)
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
module_dir <- normalizePath(dirname(args_file), mustWork = TRUE)
release_root <- normalizePath(file.path(module_dir, "..", "..", ".."), mustWork = TRUE)
source(file.path(release_root, "config", "celltype_eligibility.R"))

eqtl_root <- Sys.getenv("ONEK1K_EQTL_ROOT", unset = "")
if (!nzchar(eqtl_root) || !dir.exists(eqtl_root)) {
  stop("Set ONEK1K_EQTL_ROOT to the OneK1K single-gene cis-eQTL tarball directory.")
}
summary_root <- Sys.getenv("SC_PCQTL_PCQTL_SUMMARY_ROOT", unset = "")
if (!nzchar(summary_root)) stop("Set SC_PCQTL_PCQTL_SUMMARY_ROOT.")
overlap_root <- Sys.getenv(
  "SC_PCQTL_PCQTL_EGENE_ROOT",
  unset = file.path(dirname(normalizePath(summary_root, mustWork = FALSE)), "pcqtl_compare_saigeqtl")
)
out_dir <- file.path(overlap_root, "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

eligibility <- load_celltype_eligibility(file.path(release_root, "config", "celltype_eligibility.tsv"))
eligibility <- eligibility[include_primary == TRUE]
fdr_threshold <- 0.05

extract_celltype <- function(celltype, eqtl_celltype) {
  tar_file <- file.path(eqtl_root, sprintf("cis_%s.tar.gz", eqtl_celltype))
  if (!file.exists(tar_file)) stop("Missing OneK1K eQTL tarball: ", tar_file)

  tmp_parent <- Sys.getenv("TMPDIR", unset = tempdir())
  tmp_dir <- tempfile(sprintf("sc_pcqtl_egene_%s_", celltype), tmpdir = tmp_parent)
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  status <- system2("tar", c("-xzf", shQuote(tar_file), "-C", shQuote(tmp_dir)))
  if (!identical(status, 0L)) stop("Failed to extract ", tar_file)
  files <- list.files(tmp_dir, pattern = "[.]singleVar[.]txt$", recursive = TRUE, full.names = TRUE)
  if (!length(files)) stop("No single-gene cis-eQTL files found in ", tar_file)

  suffix <- sprintf("_%s_count_saigeqtl_cis_window_1000000.singleVar.txt", eqtl_celltype)
  rows <- lapply(files, function(path) {
    file_name <- basename(path)
    if (!endsWith(file_name, suffix)) return(NULL)
    gene <- substr(file_name, 1L, nchar(file_name) - nchar(suffix))
    x <- tryCatch(
      fread(path, select = "p.value", showProgress = FALSE),
      error = function(e) NULL
    )
    if (is.null(x) || !nrow(x) || !"p.value" %in% names(x)) return(NULL)
    p <- as.numeric(x$p.value)
    q <- p.adjust(p, method = "BH")
    data.table(
      celltype = celltype,
      gene = gene,
      n_snps = length(p),
      min_p = if (any(is.finite(p))) min(p, na.rm = TRUE) else NA_real_,
      min_fdr = if (any(is.finite(q))) min(q, na.rm = TRUE) else NA_real_,
      is_egene = any(is.finite(q) & q < fdr_threshold),
      fdr_method = "BH_within_gene_cis_snps"
    )
  })
  result <- rbindlist(rows, fill = TRUE)
  if (!nrow(result)) stop("No readable eQTL results in ", tar_file)
  message(
    celltype, ": ", nrow(result), " genes; ",
    result[is_egene == TRUE, .N], " eGenes"
  )
  result
}

all_egenes <- rbindlist(lapply(seq_len(nrow(eligibility)), function(i) {
  extract_celltype(eligibility$celltype[i], eligibility$eqtl_celltype[i])
}), fill = TRUE)
setorder(all_egenes, celltype, gene)

observed <- sort(unique(all_egenes$celltype))
expected <- sort(eligibility$celltype)
if (!identical(observed, expected)) {
  stop("eGene extraction did not return the complete primary cell-type set")
}

summary_ct <- all_egenes[, .(
  n_genes_tested = .N,
  n_egenes = sum(is_egene == TRUE),
  pct_egenes = round(100 * mean(is_egene == TRUE), 2),
  median_p = median(min_p, na.rm = TRUE)
), by = celltype][order(-n_egenes)]

fwrite(all_egenes, file.path(out_dir, "all_egenes_combined.tsv"), sep = "\t")
fwrite(summary_ct, file.path(out_dir, "egene_summary_by_celltype.tsv"), sep = "\t")
message("Wrote primary-analysis eGene summaries to ", out_dir)
