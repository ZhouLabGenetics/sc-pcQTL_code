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
window_kb <- as.integer(get_arg("--window-kb", "250"))

out_dir <- file.path(ROOT_DIR, "results/coloc/qtl_gwas_abf_raw", sprintf("window_%skb", window_kb))
chunk_dir <- file.path(out_dir, "chunks")
files <- list.files(chunk_dir, pattern = sprintf("^qtl_gwas_abf_raw_.*_window_%skb\\.tsv$", window_kb), full.names = TRUE)
if (!length(files)) stop("No ABF raw chunk files found in ", chunk_dir)

dt <- rbindlist(lapply(files, fread), fill = TRUE)
if (!"celltype" %in% names(dt)) stop("ABF chunk output lacks celltype")
n_before_eligibility <- nrow(dt)
dt <- dt[celltype %chin% names(CELLTYPE_EQTL_MAP)]
if (!nrow(dt)) stop("No primary-analysis cell types remain after merging ABF chunks")
if (nrow(dt) < n_before_eligibility) {
  message("Excluded ", n_before_eligibility - nrow(dt), " stale ABF rows below the cell-count threshold")
}
fwrite(dt, file.path(out_dir, "qtl_gwas_abf_raw_coloc_summary.tsv"), sep = "\t", quote = FALSE)

ok <- dt[status == "ok" & n_shared_snps >= 50]
best_pheno <- function(type_label, pheno, type, pp) {
  idx <- which(type == type_label & is.finite(pp))
  if (!length(idx)) return(NA_character_)
  as.character(pheno[idx[which.max(pp[idx])]])
}
max_pp <- function(type_label, type, pp) {
  vals <- pp[type == type_label & is.finite(pp)]
  if (!length(vals)) return(NA_real_)
  as.numeric(max(vals))
}
class_dt <- ok[, .(
  max_eQTL_PPH4 = max_pp("eQTL", qtl_type, PP.H4.abf),
  max_pcQTL_PPH4 = max_pp("pcQTL", qtl_type, PP.H4.abf),
  any_eQTL = any(qtl_type == "eQTL" & coloc_pass, na.rm = TRUE),
  any_pcQTL = any(qtl_type == "pcQTL" & coloc_pass, na.rm = TRUE),
  any_eQTL_h4_0_7 = any(qtl_type == "eQTL" & coloc_pass_h4_0_7, na.rm = TRUE),
  any_pcQTL_h4_0_7 = any(qtl_type == "pcQTL" & coloc_pass_h4_0_7, na.rm = TRUE),
  any_eQTL_h4_0_8 = any(qtl_type == "eQTL" & coloc_pass_h4_0_8, na.rm = TRUE),
  any_pcQTL_h4_0_8 = any(qtl_type == "pcQTL" & coloc_pass_h4_0_8, na.rm = TRUE),
  best_eQTL_phenotype = best_pheno("eQTL", phenotype_id, qtl_type, PP.H4.abf),
  best_pcQTL_phenotype = best_pheno("pcQTL", phenotype_id, qtl_type, PP.H4.abf),
  n_ok_tests = as.integer(.N),
  n_shared_snps_median = as.numeric(median(n_shared_snps, na.rm = TRUE))
), by = .(phenocode, phenotype, category, celltype, cluster_id, region_chr, region_start, region_end, window_kb)]
class_dt[!is.finite(max_eQTL_PPH4), `:=`(max_eQTL_PPH4 = NA_real_, best_eQTL_phenotype = NA_character_)]
class_dt[!is.finite(max_pcQTL_PPH4), `:=`(max_pcQTL_PPH4 = NA_real_, best_pcQTL_phenotype = NA_character_)]
class_dt[, coloc_class := fifelse(any_eQTL & any_pcQTL, "shared",
                           fifelse(any_eQTL, "eQTL_only",
                           fifelse(any_pcQTL, "pcQTL_specific", "none")))]
class_dt[, coloc_class_h4_0_7 := fifelse(any_eQTL_h4_0_7 & any_pcQTL_h4_0_7, "shared",
                                  fifelse(any_eQTL_h4_0_7, "eQTL_only",
                                  fifelse(any_pcQTL_h4_0_7, "pcQTL_specific", "none")))]
class_dt[, coloc_class_h4_0_8 := fifelse(any_eQTL_h4_0_8 & any_pcQTL_h4_0_8, "shared",
                                  fifelse(any_eQTL_h4_0_8, "eQTL_only",
                                  fifelse(any_pcQTL_h4_0_8, "pcQTL_specific", "none")))]
fwrite(class_dt, file.path(out_dir, "qtl_gwas_abf_raw_coloc_classes.tsv"), sep = "\t", quote = FALSE)

top <- ok[order(-PP.H4.abf)][seq_len(min(.N, 500L))]
fwrite(top, file.path(out_dir, "top_abf_raw_coloc_candidates.tsv"), sep = "\t", quote = FALSE)

qc <- rbindlist(list(
  dt[, .N, by = status][, .(metric = "status", value = status, N)],
  dt[, .N, by = .(value = paste(qtl_type, status, sep = ":"))][, metric := "qtl_type_status"][, .(metric, value, N)],
  class_dt[, .N, by = coloc_class][, .(metric = "coloc_class", value = coloc_class, N)],
  class_dt[, .N, by = coloc_class_h4_0_7][, .(metric = "coloc_class_h4_0_7", value = coloc_class_h4_0_7, N)],
  class_dt[, .N, by = coloc_class_h4_0_8][, .(metric = "coloc_class_h4_0_8", value = coloc_class_h4_0_8, N)],
  data.table(metric = "tests", value = "ok_n_shared_ge_50", N = nrow(ok)),
  data.table(metric = "tests", value = "coloc_pass", N = sum(ok$coloc_pass, na.rm = TRUE)),
  data.table(metric = "tests", value = "coloc_pass_h4_0_7", N = sum(ok$coloc_pass_h4_0_7, na.rm = TRUE)),
  data.table(metric = "tests", value = "coloc_pass_h4_0_8", N = sum(ok$coloc_pass_h4_0_8, na.rm = TRUE))
), fill = TRUE)
fwrite(qc, file.path(out_dir, "qtl_gwas_abf_raw_qc.tsv"), sep = "\t", quote = FALSE)

cat("Merged", length(files), "endpoint chunks and", nrow(dt), "ABF raw tests for window", window_kb, "kb\n")
