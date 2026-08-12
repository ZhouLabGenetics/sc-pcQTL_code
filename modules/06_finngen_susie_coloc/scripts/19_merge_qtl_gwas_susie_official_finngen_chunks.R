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

out_dir <- get_arg("--out-dir", file.path(ROOT_DIR, "results/coloc/qtl_gwas_susie_official_finngen"))
chunk_dir <- file.path(out_dir, "chunks")
files <- list.files(chunk_dir, pattern = "^qtl_gwas_susie_official_finngen_coloc_chunk_\\d+\\.tsv$", full.names = TRUE)
if (!length(files)) stop("No SuSiE official FinnGen coloc chunk files found in ", chunk_dir)

dt <- rbindlist(lapply(files, fread), fill = TRUE)
if (!"celltype" %in% names(dt)) stop("SuSiE coloc chunk output lacks celltype")
n_before_eligibility <- nrow(dt)
dt <- dt[celltype %chin% names(CELLTYPE_EQTL_MAP)]
if (!nrow(dt)) stop("No primary-analysis cell types remain after merging SuSiE coloc chunks")
if (nrow(dt) < n_before_eligibility) {
  message("Excluded ", n_before_eligibility - nrow(dt), " stale SuSiE coloc rows below the cell-count threshold")
}
dt[, `:=`(
  pph4 = as.numeric(pph4),
  PP.H3.abf = as.numeric(PP.H3.abf),
  PP.H4.abf = as.numeric(PP.H4.abf),
  n_shared_variants = as.integer(n_shared_variants),
  n_coloc_variants = as.integer(n_coloc_variants),
  qtl_shared_alpha_mass = as.numeric(qtl_shared_alpha_mass),
  gwas_shared_alpha_mass = as.numeric(gwas_shared_alpha_mass),
  shared_mass_pass = as.logical(shared_mass_pass)
)]
dt[, coloc_pass_main := status == "ok" & is.finite(pph4) & pph4 > COLOC_PPH4_CUTOFF]
dt[, coloc_pass_h4_0_7 := status == "ok" & is.finite(pph4) & pph4 > 0.7]
dt[, coloc_pass_h4_0_8 := status == "ok" & is.finite(pph4) & pph4 > COLOC_PPH4_SENSITIVITY_CUTOFF]
dt[, coloc_pass_main_shared_mass_qc := coloc_pass_main & shared_mass_pass == TRUE]
dt[, coloc_pass_h4_0_7_shared_mass_qc := coloc_pass_h4_0_7 & shared_mass_pass == TRUE]
dt[, coloc_pass_h4_0_8_shared_mass_qc := coloc_pass_h4_0_8 & shared_mass_pass == TRUE]
fwrite(dt, file.path(out_dir, "qtl_gwas_susie_official_finngen_coloc_summary.tsv"), sep = "\t", quote = FALSE)

classify_group <- function(any_eqtl, any_pcqtl) {
  fifelse(any_eqtl & any_pcqtl, "shared",
    fifelse(any_eqtl & !any_pcqtl, "eQTL_only",
      fifelse(!any_eqtl & any_pcqtl, "pcQTL_specific", "none")
    )
  )
}

best_pheno <- function(pheno, pph4, type, wanted_type) {
  idx <- which(type == wanted_type & is.finite(pph4))
  if (!length(idx)) return(NA_character_)
  pheno[idx[which.max(pph4[idx])]]
}

class_dt <- dt[, .(
  phenotype = phenotype[1],
  n_tests = .N,
  n_ok_tests = sum(status == "ok", na.rm = TRUE),
  n_shared_variants_median = as.numeric(median(n_shared_variants, na.rm = TRUE)),
  max_eQTL_PPH4 = suppressWarnings(max(pph4[qtl_type == "eQTL"], na.rm = TRUE)),
  max_pcQTL_PPH4 = suppressWarnings(max(pph4[qtl_type == "pcQTL"], na.rm = TRUE)),
  best_eQTL_phenotype = best_pheno(qtl_phenotype_id, pph4, qtl_type, "eQTL"),
  best_pcQTL_phenotype = best_pheno(qtl_phenotype_id, pph4, qtl_type, "pcQTL"),
  any_eQTL_main = any(qtl_type == "eQTL" & coloc_pass_main, na.rm = TRUE),
  any_pcQTL_main = any(qtl_type == "pcQTL" & coloc_pass_main, na.rm = TRUE),
  any_eQTL_h4_0_7 = any(qtl_type == "eQTL" & coloc_pass_h4_0_7, na.rm = TRUE),
  any_pcQTL_h4_0_7 = any(qtl_type == "pcQTL" & coloc_pass_h4_0_7, na.rm = TRUE),
  any_eQTL_h4_0_8 = any(qtl_type == "eQTL" & coloc_pass_h4_0_8, na.rm = TRUE),
  any_pcQTL_h4_0_8 = any(qtl_type == "pcQTL" & coloc_pass_h4_0_8, na.rm = TRUE),
  any_eQTL_main_shared_mass_qc = any(qtl_type == "eQTL" & coloc_pass_main_shared_mass_qc, na.rm = TRUE),
  any_pcQTL_main_shared_mass_qc = any(qtl_type == "pcQTL" & coloc_pass_main_shared_mass_qc, na.rm = TRUE),
  any_eQTL_h4_0_7_shared_mass_qc = any(qtl_type == "eQTL" & coloc_pass_h4_0_7_shared_mass_qc, na.rm = TRUE),
  any_pcQTL_h4_0_7_shared_mass_qc = any(qtl_type == "pcQTL" & coloc_pass_h4_0_7_shared_mass_qc, na.rm = TRUE),
  any_eQTL_h4_0_8_shared_mass_qc = any(qtl_type == "eQTL" & coloc_pass_h4_0_8_shared_mass_qc, na.rm = TRUE),
  any_pcQTL_h4_0_8_shared_mass_qc = any(qtl_type == "pcQTL" & coloc_pass_h4_0_8_shared_mass_qc, na.rm = TRUE)
), by = .(phenocode, celltype, cluster_id)]
class_dt[!is.finite(max_eQTL_PPH4), max_eQTL_PPH4 := NA_real_]
class_dt[!is.finite(max_pcQTL_PPH4), max_pcQTL_PPH4 := NA_real_]
class_dt[, coloc_class_main := classify_group(any_eQTL_main, any_pcQTL_main)]
class_dt[, coloc_class_h4_0_7 := classify_group(any_eQTL_h4_0_7, any_pcQTL_h4_0_7)]
class_dt[, coloc_class_h4_0_8 := classify_group(any_eQTL_h4_0_8, any_pcQTL_h4_0_8)]
class_dt[, coloc_class_main_shared_mass_qc := classify_group(any_eQTL_main_shared_mass_qc, any_pcQTL_main_shared_mass_qc)]
class_dt[, coloc_class_h4_0_7_shared_mass_qc := classify_group(any_eQTL_h4_0_7_shared_mass_qc, any_pcQTL_h4_0_7_shared_mass_qc)]
class_dt[, coloc_class_h4_0_8_shared_mass_qc := classify_group(any_eQTL_h4_0_8_shared_mass_qc, any_pcQTL_h4_0_8_shared_mass_qc)]
fwrite(class_dt, file.path(out_dir, "qtl_gwas_susie_official_finngen_coloc_classes.tsv"), sep = "\t", quote = FALSE)

qc_files <- list.files(chunk_dir, pattern = "^qtl_gwas_susie_official_finngen_qc_chunk_\\d+\\.tsv$", full.names = TRUE)
if (length(qc_files)) {
  liftover_qc <- rbindlist(lapply(qc_files, fread), fill = TRUE)
  if ("celltype" %in% names(liftover_qc)) {
    liftover_qc <- liftover_qc[celltype %chin% names(CELLTYPE_EQTL_MAP)]
  }
  fwrite(liftover_qc, file.path(out_dir, "qtl_gwas_susie_official_finngen_liftover_qc.tsv"), sep = "\t", quote = FALSE)
}

summarize_threshold <- function(class_col, label) {
  cls <- class_dt[[class_col]]
  e_groups <- class_dt[cls %in% c("eQTL_only", "shared")]
  p_groups <- class_dt[cls == "pcQTL_specific"]
  data.table(
    threshold = label,
    n_total_cluster_endpoint_groups = nrow(class_dt),
    n_eQTL_events = nrow(e_groups),
    n_pcQTL_specific_events = nrow(p_groups),
    n_shared_events = sum(cls == "shared", na.rm = TRUE),
    n_eQTL_only_events = sum(cls == "eQTL_only", na.rm = TRUE),
    n_none_groups = sum(cls == "none", na.rm = TRUE),
    n_unique_eQTL_hits = uniqueN(e_groups[, .(phenocode, celltype, cluster_id)]),
    n_unique_pcQTL_specific_hits = uniqueN(p_groups[, .(phenocode, celltype, cluster_id)]),
    percent_event_increase = if (nrow(e_groups) > 0) nrow(p_groups) / nrow(e_groups) else NA_real_,
    percent_unique_increase = if (nrow(e_groups) > 0) {
      uniqueN(p_groups[, .(phenocode, celltype, cluster_id)]) /
        uniqueN(e_groups[, .(phenocode, celltype, cluster_id)])
    } else NA_real_
  )
}

summary_dt <- rbindlist(list(
  summarize_threshold("coloc_class_main", "PPH4_gt_0.75"),
  summarize_threshold("coloc_class_h4_0_7", "PPH4_gt_0.70"),
  summarize_threshold("coloc_class_h4_0_8", "PPH4_gt_0.80"),
  summarize_threshold("coloc_class_main_shared_mass_qc", "PPH4_gt_0.75_shared_alpha_mass_ge_0.9_QC"),
  summarize_threshold("coloc_class_h4_0_7_shared_mass_qc", "PPH4_gt_0.70_shared_alpha_mass_ge_0.9_QC"),
  summarize_threshold("coloc_class_h4_0_8_shared_mass_qc", "PPH4_gt_0.80_shared_alpha_mass_ge_0.9_QC")
))
fwrite(summary_dt, file.path(out_dir, "qtl_gwas_susie_official_finngen_additional_hit_summary.tsv"), sep = "\t", quote = FALSE)

make_qc <- function(value, metric_label) {
  out <- data.table(value = as.character(value))[, .N, by = value]
  out[, metric := metric_label]
  setcolorder(out, c("metric", "value", "N"))
  out
}
qc <- rbindlist(list(
  make_qc(dt$status, "status"),
  make_qc(paste(dt$qtl_type, dt$status, sep = ":"), "qtl_type_status"),
  make_qc(class_dt$coloc_class_main, "coloc_class_main"),
  make_qc(class_dt$coloc_class_h4_0_7, "coloc_class_h4_0_7"),
  make_qc(class_dt$coloc_class_h4_0_8, "coloc_class_h4_0_8"),
  make_qc(class_dt$coloc_class_main_shared_mass_qc, "coloc_class_main_shared_mass_qc"),
  make_qc(class_dt$coloc_class_h4_0_7_shared_mass_qc, "coloc_class_h4_0_7_shared_mass_qc"),
  make_qc(class_dt$coloc_class_h4_0_8_shared_mass_qc, "coloc_class_h4_0_8_shared_mass_qc")
), fill = TRUE)
fwrite(qc, file.path(out_dir, "qtl_gwas_susie_official_finngen_qc.tsv"), sep = "\t", quote = FALSE)

cat("Merged", length(files), "chunk files and", nrow(dt), "SuSiE signal-pair results\n")
print(summary_dt)
