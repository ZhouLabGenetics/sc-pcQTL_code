#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
module06_code_dir <- normalizePath(file.path(dirname(args_file), "..", ".."), mustWork = TRUE)
source(file.path(module06_code_dir, "config", "config.R"))
root_dir <- ROOT_DIR

input_file <- file.path(
  root_dir,
  "results/coloc/qtl_gwas_susie_official_finngen_all_finemapped/qtl_gwas_susie_official_finngen_coloc_summary.tsv"
)
out_dir <- file.path(root_dir, "05_post_coloc_susie_official_finngen_all_finemapped/results/signal_groups")
dir_create(out_dir)
if (!file.exists(input_file)) stop("Missing SuSiE coloc summary: ", input_file)

message("Reading SuSiE official FinnGen summary: ", input_file)
dt <- fread(input_file, showProgress = TRUE)
if (!"celltype" %in% names(dt)) stop("SuSiE coloc summary lacks celltype")
dt <- dt[celltype %chin% names(CELLTYPE_EQTL_MAP)]
if (!nrow(dt)) stop("No primary-analysis cell types remain in SuSiE coloc summary")
dt[, `:=`(
  pph4 = as.numeric(pph4),
  PP.H3.abf = as.numeric(PP.H3.abf),
  PP.H4.abf = as.numeric(PP.H4.abf),
  n_shared_variants = as.integer(n_shared_variants),
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

dt[, susie_signal_group_id := paste("SUSIESG", celltype, cluster_id, phenocode, sep = "__")]
fwrite(dt, file.path(out_dir, "susie_coloc_event_table.tsv"), sep = "\t")

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

groups <- dt[, .(
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
), by = .(susie_signal_group_id, phenocode, celltype, cluster_id)]
groups[!is.finite(max_eQTL_PPH4), max_eQTL_PPH4 := NA_real_]
groups[!is.finite(max_pcQTL_PPH4), max_pcQTL_PPH4 := NA_real_]
groups[, coloc_class_main := classify_group(any_eQTL_main, any_pcQTL_main)]
groups[, coloc_class_h4_0_7 := classify_group(any_eQTL_h4_0_7, any_pcQTL_h4_0_7)]
groups[, coloc_class_h4_0_8 := classify_group(any_eQTL_h4_0_8, any_pcQTL_h4_0_8)]
groups[, coloc_class_main_shared_mass_qc := classify_group(any_eQTL_main_shared_mass_qc, any_pcQTL_main_shared_mass_qc)]
groups[, coloc_class_h4_0_7_shared_mass_qc := classify_group(any_eQTL_h4_0_7_shared_mass_qc, any_pcQTL_h4_0_7_shared_mass_qc)]
groups[, coloc_class_h4_0_8_shared_mass_qc := classify_group(any_eQTL_h4_0_8_shared_mass_qc, any_pcQTL_h4_0_8_shared_mass_qc)]
fwrite(groups, file.path(out_dir, "susie_signal_groups.tsv"), sep = "\t")

summarize_threshold <- function(class_col, label) {
  cls <- groups[[class_col]]
  e_groups <- groups[cls %in% c("eQTL_only", "shared")]
  p_groups <- groups[cls == "pcQTL_specific"]
  data.table(
    threshold = label,
    n_total_groups = nrow(groups),
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
fwrite(summary_dt, file.path(out_dir, "susie_additional_hit_summary.tsv"), sep = "\t")

by_celltype <- melt(
  groups[, .(celltype, coloc_class_main, coloc_class_h4_0_7, coloc_class_h4_0_8)],
  id.vars = "celltype",
  variable.name = "threshold",
  value.name = "coloc_class"
)[, .N, by = .(threshold, celltype, coloc_class)]
fwrite(by_celltype, file.path(out_dir, "susie_coloc_class_by_celltype.tsv"), sep = "\t")

message("Wrote grouped SuSiE outputs to ", out_dir)
print(summary_dt)
