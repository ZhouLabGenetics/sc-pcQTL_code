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
  make_option("--results_dir",
              default = file.path(get_sim_root(), "01_main_text_results", "01_null_calibration_final_sc", "results", "full_realcounts")),
  make_option("--sc_dir", default = NULL),
  make_option("--pb_dir", default = NULL),
  make_option("--setting_label", default = "default"),
  make_option("--p_thresholds", default = "0.05,0.01,0.005,0.001"),
  make_option("--output_dir",
              default = file.path(get_sim_root(), "01_main_text_results", "01_null_calibration_final_sc", "results", "full_realcounts"))
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
p_thresholds <- as.numeric(strsplit(opt$p_thresholds, ",")[[1]])
p_thresholds <- unique(p_thresholds[!is.na(p_thresholds)])
p_thresholds <- p_thresholds[order(-p_thresholds)]

SC_DIR <- if (is.null(opt$sc_dir)) file.path(opt$results_dir, "sc_hurdle") else opt$sc_dir
PB_DIR <- if (is.null(opt$pb_dir)) file.path(opt$results_dir, "pb_spearman") else opt$pb_dir

message("Loading directional SC hurdle outputs...")
count_file <- file.path(SC_DIR, "sc_hurdle_count_pvalues.tsv")
zero_file <- file.path(SC_DIR, "sc_hurdle_zero_pvalues.tsv")
if (!file.exists(count_file) || !file.exists(zero_file)) {
  stop("Directional SC p-value files are required; pair-mean matrices are not valid final inputs")
}
pair_mask <- load_pair_mask_long(file.path(SC_DIR, "pair_filter.tsv.gz"))
sc_count_long <- fread(count_file)
sc_zero_long <- fread(zero_file)
sc_pairs <- merge(sc_count_long, sc_zero_long,
                  by = c("Gene1", "Gene2"), all = TRUE,
                  suffixes = c("_count", "_zero"))
eligible_undir <- pair_mask[eligible == TRUE, .(
  key1 = pmin(Gene1, Gene2),
  key2 = pmax(Gene1, Gene2)
)]
sc_pairs[, key1 := pmin(Gene1, Gene2)]
sc_pairs[, key2 := pmax(Gene1, Gene2)]
sc_pairs <- merge(sc_pairs, unique(eligible_undir), by = c("key1", "key2"))
sc_pairs[, c("key1", "key2") := NULL]

message("Loading PB outputs...")
pb_long <- load_matrix_long(file.path(PB_DIR, "pvalue_matrix.tsv.gz"))
pb_long <- merge(pb_long, pair_mask[, .(Gene1, Gene2, eligible)], by = c("Gene1", "Gene2"))
pb_long <- pb_long[eligible == TRUE]
pb_long <- rbind(
  pb_long[, .(Gene1, Gene2, Pvalue)],
  pb_long[, .(Gene1 = Gene2, Gene2 = Gene1, Pvalue)],
  use.names = TRUE
)

calc_sc_fpr <- function(dt, thresholds) {
  res <- data.table()
  for (pth in thresholds) {
    flag_count <- dt$Pvalue_count < pth
    flag_zero <- dt$Pvalue_zero < pth
    res <- rbind(res, data.table(
      p_threshold = pth,
      count_fpr = mean(flag_count, na.rm = TRUE),
      zero_fpr = mean(flag_zero, na.rm = TRUE),
      combined_fpr = mean(flag_count | flag_zero, na.rm = TRUE)
    ))
  }
  res
}

calc_pb_fpr <- function(dt, thresholds) {
  res <- data.table()
  for (pth in thresholds) {
    flag <- dt$Pvalue < pth
    res <- rbind(res, data.table(
      p_threshold = pth,
      fpr = mean(flag, na.rm = TRUE)
    ))
  }
  res
}

sc_fpr <- calc_sc_fpr(sc_pairs, p_thresholds)
pb_fpr <- calc_pb_fpr(pb_long, p_thresholds)

fpr_dt <- data.table(
  method = c(rep("SC_count", nrow(sc_fpr)), rep("SC_zero", nrow(sc_fpr)), rep("SC_combined", nrow(sc_fpr)), rep("PB", nrow(pb_fpr))),
  p_threshold = c(sc_fpr$p_threshold, sc_fpr$p_threshold, sc_fpr$p_threshold, pb_fpr$p_threshold),
  fpr = c(sc_fpr$count_fpr, sc_fpr$zero_fpr, sc_fpr$combined_fpr, pb_fpr$fpr)
)
fpr_dt[, setting := opt$setting_label]
fwrite(fpr_dt, file.path(opt$output_dir, "fpr_comparison.tsv"), sep = "\t")

sink(file.path(opt$output_dir, "fpr_comparison_report.txt"))
cat("=== FPR comparison (real-data permutation-based null) ===\n")
cat(sprintf("Setting: %s\n", opt$setting_label))
cat(sprintf("SC dir: %s\nPB dir: %s\n\n", SC_DIR, PB_DIR))
cat(sprintf("SC directional tests: %d\n", nrow(sc_pairs)))
cat(sprintf("PB directional-equivalent records: %d\n\n", nrow(pb_long)))
print(fpr_dt)
sink()

message("FPR comparison finished. Outputs in ", opt$output_dir)
