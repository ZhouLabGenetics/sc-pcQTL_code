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
  make_option("--truth_file", help = "truth_pairs.tsv"),
  make_option("--sc_null_dir", help = "SC hurdle outputs for null dataset"),
  make_option("--sc_signal_dir", help = "SC hurdle outputs for signal dataset"),
  make_option("--pb_null_dir", help = "PB outputs for null dataset"),
  make_option("--pb_signal_dir", help = "PB outputs for signal dataset"),
  make_option("--p_thresholds", default = "0.05,0.01,0.005,0.001"),
  make_option("--output_dir", default = "results/power_full/compare"),
  make_option("--label", default = "signal")
)

opt <- parse_args(OptionParser(option_list = option_list))
required <- c("truth_file", "sc_null_dir", "sc_signal_dir",
              "pb_null_dir", "pb_signal_dir")
for (nm in required) {
  val <- opt[[nm]]
  if (is.null(val) || !file.exists(val)) {
    stop(sprintf("Missing or invalid --%s", nm))
  }
}
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
p_thresholds <- as.numeric(strsplit(opt$p_thresholds, ",")[[1]])
p_thresholds <- unique(p_thresholds[!is.na(p_thresholds)])
p_thresholds <- p_thresholds[order(-p_thresholds)]

truth <- fread(opt$truth_file)
setnames(truth, names(truth), c("Gene1", "Gene2", "Label"))
truth_dir <- rbind(
  truth,
  truth[, .(Gene1 = Gene2, Gene2 = Gene1, Label)],
  use.names = TRUE
)

pb_null <- load_matrix_long(file.path(opt$pb_null_dir, "pvalue_matrix.tsv.gz"))
pb_signal <- load_matrix_long(file.path(opt$pb_signal_dir, "pvalue_matrix.tsv.gz"))

load_sc_directional <- function(sc_dir) {
  count_file <- file.path(sc_dir, "sc_hurdle_count_pvalues.tsv")
  zero_file <- file.path(sc_dir, "sc_hurdle_zero_pvalues.tsv")
  mask <- load_pair_mask_long(file.path(sc_dir, "pair_filter.tsv.gz"))
  if (!file.exists(count_file) || !file.exists(zero_file)) {
    stop("Directional SC p-value files are required; pair-mean matrices are not valid final inputs: ", sc_dir)
  }
  dt_count <- fread(count_file)
  dt_zero <- fread(zero_file)
  dt <- merge(dt_count, dt_zero,
              by = c("Gene1", "Gene2"),
              all = TRUE, suffixes = c("_count", "_zero"))
  eligible_undir <- mask[eligible == TRUE, .(
    key1 = pmin(Gene1, Gene2),
    key2 = pmax(Gene1, Gene2)
  )]
  dt[, key1 := pmin(Gene1, Gene2)]
  dt[, key2 := pmax(Gene1, Gene2)]
  dt <- merge(dt, unique(eligible_undir), by = c("key1", "key2"))
  dt[, c("key1", "key2") := NULL]
  list(dt = dt, mask = mask)
}

sc_null_obj <- load_sc_directional(opt$sc_null_dir)
sc_signal_obj <- load_sc_directional(opt$sc_signal_dir)
sc_null <- sc_null_obj$dt
sc_signal <- sc_signal_obj$dt
mask_null <- sc_null_obj$mask
mask_signal <- sc_signal_obj$mask

annotate_truth <- function(dt) merge(dt, truth_dir, by = c("Gene1", "Gene2"), all.x = TRUE)
sc_null <- annotate_truth(sc_null)
sc_signal <- annotate_truth(sc_signal)
pb_null <- merge(pb_null, mask_null[, .(Gene1, Gene2, eligible)], by = c("Gene1", "Gene2"))[eligible == TRUE]
pb_signal <- merge(pb_signal, mask_signal[, .(Gene1, Gene2, eligible)], by = c("Gene1", "Gene2"))[eligible == TRUE]
pb_null <- rbind(pb_null[, .(Gene1, Gene2, Pvalue)],
                 pb_null[, .(Gene1 = Gene2, Gene2 = Gene1, Pvalue)],
                 use.names = TRUE)
pb_signal <- rbind(pb_signal[, .(Gene1, Gene2, Pvalue)],
                   pb_signal[, .(Gene1 = Gene2, Gene2 = Gene1, Pvalue)],
                   use.names = TRUE)
pb_null <- annotate_truth(pb_null)
pb_signal <- annotate_truth(pb_signal)

calc_metrics <- function(dt, thresholds) {
  res <- data.table()
  for (pth in thresholds) {
    dt[, sig_count := (Pvalue_count < pth)]
    dt[, sig_zero := (Pvalue_zero < pth)]
    dt[, sig_combined := (sig_count | sig_zero)]
    subset_null <- dt[Label == "null"]
    subset_signal <- dt[Label == "signal"]
    res <- rbind(res, data.table(
      p_threshold = pth,
      count_fpr = mean(subset_null$sig_count, na.rm = TRUE),
      zero_fpr = mean(subset_null$sig_zero, na.rm = TRUE),
      combined_fpr = mean(subset_null$sig_combined, na.rm = TRUE),
      count_power = mean(subset_signal$sig_count, na.rm = TRUE),
      zero_power = mean(subset_signal$sig_zero, na.rm = TRUE),
      combined_power = mean(subset_signal$sig_combined, na.rm = TRUE)
    ))
  }
  res
}

calc_pb <- function(dt, thresholds) {
  res <- data.table()
  for (pth in thresholds) {
    dt[, sig := (Pvalue < pth)]
    subset_null <- dt[Label == "null"]
    subset_signal <- dt[Label == "signal"]
    res <- rbind(res, data.table(
      p_threshold = pth,
      fpr = mean(subset_null$sig, na.rm = TRUE),
      power = mean(subset_signal$sig, na.rm = TRUE)
    ))
  }
  res
}

sc_null_metrics <- calc_metrics(sc_null, p_thresholds)
sc_signal_metrics <- calc_metrics(sc_signal, p_thresholds)
pb_null_metrics <- calc_pb(pb_null, p_thresholds)
pb_signal_metrics <- calc_pb(pb_signal, p_thresholds)

fpr_dt <- data.table(
  method = c(rep("SC_count", nrow(sc_null_metrics)),
             rep("SC_zero", nrow(sc_null_metrics)),
             rep("SC_combined", nrow(sc_null_metrics)),
             rep("PB", nrow(pb_null_metrics))),
  p_threshold = c(sc_null_metrics$p_threshold,
                  sc_null_metrics$p_threshold,
                  sc_null_metrics$p_threshold,
                  pb_null_metrics$p_threshold),
  fpr = c(sc_null_metrics$count_fpr,
          sc_null_metrics$zero_fpr,
          sc_null_metrics$combined_fpr,
          pb_null_metrics$fpr)
)

power_dt <- data.table(
  method = c(rep("SC_count", nrow(sc_signal_metrics)),
             rep("SC_zero", nrow(sc_signal_metrics)),
             rep("SC_combined", nrow(sc_signal_metrics)),
             rep("PB", nrow(pb_signal_metrics))),
  p_threshold = c(sc_signal_metrics$p_threshold,
                  sc_signal_metrics$p_threshold,
                  sc_signal_metrics$p_threshold,
                  pb_signal_metrics$p_threshold),
  power = c(sc_signal_metrics$count_power,
            sc_signal_metrics$zero_power,
            sc_signal_metrics$combined_power,
            pb_signal_metrics$power)
)

fwrite(fpr_dt, file.path(opt$output_dir, "fpr_comparison.tsv"), sep = "\t")
fwrite(power_dt, file.path(opt$output_dir, "power_comparison.tsv"), sep = "\t")

sink(file.path(opt$output_dir, "power_summary.txt"))
cat("=== Power comparison summary ===\n")
cat(sprintf("Label: %s\n", opt$label))
cat(sprintf("Eligible null directional tests shared by SC/PB: %d\n", nrow(sc_null)))
cat(sprintf("Eligible signal directional tests shared by SC/PB: %d\n\n", nrow(sc_signal)))
cat("FPR (null dataset):\n")
print(fpr_dt)
cat("\nPower (signal dataset):\n")
print(power_dt)
sink()

message("Results saved to ", opt$output_dir)
