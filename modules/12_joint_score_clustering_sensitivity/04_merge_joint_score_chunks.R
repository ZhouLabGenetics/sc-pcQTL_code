#!/usr/bin/env Rscript

.local_lib <- Sys.getenv("SC_PCQTL_R_LIBS", unset = "")
if (nzchar(.local_lib)) .libPaths(c(.local_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args_all[grep("--file=", args_all)[1]]))
module_dir <- dirname(script_path)
source(file.path(module_dir, "load_config.R"))
cfg <- load_joint_score_config(module_dir)

manifest <- fread(file.path(cfg$stage_dir, "joint_score_task_manifest.tsv"))
summary_files <- file.path(cfg$chunk_dir, manifest$chunk_id, "chunk_summary.tsv")
complete_files <- file.path(cfg$chunk_dir, manifest$chunk_id, "COMPLETE")
if (any(!file.exists(summary_files)) || any(!file.exists(complete_files))) {
  missing <- manifest$chunk_id[!file.exists(summary_files) | !file.exists(complete_files)]
  stop("Incomplete joint-score chunks: ", paste(head(missing, 30), collapse = ", "))
}

all_summary <- rbindlist(lapply(summary_files, fread), fill = TRUE)
fwrite(all_summary, file.path(cfg$stage_dir, "all_chunk_summaries.tsv"), sep = "\t")

for (chr in 1:22) {
  chr_manifest <- manifest[chromosome == chr]
  if (!nrow(chr_manifest)) next
  sig_files <- file.path(cfg$chunk_dir, chr_manifest$chunk_id, "significant_directions.tsv")
  sig <- rbindlist(lapply(sig_files, function(path) {
    value <- fread(path)
    if (!nrow(value)) return(data.table())
    value
  }), fill = TRUE)

  if (nrow(sig)) {
    setorder(sig, Gene1, Gene2, -neglog10_p_joint)
    pairs <- sig[, .SD[1L], by = .(Gene1, Gene2)]
    direction_counts <- sig[, .(n_significant_directions = .N), by = .(Gene1, Gene2)]
    pairs <- merge(pairs, direction_counts, by = c("Gene1", "Gene2"), sort = FALSE)
    setnames(
      pairs,
      c("p_joint", "stat_joint", "stat_count", "stat_detection", "direction"),
      c("Pvalue_joint", "Statistic_joint", "Statistic_count", "Statistic_detection", "Best_direction")
    )
  } else {
    pairs <- data.table(
      Gene1 = character(), Gene2 = character(), Pvalue_joint = numeric(),
      Statistic_joint = numeric(), Statistic_count = numeric(),
      Statistic_detection = numeric(), Best_direction = character(),
      n_significant_directions = integer()
    )
  }

  out_dir <- file.path(cfg$assoc_dir, sprintf("chr%d", chr))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(pairs, file.path(out_dir, sprintf("chr%d_significant_pairs.tsv", chr)), sep = "\t")
  fwrite(data.table(
    Chromosome = chr,
    NumChunks = nrow(chr_manifest),
    DirectionalTests = all_summary[chromosome == chr, sum(n_directional_tests)],
    CompleteJointTests = all_summary[chromosome == chr, sum(n_complete_joint_tests)],
    SignificantDirections = all_summary[chromosome == chr, sum(n_significant_directions)],
    SignificantPairs = nrow(pairs)
  ), file.path(out_dir, sprintf("chr%d_summary.tsv", chr)), sep = "\t")
}
writeLines("OK", file.path(cfg$stage_dir, "MERGE_COMPLETE"))
message(sprintf("[%s] Joint-score chunks merged", cfg$celltype))
