#!/usr/bin/env Rscript

.local_lib <- Sys.getenv("SC_PCQTL_R_LIBS", unset = "")
if (nzchar(.local_lib)) .libPaths(c(.local_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args_all[grep("--file=", args_all)[1]]))
module_dir <- dirname(script_path)
source(file.path(module_dir, "load_config.R"))
cfg <- load_joint_score_config(module_dir)

args <- commandArgs(trailingOnly = TRUE)
chunk_size <- if (length(args)) {
  as.integer(args[1])
} else {
  as.integer(Sys.getenv("SC_PCQTL_JOINT_SCORE_RESPONSE_CHUNK", "50"))
}
if (!is.finite(chunk_size) || chunk_size < 1L) stop("Invalid response chunk size")
if (!file.exists(file.path(cfg$stage_dir, "STAGE_COMPLETE"))) stop("Run 01_stage_inputs.R first")

genes <- fread(file.path(cfg$stage_dir, "filtered_autosomal_genes.tsv"))
rows <- list()
k <- 0L
for (chr in 1:22) {
  n <- genes[chr_numeric == chr, .N]
  if (!n) next
  for (start in seq.int(1L, n, by = chunk_size)) {
    k <- k + 1L
    stop_at <- min(n, start + chunk_size - 1L)
    rows[[k]] <- data.table(
      task_id = k,
      chromosome = chr,
      response_start = start,
      response_end = stop_at,
      n_responses = stop_at - start + 1L,
      chunk_id = sprintf("chr%d_resp%05d_%05d", chr, start, stop_at)
    )
  }
}
manifest <- rbindlist(rows)
fwrite(manifest, file.path(cfg$stage_dir, "joint_score_task_manifest.tsv"), sep = "\t")
message(sprintf("[%s] Wrote %d tasks at up to %d responses per task", cfg$celltype, nrow(manifest), chunk_size))
