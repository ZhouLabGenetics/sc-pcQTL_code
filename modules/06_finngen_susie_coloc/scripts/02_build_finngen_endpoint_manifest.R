#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))

script <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1])
root <- normalizePath(file.path(dirname(script), ".."))
source(file.path(root, "config", "config.R"))

dir_create(file.path(ROOT_DIR, "manifests"))
dir_create(file.path(FINNGEN_ROOT, "metadata"))

manifest_file <- file.path(FINNGEN_ROOT, "metadata", "finngen_R12_manifest.tsv")
if (!file.exists(manifest_file)) {
  status <- system2(
    "curl",
    c("-L", "--fail", "--show-error", "--output", shQuote(manifest_file), shQuote(FINNGEN_R12_MANIFEST_URL))
  )
  if (status != 0L) stop("Failed to download FinnGen manifest: ", FINNGEN_R12_MANIFEST_URL)
}

fg <- fread(manifest_file)
setnames(fg, gsub("^#", "", names(fg)))
fg[, text_for_filter := tolower(paste(phenocode, phenotype, category, sep = " | "))]
fg[, priority_reason := fifelse(grepl(FINNGEN_PRIORITY_REGEX, text_for_filter, perl = TRUE),
                                "immune_inflammation_blood_infection_keyword_or_category", "out_of_scope")]
fg[, analysis_tier := fifelse(priority_reason != "out_of_scope" & num_cases >= MAIN_MIN_CASES, "main",
                       fifelse(priority_reason != "out_of_scope", "excluded_low_cases", "excluded_out_of_scope"))]
fg[, local_gz := file.path(FINNGEN_ROOT, basename(path_https))]
fg[, local_tbi := paste0(local_gz, ".tbi")]
fg[, downloaded := file.exists(local_gz)]
fg[, downloaded_tbi := file.exists(local_tbi)]

fwrite(fg, file.path(ROOT_DIR, "manifests/finngen_all_endpoints_with_selection.tsv"), sep = "\t", quote = FALSE)
fwrite(fg[analysis_tier == "main"], file.path(ROOT_DIR, "manifests/finngen_main_endpoints.tsv"), sep = "\t", quote = FALSE)

qc <- fg[, .N, by = .(analysis_tier)][order(analysis_tier)]
fwrite(qc, file.path(ROOT_DIR, "manifests/finngen_endpoint_selection_qc.tsv"), sep = "\t", quote = FALSE)

message("FinnGen endpoint manifest written.")
print(qc)
