#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

script <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1])
root <- normalizePath(file.path(dirname(script), ".."))
source(file.path(root, "config", "config.R"))

dir_create <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

download_gcs_object <- function(object_name, out_file) {
  if (file.exists(out_file) && file.info(out_file)$size > 0) return(out_file)
  dir_create(dirname(out_file))
  url <- paste0(
    "https://storage.googleapis.com/download/storage/v1/b/finngen-public-data-r12/o/",
    URLencode(object_name, reserved = TRUE),
    "?alt=media"
  )
  tmp <- paste0(out_file, ".tmp.", Sys.getpid())
  status <- system2("curl", c("-k", "-L", "--fail", "--retry", "3", "-o", tmp, url))
  if (status != 0L || !file.exists(tmp) || file.info(tmp)$size == 0) {
    unlink(tmp)
    stop("Failed to download GCS object: ", object_name)
  }
  file.rename(tmp, out_file)
  out_file
}

list_gcs_objects <- function(prefix) {
  base <- "https://storage.googleapis.com/storage/v1/b/finngen-public-data-r12/o"
  token <- ""
  rows <- character()
  repeat {
    url <- paste0(base, "?prefix=", URLencode(prefix, reserved = TRUE), "&maxResults=1000")
    if (nzchar(token)) url <- paste0(url, "&pageToken=", URLencode(token, reserved = TRUE))
    x <- fromJSON(url)
    if (!is.null(x$items) && nrow(x$items)) rows <- c(rows, x$items$name)
    if (is.null(x$nextPageToken) || !nzchar(x$nextPageToken)) break
    token <- x$nextPageToken
  }
  data.table(name = rows)
}

parse_finngen_id <- function(path, suffix_re) {
  sub("^finngen_R12_", "", sub(suffix_re, "", basename(path)))
}

fg_root <- FINNGEN_ROOT
lab_root <- file.path(fg_root, "lab_values")
doc_dir <- file.path(lab_root, "analysis_documentation")
dir_create(doc_dir)

core_manifest_file <- file.path(fg_root, "finngen_R12_manifest.tsv")
download_gcs_object("summary_stats/finngen_R12_manifest.tsv", core_manifest_file)

lab_summary_doc <- file.path(doc_dir, "Kanta_labs_GWAS_results_v2_summary.txt")
download_gcs_object("lab_values/analysis_documentation/Kanta_labs_GWAS_results_v2_summary.txt", lab_summary_doc)
download_gcs_object(
  "lab_values/analysis_documentation/Kanta_labs_GWAS_results_v2_summary.xlsx",
  file.path(doc_dir, "Kanta_labs_GWAS_results_v2_summary.xlsx")
)
download_gcs_object(
  "lab_values/analysis_documentation/FinnGen_kanta_lab_documentation_v2.pdf",
  file.path(doc_dir, "FinnGen_kanta_lab_documentation_v2.pdf")
)

core <- fread(core_manifest_file)
setnames(core, gsub("^#", "", names(core)))

core_susie_objects <- list_gcs_objects("finemap/full/susie/")
core_susie_objects[, file := basename(name)]
core_susie_objects[, kind := fifelse(
  grepl("[.]SUSIE[.]cred[.]bgz$", file), "cred",
  fifelse(grepl("[.]SUSIE[.]snp[.]bgz$", file), "snp",
          fifelse(grepl("[.]SUSIE[.]snp[.]bgz[.]tbi$", file), "snp_tbi",
                  fifelse(grepl("[.]SUSIE[.]cred_99[.]bgz$", file), "cred_99", "other")))
)]
core_susie_objects[, phenocode := sub("^finngen_R12_", "", sub("[.]SUSIE[.].*$", "", file))]
fwrite(core_susie_objects, file.path(fg_root, "finngen_R12_official_susie_bucket_objects.tsv"), sep = "\t", quote = FALSE)

core_susie <- dcast(
  unique(core_susie_objects[kind %in% c("cred", "cred_99", "snp", "snp_tbi"), .(phenocode, kind)]),
  phenocode ~ kind,
  fun.aggregate = length
)
for (cc in c("cred", "cred_99", "snp", "snp_tbi")) if (!cc %in% names(core_susie)) core_susie[, (cc) := 0L]

core_out <- core[, .(
  phenocode = as.character(phenocode),
  phenotype,
  category,
  endpoint_source = "core_disease",
  trait_type = "cc",
  analysis_type = "case_control",
  unit = NA_character_,
  num_cases = as.integer(num_cases),
  num_controls = as.integer(num_controls),
  n_total = as.integer(num_cases + num_controls),
  path_bucket,
  path_https
)]
core_out <- merge(core_out, core_susie, by = "phenocode", all.x = TRUE)
for (cc in c("cred", "cred_99", "snp", "snp_tbi")) core_out[is.na(get(cc)), (cc) := 0L]
core_out[, `:=`(
  official_susie_in_bucket = cred > 0 & snp > 0 & snp_tbi > 0,
  official_susie_has_cred = cred > 0,
  official_susie_has_cred_99 = cred_99 > 0,
  official_susie_has_snp = snp > 0,
  official_susie_has_tbi = snp_tbi > 0,
  susie_cred_url = paste0("https://storage.googleapis.com/download/storage/v1/b/finngen-public-data-r12/o/finemap%2Ffull%2Fsusie%2Ffinngen_R12_", phenocode, ".SUSIE.cred.bgz?alt=media"),
  susie_snp_url = paste0("https://storage.googleapis.com/download/storage/v1/b/finngen-public-data-r12/o/finemap%2Ffull%2Fsusie%2Ffinngen_R12_", phenocode, ".SUSIE.snp.bgz?alt=media"),
  susie_snp_tbi_url = paste0("https://storage.googleapis.com/download/storage/v1/b/finngen-public-data-r12/o/finemap%2Ffull%2Fsusie%2Ffinngen_R12_", phenocode, ".SUSIE.snp.bgz.tbi?alt=media"),
  local_gz = file.path(fg_root, paste0("finngen_R12_", phenocode, ".gz")),
  local_tbi = file.path(fg_root, paste0("finngen_R12_", phenocode, ".gz.tbi")),
  local_susie_cred = file.path(fg_root, "finemapping/r12_susie", paste0("finngen_R12_", phenocode, ".SUSIE.cred.bgz")),
  local_susie_snp = file.path(fg_root, "finemapping/r12_susie", paste0("finngen_R12_", phenocode, ".SUSIE.snp.bgz")),
  local_susie_tbi = file.path(fg_root, "finemapping/r12_susie", paste0("finngen_R12_", phenocode, ".SUSIE.snp.bgz.tbi"))
)]
core_out[, `:=`(
  local_gz_exists = file.exists(local_gz),
  local_tbi_exists = file.exists(local_tbi),
  local_susie_complete = file.exists(local_susie_cred) & file.exists(local_susie_snp) & file.exists(local_susie_tbi)
)]

lab_doc <- fread(lab_summary_doc)
setnames(lab_doc, "OMOPID", "phenocode")
lab_doc[, phenocode := as.character(phenocode)]

lab_summary_objects <- list_gcs_objects("lab_values/summary_stats/")
lab_susie_objects <- list_gcs_objects("lab_values/finemap/full/susie/")
fwrite(lab_summary_objects, file.path(lab_root, "lab_values_summary_stats_bucket_objects.tsv"), sep = "\t", quote = FALSE)
fwrite(lab_susie_objects, file.path(lab_root, "lab_values_official_susie_bucket_objects.tsv"), sep = "\t", quote = FALSE)

lab_summary <- lab_summary_objects[grepl("finngen_R12_[0-9]+[.]gz$", name)]
lab_summary[, phenocode := parse_finngen_id(name, "[.]gz$")]
lab_summary_tbi <- lab_summary_objects[grepl("finngen_R12_[0-9]+[.]gz[.]tbi$", name)]
lab_summary_tbi[, phenocode := parse_finngen_id(name, "[.]gz[.]tbi$")]

lab_susie <- copy(lab_susie_objects)
lab_susie[, file := basename(name)]
lab_susie[, kind := fifelse(
  grepl("[.]SUSIE[.]cred[.]bgz$", file), "cred",
  fifelse(grepl("[.]SUSIE[.]snp[.]bgz$", file), "snp",
          fifelse(grepl("[.]SUSIE[.]snp[.]bgz[.]tbi$", file), "snp_tbi",
                  fifelse(grepl("[.]SUSIE[.]cred_99[.]bgz$", file), "cred_99", "other")))
)]
lab_susie[, phenocode := sub("^finngen_R12_", "", sub("[.]SUSIE[.].*$", "", file))]
lab_susie_wide <- dcast(
  unique(lab_susie[kind %in% c("cred", "cred_99", "snp", "snp_tbi"), .(phenocode, kind)]),
  phenocode ~ kind,
  fun.aggregate = length
)
for (cc in c("cred", "cred_99", "snp", "snp_tbi")) if (!cc %in% names(lab_susie_wide)) lab_susie_wide[, (cc) := 0L]

lab_out <- lab_doc[, .(
  phenocode,
  phenotype = phenostring,
  category = "Kanta laboratory values",
  endpoint_source = "lab_values",
  trait_type = fifelse(grepl("quantitative", AnalysisType, ignore.case = TRUE), "quant", "binary_or_other"),
  analysis_type = AnalysisType,
  unit,
  num_cases = as.integer(N_cases),
  num_controls = as.integer(N_controls),
  n_total = as.integer(N_total),
  path_bucket = paste0("gs://finngen-public-data-r12/lab_values/summary_stats/finngen_R12_", phenocode, ".gz"),
  path_https = paste0("https://storage.googleapis.com/download/storage/v1/b/finngen-public-data-r12/o/lab_values%2Fsummary_stats%2Ffinngen_R12_", phenocode, ".gz?alt=media")
)]
lab_out <- merge(lab_out, unique(lab_summary[, .(phenocode, has_summary_gz = TRUE)]), by = "phenocode", all.x = TRUE)
lab_out <- merge(lab_out, unique(lab_summary_tbi[, .(phenocode, has_summary_tbi = TRUE)]), by = "phenocode", all.x = TRUE)
lab_out <- merge(lab_out, lab_susie_wide, by = "phenocode", all.x = TRUE)
for (cc in c("has_summary_gz", "has_summary_tbi")) lab_out[is.na(get(cc)), (cc) := FALSE]
for (cc in c("cred", "cred_99", "snp", "snp_tbi")) lab_out[is.na(get(cc)), (cc) := 0L]
lab_out[, `:=`(
  official_susie_in_bucket = cred > 0 & snp > 0 & snp_tbi > 0,
  official_susie_has_cred = cred > 0,
  official_susie_has_cred_99 = cred_99 > 0,
  official_susie_has_snp = snp > 0,
  official_susie_has_tbi = snp_tbi > 0,
  susie_cred_url = paste0("https://storage.googleapis.com/download/storage/v1/b/finngen-public-data-r12/o/lab_values%2Ffinemap%2Ffull%2Fsusie%2Ffinngen_R12_", phenocode, ".SUSIE.cred.bgz?alt=media"),
  susie_snp_url = paste0("https://storage.googleapis.com/download/storage/v1/b/finngen-public-data-r12/o/lab_values%2Ffinemap%2Ffull%2Fsusie%2Ffinngen_R12_", phenocode, ".SUSIE.snp.bgz?alt=media"),
  susie_snp_tbi_url = paste0("https://storage.googleapis.com/download/storage/v1/b/finngen-public-data-r12/o/lab_values%2Ffinemap%2Ffull%2Fsusie%2Ffinngen_R12_", phenocode, ".SUSIE.snp.bgz.tbi?alt=media"),
  local_gz = file.path(fg_root, "lab_values/summary_stats", paste0("finngen_R12_", phenocode, ".gz")),
  local_tbi = file.path(fg_root, "lab_values/summary_stats", paste0("finngen_R12_", phenocode, ".gz.tbi")),
  local_susie_cred = file.path(fg_root, "lab_values/finemapping/r12_susie", paste0("finngen_R12_", phenocode, ".SUSIE.cred.bgz")),
  local_susie_snp = file.path(fg_root, "lab_values/finemapping/r12_susie", paste0("finngen_R12_", phenocode, ".SUSIE.snp.bgz")),
  local_susie_tbi = file.path(fg_root, "lab_values/finemapping/r12_susie", paste0("finngen_R12_", phenocode, ".SUSIE.snp.bgz.tbi"))
)]
lab_out[, `:=`(
  local_gz_exists = file.exists(local_gz),
  local_tbi_exists = file.exists(local_tbi),
  local_susie_complete = file.exists(local_susie_cred) & file.exists(local_susie_snp) & file.exists(local_susie_tbi)
)]

common_cols <- intersect(names(core_out), names(lab_out))
all_sources <- rbindlist(list(core_out[, ..common_cols], lab_out[, ..common_cols]), fill = TRUE)
setorder(all_sources, endpoint_source, phenocode)

dir_create(file.path(ROOT_DIR, "manifests"))
fwrite(core_out, file.path(ROOT_DIR, "manifests/finngen_core_disease_phenotype_standard.tsv"), sep = "\t", quote = FALSE)
fwrite(lab_out, file.path(ROOT_DIR, "manifests/finngen_lab_values_phenotype_standard.tsv"), sep = "\t", quote = FALSE)
fwrite(all_sources, file.path(ROOT_DIR, "manifests/finngen_coloc_gwas_phenotype_standard_all_sources.tsv"), sep = "\t", quote = FALSE)

qc <- rbindlist(list(
  data.table(metric = "core_manifest_rows", value = "n", N = nrow(core_out)),
  data.table(metric = "core_official_susie_complete_in_bucket", value = "n", N = core_out[official_susie_in_bucket == TRUE, .N]),
  data.table(metric = "lab_manifest_rows", value = "n", N = nrow(lab_out)),
  data.table(metric = "lab_summary_stats_gz_in_bucket", value = "n", N = lab_out[has_summary_gz == TRUE, .N]),
  data.table(metric = "lab_official_susie_complete_in_bucket", value = "n", N = lab_out[official_susie_in_bucket == TRUE, .N]),
  data.table(metric = "all_source_rows", value = "n", N = nrow(all_sources)),
  data.table(metric = "all_source_official_susie_complete_in_bucket", value = "n", N = all_sources[official_susie_in_bucket == TRUE, .N])
), fill = TRUE)
fwrite(qc, file.path(ROOT_DIR, "manifests/finngen_coloc_gwas_phenotype_standard_all_sources_qc.tsv"), sep = "\t", quote = FALSE)
print(qc)
