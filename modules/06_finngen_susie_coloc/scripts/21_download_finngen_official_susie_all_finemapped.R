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

manifest <- get_arg("--manifest", file.path(ROOT_DIR, "manifests/finngen_coloc_gwas_phenotype_standard_all_sources.tsv"))
chunk_index <- as.integer(get_arg("--chunk-index", NA_character_))
n_chunks <- as.integer(get_arg("--n-chunks", NA_character_))
out_dir <- get_arg("--out-dir", file.path(ROOT_DIR, "logs/finngen_official_susie_all_finemapped_download"))
max_phenotypes <- as.integer(get_arg("--max-phenotypes", NA_character_))
phenocode_filter <- get_arg("--phenocode", NA_character_)

if (!file.exists(manifest)) stop("Missing all-source FinnGen phenotype manifest: ", manifest)
dir_create(out_dir)

endpoints <- fread(manifest)
required <- c(
  "phenocode", "endpoint_source", "official_susie_in_bucket",
  "susie_cred_url", "susie_snp_url", "susie_snp_tbi_url",
  "local_susie_cred", "local_susie_snp", "local_susie_tbi"
)
miss <- setdiff(required, names(endpoints))
if (length(miss)) stop("Missing manifest columns: ", paste(miss, collapse = ", "))

endpoints <- endpoints[official_susie_in_bucket == TRUE]
if (!is.na(phenocode_filter)) endpoints <- endpoints[phenocode == phenocode_filter]
setorder(endpoints, endpoint_source, phenocode)
if (!is.na(max_phenotypes)) endpoints <- head(endpoints, max_phenotypes)

if (!is.na(chunk_index) || !is.na(n_chunks)) {
  if (is.na(chunk_index) || is.na(n_chunks) || chunk_index < 1L || n_chunks < 1L || chunk_index > n_chunks) {
    stop("--chunk-index and --n-chunks must be valid positive integers")
  }
  chunk_size <- ceiling(nrow(endpoints) / n_chunks)
  start_i <- (chunk_index - 1L) * chunk_size + 1L
  end_i <- min(chunk_index * chunk_size, nrow(endpoints))
  endpoints <- if (start_i <= nrow(endpoints)) endpoints[start_i:end_i] else endpoints[0]
}

download_one <- function(url, out_file) {
  if (file.exists(out_file) && file.info(out_file)$size > 0) return("exists")
  dir_create(dirname(out_file))
  tmp <- paste0(out_file, ".tmp.", Sys.getpid())
  curl_bin <- if (file.exists("/usr/bin/curl")) "/usr/bin/curl" else unname(Sys.which("curl"))
  wget_bin <- if (file.exists("/usr/bin/wget")) "/usr/bin/wget" else unname(Sys.which("wget"))
  if (nzchar(curl_bin)) {
    status <- system2(curl_bin, c("-sS", "-k", "-L", "--fail", "--retry", "4", "--retry-delay", "10", "-o", tmp, url))
  } else if (nzchar(wget_bin)) {
    status <- system2(wget_bin, c("--quiet", "--no-check-certificate", "--tries=4", "-O", tmp, url))
  } else {
    status <- tryCatch({
      download.file(url, tmp, mode = "wb", quiet = TRUE)
      0L
    }, error = function(e) 1L)
  }
  if (status != 0L || !file.exists(tmp) || file.info(tmp)$size == 0) {
    unlink(tmp)
    return("failed")
  }
  file.rename(tmp, out_file)
  "downloaded"
}

rows <- list()
row_i <- 0L
for (i in seq_len(nrow(endpoints))) {
  ep <- endpoints[i]
  specs <- data.table(
    file_type = c("cred", "snp", "snp_tbi"),
    url = c(ep$susie_cred_url, ep$susie_snp_url, ep$susie_snp_tbi_url),
    local_file = c(ep$local_susie_cred, ep$local_susie_snp, ep$local_susie_tbi)
  )
  for (j in seq_len(nrow(specs))) {
    status <- download_one(specs$url[j], specs$local_file[j])
    row_i <- row_i + 1L
    rows[[row_i]] <- data.table(
      phenocode = ep$phenocode,
      endpoint_source = ep$endpoint_source,
      file_type = specs$file_type[j],
      status = status,
      local_file = specs$local_file[j],
      bytes = if (file.exists(specs$local_file[j])) file.info(specs$local_file[j])$size else NA_real_
    )
  }
  message(sprintf("[%d/%d] %s (%s)", i, nrow(endpoints), ep$phenocode, ep$endpoint_source))
}

dt <- if (length(rows)) rbindlist(rows, fill = TRUE) else data.table()
suffix <- if (!is.na(chunk_index)) sprintf("_chunk_%03d", chunk_index) else ""
fwrite(dt, file.path(out_dir, paste0("finngen_official_susie_all_finemapped_download", suffix, ".tsv")), sep = "\t", quote = FALSE)
print(dt[, .N, by = .(endpoint_source, file_type, status)][order(endpoint_source, file_type, status)])
