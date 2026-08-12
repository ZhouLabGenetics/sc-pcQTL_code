#!/usr/bin/env Rscript

.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(fasthurdle)
})

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- normalizePath(dirname(args_file), mustWork = TRUE)
count_file <- Sys.getenv("SC_PCQTL_GIMAP_COUNT_FILE", unset = "")
if (!nzchar(count_file) || !file.exists(count_file)) {
  stop("Set SC_PCQTL_GIMAP_COUNT_FILE to the final CD8_NC count/covariate matrix.")
}
library_cache <- Sys.getenv("SC_PCQTL_GIMAP_LIBRARY_CACHE", unset = "")
manuscript_root <- Sys.getenv("SC_PCQTL_MANUSCRIPT_ROOT", unset = "")
out_default <- if (nzchar(manuscript_root)) file.path(manuscript_root, "data", "gimap_pairwise_hurdle_matrix.tsv") else ""
out_file <- Sys.getenv("SC_PCQTL_GIMAP_HURDLE_MATRIX", unset = out_default)
if (!nzchar(out_file)) stop("Set SC_PCQTL_GIMAP_HURDLE_MATRIX or SC_PCQTL_MANUSCRIPT_ROOT.")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

genes <- c("GIMAP7", "GIMAP4", "GIMAP5", "GIMAP1", "GIMAP2", "GIMAP6", "GIMAP8")
header <- names(fread(count_file, nrows = 0))
id_col <- intersect(c("CellID", "barcode"), header)[1]
if (is.na(id_col)) stop("Count matrix needs a CellID or barcode column.")
required_covars <- c("age", "sex", paste0("pc", 1:6), "pf1", "pf2")
missing <- setdiff(c(required_covars, genes), header)
if (length(missing)) stop("Missing required GIMAP inputs: ", paste(missing, collapse = ", "))

select <- intersect(unique(c(id_col, required_covars, "log_total_read_counts", genes)), header)
dt <- fread(count_file, select = select)
if (!"log_total_read_counts" %in% names(dt)) {
  if (!nzchar(library_cache) || !file.exists(library_cache)) {
    stop("log_total_read_counts is absent; set SC_PCQTL_GIMAP_LIBRARY_CACHE.")
  }
  lib <- fread(library_cache, select = c(id_col, "log_total_read_counts"))
  dt <- merge(dt, lib, by = id_col, all.x = TRUE, sort = FALSE)
}
if (any(!is.finite(as.numeric(dt$log_total_read_counts)))) stop("Non-finite log_total_read_counts values.")

for (name in c("age", paste0("pc", 1:6), "pf1", "pf2", "log_total_read_counts")) {
  set(dt, j = name, value = as.numeric(dt[[name]]))
}
dt[, sex := factor(sex)]

active <- function(data, names) {
  names[vapply(names, function(name) uniqueN(data[[name]][!is.na(data[[name]])]) > 1L, logical(1))]
}
extract_z <- function(component, term) {
  if (is.null(component) || !term %in% rownames(component)) return(NA_real_)
  z_col <- grep("^z", colnames(component), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(z_col)) return(NA_real_)
  as.numeric(component[term, z_col])
}
fit_direction <- function(response, predictor) {
  model_dt <- data.frame(
    response = dt[[response]], predictor = dt[[predictor]],
    dt[, c(required_covars, "log_total_read_counts"), with = FALSE],
    check.names = FALSE
  )
  base_covars <- active(model_dt, required_covars)
  count_rhs <- paste(c("predictor", base_covars), collapse = " + ")
  zero_rhs <- paste(c("predictor", "log_total_read_counts", base_covars), collapse = " + ")
  formula <- as.formula(sprintf("response ~ %s + offset(log_total_read_counts) | %s", count_rhs, zero_rhs))
  fit <- tryCatch(
    fasthurdle(formula, data = model_dt, dist = "poisson", zero.dist = "binomial"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(c(count = NA_real_, zero = NA_real_))
  coef <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(coef)) return(c(count = NA_real_, zero = NA_real_))
  c(count = extract_z(coef$count, "predictor"), zero = extract_z(coef$zero, "predictor"))
}

result <- CJ(gene_row = genes, gene_col = genes, unique = TRUE)
result[, `:=`(z_count_abs = NA_real_, z_zero_abs = NA_real_)]
for (i in seq_len(length(genes) - 1L)) {
  for (j in (i + 1L):length(genes)) {
    forward <- fit_direction(genes[i], genes[j])
    reverse <- fit_direction(genes[j], genes[i])
    z_count <- suppressWarnings(max(abs(c(forward["count"], reverse["count"])), na.rm = TRUE))
    z_zero <- suppressWarnings(max(abs(c(forward["zero"], reverse["zero"])), na.rm = TRUE))
    if (!is.finite(z_count)) z_count <- NA_real_
    if (!is.finite(z_zero)) z_zero <- NA_real_
    result[
      (gene_row == genes[i] & gene_col == genes[j]) | (gene_row == genes[j] & gene_col == genes[i]),
      `:=`(z_count_abs = z_count, z_zero_abs = z_zero)
    ]
  }
}

fwrite(result, out_file, sep = "\t")
message("Wrote final GIMAP hurdle matrix to ", out_file)
