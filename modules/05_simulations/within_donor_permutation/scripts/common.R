prepend_project_libraries <- function() {
  configured <- Sys.getenv("SC_PCQTL_R_LIBS", unset = "")
  if (!nzchar(configured)) return(invisible(NULL))
  paths <- strsplit(configured, .Platform$path.sep, fixed = TRUE)[[1L]]
  paths <- paths[nzchar(paths) & dir.exists(paths)]
  if (length(paths)) .libPaths(unique(c(paths, .libPaths())))
  invisible(NULL)
}

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args)) return(list())
  if (length(args) %% 2L != 0L) stop("Arguments must be supplied as --name value pairs")
  out <- list()
  for (i in seq.int(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[[i]])
    if (identical(key, args[[i]])) stop("Invalid argument: ", args[[i]])
    out[[key]] <- args[[i + 1L]]
  }
  out
}

require_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) stop("Missing required --", name)
  value
}

integer_arg <- function(args, name, default) {
  if (is.null(args[[name]])) return(as.integer(default))
  value <- suppressWarnings(as.integer(args[[name]]))
  if (is.na(value)) stop("--", name, " must be an integer")
  value
}

numeric_arg <- function(args, name, default) {
  if (is.null(args[[name]])) return(as.numeric(default))
  value <- suppressWarnings(as.numeric(args[[name]]))
  if (is.na(value)) stop("--", name, " must be numeric")
  value
}

atomic_write_table <- function(x, path, sep = "\t") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  con <- if (grepl("\\.gz$", path)) gzfile(tmp, open = "wt") else tmp
  write.table(
    x, con, sep = sep, quote = FALSE, row.names = FALSE, col.names = TRUE,
    na = "NA"
  )
  if (inherits(con, "connection")) close(con)
  if (!file.rename(tmp, path)) stop("Failed to move output into place: ", path)
  invisible(path)
}

atomic_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(x, tmp, compress = TRUE)
  if (!file.rename(tmp, path)) stop("Failed to move output into place: ", path)
  invisible(path)
}

safe_fraction <- function(x) {
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

extract_term_p <- function(coef_table, term) {
  if (is.null(coef_table) || !(term %in% rownames(coef_table))) return(NA_real_)
  p_col <- intersect(c("Pr(>|z|)", "Pr(>|t|)"), colnames(coef_table))
  if (!length(p_col)) return(NA_real_)
  as.numeric(coef_table[term, p_col[[1L]]])
}
