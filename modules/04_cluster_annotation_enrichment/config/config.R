suppressPackageStartupMessages(library(data.table))

.config_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) "")
.module_default <- if (nzchar(.config_file)) dirname(dirname(.config_file)) else getwd()
.release_root <- normalizePath(file.path(.module_default, "..", ".."), mustWork = TRUE)
source(file.path(.release_root, "config", "celltype_eligibility.R"))

ROOT_DIR <- normalizePath(
  Sys.getenv("SC_PCQTL_CLUSTER_ENRICHMENT_ROOT", unset = .module_default),
  mustWork = FALSE
)
ADD_COV_ROOT <- normalizePath(
  Sys.getenv(
    "SC_PCQTL_UPSTREAM_ROOT",
    unset = file.path(ROOT_DIR, "cache", "01_upstream_main_pipeline_add_cov")
  ),
  mustWork = FALSE
)
PB_INPUT_DIR <- normalizePath(
  Sys.getenv("SC_PCQTL_PB_INPUT_DIR", unset = file.path(ROOT_DIR, "cache", "pb_inputs")),
  mustWork = FALSE
)

CELLTYPES <- primary_celltypes(file.path(.release_root, "config", "celltype_eligibility.tsv"))

NON_GENE_COLS <- c(
  "barcode", "individual", "sex", paste0("pc", 1:6), "age", "pf1", "pf2",
  "cell_type", paste0("MOFA", 1:10),
  "cell_read_counts", "log_cell_read_counts", "log2_cell_read_counts",
  "total_read_counts", "log_total_read_counts",
  "phase", "G1.Score", "S.Score", "G2M.Score", "pseudotime", "pseudotime2"
)

pb_input_file <- function(celltype) {
  file.path(PB_INPUT_DIR, sprintf("%s_pb_input.tsv.gz", celltype))
}

MAIN_CLUSTER_SIZES <- 2:5
CLUSTER_THRESHOLD <- 0.70
BONF_ALPHA <- 0.05
PROMOTER_WINDOW_BP <- 1000L
BOUNDARY_WINDOW_BP <- 50000L

dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

std_chr <- function(x) sub("^chr", "", as.character(x))

safe_fread <- function(path, ...) {
  if (!file.exists(path) || file.info(path)$size == 0) return(data.table())
  if (grepl("\\.gz$", path)) return(fread(cmd = paste("zcat", shQuote(path)), ...))
  fread(path, ...)
}

fread_maybe_gz <- function(path, ...) {
  if (grepl("\\.gz$", path)) return(fread(cmd = paste("zcat", shQuote(path)), ...))
  fread(path, ...)
}
