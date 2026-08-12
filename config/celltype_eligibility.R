SC_PCQTL_MIN_CELLTYPE_CELLS <- 10000L
.SC_PCQTL_ELIGIBILITY_HELPER <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) "")
.SC_PCQTL_DEFAULT_CELLTYPE_MANIFEST <- if (nzchar(.SC_PCQTL_ELIGIBILITY_HELPER)) {
  file.path(dirname(.SC_PCQTL_ELIGIBILITY_HELPER), "celltype_eligibility.tsv")
} else {
  ""
}

load_celltype_eligibility <- function(manifest = Sys.getenv("SC_PCQTL_CELLTYPE_MANIFEST", unset = "")) {
  if (!nzchar(manifest)) {
    if (!nzchar(.SC_PCQTL_DEFAULT_CELLTYPE_MANIFEST)) stop("Set SC_PCQTL_CELLTYPE_MANIFEST explicitly.")
    manifest <- .SC_PCQTL_DEFAULT_CELLTYPE_MANIFEST
  }
  eligibility <- data.table::fread(manifest)
  required <- c("celltype", "eqtl_celltype", "display_label", "n_cells", "include_primary")
  missing <- setdiff(required, names(eligibility))
  if (length(missing)) stop("Cell-type eligibility manifest is missing: ", paste(missing, collapse = ", "))
  eligibility[, include_primary := as.logical(include_primary)]
  if (any(eligibility$include_primary & eligibility$n_cells < SC_PCQTL_MIN_CELLTYPE_CELLS)) {
    stop("Eligibility manifest includes a cell type below the 10,000-cell threshold.")
  }
  if (any(!eligibility$include_primary & eligibility$n_cells >= SC_PCQTL_MIN_CELLTYPE_CELLS)) {
    stop("Eligibility manifest excludes a cell type that meets the 10,000-cell threshold.")
  }
  eligibility
}

primary_celltypes <- function(manifest = Sys.getenv("SC_PCQTL_CELLTYPE_MANIFEST", unset = "")) {
  load_celltype_eligibility(manifest)[include_primary == TRUE, celltype]
}

canonical_celltype_map <- function(manifest = Sys.getenv("SC_PCQTL_CELLTYPE_MANIFEST", unset = ""),
                                   primary_only = FALSE) {
  eligibility <- load_celltype_eligibility(manifest)
  if (isTRUE(primary_only)) eligibility <- eligibility[include_primary == TRUE]
  stats::setNames(eligibility$eqtl_celltype, eligibility$celltype)
}

canonical_celltype_labels <- function(celltypes,
                                      manifest = Sys.getenv("SC_PCQTL_CELLTYPE_MANIFEST", unset = "")) {
  labels <- unname(canonical_celltype_map(manifest)[as.character(celltypes)])
  if (anyNA(labels)) {
    missing <- unique(as.character(celltypes)[is.na(labels)])
    stop("No canonical OneK1K label for cell type(s): ", paste(missing, collapse = ", "))
  }
  labels
}
