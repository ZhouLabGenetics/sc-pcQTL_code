# Copy this file into each cell-type pcQTL working directory or set
# SC_PCQTL_CELL_TYPE/SC_PCQTL_UPSTREAM_CELLTYPES_ROOT before running the module.
CELL_TYPE         <- Sys.getenv("SC_PCQTL_CELL_TYPE", unset = "")
if (!nzchar(CELL_TYPE)) stop("Set SC_PCQTL_CELL_TYPE explicitly.")
.config_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) "")
.module_dir <- if (nzchar(.config_file)) dirname(.config_file) else getwd()
.release_root <- normalizePath(file.path(.module_dir, "..", ".."), mustWork = TRUE)
source(file.path(.release_root, "config", "celltype_eligibility.R"))
.eligibility <- load_celltype_eligibility(file.path(.release_root, "config", "celltype_eligibility.tsv"))
.celltype_status <- .eligibility[celltype == CELL_TYPE]
if (!nrow(.celltype_status)) stop("Cell type is absent from the eligibility manifest: ", CELL_TYPE)
if (!.celltype_status$include_primary[[1]]) {
  stop("Cell type ", CELL_TYPE, " is excluded from primary analysis: ", .celltype_status$exclusion_reason[[1]])
}
DATA_ROOT         <- Sys.getenv("COQTL_DATA_ROOT", unset = "data")
UPSTREAM_ROOT     <- Sys.getenv("SC_PCQTL_UPSTREAM_ROOT", unset = normalizePath(file.path(getwd(), ".."), mustWork = FALSE))
CELLTYPES_ROOT    <- Sys.getenv("SC_PCQTL_UPSTREAM_CELLTYPES_ROOT", unset = file.path(UPSTREAM_ROOT, "celltypes"))
PCQTL_DIR         <- Sys.getenv("SC_PCQTL_PCQTL_DIR", unset = getwd())
COUNT_FILE        <- Sys.getenv("COQTL_COUNT_FILE", unset = file.path(DATA_ROOT, sprintf("%s_readcounts.tsv", CELL_TYPE)))
CLUSTER_FILE      <- Sys.getenv("SC_PCQTL_CLUSTER_FILE", unset = file.path(CELLTYPES_ROOT, CELL_TYPE, "cluster_identification/results/method2_sc_hurdle/merged_clusters/cluster_summary.tsv"))
CLUSTER_GENE_FILE <- Sys.getenv("SC_PCQTL_CLUSTER_GENE_FILE", unset = file.path(CELLTYPES_ROOT, CELL_TYPE, "cluster_identification/results/method2_sc_hurdle/merged_clusters/cluster_gene_assignments.tsv"))
