# Copy to config.R and edit paths before running this module.

CELL_TYPE <- Sys.getenv("SC_PCQTL_CELL_TYPE", unset = "")
if (!nzchar(CELL_TYPE)) stop("Set SC_PCQTL_CELL_TYPE explicitly.")
DATA_ROOT <- Sys.getenv("COQTL_DATA_ROOT", unset = "data")
WORK_DIR <- Sys.getenv("SC_PCQTL_HURDLE_WORK_DIR", unset = getwd())

COUNT_FILE <- Sys.getenv("COQTL_COUNT_FILE", unset = file.path(DATA_ROOT, sprintf("%s_readcounts.tsv", CELL_TYPE)))
GENE_INFO_FILE <- Sys.getenv("COQTL_GENE_INFO_FILE", unset = file.path(DATA_ROOT, "gene_info_with_location.tsv"))
MIN_CELLTYPE_CELLS <- 10000L

RESULTS_ROOT <- file.path(WORK_DIR, "results")
METHOD2_RESULTS <- file.path(RESULTS_ROOT, "method2_sc_hurdle")
LOG_DIR <- file.path(WORK_DIR, "logs")
