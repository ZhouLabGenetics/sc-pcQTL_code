get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_idx <- grep("--file=", cmd_args)
  if (length(file_idx) > 0) {
    return(normalizePath(sub("--file=", "", cmd_args[file_idx])))
  }
  # fallback for interactive sessions
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
  normalizePath(".")
}

load_config <- function(script_dir = NULL) {
  if (is.null(script_dir)) {
    script_dir <- dirname(get_script_path())
  }
  config_path <- file.path(script_dir, "config.R")
  if (!file.exists(config_path)) {
    stop("Missing config file: ", config_path)
  }

  cfg_env <- new.env()
  sys.source(config_path, envir = cfg_env)

  required <- c(
    "CELL_TYPE", "DATA_ROOT", "WORK_DIR",
    "COUNT_FILE", "GENE_INFO_FILE", "RESULTS_ROOT",
    "METHOD2_RESULTS", "LOG_DIR", "MIN_CELLTYPE_CELLS"
  )
  missing <- required[!vapply(required, exists, logical(1), envir = cfg_env)]
  if (length(missing)) {
    stop("Config file missing definitions: ", paste(missing, collapse = ", "))
  }

  cfg <- list(
    cell_type = cfg_env$CELL_TYPE,
    data_root = cfg_env$DATA_ROOT,
    work_dir = cfg_env$WORK_DIR,
    count_file = cfg_env$COUNT_FILE,
    gene_info_file = cfg_env$GENE_INFO_FILE,
    results_root = cfg_env$RESULTS_ROOT,
    method2_results = cfg_env$METHOD2_RESULTS,
    logs_dir = cfg_env$LOG_DIR,
    min_celltype_cells = cfg_env$MIN_CELLTYPE_CELLS
  )

  release_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
  eligibility_env <- new.env(parent = baseenv())
  sys.source(file.path(release_root, "config", "celltype_eligibility.R"), envir = eligibility_env)
  eligibility_manifest <- Sys.getenv(
    "SC_PCQTL_CELLTYPE_MANIFEST",
    unset = file.path(release_root, "config", "celltype_eligibility.tsv")
  )
  eligibility <- eligibility_env$load_celltype_eligibility(
    eligibility_manifest
  )
  if (!identical(as.integer(cfg$min_celltype_cells), eligibility_env$SC_PCQTL_MIN_CELLTYPE_CELLS)) {
    stop(
      "MIN_CELLTYPE_CELLS must match the central eligibility threshold: ",
      eligibility_env$SC_PCQTL_MIN_CELLTYPE_CELLS
    )
  }
  status <- eligibility[celltype == cfg$cell_type]
  if (!nrow(status)) stop("Cell type is absent from the eligibility manifest: ", cfg$cell_type)
  if (!status$include_primary[[1]]) {
    stop("Cell type ", cfg$cell_type, " is excluded from primary analysis: ", status$exclusion_reason[[1]])
  }

  cfg$chunk_dir <- file.path(cfg$method2_results, "gene_associations_chunked")
  cfg$assoc_dir <- file.path(cfg$method2_results, "gene_associations")
  cfg$clusters_dir <- file.path(cfg$method2_results, "clusters_by_chr")
  cfg
}
