get_joint_score_module_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("--file=", args)
  if (length(hit)) return(dirname(normalizePath(sub("--file=", "", args[hit[1]]))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  normalizePath(".")
}

required_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Set ", name, " explicitly")
  value
}

load_joint_score_config <- function(module_dir = get_joint_score_module_dir()) {
  release_root <- normalizePath(file.path(module_dir, "..", ".."), mustWork = TRUE)
  manifest <- Sys.getenv(
    "SC_PCQTL_CELLTYPE_MANIFEST",
    unset = file.path(release_root, "config", "celltype_eligibility.tsv")
  )
  if (!nzchar(manifest)) manifest <- file.path(release_root, "config", "celltype_eligibility.tsv")

  eligibility_env <- new.env(parent = baseenv())
  sys.source(file.path(release_root, "config", "celltype_eligibility.R"), envir = eligibility_env)
  eligibility <- eligibility_env$load_celltype_eligibility(manifest)

  celltype_id <- required_env("SC_PCQTL_CELL_TYPE")
  status <- eligibility[celltype == celltype_id]
  if (!nrow(status)) stop("Cell type is absent from the eligibility manifest: ", celltype_id)
  if (!status$include_primary[[1]]) {
    stop("Cell type is excluded from primary analysis: ", status$exclusion_reason[[1]])
  }

  work_dir <- required_env("SC_PCQTL_JOINT_SCORE_WORK_DIR")
  method_dir <- file.path(work_dir, "results", "method2_sc_hurdle")
  cfg <- list(
    module_dir = module_dir,
    release_root = release_root,
    celltype = celltype_id,
    count_file = required_env("COQTL_COUNT_FILE"),
    gene_info_file = required_env("COQTL_GENE_INFO_FILE"),
    primary_filtered_genes_file = required_env("SC_PCQTL_PRIMARY_FILTERED_GENES_FILE"),
    library_size_cache = Sys.getenv("SC_PCQTL_LIBRARY_SIZE_CACHE", unset = ""),
    work_dir = work_dir,
    method_dir = method_dir,
    stage_dir = file.path(method_dir, "joint_score_stage"),
    chunk_dir = file.path(method_dir, "gene_associations_chunked"),
    assoc_dir = file.path(method_dir, "gene_associations"),
    clusters_dir = file.path(method_dir, "clusters_by_chr"),
    merged_dir = file.path(method_dir, "merged_clusters"),
    min_celltype_cells = eligibility_env$SC_PCQTL_MIN_CELLTYPE_CELLS
  )
  cfg
}
