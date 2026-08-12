get_sim_root <- function() {
  Sys.getenv(
    "COQTL_SIM_ROOT",
    unset = normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
  )
}

get_data_root <- function() {
  Sys.getenv(
    "COQTL_SIM_DATA_ROOT",
    unset = file.path(get_sim_root(), "data")
  )
}

get_raw_counts_file <- function() {
  Sys.getenv(
    "COQTL_RAW_COUNTS_FILE",
    unset = ""
  )
}

get_gene_info_file <- function() {
  Sys.getenv(
    "COQTL_GENE_INFO_FILE",
    unset = ""
  )
}

load_matrix_long <- function(path, value_name = "Pvalue") {
  if (!file.exists(path)) {
    stop("Missing matrix: ", path)
  }
  dt <- data.table::fread(path)
  gene_order <- dt$Gene
  long <- data.table::melt(
    dt,
    id.vars = "Gene",
    variable.name = "Gene2",
    value.name = value_name
  )
  data.table::setnames(long, "Gene", "Gene1")
  long <- long[Gene1 != Gene2]
  long[, Gene1 := factor(Gene1, levels = gene_order)]
  long[, Gene2 := factor(Gene2, levels = gene_order)]
  long <- long[as.integer(Gene1) < as.integer(Gene2)]
  long[, Gene1 := as.character(Gene1)]
  long[, Gene2 := as.character(Gene2)]
  long[]
}

load_pair_mask_long <- function(path) {
  mask <- load_matrix_long(path, value_name = "pair_status")
  mask[, eligible := !is.na(pair_status) & pair_status == "tested"]
  mask[]
}
