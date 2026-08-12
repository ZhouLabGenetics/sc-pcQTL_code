#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
})

script <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1])
root <- normalizePath(file.path(dirname(script), ".."))
source(file.path(root, "config", "config.R"))
source(file.path(root, "scripts", "common.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  hit <- which(args == flag)
  if (!length(hit) || hit == length(args)) return(default)
  args[hit + 1L]
}
task_file <- get_arg("--tasks", file.path(ROOT_DIR, "manifests/qtl_finemap_tasks.tsv"))
task_index <- as.integer(get_arg("--task-index", NA_character_))

tasks <- fread(task_file)
tasks <- tasks[status == "ready"]
if (!is.na(task_index)) tasks <- tasks[task_index]
if (!nrow(tasks)) stop("No ready QTL fine-mapping tasks")

dir_create(file.path(ROOT_DIR, "results/fine_mapping/qtl"))
status_rows <- list()

read_qtl_stats <- function(row) {
  if (row$phenotype_type == "pcQTL") {
    fread(row$sumstats_path)
  } else {
    cmd <- sprintf("tar -xOzf %s %s 2>/dev/null", shQuote(row$sumstats_path), shQuote(row$sumstats_member))
    dt <- fread(cmd = cmd)
    if (!ncol(dt)) {
      stop("Missing or empty eQTL tar member: ", row$sumstats_member)
    }
    dt
  }
}

for (i in seq_len(nrow(tasks))) {
  row <- tasks[i]
  out <- row$output_file
  out_rds <- sub("\\.credible_sets\\.tsv$", ".susie.rds", out)
  dir_create(dirname(out))

  status <- "ok"
  msg <- ""
  n_cs <- 0L
  tryCatch({
    if (!file.exists(row$ld_matrix) || !file.exists(row$ld_variants)) stop("missing LD files")
    region <- read_region_file(row$region_file)
    raw <- read_qtl_stats(row)
    stats <- normalize_sumstats(raw, "qtl")
    stats <- stats[chr == region$chr[1] & pos >= region$start[1] & pos <= region$end[1]]

    vars <- fread(row$ld_variants)
    setnames(vars, c("chr", "variant_id", "cm", "pos", "other_allele", "effect_allele"))
    vars[, chr := sub("^chr", "", as.character(chr))]
    vars[, variant_key := standard_variant_id(chr, pos, effect_allele, other_allele)]

    stats <- harmonize_to_reference(stats, vars[, .(chr, pos, effect_allele, other_allele)])
    stats <- merge(vars[, .(variant_key, ld_order = .I)], stats, by = "variant_key")
    setorder(stats, ld_order)

    ld <- read_ld_matrix(row$ld_matrix)
    keep <- stats$ld_order
    ld <- ld[keep, keep, drop = FALSE]
    # SAIGE-QTL N is the reported donor sample size, not the cell-row count.
    fm <- run_susie_from_sumstats(stats, ld, SUSIE_L, SUSIE_COVERAGE, CS_MIN_ABS_CORR)
    cs_dt <- extract_susie_cs(fm, row$phenotype_id, row$phenotype_type,
                              paste(row$celltype, row$cluster_id, sep = "__"))
    if (nrow(cs_dt)) {
      cs_dt[, `:=`(celltype = row$celltype, cluster_id = row$cluster_id, ld_source = "OneK1K_donor_genotype")]
    } else {
      cs_dt <- data.table(
        celltype = character(), cluster_id = character(), phenotype_type = character(),
        phenotype_id = character(), credible_set_id = character(), chr = character(),
        pos = integer(), effect_allele = character(), other_allele = character(),
        pip = numeric(), beta = numeric(), se = numeric(), pvalue = numeric(),
        af = numeric(), n = integer(), ld_source = character()
      )
    }
    n_cs <- uniqueN(cs_dt$credible_set_id)
    fwrite(cs_dt, out, sep = "\t", quote = FALSE)
    saveRDS(list(fit = fm$fit, stats = fm$stats, metadata = as.list(row)), out_rds)
  }, error = function(e) {
    status <<- "failed"
    msg <<- conditionMessage(e)
    fwrite(data.table(), out, sep = "\t")
  })

  status_rows[[i]] <- data.table(
    celltype = row$celltype,
    cluster_id = row$cluster_id,
    phenotype_type = row$phenotype_type,
    phenotype_id = row$phenotype_id,
    status = status,
    message = msg,
    n_cs = n_cs,
    output_file = out,
    susie_rds = out_rds
  )
  message(sprintf("[%d/%d] %s %s %s %s", i, nrow(tasks), status, row$celltype, row$cluster_id, row$phenotype_id))
}

status_dt <- rbindlist(status_rows, fill = TRUE)
status_file <- if (is.na(task_index)) {
  file.path(ROOT_DIR, "results/fine_mapping/qtl_finemap_status.tsv")
} else {
  file.path(ROOT_DIR, "results/fine_mapping", sprintf("qtl_finemap_status_task_%s.tsv", task_index))
}
dir_create(dirname(status_file))
fwrite(status_dt, status_file, sep = "\t", quote = FALSE)
