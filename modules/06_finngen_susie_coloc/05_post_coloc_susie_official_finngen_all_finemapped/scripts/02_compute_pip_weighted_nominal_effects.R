#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
module06_code_dir <- normalizePath(file.path(dirname(args_file), "..", ".."), mustWork = TRUE)
source(file.path(module06_code_dir, "config", "config.R"))
root_dir <- ROOT_DIR

post_dir <- file.path(root_dir, "05_post_coloc_susie_official_finngen_all_finemapped")
out_dir <- file.path(post_dir, "results/gene_effects")
cache_dir <- file.path(post_dir, "cache")
dir_create(out_dir)

cs_files <- list.files(
  file.path(root_dir, "results/fine_mapping/qtl"),
  pattern = "\\.credible_sets\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(cs_files)) stop("No QTL credible set files found")

message("Reading QTL credible sets")
cs <- rbindlist(lapply(cs_files, function(f) {
  x <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x) || !"credible_set_id" %in% names(x)) return(NULL)
  x
}), fill = TRUE)
if (!nrow(cs)) stop("No non-empty credible sets found")

needed <- c("celltype", "cluster_id", "phenotype_type", "phenotype_id",
            "credible_set_id", "chr", "pos", "effect_allele", "other_allele", "pip")
missing_needed <- setdiff(needed, names(cs))
if (length(missing_needed)) stop("Missing credible-set columns: ", paste(missing_needed, collapse = ", "))
cs <- cs[, ..needed]
cs <- cs[celltype %chin% names(CELLTYPE_EQTL_MAP)]
if (!nrow(cs)) stop("No primary-analysis credible sets remain after eligibility filtering")
cs[, `:=`(
  chr = as.character(chr),
  pos = as.integer(pos),
  effect_allele = toupper(as.character(effect_allele)),
  other_allele = toupper(as.character(other_allele)),
  pip = as.numeric(pip),
  phenotype_type = fifelse(phenotype_type == "pcQTL", "pcQTL", "eQTL")
)]
cs <- cs[is.finite(pip) & pip > 0]
setkey(cs, celltype, cluster_id)

message("Reading cluster gene assignments")
assign_files <- list.files(file.path(cache_dir, "cluster_gene_assignments"),
                           pattern = "\\.cluster_gene_assignments\\.tsv$", full.names = TRUE)
if (!length(assign_files)) stop("No cached cluster gene assignments found. Run 00_stage_upstream_inputs.sh first.")
gene_map <- rbindlist(lapply(assign_files, function(f) {
  ct <- sub("\\.cluster_gene_assignments\\.tsv$", "", basename(f))
  x <- fread(f, showProgress = FALSE)
  x[, celltype := ct]
  x
}), fill = TRUE)
if (!all(c("celltype", "cluster_id", "gene_name") %in% names(gene_map))) {
  stop("Cluster gene map is missing required columns")
}
gene_map <- unique(gene_map[, .(celltype, cluster_id, gene_name)])
gene_map <- gene_map[celltype %chin% names(CELLTYPE_EQTL_MAP)]
setkey(gene_map, celltype, cluster_id)

read_nominal_gene <- function(celltype, gene) {
  cached <- file.path(cache_dir, "nominal_eqtl", celltype, paste0(gene, ".singleVar.txt"))
  if (!file.exists(cached)) return(NULL)
  x <- tryCatch(fread(cached, showProgress = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  req <- c("CHR", "POS", "Allele1", "Allele2", "BETA")
  if (!all(req %in% names(x))) return(NULL)
  x <- x[, .(
    chr = as.character(CHR),
    pos = as.integer(POS),
    nominal_other = toupper(as.character(Allele1)),
    nominal_effect = toupper(as.character(Allele2)),
    nominal_beta = as.numeric(BETA)
  )]
  setkey(x, chr, pos)
  x
}

safe_cv <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  m <- mean(v)
  if (!is.finite(m) || m == 0) return(NA_real_)
  sd(v) / m
}

clusters <- unique(cs[, .(celltype, cluster_id)])
setorder(clusters, celltype, cluster_id)
gene_effects <- list()
missing_inputs <- list()
ri <- 0L
mi <- 0L

for (i in seq_len(nrow(clusters))) {
  ct <- clusters$celltype[i]
  cl <- clusters$cluster_id[i]
  if (i %% 100 == 0) message("Processed clusters: ", i, " / ", nrow(clusters))
  cs_sub <- cs[.(ct, cl)]
  genes <- unique(gene_map[.(ct, cl), gene_name])
  genes <- genes[!is.na(genes)]
  if (!length(genes) || !nrow(cs_sub)) next
  for (gene in genes) {
    nominal <- read_nominal_gene(ct, gene)
    if (is.null(nominal)) {
      mi <- mi + 1L
      missing_inputs[[mi]] <- data.table(celltype = ct, cluster_id = cl, gene_id = gene)
      next
    }
    merged <- merge(cs_sub, nominal, by = c("chr", "pos"), allow.cartesian = TRUE)
    if (!nrow(merged)) next
    merged[, allele_match := nominal_effect == effect_allele & nominal_other == other_allele]
    merged[, allele_flip := nominal_effect == other_allele & nominal_other == effect_allele]
    merged <- merged[allele_match | allele_flip]
    if (!nrow(merged)) next
    merged[, aligned_beta := fifelse(allele_match, nominal_beta, -nominal_beta)]
    eff <- merged[, .(
      pip_sum = sum(pip, na.rm = TRUE),
      pip_weighted_nominal_beta = sum(pip * aligned_beta, na.rm = TRUE) / sum(pip, na.rm = TRUE),
      n_matched_variants = uniqueN(paste(chr, pos, effect_allele, other_allele, sep = ":")),
      n_cs_variants = uniqueN(cs_sub[credible_set_id == .BY$credible_set_id, paste(chr, pos, effect_allele, other_allele, sep = ":")])
    ), by = .(celltype, cluster_id, qtl_type = phenotype_type, phenotype_id, credible_set_id)]
    if (nrow(eff)) {
      eff[, gene_id := gene]
      ri <- ri + 1L
      gene_effects[[ri]] <- eff
    }
  }
}

gene_effects <- rbindlist(gene_effects, fill = TRUE)
if (!nrow(gene_effects)) stop("No gene-level effects could be calculated")
gene_effects[, abs_pip_weighted_nominal_beta := abs(pip_weighted_nominal_beta)]
setcolorder(gene_effects, c(
  "celltype", "cluster_id", "qtl_type", "phenotype_id", "credible_set_id",
  "gene_id", "pip_sum", "pip_weighted_nominal_beta",
  "abs_pip_weighted_nominal_beta", "n_matched_variants", "n_cs_variants"
))
fwrite(gene_effects, file.path(out_dir, "pip_weighted_nominal_gene_effects.tsv"), sep = "\t")

summary_dt <- gene_effects[, {
  abs_eff <- abs_pip_weighted_nominal_beta[is.finite(abs_pip_weighted_nominal_beta)]
  sorted <- sort(abs_eff, decreasing = TRUE)
  .(
    n_genes = .N,
    n_genes_with_effect = sum(is.finite(abs_pip_weighted_nominal_beta)),
    pip_sum_median = median(pip_sum, na.rm = TRUE),
    max_abs_effect = if (length(sorted)) sorted[1] else NA_real_,
    mean_abs_effect = if (length(abs_eff)) mean(abs_eff) else NA_real_,
    sd_abs_effect = if (length(abs_eff) > 1) sd(abs_eff) else NA_real_,
    cv_effect = safe_cv(abs_eff),
    second_largest_abs_effect = if (length(sorted) >= 2) sorted[2] else NA_real_,
    fraction_largest = if (length(sorted) && sum(abs_eff) > 0) sorted[1] / sum(abs_eff) else NA_real_
  )
}, by = .(celltype, cluster_id, qtl_type, phenotype_id, credible_set_id)]
fwrite(summary_dt, file.path(out_dir, "qtl_cs_nominal_effect_summary.tsv"), sep = "\t")

missing_dt <- if (length(missing_inputs)) rbindlist(missing_inputs, fill = TRUE) else data.table()
if (!nrow(missing_dt)) missing_dt <- data.table(celltype = character(), cluster_id = character(), gene_id = character())
fwrite(missing_dt, file.path(out_dir, "gene_effect_missing_inputs.tsv"), sep = "\t")

message("Wrote gene-level effect outputs to ", out_dir)
message("Gene-effect rows: ", nrow(gene_effects), "; CS summary rows: ", nrow(summary_dt))
