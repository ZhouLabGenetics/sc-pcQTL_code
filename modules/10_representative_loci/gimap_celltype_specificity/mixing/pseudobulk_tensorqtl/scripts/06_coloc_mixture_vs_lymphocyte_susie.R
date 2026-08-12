#!/usr/bin/env Rscript

# Fine-map each GIMAP pseudobulk readout with OneK1K donor LD, lift the QTL
# variants from hg19 to hg38, and compare it with the FinnGen R12 official
# SuSiE signal for lymphocyte count (phenocode 3019198) using coloc.susie.

require_dir <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value) || !dir.exists(value)) stop("Set ", name, " to an existing directory.")
  normalizePath(value, mustWork = TRUE)
}

code_root <- require_dir("SC_PCQTL_CODE_ROOT")
formal_root <- require_dir("SC_PCQTL_FORMAL_COLOC_ROOT")
mix_root <- require_dir("SC_PCQTL_GIMAP_MIXING_ROOT")
local_lib <- Sys.getenv("SC_PCQTL_R_LIBS", unset = "")
if (nzchar(local_lib)) .libPaths(c(local_lib, .libPaths()))

source(file.path(code_root, "modules/06_finngen_susie_coloc/config/config.R"))
source(file.path(code_root, "modules/06_finngen_susie_coloc/scripts/common.R"))
suppressPackageStartupMessages({
  library(coloc)
  library(susieR)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) default else args[hit[1] + 1L]
}

coloc_dir <- file.path(mix_root, "coloc_susie")
dir.create(coloc_dir, recursive = TRUE, showWarnings = FALSE)
sumstats_path <- get_arg("--sumstats", file.path(coloc_dir, "mix_all_sumstats.tsv"))
out_path <- get_arg("--out", file.path(coloc_dir, "mix_coloc_susie_summary.tsv"))
scenario_filter <- get_arg("--scenarios", "")
gwas_rds_path <- Sys.getenv(
  "SC_PCQTL_GIMAP_GWAS_SUSIE_RDS",
  unset = file.path(
    formal_root,
    "results/fine_mapping/gwas_all_finemapped/rds/3019198/cd8_nc/SC_chr7_cluster_001/gwas.official_finngen_susie.rds"
  )
)
ld_matrix_path <- file.path(coloc_dir, "gimap_union.ld.gz")
ld_vars_path <- file.path(coloc_dir, "gimap_union.ld_variants.tsv")
lab_manifest <- file.path(formal_root, "manifests/finngen_lab_values_phenotype_standard.tsv")

required_files <- c(sumstats_path, gwas_rds_path, ld_matrix_path, ld_vars_path, lab_manifest)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing GIMAP coloc inputs: ", paste(missing_files, collapse = ", "))

lift_cache <- new.env(parent = emptyenv())
load_liftover_posmap <- function(chr) {
  chr <- sub("^chr", "", as.character(chr))
  if (exists(chr, envir = lift_cache, inherits = FALSE)) return(get(chr, envir = lift_cache))
  path <- sprintf(ONEK1K_QTL_LIFTOVER_POSMAP_PATTERN, chr)
  if (!file.exists(path)) stop("Missing OneK1K QTL liftOver posmap: ", path)
  map <- fread(path, header = FALSE, col.names = c("old_id", "pos_hg38"))
  map[, `:=`(
    chr = sub("^chr", "", sub(":.*$", "", old_id)),
    pos_hg19 = as.integer(sub("^.*:", "", old_id)),
    pos_hg38 = as.integer(pos_hg38)
  )]
  map <- unique(
    map[is.finite(pos_hg19) & is.finite(pos_hg38), .(chr, pos_hg19, pos_hg38)],
    by = c("chr", "pos_hg19")
  )
  setkey(map, chr, pos_hg19)
  assign(chr, map, envir = lift_cache)
  map
}

lift_qtl_stats_to_hg38 <- function(qstats) {
  qstats <- as.data.table(copy(qstats))
  qstats[, `:=`(
    q_idx = .I,
    chr = sub("^chr", "", as.character(chr)),
    pos_hg19 = as.integer(pos),
    effect_allele = toupper(as.character(effect_allele)),
    other_allele = toupper(as.character(other_allele))
  )]
  maps <- rbindlist(lapply(unique(qstats$chr), load_liftover_posmap), fill = TRUE)
  lifted <- merge(qstats, maps, by = c("chr", "pos_hg19"), all = FALSE)
  if (!nrow(lifted)) return(lifted)
  lifted[, `:=`(
    pos = as.integer(pos_hg38),
    marker_id_hg19 = marker_id,
    marker_id = standard_variant_id(chr, pos_hg38, effect_allele, other_allele)
  )]
  lifted[, pos_hg38 := NULL]
  unique(lifted, by = c("q_idx", "chr", "pos", "effect_allele", "other_allele"))
}

shared_variant_table <- function(q, g) {
  qstats <- lift_qtl_stats_to_hg38(q$stats)
  gstats <- as.data.table(copy(g$stats))
  gstats[, `:=`(
    g_idx = .I,
    chr = sub("^chr", "", as.character(chr)),
    pos = as.integer(pos),
    effect_allele = toupper(as.character(effect_allele)),
    other_allele = toupper(as.character(other_allele))
  )]
  if (!nrow(qstats) || !nrow(gstats)) return(data.table())
  shared <- merge(
    qstats[, .(chr, pos, q_idx, q_effect = effect_allele, q_other = other_allele)],
    gstats[, .(chr, pos, g_idx, g_effect = effect_allele, g_other = other_allele)],
    by = c("chr", "pos"), allow.cartesian = TRUE
  )
  if (!nrow(shared)) return(data.table())
  shared[, allele_matched := q_effect == g_effect & q_other == g_other]
  shared[, allele_flipped := q_effect == g_other & q_other == g_effect]
  shared <- shared[allele_matched | allele_flipped]
  setorder(shared, chr, pos, g_effect, g_other)
  shared <- unique(shared, by = c("q_idx", "g_idx"))
  shared <- shared[!duplicated(q_idx) & !duplicated(g_idx)]
  shared[, snp := standard_variant_id(chr, pos, g_effect, g_other)]
  shared
}

set_susie_variant_names <- function(fit, snp_names) {
  snp_names <- make.unique(as.character(snp_names))
  for (field in c("lbf_variable", "alpha", "mu", "mu2")) {
    if (!is.null(fit[[field]]) && ncol(fit[[field]]) == length(snp_names)) {
      colnames(fit[[field]]) <- snp_names
    }
  }
  if (!is.null(fit$pip) && length(fit$pip) == length(snp_names)) names(fit$pip) <- snp_names
  fit
}

compact_susie_cs <- function(fit) {
  if (is.null(fit$sets$cs) || !length(fit$sets$cs) || is.null(fit$lbf_variable)) {
    return(list(fit = fit, original_signal_index = integer()))
  }
  original <- as.integer(fit$sets$cs_index %||% seq_along(fit$sets$cs))
  keep <- lengths(fit$sets$cs) > 0L & is.finite(original) &
    original >= 1L & original <= nrow(fit$lbf_variable)
  if (!any(keep)) return(list(fit = fit, original_signal_index = integer()))
  original <- original[keep]
  fit$sets$cs <- fit$sets$cs[keep]
  if (!is.null(fit$sets$coverage)) fit$sets$coverage <- fit$sets$coverage[keep]
  if (!is.null(fit$sets$purity)) fit$sets$purity <- fit$sets$purity[keep, , drop = FALSE]
  for (field in c("lbf_variable", "alpha", "mu", "mu2")) {
    if (!is.null(fit[[field]])) fit[[field]] <- fit[[field]][original, , drop = FALSE]
  }
  fit$sets$cs_index <- seq_along(original)
  list(fit = fit, original_signal_index = original)
}

lab <- fread(lab_manifest)
gwas_n <- suppressWarnings(as.numeric(lab[phenocode == "3019198"]$n_total[1]))
if (!is.finite(gwas_n) || gwas_n <= 0) stop("Missing FinnGen 3019198 sample size.")
gwas_rds <- readRDS(gwas_rds_path)
gwas <- list(fit = gwas_rds$fit, stats = as.data.table(gwas_rds$stats))

variants <- fread(ld_vars_path)
required_variant_cols <- c("chr", "variant_id", "cm", "pos", "other_allele", "effect_allele")
if (!all(required_variant_cols %in% names(variants))) stop("Invalid LD variant map: ", ld_vars_path)
variants[, chr := sub("^chr", "", as.character(chr))]
variants[, variant_key := standard_variant_id(chr, pos, effect_allele, other_allele)]
variants[, ld_order := .I]
ld_full <- read_ld_matrix(ld_matrix_path)
if (nrow(ld_full) != nrow(variants) || ncol(ld_full) != nrow(variants)) {
  stop("LD dimensions do not match the variant map.")
}

all_qtl <- fread(sumstats_path)
required_qtl_cols <- c("scenario", "qtl_type", "readout", "variant_id", "af", "slope", "slope_se", "qtl_n")
missing_qtl_cols <- setdiff(required_qtl_cols, names(all_qtl))
if (length(missing_qtl_cols)) stop("Missing QTL columns: ", paste(missing_qtl_cols, collapse = ", "))
if (nzchar(scenario_filter)) all_qtl <- all_qtl[scenario %in% strsplit(scenario_filter, ",")[[1]]]

coloc_one <- function(d) {
  stats <- merge(
    d,
    variants[, .(variant_id, chr, pos, tensorqtl_effect = other_allele, tensorqtl_other = effect_allele)],
    by = "variant_id"
  )
  stats <- stats[is.finite(slope) & is.finite(slope_se) & slope_se > 0 & af > 0 & af < 1]
  qtl_n <- unique(as.integer(stats$qtl_n))
  base <- list(
    pph4_susie = NA_real_, pph4_abf = NA_real_, n_shared = 0L,
    n_qtl_cs = 0L, qtl_n = if (length(qtl_n) == 1L) qtl_n else NA_integer_, status = "no_test"
  )
  if (length(qtl_n) != 1L || !is.finite(qtl_n) || qtl_n <= 0) {
    base$status <- "invalid_qtl_n"
    return(base)
  }
  if (nrow(stats) < 5L) return(base)

  # TensorQTL's slope is aligned to PLINK A1 (BIM column 5). The LD variant
  # map names this allele tensorqtl_effect and harmonizes it to BIM column 6.
  qstats <- data.table(
    chr = stats$chr,
    pos = as.integer(stats$pos),
    marker_id = stats$variant_id,
    effect_allele = stats$tensorqtl_effect,
    other_allele = stats$tensorqtl_other,
    beta = stats$slope,
    se = stats$slope_se,
    pvalue = 2 * pnorm(-abs(stats$slope / stats$slope_se)),
    af = stats$af,
    n = qtl_n
  )
  qstats <- harmonize_to_reference(
    qstats,
    variants[, .(chr, pos, effect_allele, other_allele)]
  )
  qstats <- merge(variants[, .(variant_key, ld_order)], qstats, by = "variant_key")
  setorder(qstats, ld_order)
  if (nrow(qstats) < 5L) return(base)
  ld_subset <- ld_full[qstats$ld_order, qstats$ld_order, drop = FALSE]
  fine_map <- tryCatch(
    run_susie_from_sumstats(qstats, ld_subset, SUSIE_L, SUSIE_COVERAGE, CS_MIN_ABS_CORR),
    error = function(e) NULL
  )
  if (is.null(fine_map)) {
    base$status <- "susie_failed"
    return(base)
  }
  qtl <- list(fit = fine_map$fit, stats = fine_map$stats)
  base$n_qtl_cs <- length(fine_map$cs$cs %||% list())
  shared <- shared_variant_table(qtl, gwas)
  base$n_shared <- nrow(shared)
  if (nrow(shared) < 2L) {
    base$status <- "few_shared"
    return(base)
  }

  q_shared <- qtl$stats[shared$q_idx]
  g_shared <- gwas$stats[shared$g_idx]
  base$pph4_abf <- tryCatch({
    dataset_qtl <- list(
      beta = q_shared$beta, varbeta = q_shared$se^2, snp = shared$snp,
      MAF = pmin(q_shared$af, 1 - q_shared$af), N = qtl_n, type = "quant"
    )
    dataset_gwas <- list(
      beta = g_shared$beta, varbeta = g_shared$se^2, snp = shared$snp,
      MAF = pmin(g_shared$af, 1 - g_shared$af), N = gwas_n, type = "quant"
    )
    invisible(capture.output(abf <- suppressMessages(coloc.abf(dataset_qtl, dataset_gwas))))
    as.numeric(abf$summary["PP.H4.abf"])
  }, error = function(e) NA_real_)

  q_fit <- compact_susie_cs(set_susie_variant_names(
    subset_susie_to_variants(qtl$fit, shared$q_idx), shared$snp
  ))
  g_fit <- compact_susie_cs(set_susie_variant_names(
    subset_susie_to_variants(gwas$fit, shared$g_idx), shared$snp
  ))
  if (!length(q_fit$original_signal_index) || !length(g_fit$original_signal_index)) {
    base$status <- "no_cs_after_share"
    return(base)
  }
  coloc_result <- tryCatch(
    coloc.susie(q_fit$fit, g_fit$fit, trim_by_posterior = FALSE),
    error = function(e) NULL
  )
  if (is.null(coloc_result) || !("summary" %in% names(coloc_result))) {
    base$status <- "coloc_no_pair"
    return(base)
  }
  coloc_summary <- as.data.table(coloc_result$summary)
  if (!nrow(coloc_summary)) {
    base$status <- "coloc_no_pair"
    return(base)
  }
  pph4_cols <- grep("PP\\.H4|PPH4", names(coloc_summary), value = TRUE)
  if (!length(pph4_cols)) {
    base$status <- "coloc_missing_pph4"
    return(base)
  }
  pph4 <- suppressWarnings(as.numeric(coloc_summary[[pph4_cols[1L]]]))
  if (!any(is.finite(pph4))) {
    base$status <- "coloc_nonfinite_pph4"
    return(base)
  }
  base$pph4_susie <- max(pph4[is.finite(pph4)])
  base$status <- "ok"
  base
}

result <- all_qtl[, {
  message(sprintf("[%s] %s %s", scenario[1], qtl_type[1], readout[1]))
  as.data.table(coloc_one(.SD))
}, by = .(scenario, qtl_type, readout)]
fwrite(result, out_path, sep = "\t")
message("Wrote ", nrow(result), " GIMAP mixing coloc rows to ", out_path)
