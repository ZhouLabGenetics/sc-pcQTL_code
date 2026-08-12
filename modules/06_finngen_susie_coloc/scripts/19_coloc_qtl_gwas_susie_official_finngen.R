#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

script <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1])
root <- normalizePath(file.path(dirname(script), ".."))
source(file.path(root, "config", "config.R"))
source(file.path(root, "scripts", "common.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}
has_flag <- function(flag) flag %in% args

pair_manifest <- get_arg("--pair-manifest", file.path(ROOT_DIR, "manifests/qtl_gwas_susie_official_finngen_pairs.tsv"))
gwas_manifest <- get_arg("--gwas-manifest", file.path(ROOT_DIR, "results/fine_mapping/gwas/finngen_official_cs_by_cluster.tsv"))
chunk_index <- as.integer(get_arg("--chunk-index", NA_character_))
n_chunks <- as.integer(get_arg("--n-chunks", NA_character_))
max_pairs <- as.integer(get_arg("--max-pairs", NA_character_))
out_dir <- get_arg("--out-dir", file.path(ROOT_DIR, "results/coloc/qtl_gwas_susie_official_finngen"))
primary_ct <- names(CELLTYPE_EQTL_MAP)

set_susie_variant_names <- function(fit, snp_names) {
  snp_names <- make.unique(as.character(snp_names))
  if (!is.null(fit$lbf_variable) && ncol(fit$lbf_variable) == length(snp_names)) colnames(fit$lbf_variable) <- snp_names
  if (!is.null(fit$alpha) && ncol(fit$alpha) == length(snp_names)) colnames(fit$alpha) <- snp_names
  if (!is.null(fit$mu) && ncol(fit$mu) == length(snp_names)) colnames(fit$mu) <- snp_names
  if (!is.null(fit$mu2) && ncol(fit$mu2) == length(snp_names)) colnames(fit$mu2) <- snp_names
  if (!is.null(fit$pip) && length(fit$pip) == length(snp_names)) names(fit$pip) <- snp_names
  fit
}

compact_susie_cs <- function(fit) {
  if (is.null(fit$sets$cs) || !length(fit$sets$cs) || is.null(fit$lbf_variable)) {
    return(list(fit = fit, original_signal_index = integer()))
  }
  orig_idx <- fit$sets$cs_index %||% seq_along(fit$sets$cs)
  orig_idx <- as.integer(orig_idx)
  keep <- lengths(fit$sets$cs) > 0L & is.finite(orig_idx) & orig_idx >= 1L & orig_idx <= nrow(fit$lbf_variable)
  if (!any(keep)) return(list(fit = fit, original_signal_index = integer()))
  orig_idx <- orig_idx[keep]
  fit$sets$cs <- fit$sets$cs[keep]
  if (!is.null(fit$sets$coverage)) fit$sets$coverage <- fit$sets$coverage[keep]
  if (!is.null(fit$sets$purity)) fit$sets$purity <- fit$sets$purity[keep, , drop = FALSE]
  if (!is.null(fit$lbf_variable)) fit$lbf_variable <- fit$lbf_variable[orig_idx, , drop = FALSE]
  if (!is.null(fit$alpha)) fit$alpha <- fit$alpha[orig_idx, , drop = FALSE]
  if (!is.null(fit$mu)) fit$mu <- fit$mu[orig_idx, , drop = FALSE]
  if (!is.null(fit$mu2)) fit$mu2 <- fit$mu2[orig_idx, , drop = FALSE]
  fit$sets$cs_index <- seq_along(orig_idx)
  list(fit = fit, original_signal_index = orig_idx)
}

alpha_mass <- function(fit, signal_idx, keep_idx) {
  if (is.null(fit$alpha) || signal_idx > nrow(fit$alpha)) return(NA_real_)
  denom <- sum(fit$alpha[signal_idx, ], na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  sum(fit$alpha[signal_idx, keep_idx], na.rm = TRUE) / denom
}

pip_mass <- function(fit, keep_idx) {
  if (is.null(fit$pip)) return(NA_real_)
  denom <- sum(fit$pip, na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  sum(fit$pip[keep_idx], na.rm = TRUE) / denom
}

lift_cache <- new.env(parent = emptyenv())
load_liftover_posmap <- function(chr) {
  chr <- sub("^chr", "", as.character(chr))
  if (exists(chr, envir = lift_cache, inherits = FALSE)) return(get(chr, envir = lift_cache))
  path <- sprintf(ONEK1K_QTL_LIFTOVER_POSMAP_PATTERN, chr)
  if (!file.exists(path)) stop("Missing OneK1K QTL liftOver posmap: ", path)
  dt <- fread(path, header = FALSE, col.names = c("old_id", "pos_hg38"))
  dt[, `:=`(
    chr = sub("^chr", "", sub(":.*$", "", old_id)),
    pos_hg19 = as.integer(sub("^.*:", "", old_id)),
    pos_hg38 = as.integer(pos_hg38)
  )]
  dt <- unique(dt[is.finite(pos_hg19) & is.finite(pos_hg38), .(chr, pos_hg19, pos_hg38)], by = c("chr", "pos_hg19"))
  setkey(dt, chr, pos_hg19)
  assign(chr, dt, envir = lift_cache)
  dt
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
  lifted <- merge(qstats, maps, by = c("chr", "pos_hg19"), all.x = FALSE, all.y = FALSE)
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
  qstats0 <- as.data.table(copy(q$stats))
  gstats <- as.data.table(copy(g$stats))
  qstats <- lift_qtl_stats_to_hg38(qstats0)
  gstats[, `:=`(
    g_idx = .I,
    chr = sub("^chr", "", as.character(chr)),
    pos = as.integer(pos),
    effect_allele = toupper(as.character(effect_allele)),
    other_allele = toupper(as.character(other_allele))
  )]
  if (!nrow(qstats) || !nrow(gstats)) {
    return(list(shared = data.table(), qc = data.table(
      n_qtl_variants_before_liftover = nrow(qstats0),
      n_qtl_variants_after_liftover = nrow(qstats),
      n_gwas_variants = nrow(gstats),
      n_position_overlaps = 0L,
      n_allele_matched = 0L,
      n_allele_flipped = 0L,
      n_shared_variants = 0L
    )))
  }
  m <- merge(
    qstats[, .(chr, pos, q_idx, q_effect = effect_allele, q_other = other_allele, q_pos_hg19 = pos_hg19)],
    gstats[, .(chr, pos, g_idx, g_effect = effect_allele, g_other = other_allele)],
    by = c("chr", "pos"),
    allow.cartesian = TRUE
  )
  if (!nrow(m)) {
    return(list(shared = data.table(), qc = data.table(
      n_qtl_variants_before_liftover = nrow(qstats0),
      n_qtl_variants_after_liftover = nrow(qstats),
      n_gwas_variants = nrow(gstats),
      n_position_overlaps = 0L,
      n_allele_matched = 0L,
      n_allele_flipped = 0L,
      n_shared_variants = 0L
    )))
  }
  m[, allele_matched := q_effect == g_effect & q_other == g_other]
  m[, allele_flipped := q_effect == g_other & q_other == g_effect]
  allele_matched_n <- sum(m$allele_matched, na.rm = TRUE)
  allele_flipped_n <- sum(m$allele_flipped, na.rm = TRUE)
  shared <- m[allele_matched | allele_flipped]
  setorder(shared, chr, pos, g_effect, g_other)
  shared <- unique(shared, by = c("q_idx", "g_idx"))
  shared <- shared[!duplicated(q_idx) & !duplicated(g_idx)]
  shared[, snp := standard_variant_id(chr, pos, g_effect, g_other)]
  list(shared = shared, qc = data.table(
    n_qtl_variants_before_liftover = nrow(qstats0),
    n_qtl_variants_after_liftover = nrow(qstats),
    n_gwas_variants = nrow(gstats),
    n_position_overlaps = nrow(m),
    n_allele_matched = allele_matched_n,
    n_allele_flipped = allele_flipped_n,
    n_shared_variants = nrow(shared)
  ))
}

build_pairs <- function() {
  qtl_status <- fread(file.path(ROOT_DIR, "results/fine_mapping/qtl_finemap_status.tsv"))
  qtl_status <- qtl_status[celltype %chin% primary_ct & status == "ok" & n_cs > 0 & file.exists(susie_rds)]
  if (!file.exists(gwas_manifest)) stop("Missing FinnGen official SuSiE manifest: ", gwas_manifest)
  gwas_meta <- unique(fread(
    gwas_manifest,
    select = c("phenocode", "phenotype", "celltype", "cluster_id",
               "gwas_susie_rds", "ld_source", "gwas_finemap_source")
  ))
  setnames(gwas_meta, "ld_source", "gwas_ld_source", skip_absent = TRUE)
  gwas_meta <- gwas_meta[celltype %chin% primary_ct & file.exists(gwas_susie_rds)]
  pairs <- merge(qtl_status, gwas_meta, by = c("celltype", "cluster_id"), allow.cartesian = TRUE)
  setorder(pairs, phenocode, celltype, cluster_id, phenotype_type, phenotype_id)
  pairs
}

if (file.exists(pair_manifest) && !has_flag("--force-rebuild-pairs")) {
  pairs <- fread(pair_manifest)
  if (!"celltype" %in% names(pairs)) stop("Existing QTL-GWAS pair manifest lacks celltype")
  n_pairs_before <- nrow(pairs)
  pairs <- pairs[celltype %chin% primary_ct]
  if (nrow(pairs) < n_pairs_before) {
    message("Purged ", n_pairs_before - nrow(pairs), " excluded-celltype rows from existing pair manifest")
    fwrite(pairs, pair_manifest, sep = "\t", quote = FALSE)
  }
} else {
  pairs <- build_pairs()
  dir_create(dirname(pair_manifest))
  fwrite(pairs, pair_manifest, sep = "\t", quote = FALSE)
}
if (!nrow(pairs)) stop("No primary-analysis QTL-GWAS pairs remain")
if (!is.na(max_pairs)) pairs <- head(pairs, max_pairs)

if (has_flag("--build-pair-manifest-only")) {
  message("QTL-GWAS SuSiE official FinnGen pairs: ", nrow(pairs))
  quit(save = "no")
}

if (!is.na(chunk_index) || !is.na(n_chunks)) {
  if (is.na(chunk_index) || is.na(n_chunks) || chunk_index < 1L || n_chunks < 1L || chunk_index > n_chunks) {
    stop("--chunk-index and --n-chunks must be valid positive integers")
  }
  chunk_size <- ceiling(nrow(pairs) / n_chunks)
  start_i <- (chunk_index - 1L) * chunk_size + 1L
  end_i <- min(chunk_index * chunk_size, nrow(pairs))
  pairs <- if (start_i <= nrow(pairs)) pairs[start_i:end_i] else pairs[0]
}

dir_create(out_dir)
rows <- list()
qc_rows <- list()
row_i <- 0L

for (i in seq_len(nrow(pairs))) {
  row <- pairs[i]
  base_row <- data.table(
    phenocode = row$phenocode,
    phenotype = row$phenotype,
    celltype = row$celltype,
    cluster_id = row$cluster_id,
    qtl_type = row$phenotype_type,
    qtl_phenotype_id = row$phenotype_id,
    qtl_susie_rds = row$susie_rds,
    gwas_susie_rds = row$gwas_susie_rds,
    gwas_finemap_source = row$gwas_finemap_source %||% "FinnGen_R12_official_SuSiE",
    gwas_ld_source = row$gwas_ld_source %||% "FinnGen_R12_official_in_sample_LD"
  )
  tryCatch({
    q <- readRDS(row$susie_rds)
    g <- readRDS(row$gwas_susie_rds)
    sv <- shared_variant_table(q, g)
    shared <- sv$shared
    qc <- cbind(base_row, sv$qc)
    qc_rows[[length(qc_rows) + 1L]] <- qc

    if (nrow(shared) < 2L) {
      row_i <- row_i + 1L
      rows[[row_i]] <- cbind(base_row, data.table(
        status = "no_test",
        message = "fewer than two allele-compatible shared variants after QTL hg19-to-hg38 liftOver",
        qtl_signal_index = NA_integer_,
        gwas_signal_index = NA_integer_,
        n_shared_variants = nrow(shared),
        n_coloc_variants = 0L,
        qtl_shared_alpha_mass = NA_real_,
        gwas_shared_alpha_mass = NA_real_,
        qtl_shared_pip_mass = NA_real_,
        gwas_shared_pip_mass = NA_real_,
        PP.H0.abf = NA_real_, PP.H1.abf = NA_real_, PP.H2.abf = NA_real_,
        PP.H3.abf = NA_real_, PP.H4.abf = NA_real_,
        pph4 = NA_real_,
        shared_mass_pass = FALSE,
        colocalized = FALSE,
        colocalized_h4_0_7 = FALSE,
        colocalized_h4_0_8 = FALSE
      ))
      next
    }

    q_shared_fit <- subset_susie_to_variants(q$fit, shared$q_idx)
    g_shared_fit <- subset_susie_to_variants(g$fit, shared$g_idx)
    q_shared_fit <- set_susie_variant_names(q_shared_fit, shared$snp)
    g_shared_fit <- set_susie_variant_names(g_shared_fit, shared$snp)
    q_compact <- compact_susie_cs(q_shared_fit)
    g_compact <- compact_susie_cs(g_shared_fit)

    if (!length(q_compact$original_signal_index) || !length(g_compact$original_signal_index)) {
      row_i <- row_i + 1L
      rows[[row_i]] <- cbind(base_row, data.table(
        status = "no_test",
        message = "no non-empty QTL or GWAS SuSiE CS after shared variant subsetting",
        qtl_signal_index = NA_integer_,
        gwas_signal_index = NA_integer_,
        n_shared_variants = nrow(shared),
        n_coloc_variants = 0L,
        qtl_shared_alpha_mass = NA_real_,
        gwas_shared_alpha_mass = NA_real_,
        qtl_shared_pip_mass = pip_mass(q$fit, shared$q_idx),
        gwas_shared_pip_mass = pip_mass(g$fit, shared$g_idx),
        PP.H0.abf = NA_real_, PP.H1.abf = NA_real_, PP.H2.abf = NA_real_,
        PP.H3.abf = NA_real_, PP.H4.abf = NA_real_,
        pph4 = NA_real_,
        shared_mass_pass = FALSE,
        colocalized = FALSE,
        colocalized_h4_0_7 = FALSE,
        colocalized_h4_0_8 = FALSE
      ))
      next
    }

    cr <- coloc.susie(q_compact$fit, g_compact$fit, trim_by_posterior = FALSE)
    res <- if ("summary" %in% names(cr)) as.data.table(cr$summary) else as.data.table(cr)
    if (!nrow(res) || !"idx1" %in% names(res) || !"idx2" %in% names(res)) {
      row_i <- row_i + 1L
      rows[[row_i]] <- cbind(base_row, data.table(
        status = "no_test",
        message = "coloc.susie returned no signal-pair summary",
        qtl_signal_index = NA_integer_,
        gwas_signal_index = NA_integer_,
        n_shared_variants = nrow(shared),
        n_coloc_variants = nrow(shared),
        qtl_shared_alpha_mass = NA_real_,
        gwas_shared_alpha_mass = NA_real_,
        qtl_shared_pip_mass = pip_mass(q$fit, shared$q_idx),
        gwas_shared_pip_mass = pip_mass(g$fit, shared$g_idx),
        PP.H0.abf = NA_real_, PP.H1.abf = NA_real_, PP.H2.abf = NA_real_,
        PP.H3.abf = NA_real_, PP.H4.abf = NA_real_,
        pph4 = NA_real_,
        shared_mass_pass = FALSE,
        colocalized = FALSE,
        colocalized_h4_0_7 = FALSE,
        colocalized_h4_0_8 = FALSE
      ))
      next
    }

    pp_col <- grep("PP\\.H4|PPH4", names(res), value = TRUE)[1]
    if (is.na(pp_col)) stop("Cannot find PPH4 column in coloc.susie result")
    for (j in seq_len(nrow(res))) {
      qsig <- q_compact$original_signal_index[as.integer(res$idx1[j])]
      gsig <- g_compact$original_signal_index[as.integer(res$idx2[j])]
      q_alpha_mass <- alpha_mass(q$fit, qsig, shared$q_idx)
      g_alpha_mass <- alpha_mass(g$fit, gsig, shared$g_idx)
      shared_pass <- isTRUE(q_alpha_mass >= COLOC_SHARED_ALPHA_MASS_MIN) &&
        isTRUE(g_alpha_mass >= COLOC_SHARED_ALPHA_MASS_MIN)
      pph4 <- as.numeric(res[[pp_col]][j])
      row_i <- row_i + 1L
      rows[[row_i]] <- cbind(base_row, data.table(
        status = "ok",
        message = "",
        qtl_signal_index = qsig,
        gwas_signal_index = gsig,
        n_shared_variants = nrow(shared),
        n_coloc_variants = as.integer(res$nsnps[j] %||% nrow(shared)),
        qtl_shared_alpha_mass = q_alpha_mass,
        gwas_shared_alpha_mass = g_alpha_mass,
        qtl_shared_pip_mass = pip_mass(q$fit, shared$q_idx),
        gwas_shared_pip_mass = pip_mass(g$fit, shared$g_idx),
        PP.H0.abf = as.numeric(res$PP.H0.abf[j] %||% NA_real_),
        PP.H1.abf = as.numeric(res$PP.H1.abf[j] %||% NA_real_),
        PP.H2.abf = as.numeric(res$PP.H2.abf[j] %||% NA_real_),
        PP.H3.abf = as.numeric(res$PP.H3.abf[j] %||% NA_real_),
        PP.H4.abf = pph4,
        pph4 = pph4,
        shared_mass_pass = shared_pass,
        colocalized = is.finite(pph4) && pph4 > COLOC_PPH4_CUTOFF,
        colocalized_h4_0_7 = is.finite(pph4) && pph4 > 0.7,
        colocalized_h4_0_8 = is.finite(pph4) && pph4 > COLOC_PPH4_SENSITIVITY_CUTOFF,
        colocalized_shared_mass_qc = shared_pass && is.finite(pph4) && pph4 > COLOC_PPH4_CUTOFF,
        colocalized_h4_0_7_shared_mass_qc = shared_pass && is.finite(pph4) && pph4 > 0.7,
        colocalized_h4_0_8_shared_mass_qc = shared_pass && is.finite(pph4) && pph4 > COLOC_PPH4_SENSITIVITY_CUTOFF
      ))
    }
  }, error = function(e) {
    row_i <<- row_i + 1L
    rows[[row_i]] <<- cbind(base_row, data.table(
      status = "failed",
      message = conditionMessage(e),
      qtl_signal_index = NA_integer_,
      gwas_signal_index = NA_integer_,
      n_shared_variants = NA_integer_,
      n_coloc_variants = NA_integer_,
      qtl_shared_alpha_mass = NA_real_,
      gwas_shared_alpha_mass = NA_real_,
      qtl_shared_pip_mass = NA_real_,
      gwas_shared_pip_mass = NA_real_,
      PP.H0.abf = NA_real_, PP.H1.abf = NA_real_, PP.H2.abf = NA_real_,
      PP.H3.abf = NA_real_, PP.H4.abf = NA_real_,
      pph4 = NA_real_,
      shared_mass_pass = FALSE,
      colocalized = FALSE,
      colocalized_h4_0_7 = FALSE,
      colocalized_h4_0_8 = FALSE
    ))
  })
  if (i %% 100 == 0) message(sprintf("[%d/%d] SuSiE official FinnGen coloc pairs processed", i, nrow(pairs)))
}

dt <- rbindlist(rows, fill = TRUE)
qc_dt <- rbindlist(qc_rows, fill = TRUE)
if (!is.na(chunk_index)) {
  chunk_dir <- file.path(out_dir, "chunks")
  dir_create(chunk_dir)
  fwrite(dt, file.path(chunk_dir, sprintf("qtl_gwas_susie_official_finngen_coloc_chunk_%03d.tsv", chunk_index)),
         sep = "\t", quote = FALSE)
  fwrite(qc_dt, file.path(chunk_dir, sprintf("qtl_gwas_susie_official_finngen_qc_chunk_%03d.tsv", chunk_index)),
         sep = "\t", quote = FALSE)
} else {
  fwrite(dt, file.path(out_dir, "qtl_gwas_susie_official_finngen_coloc_summary.tsv"), sep = "\t", quote = FALSE)
  fwrite(qc_dt, file.path(out_dir, "qtl_gwas_susie_official_finngen_liftover_qc.tsv"), sep = "\t", quote = FALSE)
}
