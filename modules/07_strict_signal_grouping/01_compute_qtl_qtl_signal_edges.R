#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (is.na(args_file) || !nzchar(args_file)) args_file <- getwd()
module06_code_dir <- normalizePath(file.path(dirname(args_file), "..", "06_finngen_susie_coloc"), mustWork = TRUE)
source(file.path(module06_code_dir, "config", "config.R"))
source(file.path(module06_code_dir, "scripts", "common.R"))
root_dir <- ROOT_DIR

out_dir <- file.path(root_dir, "06_strict_signal_grouping/results/qtl_qtl")
dir_create(out_dir)

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
  fit$lbf_variable <- fit$lbf_variable[orig_idx, , drop = FALSE]
  if (!is.null(fit$alpha)) fit$alpha <- fit$alpha[orig_idx, , drop = FALSE]
  if (!is.null(fit$mu)) fit$mu <- fit$mu[orig_idx, , drop = FALSE]
  if (!is.null(fit$mu2)) fit$mu2 <- fit$mu2[orig_idx, , drop = FALSE]
  fit$sets$cs_index <- seq_along(orig_idx)
  list(fit = fit, original_signal_index = orig_idx)
}

shared_qtl_variants <- function(a, b) {
  ast <- as.data.table(copy(a$stats))
  bst <- as.data.table(copy(b$stats))
  ast[, a_idx := .I]
  bst[, b_idx := .I]
  key_cols <- intersect(c("variant_key", "chr", "pos", "effect_allele", "other_allele"), names(ast))
  key_cols <- intersect(key_cols, names(bst))
  if ("variant_key" %in% key_cols) {
    m <- merge(ast[, .(variant_key, a_idx)], bst[, .(variant_key, b_idx)], by = "variant_key")
    if (nrow(m)) {
      m <- unique(m, by = c("a_idx", "b_idx", "variant_key"))
      return(m[, .(snp = variant_key, a_idx, b_idx)])
    }
  }
  need <- c("chr", "pos", "effect_allele", "other_allele")
  if (!all(need %in% key_cols)) return(data.table(snp = character(), a_idx = integer(), b_idx = integer()))
  m <- merge(
    ast[, .(chr, pos, effect_allele, other_allele, a_idx)],
    bst[, .(chr, pos, effect_allele, other_allele, b_idx)],
    by = c("chr", "pos", "effect_allele", "other_allele")
  )
  if (!nrow(m)) return(data.table(snp = character(), a_idx = integer(), b_idx = integer()))
  m[, snp := standard_variant_id(chr, pos, effect_allele, other_allele)]
  unique(m[, .(snp, a_idx, b_idx)], by = c("a_idx", "b_idx", "snp"))
}

status <- fread(file.path(root_dir, "results/fine_mapping/qtl_finemap_status.tsv"))
status <- status[celltype %chin% names(CELLTYPE_EQTL_MAP)]
status <- status[status == "ok" & n_cs > 0 & file.exists(susie_rds)]
pc <- status[phenotype_type == "pcQTL"]
eg <- status[phenotype_type == "eQTL"]
pairs <- merge(
  pc[, .(celltype, cluster_id, pc_id = phenotype_id, pc_susie_rds = susie_rds)],
  eg[, .(celltype, cluster_id, gene_id = phenotype_id, eqtl_susie_rds = susie_rds)],
  by = c("celltype", "cluster_id"),
  allow.cartesian = TRUE
)
setorder(pairs, celltype, cluster_id, pc_id, gene_id)

process_pair <- function(row) {
  base <- row[, .(celltype, cluster_id, pc_id, gene_id, pc_susie_rds, eqtl_susie_rds)]
  tryCatch({
    a <- readRDS(row$pc_susie_rds)
    b <- readRDS(row$eqtl_susie_rds)
    shared <- shared_qtl_variants(a, b)
    if (nrow(shared) < 2L) {
      return(cbind(base, data.table(
        status = "no_test",
        message = "fewer than two shared QTL variants",
        pc_signal_index = NA_integer_,
        eqtl_signal_index = NA_integer_,
        n_shared_variants = nrow(shared),
        n_coloc_variants = 0L,
        PP.H0.abf = NA_real_, PP.H1.abf = NA_real_, PP.H2.abf = NA_real_,
        PP.H3.abf = NA_real_, PP.H4.abf = NA_real_, pph4 = NA_real_,
        colocalized = FALSE
      )))
    }
    a_fit <- set_susie_variant_names(subset_susie_to_variants(a$fit, shared$a_idx), shared$snp)
    b_fit <- set_susie_variant_names(subset_susie_to_variants(b$fit, shared$b_idx), shared$snp)
    a_cmp <- compact_susie_cs(a_fit)
    b_cmp <- compact_susie_cs(b_fit)
    if (!length(a_cmp$original_signal_index) || !length(b_cmp$original_signal_index)) {
      return(cbind(base, data.table(
        status = "no_test",
        message = "no non-empty QTL CS after shared variant subsetting",
        pc_signal_index = NA_integer_,
        eqtl_signal_index = NA_integer_,
        n_shared_variants = nrow(shared),
        n_coloc_variants = 0L,
        PP.H0.abf = NA_real_, PP.H1.abf = NA_real_, PP.H2.abf = NA_real_,
        PP.H3.abf = NA_real_, PP.H4.abf = NA_real_, pph4 = NA_real_,
        colocalized = FALSE
      )))
    }
    cr <- coloc.susie(a_cmp$fit, b_cmp$fit, trim_by_posterior = FALSE)
    res <- if ("summary" %in% names(cr)) as.data.table(cr$summary) else as.data.table(cr)
    if (!nrow(res) || !"idx1" %in% names(res) || !"idx2" %in% names(res)) {
      return(cbind(base, data.table(
        status = "no_test",
        message = "coloc.susie returned no signal-pair summary",
        pc_signal_index = NA_integer_,
        eqtl_signal_index = NA_integer_,
        n_shared_variants = nrow(shared),
        n_coloc_variants = nrow(shared),
        PP.H0.abf = NA_real_, PP.H1.abf = NA_real_, PP.H2.abf = NA_real_,
        PP.H3.abf = NA_real_, PP.H4.abf = NA_real_, pph4 = NA_real_,
        colocalized = FALSE
      )))
    }
    pp_col <- grep("PP\\.H4|PPH4", names(res), value = TRUE)[1]
    if (is.na(pp_col)) stop("Cannot find PPH4 column in QTL-QTL coloc result")
    out <- data.table(
      status = "ok",
      message = "",
      pc_signal_index = a_cmp$original_signal_index[as.integer(res$idx1)],
      eqtl_signal_index = b_cmp$original_signal_index[as.integer(res$idx2)],
      n_shared_variants = nrow(shared),
      n_coloc_variants = as.integer(res$nsnps %||% nrow(shared)),
      PP.H0.abf = as.numeric(res$PP.H0.abf %||% NA_real_),
      PP.H1.abf = as.numeric(res$PP.H1.abf %||% NA_real_),
      PP.H2.abf = as.numeric(res$PP.H2.abf %||% NA_real_),
      PP.H3.abf = as.numeric(res$PP.H3.abf %||% NA_real_),
      PP.H4.abf = as.numeric(res[[pp_col]]),
      pph4 = as.numeric(res[[pp_col]])
    )
    out[, colocalized := is.finite(pph4) & pph4 > COLOC_PPH4_CUTOFF]
    cbind(base, out)
  }, error = function(e) {
    cbind(base, data.table(
      status = "failed",
      message = conditionMessage(e),
      pc_signal_index = NA_integer_,
      eqtl_signal_index = NA_integer_,
      n_shared_variants = NA_integer_,
      n_coloc_variants = NA_integer_,
      PP.H0.abf = NA_real_, PP.H1.abf = NA_real_, PP.H2.abf = NA_real_,
      PP.H3.abf = NA_real_, PP.H4.abf = NA_real_, pph4 = NA_real_,
      colocalized = FALSE
    ))
  })
}

rows <- vector("list", nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  rows[[i]] <- process_pair(pairs[i])
  if (i %% 100L == 0L) message(sprintf("Processed QTL-QTL pairs: %d / %d", i, nrow(pairs)))
}

edge_dt <- rbindlist(rows, fill = TRUE)
fwrite(edge_dt, file.path(out_dir, "eqtl_pcqtl_signal_coloc_edges.tsv"), sep = "\t")

summary_dt <- edge_dt[, .N, by = .(status)][order(-N)]
summary_dt <- rbind(
  summary_dt,
  data.table(status = "colocalized_pph4_gt_0.75", N = edge_dt[status == "ok" & pph4 > COLOC_PPH4_CUTOFF, .N]),
  fill = TRUE
)
fwrite(summary_dt, file.path(out_dir, "eqtl_pcqtl_signal_coloc_edges_qc.tsv"), sep = "\t")
print(summary_dt)
