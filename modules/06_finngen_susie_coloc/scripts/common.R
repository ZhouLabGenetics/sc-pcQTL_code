suppressPackageStartupMessages({
  library(data.table)
})

load_config <- function() {
  config_path <- file.path(dirname(dirname(normalizePath(sys.frame(1)$ofile %||% "."))), "config", "config.R")
  if (!file.exists(config_path)) {
    config_path <- file.path(getwd(), "config", "config.R")
  }
  source(config_path, local = FALSE)
  invisible(TRUE)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x

messagef <- function(...) message(sprintf(...))

read_region_file <- function(path) {
  x <- fread(path, header = FALSE, col.names = c("chr", "start", "end"))
  x[, chr := sub("^chr", "", as.character(chr))]
  x
}

standard_variant_id <- function(chr, pos, a1 = NULL, a2 = NULL) {
  chr <- sub("^chr", "", as.character(chr))
  pos <- as.integer(pos)
  if (is.null(a1) || is.null(a2)) {
    paste(chr, pos, sep = ":")
  } else {
    paste0(chr, ":", pos, "_", a1, "/", a2)
  }
}

normalize_sumstats <- function(dt, source = c("qtl", "gwas")) {
  source <- match.arg(source)
  nm <- names(dt)
  setnames(dt, nm, gsub("^#", "", nm))

  if (source == "qtl") {
    required <- c("CHR", "POS", "MarkerID", "Allele1", "Allele2", "BETA", "SE", "p.value", "N")
    miss <- setdiff(required, names(dt))
    if (length(miss)) stop("Missing QTL columns: ", paste(miss, collapse = ", "))
    out <- dt[, .(
      chr = sub("^chr", "", as.character(CHR)),
      pos = as.integer(POS),
      marker_id = as.character(MarkerID),
      effect_allele = as.character(Allele2),
      other_allele = as.character(Allele1),
      beta = as.numeric(BETA),
      se = as.numeric(SE),
      pvalue = as.numeric(p.value),
      af = as.numeric(AF_Allele2),
      n = as.integer(N)
    )]
  } else {
    required <- c("chrom", "pos", "ref", "alt", "pval", "beta", "sebeta")
    miss <- setdiff(required, names(dt))
    if (length(miss)) stop("Missing FinnGen columns: ", paste(miss, collapse = ", "))
    out <- dt[, .(
      chr = sub("^chr", "", as.character(chrom)),
      pos = as.integer(pos),
      marker_id = if ("rsids" %in% names(dt)) as.character(rsids) else paste(chrom, pos, sep = ":"),
      effect_allele = as.character(alt),
      other_allele = as.character(ref),
      beta = as.numeric(beta),
      se = as.numeric(sebeta),
      pvalue = as.numeric(pval),
      af = as.numeric(af_alt),
      n = NA_integer_
    )]
  }
  out[, variant_key := standard_variant_id(chr, pos)]
  out <- out[is.finite(beta) & is.finite(se) & se > 0 & is.finite(pvalue)]
  unique(out, by = c("chr", "pos", "effect_allele", "other_allele"))
}

harmonize_to_reference <- function(stats, ref_variants) {
  required <- c("chr", "pos", "effect_allele", "other_allele")
  miss <- setdiff(required, names(ref_variants))
  if (length(miss)) stop("Missing reference columns: ", paste(miss, collapse = ", "))
  ref <- copy(ref_variants)
  stats <- copy(stats)
  stats[, stat_row := .I]
  ref[, `:=`(
    ref_row = .I,
    ref_effect = effect_allele,
    ref_other = other_allele,
    ref_key = standard_variant_id(chr, pos, effect_allele, other_allele)
  )]
  m <- merge(stats, ref[, .(chr, pos, ref_row, ref_effect, ref_other, ref_key)],
             by = c("chr", "pos"), allow.cartesian = TRUE)
  same <- m$effect_allele == m$ref_effect & m$other_allele == m$ref_other
  flip <- m$effect_allele == m$ref_other & m$other_allele == m$ref_effect
  ambig <- paste0(m$effect_allele, m$other_allele) %in% c("AT", "TA", "CG", "GC")
  m[, matched_orientation := same]
  m[, flipped_orientation := flip]
  m[, ambiguous_snp := ambig]
  m <- m[matched_orientation | flipped_orientation]
  m[flipped_orientation == TRUE, beta := -beta]
  m[, `:=`(effect_allele = ref_effect, other_allele = ref_other, variant_key = ref_key)]
  m[, allele_status := fifelse(ambiguous_snp, "ambiguous_retained_after_exact_match", fifelse(flipped_orientation, "flipped", "matched"))]
  setorder(m, stat_row, ref_row)
  m <- unique(m, by = c("stat_row", "ref_row", "variant_key"))
  m[, c("stat_row", "ref_row", "ref_effect", "ref_other", "ref_key", "matched_orientation", "flipped_orientation", "ambiguous_snp") := NULL]
  m[]
}

read_ld_matrix <- function(path) {
  if (grepl("\\.gz$", path)) {
    as.matrix(fread(cmd = paste("zcat", shQuote(path)), header = FALSE))
  } else {
    as.matrix(fread(path, header = FALSE))
  }
}

run_susie_from_sumstats <- function(stats, ld, L, coverage, min_abs_corr) {
  if (!requireNamespace("susieR", quietly = TRUE)) stop("susieR is not installed")
  z <- stats$beta / stats$se
  keep <- is.finite(z) & complete.cases(ld)
  stats <- stats[keep]
  ld <- ld[keep, keep, drop = FALSE]
  if (nrow(stats) < 2L) stop("Need at least two variants after filtering")
  sample_n <- unique(stats$n[is.finite(stats$n) & stats$n > 0])
  if (length(sample_n) != 1L) {
    stop("Expected exactly one positive summary-statistic sample size per phenotype")
  }
  fit <- susieR::susie_rss(z = z, R = ld, n = sample_n,
                           L = L, coverage = coverage, estimate_residual_variance = FALSE)
  cs <- susieR::susie_get_cs(fit, Xcorr = ld, min_abs_corr = min_abs_corr, coverage = coverage)
  list(fit = fit, cs = cs, stats = stats)
}

subset_susie_to_variants <- function(susie_fit, keep_idx) {
  old_to_new <- rep(NA_integer_, ncol(susie_fit$lbf_variable))
  old_to_new[keep_idx] <- seq_along(keep_idx)
  fit <- susie_fit
  if (!is.null(fit$lbf_variable)) fit$lbf_variable <- fit$lbf_variable[, keep_idx, drop = FALSE]
  if (!is.null(fit$alpha)) fit$alpha <- fit$alpha[, keep_idx, drop = FALSE]
  if (!is.null(fit$mu)) fit$mu <- fit$mu[, keep_idx, drop = FALSE]
  if (!is.null(fit$mu2)) fit$mu2 <- fit$mu2[, keep_idx, drop = FALSE]
  if (!is.null(fit$Xr)) fit$Xr <- fit$Xr[keep_idx]
  if (!is.null(fit$pip)) fit$pip <- fit$pip[keep_idx]
  if (!is.null(fit$sets$cs) && length(fit$sets$cs)) {
    fit$sets$cs <- lapply(fit$sets$cs, function(x) {
      y <- old_to_new[x]
      y[!is.na(y)]
    })
    keep_cs <- lengths(fit$sets$cs) > 0L
    fit$sets$cs <- fit$sets$cs[keep_cs]
    if (!is.null(fit$sets$cs_index)) fit$sets$cs_index <- fit$sets$cs_index[keep_cs]
    if (!is.null(fit$sets$coverage)) fit$sets$coverage <- fit$sets$coverage[keep_cs]
    if (!is.null(fit$sets$purity)) fit$sets$purity <- fit$sets$purity[keep_cs, , drop = FALSE]
  }
  fit
}

extract_susie_cs <- function(fm, phenotype_id, phenotype_type, context) {
  stats <- fm$stats
  fit <- fm$fit
  cs <- fm$cs$cs
  if (is.null(cs) || length(cs) == 0L) {
    return(data.table())
  }
  rbindlist(lapply(seq_along(cs), function(i) {
    idx <- cs[[i]]
    d <- stats[idx]
    d[, `:=`(
      phenotype_id = phenotype_id,
      phenotype_type = phenotype_type,
      credible_set_id = paste(context, phenotype_type, phenotype_id, paste0("CS", i), sep = "__"),
      cs_index = i,
      pip = fit$pip[idx],
      cs_size = length(idx),
      lead_pip = max(fit$pip[idx], na.rm = TRUE)
    )]
    lead <- d[which.max(pip)]
    d[, `:=`(
      lead_variant = standard_variant_id(lead$chr, lead$pos, lead$effect_allele, lead$other_allele),
      lead_chr = lead$chr,
      lead_pos = lead$pos
    )]
    d
  }), fill = TRUE)
}
