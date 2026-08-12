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

endpoint_index <- as.integer(get_arg("--endpoint-index", Sys.getenv("SLURM_ARRAY_TASK_ID", NA_character_)))
window_kb <- as.integer(get_arg("--window-kb", "250"))
min_shared <- as.integer(get_arg("--min-shared", "50"))

if (is.na(endpoint_index)) stop("Provide --endpoint-index or SLURM_ARRAY_TASK_ID")

out_dir <- file.path(ROOT_DIR, "results/coloc/qtl_gwas_abf_raw", sprintf("window_%skb", window_kb))
chunk_dir <- file.path(out_dir, "chunks")
dir_create(chunk_dir)

endpoint_tasks <- fread(file.path(ROOT_DIR, "manifests/qtl_gwas_abf_raw_endpoint_tasks.tsv"))
endpoint <- endpoint_tasks[endpoint_task_index == endpoint_index]
if (nrow(endpoint) != 1L) stop("Expected exactly one endpoint row for index ", endpoint_index)
linked_gwas <- file.path("data/finngen_raw_links", basename(endpoint$local_gz))
gwas_path <- if (file.exists(linked_gwas)) linked_gwas else endpoint$resolved_local_gz

qtl_leads <- fread(file.path(ROOT_DIR, "manifests/qtl_gwas_abf_raw_qtl_leads.tsv"))
qtl_stats <- readRDS(file.path(ROOT_DIR, "results/coloc/qtl_gwas_abf_raw/cache/qtl_abf_stats.rds"))
setDT(qtl_stats)
qtl_leads <- qtl_leads[celltype %chin% names(CELLTYPE_EQTL_MAP)]
qtl_stats <- qtl_stats[celltype %chin% names(CELLTYPE_EQTL_MAP)]
if (!nrow(qtl_leads) || !nrow(qtl_stats)) stop("No primary-analysis QTL inputs remain")
setkey(qtl_stats, qtl_id, chr, pos)

is_ambiguous <- function(a, b) paste0(a, b) %in% c("AT", "TA", "CG", "GC")

read_finngen_raw <- function(path) {
  raw <- fread(
    cmd = paste("zcat", shQuote(path)),
    select = c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"),
    showProgress = FALSE
  )
  setnames(raw, names(raw), gsub("^#", "", names(raw)))
  stats <- raw[, .(
    chr = sub("^chr", "", as.character(chrom)),
    pos = as.integer(pos),
    marker_id = as.character(rsids),
    effect_allele = as.character(alt),
    other_allele = as.character(ref),
    gwas_beta = as.numeric(beta),
    gwas_se = as.numeric(sebeta),
    gwas_pvalue = as.numeric(pval),
    gwas_maf = pmin(as.numeric(af_alt), 1 - as.numeric(af_alt))
  )]
  rm(raw)
  gc()
  stats <- stats[is.finite(gwas_beta) & is.finite(gwas_se) & gwas_se > 0 &
                   is.finite(gwas_pvalue) & is.finite(gwas_maf) & gwas_maf > 0 & gwas_maf < 1]
  unique(stats, by = c("chr", "pos", "effect_allele", "other_allele"))
}

harmonize_region <- function(qtl, gwas) {
  q <- copy(qtl)
  g <- copy(gwas)
  setnames(q, c("effect_allele", "other_allele", "beta", "se", "pvalue"),
           c("q_effect", "q_other", "qtl_beta", "qtl_se", "qtl_pvalue"))
  setnames(g, c("effect_allele", "other_allele"),
           c("g_effect", "g_other"))
  m <- merge(
    q[, .(chr, pos, q_effect, q_other, qtl_beta, qtl_se, qtl_pvalue, qtl_maf, qtl_n)],
    g[, .(chr, pos, marker_id, g_effect, g_other, gwas_beta, gwas_se, gwas_pvalue, gwas_maf)],
    by = c("chr", "pos"),
    allow.cartesian = TRUE
  )
  if (!nrow(m)) return(list(dt = m, allele_flipped_n = 0L, ambiguous_removed_n = 0L, mismatch_removed_n = 0L, duplicated_removed_n = 0L))
  same <- m$q_effect == m$g_effect & m$q_other == m$g_other
  flip <- m$q_effect == m$g_other & m$q_other == m$g_effect
  mismatch_removed_n <- sum(!(same | flip))
  m <- m[same | flip]
  allele_flipped_n <- sum(flip[same | flip])
  if (nrow(m)) m[flip[same | flip], gwas_beta := -gwas_beta]
  ambiguous <- is_ambiguous(m$q_effect, m$q_other)
  ambiguous_removed_n <- sum(ambiguous, na.rm = TRUE)
  m <- m[!ambiguous]
  if (!nrow(m)) return(list(dt = m, allele_flipped_n = allele_flipped_n, ambiguous_removed_n = ambiguous_removed_n, mismatch_removed_n = mismatch_removed_n, duplicated_removed_n = 0L))
  m[, snp := standard_variant_id(chr, pos, q_effect, q_other)]
  dup <- duplicated(m$snp) | duplicated(m$snp, fromLast = TRUE)
  duplicated_removed_n <- sum(dup)
  m <- m[!dup]
  list(dt = m, allele_flipped_n = allele_flipped_n, ambiguous_removed_n = ambiguous_removed_n,
       mismatch_removed_n = mismatch_removed_n, duplicated_removed_n = duplicated_removed_n)
}

make_status_row <- function(row, status, message, region_start, region_end, n_shared = NA_integer_,
                            allele_flipped_n = NA_integer_, ambiguous_removed_n = NA_integer_,
                            mismatch_removed_n = NA_integer_, duplicated_removed_n = NA_integer_) {
  data.table(
    phenocode = endpoint$phenocode,
    phenotype = endpoint$phenotype,
    category = endpoint$category,
    celltype = row$celltype,
    cluster_id = row$cluster_id,
    qtl_type = row$qtl_type,
    phenotype_id = row$phenotype_id,
    qtl_id = row$qtl_id,
    region_chr = row$lead_chr,
    region_start = region_start,
    region_end = region_end,
    window_kb = window_kb,
    n_shared_snps = n_shared,
    qtl_lead_variant = row$qtl_lead_variant,
    qtl_lead_p = row$qtl_lead_p,
    gwas_min_p_variant = NA_character_,
    qtl_min_p = NA_real_,
    gwas_min_p = NA_real_,
    qtl_N = row$qtl_n,
    gwas_N = as.integer(endpoint$num_cases + endpoint$num_controls),
    n_cases = as.integer(endpoint$num_cases),
    n_controls = as.integer(endpoint$num_controls),
    allele_flipped_n = allele_flipped_n,
    ambiguous_removed_n = ambiguous_removed_n,
    mismatch_removed_n = mismatch_removed_n,
    duplicated_removed_n = duplicated_removed_n,
    PP.H0.abf = NA_real_,
    PP.H1.abf = NA_real_,
    PP.H2.abf = NA_real_,
    PP.H3.abf = NA_real_,
    PP.H4.abf = NA_real_,
    PPH4_gt_PPH3 = FALSE,
    coloc_pass = FALSE,
    coloc_pass_h4_0_7 = FALSE,
    coloc_pass_h4_0_8 = FALSE,
    status = status,
    message = message
  )
}

rows <- list()
if (!file.exists(gwas_path)) {
  for (i in seq_len(nrow(qtl_leads))) {
    r <- qtl_leads[i]
    rows[[i]] <- make_status_row(r, "missing_gwas", "FinnGen raw GWAS file is missing",
                                 r$lead_pos - window_kb * 1000L, r$lead_pos + window_kb * 1000L)
  }
} else {
  message("Reading FinnGen raw GWAS: ", endpoint$phenocode)
  gwas_stats <- read_finngen_raw(gwas_path)
  setkey(gwas_stats, chr, pos)
  row_i <- 0L
  hla_chr <- "6"
  hla_start <- 25000000L
  hla_end <- 34000000L
  for (i in seq_len(nrow(qtl_leads))) {
    qrow <- qtl_leads[i]
    start <- max(1L, as.integer(qrow$lead_pos - window_kb * 1000L))
    end <- as.integer(qrow$lead_pos + window_kb * 1000L)
    if (as.character(qrow$lead_chr) == hla_chr && end >= hla_start && start <= hla_end) {
      row_i <- row_i + 1L
      rows[[row_i]] <- make_status_row(qrow, "hla_excluded", "Region overlaps HLA", start, end)
      next
    }
    qtl_chr <- qtl_stats[J(qrow$qtl_id, as.character(qrow$lead_chr)), nomatch = 0]
    qtl_region <- qtl_chr[pos >= start & pos <= end]
    gwas_chr <- gwas_stats[J(as.character(qrow$lead_chr)), nomatch = 0]
    gwas_region <- gwas_chr[pos >= start & pos <= end]
    hm <- harmonize_region(qtl_region, gwas_region)
    shared <- hm$dt
    if (!nrow(shared)) {
      row_i <- row_i + 1L
      rows[[row_i]] <- make_status_row(qrow, "no_shared_snps", "No allele-compatible shared SNPs after filtering",
                                       start, end, 0L, hm$allele_flipped_n, hm$ambiguous_removed_n,
                                       hm$mismatch_removed_n, hm$duplicated_removed_n)
      next
    }
    if (nrow(shared) < min_shared) {
      row_i <- row_i + 1L
      rows[[row_i]] <- make_status_row(qrow, "low_shared_snps", sprintf("n_shared_snps < %d", min_shared),
                                       start, end, nrow(shared), hm$allele_flipped_n, hm$ambiguous_removed_n,
                                       hm$mismatch_removed_n, hm$duplicated_removed_n)
      next
    }
    shared <- shared[is.finite(qtl_beta) & is.finite(qtl_se) & qtl_se > 0 &
                       is.finite(gwas_beta) & is.finite(gwas_se) & gwas_se > 0 &
                       is.finite(qtl_maf) & qtl_maf > 0 & qtl_maf < 1 &
                       is.finite(gwas_maf) & gwas_maf > 0 & gwas_maf < 1]
    if (nrow(shared) < min_shared) {
      row_i <- row_i + 1L
      rows[[row_i]] <- make_status_row(qrow, "low_shared_snps", sprintf("n_shared_snps < %d after finite-value filtering", min_shared),
                                       start, end, nrow(shared), hm$allele_flipped_n, hm$ambiguous_removed_n,
                                       hm$mismatch_removed_n, hm$duplicated_removed_n)
      next
    }
    qtl_dataset <- list(
      beta = shared$qtl_beta,
      varbeta = shared$qtl_se^2,
      snp = shared$snp,
      position = shared$pos,
      MAF = shared$qtl_maf,
      N = suppressWarnings(max(shared$qtl_n, na.rm = TRUE)),
      type = "quant"
    )
    gwas_n <- as.integer(endpoint$num_cases + endpoint$num_controls)
    gwas_dataset <- list(
      beta = shared$gwas_beta,
      varbeta = shared$gwas_se^2,
      snp = shared$snp,
      position = shared$pos,
      MAF = shared$gwas_maf,
      N = gwas_n,
      type = "cc",
      s = as.numeric(endpoint$num_cases) / gwas_n
    )
    status <- "ok"
    msg <- ""
    pp <- rep(NA_real_, 5)
    names(pp) <- paste0("PP.H", 0:4, ".abf")
    tryCatch({
      invisible(capture.output({
        res <- suppressMessages(coloc::coloc.abf(qtl_dataset, gwas_dataset))
      }))
      s <- as.list(res$summary)
      pp[] <- as.numeric(unlist(s[names(pp)]))
    }, error = function(e) {
      status <<- "failed"
      msg <<- conditionMessage(e)
    })
    qmin <- shared[which.min(qtl_pvalue)]
    gmin <- shared[which.min(gwas_pvalue)]
    row_i <- row_i + 1L
    out <- make_status_row(qrow, status, msg, start, end, nrow(shared), hm$allele_flipped_n,
                           hm$ambiguous_removed_n, hm$mismatch_removed_n, hm$duplicated_removed_n)
    out[, `:=`(
      gwas_min_p_variant = standard_variant_id(gmin$chr, gmin$pos, gmin$q_effect, gmin$q_other),
      qtl_min_p = as.numeric(qmin$qtl_pvalue),
      gwas_min_p = as.numeric(gmin$gwas_pvalue),
      PP.H0.abf = pp["PP.H0.abf"],
      PP.H1.abf = pp["PP.H1.abf"],
      PP.H2.abf = pp["PP.H2.abf"],
      PP.H3.abf = pp["PP.H3.abf"],
      PP.H4.abf = pp["PP.H4.abf"]
    )]
    out[, `:=`(
      PPH4_gt_PPH3 = !is.na(PP.H4.abf) & !is.na(PP.H3.abf) & PP.H4.abf > PP.H3.abf,
      coloc_pass = status == "ok" & !is.na(PP.H4.abf) & PP.H4.abf >= 0.75 & PP.H4.abf > PP.H3.abf,
      coloc_pass_h4_0_7 = status == "ok" & !is.na(PP.H4.abf) & PP.H4.abf >= 0.7 & PP.H4.abf > PP.H3.abf,
      coloc_pass_h4_0_8 = status == "ok" & !is.na(PP.H4.abf) & PP.H4.abf >= 0.8 & PP.H4.abf > PP.H3.abf
    )]
    rows[[row_i]] <- out
  }
}

dt <- rbindlist(rows, fill = TRUE)
out_file <- file.path(chunk_dir, sprintf("qtl_gwas_abf_raw_%s_window_%skb.tsv", endpoint$phenocode, window_kb))
fwrite(dt, out_file, sep = "\t", quote = FALSE)
message("Wrote ", nrow(dt), " ABF raw rows to ", out_file)
