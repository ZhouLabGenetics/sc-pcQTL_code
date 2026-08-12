#!/usr/bin/env Rscript

.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

.this_dir <- local({
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
})
root <- Sys.getenv("SC_PCQTL_FORMAL_COLOC_ROOT", unset = normalizePath(file.path(getwd(), "..", "06_finngen_susie_coloc"), mustWork = FALSE))
work_dir <- Sys.getenv(
  "SC_PCQTL_LOCUSZOOM_WORK_ROOT",
  unset = file.path(root, "10_publication_locuszoom_redesign")
)
source(file.path(dirname(.this_dir), "06_finngen_susie_coloc/config/config.R"))

data_dir <- file.path(work_dir, "data", "gimap_finngen_style")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

example <- list(
  celltype = "cd8_nc",
  cluster_id = "SC_chr7_cluster_001",
  qtl_type = "pcQTL",
  qtl_phenotype_id = "PC3",
  phenocode = "3019198"
)

qtl_rds <- file.path(root, "results/fine_mapping/qtl/cd8_nc/SC_chr7_cluster_001/pcQTL/PC3.susie.rds")
gwas_rds <- file.path(root, "results/fine_mapping/gwas_all_finemapped/rds/3019198/cd8_nc/SC_chr7_cluster_001/gwas.official_finngen_susie.rds")
coloc_file <- file.path(root, "results/coloc/qtl_gwas_susie_official_finngen_all_finemapped/qtl_gwas_susie_official_finngen_coloc_summary.tsv")
locuszoom_data_root <- Sys.getenv("SC_PCQTL_LOCUSZOOM_DATA_ROOT", unset = "")
gene_json <- file.path(locuszoom_data_root, "candidate__adjusted_only__adjusted__cd8_nc__SC_chr7_cluster_001__PC3_genes.json")
meta_json <- file.path(locuszoom_data_root, "candidate__adjusted_only__adjusted__cd8_nc__SC_chr7_cluster_001__PC3_meta.json")
effect_file <- file.path(root, "09_mechanistic_celltype_analysis/results/gene_effects/strict_pip_weighted_nominal_gene_effects.tsv")

stopifnot(file.exists(qtl_rds), file.exists(gwas_rds), file.exists(coloc_file))

standard_variant_id <- function(chr, pos, ea = NULL, oa = NULL) {
  chr <- sub("^chr", "", as.character(chr))
  if (is.null(ea) || is.null(oa)) return(paste(chr, as.integer(pos), sep = ":"))
  paste0(chr, ":", as.integer(pos), "_", toupper(as.character(ea)), "/", toupper(as.character(oa)))
}

display_ld_key <- function(pos, ea, oa) {
  ea <- toupper(as.character(ea))
  oa <- toupper(as.character(oa))
  allele_pair <- ifelse(ea <= oa, paste0(ea, "/", oa), paste0(oa, "/", ea))
  paste(as.integer(pos), allele_pair, sep = ":")
}

locuszoom_ld_levels <- c(
  "LD reference", "0.0-0.2", "0.2-0.4", "0.4-0.6",
  "0.6-0.8", "0.8-1.0", "No LD"
)
locuszoom_ld_group <- function(r2, is_ref) {
  r2 <- suppressWarnings(as.numeric(r2))
  out <- rep("No LD", length(r2))
  ok <- is.finite(r2)
  out[ok] <- as.character(cut(
    pmax(pmin(r2[ok], 1), 0),
    breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
    labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"),
    right = FALSE
  ))
  out[is_ref %in% TRUE] <- "LD reference"
  factor(out, levels = locuszoom_ld_levels)
}

apply_unified_display_ld <- function(gstats, qstats) {
  qref <- unique(
    qstats[is.finite(lead_r2), .(
      display_ld_key = display_ld_key(pos, effect_allele, other_allele),
      unified_lead_r2 = pmax(pmin(as.numeric(lead_r2), 1), 0)
    )],
    by = "display_ld_key"
  )
  gstats[, display_ld_key := display_ld_key(pos, effect_allele, other_allele)]
  gstats <- merge(gstats, qref, by = "display_ld_key", all.x = TRUE, sort = FALSE)
  gstats[, `:=`(
    lead_r2 = fifelse(is.finite(unified_lead_r2), unified_lead_r2, NA_real_),
    ld_source_for_plot = fifelse(
      is.finite(unified_lead_r2),
      "intermediate QTL LD before 1000G display override",
      "missing intermediate QTL LD before 1000G display override"
    )
  )]
  gstats[, c("display_ld_key", "unified_lead_r2") := NULL]
  gstats[]
}

neglog10p <- function(p) {
  p <- as.numeric(p)
  p[!is.finite(p) | p <= 0] <- .Machine$double.xmin
  -log10(p)
}

read_liftover_posmap <- function(chr) {
  path <- sprintf(ONEK1K_QTL_LIFTOVER_POSMAP_PATTERN, sub("^chr", "", as.character(chr)))
  if (!file.exists(path)) stop("Missing liftOver position map: ", path)
  dt <- fread(path, header = FALSE, col.names = c("old_id", "pos_hg38"))
  dt[, `:=`(
    chr = sub("^chr", "", sub(":.*$", "", old_id)),
    pos_hg19 = as.integer(sub("^.*:", "", old_id)),
    pos_hg38 = as.integer(pos_hg38)
  )]
  unique(dt[is.finite(pos_hg19) & is.finite(pos_hg38), .(chr, pos_hg19, pos_hg38)],
         by = c("chr", "pos_hg19"))
}

liftover_qtl_positions <- function(dt, row_col = "q_idx") {
  liftover_bin <- file.path(root, "resources/liftover_tools/liftOver")
  chain_file <- file.path(root, "resources/liftover_tools/hg19ToHg38.over.chain.gz")
  out <- data.table(
    row_id = dt[[row_col]],
    chr = sub("^chr", "", as.character(dt$chr)),
    pos_hg19 = as.integer(dt$pos_hg19),
    pos_hg38 = NA_integer_,
    liftover_source = NA_character_
  )
  if (file.exists(liftover_bin) && file.exists(chain_file)) {
    tmp_dir <- tempfile("qtl_liftover_")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    bed_in <- file.path(tmp_dir, "qtl.hg19.bed")
    bed_out <- file.path(tmp_dir, "qtl.hg38.bed")
    bed_unmapped <- file.path(tmp_dir, "qtl.unmapped.bed")
    bed <- out[is.finite(pos_hg19), .(
      chrom = paste0("chr", chr),
      start = pmax(0L, pos_hg19 - 1L),
      end = pos_hg19,
      name = as.character(row_id)
    )]
    fwrite(bed, bed_in, sep = "\t", col.names = FALSE)
    cmd_status <- system2(liftover_bin, c(bed_in, chain_file, bed_out, bed_unmapped),
                          stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(cmd_status, "status")) && attr(cmd_status, "status") != 0) {
      stop("UCSC liftOver failed: ", paste(cmd_status, collapse = "\n"))
    }
    if (file.exists(bed_out) && file.info(bed_out)$size > 0) {
      mapped <- fread(
        bed_out,
        header = FALSE,
        col.names = c("chrom", "start", "end", "name")
      )
      mapped[, `:=`(
        row_id = suppressWarnings(as.integer(name)),
        pos_hg38 = as.integer(end)
      )]
      mapped <- unique(mapped[is.finite(row_id) & is.finite(pos_hg38), .(row_id, pos_hg38)], by = "row_id")
      out[mapped, on = "row_id", `:=`(pos_hg38 = i.pos_hg38, liftover_source = "ucsc_hg19ToHg38_chain")]
    }
  } else {
    posmap <- read_liftover_posmap(unique(out$chr))
    out <- merge(out, posmap, by = c("chr", "pos_hg19"), all.x = TRUE, suffixes = c("", "_posmap"))
    out[is.na(pos_hg38) & is.finite(pos_hg38_posmap), `:=`(
      pos_hg38 = pos_hg38_posmap,
      liftover_source = "precomputed_posmap"
    )]
    out[, pos_hg38_posmap := NULL]
  }
  if (any(!is.finite(out$pos_hg38))) {
    warning(sprintf("Dropped %d QTL variants without hg38 liftOver coordinates.", sum(!is.finite(out$pos_hg38))))
  }
  setnames(out, "row_id", row_col)
  out[is.finite(pos_hg38)]
}

source(file.path(.this_dir, "eur_ld_display_helpers.R"))

get_cs_indices <- function(fit, signal_index) {
  cs_idx <- suppressWarnings(as.integer(fit$sets$cs_index))
  cs_pos <- match(as.integer(signal_index), cs_idx)
  if (is.na(cs_pos)) {
    if (length(fit$sets$cs) < signal_index) return(integer())
    return(as.integer(fit$sets$cs[[signal_index]]))
  }
  as.integer(fit$sets$cs[[cs_pos]])
}

read_finngen_region_snp <- function(snp_file, chr, start, end) {
  if (!file.exists(snp_file)) return(data.table())
  cols <- c(
    "region", "v", "rsid", "chromosome", "position", "allele1", "allele2",
    "p", "prob", "cs", "lead_r2", paste0("alpha", 1:10)
  )
  dt <- fread(cmd = paste("zcat", shQuote(snp_file)), select = cols, showProgress = FALSE)
  dt[chromosome == paste0("chr", sub("^chr", "", chr)) & position >= start & position <= end]
}

read_ld_matrix <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (grepl("\\.gz$", path)) {
    as.matrix(fread(cmd = paste("zcat", shQuote(path)), header = FALSE, showProgress = FALSE))
  } else {
    as.matrix(fread(path, header = FALSE, showProgress = FALSE))
  }
}

compact_gene_lanes <- function(genes, min_gap = 60000L) {
  genes <- copy(genes)
  setorder(genes, xmin, xmax)
  lane_end <- numeric()
  genes[, lane := 1L]
  for (i in seq_len(nrow(genes))) {
    placed <- FALSE
    if (length(lane_end)) {
      for (j in seq_along(lane_end)) {
        if (genes$xmin[i] > lane_end[j] + min_gap) {
          genes$lane[i] <- j
          lane_end[j] <- genes$xmax[i]
          placed <- TRUE
          break
        }
      }
    }
    if (!placed) {
      lane_end <- c(lane_end, genes$xmax[i])
      genes$lane[i] <- length(lane_end)
    }
  }
  genes[]
}

q <- readRDS(qtl_rds)
g <- readRDS(gwas_rds)

coloc <- fread(coloc_file)
coloc_row <- coloc[
  as.character(phenocode) == example$phenocode &
    celltype == example$celltype &
    cluster_id == example$cluster_id &
    qtl_type == example$qtl_type &
    qtl_phenotype_id == example$qtl_phenotype_id &
    status == "ok"
][which.max(pph4)]
if (nrow(coloc_row) != 1) stop("Could not identify one formal coloc row for the example.")

q_signal <- as.integer(coloc_row$qtl_signal_index)
g_signal <- as.integer(coloc_row$gwas_signal_index)

qstats <- as.data.table(copy(q$stats))
qstats[, `:=`(
  q_idx = .I,
  chr = sub("^chr", "", as.character(chr)),
  pos_hg19 = as.integer(pos),
  effect_allele = toupper(as.character(effect_allele)),
  other_allele = toupper(as.character(other_allele))
)]
q_lift <- liftover_qtl_positions(qstats[, .(q_idx, chr, pos_hg19)])
qstats <- merge(qstats, q_lift[, .(q_idx, pos_hg38, liftover_source)], by = "q_idx", all.x = FALSE, all.y = FALSE)
qstats[, `:=`(
  pos = as.integer(pos_hg38),
  mlog10p = neglog10p(pvalue),
  pip = as.numeric(q$fit$pip[q_idx]),
  alpha = as.numeric(q$fit$alpha[q_signal, q_idx]),
  variant_label = standard_variant_id(chr, pos_hg38, effect_allele, other_allele)
)]

gstats <- as.data.table(copy(g$stats))
gstats[, `:=`(
  g_idx = .I,
  chr = sub("^chr", "", as.character(chr)),
  pos = as.integer(pos),
  effect_allele = toupper(as.character(effect_allele)),
  other_allele = toupper(as.character(other_allele)),
  mlog10p = neglog10p(pvalue),
  pip = as.numeric(g$fit$pip),
  alpha = as.numeric(g$fit$alpha[g_signal, ]),
  variant_label = standard_variant_id(chr, pos, effect_allele, other_allele)
)]

q_cs <- get_cs_indices(q$fit, q_signal)
g_cs <- get_cs_indices(g$fit, g_signal)
qstats[, in_cs := q_idx %in% q_cs]
gstats[, in_cs := g_idx %in% g_cs]

q_lead <- qstats[which.max(alpha)]
g_lead <- gstats[which.max(alpha)]

shared <- merge(
  qstats[, .(chr, pos, q_idx, q_effect = effect_allele, q_other = other_allele)],
  gstats[, .(chr, pos, g_idx, g_effect = effect_allele, g_other = other_allele)],
  by = c("chr", "pos"),
  allow.cartesian = TRUE
)
shared[, allele_matched := q_effect == g_effect & q_other == g_other]
shared[, allele_flipped := q_effect == g_other & q_other == g_effect]
shared <- unique(shared[allele_matched | allele_flipped], by = c("q_idx", "g_idx"))
qstats[, shared := q_idx %in% shared$q_idx]
gstats[, shared := g_idx %in% shared$g_idx]

ld <- read_ld_matrix(q$metadata$ld_matrix)
if (!is.null(ld) && q_lead$ld_order <= nrow(ld)) {
  qstats[, lead_r2 := as.numeric(ld[as.integer(q_lead$ld_order), as.integer(ld_order)]^2)]
} else {
  qstats[, lead_r2 := NA_real_]
}
qstats[q_idx == q_lead$q_idx, lead_r2 := 1]
qstats[, ld_source_for_plot := fifelse(
  is.finite(lead_r2),
  "intermediate QTL LD before 1000G display override",
  "missing intermediate QTL LD before 1000G display override"
)]

fg_region <- read_finngen_region_snp(g$metadata$snp_file, "7", min(gstats$pos), max(gstats$pos))
if (nrow(fg_region)) {
  fg_region[, `:=`(
    chr = sub("^chr", "", chromosome),
    effect_allele = toupper(allele1),
    other_allele = toupper(allele2),
    fg_variant_label = standard_variant_id(chromosome, position, allele1, allele2)
  )]
  fg_region <- fg_region[, .(
    chr,
    pos = as.integer(position),
    fg_effect_allele = effect_allele,
    fg_other_allele = other_allele,
    fg_prob = as.numeric(prob),
    fg_cs = as.integer(cs),
    fg_lead_r2 = as.numeric(lead_r2),
    fg_alpha = as.numeric(get(paste0("alpha", g_signal)))
  )]
  fg_join <- merge(
    gstats[, .(g_idx, chr, pos, effect_allele, other_allele)],
    fg_region,
    by = c("chr", "pos"),
    allow.cartesian = TRUE
  )
  fg_join[, allele_matched := effect_allele == fg_effect_allele & other_allele == fg_other_allele]
  fg_join[, allele_flipped := effect_allele == fg_other_allele & other_allele == fg_effect_allele]
  fg_join <- fg_join[allele_matched | allele_flipped]
  fg_join[, match_rank := fifelse(!is.na(fg_lead_r2), 0L, fifelse(!is.na(fg_prob), 1L, 2L))]
  setorder(fg_join, g_idx, match_rank)
  fg_join <- unique(fg_join, by = "g_idx")
  gstats <- merge(
    gstats,
    fg_join[, .(g_idx, fg_prob, fg_cs, fg_lead_r2, fg_alpha, fg_allele_flipped = allele_flipped)],
    by = "g_idx",
    all.x = TRUE,
    sort = FALSE
  )
  gstats[, lead_r2 := fg_lead_r2]
  gstats[is.na(lead_r2) & g_idx == g_lead$g_idx, lead_r2 := 1]
} else {
  gstats[, lead_r2 := NA_real_]
}
gstats[, ld_source_for_plot := fifelse(!is.na(fg_lead_r2), "FinnGen official lead_r2", NA_character_)]
shared_ld <- merge(
  shared[, .(q_idx, g_idx)],
  qstats[, .(q_idx, qtl_lead_r2_for_shared = lead_r2)],
  by = "q_idx",
  all.x = TRUE,
  sort = FALSE
)
shared_ld <- unique(shared_ld[is.finite(qtl_lead_r2_for_shared), .(g_idx, qtl_lead_r2_for_shared)], by = "g_idx")
if (nrow(shared_ld)) {
  gstats <- merge(gstats, shared_ld, by = "g_idx", all.x = TRUE, sort = FALSE)
  gstats[
    is.na(lead_r2) & shared == TRUE & is.finite(qtl_lead_r2_for_shared),
    `:=`(
      lead_r2 = qtl_lead_r2_for_shared,
      ld_source_for_plot = "OneK1K LD fallback for shared SNP"
    )
  ]
} else {
  gstats[, qtl_lead_r2_for_shared := NA_real_]
}
gstats <- apply_unified_display_ld(gstats, qstats)
eur_ld <- add_1000g_eur_display_ld(
  qstats = qstats,
  gstats = gstats,
  q_lead = q_lead,
  g_lead = g_lead,
  chr = "7",
  formal_root = root,
  out_dir = file.path(data_dir, "eur_1000g_display_ld")
)
qstats <- eur_ld$qstats
gstats <- eur_ld$gstats
message(sprintf(
  "Using 1000G EUR Phase3 display LD anchor %s (%s reference SNPs in window).",
  ifelse(is.na(eur_ld$anchor_snp), "NA", eur_ld$anchor_snp),
  format(eur_ld$ref_n, big.mark = ",")
))

plot_start <- min(qstats$pos, na.rm = TRUE)
plot_end <- max(qstats$pos, na.rm = TRUE)

qplot_dt <- qstats[pos >= plot_start & pos <= plot_end]
gplot_dt <- gstats[pos >= plot_start & pos <= plot_end]
if (nrow(qplot_dt) == 0 || nrow(gplot_dt) == 0) stop("Empty plot range after hg38 alignment.")

gene_dt <- as.data.table(fromJSON(gene_json))
gene_dt[, gene_row := .I]
gene_start_lift <- liftover_qtl_positions(
  gene_dt[, .(gene_row, chr = chrom, pos_hg19 = as.integer(start))],
  row_col = "gene_row"
)
setnames(gene_start_lift, "pos_hg38", "start_hg38")
gene_end_lift <- liftover_qtl_positions(
  gene_dt[, .(gene_row, chr = chrom, pos_hg19 = as.integer(end))],
  row_col = "gene_row"
)
setnames(gene_end_lift, "pos_hg38", "end_hg38")
gene_dt <- merge(gene_dt, gene_start_lift[, .(gene_row, start_hg38)], by = "gene_row", all.x = TRUE, sort = FALSE)
gene_dt <- merge(gene_dt, gene_end_lift[, .(gene_row, end_hg38)], by = "gene_row", all.x = TRUE, sort = FALSE)
gene_dt[, `:=`(
  start_hg19 = as.integer(start),
  end_hg19 = as.integer(end)
)]
gene_dt <- gene_dt[is.finite(start_hg38) & is.finite(end_hg38)]
gene_dt[, `:=`(
  xmin = pmin(start_hg38, end_hg38),
  xmax = pmax(start_hg38, end_hg38),
  strand_label = fifelse(strand == "-", "<", ">")
)]

meta <- fromJSON(meta_json)
loading <- data.table(
  gene_name = names(meta$cluster_gene_loadings),
  pc_loading = as.numeric(meta$cluster_gene_loadings)
)
effects <- fread(effect_file)
effects <- effects[
  celltype == example$celltype &
    cluster_id == example$cluster_id &
    qtl_type == example$qtl_type &
    phenotype_id == example$qtl_phenotype_id,
  .(gene_name = gene_id, pip_weighted_nominal_beta, abs_pip_weighted_nominal_beta)
]
gene_dt <- merge(gene_dt, loading, by = "gene_name", all.x = TRUE)
gene_dt <- merge(gene_dt, effects, by = "gene_name", all.x = TRUE)
gene_dt <- compact_gene_lanes(gene_dt)

qplot_dt[, point_role := fcase(
  q_idx == q_lead$q_idx, "Lead / target",
  in_cs & shared, "Shared CS",
  in_cs, "QTL CS",
  shared, "Shared",
  default = "Other"
)]
gplot_dt[, point_role := fcase(
  g_idx == g_lead$g_idx, "Lead / target",
  in_cs & shared, "Shared CS",
  in_cs, "GWAS CS",
  shared, "Shared",
  default = "Other"
)]

qtl_ref_key <- qplot_dt[point_role == "Lead / target"][1, display_ld_key(pos, effect_allele, other_allele)]
add_locuszoom_ld_group <- function(dt) {
  dt[, ld_color_group := locuszoom_ld_group(
    lead_r2,
    display_ld_key(pos, effect_allele, other_allele) == qtl_ref_key
  )]
  dt
}
qplot_dt <- add_locuszoom_ld_group(qplot_dt)
gplot_dt <- add_locuszoom_ld_group(gplot_dt)

q_out <- copy(qplot_dt)
g_out <- copy(gplot_dt)
gene_out <- copy(gene_dt)
list_cols <- names(gene_out)[vapply(gene_out, is.list, logical(1))]
if (length(list_cols)) gene_out[, (list_cols) := NULL]

fwrite(q_out, file.path(data_dir, "gimap_pcqtl_pc3_plot_variants_hg38.tsv"), sep = "\t", quote = FALSE)
fwrite(g_out, file.path(data_dir, "gimap_finngen_3019198_plot_variants_hg38.tsv"), sep = "\t", quote = FALSE)
fwrite(gene_out, file.path(data_dir, "gimap_gene_track_hg38_approx.tsv"), sep = "\t", quote = FALSE)
fwrite(coloc_row, file.path(data_dir, "gimap_formal_coloc_row.tsv"), sep = "\t", quote = FALSE)

message("Wrote final GIMAP plot-data tables to ", data_dir)
