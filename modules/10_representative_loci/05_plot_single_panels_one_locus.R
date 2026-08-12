#!/usr/bin/env Rscript

.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
module_dir <- normalizePath(dirname(args_file), mustWork = TRUE)
root <- Sys.getenv("SC_PCQTL_FORMAL_COLOC_ROOT", unset = "")
if (!nzchar(root)) stop("Set SC_PCQTL_FORMAL_COLOC_ROOT to the formal colocalization results root.")
root <- normalizePath(root, mustWork = TRUE)
work_dir <- normalizePath(
  Sys.getenv("SC_PCQTL_LOCUSZOOM_WORK_ROOT", unset = file.path(root, "10_publication_locuszoom_redesign")),
  mustWork = TRUE
)
source(file.path(dirname(module_dir), "06_finngen_susie_coloc/config/config.R"))

manifest_file <- file.path(work_dir, "data/publication_locus_manifest/publication_locus_manifest.tsv")
gencode_file <- file.path(work_dir, "data/gencode/gencode.v38.gene_coordinates.tsv")
gencode_gtf <- file.path(root, "resources/gencode/gencode.v38.annotation.gtf.gz")
effect_file <- file.path(root, "09_mechanistic_celltype_analysis/results/gene_effects/strict_pip_weighted_nominal_gene_effects.tsv")
coloc_file <- file.path(root, "results/coloc/qtl_gwas_susie_official_finngen_all_finemapped/qtl_gwas_susie_official_finngen_coloc_summary.tsv")

figure_id <- Sys.getenv("FIGURE_ID", "high_confidence_shared_mass_pass_cd8_et_SC_chr11_cluster_002_PC2_3009542")
out_root <- file.path(work_dir, "figures/single_panels", figure_id)
data_root <- file.path(work_dir, "data/single_panels", figure_id)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(data_root, recursive = TRUE, showWarnings = FALSE)

manifest <- fread(manifest_file)
target_figure_id <- figure_id
row_dt <- manifest[figure_id == target_figure_id]
if (nrow(row_dt) != 1) stop("Expected exactly one manifest row for FIGURE_ID=", figure_id)
row <- as.list(row_dt[1])

gencode <- fread(gencode_file)
effects_all <- fread(effect_file)
coloc_all <- fread(coloc_file)

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
locuszoom_ld_palette <- c(
  "LD reference" = "#9632b8",
  "0.0-0.2" = "#357ebd",
  "0.2-0.4" = "#46b8da",
  "0.4-0.6" = "#5cb85c",
  "0.6-0.8" = "#eea236",
  "0.8-1.0" = "#d43f3a",
  "No LD" = "#B8B8B8"
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

short_allele <- function(x, max_len = 9L) {
  x <- toupper(as.character(x))
  ifelse(nchar(x) > max_len, paste0(substr(x, 1, max_len), "..."), x)
}

shorten_text <- function(x, max_chars = 48L) {
  x <- gsub("\\s+", " ", as.character(x))
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3L), "..."), x)
}

clean_trait_label <- function(x) {
  x <- gsub("\\s+", " ", as.character(x))
  x <- sub("\\s*\\[[^]]+\\].*$", "", x)
  x <- sub(", definitions combined$", "", x)
  x <- sub(", sample-wise median$", "", x)
  x
}

lead_label <- function(dt) {
  sprintf(
    "chr%s:%.3f Mb (PIP=%.2g)",
    sub("^chr", "", dt$chr), dt$pos / 1e6,
    dt$pip
  )
}

neglog10p <- function(p) {
  p <- as.numeric(p)
  p[!is.finite(p) | p <= 0] <- .Machine$double.xmin
  -log10(p)
}

extract_attr <- function(x, key) {
  pat <- paste0(key, " \"([^\"]+)\"")
  m <- regexec(pat, x, perl = TRUE)
  v <- regmatches(x, m)
  out <- rep(NA_character_, length(x))
  hit <- lengths(v) >= 2
  if (any(hit)) out[hit] <- vapply(v[hit], function(z) z[2], character(1))
  out
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

source(file.path(module_dir, "eur_ld_display_helpers.R"))

get_cs_indices <- function(fit, signal_index) {
  signal_index <- as.integer(signal_index)
  cs_idx <- suppressWarnings(as.integer(fit$sets$cs_index))
  cs_pos <- match(signal_index, cs_idx)
  if (!is.na(cs_pos)) return(as.integer(fit$sets$cs[[cs_pos]]))
  if (length(fit$sets$cs) >= signal_index) return(as.integer(fit$sets$cs[[signal_index]]))
  integer()
}

read_ld_matrix <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (grepl("\\.gz$", path)) {
    as.matrix(fread(cmd = paste("zcat", shQuote(path)), header = FALSE, showProgress = FALSE))
  } else {
    as.matrix(fread(path, header = FALSE, showProgress = FALSE))
  }
}

read_finngen_region_snp <- function(snp_file, chr, start, end, signal_index) {
  alpha_col <- paste0("alpha", as.integer(signal_index))
  cols <- c("region", "chromosome", "position", "allele1", "allele2", "p", "prob", "cs", "lead_r2", alpha_col)
  dt <- fread(cmd = paste("zcat", shQuote(snp_file)), select = cols, showProgress = FALSE)
  if (!alpha_col %in% names(dt)) dt[, (alpha_col) := NA_real_]
  setnames(dt, alpha_col, "fg_alpha")
  dt[chromosome == paste0("chr", sub("^chr", "", chr)) & position >= start & position <= end]
}

compact_gene_lanes <- function(genes, min_gap = 60000L) {
  genes <- copy(genes)
  setorder(genes, start, end)
  lane_end <- numeric()
  genes[, lane := 1L]
  for (i in seq_len(nrow(genes))) {
    placed <- FALSE
    if (length(lane_end)) {
      for (j in seq_along(lane_end)) {
        if (genes$start[i] > lane_end[j] + min_gap) {
          genes$lane[i] <- j
          lane_end[j] <- genes$end[i]
          placed <- TRUE
          break
        }
      }
    }
    if (!placed) {
      lane_end <- c(lane_end, genes$end[i])
      genes$lane[i] <- length(lane_end)
    }
  }
  genes[]
}

read_gencode_exons <- function(gtf_file, chr, start, end, gene_names) {
  if (!file.exists(gtf_file) || !length(gene_names)) return(data.table())
  cmd <- sprintf(
    "zcat %s | awk -v chr=%s -v s=%d -v e=%d 'BEGIN{FS=OFS=\"\\t\"} $1==chr && $3==\"exon\" && $4<=e && $5>=s {print}'",
    shQuote(gtf_file), shQuote(paste0("chr", sub("^chr", "", chr))), as.integer(start), as.integer(end)
  )
  dt <- tryCatch(
    fread(
      cmd = cmd,
      header = FALSE,
      col.names = c("chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
      showProgress = FALSE
    ),
    error = function(e) data.table()
  )
  if (!nrow(dt)) return(dt)
  dt[, `:=`(
    gene_name = extract_attr(attribute, "gene_name"),
    transcript_id = sub("\\..*$", "", extract_attr(attribute, "transcript_id")),
    start = as.integer(start),
    end = as.integer(end)
  )]
  dt <- dt[gene_name %chin% gene_names]
  if (!nrow(dt)) return(dt)

  # Collapse all transcript exons into non-overlapping exon blocks per gene for a compact publication track.
  setorder(dt, gene_name, start, end)
  dt[, prev_max_end := shift(cummax(end)), by = gene_name]
  dt[, exon_group := cumsum(is.na(prev_max_end) | start > prev_max_end + 1L), by = gene_name]
  dt[, .(start = min(start), end = max(end)), by = .(gene_name, exon_group)]
}

theme_panel <- function(base_size = 5.2) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(color = "#111827"),
      axis.text = element_text(color = "#111827"),
      axis.line = element_line(linewidth = 0.22, color = "#111827"),
      axis.ticks = element_line(linewidth = 0.20, color = "#111827"),
      axis.ticks.length = unit(1.2, "mm"),
      panel.grid.major.y = element_line(color = "#e5e7eb", linewidth = 0.16),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 0.4),
      plot.subtitle = element_text(color = "#4b5563", size = base_size - 0.6),
      legend.title = element_text(size = base_size - 0.9),
      legend.text = element_text(size = base_size - 1.0),
      plot.margin = margin(2, 2, 2, 2)
    )
}

make_x_scale <- function(start, end) {
  scale_x_continuous(
    limits = c(start, end),
    labels = function(x) sprintf("%.3f", x / 1e6),
    expand = expansion(mult = c(0.012, 0.012))
  )
}

lead_annotation_layer <- function(lead_dt, plot_start, plot_end, y_max, prefer_right = TRUE) {
  width <- plot_end - plot_start
  lead_dt <- copy(lead_dt)
  lead_dt[, label := lead_label(.SD)]
  right_x <- lead_dt$pos + 0.12 * width
  left_x <- lead_dt$pos - 0.12 * width
  use_right <- prefer_right && right_x < (plot_end - 0.06 * width)
  if (!use_right && left_x < (plot_start + 0.06 * width)) use_right <- TRUE
  lead_dt[, `:=`(
    label_x = if (use_right) right_x else left_x,
    label_y = y_max * 1.05,
    label_hjust = if (use_right) 0 else 1
  )]
  list(
    geom_segment(
      data = lead_dt,
      aes(x = pos, y = mlog10p, xend = label_x, yend = label_y),
      inherit.aes = FALSE,
      linewidth = 0.16,
      color = "#52525b",
      alpha = 0.75
    ),
    geom_text(
      data = lead_dt,
      aes(x = label_x, y = label_y, label = label, hjust = label_hjust),
      inherit.aes = FALSE,
      size = 1.35,
      fontface = "bold",
      color = "#111827"
    )
  )
}

save_panel <- function(plot, name, width, height) {
  png_file <- file.path(out_root, paste0(name, ".png"))
  pdf_file <- file.path(out_root, paste0(name, ".pdf"))
  svg_file <- file.path(out_root, paste0(name, ".svg"))
  ggsave(png_file, plot, width = width, height = height, dpi = 450, bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  svg(svg_file, width = width, height = height, bg = "white")
  print(plot)
  dev.off()
  data.table(panel = name, png = png_file, pdf = pdf_file, svg = svg_file)
}

q <- readRDS(row$qtl_susie_rds)
g <- readRDS(row$gwas_susie_rds)
q_signal <- as.integer(row$qtl_signal_index)
g_signal <- as.integer(row$gwas_signal_index)

qchr <- sub("^chr", "", as.character(q$stats$chr[1]))

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
qtl_lead_pos <- q_lead$pos
gwas_lead_pos <- g_lead$pos

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
if (!is.null(ld) && is.finite(q_lead$ld_order) && q_lead$ld_order <= nrow(ld)) {
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

loading_file <- file.path(PCQTL_ROOT, "celltypes", row$celltype, "pcQTL", "step2_pca", row$cluster_id, "gene_loadings.tsv")
loading <- fread(loading_file)
pc_col <- row$qtl_phenotype_id
loading <- loading[, .(gene_name = Gene, pc_loading = as.numeric(get(pc_col)))]
effect <- effects_all[
  celltype == row$celltype &
    cluster_id == row$cluster_id &
    qtl_type == "pcQTL" &
    phenotype_id == row$qtl_phenotype_id,
  .(gene_name = gene_id, credible_set_id, pip_weighted_nominal_beta, abs_pip_weighted_nominal_beta)
]
effect_target <- effect[grepl(paste0("__CS", q_signal, "$"), credible_set_id)]
if (nrow(effect_target)) effect <- effect_target
effect <- unique(effect[order(-abs(pip_weighted_nominal_beta))], by = "gene_name")
effect[, credible_set_id := NULL]

gene_names <- unique(c(loading$gene_name, effect$gene_name))
gene_dt <- gencode[gene_name %chin% gene_names & chrom == qchr]
gene_dt <- merge(gene_dt, loading, by = "gene_name", all.x = TRUE)
gene_dt <- merge(gene_dt, effect, by = "gene_name", all.x = TRUE)
gene_dt <- gene_dt[order(fifelse(gene_type == "protein_coding", 0L, 1L), start, end)]
gene_dt <- unique(gene_dt, by = "gene_name")
gene_dt[, `:=`(strand_label = fifelse(strand == "-", "<", ">"), xmin = start, xmax = end)]
gene_dt <- compact_gene_lanes(gene_dt)

assoc_start <- min(qstats$pos, na.rm = TRUE)
assoc_end <- max(qstats$pos, na.rm = TRUE)
gene_start <- assoc_start
gene_end <- assoc_end

fg_region <- read_finngen_region_snp(g$metadata$snp_file, qchr, assoc_start, assoc_end, g_signal)
fg_region[, `:=`(
  chr = sub("^chr", "", chromosome),
  effect_allele = toupper(allele1),
  other_allele = toupper(allele2)
)]
fg_region <- fg_region[, .(
  chr,
  pos = as.integer(position),
  fg_effect_allele = effect_allele,
  fg_other_allele = other_allele,
  fg_prob = as.numeric(prob),
  fg_cs = as.integer(cs),
  fg_lead_r2 = as.numeric(lead_r2),
  fg_alpha = as.numeric(fg_alpha)
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
  chr = qchr,
  formal_root = root,
  out_dir = file.path(data_root, "eur_1000g_display_ld")
)
qstats <- eur_ld$qstats
gstats <- eur_ld$gstats
message(sprintf(
  "Using 1000G EUR Phase3 display LD anchor %s (%s reference SNPs in window).",
  ifelse(is.na(eur_ld$anchor_snp), "NA", eur_ld$anchor_snp),
  format(eur_ld$ref_n, big.mark = ",")
))

qbg_dt <- qstats[pos >= assoc_start & pos <= assoc_end]
gbg_dt <- gstats[pos >= assoc_start & pos <= assoc_end]
qplot_dt <- qbg_dt
gplot_dt <- gbg_dt
q_lead_panel <- qplot_dt[q_idx == q_lead$q_idx][1]
if (!nrow(q_lead_panel) || !is.finite(q_lead_panel$pos)) q_lead_panel <- qplot_dt[which.max(alpha)]
g_lead_panel <- gplot_dt[g_idx == g_lead$g_idx][1]
if (!nrow(g_lead_panel) || !is.finite(g_lead_panel$pos)) g_lead_panel <- gplot_dt[which.max(alpha)]
qtl_lead_pos <- q_lead_panel$pos
gwas_lead_pos <- g_lead_panel$pos
qtl_ref_key <- display_ld_key(q_lead_panel$pos, q_lead_panel$effect_allele, q_lead_panel$other_allele)
add_locuszoom_ld_group <- function(dt) {
  dt[, ld_color_group := locuszoom_ld_group(
    lead_r2,
    display_ld_key(pos, effect_allele, other_allele) == qtl_ref_key
  )]
  dt
}
qplot_dt <- add_locuszoom_ld_group(qplot_dt)
gplot_dt <- add_locuszoom_ld_group(gplot_dt)
qbg_dt <- add_locuszoom_ld_group(qbg_dt)
gbg_dt <- add_locuszoom_ld_group(gbg_dt)
gene_plot_dt <- gene_dt[start <= gene_end & end >= gene_start]
exon_plot_dt <- read_gencode_exons(gencode_gtf, qchr, gene_start, gene_end, gene_plot_dt$gene_name)
if (nrow(exon_plot_dt)) {
  exon_plot_dt <- merge(
    exon_plot_dt,
    gene_plot_dt[, .(gene_name, lane, strand)],
    by = "gene_name",
    all.x = TRUE,
    sort = FALSE
  )
  exon_plot_dt <- exon_plot_dt[!is.na(lane)]
  exon_plot_dt[, `:=`(
    plot_start = pmax(start, gene_start),
    plot_end = pmin(end, gene_end)
  )]
}
if (nrow(gene_plot_dt)) {
  gene_width <- gene_end - gene_start
  gene_plot_dt[, `:=`(
    plot_start = pmax(start, gene_start),
    plot_end = pmin(end, gene_end),
    label_x = pmin(pmax((start + end) / 2, gene_start + 0.018 * gene_width), gene_end - 0.018 * gene_width),
    tss_x = pmin(pmax(ifelse(strand == "+", start, end), gene_start), gene_end),
    strand_x = pmin(pmax(ifelse(strand == "+", end, start), gene_start + 0.010 * gene_width), gene_end - 0.010 * gene_width)
  )]
}

qplot_dt[, point_role := fcase(
  q_idx == q_lead$q_idx, "Lead",
  in_cs & shared, "Shared CS",
  in_cs, "QTL CS",
  shared, "Shared",
  default = "Other"
)]
gplot_dt[, point_role := fcase(
  g_idx == g_lead$g_idx, "Lead",
  in_cs & shared, "Shared CS",
  in_cs, "GWAS CS",
  shared, "Shared",
  default = "Other"
)]

assoc_x_scale <- make_x_scale(assoc_start, assoc_end)
gene_x_scale <- make_x_scale(gene_start, gene_end)
ld_scale <- scale_color_manual(
  values = locuszoom_ld_palette,
  breaks = locuszoom_ld_levels,
  drop = FALSE,
  name = expression(display~LD~italic(r)^2),
  guide = guide_legend(override.aes = list(size = 1.9, alpha = 1), title.position = "top")
)

association_panel <- function(dt, lead_dt, title, subtitle, ylab, threshold = NULL, bg_dt = NULL) {
  y_max <- max(c(dt$mlog10p, if (!is.null(bg_dt)) bg_dt$mlog10p else NA_real_), na.rm = TRUE)
  y_top <- y_max * 1.22
  bg_plot_dt <- if (!is.null(bg_dt) && nrow(bg_dt)) bg_dt[is.na(shared) | shared != TRUE] else NULL
  ref_dt <- dt[ld_color_group == "LD reference"]
  marker_dt <- if (nrow(ref_dt)) ref_dt[1] else lead_dt
  p <- ggplot() +
    {if (!is.null(bg_plot_dt) && nrow(bg_plot_dt)) geom_point(
      data = bg_plot_dt,
      aes(pos, mlog10p, color = ld_color_group),
      inherit.aes = FALSE,
      size = 0.46,
      alpha = 0.50,
      stroke = 0
    )} +
    geom_point(data = dt, aes(pos, mlog10p, color = ld_color_group), size = 0.95, alpha = 0.93, stroke = 0) +
    geom_point(data = dt[in_cs == TRUE], aes(pos, mlog10p), shape = 21, fill = NA, color = "#111827", size = 1.28, stroke = 0.20) +
    geom_vline(xintercept = gwas_lead_pos, linetype = "dashed", color = "#9ca3af", linewidth = 0.24) +
    geom_vline(xintercept = qtl_lead_pos, linetype = "dashed", color = "#9ca3af", linewidth = 0.24) +
    geom_point(data = marker_dt, aes(pos, mlog10p), shape = 23, size = 1.95, fill = "#9632b8", color = "#9632b8", stroke = 0.22) +
    lead_annotation_layer(marker_dt, assoc_start, assoc_end, y_max) +
    annotate("text", x = assoc_start + 0.012 * (assoc_end - assoc_start), y = y_top * 0.97,
             label = title, hjust = 0, vjust = 1, size = 1.62, fontface = "bold", color = "#111827") +
    assoc_x_scale + ld_scale +
    scale_y_continuous(limits = c(0, y_top), expand = expansion(mult = c(0, 0.015))) +
    labs(x = paste0("Chromosome ", qchr, " position (Mb, hg38)"), y = ylab) +
    theme_panel(5.2) +
    theme(
      legend.position = c(0.975, 0.76),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = scales::alpha("white", 0.78), color = NA),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.key.width = unit(1.9, "mm"),
      legend.key.height = unit(1.9, "mm")
    )
  if (!is.null(threshold)) p <- p + geom_hline(yintercept = threshold, linetype = "dashed", color = "#7a7a7a", linewidth = 0.22)
  p
}

p_gwas <- association_panel(
  gplot_dt, g_lead_panel,
  paste0(shorten_text(clean_trait_label(row$phenotype), 34), " (FinnGen R12)"),
  paste0(row$phenocode, ": ", row$phenotype, " | official SuSiE signal ", g_signal),
  expression(-log[10](italic(P)[GWAS])),
  -log10(5e-8),
  bg_dt = gbg_dt
)

p_qtl <- association_panel(
  qplot_dt, q_lead_panel,
  paste0(row$qtl_phenotype_id, " pcQTL (", row$celltype, ")"),
  paste0(row$celltype, " ", row$cluster_id, " ", row$qtl_phenotype_id, " | OneK1K LD SuSiE signal ", q_signal),
  expression(-log[10](italic(P)[pcQTL])),
  -log10(5e-8),
  bg_dt = qbg_dt
)

p_gene <- ggplot() +
  geom_segment(data = gene_plot_dt, aes(x = plot_start, xend = plot_end, y = lane, yend = lane), color = "#475569", linewidth = 0.20) +
  {if (nrow(exon_plot_dt)) geom_rect(
    data = exon_plot_dt,
    aes(xmin = plot_start, xmax = plot_end, ymin = lane - 0.085, ymax = lane + 0.085, fill = strand),
    color = NA,
    alpha = 0.95
  ) else geom_rect(
    data = gene_plot_dt,
    aes(xmin = plot_start, xmax = plot_end, ymin = lane - 0.075, ymax = lane + 0.075, fill = strand),
    color = NA,
    alpha = 0.8
  )} +
  geom_segment(data = gene_plot_dt, aes(x = tss_x, xend = tss_x,
                   y = lane - 0.16, yend = lane + 0.16),
               color = "#111827", linewidth = 0.20) +
  geom_text(data = gene_plot_dt, aes(x = label_x, y = lane + 0.24, label = gene_name), size = 1.85, fontface = "bold", color = "#111827") +
  geom_text(data = gene_plot_dt, aes(x = strand_x, y = lane - 0.19, label = strand_label), size = 1.55, color = "#334155") +
  geom_vline(xintercept = gwas_lead_pos, linetype = "dashed", color = "#9ca3af", linewidth = 0.22) +
  geom_vline(xintercept = qtl_lead_pos, linetype = "dashed", color = "#9ca3af", linewidth = 0.22) +
  scale_fill_manual(values = c("+" = "#64748b", "-" = "#94a3b8"), guide = "none") +
  gene_x_scale +
  scale_y_continuous(NULL, breaks = NULL, expand = expansion(mult = c(0.12, 0.28))) +
  labs(x = paste0("Chromosome ", qchr, " position (Mb, hg38)")) +
  theme_panel(5.2) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())

coloc_compare <- coloc_all[
  as.character(phenocode) == as.character(row$phenocode) &
    celltype == row$celltype &
    cluster_id == row$cluster_id &
    status == "ok" &
    qtl_type %chin% c("pcQTL", "eQTL")
]
coloc_compare <- coloc_compare[, .SD[which.max(pph4)], by = .(qtl_type, qtl_phenotype_id)]
coloc_compare[, label := paste(qtl_type, qtl_phenotype_id)]
setorder(coloc_compare, qtl_type, -pph4)
coloc_compare[, label := factor(label, levels = unique(rev(label)))]
p_coloc <- ggplot(coloc_compare, aes(label, pph4, fill = qtl_type)) +
  geom_col(width = 0.62) +
  geom_hline(yintercept = COLOC_PPH4_CUTOFF, linetype = "dashed", linewidth = 0.22, color = "#111827") +
  coord_flip() +
  scale_fill_manual(values = c("pcQTL" = "#7f1d1d", "eQTL" = "#2563eb"), name = "QTL type") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(
    title = "Same-cluster QTL-GWAS colocalization",
    subtitle = sprintf("pcQTL PP.H4=%.3f; shared-mass QC=%s", row$formal_pph4, ifelse(row$shared_mass_pass, "pass", "fail")),
    x = NULL,
    y = "PP.H4"
  ) +
  theme_panel(5.2) +
  theme(
    legend.position = "top",
    legend.key.size = unit(2.6, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.title = element_text(size = 6.4, face = "bold"),
    plot.subtitle = element_text(size = 5.0, color = "#4b5563")
  )

effect_dt <- merge(loading, effect, by = "gene_name", all = TRUE)
effect_dt <- unique(effect_dt[!is.na(gene_name)], by = "gene_name")
effect_dt[, gene_name := factor(gene_name, levels = unique(effect_dt[order(pc_loading)]$gene_name))]
effect_long <- melt(effect_dt[, .(gene_name, pc_loading, pip_weighted_nominal_beta)], id.vars = "gene_name", variable.name = "metric", value.name = "value")
effect_long[, metric := factor(metric, levels = c("pc_loading", "pip_weighted_nominal_beta"), labels = c(paste0(row$qtl_phenotype_id, " loading"), "PIP-weighted nominal expression effect"))]
p_effect <- ggplot(effect_long, aes(gene_name, value, fill = value > 0)) +
  geom_hline(yintercept = 0, color = "#111827", linewidth = 0.20) +
  geom_col(width = 0.62) +
  coord_flip() +
  facet_grid(. ~ metric, scales = "free_x") +
  scale_fill_manual(values = c("TRUE" = "#0f766e", "FALSE" = "#be123c"), guide = "none") +
  labs(title = "PC loading and PIP-weighted expression effect", x = NULL, y = NULL) +
  theme_panel(5.2) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 5.5),
    plot.title = element_text(size = 6.4, face = "bold")
  )

panels <- rbindlist(list(
  save_panel(p_gwas, "01_gwas_association", 6.20, 1.12),
  save_panel(p_qtl, "02_pcqtl_association", 6.20, 1.12),
  save_panel(p_gene, "04_gene_track", 6.20, 0.70),
  save_panel(p_coloc, "05_coloc_pph4_comparison", 3.10, 1.85),
  save_panel(p_effect, "06_pc_loading_gene_effect", 3.40, 1.90)
), fill = TRUE)

fwrite(panels, file.path(data_root, "single_panel_outputs.tsv"), sep = "\t", quote = FALSE)
fwrite(qplot_dt, file.path(data_root, "qtl_plot_variants_hg38.tsv"), sep = "\t", quote = FALSE)
fwrite(gplot_dt, file.path(data_root, "gwas_plot_variants_hg38.tsv"), sep = "\t", quote = FALSE)
fwrite(qbg_dt, file.path(data_root, "qtl_plot_background_variants_hg38.tsv"), sep = "\t", quote = FALSE)
fwrite(gbg_dt, file.path(data_root, "gwas_plot_background_variants_hg38.tsv"), sep = "\t", quote = FALSE)
fwrite(gene_plot_dt, file.path(data_root, "gene_track_gencode_v38.tsv"), sep = "\t", quote = FALSE)
fwrite(exon_plot_dt, file.path(data_root, "gene_track_gencode_v38_exons.tsv"), sep = "\t", quote = FALSE)
fwrite(coloc_compare, file.path(data_root, "same_cluster_coloc_comparison.tsv"), sep = "\t", quote = FALSE)

readme <- c(
  paste0("# Single Panels: ", figure_id),
  "",
  "These are deliberately separate panels. They should be inspected and refined individually before any multi-panel composition.",
  "",
  "- `01_gwas_association`: FinnGen regional GWAS association.",
  "- `02_pcqtl_association`: OneK1K pcQTL regional association.",
  "- `04_gene_track`: GENCODE v38 gene model track.",
  "- `05_coloc_pph4_comparison`: same-cluster eQTL vs pcQTL PP.H4.",
  "- `06_pc_loading_gene_effect`: PC loading and PIP-weighted nominal gene effect.",
  "",
  "Grey dashed lines mark FinnGen and OneK1K lead variants. Open circles mark SuSiE credible-set variants.",
  "Display LD color uses one unified 1000G EUR Phase3 PLINK reference anchored to the pcQTL lead across GWAS and QTL association panels.",
  "Foreground points show allele-compatible shared SNPs used for QTL-GWAS colocalization.",
  "Association panels span the full OneK1K cisQTL fine-mapping window."
)
writeLines(readme, file.path(out_root, "README.md"))
message("Wrote single panels to: ", out_root)
