#!/usr/bin/env Rscript

.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(gridExtra)
  library(scales)
})

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
module_dir <- normalizePath(dirname(args_file), mustWork = TRUE)
release_root <- normalizePath(file.path(module_dir, "..", ".."), mustWork = TRUE)
source(file.path(release_root, "config", "celltype_eligibility.R"))
formal_root <- Sys.getenv("SC_PCQTL_FORMAL_COLOC_ROOT", unset = "")
if (!nzchar(formal_root)) stop("Set SC_PCQTL_FORMAL_COLOC_ROOT to the formal colocalization results root.")
formal_root <- normalizePath(formal_root, mustWork = TRUE)
workflow_root <- normalizePath(
  Sys.getenv("SC_PCQTL_LOCUSZOOM_WORK_ROOT", unset = file.path(formal_root, "10_publication_locuszoom_redesign")),
  mustWork = TRUE
)
single_data_dir <- file.path(workflow_root, "data/single_panels")
manifest_file <- file.path(workflow_root, "data/publication_locus_manifest/publication_locus_manifest.tsv")
section4_dir <- Sys.getenv("SC_PCQTL_LOCUS_OUTPUT_DIR", unset = "")
if (!nzchar(section4_dir)) stop("Set SC_PCQTL_LOCUS_OUTPUT_DIR to the final figure-output directory.")
section4_dir <- normalizePath(section4_dir, mustWork = FALSE)
preview_dir <- file.path(section4_dir, "previews")
dir.create(section4_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)

# Build directional gene-arrow polygons (rectangle + strand-pointing head); y may be per-gene (lane).
gene_arrow_poly <- function(name, x0, x1, strand, start, end, y = 1, h = 0.16,
                            min_frac = 0.007, head_frac = 0.42) {
  win <- end - start
  minw <- win * min_frac
  yv <- rep_len(as.numeric(y), length(name))
  parts <- lapply(seq_along(name), function(i) {
    a <- as.numeric(x0[i]); b <- as.numeric(x1[i])
    if (!is.finite(a) || !is.finite(b)) return(NULL)
    if (b < a) { tmp <- a; a <- b; b <- tmp }
    w <- b - a
    if (w < minw) { m <- (a + b) / 2; a <- m - minw / 2; b <- m + minw / 2; w <- minw }
    head <- min(head_frac * w, minw * 0.6)
    plus <- !identical(as.character(strand[i]), "-")
    px <- if (plus) c(a, a, b - head, b, b - head) else c(b, b, a + head, a, a + head)
    yy <- yv[i]
    py <- c(yy - h, yy + h, yy + h, yy, yy - h)
    data.table(gene_name = name[i], px = px, py = py)
  })
  rbindlist(parts)
}

# Pack genes into stacked rows so neither glyphs nor labels overlap within a row.
assign_gene_lanes <- function(mid, x0, x1, nlab, win, char_frac = 0.0072, gap_frac = 0.004) {
  labw <- nlab * win * char_frac
  ext_lo <- pmin(as.numeric(x0), mid - labw / 2)
  ext_hi <- pmax(as.numeric(x1), mid + labw / 2)
  gap <- win * gap_frac
  lane <- integer(length(mid))
  lane_right <- numeric(0)
  for (i in order(mid)) {
    L <- which(ext_lo[i] > (lane_right + gap))
    if (length(L)) { L <- L[1]; lane[i] <- L; lane_right[L] <- ext_hi[i] }
    else { lane_right <- c(lane_right, ext_hi[i]); lane[i] <- length(lane_right) }
  }
  lane
}

# Load all GENCODE v38 genes overlapping a window, flagging cluster membership for red/blue colouring.
.GENCODE_GENES <- file.path(workflow_root, "data/gencode/gencode.v38.gene_coordinates.tsv")
load_region_genes <- function(chr, lo, hi, cluster_names, types = c("protein_coding")) {
  if (!file.exists(.GENCODE_GENES)) return(data.table())
  g <- fread(.GENCODE_GENES)
  cc <- sub("^chr", "", as.character(chr))
  g <- g[chrom == cc & end >= lo & start <= hi]
  if (length(types)) g <- g[gene_type %chin% types | gene_name %chin% cluster_names]
  g <- unique(g, by = "gene_name")
  g[, is_cluster := gene_name %chin% cluster_names]
  g[, gene_mid := (as.numeric(start) + as.numeric(end)) / 2]
  g[]
}

# Collapsed exon boxes for genes in a window (LocusZoom gene-model glyph: line + exon boxes).
.GENCODE_EXONS <- file.path(workflow_root, "data/gencode/gencode.v38.protein_coding_exons.tsv")
load_region_exons <- function(chr, lo, hi, gene_names) {
  if (!file.exists(.GENCODE_EXONS)) return(data.table(gene_name = character(), exon_start = numeric(), exon_end = numeric()))
  e <- fread(.GENCODE_EXONS)
  cc <- sub("^chr", "", as.character(chr))
  e[chrom == cc & gene_name %chin% gene_names & exon_end >= lo & exon_start <= hi,
    .(gene_name, exon_start = as.numeric(exon_start), exon_end = as.numeric(exon_end))]
}

figure_map <- data.table(
  figure_id = c(
    "high_confidence_shared_mass_pass_cd8_et_SC_chr11_cluster_002_PC2_3009542",
    "high_confidence_shared_mass_pass_cd8_et_SC_chr19_cluster_024_PC2_3007461",
    "high_confidence_shared_mass_pass_b_in_SC_chr17_cluster_004_PC2_T2D"
  ),
  output_base = c(
    "chr11_cd8et_ascl2_c11orf21_tspan32_hematocrit_locus_block",
    "chr19_cd8et_eif3k_actn4_platelets_locus_block",
    "chr17_bin_evi2b_evi2a_t2d_locus_block"
  )
)

manifest <- fread(manifest_file)
manifest <- merge(figure_map, manifest, by = "figure_id", all.x = TRUE, sort = FALSE)
if (anyNA(manifest$celltype)) stop("Some supplement locus figure IDs are missing from the workflow manifest.")
.fig_only <- Sys.getenv("FIG_ONLY", "")
if (nzchar(.fig_only)) {
  .sel <- trimws(strsplit(.fig_only, ",", fixed = TRUE)[[1]])
  manifest <- manifest[output_base %chin% .sel | figure_id %chin% .sel]
  if (!nrow(manifest)) stop("FIG_ONLY matched no rows: ", .fig_only)
}

celltype_label <- function(x) {
  canonical_celltype_labels(x)
}

clean_trait_label <- function(x) {
  x <- gsub("\\s+", " ", as.character(x))
  x <- sub("\\s*\\[[^]]+\\].*$", "", x)
  x <- sub(", definitions combined$", "", x)
  x <- sub(", sample-wise median$", "", x)
  x
}

shorten_text <- function(x, max_chars = 78L) {
  x <- gsub("\\s+", " ", as.character(x))
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3L), "..."), x)
}

neglog10p <- function(p) {
  p <- as.numeric(p)
  p[!is.finite(p) | p <= 0] <- .Machine$double.xmin
  -log10(p)
}

display_ld_key <- function(pos, ea, oa) {
  ea <- toupper(as.character(ea))
  oa <- toupper(as.character(oa))
  allele_pair <- ifelse(ea <= oa, paste0(ea, "/", oa), paste0(oa, "/", ea))
  paste(as.integer(pos), allele_pair, sep = ":")
}

standard_variant_id <- function(chr, pos, ea, oa) {
  paste0(sub("^chr", "", as.character(chr)), ":", as.integer(pos), "_",
         toupper(as.character(ea)), "/", toupper(as.character(oa)))
}

short_variant_label <- function(chr, pos, ea, oa) {
  sprintf("chr%s:%.3f Mb\n%s/%s", sub("^chr", "", chr), as.numeric(pos) / 1e6,
          toupper(as.character(ea)), toupper(as.character(oa)))
}

fmt_num <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

fmt_int <- function(x) {
  x <- suppressWarnings(as.integer(x))
  ifelse(is.finite(x), format(x, big.mark = ",", scientific = FALSE), "NA")
}

get_cs_indices <- function(fit, signal_index) {
  signal_index <- as.integer(signal_index)
  cs_idx <- suppressWarnings(as.integer(fit$sets$cs_index))
  cs_pos <- match(signal_index, cs_idx)
  if (!is.na(cs_pos)) return(as.integer(fit$sets$cs[[cs_pos]]))
  if (length(fit$sets$cs) >= signal_index) return(as.integer(fit$sets$cs[[signal_index]]))
  integer()
}

theme_locus <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size, color = "#17202a"),
      axis.text = element_text(size = base_size - 1.0, color = "#17202a"),
      axis.line = element_line(linewidth = 0.28, color = "#17202a"),
      axis.ticks = element_line(linewidth = 0.28, color = "#17202a"),
      panel.grid.major.y = element_line(color = "#e5e7eb", linewidth = 0.22),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 0.6),
      plot.subtitle = element_text(hjust = 0, size = base_size - 1.5, color = "#475569", margin = margin(t = 1, b = 1)),
      plot.tag = element_text(face = "bold", size = base_size + 1.2, color = "#111827"),
      plot.tag.position = "topleft",
      legend.title = element_text(size = base_size - 0.8),
      legend.text = element_text(size = base_size - 1.0),
      legend.key.height = unit(0.12, "in"),
      legend.key.width = unit(0.18, "in"),
      plot.margin = margin(2, 8, 2, 4)
    )
}

locuszoom_ld_levels <- c("LD reference", "0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0", "No LD")
locuszoom_ld_palette <- c(
  "LD reference" = "#9632b8", "0.0-0.2" = "#357ebd", "0.2-0.4" = "#46b8da",
  "0.4-0.6" = "#5cb85c", "0.6-0.8" = "#eea236", "0.8-1.0" = "#d43f3a", "No LD" = "#B8B8B8"
)
ld_rank <- c("No LD" = 0, "0.0-0.2" = 1, "0.2-0.4" = 2, "0.4-0.6" = 3,
             "0.6-0.8" = 4, "0.8-1.0" = 5, "LD reference" = 6)
locuszoom_ld_group <- function(r2, is_ref) {
  r2 <- suppressWarnings(as.numeric(r2))
  out <- rep("No LD", length(r2))
  ok <- is.finite(r2)
  out[ok] <- as.character(cut(
    pmax(pmin(r2[ok], 1), 0),
    breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
    labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"), right = FALSE
  ))
  out[is_ref %in% TRUE] <- "LD reference"
  factor(out, levels = locuszoom_ld_levels)
}
ld_scale <- function() {
  scale_color_manual(
    values = locuszoom_ld_palette,
    breaks = c("0.8-1.0", "0.6-0.8", "0.4-0.6", "0.2-0.4", "0.0-0.2"),
    labels = c("0.8 - 1.0", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0.0 - 0.2"),
    name = expression(italic(r)^2), drop = FALSE, na.value = "#B8B8B8",
    guide = guide_legend(override.aes = list(shape = 15, size = 2.7, alpha = 1, stroke = 0))
  )
}

make_x_scale <- function(start, end) {
  pad <- (end - start) * 0.015
  scale_x_continuous(
    limits = c(start - pad, end + pad),
    labels = function(x) sprintf("%.2f", x / 1e6),
    expand = expansion(mult = c(0, 0))
  )
}

add_display_ld <- function(dt, ld_map, qtl_ref_key) {
  dt[, display_ld_key := display_ld_key(pos, effect_allele, other_allele)]
  dt <- merge(dt, ld_map, by = "display_ld_key", all.x = TRUE, sort = FALSE)
  dt[, lead_r2_plot := pmax(pmin(as.numeric(display_lead_r2), 1), 0)]
  dt[display_ld_key == qtl_ref_key, lead_r2_plot := 1]
  dt[, is_ld_reference := display_ld_key == qtl_ref_key]
  dt[, ld_color_group := locuszoom_ld_group(display_lead_r2, is_ld_reference)]
  dt[, display_lead_r2 := NULL]
  dt[]
}

association_panel <- function(dt, title, subtitle, ylab, start, end, chr, qtl_lead_pos, gwas_lead_pos,
                              panel_label, lead_fill = "#7f1d1d", show_x = FALSE, show_ld_legend = FALSE,
                              show_lead_label = FALSE, lead_label_prefix = "Lead") {
  dt <- copy(dt)
  dt <- dt[pos >= start & pos <= end]
  dt[, .ld_ord := ld_rank[as.character(ld_color_group)]]
  dt[is.na(.ld_ord), .ld_ord := 0L]
  setorder(dt, .ld_ord)            # low LD underneath, high LD / reference on top
  y_max <- max(dt$mlog10p, na.rm = TRUE)
  if (!is.finite(y_max)) y_max <- 1
  ref_dt <- dt[is_ld_reference == TRUE][1]
  lead_dt <- dt[which.max(mlog10p)][1]
  label_dt <- if (nrow(lead_dt) && !is.na(lead_dt$pos)) copy(lead_dt) else copy(ref_dt)
  if (nrow(label_dt) && !is.na(label_dt$pos)) {
    pw <- end - start
    .rsid <- if ("eur_ref_snp" %in% names(label_dt) && !is.na(label_dt$eur_ref_snp[1]) && nzchar(as.character(label_dt$eur_ref_snp[1])))
               as.character(label_dt$eur_ref_snp[1])
             else paste0("chr", sub("^chr", "", as.character(chr)), ":", sprintf("%.3f", label_dt$pos[1] / 1e6), " Mb")
    .pv <- if ("pvalue" %in% names(label_dt)) suppressWarnings(as.numeric(label_dt$pvalue[1])) else NA_real_
    if (!is.finite(.pv) || .pv <= 0) .pv <- 10^(-suppressWarnings(as.numeric(label_dt$mlog10p[1])))
    .p_str <- if (is.finite(.pv) && .pv > 0) formatC(.pv, format = "e", digits = 1) else "NA"
    label_dt[, callout := paste0(.rsid, " (P = ", .p_str, ")")]
    label_dt[, label_x := pmin(pmax(pos, start + pw * 0.30), end - pw * 0.15)]
    label_dt[, label_y := y_max * 1.02]
  }
  top_expand <- if (show_ld_legend) 0.28 else if (show_lead_label) 0.22 else 0.16
  ggplot(dt, aes(x = pos, y = mlog10p)) +
    geom_vline(xintercept = qtl_lead_pos, linetype = "dashed", color = "#9aa3af", linewidth = 0.26) +
    geom_vline(xintercept = gwas_lead_pos, linetype = "dashed", color = "#c9b6d8", linewidth = 0.26) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "#9ca3af", linewidth = 0.28) +
    geom_point(data = dt[shared != TRUE], aes(color = ld_color_group), size = 0.90, alpha = 0.66, stroke = 0) +
    geom_point(data = dt[shared == TRUE], aes(color = ld_color_group), size = 1.50, alpha = 0.95, stroke = 0) +
    geom_point(data = dt[in_cs == TRUE], aes(x = pos, y = mlog10p),
               shape = 21, fill = NA, color = "#374151", size = 1.70, stroke = 0.22, inherit.aes = FALSE) +
    {if (nrow(lead_dt) && !is.na(lead_dt$pos)) geom_point(data = lead_dt, aes(x = pos, y = mlog10p),
               shape = 24, size = 2.0, fill = lead_fill, color = "white", stroke = 0.30, inherit.aes = FALSE)} +
    {if (nrow(ref_dt) && !is.na(ref_dt$pos)) geom_point(data = ref_dt, aes(x = pos, y = mlog10p),
               shape = 23, size = 2.6, fill = "#9632b8", color = "white", stroke = 0.38, inherit.aes = FALSE)} +
    {if (show_lead_label && nrow(label_dt) && !is.na(label_dt$pos)) geom_segment(
       data = label_dt, aes(x = pos, y = mlog10p, xend = label_x, yend = label_y),
       color = "#94a3b8", linewidth = 0.24, inherit.aes = FALSE)} +
    {if (show_lead_label && nrow(label_dt) && !is.na(label_dt$pos)) geom_text(
       data = label_dt, aes(x = label_x, y = label_y, label = callout),
       vjust = 0, hjust = 0.5, size = 2.75, lineheight = 0.86, color = "#1f2937", inherit.aes = FALSE)} +
    annotate("text", x = start, y = Inf, label = title, hjust = 0, vjust = 1.4, size = 3.10, color = "#111827") +
    ld_scale() +
    make_x_scale(start, end) +
    scale_y_continuous(expand = expansion(mult = c(0.04, top_expand))) +
    labs(tag = panel_label, x = if (show_x) paste0("Chromosome ", chr, " position (Mb, hg38)") else NULL, y = ylab) +
    theme_locus(10.6) +
    theme(
      panel.grid = element_blank(),
      legend.position = if (show_ld_legend) "inside" else "none",
      legend.position.inside = c(0.997, 0.99),
      legend.justification = c(1, 1),
      legend.direction = "vertical",
      legend.title = element_text(size = 9.6),
      legend.text = element_text(size = 8.4),
      legend.key.width = unit(0.11, "in"),
      legend.key.height = unit(0.11, "in"),
      legend.key.spacing.y = unit(0.004, "in"),
      legend.background = element_rect(fill = scales::alpha("white", 0.72), color = NA),
      legend.margin = margin(1, 2, 1, 2),
      legend.key = element_blank(),
      axis.text.x = if (show_x) element_text(color = "#17202a") else element_blank(),
      axis.ticks.x = if (show_x) element_line(linewidth = 0.28, color = "#17202a") else element_blank(),
      axis.line.x = element_line(linewidth = 0.28, color = "#17202a")
    )
}

gene_panel <- function(gene_dt, start, end, chr, qtl_lead_pos, gwas_lead_pos, panel_label) {
  window_start <- as.numeric(start)
  window_end <- as.numeric(end)
  win <- window_end - window_start
  # LocusZoom-style regional gene track: ALL genes in the window, cluster genes red, others blue.
  rg <- load_region_genes(chr, window_start, window_end, gene_dt$gene_name)
  if (!nrow(rg)) {
    rg <- copy(gene_dt)[, .(gene_name, start = window_start, end = window_end, strand,
                            gene_mid = (window_start + window_end) / 2, is_cluster = TRUE)]
  }
  rg <- rg[!grepl("-", gene_name)]                              # drop read-through annotations
  # GENCODE genes may overlap the window while extending beyond it. Clip only the
  # drawing coordinates so ggplot does not discard boundary-spanning segments.
  rg[, draw_start := pmax(as.numeric(get("start")), window_start)]
  rg[, draw_end := pmin(as.numeric(get("end")), window_end)]
  rg[, gene_mid := (draw_start + draw_end) / 2]
  rg[, dir_label := fifelse(strand == "-", paste0("←", gene_name), paste0(gene_name, "→"))]
  # glyph-extent packing keeps the track to a few rows; only cluster genes are named (others blue, unlabelled)
  rg[, lane := assign_gene_lanes(gene_mid, draw_start, draw_end, 0, win, gap_frac = 0.005)]
  n_lane <- max(rg$lane)
  rex <- load_region_exons(chr, window_start, window_end, rg$gene_name)
  rex[, exon_start := pmax(exon_start, window_start)]
  rex[, exon_end := pmin(exon_end, window_end)]
  rex <- merge(rex, rg[, .(gene_name, lane, is_cluster)], by = "gene_name", sort = FALSE)
  # LocusZoom gene-model glyph: gene line + exon boxes + directional names; cluster red, others blue
  p <- ggplot() +
    geom_vline(xintercept = qtl_lead_pos, linetype = "dashed", color = "#9aa3af", linewidth = 0.26) +
    geom_vline(xintercept = gwas_lead_pos, linetype = "dashed", color = "#c9b6d8", linewidth = 0.26) +
    geom_segment(data = rg, aes(x = draw_start, xend = draw_end, y = lane, yend = lane, color = is_cluster),
                 linewidth = 0.42) +
    geom_rect(data = rex, aes(xmin = exon_start, xmax = exon_end, ymin = lane - 0.16, ymax = lane + 0.16,
                              fill = is_cluster), color = NA) +
    ggrepel::geom_text_repel(data = rg[is_cluster == TRUE],
        aes(x = gene_mid, y = lane + 0.20, label = dir_label),
        color = "#B22222", fontface = "bold.italic", size = 2.3,
        segment.color = "#cf9a9a", segment.size = 0.18, min.segment.length = 0,
        box.padding = 0.10, point.padding = 0.10, direction = "both",
        ylim = c(n_lane + 0.35, NA), max.overlaps = Inf, seed = 42) +
    scale_color_manual(values = c("TRUE" = "#B22222", "FALSE" = "#3a51a3"), guide = "none") +
    scale_fill_manual(values = c("TRUE" = "#B22222", "FALSE" = "#5b78c0"), guide = "none") +
    make_x_scale(start, end) +
    scale_y_continuous(NULL, breaks = NULL, limits = c(0.5, n_lane + 1.35), expand = expansion(mult = c(0.05, 0.04))) +
    labs(tag = panel_label, x = paste0("Chromosome ", chr, " position (Mb, hg38)")) +
    theme_locus(10.6) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(color = "#17202a"),
      axis.ticks.x = element_line(linewidth = 0.28, color = "#17202a"),
      axis.line.x = element_line(linewidth = 0.28, color = "#17202a"),
      axis.title.x = element_text(size = 10.7, color = "#17202a")
    )
  attr(p, "n_lane") <- n_lane
  p
}

coloc_panel <- function(dt, panel_label) {
  plot_dt <- copy(dt[qtl_type %chin% c("pcQTL", "eQTL") & status == "ok"])
  plot_dt[, pph4 := suppressWarnings(as.numeric(pph4))]
  plot_dt <- plot_dt[, .SD[which.max(pph4)], by = .(qtl_type, qtl_phenotype_id)]
  setorder(plot_dt, qtl_type, -pph4)
  plot_dt[, label := factor(paste(qtl_type, qtl_phenotype_id), levels = unique(rev(paste(qtl_type, qtl_phenotype_id))))]
  ggplot(plot_dt, aes(x = label, y = pph4, fill = qtl_type)) +
    geom_col(width = 0.68) +
    geom_hline(yintercept = 0.75, linetype = "dashed", linewidth = 0.33, color = "#111827") +
    coord_flip() +
    scale_fill_manual(values = c("pcQTL" = "#7f1d1d", "eQTL" = "#2563eb"), name = "QTL type") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(title = "Same-cluster PPH4", tag = panel_label, x = NULL, y = "PP.H4") +
    theme_locus(12.0) +
    theme(legend.position = "top")
}

effect_panel <- function(gene_dt, pc_id, panel_label) {
  effect_dt <- unique(gene_dt[!is.na(gene_name), .(gene_name, pc_loading, pip_weighted_nominal_beta)], by = "gene_name")
  effect_dt <- effect_dt[is.finite(pc_loading) | is.finite(pip_weighted_nominal_beta)]
  effect_dt[, order_value := fifelse(is.finite(pc_loading), pc_loading, 0)]
  setorder(effect_dt, order_value)
  effect_dt[, gene_name := factor(gene_name, levels = unique(gene_name))]
  effect_long <- melt(
    effect_dt[, .(gene_name, pc_loading, pip_weighted_nominal_beta)],
    id.vars = "gene_name",
    variable.name = "metric",
    value.name = "value"
  )
  effect_long[, metric := factor(
    metric,
    levels = c("pc_loading", "pip_weighted_nominal_beta"),
    labels = c(paste0(pc_id, " loading"), "PIP-weighted effect")
  )]
  ggplot(effect_long, aes(x = gene_name, y = value, fill = value > 0)) +
    geom_hline(yintercept = 0, color = "#111827", linewidth = 0.24) +
    geom_col(width = 0.70) +
    coord_flip() +
    facet_grid(. ~ metric, scales = "free_x") +
    scale_fill_manual(values = c("TRUE" = "#0f766e", "FALSE" = "#be123c"), guide = "none") +
    scale_y_continuous(
      breaks = function(x) pretty(x, n = 3),
      labels = function(x) format(signif(x, 2), trim = TRUE, scientific = FALSE),
      guide = guide_axis(check.overlap = TRUE)
    ) +
    labs(title = "Loading and nominal effect", tag = panel_label, x = NULL, y = NULL) +
    theme_locus(12.0) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.position = "none")
}

prepare_eqtl_plot_data <- function(gene_name, row, gwas_dt, qtl_ld_map, qtl_ref_key, start, end) {
  rds_path <- file.path(formal_root, "results/fine_mapping/qtl", row$celltype, row$cluster_id, "eQTL", paste0(gene_name, ".susie.rds"))
  if (!file.exists(rds_path)) stop("Missing eQTL SuSiE RDS: ", rds_path)
  e <- readRDS(rds_path)
  estats <- as.data.table(copy(e$stats))
  estats[, `:=`(
    e_idx = .I,
    chr = sub("^chr", "", as.character(chr)),
    pos_hg19 = as.integer(pos),
    effect_allele = toupper(as.character(effect_allele)),
    other_allele = toupper(as.character(other_allele))
  )]
  pos_map <- unique(qtl_ld_map[, .(variant_key, pos_hg38 = pos)], by = "variant_key")
  estats <- merge(estats, pos_map, by = "variant_key", all.x = FALSE, all.y = FALSE, sort = FALSE)
  estats[, pos := as.integer(pos_hg38)]

  coloc_gene <- fread(file.path(single_data_dir, row$figure_id, "same_cluster_coloc_comparison.tsv"))[
    qtl_type == "eQTL" & qtl_phenotype_id == gene_name & status == "ok"
  ]
  coloc_gene[, pph4_num := suppressWarnings(as.numeric(pph4))]
  if (nrow(coloc_gene)) coloc_gene <- coloc_gene[which.max(pph4_num)]
  signal_index <- if (nrow(coloc_gene)) suppressWarnings(as.integer(coloc_gene$qtl_signal_index)) else NA_integer_
  if (!is.finite(signal_index) || signal_index < 1 || signal_index > nrow(e$fit$alpha)) {
    signal_index <- if (!is.null(e$fit$sets$cs_index) && length(e$fit$sets$cs_index)) as.integer(e$fit$sets$cs_index[1]) else 1L
  }
  signal_index <- max(1L, min(signal_index, nrow(e$fit$alpha)))
  e_cs <- get_cs_indices(e$fit, signal_index)
  estats[, `:=`(
    mlog10p = neglog10p(pvalue),
    pip = as.numeric(e$fit$pip[e_idx]),
    alpha = as.numeric(e$fit$alpha[signal_index, e_idx]),
    in_cs = e_idx %in% e_cs,
    gene_name = gene_name
  )]
  shared <- merge(
    estats[, .(chr, pos, e_idx, e_effect = effect_allele, e_other = other_allele)],
    gwas_dt[, .(chr = as.character(chr), pos, g_effect = effect_allele, g_other = other_allele)],
    by = c("chr", "pos"),
    allow.cartesian = TRUE
  )
  shared[, allele_matched := e_effect == g_effect & e_other == g_other]
  shared[, allele_flipped := e_effect == g_other & e_other == g_effect]
  shared <- unique(shared[allele_matched | allele_flipped], by = "e_idx")
  estats[, shared := e_idx %in% shared$e_idx]
  estats <- add_display_ld(estats, unique(qtl_ld_map[, .(
    display_ld_key = display_ld_key(pos, effect_allele, other_allele),
    display_lead_r2 = lead_r2_plot,
    eur_ref_snp = as.character(eur_ref_snp)
  )], by = "display_ld_key"), qtl_ref_key)
  estats <- estats[pos >= start & pos <= end]
  estats[, `:=`(
    pph4 = if (nrow(coloc_gene)) suppressWarnings(as.numeric(coloc_gene$pph4)) else NA_real_,
    signal_label = if (nrow(coloc_gene)) paste0("SuSiE signal ", signal_index) else "no tested CS"
  )]
  estats[]
}

save_outputs <- function(grob, base, width, height) {
  .outd <- Sys.getenv("OUT_DIR", "")
  pdf_file <- file.path(if (nzchar(.outd)) .outd else section4_dir, paste0(base, ".pdf"))
  png_file <- file.path(if (nzchar(.outd)) .outd else preview_dir, paste0(base, ".png"))
  cairo_pdf(pdf_file, width = width, height = height, bg = "white")
  grid.draw(grob)
  dev.off()
  png(png_file, width = width, height = height, units = "in", res = 360, bg = "white")
  grid.draw(grob)
  dev.off()
  invisible(c(pdf = pdf_file, preview_png = png_file))
}

for (i in seq_len(nrow(manifest))) {
  row <- as.list(manifest[i])
  data_dir <- file.path(single_data_dir, row$figure_id)
  message("Rendering supplement locus block ", i, "/", nrow(manifest), ": ", row$figure_id)
  qdt <- fread(file.path(data_dir, "qtl_plot_variants_hg38.tsv"))
  gdt <- fread(file.path(data_dir, "gwas_plot_variants_hg38.tsv"))
  genes <- fread(file.path(data_dir, "gene_track_gencode_v38.tsv"))
  coloc <- fread(file.path(data_dir, "same_cluster_coloc_comparison.tsv"))
  qdt[, chr := as.character(chr)]
  gdt[, chr := as.character(chr)]
  chr <- as.character(qdt$chr[1])
  plot_start <- min(c(qdt$pos, gdt$pos), na.rm = TRUE)
  plot_end <- max(c(qdt$pos, gdt$pos), na.rm = TRUE)
  q_lead <- qdt[which.max(alpha)]
  g_lead <- gdt[which.max(alpha)]
  qtl_lead_pos <- q_lead$pos
  gwas_lead_pos <- g_lead$pos
  qtl_ref_key <- display_ld_key(q_lead$pos, q_lead$effect_allele, q_lead$other_allele)

  qdt[, lead_r2_plot := pmax(pmin(as.numeric(lead_r2), 1), 0)]
  gdt[, lead_r2_plot := pmax(pmin(as.numeric(lead_r2), 1), 0)]
  qdt[display_ld_key(pos, effect_allele, other_allele) == qtl_ref_key, lead_r2_plot := 1]
  gdt[display_ld_key(pos, effect_allele, other_allele) == qtl_ref_key, lead_r2_plot := 1]
  qdt[, is_ld_reference := display_ld_key(pos, effect_allele, other_allele) == qtl_ref_key]
  gdt[, is_ld_reference := display_ld_key(pos, effect_allele, other_allele) == qtl_ref_key]
  qdt[, `:=`(shared = as.logical(shared), in_cs = as.logical(in_cs))]
  gdt[, `:=`(shared = as.logical(shared), in_cs = as.logical(in_cs))]
  qdt[, ld_color_group := locuszoom_ld_group(lead_r2, is_ld_reference)]
  gdt[, ld_color_group := locuszoom_ld_group(lead_r2, is_ld_reference)]

  # The final manuscript uses the complete shared OneK1K cis-QTL range.
  ld_map <- unique(qdt[, .(
    display_ld_key = display_ld_key(pos, effect_allele, other_allele),
    display_lead_r2 = lead_r2_plot
  )], by = "display_ld_key")

  gene_names <- genes[order(start), unique(gene_name)]
  eqtl_data <- setNames(lapply(gene_names, prepare_eqtl_plot_data, row = row, gwas_dt = gdt,
                               qtl_ld_map = qdt, qtl_ref_key = qtl_ref_key,
                               start = plot_start, end = plot_end), gene_names)

  trait <- clean_trait_label(row$phenotype)
  ct <- celltype_label(row$celltype)
  p_gwas <- association_panel(
    gdt,
    paste0(shorten_text(sub("\\s*\\(.*$", "", trait), 46L), " (FinnGen R12)"),
    NULL,
    expression(-log[10](italic(P))),
    plot_start, plot_end, chr, qtl_lead_pos, gwas_lead_pos,
    panel_label = "a",
    lead_fill = "#7f1d1d",
    show_x = FALSE,
    show_ld_legend = TRUE,
    show_lead_label = TRUE,
    lead_label_prefix = "FinnGen lead"
  )
  p_pcqtl <- association_panel(
    qdt,
    paste0(row$qtl_phenotype_id, " pcQTL (", ct, ")"),
    NULL,
    expression(-log[10](italic(P))),
    plot_start, plot_end, chr, qtl_lead_pos, gwas_lead_pos,
    panel_label = "b",
    lead_fill = "#1d4ed8",
    show_x = FALSE,
    show_ld_legend = FALSE,
    show_lead_label = TRUE,
    lead_label_prefix = "pcQTL lead"
  )
  p_gene <- gene_panel(genes, plot_start, plot_end, chr, qtl_lead_pos, gwas_lead_pos, "c")
  gene_n_lane <- attr(p_gene, "n_lane"); if (is.null(gene_n_lane)) gene_n_lane <- 1L
  gt_h <- 0.14 * gene_n_lane + 0.78
  p_eqtl <- lapply(seq_along(gene_names), function(j) {
    gene <- gene_names[[j]]
    dt <- eqtl_data[[gene]]
    subtitle <- sprintf("%s | %s | PP.H4 %s",
                        ct,
                        dt$signal_label[1],
                        ifelse(is.finite(dt$pph4[1]), sprintf("%.3g", dt$pph4[1]), "NA"))
    association_panel(
      dt,
      paste0(gene, " eQTL (", ct, ")"),
      NULL,
      expression(-log[10](italic(P))),
      plot_start, plot_end, chr, qtl_lead_pos, gwas_lead_pos,
      panel_label = letters[3L + j],
      lead_fill = "#0f766e",
      show_x = j == length(gene_names),
      show_ld_legend = FALSE,
      show_lead_label = TRUE
    )
  })
  p_coloc <- coloc_panel(coloc, letters[4L + length(gene_names)])
  p_effect <- effect_panel(genes, row$qtl_phenotype_id, letters[5L + length(gene_names)])

  # unify left-column gtable widths so the gene track aligns with the LocusZoom panels' x-axis
  left_plots <- c(list(p_gwas, p_pcqtl, p_gene), p_eqtl)
  left_grobs <- lapply(left_plots, ggplotGrob)
  .maxw <- do.call(grid::unit.pmax, lapply(left_grobs, function(g) g$widths))
  left_grobs <- lapply(left_grobs, function(g) { g$widths <- .maxw; g })
  left_stack <- arrangeGrob(
    grobs = left_grobs,
    ncol = 1,
    heights = c(1, 1, gt_h, rep(1, length(p_eqtl) - 1), 1.28)
  )
  right_stack <- arrangeGrob(
    p_coloc,
    p_effect,
    ncol = 1,
    heights = c(1.25, 1.55)
  )
  figure <- arrangeGrob(
    arrangeGrob(left_stack, right_stack, ncol = 2, widths = c(2.05, 1.05)),
    ncol = 1
  )
  height <- max(7.0, round((2 + gt_h + (length(p_eqtl) - 1) + 1.28) * 1.25, 2))
  save_outputs(figure, row$output_base, width = 10.6, height = height)
}
