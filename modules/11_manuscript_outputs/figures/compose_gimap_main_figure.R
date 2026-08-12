#!/usr/bin/env Rscript

.local_libs <- c(
  Sys.getenv("SC_PCQTL_R_LIBS", unset = ""),
  Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = "")
)
.local_libs <- unlist(strsplit(.local_libs[nzchar(.local_libs)], .Platform$path.sep, fixed = TRUE), use.names = FALSE)
if (length(.local_libs)) .libPaths(c(.local_libs, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(gridExtra)
  library(scales)
})

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- normalizePath(dirname(args_file), mustWork = TRUE)
release_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
source(file.path(release_root, "config", "celltype_eligibility.R"))
formal_root <- Sys.getenv("COQTL_FORMAL_COLOC_ROOT", Sys.getenv("SC_PCQTL_FORMAL_COLOC_ROOT", ""))
if (!dir.exists(formal_root)) {
  stop("Missing formal colocalization workflow root. Set COQTL_FORMAL_COLOC_ROOT or SC_PCQTL_FORMAL_COLOC_ROOT.")
}
formal_root <- normalizePath(formal_root)

manuscript_root <- Sys.getenv("SC_PCQTL_MANUSCRIPT_ROOT", unset = "")
figure_root_override <- Sys.getenv("SC_PCQTL_FIGURE_OUTPUT_ROOT", unset = "")
figure_root <- if (nzchar(figure_root_override)) {
  figure_root_override
} else if (nzchar(manuscript_root)) {
  file.path(manuscript_root, "figures")
} else {
  stop("Set SC_PCQTL_MANUSCRIPT_ROOT or SC_PCQTL_FIGURE_OUTPUT_ROOT.")
}
section4_dir <- Sys.getenv("SC_PCQTL_LOCUS_OUTPUT_DIR", unset = file.path(figure_root, "section4"))
section4_dir <- normalizePath(section4_dir, mustWork = FALSE)
dir.create(section4_dir, recursive = TRUE, showWarnings = FALSE)
locus_work_root <- normalizePath(
  Sys.getenv(
    "SC_PCQTL_LOCUSZOOM_WORK_ROOT",
    unset = file.path(formal_root, "10_publication_locuszoom_redesign")
  ),
  mustWork = FALSE
)

# Build directional gene-arrow polygons (rectangle + strand-pointing head) for a clean gene track.
# y may be a scalar or a per-gene vector (lane index) for LocusZoom-style stacked rows.
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

# Pack genes into stacked rows (LocusZoom-style) so neither glyphs nor labels overlap within a row.
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

# Load all GENCODE v38 genes overlapping a window (LocusZoom-style regional gene context),
# flagging cluster membership for red/blue colouring.
.GENCODE_GENES <- file.path(locus_work_root, "data", "gencode", "gencode.v38.gene_coordinates.tsv")
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
.GENCODE_EXONS <- file.path(locus_work_root, "data", "gencode", "gencode.v38.protein_coding_exons.tsv")
load_region_exons <- function(chr, lo, hi, gene_names) {
  if (!file.exists(.GENCODE_EXONS)) return(data.table(gene_name = character(), exon_start = numeric(), exon_end = numeric()))
  e <- fread(.GENCODE_EXONS)
  cc <- sub("^chr", "", as.character(chr))
  e[chrom == cc & gene_name %chin% gene_names & exon_end >= lo & exon_start <= hi,
    .(gene_name, exon_start = as.numeric(exon_start), exon_end = as.numeric(exon_end))]
}

make_gimap_vertical <- function() {
  source(file.path(formal_root, "config", "config.R"))

  data_dir <- file.path(
    locus_work_root,
    "data/gimap_finngen_style"
  )
  q_file <- file.path(data_dir, "gimap_pcqtl_pc3_plot_variants_hg38.tsv")
  g_file <- file.path(data_dir, "gimap_finngen_3019198_plot_variants_hg38.tsv")
  gene_file <- file.path(data_dir, "gimap_gene_track_hg38_approx.tsv")
  needed <- c(q_file, g_file, gene_file)
  if (!all(file.exists(needed))) {
    stop("Missing GIMAP plot-data TSVs. Run modules/10_representative_loci/run_final_locus_figures.sh first.")
  }

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
    clipped <- pmax(pmin(r2[ok], 1), 0)
    out[ok] <- as.character(cut(
      clipped,
      breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
      labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"),
      right = FALSE
    ))
    out[is_ref %in% TRUE] <- "LD reference"
    factor(out, levels = locuszoom_ld_levels)
  }
  neglog10p <- function(p) {
    p <- as.numeric(p)
    p[!is.finite(p) | p <= 0] <- .Machine$double.xmin
    -log10(p)
  }
  read_liftover_posmap <- function(chr) {
    path <- sprintf(ONEK1K_QTL_LIFTOVER_POSMAP_PATTERN, sub("^chr", "", as.character(chr)))
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
    liftover_bin <- file.path(formal_root, "resources", "liftover_tools", "liftOver")
    chain_file <- file.path(formal_root, "resources", "liftover_tools", "hg19ToHg38.over.chain.gz")
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
    }
    if (any(!is.finite(out$pos_hg38))) {
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
  read_ld_matrix <- function(path) {
    if (!file.exists(path)) return(NULL)
    if (grepl("\\.gz$", path)) {
      as.matrix(fread(cmd = paste("zcat", shQuote(path)), header = FALSE, showProgress = FALSE))
    } else {
      as.matrix(fread(path, header = FALSE, showProgress = FALSE))
    }
  }
  get_cs_indices <- function(fit, signal_index) {
    cs_idx <- suppressWarnings(as.integer(fit$sets$cs_index))
    cs_pos <- match(as.integer(signal_index), cs_idx)
    if (is.na(cs_pos)) {
      if (length(fit$sets$cs) < signal_index) return(integer())
      return(as.integer(fit$sets$cs[[signal_index]]))
    }
    as.integer(fit$sets$cs[[cs_pos]])
  }
  short_allele <- function(x, max_len = 8L) {
    x <- toupper(as.character(x))
    ifelse(nchar(x) > max_len, paste0(substr(x, 1, max_len), "..."), x)
  }
  short_variant_label <- function(chr, pos, ea, oa) {
    sprintf("chr%s:%.3f Mb\n%s/%s", sub("^chr", "", chr), as.numeric(pos) / 1e6,
            short_allele(ea), short_allele(oa))
  }
  lead_callout_label <- function(dt, prefix) {
    if (nrow(dt) == 0) return(dt)
    rsid <- if ("eur_ref_snp" %in% names(dt) && !is.na(dt$eur_ref_snp[1]) && nzchar(as.character(dt$eur_ref_snp[1])))
              as.character(dt$eur_ref_snp[1])
            else paste0("chr", sub("^chr", "", as.character(dt$chr[1])), ":",
                        sprintf("%.3f", as.numeric(dt$pos[1]) / 1e6), " Mb")
    pv <- suppressWarnings(as.numeric(dt$pvalue[1]))
    if (!is.finite(pv) || pv <= 0) pv <- 10^(-suppressWarnings(as.numeric(dt$mlog10p[1])))
    p_str <- if (is.finite(pv) && pv > 0) formatC(pv, format = "e", digits = 1) else "NA"
    dt[, callout_label := paste0(rsid, " (P = ", p_str, ")")]
    dt
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
        legend.title = element_text(size = base_size - 0.8),
        legend.text = element_text(size = base_size - 1.0),
        legend.key.height = unit(0.12, "in"),
        legend.key.width = unit(0.18, "in"),
        plot.tag = element_text(face = "bold", size = base_size + 1.2, color = "#111827"),
        plot.tag.position = "topleft",
        plot.margin = margin(2, 16, 2, 4)
      )
  }
  celltype_label <- canonical_celltype_labels("cd8_nc")
  coloc_row_file <- file.path(data_dir, "gimap_formal_coloc_row.tsv")
  if (!file.exists(coloc_row_file)) stop("Missing GIMAP formal coloc row: ", coloc_row_file)
  gimap_coloc_row <- fread(coloc_row_file)[1]
  fmt_num <- function(x, digits = 3) {
    x <- suppressWarnings(as.numeric(x))
    ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
  }
  fmt_int <- function(x) {
    x <- suppressWarnings(as.integer(x))
    ifelse(is.finite(x), format(x, big.mark = ",", scientific = FALSE), "NA")
  }
  clean_endpoint <- function(x) {
    x <- sub(", sample-wise median$", "", as.character(x))
    x <- gsub("\\[#/volume\\]", "count", x)
    x
  }

  qplot_dt <- fread(q_file)
  gplot_dt <- fread(g_file)
  gene_dt <- fread(gene_file)
  qplot_dt[, `:=`(shared = as.logical(shared), in_cs = as.logical(in_cs))]
  gplot_dt[, `:=`(shared = as.logical(shared), in_cs = as.logical(in_cs))]
  eur_ld_source <- "1000G EUR Phase3 PLINK display reference"
  eur_ld_missing_source <- "not in 1000G EUR Phase3 display reference"
  gimap_eur_ld_map <- rbindlist(list(
    qplot_dt[ld_source_for_plot == eur_ld_source & is.finite(lead_r2), .(
      display_ld_key = display_ld_key(pos, effect_allele, other_allele),
      eur_lead_r2 = as.numeric(lead_r2),
      eur_ref_snp = as.character(eur_ref_snp)
    )],
    gplot_dt[ld_source_for_plot == eur_ld_source & is.finite(lead_r2), .(
      display_ld_key = display_ld_key(pos, effect_allele, other_allele),
      eur_lead_r2 = as.numeric(lead_r2),
      eur_ref_snp = as.character(eur_ref_snp)
    )]
  ), fill = TRUE)
  gimap_eur_ld_map <- unique(gimap_eur_ld_map[order(-eur_lead_r2)], by = "display_ld_key")

  qtl_lead_pos <- qplot_dt[point_role == "Lead / target"][1, pos]
  gwas_lead_pos <- gplot_dt[point_role == "Lead / target"][1, pos]
  qtl_ref_key <- qplot_dt[point_role == "Lead / target"][1, display_ld_key(pos, effect_allele, other_allele)]
  add_locuszoom_ld_group <- function(dt) {
    dt[, `:=`(
      lead_r2_plot = pmax(pmin(as.numeric(lead_r2), 1), 0),
      is_ld_reference = display_ld_key(pos, effect_allele, other_allele) == qtl_ref_key
    )]
    dt[!is.finite(lead_r2_plot), lead_r2_plot := NA_real_]
    dt[is_ld_reference == TRUE, lead_r2_plot := 1]
    dt[, ld_color_group := locuszoom_ld_group(
      lead_r2,
      is_ld_reference
    )]
    dt
  }
  qplot_dt <- add_locuszoom_ld_group(qplot_dt)
  gplot_dt <- add_locuszoom_ld_group(gplot_dt)
  plot_start <- min(c(qplot_dt$pos, gplot_dt$pos), na.rm = TRUE)
  plot_end <- max(c(qplot_dt$pos, gplot_dt$pos), na.rm = TRUE)

  coloc_file <- file.path(
    formal_root,
    "results/coloc/qtl_gwas_susie_official_finngen_all_finemapped/qtl_gwas_susie_official_finngen_coloc_summary.tsv"
  )
  coloc <- fread(coloc_file)
  eqtl_coloc <- coloc[
    as.character(phenocode) == "3019198" &
      celltype == "cd8_nc" &
      cluster_id == "SC_chr7_cluster_001" &
      qtl_type == "eQTL" &
      status == "ok"
  ]
  eqtl_coloc[, pph4_num := suppressWarnings(as.numeric(pph4))]
  eqtl_best_by_gene <- eqtl_coloc[, .SD[which.max(pph4_num)], by = qtl_phenotype_id]
  eqtl_rds_dir <- file.path(
    formal_root,
    "results/fine_mapping/qtl/cd8_nc/SC_chr7_cluster_001/eQTL"
  )

  prepare_eqtl_plot_data <- function(gene_name) {
    coloc_row <- eqtl_best_by_gene[qtl_phenotype_id == gene_name][1]
    rds_path <- if (nrow(coloc_row) == 1 && nzchar(coloc_row$qtl_susie_rds) && file.exists(coloc_row$qtl_susie_rds)) {
      coloc_row$qtl_susie_rds
    } else {
      file.path(eqtl_rds_dir, paste0(gene_name, ".susie.rds"))
    }
    if (!file.exists(rds_path)) stop("Missing GIMAP eQTL SuSiE RDS: ", rds_path)
    e <- readRDS(rds_path)
    estats <- as.data.table(copy(e$stats))
    estats[, `:=`(
      e_idx = .I,
      chr = sub("^chr", "", as.character(chr)),
      pos_hg19 = as.integer(pos),
      effect_allele = toupper(as.character(effect_allele)),
      other_allele = toupper(as.character(other_allele))
    )]
    e_lift <- liftover_qtl_positions(estats[, .(e_idx, chr, pos_hg19)], row_col = "e_idx")
    estats <- merge(estats, e_lift[, .(e_idx, pos_hg38, liftover_source)], by = "e_idx", all.x = FALSE, all.y = FALSE)

    signal_index <- suppressWarnings(as.integer(coloc_row$qtl_signal_index))
    signal_is_tested <- nrow(coloc_row) == 1 && is.finite(signal_index) && signal_index >= 1 &&
      !is.null(e$fit$alpha) && signal_index <= nrow(e$fit$alpha)
    alpha_index <- if (signal_is_tested) {
      signal_index
    } else if (!is.null(e$fit$sets$cs_index) && length(e$fit$sets$cs_index) > 0) {
      suppressWarnings(as.integer(e$fit$sets$cs_index[1]))
    } else {
      1L
    }
    alpha_index <- max(1L, min(alpha_index, nrow(e$fit$alpha)))

    e_cs <- if (signal_is_tested || (!is.null(e$fit$sets$cs_index) && length(e$fit$sets$cs_index) > 0)) {
      get_cs_indices(e$fit, alpha_index)
    } else {
      integer()
    }
    estats[, `:=`(
      pos = as.integer(pos_hg38),
      mlog10p = neglog10p(pvalue),
      pip = as.numeric(e$fit$pip[e_idx]),
      alpha = as.numeric(e$fit$alpha[alpha_index, e_idx]),
      variant_label = standard_variant_id(chr, pos_hg38, effect_allele, other_allele),
      gene_name = gene_name,
      in_cs = e_idx %in% e_cs
    )]
    e_lead <- estats[which.min(pvalue)]

    e_shared <- merge(
      estats[, .(chr, pos, e_idx, e_effect = effect_allele, e_other = other_allele)],
      gplot_dt[, .(chr = as.character(chr), pos, g_effect = effect_allele, g_other = other_allele)],
      by = c("chr", "pos"),
      allow.cartesian = TRUE
    )
    e_shared[, allele_matched := e_effect == g_effect & e_other == g_other]
    e_shared[, allele_flipped := e_effect == g_other & e_other == g_effect]
    e_shared <- unique(e_shared[allele_matched | allele_flipped], by = "e_idx")
    estats[, shared := e_idx %in% e_shared$e_idx]

    estats[, display_ld_key := display_ld_key(pos, effect_allele, other_allele)]
    estats <- merge(estats, gimap_eur_ld_map, by = "display_ld_key", all.x = TRUE, sort = FALSE)
    estats[, `:=`(
      lead_r2 = fifelse(is.finite(eur_lead_r2), pmax(pmin(eur_lead_r2, 1), 0), NA_real_),
      ld_source_for_plot = fifelse(is.finite(eur_lead_r2), eur_ld_source, eur_ld_missing_source)
    )]
    estats[, c("display_ld_key", "eur_lead_r2") := NULL]
    estats <- add_locuszoom_ld_group(estats)
    estats[, point_role := fcase(
      e_idx == e_lead$e_idx, "Lead / target",
      in_cs & shared, "Shared CS",
      in_cs, "eQTL CS",
      shared, "Shared",
      default = "Other"
    )]
    plot_dt <- estats[pos >= plot_start & pos <= plot_end]
    plot_dt[, `:=`(
      signal_is_tested = signal_is_tested,
      pph4 = if (nrow(coloc_row) == 1) suppressWarnings(as.numeric(coloc_row$pph4)) else NA_real_,
      signal_label = if (signal_is_tested) paste0("SuSiE signal ", signal_index) else "no coloc-tested CS"
    )]
    plot_dt
  }

  eqtl_genes <- gene_dt[order(start_hg38), unique(gene_name)]
  eqtl_plot_data <- setNames(lapply(eqtl_genes, prepare_eqtl_plot_data), eqtl_genes)
  best_eqtl_gene <- eqtl_best_by_gene[which.max(pph4_num), qtl_phenotype_id]
  best_eqtl_plot <- eqtl_plot_data[[best_eqtl_gene]]

  # ---- display window: full cis region by default; WIN_START/WIN_END (Mb) can focus it ----
  .ws <- Sys.getenv("WIN_START", ""); .we <- Sys.getenv("WIN_END", "")
  if (nzchar(.ws) && nzchar(.we)) {
    win_start <- suppressWarnings(as.numeric(.ws)) * 1e6
    win_end   <- suppressWarnings(as.numeric(.we)) * 1e6
    if (is.finite(win_start) && is.finite(win_end) && win_end > win_start) {
      plot_start <- win_start; plot_end <- win_end
    }
  }

  ld_rank <- c("No LD" = 0, "0.0-0.2" = 1, "0.2-0.4" = 2, "0.4-0.6" = 3,
               "0.6-0.8" = 4, "0.8-1.0" = 5, "LD reference" = 6)
  ld_scale <- scale_color_manual(
    values = locuszoom_ld_palette,
    breaks = c("0.8-1.0", "0.6-0.8", "0.4-0.6", "0.2-0.4", "0.0-0.2"),
    labels = c("0.8 - 1.0", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0.0 - 0.2"),
    name = expression(italic(r)^2),
    drop = FALSE,
    na.value = "#B8B8B8",
    guide = guide_legend(
      override.aes = list(shape = 15, size = 2.7, alpha = 1, stroke = 0)
    )
  )
  plot_pad <- (plot_end - plot_start) * 0.015
  x_scale <- scale_x_continuous(
    limits = c(plot_start - plot_pad, plot_end + plot_pad),
    labels = function(x) sprintf("%.2f", x / 1e6),
    expand = expansion(mult = c(0, 0))
  )
  make_assoc_simple <- function(dt, title, subtitle, ylab, lead_fill,
                                show_x = FALSE, show_lead_label = TRUE,
                                show_ld_legend = FALSE, panel_note = NULL,
                                base_size = 11.2, panel_label = NULL,
                                lead_label_prefix = "Lead SNP", y_top = NULL) {
    dt <- copy(dt)
    dt <- dt[pos >= plot_start & pos <= plot_end]
    dt[, .ld_ord := ld_rank[as.character(ld_color_group)]]
    dt[is.na(.ld_ord), .ld_ord := 0L]
    setorder(dt, .ld_ord)            # low LD underneath, high LD / reference on top
    lead_dt <- dt[point_role == "Lead / target"][1]
    ref_dt <- dt[is_ld_reference == TRUE]
    if (nrow(ref_dt)) ref_dt <- ref_dt[1]
    label_dt <- if (nrow(lead_dt)) copy(lead_dt) else copy(ref_dt)
    y_max <- max(dt$mlog10p, na.rm = TRUE)
    if (!is.finite(y_max)) y_max <- 1
    if (nrow(label_dt)) {
      label_dt <- lead_callout_label(label_dt, lead_label_prefix)
      panel_width <- plot_end - plot_start
      # right-align the callout near the right edge so it clears the top-left in-panel
      # title (and the top-right LD legend on the first panel)
      label_dt[, label_x := plot_end - panel_width * 0.15]
      label_dt[, label_y := y_max + 0.20]
      label_dt[, label_hjust := 1]
    }
    top_expand <- if (show_ld_legend) 0.28 else if (show_lead_label) 0.22 else 0.16

    ggplot(dt, aes(x = pos, y = mlog10p)) +
      # reference lines: pcQTL lead (anchor) and FinnGen lead positions, plus genome-wide significance
      geom_vline(xintercept = qtl_lead_pos, linetype = "dashed", color = "#9aa3af", linewidth = 0.28) +
      geom_vline(xintercept = gwas_lead_pos, linetype = "dashed", color = "#c9b6d8", linewidth = 0.28) +
      geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "#9ca3af", linewidth = 0.3) +
      # association points coloured by discrete LD bin; allele-compatible shared variants emphasised
      geom_point(data = dt[shared != TRUE], aes(color = ld_color_group),
                 size = 0.95, alpha = 0.68, stroke = 0) +
      geom_point(data = dt[shared == TRUE], aes(color = ld_color_group),
                 size = 1.55, alpha = 0.95, stroke = 0) +
      # credible-set members ringed
      geom_point(data = dt[in_cs == TRUE], aes(x = pos, y = mlog10p),
                 shape = 21, fill = NA, color = "#374151", size = 1.75, stroke = 0.22,
                 inherit.aes = FALSE) +
      # panel-specific lead (filled triangle)
      {if (nrow(lead_dt)) geom_point(data = lead_dt, aes(x = pos, y = mlog10p),
                                     shape = 24, size = 2.1, fill = lead_fill, color = "white", stroke = 0.30,
                                     inherit.aes = FALSE)} +
      # pcQTL lead / LD anchor (purple diamond)
      {if (nrow(ref_dt)) geom_point(data = ref_dt, aes(x = pos, y = mlog10p),
                                    shape = 23, size = 2.7, fill = "#9632b8", color = "white", stroke = 0.38,
                                    inherit.aes = FALSE)} +
      # lead callout (plain text + thin leader) on key panels only
      {if (show_lead_label && nrow(label_dt)) geom_segment(
        data = label_dt, aes(x = pos, y = mlog10p, xend = label_x, yend = label_y),
        color = "#94a3b8", linewidth = 0.24, inherit.aes = FALSE)} +
      {if (show_lead_label && nrow(label_dt)) geom_text(
        data = label_dt, aes(x = label_x, y = label_y, label = callout_label, hjust = label_hjust),
        vjust = 0, size = 3.9, lineheight = 0.88, color = "#1f2937", inherit.aes = FALSE)} +
      # in-panel title, top-left (Fingen style)
      annotate("text", x = plot_start, y = Inf, label = title,
               hjust = 0, vjust = 1.4, size = base_size * 0.31, color = "#111827") +
      ld_scale + x_scale +
      # Shared-axis panels use a COMMON tick interval (every 5) so that, with the
      # proportional panel heights, the physical spacing per -log10(P) unit is identical in
      # both groups: the GWAS/pcQTL panels are simply taller versions of the same scale.
      (if (!is.null(y_top))
         scale_y_continuous(limits = c(0, y_top), breaks = seq(0, y_top, by = 5),
                            expand = expansion(mult = c(0.02, 0.06)))
       else
         scale_y_continuous(expand = expansion(mult = c(0.04, top_expand)))) +
      labs(tag = panel_label,
           x = if (show_x) "Chromosome 7 position (Mb, hg38)" else NULL, y = ylab) +
      theme_locus(base_size) +
      theme(
        panel.grid = element_blank(),
        legend.position = if (show_ld_legend) "inside" else "none",
        legend.position.inside = c(0.997, 0.99),
        legend.justification = c(1, 1),
        legend.direction = "vertical",
        legend.title = element_text(size = base_size - 0.3),
        legend.text = element_text(size = base_size - 1.6),
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

  # Two shared y-axes so signal strength is directly comparable WITHOUT crushing the
  # weak panels onto a single huge scale: the GWAS + pcQTL panels share one axis (large
  # range), and the same-cluster eQTL panels share a second axis (small range). Panel
  # HEIGHTS are then made proportional to each shared range (constant -log10(P) per inch),
  # so the strong panels are physically taller and the weak panels stay short -- no empty
  # headroom, no compression. The proportional heights are computed at assembly below.
  in_win <- function(d) d[pos >= plot_start & pos <= plot_end, mlog10p]
  g1_max <- max(c(in_win(gplot_dt), in_win(qplot_dt)), na.rm = TRUE)        # GWAS + pcQTL
  g2_max <- max(unlist(lapply(eqtl_plot_data, in_win)), na.rm = TRUE)       # same-cluster eQTLs
  if (!is.finite(g1_max) || g1_max <= 0) g1_max <- 1
  if (!is.finite(g2_max) || g2_max <= 0) g2_max <- 1
  # Reserve a fixed vertical band for lead-variant callouts. A proportional
  # expansion alone can clip long labels in the shorter eQTL panels.
  g1_top <- g1_max + 2.4
  g2_top <- g2_max + 2.4
  p_gwas2 <- make_assoc_simple(
    gplot_dt,
    "Lymphocyte count (FinnGen R12)",
    NULL,
    expression(-log[10](italic(P))),
    "#7f1d1d",
    show_x = FALSE,
    show_ld_legend = TRUE,
    panel_note = NULL,
    base_size = 11.2,
    panel_label = "a",
    lead_label_prefix = "FinnGen lead",
    y_top = g1_top
  )
  p_pcqtl2 <- make_assoc_simple(
    qplot_dt,
    paste0("GIMAP PC3 pcQTL (", celltype_label, ")"),
    NULL,
    expression(-log[10](italic(P))),
    "#1d4ed8",
    show_x = FALSE,
    panel_note = NULL,
    base_size = 11.2,
    panel_label = "b",
    lead_label_prefix = "pcQTL lead",
    y_top = g1_top
  )
  p_eqtl_list <- lapply(seq_along(eqtl_genes), function(i) {
    gene <- eqtl_genes[[i]]
    dt <- eqtl_plot_data[[gene]]
    subtitle <- sprintf("%s | %s | PP.H4 %s",
                        celltype_label,
                        dt$signal_label[1],
                        ifelse(is.finite(dt$pph4[1]), sprintf("%.3g", dt$pph4[1]), "NA"))
    make_assoc_simple(
      dt,
      paste0(gene, " eQTL (", celltype_label, ")"),
      NULL,
      expression(-log[10](italic(P))),
      "#0f766e",
      show_x = i == length(eqtl_genes),
      show_lead_label = TRUE,
      base_size = 11.2,
      panel_label = letters[3L + i],
      lead_label_prefix = "eQTL lead",
      y_top = g2_top
    )
  })

  hurdle_matrix_tsv <- Sys.getenv(
    "SC_PCQTL_GIMAP_HURDLE_MATRIX",
    unset = if (nzchar(manuscript_root)) file.path(manuscript_root, "data", "gimap_pairwise_hurdle_matrix.tsv") else ""
  )
  if (!file.exists(hurdle_matrix_tsv)) stop("Missing GIMAP hurdle matrix: ", hurdle_matrix_tsv)
  hurdle_genes <- c("GIMAP7", "GIMAP4", "GIMAP5", "GIMAP1", "GIMAP2", "GIMAP6", "GIMAP8")
  hurdle_dt <- fread(hurdle_matrix_tsv)
  hurdle_z_lim <- ceiling(max(c(hurdle_dt$z_count_abs, hurdle_dt$z_zero_abs), na.rm = TRUE))
  if (!is.finite(hurdle_z_lim) || hurdle_z_lim <= 0) hurdle_z_lim <- 1
  make_hurdle_heatmap <- function(dt, value_col, title, colors, z_lim,
                                  show_x = FALSE, panel_label = NULL) {
    plot_dt <- copy(dt)
    plot_dt[, `:=`(
      gene_row_idx = match(gene_row, hurdle_genes),
      gene_col_idx = match(gene_col, hurdle_genes)
    )]
    plot_dt <- plot_dt[gene_row_idx >= gene_col_idx]
    plot_dt[, `:=`(
      gene_row = factor(gene_row, levels = rev(hurdle_genes)),
      gene_col = factor(gene_col, levels = hurdle_genes),
      fill_value = fifelse(gene_row == gene_col, NA_real_, as.numeric(get(value_col))),
      is_diag = gene_row == gene_col
    )]
    support_n <- sum(!is.na(plot_dt$fill_value) & plot_dt$fill_value > 0, na.rm = TRUE)
    support_total <- choose(length(hurdle_genes), 2)
    support_frac <- sprintf("%.1f", 100 * support_n / support_total)
    ggplot(plot_dt, aes(gene_col, gene_row)) +
      geom_tile(aes(fill = fill_value), color = "white", linewidth = 0.45, na.rm = FALSE) +
      geom_tile(
        data = plot_dt[is_diag == TRUE],
        fill = "#d9d9d9",
        color = "white",
        linewidth = 0.45
      ) +
      scale_fill_gradientn(
        colors = colors,
        limits = c(0, z_lim),
        breaks = pretty(c(0, z_lim), n = 5),
        oob = scales::squish,
        na.value = "#f5f5f5",
        name = expression("|z|")
      ) +
      coord_fixed() +
      labs(
        title = title,
        tag = panel_label,
        x = NULL,
        y = NULL,
        caption = paste0(
          support_n, "/", support_total, " pairs with non-missing |z|"
        )
      ) +
      guides(fill = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        direction = "vertical",
        barwidth = unit(0.38, "cm"),
        barheight = unit(3.8, "cm"),
        ticks.colour = "#4d4d4d"
      )) +
      theme_minimal(base_size = 12.0) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = if (show_x) {
          element_text(angle = 45, hjust = 1, face = "italic", colour = "black", size = 10.5)
        } else {
          element_blank()
        },
        axis.text.y = element_text(face = "italic", colour = "black", size = 10.5),
        axis.title = element_blank(),
        plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
        plot.tag = element_text(face = "bold", size = 10.4, color = "#111827"),
        plot.tag.position = "topleft",
        legend.position = "right",
        legend.title = element_text(size = 10.0),
        legend.text = element_text(size = 9.0),
        plot.caption = element_text(hjust = 0, size = 9.3, margin = margin(t = 6)),
        plot.margin = margin(2, 3, 2, 3)
      )
  }
  p_hurdle_count <- make_hurdle_heatmap(
    hurdle_dt,
    "z_count_abs",
    "Count component",
    c("#fff5eb", "#fdd0a2", "#fd8d3c", "#d94801", "#7f2704"),
    hurdle_z_lim,
    show_x = FALSE,
    panel_label = letters[4L + length(eqtl_genes)]
  )
  p_hurdle_detect <- make_hurdle_heatmap(
    hurdle_dt,
    "z_zero_abs",
    "Detection component",
    c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    hurdle_z_lim,
    show_x = TRUE
  )
  p_hurdle2 <- arrangeGrob(
    p_hurdle_count,
    p_hurdle_detect,
    ncol = 1,
    heights = c(1, 1.08)
  )

  gene_dt[, gene_mid := (xmin + xmax) / 2]
  # LocusZoom gene-model track: gene body line + exon boxes + directional gene names; cluster red, others blue.
  reg_genes <- load_region_genes(7, plot_start, plot_end, gene_dt$gene_name)
  if (!nrow(reg_genes)) reg_genes <- gene_dt[, .(gene_name, start = xmin, end = xmax, strand,
                                                 is_cluster = TRUE, gene_mid)]
  reg_genes <- reg_genes[!grepl("-", gene_name)]               # drop read-through annotations
  reg_genes[, `:=`(
    draw_start = pmax(as.numeric(start), plot_start),
    draw_end = pmin(as.numeric(end), plot_end)
  )]
  reg_genes[, gene_mid := (draw_start + draw_end) / 2]
  reg_genes[, dir_label := fifelse(strand == "-", paste0("←", gene_name), paste0(gene_name, "→"))]
  # Pack gene glyphs and cluster-gene labels independently. Explicit label lanes
  # keep the seven GIMAP names legible and reproducible without repel collisions.
  reg_genes[, lane := assign_gene_lanes(gene_mid, draw_start, draw_end, 0, plot_end - plot_start, gap_frac = 0.005)]
  gene_n_lane <- max(reg_genes$lane)
  cluster_label_dt <- copy(reg_genes[is_cluster == TRUE])
  cluster_label_dt[, label_lane := assign_gene_lanes(
    gene_mid, draw_start, draw_end, nchar(dir_label), plot_end - plot_start,
    char_frac = 0.013, gap_frac = 0.004
  )]
  cluster_label_dt[, label_y := gene_n_lane + 0.55 + (label_lane - 1) * 0.52]
  gene_label_top <- max(cluster_label_dt$label_y, na.rm = TRUE)
  gt_h <- 0.14 * gene_n_lane + 0.32 * max(cluster_label_dt$label_lane) + 0.78
  gene_exons <- load_region_exons(7, plot_start, plot_end, reg_genes$gene_name)
  gene_exons <- merge(gene_exons, reg_genes[, .(gene_name, lane, is_cluster)], by = "gene_name", sort = FALSE)
  gene_exons[, `:=`(
    exon_start = pmax(exon_start, plot_start),
    exon_end = pmin(exon_end, plot_end)
  )]
  p_gene2 <- ggplot() +
    geom_vline(xintercept = qtl_lead_pos, linetype = "dashed", color = "#9aa3af", linewidth = 0.28) +
    geom_vline(xintercept = gwas_lead_pos, linetype = "dashed", color = "#c9b6d8", linewidth = 0.28) +
    geom_segment(data = reg_genes, aes(x = draw_start, xend = draw_end, y = lane, yend = lane, color = is_cluster),
                 linewidth = 0.42) +
    geom_rect(data = gene_exons, aes(xmin = exon_start, xmax = exon_end, ymin = lane - 0.16, ymax = lane + 0.16,
                                     fill = is_cluster), color = NA) +
    geom_segment(data = cluster_label_dt,
        aes(x = gene_mid, xend = gene_mid, y = lane + 0.20, yend = label_y - 0.10),
        color = "#cf9a9a", linewidth = 0.18) +
    geom_text(data = cluster_label_dt,
        aes(x = gene_mid, y = label_y, label = dir_label),
        color = "#B22222", fontface = "bold.italic", size = 2.45) +
    scale_color_manual(values = c("TRUE" = "#B22222", "FALSE" = "#3a51a3"), guide = "none") +
    scale_fill_manual(values = c("TRUE" = "#B22222", "FALSE" = "#5b78c0"), guide = "none") +
    x_scale +
    scale_y_continuous(NULL, breaks = NULL, limits = c(0.5, gene_label_top + 0.35),
                       expand = expansion(mult = c(0.05, 0.04))) +
    labs(tag = "c", x = "Chromosome 7 position (Mb, hg38)") +
    theme_locus(11.2) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(color = "#17202a"),
      axis.ticks.x = element_line(linewidth = 0.28, color = "#17202a"),
      axis.line.x = element_line(linewidth = 0.28, color = "#17202a"),
      axis.title.x = element_text(size = 11.2, color = "#17202a")
    )

  coloc_compare <- coloc[
    as.character(phenocode) == "3019198" &
      celltype == "cd8_nc" &
      cluster_id == "SC_chr7_cluster_001" &
      status == "ok" &
      qtl_type %chin% c("pcQTL", "eQTL")
  ]
  coloc_compare <- coloc_compare[, .SD[which.max(pph4)], by = .(qtl_type, qtl_phenotype_id)]
  coloc_compare <- coloc_compare[, .(qtl_type, gene = qtl_phenotype_id,
                                     pph4 = suppressWarnings(as.numeric(pph4)), tested = TRUE)]
  # Same-cluster single-gene eQTLs with no SuSiE credible set are not colocalization-tested
  # (no PP.H4). Show every cluster gene by marking these explicitly as "no credible set".
  .cs_n <- function(g) {
    f <- file.path(eqtl_rds_dir, paste0(g, ".credible_sets.tsv"))
    if (!file.exists(f)) return(NA_integer_)
    suppressWarnings(tryCatch(nrow(fread(f)), error = function(e) NA_integer_))
  }
  no_cs_genes <- setdiff(eqtl_genes, coloc_compare[qtl_type == "eQTL", gene])
  no_cs_genes <- no_cs_genes[vapply(no_cs_genes, function(g) isTRUE(.cs_n(g) == 0L), logical(1))]
  if (length(no_cs_genes)) {
    coloc_compare <- rbind(coloc_compare,
      data.table(qtl_type = "eQTL", gene = no_cs_genes, pph4 = 0, tested = FALSE))
  }
  coloc_compare[, label := paste(qtl_type, gene)]
  setorder(coloc_compare, qtl_type, -pph4)
  coloc_compare[, label := factor(label, levels = rev(label))]
  no_cs_dt <- coloc_compare[tested == FALSE]
  no_cs_layers <- if (nrow(no_cs_dt)) list(
    geom_point(data = no_cs_dt, aes(x = label, y = 0), shape = 4, size = 2.8,
               stroke = 1.0, color = "#9ca3af", inherit.aes = FALSE),
    geom_text(data = no_cs_dt, aes(x = label, y = 0.05, label = "no credible set"),
              hjust = 0, size = 4.0, fontface = "italic", color = "#6b7280", inherit.aes = FALSE)
  ) else NULL
  p_coloc2 <- ggplot(coloc_compare, aes(x = label, y = pph4, fill = qtl_type)) +
    geom_col(width = 0.68) +
    geom_hline(yintercept = 0.75, linetype = "dashed", linewidth = 0.33, color = "#111827") +
    no_cs_layers +
    coord_flip() +
    scale_fill_manual(values = c("pcQTL" = "#7f1d1d", "eQTL" = "#2563eb"), name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(title = "Same-cluster PPH4", tag = letters[5L + length(eqtl_genes)], x = NULL, y = "PP.H4") +
    theme_locus(11.3) +
    theme(legend.position = "top")

  effect_dt <- gene_dt[, .(gene_name, pc_loading, pip_weighted_nominal_beta)]
  effect_dt[, gene_name := factor(gene_name, levels = gene_dt[order(pc_loading)]$gene_name)]
  effect_long <- melt(effect_dt, id.vars = "gene_name", variable.name = "metric", value.name = "value")
  effect_long[, metric := factor(metric,
                                 levels = c("pc_loading", "pip_weighted_nominal_beta"),
                                 labels = c("PC3 loading", "PIP effect"))]
  p_effect2 <- ggplot(effect_long, aes(x = gene_name, y = value, fill = value > 0)) +
    geom_hline(yintercept = 0, color = "#111827", linewidth = 0.24) +
    geom_col(width = 0.7) +
    coord_flip() +
    facet_grid(. ~ metric, scales = "free_x") +
    scale_fill_manual(values = c("TRUE" = "#0f766e", "FALSE" = "#be123c"), guide = "none") +
    scale_y_continuous(
      breaks = function(x) pretty(x, n = 3),
      labels = function(x) format(signif(x, 2), trim = TRUE, scientific = FALSE),
      guide = guide_axis(check.overlap = TRUE)
    ) +
    labs(title = "Loading and nominal effect", tag = letters[6L + length(eqtl_genes)], x = NULL, y = NULL) +
    theme_locus(11.3) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.position = "none")

  # convert to gtables and unify column widths so every panel's plotting area (and
  # therefore the genomic x-axis) aligns vertically -- the gene track lines up with the LocusZoom panels
  left_plots <- c(list(p_gwas2, p_pcqtl2, p_gene2), p_eqtl_list)
  left_grobs <- lapply(left_plots, ggplotGrob)
  .maxw <- do.call(grid::unit.pmax, lapply(left_grobs, function(g) g$widths))
  left_grobs <- lapply(left_grobs, function(g) { g$widths <- .maxw; g })
  # Panel heights proportional to each group's shared y-range: a fixed per-panel
  # decoration allowance (title + margins, ~constant) plus a term proportional to the
  # shared axis top, using the SAME -log10(P)-per-inch slope for both groups. This makes
  # the GWAS/pcQTL panels physically taller in proportion to their larger signal while the
  # eQTL panels stay short and uniform -- the visual height of a given -log10(P) is the
  # same everywhere, so strengths are directly comparable with no empty headroom.
  # The decoration allowances below are calibrated against the actual rendered title/axis
  # heights so that, after subtracting decoration, each panel's DATA REGION is exactly
  # .u * (its shared axis range): the -log10(P)-per-inch scale (hence the tick spacing) is
  # then identical across both groups -- the GWAS/pcQTL panels are taller in exact
  # proportion to their larger range, not merely approximately.
  .u <- 0.075                                   # plotting-area height per -log10(P) unit
  h_assoc  <- 0.62 + .u * g1_top                # GWAS / pcQTL (shared large axis)
  h_eqtl   <- max(0.92, 0.46 + .u * g2_top)     # each eQTL (shared small axis), legible floor
  h_eqtl_last <- h_eqtl + 0.48                  # bottom eQTL panel also carries the x-axis
  left_h <- c(h_assoc, h_assoc, gt_h,
              rep(h_eqtl, length(p_eqtl_list) - 1), h_eqtl_last)
  left_stack <- arrangeGrob(
    grobs = left_grobs,
    ncol = 1,
    heights = left_h
  )
  right_stack <- arrangeGrob(
    p_hurdle2, p_coloc2, p_effect2,
    ncol = 1,
    heights = c(2.8, 1.85, 1.8)
  )
  figure <- arrangeGrob(
    arrangeGrob(left_stack, right_stack, ncol = 2, widths = c(2.38, 1.28)),
    ncol = 1
  )

  base <- file.path(section4_dir, "gimap_3019198_susie_locus_vertical")
  # scale total height so the taller (multi-lane) gene track does not compress the other panels
  left_units <- sum(left_h)
  # fill (most of) a manuscript page at width=\textwidth (caption leaves ~1.5in below)
  fig_h <- round(left_units * 1.36, 2)
  # Render at a smaller physical canvas (aspect preserved) so every text element is
  # larger relative to the page when the figure is placed at \textwidth.
  sv_w <- 9.0; sv_h <- round(fig_h * sv_w / 13.2, 2)
  # Keep the figure at its full proportional height (no compression of the two-group
  # heights), but bound the displayed image (at width = \textwidth = 6.5in) to just under a
  # page so the image sits alone on a dedicated float page; the detailed legend then flows
  # onto the following page (\ContinuedFloat in main.tex). Width (hence text size) is fixed.
  disp_h <- sv_h / sv_w * 6.5
  if (disp_h > 8.5) sv_h <- round(8.5 / 6.5 * sv_w, 2)
  ggsave(paste0(base, ".pdf"), figure, width = sv_w, height = sv_h, device = cairo_pdf, bg = "white")
  ggsave(paste0(base, ".png"), figure, width = sv_w, height = sv_h, dpi = 450, bg = "white")
  svg(paste0(base, ".svg"), width = sv_w, height = sv_h, bg = "white")
  grid.draw(figure)
  dev.off()
  message("Wrote ", paste0(base, ".pdf"))
  message("Wrote ", paste0(base, ".png"))
  message("Wrote ", paste0(base, ".svg"))
}

fig_target <- Sys.getenv("FIG_TARGET", "gimap")
if (!fig_target %in% c("all", "gimap")) {
  stop("FIG_TARGET must be 'gimap' or 'all'; the strict SuSiE scatter is produced by compose_fig3_cluster_and_coloc.R")
}
make_gimap_vertical()
