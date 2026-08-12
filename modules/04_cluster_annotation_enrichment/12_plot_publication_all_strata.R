#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("config/config.R")

res <- fread(file.path(ROOT_DIR, "results", "enrichment", "enrichment_by_method.tsv"))
out_dir <- Sys.getenv(
  "SC_PCQTL_CLUSTER_ENRICHMENT_FIGURE_OUTPUT_DIR",
  unset = file.path(ROOT_DIR, "figures", "publication")
)
dir_create(out_dir)

ann_labels <- c(
  has_paralog = "Paralogs",
  has_shared_go_bp = "Same GO term",
  has_same_strand_overlap = "Overlapping\n(same strand)",
  has_shared_same_strand_promoter = "Shared promoter\n(same strand)",
  has_shared_abc_enhancer = "Shared enhancer",
  has_shared_opposite_strand_promoter = "Shared promoter\n(opposite strand)",
  has_opposite_strand_overlap = "Overlapping\n(opposite strand)",
  cross_ctcf_peak = "Contains CTCF",
  cross_tad_boundary = "Cross TAD"
)
ann_order <- names(ann_labels)

stratum_labels <- c(
  positive_only = "Positive only",
  mixed = "Mixed",
  negative_only = "Negative only"
)
stratum_colors <- c(
  positive_only = "#B83A4B",
  mixed = "#6F6F6F",
  negative_only = "#69AED1"
)
stratum_offsets <- c(
  positive_only = 0.22,
  mixed = 0,
  negative_only = -0.22
)

sig_label <- function(p) {
  fifelse(is.na(p), "",
    fifelse(p < 1e-3, "***",
      fifelse(p < 1e-2, "**",
        fifelse(p < 5e-2, "*", "ns"))))
}

prepare_plot_dt <- function() {
  dt <- copy(res[
    method == "add_cov_sc_hurdle" &
      stratum %in% names(stratum_labels) &
      annotation %in% ann_order
  ])
  dt[, annotation_label := ann_labels[annotation]]
  dt[, annotation_label := factor(annotation_label, levels = rev(ann_labels[ann_order]))]
  dt[, y_base := as.numeric(annotation_label)]
  dt[, y_plot := y_base + stratum_offsets[stratum]]
  dt[, stratum_label := factor(stratum_labels[stratum], levels = stratum_labels)]
  dt[, color_key := stratum]
  dt[status == "ok", `:=`(
    OR_plot = pmin(pmax(OR, 0.20), 25),
    CI_low_plot = pmin(pmax(CI_low, 0.20), 25),
    CI_high_plot = pmin(pmax(CI_high, 0.20), 25),
    sig = sig_label(pvalue)
  )]
  dt[status != "ok", `:=`(
    OR_plot = 1,
    CI_low_plot = 1,
    CI_high_plot = 1,
    sig = "X"
  )]
  dt[]
}

plot_pcqtl_style <- function(output_prefix, width, height) {
  dt <- prepare_plot_dt()
  ok_dt <- dt[status == "ok"]
  skip_dt <- dt[status != "ok"]

  p <- ggplot() +
    geom_vline(xintercept = 1, linewidth = 0.28, linetype = "dashed", color = "grey45") +
    geom_segment(
      data = ok_dt,
      aes(x = CI_low_plot, xend = CI_high_plot, y = y_plot, yend = y_plot, color = color_key),
      linewidth = 0.42,
      alpha = 0.85
    ) +
    geom_point(
      data = ok_dt,
      aes(x = OR_plot, y = y_plot, color = color_key),
      size = 1.85,
      stroke = 0.35
    ) +
    geom_point(
      data = skip_dt,
      aes(x = 1, y = y_plot),
      shape = 4,
      size = 2.0,
      stroke = 0.55,
      color = "black"
    ) +
    geom_text(
      data = ok_dt,
      aes(x = pmin(CI_high_plot * 1.18, 25), y = y_plot, label = sig, color = color_key),
      size = 2.4,
      hjust = 0,
      show.legend = FALSE
    ) +
    scale_color_manual(values = stratum_colors, labels = stratum_labels, breaks = names(stratum_labels), name = "Cluster correlation type") +
    scale_x_log10(
      breaks = c(0.25, 0.5, 1, 2, 4, 8, 16),
      labels = c("0.25", "0.5", "1", "2", "4", "8", "16"),
      limits = c(0.20, 25)
    ) +
    scale_y_continuous(
      breaks = seq_along(levels(dt$annotation_label)),
      labels = levels(dt$annotation_label),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    labs(
      x = "Odds ratio correlated vs uncorrelated neighboring clusters (log scale)",
      y = NULL,
      caption = "X: skipped because the minimum expected cell count was < 1. Significance: ns p >= 0.05, * p < 0.05, ** p < 0.01, *** p < 0.001."
    ) +
    theme_classic(base_size = 8.5) +
    theme(
      axis.text.y = element_text(size = 7.2, color = "black"),
      axis.text.x = element_text(size = 7.0, color = "black"),
      axis.title.x = element_text(size = 8.5),
      legend.position = "top",
      legend.title = element_text(size = 7.8),
      legend.text = element_text(size = 7.5),
      legend.key.width = unit(0.9, "lines"),
      plot.caption = element_text(size = 6.4, color = "grey30", hjust = 0),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 8),
      panel.spacing.x = unit(1.0, "lines")
    )

  ggsave(file.path(out_dir, paste0(output_prefix, ".pdf")), p, width = width, height = height, useDingbats = FALSE)
  ggsave(file.path(out_dir, paste0(output_prefix, ".png")), p, width = width, height = height, dpi = 450)
}

plot_pcqtl_style("cluster_annotation_enrichment_add_cov", width = 6.8, height = 4.6)

message("Wrote final sc-pcQTL cluster-enrichment figure to ", out_dir)
