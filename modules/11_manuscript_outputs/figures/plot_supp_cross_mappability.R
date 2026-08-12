#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

extra_libs <- Sys.getenv("SC_PCQTL_R_LIBS", unset = "")
if (nzchar(extra_libs)) {
  .libPaths(c(strsplit(extra_libs, .Platform$path.sep, fixed = TRUE)[[1]], .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

crossmap_root <- Sys.getenv("SC_PCQTL_CROSSMAP_ROOT", unset = "")
manuscript_root <- Sys.getenv("SC_PCQTL_MANUSCRIPT_ROOT", unset = "")
figure_root <- Sys.getenv("SC_PCQTL_FIGURE_OUTPUT_ROOT", unset = "")

if (!nzchar(crossmap_root) || !dir.exists(crossmap_root)) {
  stop("Set SC_PCQTL_CROSSMAP_ROOT to the fixed cross-mappability workflow root.")
}
if (!nzchar(figure_root)) {
  if (!nzchar(manuscript_root) || !dir.exists(manuscript_root)) {
    stop("Set SC_PCQTL_MANUSCRIPT_ROOT or SC_PCQTL_FIGURE_OUTPUT_ROOT.")
  }
  figure_root <- file.path(manuscript_root, "figures")
}

table_dir <- file.path(crossmap_root, "results", "tables")
plot_dir <- file.path(figure_root, "section5")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

distribution_file <- file.path(
  table_dir,
  "crossmap_genes_per_cluster_distribution.tsv"
)
credible_set_file <- file.path(
  table_dir,
  "pcqtl_credible_set_by_celltype_crossmap_status.tsv"
)
strict_group_file <- file.path(
  table_dir,
  "strict_group_class_counts_crossmap_exclusion.tsv"
)
if (
  !file.exists(distribution_file) ||
  !file.exists(credible_set_file) ||
  !file.exists(strict_group_file)
) {
  stop("Required cross-mappability source table is missing.")
}

distribution <- fread(distribution_file)
credible_sets <- fread(credible_set_file)
strict_groups <- fread(strict_group_file)
credible_sets <- credible_sets[
  crossmap_status %in% c("crossmap_clean", "crossmap_flagged")
]
stopifnot(
  nrow(credible_sets) == 20L,
  sum(credible_sets$n_selected_pc_phenotypes) == 4430L,
  all(table(credible_sets$crossmap_status) == 10L)
)
credible_sets[, status := factor(
  crossmap_status,
  levels = c("crossmap_clean", "crossmap_flagged"),
  labels = c("Absent", "Present")
)]

clean_fraction <- credible_sets[
  crossmap_status == "crossmap_clean",
  credible_set_fraction
]
cross_mappable_fraction <- credible_sets[
  crossmap_status == "crossmap_flagged",
  credible_set_fraction
]
independent_test <- t.test(
  clean_fraction,
  cross_mappable_fraction,
  paired = FALSE
)

theme_crossmap <- function() {
  theme_classic(base_size = 7.5, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.34, colour = "black"),
      axis.ticks = element_line(linewidth = 0.34, colour = "black"),
      axis.text = element_text(colour = "black", size = 6.7),
      axis.title = element_text(colour = "black", size = 7.5),
      panel.grid = element_blank(),
      plot.tag = element_text(face = "bold", size = 9.5, colour = "black"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(6, 5, 4, 5)
    )
}

panel_a <- ggplot(
  distribution,
  aes(x = n_crossmap_genes_gt100, y = N)
) +
  geom_col(
    width = 0.76,
    fill = "#6F99B8",
    colour = "black",
    linewidth = 0.35
  ) +
  scale_x_continuous(
    limits = c(-0.5, 10),
    breaks = seq(0, 10, by = 2),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    breaks = seq(0, 2000, by = 500),
    expand = expansion(mult = c(0, 0.055))
  ) +
  labs(
    x = "Genes per cluster with\ncross-mappability score >100",
    y = "Gene clusters",
    tag = "A"
  ) +
  theme_crossmap() +
  theme(plot.margin = margin(6, 8, 4, 5))

bracket_y <- 0.90
panel_b <- ggplot(
  credible_sets,
  aes(x = status, y = credible_set_fraction, fill = status)
) +
  geom_boxplot(
    width = 0.54,
    colour = "black",
    linewidth = 0.45,
    outlier.shape = 16,
    outlier.size = 1.5,
    outlier.colour = "black"
  ) +
  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = bracket_y,
    yend = bracket_y,
    linewidth = 0.4
  ) +
  annotate(
    "segment",
    x = c(1, 2),
    xend = c(1, 2),
    y = bracket_y,
    yend = bracket_y - 0.035,
    linewidth = 0.4
  ) +
  annotate(
    "text",
    x = 1.5,
    y = bracket_y + 0.045,
    label = sprintf("italic(p) == %.3f", independent_test$p.value),
    parse = TRUE,
    size = 2.7
  ) +
  scale_fill_manual(
    values = c(
      "Absent" = "#D0D0D0",
      "Present" = "#7FA6C9"
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = function(x) sprintf("%.1f", x),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Cluster contains a pair with\ncross-mappability score >100",
    y = "Fraction with\nSuSiE credible set",
    tag = "B"
  ) +
  theme_crossmap() +
  theme(
    axis.text.x = element_text(size = 6.3, lineheight = 0.94),
    plot.margin = margin(6, 5, 4, 7)
  )

class_palette <- c(
  eQTL_only = "#4C78A8",
  shared = "#52A89D",
  pcQTL_specific = "#E28E2C"
)
strict_plot_data <- strict_groups[
  abs(threshold - 0.75) < 1e-9 &
    scope %in% c("all_groups", "exclude_known_flagged")
]
stopifnot(
  nrow(strict_plot_data) == 6L,
  sum(strict_plot_data[scope == "all_groups", N]) == 394L
)
strict_plot_data[, scope_label := factor(
  scope,
  levels = c("all_groups", "exclude_known_flagged"),
  labels = c("All clusters", "Cross-mappable\nclusters excluded")
)]
strict_plot_data[, coloc_class := factor(
  coloc_class,
  levels = c("eQTL_only", "shared", "pcQTL_specific")
)]

panel_c <- ggplot(
  strict_plot_data,
  aes(x = scope_label, y = N, fill = coloc_class)
) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.25) +
  geom_text(
    aes(label = N),
    position = position_stack(vjust = 0.5),
    size = 2.15,
    colour = "white"
  ) +
  scale_fill_manual(
    values = class_palette,
    breaks = c("eQTL_only", "shared", "pcQTL_specific"),
    labels = c("eQTL only", "Shared", "pcQTL specific")
  ) +
  scale_y_continuous(
    limits = c(0, 420),
    breaks = c(0, 200, 400),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Signal groups",
    fill = NULL,
    tag = "C"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_crossmap() +
  theme(
    axis.text.x = element_text(size = 6.3, lineheight = 0.94),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 5.8),
    legend.key.size = unit(3, "mm"),
    legend.spacing.x = unit(1, "mm"),
    plot.margin = margin(6, 5, 4, 7)
  )

composite <- arrangeGrob(
  panel_a,
  panel_b,
  panel_c,
  ncol = 3,
  widths = c(1.45, 1, 1.25)
)

width_in <- 178 / 25.4
height_in <- 58 / 25.4
stem <- file.path(plot_dir, "crossmap_paper_reference_composite")

cairo_pdf(
  paste0(stem, ".pdf"),
  width = width_in,
  height = height_in,
  family = "sans"
)
grid.newpage()
grid.draw(composite)
dev.off()

png(
  paste0(stem, ".png"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  type = "cairo",
  bg = "white"
)
grid.newpage()
grid.draw(composite)
dev.off()

svg(
  paste0(stem, ".svg"),
  width = width_in,
  height = height_in,
  family = "sans",
  bg = "white"
)
grid.newpage()
grid.draw(composite)
dev.off()

message(
  "Saved Supplementary Figure S4 A-C; independent Welch t-test p = ",
  signif(independent_test$p.value, 6)
)
