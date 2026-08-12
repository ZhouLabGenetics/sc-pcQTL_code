#!/usr/bin/env Rscript

# Plot-only reproduction of the exploratory MIF--MTFP1 diagnostic.
# Displayed p-values are fixed diagnostic annotations; this script rebuilds the plot only.

lib_dirs <- c(Sys.getenv("SC_PCQTL_R_LIBS"), Sys.getenv("SC_PCQTL_SHARED_R_LIBS"))
lib_dirs <- lib_dirs[nzchar(lib_dirs) & dir.exists(lib_dirs)]
if (length(lib_dirs)) .libPaths(c(lib_dirs, .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

gene_a <- "MIF"
gene_b <- "MTFP1"
pb_p <- 0.9981673
count_p <- 0.02803122
detection_p <- 1.964622e-11
sample_seed <- 20260622L

count_file <- Sys.getenv("COQTL_RAW_COUNTS_FILE")
manuscript_root <- Sys.getenv("SC_PCQTL_MANUSCRIPT_ROOT")
if (!file.exists(count_file)) stop("Set COQTL_RAW_COUNTS_FILE to the B_IN count matrix.")
if (!dir.exists(manuscript_root)) stop("Set SC_PCQTL_MANUSCRIPT_ROOT.")

counts <- fread(count_file, select = c("CellID", gene_a, gene_b))
counts[, donor_id := sub("_.*$", "", CellID)]
if (nrow(counts) != 124968L || uniqueN(counts$donor_id) != 982L) {
  stop("Unexpected exploratory B_IN input dimensions.")
}
pseudobulk <- counts[, lapply(.SD, mean), by = donor_id, .SDcols = c(gene_a, gene_b)]

format_p <- function(x) format(x, digits = 2, scientific = TRUE)
with_marginals <- function(scatter, values, color) {
  x_density <- ggplot(values, aes(x)) +
    geom_density(fill = color, alpha = 0.4, color = NA) + theme_void()
  y_density <- ggplot(values, aes(y)) +
    geom_density(fill = color, alpha = 0.4, color = NA) + coord_flip() + theme_void()
  plot_grid(
    plot_grid(x_density, NULL, rel_widths = c(4, 1)),
    plot_grid(scatter, y_density, rel_widths = c(4, 1)),
    ncol = 1, rel_heights = c(1, 4)
  )
}

pb <- pseudobulk[, .(x = get(gene_a), y = get(gene_b))]
pb_rho <- cor(pb$x, pb$y, method = "spearman")
pb_plot <- ggplot(pb, aes(x, y)) +
  geom_point(alpha = 0.5, color = "#377EB8", size = 1.1) +
  theme_bw(base_size = 13) +
  labs(
    title = sprintf("Pseudobulk: %s vs %s", gene_a, gene_b),
    subtitle = sprintf("Spearman rho = %.2f\nP = %s", pb_rho, format_p(pb_p)),
    x = "", y = ""
  ) +
  theme(plot.title = element_text(size = 15, face = "bold"),
        plot.subtitle = element_text(size = 12), axis.text = element_text(size = 11),
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))

sc <- counts[, .(x = get(gene_a), y = get(gene_b))]
set.seed(sample_seed)
sc <- sc[sample(.N, 5000L)]
sc_rho <- cor(sc$x, sc$y, method = "spearman")
sc_plot <- ggplot(sc[, .N, by = .(x, y)], aes(x, y, size = N)) +
  geom_point(alpha = 0.6, color = "#E41A1C") +
  scale_size_continuous(range = c(1, 5)) +
  theme_bw(base_size = 13) +
  labs(
    title = sprintf("Single-cell: %s vs %s", gene_a, gene_b),
    subtitle = sprintf("Spearman rho = %.2f\nCount P = %s\nDetection P = %s",
                       sc_rho, format_p(count_p), format_p(detection_p)),
    x = "", y = "", size = "Cells"
  ) +
  theme(plot.title = element_text(size = 15, face = "bold"),
        plot.subtitle = element_text(size = 12), axis.text = element_text(size = 11),
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))

output <- file.path(manuscript_root, "figures", "section2", "pair_compare_MIF_MTFP1.pdf")
ggsave(output, plot_grid(with_marginals(pb_plot, pb, "#377EB8"),
                         with_marginals(sc_plot, sc, "#E41A1C"), ncol = 2),
       width = 11, height = 5.8, device = cairo_pdf, bg = "white")
message("Wrote ", output)
