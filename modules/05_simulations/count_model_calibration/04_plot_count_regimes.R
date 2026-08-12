#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: Rscript 04_plot_count_regimes.R RESULTS_DIR OUTPUT_PDF",
    call. = FALSE
  )
}

results_dir <- args[[1L]]
output_pdf <- args[[2L]]
calibration_file <- file.path(results_dir, "controlled_count_regime_calibration.tsv")
if (!file.exists(calibration_file)) stop("Missing calibration table: ", calibration_file)
replicate_file <- file.path(results_dir, "controlled_count_regime_replicates.tsv")
if (!file.exists(replicate_file)) stop("Missing replicate table: ", replicate_file)

calibration <- fread(calibration_file)
calibration[, regime := factor(
  regime,
  levels = c(
    "poisson_generated",
    "shifted_poisson_generated",
    "negative_binomial_generated"
  ),
  labels = c(
    "Zero-truncated Poisson",
    "Shifted Poisson",
    "Negative-Binomial"
  )
)]
calibration[, model := factor(model, levels = c("Poisson", "Negative-Binomial"))]

palette <- c("Poisson" = "#0072B2", "Negative-Binomial" = "#D55E00")

plot <- ggplot(
  calibration,
  aes(nominal_alpha, rejection_ratio, color = model, shape = model, group = model)
) +
  geom_hline(yintercept = 1, linewidth = 0.45, linetype = "dashed", color = "grey45") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2, stroke = 0.7) +
  facet_wrap(~regime, nrow = 1) +
  scale_x_log10(
    breaks = c(0.001, 0.005, 0.01, 0.05),
    labels = c("0.001", "0.005", "0.01", "0.05")
  ) +
  scale_color_manual(values = palette) +
  scale_shape_manual(values = c("Poisson" = 16, "Negative-Binomial" = 17)) +
  labs(
    x = "Nominal alpha",
    y = "Empirical rejection rate / nominal alpha",
    color = "Fitted count model",
    shape = "Fitted count model"
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8, color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8),
    panel.spacing.x = grid::unit(1.1, "lines"),
    plot.margin = margin(4, 6, 4, 4)
  )

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(
  output_pdf, plot, width = 7.6, height = 3.0, units = "in",
  device = grDevices::pdf, useDingbats = FALSE
)
ggsave(sub("\\.pdf$", ".png", output_pdf), plot, width = 7.6, height = 3.0, units = "in", dpi = 300)

replicates <- fread(replicate_file)
theta <- replicates[model == "Negative-Binomial" & is.finite(theta) & theta > 0]
theta[, regime := factor(
  regime,
  levels = c(
    "poisson_generated",
    "shifted_poisson_generated",
    "negative_binomial_generated"
  ),
  labels = c("ZT Poisson", "Shifted Poisson", "Negative-Binomial")
)]

theta_plot <- ggplot(theta, aes(regime, theta, fill = regime)) +
  geom_hline(yintercept = 2, linewidth = 0.45, linetype = "dashed", color = "grey45") +
  geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.55) +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000),
    labels = scales::label_log()
  ) +
  scale_fill_manual(
    values = c(
      "ZT Poisson" = "#56B4E9",
      "Shifted Poisson" = "#009E73",
      "Negative-Binomial" = "#E69F00"
    ),
    guide = "none"
  ) +
  coord_cartesian(ylim = c(1, 1000000)) +
  labs(
    x = "Generating count distribution",
    y = expression("Estimated Negative-Binomial dispersion " * hat(theta))
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8, color = "black"),
    plot.margin = margin(4, 6, 4, 4)
  )

theta_pdf <- sub("\\.pdf$", "_theta.pdf", output_pdf)
ggsave(
  theta_pdf, theta_plot, width = 3.6, height = 3.0, units = "in",
  device = grDevices::pdf, useDingbats = FALSE
)
ggsave(sub("\\.pdf$", ".png", theta_pdf), theta_plot, width = 3.6, height = 3.0, units = "in", dpi = 300)

message(
  "Controlled count-regime figures written to ",
  normalizePath(output_pdf), " and ", normalizePath(theta_pdf)
)
