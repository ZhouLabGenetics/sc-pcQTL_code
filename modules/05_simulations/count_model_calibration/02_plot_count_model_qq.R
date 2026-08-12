#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop(
    paste(
      "Usage: Rscript 02_plot_count_model_qq.R",
      "RESULTS_DIR OUTPUT_PDF [SETTINGS_TSV]"
    ),
    call. = FALSE
  )
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
results_dir <- args[[1L]]
output_pdf <- args[[2L]]
settings_file <- if (length(args) == 3L) args[[3L]] else file.path(script_dir, "settings.tsv")
settings <- fread(settings_file)
settings[, score_test_label := ifelse(
  tolower(score_test) %in% c("true", "1", "yes"), "score", "wald"
)]
settings[, facet_label := sprintf(
  "count=%s\nzero=%s\n%s", count_dist, zero_dist, score_test_label
)]

load_pvalues <- function(path) {
  table <- fread(path)
  setnames(table, tolower(names(table)))
  table[!is.na(pvalue) & pvalue > 0 & pvalue <= 1, pvalue]
}

make_qq_table <- function(pvalues, setting_label, component) {
  pvalues <- sort(pvalues)
  n <- length(pvalues)
  data.table(
    setting_label = setting_label,
    component = component,
    Expected = -log10((seq_len(n) - 0.5) / n),
    Observed = -log10(pvalues)
  )
}

qq_table <- rbindlist(lapply(settings$setting_label, function(setting_label) {
  setting_dir <- file.path(results_dir, setting_label)
  count_p <- load_pvalues(file.path(setting_dir, "sc_hurdle_count_pvalues.tsv"))
  zero_p <- load_pvalues(file.path(setting_dir, "sc_hurdle_zero_pvalues.tsv"))
  rbindlist(list(
    make_qq_table(count_p, setting_label, "Count"),
    make_qq_table(zero_p, setting_label, "Zero"),
    make_qq_table(pmin(count_p, zero_p), setting_label, "Hurdle")
  ), use.names = TRUE)
}), use.names = TRUE)

qq_table <- merge(
  qq_table,
  settings[, .(setting_label, facet_label)],
  by = "setting_label",
  all.x = TRUE
)
qq_table[, component := factor(component, levels = c("Count", "Zero", "Hurdle"))]
qq_table[, facet_label := factor(facet_label, levels = settings$facet_label)]

plot <- ggplot(qq_table, aes(Expected, Observed)) +
  geom_point(size = 0.5, alpha = 0.6, color = "#1f78b4") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  facet_grid(component ~ facet_label, scales = "free") +
  labs(
    title = "fasthurdle v1.2.0 score-test null QQ comparison",
    x = "Expected -log10(p)",
    y = "Observed -log10(p)"
  ) +
  theme_bw()

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(output_pdf, plot, width = 12, height = 8)
fwrite(
  qq_table,
  file.path(dirname(output_pdf), "qq_distribution_assumptions_score_test_v12_source.tsv"),
  sep = "\t"
)

message("QQ grid written to ", normalizePath(output_pdf))
