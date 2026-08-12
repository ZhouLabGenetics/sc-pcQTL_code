#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(script_file), "common.R"))
prepend_project_libraries()
args <- parse_cli()
opt <- list(
  global_dir = require_arg(args, "global_dir"),
  within_dir = require_arg(args, "within_dir"),
  output_dir = require_arg(args, "output_dir"),
  p_thresholds = if (is.null(args$p_thresholds)) {
    "0.05,0.01,0.005,0.001"
  } else {
    args$p_thresholds
  }
)
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(opt$output_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
thresholds <- as.numeric(strsplit(opt$p_thresholds, ",", fixed = TRUE)[[1L]])

read_p <- function(path) {
  if (!file.exists(path)) stop("Missing p-value file: ", path)
  x <- read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                  check.names = FALSE, na.strings = c("NA", ""))
  if (!("Pvalue" %in% names(x))) stop("Pvalue column missing from ", path)
  x[, intersect(c("Gene1", "Gene2", "Pvalue"), names(x)), drop = FALSE]
}

load_scheme <- function(label, count_file, detection_file, pb_file) {
  count <- read_p(count_file)
  detection <- read_p(detection_file)
  pb <- read_p(pb_file)
  keys <- c("Gene1", "Gene2")
  components <- merge(
    count, detection, by = keys, suffixes = c("_count", "_detection"), all = TRUE
  )
  list(label = label, count = count, detection = detection, pb = pb,
       components = components)
}

global <- load_scheme(
  "Global permutation",
  file.path(opt$global_dir, "sc_hurdle_count_pvalues_pairlevel.tsv"),
  file.path(opt$global_dir, "sc_hurdle_zero_pvalues_pairlevel.tsv"),
  file.path(opt$global_dir, "pb_spearman_pvalues.tsv")
)
within <- load_scheme(
  "Within-donor permutation",
  file.path(opt$within_dir, "sc_hurdle", "sc_hurdle_count_pvalues.tsv"),
  file.path(opt$within_dir, "sc_hurdle", "sc_hurdle_zero_pvalues.tsv"),
  file.path(opt$within_dir, "pb_spearman", "pb_spearman_pvalues.tsv")
)

qq_frame <- function(pvalues, permutation, method) {
  pvalues <- sort(pvalues[is.finite(pvalues) & pvalues >= 0 & pvalues <= 1])
  n <- length(pvalues)
  if (!n) stop("No finite p-values for ", permutation, " / ", method)
  data.frame(
    permutation = permutation,
    method = method,
    expected = -log10((seq_len(n) - 0.5) / n),
    observed = -log10(pmax(pvalues, .Machine$double.xmin)),
    raw_pvalue = pvalues,
    stringsAsFactors = FALSE
  )
}

qq_source <- do.call(rbind, list(
  qq_frame(global$count$Pvalue, global$label, "Count component"),
  qq_frame(global$detection$Pvalue, global$label, "Detection component"),
  qq_frame(global$pb$Pvalue, global$label, "Pseudobulk Spearman"),
  qq_frame(within$count$Pvalue, within$label, "Count component"),
  qq_frame(within$detection$Pvalue, within$label, "Detection component"),
  qq_frame(within$pb$Pvalue, within$label, "Pseudobulk Spearman")
))

summarize_scheme <- function(x) {
  do.call(rbind, lapply(thresholds, function(threshold) {
    count_flag <- x$components$Pvalue_count < threshold
    detection_flag <- x$components$Pvalue_detection < threshold
    data.frame(
      permutation = x$label,
      method = c(
        "Count component", "Detection component", "Component union",
        "Pseudobulk Spearman"
      ),
      p_threshold = threshold,
      positive_fraction = c(
        safe_fraction(count_flag),
        safe_fraction(detection_flag),
        safe_fraction(count_flag | detection_flag),
        safe_fraction(x$pb$Pvalue < threshold)
      ),
      stringsAsFactors = FALSE
    )
  }))
}
fraction_summary <- rbind(summarize_scheme(global), summarize_scheme(within))
atomic_write_table(qq_source, file.path(opt$output_dir, "qq_source_global_vs_within.tsv.gz"))
atomic_write_table(
  fraction_summary,
  file.path(opt$output_dir, "threshold_fraction_comparison_global_vs_within.tsv")
)

panel_order <- data.frame(
  permutation = c(rep(global$label, 3L), rep(within$label, 3L)),
  method = rep(c("Count component", "Detection component", "Pseudobulk Spearman"), 2L),
  panel_tag = letters[1:6],
  stringsAsFactors = FALSE
)
panel_summary <- do.call(rbind, lapply(seq_len(nrow(panel_order)), function(i) {
  panel <- panel_order[i, ]
  d <- qq_source[
    qq_source$permutation == panel$permutation & qq_source$method == panel$method,
    , drop = FALSE
  ]
  threshold_row <- fraction_summary[
    fraction_summary$permutation == panel$permutation &
      fraction_summary$method == panel$method &
      abs(fraction_summary$p_threshold - 0.01) < 1e-12,
    , drop = FALSE
  ]
  data.frame(
    permutation = panel$permutation,
    method = panel$method,
    panel_tag = panel$panel_tag,
    n_finite_tests = nrow(d),
    max_expected_minus_log10_p = max(d$expected),
    max_observed_minus_log10_p = max(d$observed),
    fraction_p_lt_0.01 = threshold_row$positive_fraction,
    stringsAsFactors = FALSE
  )
}))
atomic_write_table(panel_summary, file.path(opt$output_dir, "qq_panel_summary.tsv"))

method_colors <- c(
  "Count component" = "#1F78B4",
  "Detection component" = "#33A02C",
  "Pseudobulk Spearman" = "#6A3D9A"
)

draw_row_strip <- function(label) {
  par(mar = c(0.2, 0.1, 0.2, 0.1), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rect(0, 0, 1, 1, col = "#F2F2F2", border = "#D0D0D0", lwd = 0.45)
  text(0.5, 0.5, label, srt = 90, font = 2, cex = 0.68, col = "#222222")
}

draw_panel <- function(permutation, method, panel_tag, top_row) {
  d <- qq_source[
    qq_source$permutation == permutation & qq_source$method == method,
    , drop = FALSE
  ]
  s <- panel_summary[
    panel_summary$permutation == permutation & panel_summary$method == method,
    , drop = FALSE
  ]
  x_max <- max(d$expected) * 1.025
  y_ticks <- pretty(c(0, max(d$observed) * 1.02), n = 5)
  y_ticks <- y_ticks[y_ticks >= 0]
  y_max <- max(y_ticks)
  par(
    mar = c(2.3, 2.75, if (top_row) 1.75 else 0.9, 0.55),
    mgp = c(1.55, 0.42, 0), tcl = -0.22, xaxs = "i", yaxs = "i"
  )
  plot.new()
  plot.window(xlim = c(0, x_max), ylim = c(0, y_max))
  segments(0, 0, max(d$expected), max(d$expected),
           col = "#8C8C8C", lty = 2, lwd = 0.65)
  points(
    d$expected, d$observed, pch = 16, cex = 0.23,
    col = grDevices::adjustcolor(method_colors[[method]], alpha.f = 0.56)
  )
  axis(1, at = 0:4, labels = 0:4, cex.axis = 0.70, lwd = 0.45, lwd.ticks = 0.45)
  axis(2, at = y_ticks, labels = format(y_ticks, trim = TRUE), las = 1,
       cex.axis = 0.70, lwd = 0.45, lwd.ticks = 0.45)
  box(lwd = 0.5, col = "#222222")
  mtext(panel_tag, side = 3, line = 0.25, adj = -0.075, font = 2, cex = 0.90)
  if (top_row) mtext(method, side = 3, line = 0.25, adj = 0.52, font = 2, cex = 0.78)
  text(0.11, y_max * 0.94,
       sprintf("n = %s", format(s$n_finite_tests, big.mark = ",", scientific = FALSE)),
       adj = c(0, 1), cex = 0.61, col = "#333333")
  threshold_label <- sprintf("%.2f%%", 100 * s$fraction_p_lt_0.01)
  text(0.11, y_max * 0.865,
       bquote(Pr(italic(P) < 0.01) == .(threshold_label)),
       adj = c(0, 1), cex = 0.61, col = "#333333")
}

draw_grid <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(
    oma = c(2.8, 3.15, 0.55, 0.25), family = "sans", fg = "#222222",
    col.axis = "#222222", col.lab = "#222222"
  )
  layout(matrix(seq_len(8L), nrow = 2L, byrow = TRUE),
         widths = c(0.13, 1, 1, 1), heights = c(1, 1))
  par(cex = 0.82)
  draw_row_strip(global$label)
  for (i in 1:3) {
    draw_panel(panel_order$permutation[[i]], panel_order$method[[i]],
               panel_order$panel_tag[[i]], TRUE)
  }
  draw_row_strip(within$label)
  for (i in 4:6) {
    draw_panel(panel_order$permutation[[i]], panel_order$method[[i]],
               panel_order$panel_tag[[i]], FALSE)
  }
  mtext(expression(Expected~~-log[10](italic(P))), side = 1, outer = TRUE,
        line = 1.35, cex = 0.76)
  mtext(expression(Observed~~-log[10](italic(P))), side = 2, outer = TRUE,
        line = 1.95, cex = 0.76)
}

width_mm <- 183
height_mm <- 125
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
render_atomic <- function(path, open_device) {
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  on.exit({
    if (dev.cur() > 1L) try(dev.off(), silent = TRUE)
    if (file.exists(tmp)) unlink(tmp)
  }, add = TRUE)
  open_device(tmp)
  draw_grid()
  dev.off()
  if (!file.rename(tmp, path)) stop("Failed to move rendered figure: ", path)
}

figure_base <- file.path(plot_dir, "qq_global_vs_within_donor")
render_atomic(paste0(figure_base, ".pdf"), function(path) {
  grDevices::cairo_pdf(path, width = width_in, height = height_in, family = "sans")
})
render_atomic(paste0(figure_base, ".png"), function(path) {
  grDevices::png(path, width = width_in, height = height_in, units = "in", res = 450,
                 type = "cairo", antialias = "subpixel")
})
render_atomic(paste0(figure_base, ".tiff"), function(path) {
  grDevices::tiff(path, width = width_in, height = height_in, units = "in", res = 600,
                  compression = "lzw", type = "cairo", antialias = "subpixel")
})

message("QQ diagnostic written to ", opt$output_dir)
