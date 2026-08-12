#!/usr/bin/env Rscript
# Shared helpers for composing dense multi-panel MAIN-TEXT figures.
#
# Goal: recompose several low-density single/two-panel main figures into denser
# a/b/c/d composites, and promote a few informative supplementary panels into the
# main text. Panels are REDRAWN from their frozen source tables (content
# unchanged) with a consistent theme and a controlled font-to-plot ratio so the
# composite stays legible at \textwidth. Composition uses gridExtra::arrangeGrob,
# mirroring the existing shared main-figure convention.
#
# Source this from a per-figure compose_*.R run at the manuscript repository root.

suppressWarnings({
  .lib_candidates <- c(
    Sys.getenv("SC_PCQTL_R_LIBS", unset = ""),
    Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = "")
  )
  .lib_candidates <- .lib_candidates[dir.exists(.lib_candidates)]
  if (length(.lib_candidates)) .libPaths(c(.lib_candidates, .libPaths()))
})

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(scales)
})

## Repository roots ----------------------------------------------------------
.common_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) "")
.figure_dir <- if (nzchar(.common_file)) dirname(.common_file) else getwd()
COQTL_WF <- Sys.getenv("COQTL_WORKFLOW_ROOT", unset = "")
if (!nzchar(COQTL_WF)) stop("Set COQTL_WORKFLOW_ROOT to the analysis workflow root.")
COQTL_WF <- normalizePath(COQTL_WF, mustWork = TRUE)
FORMAL <- normalizePath(
  Sys.getenv("SC_PCQTL_FORMAL_COLOC_ROOT", unset = file.path(COQTL_WF, "04_formal_colocalization")),
  mustWork = TRUE
)
.output_override <- Sys.getenv("SC_PCQTL_FIGURE_OUTPUT_ROOT", unset = "")
.manuscript_root <- Sys.getenv("SC_PCQTL_MANUSCRIPT_ROOT", unset = "")
MANUSCRIPT_ROOT <- .manuscript_root
.output_default <- if (nzchar(.manuscript_root)) {
  file.path(.manuscript_root, "figures")
} else {
  ""
}
if (!nzchar(.output_override) && !nzchar(.output_default)) {
  stop("Set SC_PCQTL_MANUSCRIPT_ROOT or SC_PCQTL_FIGURE_OUTPUT_ROOT.")
}
OUTPUT_ROOT <- normalizePath(if (nzchar(.output_override)) .output_override else .output_default,
                             mustWork = FALSE)
dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

## Shared theme for composite sub-panels ------------------------------------
## Keep text readable when the composite is rendered near manuscript width.
theme_panel <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text  = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      plot.subtitle = element_text(hjust = 0, size = base_size - 1, color = "grey25"),
      legend.title = element_text(face = "bold", size = base_size - 0.5),
      legend.text  = element_text(size = base_size - 1),
      legend.key.size = unit(0.36, "lines"),
      plot.tag = element_text(face = "bold", size = base_size + 4),
      plot.tag.position = "topleft",
      plot.margin = margin(5, 6, 3, 4)
    )
}

## Add a bold panel tag (a, b, c, ...) anchored at the top-left.
add_tag <- function(p, tag) {
  p + labs(tag = tag)
}

## Extract a single guide-box grob from a plot (for a shared legend). Robust to
## ggplot >= 3.5 guide-box naming.
extract_legend <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(grepl("guide-box", vapply(g$grobs, function(x) x$name, character(1))))
  for (i in idx) if (!inherits(g$grobs[[i]], "zeroGrob")) return(g$grobs[[i]])
  if (length(idx)) return(g$grobs[[idx[1]]])
  NULL
}

## Convert ggplots to grobs with equalized widths so stacked panels align.
align_widths <- function(plots) {
  grobs <- lapply(plots, ggplotGrob)
  maxw <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$widths))
  lapply(grobs, function(g) { g$widths <- maxw; g })
}

## Save a composite grob to PDF (cairo_pdf) + PNG preview.
save_composite <- function(grob, out_basename, width, height) {
  pdf_path <- paste0(out_basename, ".pdf")
  png_path <- paste0(out_basename, ".png")
  ggsave(pdf_path, grob, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(png_path, grob, width = width, height = height, dpi = 320, bg = "white")
  message("Wrote ", pdf_path)
  message("Wrote ", png_path)
  invisible(pdf_path)
}
