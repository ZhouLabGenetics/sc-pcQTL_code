#!/usr/bin/env Rscript
# Supplementary Figure S6: observed-shared-variant coloc.abf sensitivity vs the
# strict coloc.susie signal grouping. By default, draws the three panels in the
# final S7 caption from analysis outputs already produced by the release:
#   (a) class counts at PPH4>0.75 for the coloc.susie graph and the 250 kb ABF run;
#   (b) relative pcQTL-specific increment across PPH4 thresholds (SuSiE vs 250 kb ABF);
#   (c) max eQTL vs pcQTL PP.H4.abf for 250 kb ABF cluster-endpoint groups.
#
# Inputs (all produced upstream in the release; no hidden state):
#   - ABF classes per window: <abf-results-dir>/window_<W>kb/qtl_gwas_abf_raw_coloc_classes.tsv
#     (from scripts 12_coloc_qtl_gwas_abf_raw.R + 12_merge_qtl_gwas_abf_raw_chunks.R)
#   - coloc.susie strict signal groups: module 07 strict_signal_groups.tsv
#     (deposited as Supplementary Data; has `threshold` and `coloc_class`).
# Class labels are shared across both sources: eQTL_only / shared / pcQTL_specific.
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)  # base R; used to arrange the 2x2 panel grid without patchwork
})

script <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1])
root <- normalizePath(file.path(dirname(script), ".."))
source(file.path(root, "config", "config.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

abf_dir      <- get_arg("--abf-results-dir", file.path(ROOT_DIR, "results/coloc/qtl_gwas_abf_raw"))
susie_groups <- get_arg("--susie-groups",
                        file.path(ROOT_DIR, "06_strict_signal_grouping/results/strict_graph/strict_signal_groups.tsv"))
out_pdf      <- get_arg("--out", file.path(ROOT_DIR, "results/plots/abf_susie_sensitivity_comparison.pdf"))
windows      <- as.integer(strsplit(get_arg("--windows", "250"), ",")[[1]])
primary_kb   <- as.integer(get_arg("--primary-window", "250"))
thresholds   <- c(0.70, 0.75, 0.80)
primary_ct   <- names(CELLTYPE_EQTL_MAP)

is_primary_group <- function(x) {
  vapply(strsplit(as.character(x), ",", fixed = TRUE), function(ids) {
    ids <- trimws(ids[nzchar(ids)])
    length(ids) > 0L && all(ids %chin% primary_ct)
  }, logical(1))
}

# Unified colocalization-class palette + order (matches main-text Figure 3).
CLASS_LEVELS <- c("eQTL_only", "shared", "pcQTL_specific")
CLASS_COLS   <- c(eQTL_only = "#4C78A8", shared = "#B279A2", pcQTL_specific = "#F58518")
norm_class <- function(x) {
  x <- gsub("-", "_", tolower(as.character(x)))
  x <- sub("^eqtl_only$", "eQTL_only", x)
  x <- sub("^pcqtl_specific$", "pcQTL_specific", x)
  x <- sub("^shared$", "shared", x)
  factor(x, levels = CLASS_LEVELS)
}

abf_classes <- function(window_kb) {
  f <- file.path(abf_dir, sprintf("window_%skb", window_kb), "qtl_gwas_abf_raw_coloc_classes.tsv")
  if (!file.exists(f)) stop("Missing ABF class table: ", f)
  dt <- fread(f)
  if (!"celltype" %in% names(dt)) stop("ABF class table lacks celltype: ", f)
  dt <- dt[celltype %chin% primary_ct]
  if (!nrow(dt)) stop("No primary-analysis cell types remain in ABF class table: ", f)
  dt[, window_kb := window_kb]
  dt
}

# Threshold -> ABF coloc_class column name (12_merge writes 0.75 as `coloc_class`
# and 0.70/0.80 as `coloc_class_h4_0_7` / `coloc_class_h4_0_8`).
abf_class_col <- function(th) if (isTRUE(all.equal(th, 0.75))) "coloc_class" else
  sprintf("coloc_class_h4_0_%d", as.integer(round(th * 10)))

count_classes <- function(dt, class_col) {
  x <- dt[get(class_col) %in% CLASS_LEVELS, .N, by = .(coloc_class = get(class_col))]
  x[, coloc_class := norm_class(coloc_class)]
  x
}

increment <- function(dt, class_col) {
  n_pc <- dt[get(class_col) == "pcQTL_specific", .N]
  n_eqtl_containing <- dt[get(class_col) %in% c("eQTL_only", "shared"), .N]
  if (n_eqtl_containing == 0L) return(NA_real_)
  n_pc / n_eqtl_containing
}

## ---- Load sources ----------------------------------------------------------
abf_by_window <- rbindlist(lapply(windows, abf_classes), fill = TRUE)
abf_primary <- abf_by_window[window_kb == primary_kb]
if (!file.exists(susie_groups)) stop("Missing coloc.susie strict-group data: ", susie_groups)
susie <- fread(susie_groups)
setnames(susie, tolower(names(susie)))
stopifnot(all(c("threshold", "coloc_class") %in% names(susie)))
if (!"celltypes" %in% names(susie)) stop("Strict SuSiE group table lacks celltypes")
susie <- susie[is_primary_group(celltypes)]
if (!nrow(susie)) stop("No primary-analysis cell types remain in strict SuSiE groups")

## ---- Panel a: class counts at 0.75, SuSiE vs ABF(primary) ------------------
a_su <- count_classes(susie[abs(as.numeric(threshold) - 0.75) < 1e-8], "coloc_class")[, method := "coloc.susie"]
a_ab <- count_classes(abf_primary, "coloc_class")[, method := sprintf("coloc.abf %dkb", primary_kb)]
a_df <- rbind(a_su, a_ab)
a_df[, method := factor(method, levels = c("coloc.susie", sprintf("coloc.abf %dkb", primary_kb)))]
pa <- ggplot(a_df, aes(method, N, fill = coloc_class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72, color = "grey30", linewidth = 0.12) +
  scale_fill_manual(values = CLASS_COLS, drop = FALSE, name = NULL) +
  labs(title = "Class counts (PPH4 > 0.75)", x = NULL, y = "Cluster-endpoint groups", tag = "a") +
  theme_bw(base_size = 9) + theme(legend.position = "top", plot.title = element_text(face = "bold"),
                                  plot.tag = element_text(face = "bold", size = 13))

## ---- Panel b: ABF class counts across windows ------------------------------
b_df <- rbindlist(lapply(windows, function(w) {
  count_classes(abf_by_window[window_kb == w], "coloc_class")[, window := sprintf("%d kb", w)]
}))
b_df[, window := factor(window, levels = sprintf("%d kb", windows))]
pb <- ggplot(b_df, aes(window, N, fill = coloc_class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72, color = "grey30", linewidth = 0.12) +
  scale_fill_manual(values = CLASS_COLS, drop = FALSE, guide = "none") +
  labs(title = "ABF class counts by window", x = NULL, y = "Cluster-endpoint groups", tag = "b") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold"),
                                  plot.tag = element_text(face = "bold", size = 13))

## ---- Panel c: relative pcQTL-specific increment across thresholds ----------
c_df <- rbindlist(lapply(thresholds, function(th) {
  su <- increment(susie[abs(as.numeric(threshold) - th) < 1e-8], "coloc_class")
  ab <- increment(abf_primary, abf_class_col(th))
  data.table(threshold = th,
             value = c(su, ab),
             method = c("coloc.susie", sprintf("coloc.abf %dkb", primary_kb)))
}))
pc <- ggplot(c_df[is.finite(value)], aes(threshold, value, color = method)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_color_manual(values = c("coloc.susie" = "#333333",
                                setNames("#F58518", sprintf("coloc.abf %dkb", primary_kb))), name = NULL) +
  scale_x_continuous(breaks = thresholds) +
  labs(title = "pcQTL-specific increment", x = "PPH4 threshold",
       y = "pcQTL-specific / eQTL-containing", tag = "c") +
  theme_bw(base_size = 9) + theme(legend.position = "top", plot.title = element_text(face = "bold"),
                                  plot.tag = element_text(face = "bold", size = 13))

## ---- Panel d: ABF max eQTL vs pcQTL PP.H4.abf (primary window) --------------
d_df <- copy(abf_primary)
d_df[, class_plot := ifelse(coloc_class %in% CLASS_LEVELS, as.character(coloc_class), "none")]
d_df[, x := fifelse(is.finite(max_eQTL_PPH4), max_eQTL_PPH4, 0)]
d_df[, y := fifelse(is.finite(max_pcQTL_PPH4), max_pcQTL_PPH4, 0)]
pd <- ggplot() +
  geom_point(data = d_df[class_plot == "none"], aes(x, y), color = "grey75", size = 1.1, alpha = 0.6) +
  geom_point(data = d_df[class_plot != "none"],
             aes(x, y, color = celltype, shape = norm_class(class_plot)), size = 1.5, alpha = 0.85) +
  geom_vline(xintercept = 0.75, linetype = 2, color = "grey40") +
  geom_hline(yintercept = 0.75, linetype = 2, color = "grey40") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_shape_manual(values = c(eQTL_only = 15, shared = 17, pcQTL_specific = 16), drop = FALSE, name = "class") +
  labs(title = sprintf("ABF %d kb: max eQTL vs pcQTL PP.H4", primary_kb),
       x = "Max eQTL PP.H4.abf", y = "Max pcQTL PP.H4.abf", color = "cell type", tag = "d") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold"),
                                  plot.tag = element_text(face = "bold", size = 13),
                                  legend.key.size = unit(0.7, "lines"), legend.text = element_text(size = 6.5))

# Arrange the panels with base grid (no patchwork dependency). A single-window
# refresh omits the redundant window-comparison panel.
dir_create(dirname(out_pdf))
if (length(windows) == 1L) {
  pc <- pc + labs(tag = "b")
  pd <- pd + labs(tag = "c")
  grDevices::cairo_pdf(out_pdf, width = 12, height = 4.3)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 3, widths = unit(c(1, 1, 1.15), "null"))))
  print(pa, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pc, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(pd, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
  popViewport()
} else {
  grDevices::cairo_pdf(out_pdf, width = 11, height = 8.5)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))
  cell <- function(r, c) viewport(layout.pos.row = r, layout.pos.col = c)
  print(pa, vp = cell(1, 1)); print(pb, vp = cell(1, 2))
  print(pc, vp = cell(2, 1)); print(pd, vp = cell(2, 2))
  popViewport()
}
invisible(dev.off())
message("Wrote ABF-vs-SuSiE sensitivity comparison to ", out_pdf)
