#!/usr/bin/env Rscript
# Supplementary gene-effect figure for strict signal groups.
#   a  max |PIP-weighted nominal effect| by strict class (with BH-Wilcoxon brackets)
#   b  cross-gene effect CV by strict class (with BH-Wilcoxon brackets)
# The previously explored joint scatter and regulatory-annotation panel are not
# part of this final figure.

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- normalizePath(dirname(args_file), mustWork = TRUE)
source(file.path(script_dir, "main_figure_common.R"))
release_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
source(file.path(release_root, "config", "celltype_eligibility.R"))
primary_ct <- primary_celltypes(file.path(release_root, "config", "celltype_eligibility.tsv"))

out_dir <- file.path(OUTPUT_ROOT, "section5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ge_dir  <- file.path(FORMAL, "09_mechanistic_celltype_analysis/results/gene_effects")

class_levels <- c("pcQTL_specific", "shared", "eQTL_only", "non_coloc")
palette <- c(pcQTL_specific = "#F58518", shared = "#72B7B2", eQTL_only = "#4C78A8", non_coloc = "#9A9A9A")

## ---- a/b: gene-effect violins with BH-Wilcoxon brackets --------------------
effect_summary <- fread(file.path(ge_dir, "strict_qtl_cs_nominal_effect_summary.tsv"))
effect_summary <- effect_summary[celltype %chin% primary_ct]
if (!nrow(effect_summary)) stop("No primary-analysis gene-effect summaries remain")
test_input <- effect_summary[
  strict_coloc_class %in% class_levels & is.finite(max_abs_effect) & is.finite(cv_effect)
]
plot_dt <- copy(test_input)
plot_dt[, strict_coloc_class := factor(strict_coloc_class, levels = class_levels)]

test_pairs <- list(c("pcQTL_specific", "non_coloc"),
                   c("pcQTL_specific", "eQTL_only"),
                   c("shared", "non_coloc"))
tests <- rbindlist(lapply(test_pairs, function(pair) {
  x <- test_input[strict_coloc_class %in% pair]
  if (length(unique(x$strict_coloc_class)) < 2L) return(NULL)
  data.table(
    comparison = paste(pair, collapse = "_vs_"),
    metric = c("max_abs_effect", "cv_effect"),
    wilcox_p = c(wilcox.test(max_abs_effect ~ strict_coloc_class, data = x)$p.value,
                 wilcox.test(cv_effect ~ strict_coloc_class, data = x)$p.value)
  )
}), fill = TRUE)
if (nrow(tests)) tests[, fdr := p.adjust(wilcox_p, method = "BH"), by = metric]

format_p <- function(p) ifelse(is.na(p), "NA",
  ifelse(p < 1e-4, formatC(p, format = "e", digits = 1), paste0("p=", signif(p, 2))))

class_counts <- plot_dt[, .N, by = strict_coloc_class]
class_n <- setNames(class_counts$N, as.character(class_counts$strict_coloc_class))
class_labels <- setNames(sprintf("%s\nn=%s", class_levels, format(class_n[class_levels], big.mark = ",")),
                         class_levels)

make_annotation <- function(metric_col, log_scale = FALSE) {
  if (!nrow(tests) || !"fdr" %in% names(tests)) return(data.table())
  ann <- tests[metric == metric_col & is.finite(fdr)]
  if (!nrow(ann)) return(ann)
  parts <- tstrsplit(ann$comparison, "_vs_", fixed = TRUE)
  ann[, `:=`(group1 = parts[[1]], group2 = parts[[2]])]
  ann[, `:=`(x = match(group1, class_levels), xend = match(group2, class_levels))]
  ann <- ann[!is.na(x) & !is.na(xend)]
  if (!nrow(ann)) return(ann)
  ann[, p_label := paste0("BH ", format_p(fdr))]
  ann[, span_w := abs(xend - x)]; setorder(ann, span_w)
  vals <- plot_dt[[metric_col]]; vals <- vals[is.finite(vals) & vals > 0]
  if (log_scale) {
    top <- 10 ^ quantile(log10(vals), 0.995, na.rm = TRUE)
    mult <- c(1.25, 1.68, 2.25, 3.0)
    ann[, y := as.numeric(top) * mult[seq_len(.N)]]; ann[, y_text := y * 1.07]
  } else {
    top <- quantile(vals, 0.995, na.rm = TRUE)
    span <- diff(range(vals, na.rm = TRUE)); if (!is.finite(span) || span == 0) span <- top
    ann[, y := as.numeric(top) + span * c(0.08, 0.21, 0.35, 0.5)[seq_len(.N)]]; ann[, y_text := y + span * 0.03]
  }
  ann
}
add_pairwise_annotations <- function(p, ann) {
  p +
    geom_segment(data = ann, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.32) +
    geom_segment(data = ann, aes(x = x, xend = x, y = y, yend = y * 0.985), inherit.aes = FALSE, linewidth = 0.32) +
    geom_segment(data = ann, aes(x = xend, xend = xend, y = y, yend = y * 0.985), inherit.aes = FALSE, linewidth = 0.32) +
    geom_text(data = ann, aes(x = (x + xend) / 2, y = y_text, label = p_label), inherit.aes = FALSE, size = 2.5)
}

violin_base <- function(metric_col) {
  ggplot(plot_dt, aes(x = strict_coloc_class, y = get(metric_col), fill = strict_coloc_class)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.55, na.rm = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9, linewidth = 0.3, na.rm = TRUE) +
    geom_point(position = position_jitter(width = 0.08, height = 0), alpha = 0.12, size = 0.25, na.rm = TRUE) +
    scale_fill_manual(values = palette, drop = FALSE) +
    scale_x_discrete(labels = class_labels, drop = FALSE) +
    theme_panel(9) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 22, hjust = 1, size = 7.5))
}
ann_max <- make_annotation("max_abs_effect", log_scale = TRUE)
p_a <- violin_base("max_abs_effect") + scale_y_continuous(trans = "log10") +
  labs(title = "Maximum gene-level effect", x = NULL, y = "Max |PIP-weighted nominal beta|")
if (nrow(ann_max)) p_a <- add_pairwise_annotations(p_a, ann_max)

ann_cv <- make_annotation("cv_effect", log_scale = FALSE)
p_b <- violin_base("cv_effect") +
  labs(title = "Cross-gene effect dispersion", x = NULL,
       y = "CV of |PIP-weighted nominal beta|")
if (nrow(ann_cv)) p_b <- add_pairwise_annotations(p_b, ann_cv)

## ---- assemble --------------------------------------------------------------
p_a <- add_tag(p_a, "a"); p_b <- add_tag(p_b, "b")
figure <- arrangeGrob(p_a, p_b, ncol = 2)
save_composite(figure, file.path(out_dir, "main_mechanism_composite"),
               width = 7.0, height = 3.4)
