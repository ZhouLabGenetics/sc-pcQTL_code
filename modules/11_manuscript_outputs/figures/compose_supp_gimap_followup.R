#!/usr/bin/env Rscript
# Combined supplementary GIMAP-example figure on ONE full page (single 5-panel composite):
#   a  coloc.susie PP.H4 by cell type          (cell-type specificity)
#   b  GIMAP module completeness by cell type  (cell-type specificity)
#   c  coloc.susie PP.H4 vs lymphocyte count under CD8_NC pseudobulk dilution (single-cell vs bulk)
#   d  Cluster-PC versus gene SMR              (targeted SMR follow-up)
#   e  FinnGen lymphocyte-count GWAS-SMR       (targeted SMR follow-up)
# Final SMR panels d/e share one status palette and legend; the other panels keep
# compact top legends. The theme is harmonized and the panels fill one page. Run
# from the repository root.

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- normalizePath(dirname(args_file), mustWork = TRUE)
source(file.path(script_dir, "main_figure_common.R"))
source(file.path(script_dir, "..", "..", "..", "config", "celltype_eligibility.R"))
primary_ct <- primary_celltypes()
ct <- function(x) canonical_celltype_labels(x)

out_dir <- file.path(OUTPUT_ROOT, "section4")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ===========================================================================
## d, e: targeted SMR follow-up (cluster-PC-vs-gene and FinnGen GWAS-SMR)
## ===========================================================================
cluster_smr_tsv <- Sys.getenv("SC_PCQTL_GIMAP_CLUSTER_SMR_SUMMARY", unset = "")
gwas_smr_tsv <- Sys.getenv("SC_PCQTL_GIMAP_GWAS_SMR_SUMMARY", unset = "")
if (!nzchar(cluster_smr_tsv) || !file.exists(cluster_smr_tsv)) {
  stop("Set SC_PCQTL_GIMAP_CLUSTER_SMR_SUMMARY to the module-03 targeted FinnGen-3019198 output")
}
if (!nzchar(gwas_smr_tsv) || !file.exists(gwas_smr_tsv)) {
  stop("Set SC_PCQTL_GIMAP_GWAS_SMR_SUMMARY to the module-03 targeted FinnGen-3019198 output")
}

smr_levels  <- c("SMR + HEIDI consistent", "SMR strong, HEIDI inconsistent", "Weak/unclear")
smr_palette <- c("SMR + HEIDI consistent"       = "#1B7837",
                 "SMR strong, HEIDI inconsistent" = "#B35806",
                 "Weak/unclear"                   = "#BDBDBD")

fmt_p <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 1))
status_label <- function(x) fcase(
  x == "smr_bonferroni_heidi_pass",         "SMR + HEIDI consistent",
  x == "smr_bonferroni_heidi_inconsistent", "SMR strong, HEIDI inconsistent",
  default = "Weak/unclear"
)

theme_smr <- function(bs = 8.5) theme_panel(bs) +
  theme(axis.text.y = element_text(face = "italic"),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = bs + 0.5, margin = margin(l = 14, b = 2)),
        legend.position = "none")

make_smr_plot <- function(dt, status_col, title, tag) {
  dt <- copy(dt)
  dt[, p_smr_num := suppressWarnings(as.numeric(p_smr))]
  dt[!is.finite(p_smr_num) | p_smr_num <= 0, p_smr_num := .Machine$double.xmin]
  dt[, log10_p := -log10(p_smr_num)]
  dt[, class := factor(status_label(get(status_col)), levels = smr_levels)]
  dt[, gene := factor(gene, levels = dt[order(log10_p)]$gene)]
  dt[, lab := paste0("p=", fmt_p(p_smr_num))]
  bonf <- -log10(0.05 / nrow(dt))
  label_x <- max(dt$log10_p, na.rm = TRUE) + 0.42
  ggplot(dt, aes(log10_p, gene)) +
    geom_segment(aes(x = 0, xend = log10_p, y = gene, yend = gene),
                 linewidth = 0.35, color = "#BDBDBD") +
    geom_vline(xintercept = bonf, linetype = "dashed", linewidth = 0.32, color = "#4B5563") +
    geom_point(aes(fill = class), shape = 21, size = 2.9, color = "#111827", stroke = 0.30) +
    geom_text(aes(x = label_x, label = lab), hjust = 0, size = 2.35, color = "#111827") +
    scale_fill_manual(values = smr_palette, drop = FALSE, name = NULL) +
    scale_x_continuous(name = expression(-log[10](italic(P)[SMR])),
                       expand = expansion(mult = c(0, 0.32))) +
    labs(title = title, tag = tag, y = NULL) +
    theme_smr(8.5)
}

cluster_smr <- fread(cluster_smr_tsv)
gwas_smr <- fread(gwas_smr_tsv)
required_smr <- c("gene", "p_smr", "status")
if (!all(required_smr %in% names(cluster_smr))) {
  stop("Cluster-PC SMR summary lacks final-runner columns: ", paste(setdiff(required_smr, names(cluster_smr)), collapse = ", "))
}
if (!all(required_smr %in% names(gwas_smr))) {
  stop("GWAS SMR summary lacks final-runner columns: ", paste(setdiff(required_smr, names(gwas_smr)), collapse = ", "))
}

p_a <- make_smr_plot(cluster_smr, "status", "Cluster-PC versus gene SMR", "a")
p_b <- make_smr_plot(gwas_smr,    "status", "FinnGen lymphocyte-count GWAS-SMR", "b")

# One shared SMR-status legend for panels a + b.
smr_legend <- extract_legend(
  p_a + guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2.7))) +
    theme(legend.position = "bottom", legend.key.size = unit(0.26, "cm"),
          legend.spacing.x = unit(0.10, "cm"), legend.text = element_text(size = 7.6))
)

## ===========================================================================
## c, d, e: cell-type specificity and pseudobulk/bulk mixing
## ===========================================================================
wf <- file.path(COQTL_WF,
  "03_analysis_celltypes/02_downstream_analysis_modules_add_cov_fdr/gimap_celltype_specificity_analysis")
coloc_event_tsv <- file.path(FORMAL,
  "05_post_coloc_susie_official_finngen_all_finemapped/results/signal_groups/susie_coloc_event_table.tsv")

# Panels a/b are highlight (pass / complete-module) panels, not pcQTL-vs-eQTL
# comparisons, so they use a blue-grey scheme: highlighted cell types in the
# house blue (#4C78A8), the remainder muted grey.
fill2 <- function(true_lab, false_lab) scale_fill_manual(
  values = c(`TRUE` = "#4C78A8", `FALSE` = "#bdbdbd"), breaks = c(TRUE, FALSE),
  labels = c(true_lab, false_lab), na.value = "#bdbdbd", name = NULL)
theme_g <- function(bs = 8.5) theme_panel(bs) +
  theme(legend.position = "top", legend.key.size = unit(0.34, "lines"),
        legend.text = element_text(size = 7.1), legend.margin = margin(0, 0, 0, 0),
        plot.title = element_text(face = "bold", size = bs + 0.5))

## c: formal coloc.susie PP.H4 by cell type
gmap <- unique(fread(file.path(wf, "data/sc_gimap_best_pcqtl_by_celltype.tsv"))[
  celltype %in% primary_ct, .(celltype, cluster_id)
])
ce <- fread(coloc_event_tsv)
ce <- ce[qtl_type == "pcQTL", .(celltype, cluster_id, PP.H4.abf = as.numeric(`PP.H4.abf`),
                                pass = coloc_pass_main %in% c(TRUE, "TRUE", "True"))]
cc <- merge(ce, gmap, by = c("celltype", "cluster_id"))
cm <- cc[, .(pph4 = max(PP.H4.abf, na.rm = TRUE), pass = any(pass)), by = celltype]
cm <- merge(gmap[, .(celltype)], cm, by = "celltype", all.x = TRUE)
cm[!is.finite(pph4), pph4 := 0]; cm[is.na(pass), pass := FALSE]
cm[, lab := factor(ct(celltype), levels = rev(ct(primary_ct)))]
p_c <- ggplot(cm, aes(lab, pph4, fill = pass)) +
  geom_hline(yintercept = 0.75, linetype = "dashed", linewidth = 0.3, color = "grey45") +
  geom_col(width = 0.72, color = "grey30", linewidth = 0.12) + coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  fill2("Colocalized", "Not colocalized") +
  labs(title = "coloc.susie PP.H4 by cell type", x = NULL,
       y = "Max coloc.susie PP.H4 vs FinnGen") +
  theme_g(8.5)

## d: GIMAP module completeness by cell type
m <- fread(file.path(wf, "data/sc_gimap_module_completeness_by_celltype.tsv"))
m <- m[celltype %in% primary_ct,
  .(celltype, n = as.integer(max_gimap_genes_in_any_cluster), full = as.logical(has_full_7_gene_cluster))]
m[, lab := factor(ct(celltype), levels = rev(ct(primary_ct)))]
p_d <- ggplot(m, aes(lab, n, fill = full)) +
  geom_col(width = 0.72, color = "grey30", linewidth = 0.12) + coord_flip() +
  scale_y_continuous(breaks = 1:7, limits = c(0, 7.2), expand = expansion(mult = c(0, 0.02))) +
  fill2("Complete 7-gene module", "Partial") +
  labs(title = "GIMAP module completeness", x = NULL,
       y = "Max GIMAP genes in one cluster") +
  theme_g(8.5)

## c: coloc.susie PP.H4 with lymphocyte count under CD8_NC pseudobulk dilution.
##    Only the OneK1K dilution gradient (0-100% CD8_NC) + the observed realistic
##    composition are plotted. Data frozen by the mixing pipeline's
##    07_make_coloc_susie_publication_summary.R (SuSiE fine-map each mixture QTL with
##    OneK1K in-sample LD, hg19->hg38 liftover, coloc.susie vs FinnGen lymphocyte count 3019198).
cs <- fread(file.path(wf, "mixing/pseudobulk_tensorqtl/coloc_susie/gimap_mixing_coloc_susie_publication_summary.tsv"))
plot_lvls <- c("0%", "5%", "10%", "25%", "50%", "75%", "100%", "Observed")
cs <- cs[comp_label %in% plot_lvls]
csl <- melt(cs[, .(comp_label, pcQTL = pcQTL_pph4_susie, eQTL = eQTL_pph4_susie)],
            id.vars = "comp_label", variable.name = "qtl_type", value.name = "pph4")
csl[, comp_label := factor(comp_label, levels = plot_lvls)]
csl[, qtl_type := factor(qtl_type, levels = c("pcQTL", "eQTL"))]
obs_div <- length(plot_lvls) - 0.5   # dotted divider before the observed-composition reference
p_e <- ggplot(csl, aes(comp_label, pph4, fill = qtl_type)) +
  geom_hline(yintercept = 0.75, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  geom_vline(xintercept = obs_div, linetype = "dotted", linewidth = 0.35, color = "grey55") +
  geom_col(position = position_dodge(width = 0.72), width = 0.66, color = "grey30", linewidth = 0.12) +
  geom_text(data = csl[pph4 >= 0.1], aes(label = sprintf("%.2f", pph4)),
            position = position_dodge(width = 0.72), vjust = -0.3, size = 2.05) +
  scale_fill_manual(values = c(pcQTL = "#F58518", eQTL = "#4C78A8"),
                    labels = c("pcQTL (best cluster PC)", "eQTL (best single gene)"), name = NULL) +
  scale_y_continuous(limits = c(0, 1.03), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "coloc.susie with lymphocyte count under CD8_NC dilution",
       x = "CD8_NC fraction in pseudobulk   (Observed = realistic composition)",
       y = "coloc.susie PP.H4 vs lymphocyte count") +
  theme_g(8.5) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

## ===========================================================================
## Assemble the five-panel supplementary GIMAP example:
##   a = coloc.susie PP.H4 by cell type, b = GIMAP module completeness,
##   c = coloc.susie under CD8_NC pseudobulk dilution,
##   d = cluster-PC-vs-gene targeted SMR, e = FinnGen lymphocyte-count GWAS-SMR.
## ===========================================================================
p_c <- add_tag(p_c, "a"); p_d <- add_tag(p_d, "b"); p_e <- add_tag(p_e, "c")
p_a <- add_tag(p_a, "d"); p_b <- add_tag(p_b, "e")
row_ab <- arrangeGrob(p_c, p_d, ncol = 2)                                   # cell-type specificity
row_de <- arrangeGrob(arrangeGrob(p_a, p_b, ncol = 2), smr_legend,
                      ncol = 1, heights = c(1, 0.16))                        # targeted SMR + shared legend
gimap_supp <- arrangeGrob(row_ab, p_e, row_de, ncol = 1, heights = c(1.0, 1.0, 1.12))
save_composite(gimap_supp, file.path(out_dir, "gimap_example_supp"),
               width = 7.5, height = 8.7)
