#!/usr/bin/env Rscript
# Main cluster-landscape and pcQTL-colocalization figure. Five panels:
#   a  detected clusters per cell type            (cluster landscape)
#   b  cluster-size distribution (genes/cluster)  (cluster landscape)
#   c  FDR-significant cluster-PC phenotypes per cell type
#   d  connected-component classes by cell type (PPH4 > 0.75)
#   e  maximum single-gene eQTL vs pcQTL PP.H4 by cluster-endpoint group
# The strict-SuSiE threshold sensitivity panel is exported separately.
# Run from the manuscript repository root.

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- normalizePath(dirname(args_file), mustWork = TRUE)
source(file.path(script_dir, "main_figure_common.R"))
release_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
source(file.path(release_root, "config", "celltype_eligibility.R"))
eligibility_manifest <- file.path(release_root, "config", "celltype_eligibility.tsv")
cell_levels <- primary_celltypes(eligibility_manifest)
celltype_labels <- canonical_celltype_map(eligibility_manifest, primary_only = TRUE)
cell_label_levels <- canonical_celltype_labels(cell_levels, eligibility_manifest)

out_dir <- file.path(OUTPUT_ROOT, "section5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
strict_dir <- file.path(FORMAL, "06_strict_signal_grouping/results/strict_graph")
signal_file <- file.path(FORMAL, "05_post_coloc_susie_official_finngen_all_finemapped/results/signal_groups/susie_signal_groups.tsv")

class_colors <- c(eQTL_only = "#4C78A8", shared = "#72B7B2", pcQTL_specific = "#F58518")
class_levels <- c("eQTL_only", "shared", "pcQTL_specific")
class_labels3 <- c("eQTL-only", "shared", "pcQTL-specific")
cell_cols <- c(b_in="#1b9e77",b_mem="#66a61e",cd4_et="#7570b3",cd4_nc="#386cb0",cd8_et="#e7298a",
               cd8_nc="#e78ac3",cd8_s100b="#a6761d",mono_c="#e6ab02",mono_nc="#a6d854",nk="#1f78b4")
## ---- a, b: cluster landscape -----------------------------------------------
## Unified bar palette: cluster counts = pcQTL orange, size distribution = blue.
land_dir <- Sys.getenv(
  "SC_PCQTL_CLUSTER_LANDSCAPE_DATA",
  unset = if (nzchar(MANUSCRIPT_ROOT)) file.path(MANUSCRIPT_ROOT, "data", "cluster_landscape") else ""
)
if (!nzchar(land_dir) || !dir.exists(land_dir)) {
  stop("Set SC_PCQTL_CLUSTER_LANDSCAPE_DATA or run build_cluster_landscape_tables.R with SC_PCQTL_MANUSCRIPT_ROOT.")
}
land <- fread(file.path(land_dir, "sc_cluster_summary_by_celltype.tsv"))
land[, cell_type_label := unname(celltype_labels[cell_type])]
if (anyNA(land$cell_type_label)) stop("Missing abbreviated cell-type label in Figure 3a.")
land[, cell_type_label := factor(
  cell_type_label,
  levels = land[order(clusters)]$cell_type_label
)]
p_land <- ggplot(land, aes(clusters, cell_type_label)) +
  geom_col(width = 0.72, fill = "#F58518", color = "grey30", linewidth = 0.12) +
  geom_text(aes(label = scales::comma(clusters)), hjust = -0.18, size = 2.5, color = "black") +
  scale_x_continuous(name = "Detected clusters", labels = scales::comma,
                     expand = expansion(mult = c(0, 0.22))) +
  labs(title = "Clusters per cell type", y = NULL) +
  theme_panel(9) + theme(panel.grid.major.y = element_blank())

dist <- fread(file.path(land_dir, "sc_cluster_size_distribution.tsv"))
dist[, size_group := factor(size_group, levels = size_group)]
p_dist <- ggplot(dist, aes(size_group, n_clusters)) +
  geom_col(width = 0.70, fill = "#4C78A8", color = "grey30", linewidth = 0.12) +
  geom_text(aes(label = scales::comma(n_clusters)), vjust = -0.25, size = 2.3, color = "black") +
  scale_y_log10(name = "Clusters", labels = scales::comma, expand = expansion(mult = c(0, 0.28))) +
  scale_x_discrete(name = "Genes per cluster") +
  labs(title = "Cluster size distribution") +
  theme_panel(9) + theme(panel.grid.major.x = element_blank())

## ---- Supplementary: class counts across thresholds ------------------------
groups <- fread(file.path(strict_dir, "strict_signal_groups.tsv"))
groups <- groups[celltypes %chin% cell_levels]
cc <- groups[, .N, by = .(threshold, coloc_class)][coloc_class %in% names(class_colors)]
cc[, threshold_label := paste0("> ", threshold)]
cc[, coloc_class := factor(coloc_class, levels = class_levels)]
p_thr <- ggplot(cc, aes(threshold_label, N, fill = coloc_class)) +
  geom_col(color = "grey25", linewidth = 0.2, width = 0.78) +
  scale_fill_manual(values = class_colors, name = "Strict group class", labels = class_labels3) +
  labs(title = "Across PPH4 thresholds", x = "Colocalization threshold", y = "Connected components") +
  theme_panel(9) +
  theme(legend.position = "top", panel.grid.major.x = element_blank()) +
  guides(fill = guide_legend(nrow = 1))
save_composite(
  p_thr,
  file.path(out_dir, "supp_susie_class_counts_by_threshold"),
  width = 5.2,
  height = 3.2
)

## ---- c: significant cluster-PC phenotypes by cell type --------------------
if (!nzchar(MANUSCRIPT_ROOT)) {
  stop("Set SC_PCQTL_MANUSCRIPT_ROOT to read the verified pcQTL discovery summary.")
}
discovery_file <- file.path(MANUSCRIPT_ROOT, "data", "supp_pcqtl_discovery_by_celltype.tsv")
if (!file.exists(discovery_file)) {
  stop("Missing pcQTL discovery summary; run build_pcqtl_discovery_by_celltype.R: ", discovery_file)
}
discovery <- fread(discovery_file)
required_discovery <- c(
  "celltype",
  "n_successfully_tested_cluster_pc_phenotypes",
  "n_within_phenotype_fdr_hit_cluster_pc_phenotypes"
)
if (length(setdiff(required_discovery, names(discovery)))) {
  stop("pcQTL discovery summary is missing required columns: ",
       paste(setdiff(required_discovery, names(discovery)), collapse = ", "))
}
if (anyDuplicated(discovery$celltype)) stop("Duplicate cell types in ", discovery_file)
if (nrow(discovery) != nrow(land) ||
    !setequal(discovery$celltype, land$cell_type)) {
  stop("Figure 3 cluster and pcQTL discovery summaries do not contain the same cell types.")
}
if (sum(discovery$n_successfully_tested_cluster_pc_phenotypes) != 4353L ||
    sum(discovery$n_within_phenotype_fdr_hit_cluster_pc_phenotypes) != 2040L) {
  stop("pcQTL discovery totals differ from the 10-cell-type analysis contract (4,353 tested; 2,040 significant).")
}
discovery[, celltype_label := unname(celltype_labels[celltype])]
if (anyNA(discovery$celltype_label)) stop("Missing abbreviated cell-type label in Figure 3c.")
discovery[, celltype_plot := factor(
  celltype_label,
  levels = discovery[
    order(n_within_phenotype_fdr_hit_cluster_pc_phenotypes, celltype)
  ]$celltype_label
)]
p_sig <- ggplot(discovery, aes(
  n_within_phenotype_fdr_hit_cluster_pc_phenotypes,
  celltype_plot,
  fill = celltype
)) +
  geom_col(width = 0.72, color = "grey30", linewidth = 0.12) +
  geom_text(
    aes(label = scales::comma(n_within_phenotype_fdr_hit_cluster_pc_phenotypes)),
    hjust = -0.18,
    size = 2.5,
    color = "black"
  ) +
  scale_x_continuous(
    name = "Significant cis-pcQTLs",
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.20))
  ) +
  scale_fill_manual(values = cell_cols, guide = "none", drop = FALSE) +
  labs(y = NULL, title = "Significant cis-pcQTLs") +
  theme_panel(9) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(size = 6.5),
    axis.text.y = element_text(size = 6.5),
    axis.title.x = element_text(size = 8)
  )

## ---- d: class counts by cell type (PPH4 > 0.75) ----------------------------
mg <- groups[abs(threshold - 0.75) < 1e-9 & coloc_class %in% names(class_colors)]
mg[, coloc_class := factor(coloc_class, levels = class_levels)]
ctc <- mg[, .N, by = .(celltypes, coloc_class)]; setnames(ctc, "celltypes", "celltype")
ctc[, celltype := factor(
  celltype,
  levels = ctc[, sum(N), by = celltype][order(-V1)]$celltype
)]
p_ct <- ggplot(ctc, aes(celltype, N, fill = coloc_class)) +
  geom_col(color = "grey25", linewidth = 0.15, width = 0.82) +
  scale_fill_manual(values = class_colors, name = "Strict group class", labels = class_labels3) +
  scale_x_discrete(labels = celltype_labels) +
  labs(title = "Signal groups by cell type (PPH4 > 0.75)", x = "Cell type", y = "Connected components") +
  theme_panel(9) +
  theme(legend.position = "top", panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  guides(fill = guide_legend(nrow = 1))

## ---- e: max eQTL vs pcQTL PP.H4 scatter ------------------------------------
dt <- fread(signal_file)
dt <- dt[celltype %chin% cell_levels]
dt[, `:=`(max_eQTL_PPH4 = fifelse(is.na(as.numeric(max_eQTL_PPH4)), 0, as.numeric(max_eQTL_PPH4)),
          max_pcQTL_PPH4 = fifelse(is.na(as.numeric(max_pcQTL_PPH4)), 0, as.numeric(max_pcQTL_PPH4)))]
class_shapes <- c("pcQTL-specific" = 24, "shared" = 16, "eQTL-only" = 15)
# Non-colocalized groups (both PPH4 < 0.75) are kept as a light grey background
# layer so their presence is visible without competing with the colocalized
# points; subsampled with a fixed seed because there are ~69k of them.
bg <- dt[!coloc_class_main %chin% c("eQTL_only", "pcQTL_specific", "shared")]
set.seed(1); if (nrow(bg) > 8000) bg <- bg[sample(.N, 8000)]
fg <- dt[coloc_class_main %chin% c("eQTL_only", "pcQTL_specific", "shared")]
fg[, celltype := factor(celltype, levels = cell_levels)]
fg[, coloc_class_main := factor(coloc_class_main, levels = c("pcQTL_specific", "shared", "eQTL_only"),
                                labels = c("pcQTL-specific", "shared", "eQTL-only"))]
p_sc <- ggplot(mapping = aes(max_eQTL_PPH4, max_pcQTL_PPH4)) +
  geom_point(data = bg, color = "grey80", size = 0.65, alpha = 0.4, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, color = "#9e9e9e", linewidth = 0.35) +
  geom_vline(xintercept = 0.75, linetype = "dashed", color = "#333333", linewidth = 0.35) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "#333333", linewidth = 0.35) +
  geom_point(data = fg, aes(color = celltype, shape = coloc_class_main), size = 1.7, alpha = 0.9, stroke = 0.25) +
  scale_color_manual(values = cell_cols, name = "Cell type", labels = celltype_labels, drop = FALSE) +
  scale_shape_manual(values = class_shapes, name = "SuSiE class", drop = FALSE) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.02, 0.02))) +
  coord_fixed(ratio = 1) +
  labs(title = "Max eQTL vs pcQTL PP.H4",
       x = "Max eQTL PP.H4", y = "Max pcQTL PP.H4") +
  theme_panel(9) +
  theme(legend.position = "right", legend.key.size = unit(0.36, "lines"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 7),
        legend.spacing.y = unit(1, "pt")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 1.5)),
         shape = guide_legend(override.aes = list(size = 1.7)))

## ---- assemble (two columns) ------------------------------------------------
## Preserve the original Figure 3 arrangement: panels a-c stacked in the left
## column and panels d-e stacked in the right column. Panel c alone is replaced
## by the significant cis-pcQTL count plot.
p_land <- add_tag(p_land, "a"); p_dist <- add_tag(p_dist, "b")
p_sig  <- add_tag(p_sig,  "c"); p_ct   <- add_tag(p_ct,   "d")
p_sc   <- add_tag(p_sc,   "e")
left_col  <- arrangeGrob(p_land, p_dist, p_sig, ncol = 1, heights = c(1.35, 0.7, 1.4))
right_col <- arrangeGrob(p_ct, p_sc, ncol = 1, heights = c(2.05, 1.4))
figure <- arrangeGrob(left_col, right_col, ncol = 2, widths = c(0.34, 0.66))
save_composite(figure, file.path(out_dir, "main_pcqtl_results"), width = 7.2, height = 7.7)
