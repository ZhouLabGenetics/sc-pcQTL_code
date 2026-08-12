#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
})

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (is.na(args_file) || !nzchar(args_file)) args_file <- getwd()
module06_code_dir <- normalizePath(file.path(dirname(args_file), "..", "06_finngen_susie_coloc"), mustWork = TRUE)
source(file.path(module06_code_dir, "config", "config.R"))
root_dir <- ROOT_DIR

strict_dir <- file.path(root_dir, "06_strict_signal_grouping/results/strict_graph")
post_dir <- file.path(root_dir, "05_post_coloc_susie_official_finngen_all_finemapped")
out_dir <- file.path(root_dir, "09_mechanistic_celltype_analysis/results/gene_effects")
dir_create(out_dir)

threshold_main <- 0.75

read_loadings <- function(celltype, cluster_id, pc_id) {
  f <- file.path(post_dir, "cache", "pc_loadings", celltype, cluster_id, "gene_loadings.tsv")
  if (!file.exists(f)) return(NA_character_)
  x <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x) || !"Gene" %in% names(x) || !pc_id %in% names(x)) return(NA_character_)
  x[, loading_abs := abs(get(pc_id))]
  x <- x[order(-loading_abs)]
  n <- min(5L, nrow(x))
  paste(sprintf("%s(%.3f)", x$Gene[seq_len(n)], x[[pc_id]][seq_len(n)]), collapse = ";")
}

cs_files <- list.files(
  file.path(root_dir, "results/fine_mapping/qtl"),
  pattern = "\\.credible_sets\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(cs_files)) stop("No QTL credible set TSV files found")

message("Reading QTL credible set maps")
cs_map <- rbindlist(lapply(cs_files, function(f) {
  x <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  req <- c("celltype", "cluster_id", "phenotype_type", "phenotype_id", "cs_index", "credible_set_id")
  if (!all(req %in% names(x))) return(NULL)
  unique(x[, ..req])
}), fill = TRUE)
cs_map <- cs_map[celltype %chin% names(CELLTYPE_EQTL_MAP)]
cs_map <- unique(cs_map)
cs_map[, qtl_type := fifelse(phenotype_type == "pcQTL", "pcQTL", "eQTL")]
cs_map[, signal_index := as.integer(cs_index)]
cs_map <- unique(cs_map[, .(celltype, cluster_id, qtl_type, phenotype_id, signal_index, credible_set_id)])

message("Reading strict graph nodes and groups")
nodes <- fread(file.path(strict_dir, "strict_graph_nodes.tsv"))
groups <- fread(file.path(strict_dir, "strict_signal_groups.tsv"))
groups <- groups[abs(threshold - threshold_main) < 1e-9 & celltypes %chin% names(CELLTYPE_EQTL_MAP)]
nodes <- nodes[abs(threshold - threshold_main) < 1e-9 &
                 celltype %chin% names(CELLTYPE_EQTL_MAP) & node_type %in% c("eQTL", "pcQTL")]
nodes[, qtl_type := node_type]
nodes[, signal_index := as.integer(signal_index)]

qtl_group_membership <- merge(
  nodes[, .(signal_group_id, node_id, celltype, cluster_id, qtl_type, phenotype_id, signal_index)],
  groups[, .(signal_group_id, coloc_class, n_gwas, n_eqtl, n_pcqtl, phenocodes, gwas_lead_variants)],
  by = "signal_group_id",
  all.x = TRUE
)
qtl_group_membership <- merge(
  qtl_group_membership,
  cs_map,
  by = c("celltype", "cluster_id", "qtl_type", "phenotype_id", "signal_index"),
  all.x = TRUE
)
fwrite(qtl_group_membership, file.path(out_dir, "strict_qtl_cs_group_membership.tsv"), sep = "\t")

message("Reading existing PIP-weighted nominal gene effects")
gene_effects <- fread(file.path(post_dir, "results/gene_effects/pip_weighted_nominal_gene_effects.tsv"))
cs_summary <- fread(file.path(post_dir, "results/gene_effects/qtl_cs_nominal_effect_summary.tsv"))
gene_effects <- gene_effects[celltype %chin% names(CELLTYPE_EQTL_MAP)]
cs_summary <- cs_summary[celltype %chin% names(CELLTYPE_EQTL_MAP)]

strict_cs_class <- qtl_group_membership[!is.na(credible_set_id), .(
  strict_group_ids = paste(sort(unique(signal_group_id)), collapse = ";"),
  strict_coloc_classes = paste(sort(unique(coloc_class)), collapse = ";"),
  n_strict_groups = uniqueN(signal_group_id),
  n_strict_gwas_groups = sum(n_gwas > 0, na.rm = TRUE),
  strict_pcqtl_specific = any(coloc_class == "pcQTL_specific"),
  strict_shared = any(coloc_class == "shared"),
  strict_eqtl_only = any(coloc_class == "eQTL_only"),
  strict_phenocodes = paste(sort(unique(unlist(strsplit(phenocodes[!is.na(phenocodes)], ",", fixed = TRUE)))), collapse = ","),
  strict_gwas_lead_variants = paste(sort(unique(unlist(strsplit(gwas_lead_variants[!is.na(gwas_lead_variants)], ",", fixed = TRUE)))), collapse = ",")
), by = .(celltype, cluster_id, qtl_type, phenotype_id, credible_set_id)]

strict_cs_class[, strict_coloc_class := fcase(
  strict_pcqtl_specific, "pcQTL_specific",
  strict_shared, "shared",
  strict_eqtl_only, "eQTL_only",
  default = "non_coloc"
)]

effect_summary <- merge(
  cs_summary,
  strict_cs_class,
  by = c("celltype", "cluster_id", "qtl_type", "phenotype_id", "credible_set_id"),
  all.x = TRUE
)
effect_summary[is.na(strict_coloc_class), strict_coloc_class := "non_coloc"]
effect_summary[is.na(n_strict_groups), `:=`(
  strict_group_ids = "",
  strict_coloc_classes = "",
  n_strict_groups = 0L,
  n_strict_gwas_groups = 0L,
  strict_pcqtl_specific = FALSE,
  strict_shared = FALSE,
  strict_eqtl_only = FALSE,
  strict_phenocodes = "",
  strict_gwas_lead_variants = ""
)]
fwrite(effect_summary, file.path(out_dir, "strict_qtl_cs_nominal_effect_summary.tsv"), sep = "\t")

gene_effects_strict <- merge(
  gene_effects,
  strict_cs_class,
  by = c("celltype", "cluster_id", "qtl_type", "phenotype_id", "credible_set_id"),
  all.x = TRUE
)
gene_effects_strict[is.na(strict_coloc_class), strict_coloc_class := "non_coloc"]
fwrite(gene_effects_strict, file.path(out_dir, "strict_pip_weighted_nominal_gene_effects.tsv"), sep = "\t")

message("Building strict pcQTL-specific mechanism table")
pc_group_members <- qtl_group_membership[coloc_class == "pcQTL_specific" & qtl_type == "pcQTL" & !is.na(credible_set_id)]
pc_group_members <- unique(pc_group_members, by = c("signal_group_id", "celltype", "cluster_id", "phenotype_id", "credible_set_id"))
pc_group_members <- merge(
  pc_group_members,
  effect_summary[, .(
    celltype, cluster_id, qtl_type, phenotype_id, credible_set_id,
    n_genes, n_genes_with_effect, max_abs_effect, mean_abs_effect, sd_abs_effect,
    cv_effect, second_largest_abs_effect, fraction_largest
  )],
  by = c("celltype", "cluster_id", "qtl_type", "phenotype_id", "credible_set_id"),
  all.x = TRUE
)

top_effect_genes <- gene_effects_strict[
  order(-abs_pip_weighted_nominal_beta),
  .(top_effect_genes = paste(
    sprintf("%s(%.4g)", gene_id[seq_len(min(5L, .N))], pip_weighted_nominal_beta[seq_len(min(5L, .N))]),
    collapse = ";"
  )),
  by = .(celltype, cluster_id, qtl_type, phenotype_id, credible_set_id)
]
pc_group_members <- merge(
  pc_group_members,
  top_effect_genes,
  by = c("celltype", "cluster_id", "qtl_type", "phenotype_id", "credible_set_id"),
  all.x = TRUE
)
if (nrow(pc_group_members)) {
  pc_group_members[, top_loading_genes := mapply(read_loadings, celltype, cluster_id, phenotype_id)]
}
setorder(pc_group_members, celltype, cluster_id, signal_group_id, phenotype_id, credible_set_id)
fwrite(pc_group_members, file.path(out_dir, "strict_pcqtl_specific_mechanism_table_nominal_effect.tsv"), sep = "\t")

message("Testing distributed-effect differences")
test_input <- effect_summary[
  strict_coloc_class %in% c("pcQTL_specific", "shared", "eQTL_only", "non_coloc") &
    is.finite(max_abs_effect) & is.finite(cv_effect)
]
test_pairs <- list(
  c("pcQTL_specific", "non_coloc"),
  c("pcQTL_specific", "eQTL_only"),
  c("shared", "non_coloc")
)
tests <- rbindlist(lapply(test_pairs, function(pair) {
  x <- test_input[strict_coloc_class %in% pair]
  if (length(unique(x$strict_coloc_class)) < 2L) return(NULL)
  data.table(
    comparison = paste(pair, collapse = "_vs_"),
    metric = c("max_abs_effect", "cv_effect", "fraction_largest"),
    wilcox_p = c(
      wilcox.test(max_abs_effect ~ strict_coloc_class, data = x)$p.value,
      wilcox.test(cv_effect ~ strict_coloc_class, data = x)$p.value,
      wilcox.test(fraction_largest ~ strict_coloc_class, data = x)$p.value
    ),
    n_group1 = x[strict_coloc_class == pair[1], .N],
    n_group2 = x[strict_coloc_class == pair[2], .N],
    median_group1 = c(
      median(x[strict_coloc_class == pair[1], max_abs_effect], na.rm = TRUE),
      median(x[strict_coloc_class == pair[1], cv_effect], na.rm = TRUE),
      median(x[strict_coloc_class == pair[1], fraction_largest], na.rm = TRUE)
    ),
    median_group2 = c(
      median(x[strict_coloc_class == pair[2], max_abs_effect], na.rm = TRUE),
      median(x[strict_coloc_class == pair[2], cv_effect], na.rm = TRUE),
      median(x[strict_coloc_class == pair[2], fraction_largest], na.rm = TRUE)
    )
  )
}), fill = TRUE)
if (nrow(tests)) tests[, fdr := p.adjust(wilcox_p, method = "BH"), by = metric]
fwrite(tests, file.path(out_dir, "strict_distributed_effect_wilcoxon_tests.tsv"), sep = "\t")

message("Wrote strict gene-effect outputs to ", out_dir)
