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

post_dir <- file.path(root_dir, "05_post_coloc_susie_official_finngen_all_finemapped")
strict_effect_dir <- file.path(root_dir, "09_mechanistic_celltype_analysis/results/gene_effects")
cluster_annotation_root <- Sys.getenv(
  "SC_PCQTL_CLUSTER_ENRICHMENT_ROOT",
  unset = file.path(dirname(root_dir), "03_analysis_celltypes/04_cluster_annotation_enrichment_add_cov")
)
ann_src <- file.path(cluster_annotation_root, "annotations/processed")
if (!dir.exists(ann_src)) stop("Missing processed cluster annotations: ", ann_src)
out_dir <- file.path(root_dir, "09_mechanistic_celltype_analysis/results/regulatory_annotation")
dir_create(out_dir)

promoter_window <- 1000L
tad_boundary_window <- 50000L
threshold_main <- 0.75

std_chr_local <- function(x) sub("^chr", "", as.character(x))

overlap_points <- function(points, intervals, id_col = "row_id") {
  if (!nrow(points) || !nrow(intervals)) return(integer())
  p <- copy(points)
  p[, `:=`(chr = std_chr_local(chr), start = as.integer(pos), end = as.integer(pos))]
  p <- p[chr %in% as.character(1:22)]
  iv <- copy(intervals)
  iv[, `:=`(chr = std_chr_local(chr), start = as.integer(start), end = as.integer(end))]
  iv <- iv[chr %in% as.character(1:22) & is.finite(start) & is.finite(end) & start <= end]
  if (!nrow(p) || !nrow(iv)) return(integer())
  setkey(iv, chr, start, end)
  setkey(p, chr, start, end)
  hit <- foverlaps(p, iv, by.x = c("chr", "start", "end"), by.y = c("chr", "start", "end"),
                   type = "within", nomatch = 0L)
  unique(hit[[id_col]])
}

message("Reading QTL credible set variants")
cs_files <- list.files(
  file.path(root_dir, "results/fine_mapping/qtl"),
  pattern = "\\.credible_sets\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(cs_files)) stop("No QTL credible set TSV files found")
cs <- rbindlist(lapply(cs_files, function(f) {
  x <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  req <- c("celltype", "cluster_id", "phenotype_type", "phenotype_id", "credible_set_id",
           "cs_index", "chr", "pos", "effect_allele", "other_allele", "pip")
  if (!all(req %in% names(x))) return(NULL)
  x[, ..req]
}), fill = TRUE)
if (!nrow(cs)) stop("No credible set variants found")
cs <- cs[celltype %chin% names(CELLTYPE_EQTL_MAP)]
if (!nrow(cs)) stop("No primary-analysis credible set variants remain after eligibility filtering")
cs[, `:=`(
  row_id = .I,
  qtl_type = fifelse(phenotype_type == "pcQTL", "pcQTL", "eQTL"),
  chr = std_chr_local(chr),
  pos = as.integer(pos),
  pip = as.numeric(pip),
  signal_index = as.integer(cs_index)
)]
cs <- cs[is.finite(pos) & is.finite(pip) & pip > 0]

message("Adding strict coloc classes")
strict_class <- fread(file.path(strict_effect_dir, "strict_qtl_cs_nominal_effect_summary.tsv"))
strict_class <- strict_class[celltype %chin% names(CELLTYPE_EQTL_MAP)]
strict_class <- unique(strict_class[, .(
  celltype, cluster_id, qtl_type, phenotype_id, credible_set_id,
  strict_coloc_class, strict_group_ids, strict_phenocodes, strict_gwas_lead_variants
)])
cs <- merge(
  cs,
  strict_class,
  by = c("celltype", "cluster_id", "qtl_type", "phenotype_id", "credible_set_id"),
  all.x = TRUE
)
cs[is.na(strict_coloc_class), strict_coloc_class := "non_coloc"]
cs[, variant_id := paste0(chr, ":", pos, "_", effect_allele, "/", other_allele)]

message("Reading cluster gene assignments")
assign_files <- list.files(
  file.path(post_dir, "cache/cluster_gene_assignments"),
  pattern = "\\.cluster_gene_assignments\\.tsv$",
  full.names = TRUE
)
gene_map <- rbindlist(lapply(assign_files, function(f) {
  ct <- sub("\\.cluster_gene_assignments\\.tsv$", "", basename(f))
  x <- fread(f, showProgress = FALSE)
  x[, celltype := ct]
  x
}), fill = TRUE)
gene_map <- gene_map[celltype %chin% names(CELLTYPE_EQTL_MAP)]
gene_map <- unique(gene_map[, .(celltype, cluster_id, gene_id = gene_name)])

message("Reading regulatory annotation resources")
gene_model <- fread(file.path(ann_src, "gene_model.tsv"))
gene_model[, `:=`(
  chr = std_chr_local(chr),
  start = as.integer(start),
  end = as.integer(end),
  gene_id = gene_name,
  tss = fifelse(strand == "-", as.integer(end), as.integer(start))
)]
gene_model <- unique(gene_model[chr %in% as.character(1:22), .(gene_id, chr, start, end, strand, tss)])
gene_model[, `:=`(
  promoter_start = pmax(1L, tss - promoter_window),
  promoter_end = tss + promoter_window
)]

abc <- fread(file.path(ann_src, "abc_links.tsv"))
if (nrow(abc)) {
  abc[, `:=`(chr = std_chr_local(chr), start = as.integer(start), end = as.integer(end), gene_id = gene_name)]
  abc <- abc[chr %in% as.character(1:22)]
}
ctcf <- fread(file.path(ann_src, "ctcf_peaks.tsv"))
tad <- fread(file.path(ann_src, "tad_boundaries.tsv"))
if (nrow(tad)) {
  tad[, `:=`(
    chr = std_chr_local(chr),
    start = pmax(1L, as.integer(start) - tad_boundary_window),
    end = as.integer(end) + tad_boundary_window
  )]
}

message("Annotating gene-specific promoter/body effects")
cs_gene <- merge(
  cs[, .(row_id, celltype, cluster_id, qtl_type, phenotype_id, credible_set_id, variant_id, chr, pos, pip, strict_coloc_class)],
  gene_map,
  by = c("celltype", "cluster_id"),
  allow.cartesian = TRUE
)
cs_gene <- merge(cs_gene, gene_model, by = "gene_id", all.x = TRUE, suffixes = c("", "_gene"))
cs_gene[, `:=`(
  same_chr_gene = chr == chr_gene,
  in_gene_body = chr == chr_gene & pos >= start & pos <= end,
  in_promoter = chr == chr_gene & pos >= promoter_start & pos <= promoter_end,
  distance_to_tss = fifelse(chr == chr_gene, abs(pos - tss), NA_integer_)
)]

if (nrow(abc)) {
  points <- unique(cs_gene[, .(gene_variant_row = .I, row_id, gene_id, chr, pos)])
  p <- copy(points)
  p[, `:=`(start = pos, end = pos)]
  abc_iv <- unique(abc[, .(gene_id, chr, start, end)])
  setkey(abc_iv, chr, start, end)
  setkey(p, chr, start, end)
  abc_hit <- foverlaps(p, abc_iv, by.x = c("chr", "start", "end"), by.y = c("chr", "start", "end"),
                       type = "within", nomatch = 0L)
  abc_hit <- unique(abc_hit[gene_id == i.gene_id, .(row_id, gene_id)])
  cs_gene[, in_abc_link := FALSE]
  if (nrow(abc_hit)) {
    cs_gene[paste(row_id, gene_id) %in% paste(abc_hit$row_id, abc_hit$gene_id), in_abc_link := TRUE]
  }
} else {
  cs_gene[, in_abc_link := FALSE]
}

message("Annotating generic variant-level regulatory features")
cs[, in_ctcf_peak := row_id %in% overlap_points(cs[, .(row_id, chr, pos)], ctcf)]
cs[, near_tad_boundary := row_id %in% overlap_points(cs[, .(row_id, chr, pos)], tad)]
if (nrow(abc)) {
  cs[, in_any_abc_enhancer := row_id %in% overlap_points(cs[, .(row_id, chr, pos)], abc)]
} else {
  cs[, in_any_abc_enhancer := FALSE]
}

row_gene_any <- cs_gene[, .(
  in_cluster_gene_body = any(in_gene_body, na.rm = TRUE),
  in_cluster_gene_promoter = any(in_promoter, na.rm = TRUE),
  in_cluster_gene_abc_link = any(in_abc_link, na.rm = TRUE),
  min_distance_to_cluster_gene_tss = as.numeric(suppressWarnings(min(as.numeric(distance_to_tss), na.rm = TRUE)))
), by = row_id]
row_gene_any[!is.finite(min_distance_to_cluster_gene_tss), min_distance_to_cluster_gene_tss := NA_real_]
cs <- merge(cs, row_gene_any, by = "row_id", all.x = TRUE)
for (col in c("in_cluster_gene_body", "in_cluster_gene_promoter", "in_cluster_gene_abc_link")) {
  cs[is.na(get(col)), (col) := FALSE]
}

variant_out <- cs[, .(
  celltype, cluster_id, qtl_type, phenotype_id, credible_set_id, signal_index,
  strict_coloc_class, strict_group_ids, strict_phenocodes, strict_gwas_lead_variants,
  variant_id, chr, pos, effect_allele, other_allele, pip,
  in_cluster_gene_promoter, in_cluster_gene_body, in_cluster_gene_abc_link,
  in_any_abc_enhancer, in_ctcf_peak, near_tad_boundary,
  min_distance_to_cluster_gene_tss
)]
fwrite(variant_out, file.path(out_dir, "strict_qtl_cs_variant_regulatory_annotations.tsv"), sep = "\t")

message("Summarizing PIP-weighted annotation probabilities")
pip_prob <- function(flag, pip) {
  if (!length(pip) || sum(pip, na.rm = TRUE) <= 0) return(NA_real_)
  sum(as.numeric(flag) * pip, na.rm = TRUE) / sum(pip, na.rm = TRUE)
}
cs_summary <- cs[, .(
  n_cs_variants = .N,
  pip_sum = sum(pip, na.rm = TRUE),
  pip_weighted_promoter = pip_prob(in_cluster_gene_promoter, pip),
  pip_weighted_gene_body = pip_prob(in_cluster_gene_body, pip),
  pip_weighted_abc_link = pip_prob(in_cluster_gene_abc_link, pip),
  pip_weighted_any_abc_enhancer = pip_prob(in_any_abc_enhancer, pip),
  pip_weighted_ctcf = pip_prob(in_ctcf_peak, pip),
  pip_weighted_tad_boundary_50kb = pip_prob(near_tad_boundary, pip),
  min_distance_to_cluster_gene_tss = suppressWarnings(min(min_distance_to_cluster_gene_tss, na.rm = TRUE))
), by = .(
  celltype, cluster_id, qtl_type, phenotype_id, credible_set_id, signal_index,
  strict_coloc_class, strict_group_ids, strict_phenocodes, strict_gwas_lead_variants
)]
cs_summary[!is.finite(min_distance_to_cluster_gene_tss), min_distance_to_cluster_gene_tss := NA_real_]
fwrite(cs_summary, file.path(out_dir, "strict_qtl_cs_pip_weighted_regulatory_annotation_summary.tsv"), sep = "\t")

gene_ann <- cs_gene[, .(
  n_cs_variants = .N,
  pip_sum = sum(pip, na.rm = TRUE),
  pip_weighted_gene_promoter = pip_prob(in_promoter, pip),
  pip_weighted_gene_body = pip_prob(in_gene_body, pip),
  pip_weighted_gene_abc_link = pip_prob(in_abc_link, pip),
  min_distance_to_tss = as.numeric(suppressWarnings(min(as.numeric(distance_to_tss), na.rm = TRUE)))
), by = .(celltype, cluster_id, qtl_type, phenotype_id, credible_set_id, strict_coloc_class, gene_id)]
gene_ann[!is.finite(min_distance_to_tss), min_distance_to_tss := NA_real_]
fwrite(gene_ann, file.path(out_dir, "strict_qtl_cs_gene_pip_weighted_regulatory_annotation_summary.tsv"), sep = "\t")

message("Updating strict pcQTL-specific mechanism table with annotations")
mech <- fread(file.path(strict_effect_dir, "strict_pcqtl_specific_mechanism_table_nominal_effect.tsv"))
mech_ann <- merge(
  mech,
  cs_summary[, .(
    celltype, cluster_id, qtl_type, phenotype_id, credible_set_id,
    pip_weighted_promoter, pip_weighted_gene_body, pip_weighted_abc_link,
    pip_weighted_any_abc_enhancer, pip_weighted_ctcf, pip_weighted_tad_boundary_50kb,
    min_distance_to_cluster_gene_tss
  )],
  by = c("celltype", "cluster_id", "qtl_type", "phenotype_id", "credible_set_id"),
  all.x = TRUE
)
fwrite(mech_ann, file.path(out_dir, "strict_pcqtl_specific_mechanism_table_with_regulatory_annotation.tsv"), sep = "\t")

ann_long <- melt(
  cs_summary[strict_coloc_class %in% c("pcQTL_specific", "shared", "eQTL_only", "non_coloc")],
  id.vars = c("celltype", "cluster_id", "qtl_type", "phenotype_id", "credible_set_id", "strict_coloc_class"),
  measure.vars = c("pip_weighted_promoter", "pip_weighted_gene_body", "pip_weighted_abc_link",
                   "pip_weighted_ctcf", "pip_weighted_tad_boundary_50kb"),
  variable.name = "annotation",
  value.name = "pip_weighted_probability"
)
ann_long[, annotation := factor(annotation, levels = c(
  "pip_weighted_promoter", "pip_weighted_gene_body", "pip_weighted_abc_link",
  "pip_weighted_ctcf", "pip_weighted_tad_boundary_50kb"
), labels = c("Promoter", "Gene body", "ABC link", "CTCF", "TAD boundary +/-50kb"))]
ann_long[, strict_coloc_class := factor(strict_coloc_class, levels = c("pcQTL_specific", "shared", "eQTL_only", "non_coloc"))]

tests <- ann_long[
  strict_coloc_class %in% c("pcQTL_specific", "non_coloc") & is.finite(pip_weighted_probability),
  {
    if (length(unique(strict_coloc_class)) < 2L) {
      .(wilcox_p = NA_real_, n_pcqtl_specific = sum(strict_coloc_class == "pcQTL_specific"), n_non_coloc = sum(strict_coloc_class == "non_coloc"),
        median_pcqtl_specific = NA_real_, median_non_coloc = NA_real_)
    } else {
      .(
        wilcox_p = wilcox.test(pip_weighted_probability ~ strict_coloc_class)$p.value,
        n_pcqtl_specific = sum(strict_coloc_class == "pcQTL_specific"),
        n_non_coloc = sum(strict_coloc_class == "non_coloc"),
        median_pcqtl_specific = median(pip_weighted_probability[strict_coloc_class == "pcQTL_specific"], na.rm = TRUE),
        median_non_coloc = median(pip_weighted_probability[strict_coloc_class == "non_coloc"], na.rm = TRUE)
      )
    }
  },
  by = annotation
]
tests[, fdr := p.adjust(wilcox_p, method = "BH")]
fwrite(tests, file.path(out_dir, "strict_regulatory_annotation_pcqtl_specific_vs_non_coloc_tests.tsv"), sep = "\t")

message("Wrote strict regulatory annotation outputs to ", out_dir)
