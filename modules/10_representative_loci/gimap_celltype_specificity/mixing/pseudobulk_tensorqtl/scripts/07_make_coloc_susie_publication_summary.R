#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

mix_root <- Sys.getenv("SC_PCQTL_GIMAP_MIXING_ROOT", unset = "")
if (!nzchar(mix_root) || !dir.exists(mix_root)) {
  stop("Set SC_PCQTL_GIMAP_MIXING_ROOT to the completed GIMAP mixing directory.")
}
coloc_dir <- file.path(mix_root, "coloc_susie")
coloc_file <- file.path(coloc_dir, "mix_coloc_susie_summary.tsv")
if (!file.exists(coloc_file)) stop("Missing colocalization input table: ", coloc_file)

result <- fread(coloc_file)
result[, pph4 := fifelse(is.na(pph4_susie), 0, pph4_susie)]
best <- result[, .SD[which.max(pph4)], by = .(scenario, qtl_type)]

pph4 <- dcast(best, scenario ~ qtl_type, value.var = "pph4")
readout <- dcast(best, scenario ~ qtl_type, value.var = "readout")
setnames(readout, c("pcQTL", "eQTL"), c("pcQTL_best_readout", "eQTL_best_readout"), skip_absent = TRUE)
abf <- dcast(best, scenario ~ qtl_type, value.var = "pph4_abf")
setnames(abf, c("pcQTL", "eQTL"), c("pcQTL_abf", "eQTL_abf"), skip_absent = TRUE)
summary <- Reduce(function(x, y) merge(x, y, by = "scenario"), list(pph4, readout, abf))
setnames(summary, c("pcQTL", "eQTL"), c("pcQTL_pph4_susie", "eQTL_pph4_susie"), skip_absent = TRUE)

phenotype <- fread(file.path(mix_root, "data/gimap_pseudobulk_tensorqtl_phenotype_level_results.tsv"))
phenotype <- phenotype[phenotype_set %in% c("pc_only", "single_gene_only")]
phenotype[, qtl_type := fifelse(phenotype_set == "pc_only", "pcQTL", "eQTL")]
detection <- phenotype[is.finite(pval_nominal_cis) & pval_nominal_cis > 0,
  .(neglog10p = -log10(min(pval_nominal_cis))), by = .(scenario, qtl_type)]
detection <- dcast(detection, scenario ~ qtl_type, value.var = "neglog10p")
setnames(
  detection, c("pcQTL", "eQTL"),
  c("pcQTL_detect_neglog10p", "eQTL_detect_neglog10p"), skip_absent = TRUE
)
summary <- merge(summary, detection, by = "scenario", all.x = TRUE)

labels <- c(
  cd8_nc_weight_000_of_100 = "0%",
  cd8_nc_weight_005_of_100 = "5%",
  cd8_nc_weight_010_of_100 = "10%",
  cd8_nc_weight_025_of_100 = "25%",
  cd8_nc_weight_050_of_100 = "50%",
  cd8_nc_weight_075_of_100 = "75%",
  cd8_nc_only = "100%",
  observed_cell_count_weighted = "Observed",
  equal_celltype_weighted = "Equal mix"
)
summary <- summary[scenario %in% names(labels)]
summary[, `:=`(comp_label = labels[scenario], source = "OneK1K pseudobulk")]

output <- copy(summary)
order_labels <- c("0%", "5%", "10%", "25%", "50%", "75%", "100%", "Observed", "Equal mix")
output[, comp_order := match(comp_label, order_labels)]
setorder(output, comp_order)
setcolorder(output, c(
  "comp_order", "comp_label", "source", "scenario",
  "pcQTL_pph4_susie", "eQTL_pph4_susie",
  "pcQTL_best_readout", "eQTL_best_readout", "pcQTL_abf", "eQTL_abf"
))
out_file <- file.path(coloc_dir, "gimap_mixing_coloc_susie_publication_summary.tsv")
fwrite(output, out_file, sep = "\t")
message("Wrote ", out_file)
