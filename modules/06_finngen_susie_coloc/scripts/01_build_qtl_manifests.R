#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))

root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1] |> sub("--file=", "", x = _)), ".."), mustWork = FALSE)
if (is.na(root) || root == ".") root <- getwd()
source(file.path(root, "config", "config.R"))

celltypes <- names(CELLTYPE_EQTL_MAP)

dir_create(file.path(ROOT_DIR, "manifests"))

make_eqtl_member <- function(tar_path, gene) {
  inner <- sub("\\.tar\\.gz$", "", basename(tar_path))
  suffix <- sub("^cis_", "", inner)
  file.path(inner, paste0(gene, "_", suffix, "_count_saigeqtl_cis_window_1000000.singleVar.txt"))
}

list_tar_members <- function(tar_path) {
  if (!file.exists(tar_path)) return(character())
  out <- tryCatch(
    system2("tar", c("-tzf", tar_path), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  sub("^\\./", "", out)
}

cluster_rows <- list()
pc_rows <- list()
eqtl_rows <- list()
ld_rows <- list()

for (ct in celltypes) {
  cell_dir <- file.path(PCQTL_ROOT, "celltypes", ct)
  pcqtl_dir <- file.path(cell_dir, "pcQTL")
  cluster_file <- file.path(cell_dir, "cluster_identification/results/method2_sc_hurdle/merged_clusters/cluster_summary.tsv")
  gene_file <- file.path(cell_dir, "cluster_identification/results/method2_sc_hurdle/merged_clusters/cluster_gene_assignments.tsv")
  pc_map_file <- file.path(pcqtl_dir, "step3_saige/cluster_pc_map.tsv")
  if (!file.exists(cluster_file) || !file.exists(gene_file) || !file.exists(pc_map_file)) {
    warning("Skipping ", ct, ": missing cluster or pcQTL manifest files")
    next
  }

  clusters <- fread(cluster_file)
  genes <- fread(gene_file)
  pc_map <- fread(pc_map_file)
  clusters <- clusters[cluster_id %in% pc_map$cluster_id]
  clusters[, celltype := ct]
  clusters[, genes := as.character(genes)]
  clusters[, region_file := file.path(pcqtl_dir, "step3_saige/regions", paste0(cluster_id, "_region.txt"))]
  clusters[, region_exists := TRUE]
  # Existing SAIGE regions are generated as cluster start/end +/- 500 kb.
  # Use the same deterministic definition here instead of opening thousands of
  # tiny NFS files during manifest construction.
  clusters[, `:=`(
    region_chr = sub("^chr", "", as.character(chromosome)),
    region_start = pmax(0L, as.integer(start_position) - 500000L),
    region_end = as.integer(end_position) + 500000L
  )]
  clusters[, region_source := "pcQTL_step3_saige_regions_cluster_pm500kb"]
  cluster_rows[[ct]] <- clusters

  pc_map[, celltype := ct]
  pc_map[, pcqtl_sumstats := file.path(pcqtl_dir, "step3_saige/step2", cluster_id, PC)]
  pc_map[, sumstats_exists := TRUE]
  pc_rows[[ct]] <- pc_map[, .(
    celltype, cluster_id, phenotype_type = "pcQTL", phenotype_id = PC,
    chromosome = sub("^chr", "", as.character(chromosome)),
    pheno_file, region_file, sumstats_path = pcqtl_sumstats, sumstats_member = NA_character_,
    sumstats_exists
  )]

  tar_path <- file.path(EQTL_ROOT, CELLTYPE_EQTL_MAP[[ct]])
  tar_members <- list_tar_members(tar_path)
  ct_genes <- unique(genes[, .(celltype = ct, cluster_id, chromosome = sub("^chr", "", as.character(chromosome)), gene_name)])
  ct_genes[, `:=`(
    phenotype_type = "eQTL",
    phenotype_id = gene_name,
    region_file = file.path(pcqtl_dir, "step3_saige/regions", paste0(cluster_id, "_region.txt")),
    sumstats_path = tar_path,
    sumstats_member = if (file.exists(tar_path)) vapply(gene_name, function(g) make_eqtl_member(tar_path, g), character(1)) else NA_character_
  )]
  ct_genes[, sumstats_exists := file.exists(sumstats_path) & !is.na(sumstats_member) & sumstats_member %chin% tar_members]
  eqtl_rows[[ct]] <- ct_genes[, .(
    celltype, cluster_id, phenotype_type, phenotype_id, chromosome,
    pheno_file = NA_character_, region_file, sumstats_path, sumstats_member,
    sumstats_exists
  )]

  ld_rows[[ct]] <- unique(clusters[, .(
    celltype, cluster_id,
    chromosome = region_chr,
    region_start, region_end, region_file,
    genotype_prefix = sprintf(ONEK1K_MAF005_GENOTYPE_PREFIX, region_chr),
    ld_prefix = file.path(ROOT_DIR, "results/ld/onek1k", celltype, cluster_id),
    ld_matrix = file.path(ROOT_DIR, "results/ld/onek1k", celltype, paste0(cluster_id, ".ld.gz")),
    ld_variants = file.path(ROOT_DIR, "results/ld/onek1k", celltype, paste0(cluster_id, ".variants.tsv"))
  )])
}

cluster_manifest <- rbindlist(cluster_rows, fill = TRUE)
qtl_pheno_manifest <- rbindlist(c(pc_rows, eqtl_rows), fill = TRUE)
ld_manifest <- rbindlist(ld_rows, fill = TRUE)

fwrite(cluster_manifest, file.path(ROOT_DIR, "manifests/cluster_manifest.tsv"), sep = "\t", quote = FALSE)
fwrite(qtl_pheno_manifest, file.path(ROOT_DIR, "manifests/qtl_phenotype_manifest.tsv"), sep = "\t", quote = FALSE)
fwrite(ld_manifest, file.path(ROOT_DIR, "manifests/qtl_ld_tasks.tsv"), sep = "\t", quote = FALSE)

qtl_tasks <- merge(qtl_pheno_manifest, ld_manifest[, .(celltype, cluster_id, ld_matrix, ld_variants)],
                   by = c("celltype", "cluster_id"), all.x = TRUE)
qtl_tasks[, output_file := file.path(ROOT_DIR, "results/fine_mapping/qtl", celltype, cluster_id, phenotype_type, paste0(phenotype_id, ".credible_sets.tsv"))]
qtl_tasks[, status := fifelse(sumstats_exists, "ready", "missing_input")]
fwrite(qtl_tasks, file.path(ROOT_DIR, "manifests/qtl_finemap_tasks.tsv"), sep = "\t", quote = FALSE)

qc <- qtl_tasks[, .(
  n_tasks = .N,
  n_ready = sum(status == "ready"),
  n_missing_sumstats = sum(!sumstats_exists),
  n_missing_region = sum(!region_file %chin% cluster_manifest$region_file)
), by = .(celltype, phenotype_type)]
fwrite(qc, file.path(ROOT_DIR, "manifests/qtl_manifest_qc.tsv"), sep = "\t", quote = FALSE)

message("Wrote manifests under ", file.path(ROOT_DIR, "manifests"))
