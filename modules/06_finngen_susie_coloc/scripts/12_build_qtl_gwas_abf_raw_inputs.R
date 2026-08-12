#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))

script <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1])
root <- normalizePath(file.path(dirname(script), ".."))
source(file.path(root, "config", "config.R"))
source(file.path(root, "scripts", "common.R"))

out_dir <- file.path(ROOT_DIR, "results/coloc/qtl_gwas_abf_raw")
cache_dir <- file.path(out_dir, "cache")
dir_create(out_dir)
dir_create(cache_dir)
primary_ct <- names(CELLTYPE_EQTL_MAP)

qtl_status_all <- fread(file.path(ROOT_DIR, "results/fine_mapping/qtl_finemap_status.tsv"))
if (!"celltype" %in% names(qtl_status_all)) stop("QTL fine-mapping status lacks celltype")
qtl_status_all <- qtl_status_all[celltype %chin% primary_ct]
if (!nrow(qtl_status_all)) stop("No primary-analysis cell types remain in QTL fine-mapping status")
qtl_status_all[, qtl_id := paste(celltype, cluster_id, phenotype_type, phenotype_id, sep = "__")]

# Keep these two sets separate:
# - qtl_stats_status: all successful QTL summary statistics, used for ABF windows and plots.
# - qtl_lead_status: credible/significant QTL signals used to define ABF candidate regions.
qtl_stats_status <- qtl_status_all[status == "ok"]
qtl_lead_status <- qtl_stats_status[n_cs > 0]

qtl_stats_cache <- file.path(cache_dir, "qtl_abf_stats.rds")
qtl_leads_file <- file.path(ROOT_DIR, "manifests/qtl_gwas_abf_raw_qtl_leads.tsv")
qtl_liftover_qc_file <- file.path(out_dir, "qtl_gwas_abf_raw_liftover_qc.tsv")
if (file.exists(qtl_stats_cache) && file.exists(qtl_leads_file)) {
  qtl_stats <- as.data.table(readRDS(qtl_stats_cache))
  qtl_leads <- fread(qtl_leads_file)
  if (!"celltype" %in% names(qtl_stats) || !"celltype" %in% names(qtl_leads)) {
    stop("Existing QTL ABF cache lacks celltype; remove it and rebuild")
  }
  n_stats_before <- nrow(qtl_stats)
  n_leads_before <- nrow(qtl_leads)
  qtl_stats <- qtl_stats[celltype %chin% primary_ct]
  qtl_leads <- qtl_leads[celltype %chin% primary_ct]
  if (nrow(qtl_stats) < n_stats_before || nrow(qtl_leads) < n_leads_before) {
    message(sprintf(
      "Purged excluded cell types from existing QTL ABF cache: %d statistic rows and %d leads",
      n_stats_before - nrow(qtl_stats), n_leads_before - nrow(qtl_leads)
    ))
    saveRDS(qtl_stats, qtl_stats_cache)
    fwrite(qtl_leads, qtl_leads_file, sep = "\t", quote = FALSE)
    if (file.exists(qtl_liftover_qc_file)) {
      liftover_qc <- fread(qtl_liftover_qc_file)
      if ("celltype" %in% names(liftover_qc)) {
        liftover_qc <- liftover_qc[celltype %chin% primary_ct]
        fwrite(liftover_qc, qtl_liftover_qc_file, sep = "\t", quote = FALSE)
      }
    }
  }
  message("Reusing existing QTL ABF cache")
} else {
  load_liftover_posmap <- function(chrs) {
    rows <- lapply(sort(unique(as.character(chrs))), function(chr) {
      path <- sprintf(ONEK1K_QTL_LIFTOVER_POSMAP_PATTERN, chr)
      if (!file.exists(path)) stop("Missing OneK1K QTL liftOver posmap: ", path)
      dt <- fread(path, header = FALSE, col.names = c("old_id", "pos_hg38"))
      dt[, `:=`(
        chr = sub("^chr", "", sub(":.*$", "", old_id)),
        pos_hg19 = as.integer(sub("^.*:", "", old_id)),
        pos_hg38 = as.integer(pos_hg38)
      )]
      dt <- dt[is.finite(pos_hg19) & is.finite(pos_hg38)]
      unique(dt[, .(chr, pos_hg19, pos_hg38)], by = c("chr", "pos_hg19"))
    })
    out <- rbindlist(rows, fill = TRUE)
    setkey(out, chr, pos_hg19)
    out
  }

  lift_qtl_stats_to_hg38 <- function(stats, lift_map) {
    stats <- copy(stats)
    stats[, `:=`(
      chr = sub("^chr", "", as.character(chr)),
      pos_hg19 = as.integer(pos)
    )]
    lifted <- merge(stats, lift_map, by = c("chr", "pos_hg19"), all.x = FALSE, all.y = FALSE)
    if (!nrow(lifted)) return(lifted)
    lifted[, pos := as.integer(pos_hg38)]
    lifted[, pos_hg38 := NULL]
    lifted[, marker_id_hg19 := marker_id]
    lifted[, marker_id := standard_variant_id(chr, pos, effect_allele, other_allele)]
    lifted <- unique(lifted, by = c("qtl_id", "chr", "pos", "effect_allele", "other_allele"))
    lifted[]
  }

  lift_map <- load_liftover_posmap(qtl_stats_status$chromosome %||% as.character(1:22))
  stats_rows <- vector("list", nrow(qtl_stats_status))
  lead_rows <- vector("list", nrow(qtl_stats_status))
  liftover_qc_rows <- vector("list", nrow(qtl_stats_status))

  for (i in seq_len(nrow(qtl_stats_status))) {
    row <- qtl_stats_status[i]
    fm <- readRDS(row$susie_rds)
    stats <- as.data.table(copy(fm$stats))
    stats <- stats[is.finite(beta) & is.finite(se) & se > 0 & is.finite(pvalue)]
    stats[, `:=`(
      qtl_id = row$qtl_id,
      celltype = row$celltype,
      cluster_id = row$cluster_id,
      qtl_type = row$phenotype_type,
      phenotype_id = row$phenotype_id,
      qtl_maf = pmin(as.numeric(af), 1 - as.numeric(af)),
      qtl_n = as.integer(n)
    )]
    stats <- stats[is.finite(qtl_maf) & qtl_maf > 0 & qtl_maf < 1]
    n_before_liftover <- nrow(stats)
    stats <- lift_qtl_stats_to_hg38(stats, lift_map)
    n_after_liftover <- nrow(stats)
    liftover_qc_rows[[i]] <- data.table(
      qtl_id = row$qtl_id,
      celltype = row$celltype,
      cluster_id = row$cluster_id,
      qtl_type = row$phenotype_type,
      phenotype_id = row$phenotype_id,
      coordinate_build = ONEK1K_QTL_COORDINATE_BUILD,
      n_before_liftover = n_before_liftover,
      n_after_liftover = n_after_liftover,
      n_liftover_dropped = n_before_liftover - n_after_liftover
    )
    if (!nrow(stats)) next
    if (row$n_cs > 0) {
      lead <- stats[which.min(pvalue)]
      lead_rows[[i]] <- data.table(
        qtl_id = row$qtl_id,
        celltype = row$celltype,
        cluster_id = row$cluster_id,
        qtl_type = row$phenotype_type,
        phenotype_id = row$phenotype_id,
        susie_rds = row$susie_rds,
        lead_chr = lead$chr,
        lead_pos = as.integer(lead$pos),
        lead_pos_hg19 = as.integer(lead$pos_hg19),
        qtl_lead_variant = standard_variant_id(lead$chr, lead$pos, lead$effect_allele, lead$other_allele),
        qtl_lead_variant_hg19 = standard_variant_id(lead$chr, lead$pos_hg19, lead$effect_allele, lead$other_allele),
        qtl_lead_p = as.numeric(lead$pvalue),
        qtl_lead_beta = as.numeric(lead$beta),
        qtl_lead_se = as.numeric(lead$se),
        qtl_n = as.integer(max(stats$qtl_n, na.rm = TRUE)),
        n_qtl_variants = nrow(stats),
        coordinate_build = ONEK1K_QTL_COORDINATE_BUILD
      )
    }
    stats_rows[[i]] <- stats[, .(
      qtl_id, celltype, cluster_id, qtl_type, phenotype_id,
      chr, pos, pos_hg19, marker_id, marker_id_hg19, effect_allele, other_allele,
      beta, se, pvalue, qtl_maf, qtl_n
    )]
    if (i %% 250 == 0) message(sprintf("[%d/%d] cached QTL ABF stats", i, nrow(qtl_stats_status)))
  }

  qtl_stats <- rbindlist(stats_rows, fill = TRUE)
  qtl_leads <- rbindlist(lead_rows, fill = TRUE)
  liftover_qc <- rbindlist(liftover_qc_rows, fill = TRUE)
  setkey(qtl_stats, qtl_id, chr, pos)

  saveRDS(qtl_stats, qtl_stats_cache)
  fwrite(qtl_leads, qtl_leads_file, sep = "\t", quote = FALSE)
  fwrite(liftover_qc, qtl_liftover_qc_file, sep = "\t", quote = FALSE)
}

endpoints <- fread(file.path(ROOT_DIR, "manifests/finngen_main_endpoints.tsv"))
endpoints[, alt_local_gz := file.path(FINNGEN_ROOT, "storage.googleapis.com/finngen-public-data-r12/summary_stats/release", basename(local_gz))]
endpoints[, resolved_local_gz := fifelse(file.exists(local_gz), local_gz,
                                  fifelse(file.exists(alt_local_gz), alt_local_gz, local_gz))]
endpoints[, raw_gwas_available := file.exists(resolved_local_gz)]
endpoints[, endpoint_task_index := .I]
fwrite(endpoints, file.path(ROOT_DIR, "manifests/qtl_gwas_abf_raw_endpoint_tasks.tsv"), sep = "\t", quote = FALSE)

qc <- rbindlist(list(
  data.table(metric = "qtl_stats_successful_phenotypes", value = "n", N = nrow(qtl_stats_status)),
  data.table(metric = "qtl_with_cs", value = "n", N = nrow(qtl_lead_status)),
  data.table(metric = "qtl_leads_cached", value = "n", N = nrow(qtl_leads)),
  data.table(metric = "qtl_stats_rows", value = "n", N = nrow(qtl_stats)),
  endpoints[, .N, by = .(value = as.character(raw_gwas_available))][, metric := "raw_gwas_available"][, .(metric, value, N)]
), fill = TRUE)
fwrite(qc, file.path(out_dir, "qtl_gwas_abf_raw_input_qc.tsv"), sep = "\t", quote = FALSE)

message("Wrote QTL ABF leads: ", nrow(qtl_leads))
message("Wrote QTL ABF stats rows: ", nrow(qtl_stats))
message("Wrote endpoint tasks: ", nrow(endpoints))
