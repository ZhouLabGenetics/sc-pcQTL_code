#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))

script <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1])
root <- normalizePath(file.path(dirname(script), ".."))
source(file.path(root, "config", "config.R"))
source(file.path(root, "scripts", "common.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

phenotype_manifest <- get_arg("--phenotype-manifest", file.path(ROOT_DIR, "manifests/finngen_coloc_gwas_phenotype_standard_all_sources.tsv"))
chunk_index <- as.integer(get_arg("--chunk-index", NA_character_))
n_chunks <- as.integer(get_arg("--n-chunks", NA_character_))
out_dir <- get_arg("--out-dir", file.path(ROOT_DIR, "results/fine_mapping/gwas_all_finemapped"))
max_phenotypes <- as.integer(get_arg("--max-phenotypes", NA_character_))
phenocode_filter <- get_arg("--phenocode", NA_character_)
endpoint_source_filter <- get_arg("--endpoint-source", NA_character_)
chunk_label <- get_arg("--chunk-label", NA_character_)

read_bgz <- function(path) fread(cmd = paste("zcat", shQuote(path)))

finngen_susie_snp_columns <- function() {
  c(
    "trait", "region", "v", "rsid", "chromosome", "position",
    "allele1", "allele2", "maf", "beta", "se", "p",
    "mean", "sd", "prob", "cs", "cs_specific_prob", "low_purity",
    "lead_r2", "mean_99", "sd_99", "prob_99", "cs_99",
    "cs_specific_prob_99", "low_purity_99", "lead_r2_99",
    paste0("alpha", 1:10),
    paste0("mean", 1:10),
    paste0("sd", 1:10),
    paste0("lbf_variable", 1:10)
  )
}

read_susie_snp_header <- function(path) {
  header_line <- tryCatch(
    suppressWarnings(system2("bash", c("-c", paste("gzip -dc", shQuote(path), "| head -n 1")), stdout = TRUE, stderr = FALSE)),
    error = function(e) character()
  )
  if (length(header_line) && nzchar(header_line[1])) {
    cols <- tryCatch(
      suppressWarnings(strsplit(header_line[1], "\t", fixed = TRUE)[[1]]),
      error = function(e) character()
    )
    if (length(cols) == length(finngen_susie_snp_columns()) && identical(cols[1:6], finngen_susie_snp_columns()[1:6])) {
      return(cols)
    }
  }
  finngen_susie_snp_columns()
}

find_tabix <- function() {
  tabix <- Sys.which("tabix")
  if (nzchar(tabix)) return(unname(tabix))
  candidates <- c(
    "/apps/software/HTSlib/1.21-GCC-13.3.0/bin/tabix",
    "/usr/bin/tabix",
    "/usr/local/bin/tabix"
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) return(hit[1])
  NA_character_
}

merge_query_ranges <- function(ranges) {
  ranges <- unique(ranges[is.finite(start) & is.finite(end) & end >= start])
  if (!nrow(ranges)) return(ranges)
  ranges[, chromosome := paste0("chr", sub("^chr", "", as.character(chromosome)))]
  setorder(ranges, chromosome, start, end)
  out <- vector("list", nrow(ranges))
  k <- 0L
  for (chr in unique(ranges$chromosome)) {
    x <- ranges[chromosome == chr]
    cur_start <- x$start[1]
    cur_end <- x$end[1]
    if (nrow(x) > 1L) {
      for (i in 2:nrow(x)) {
        if (x$start[i] <= cur_end + 1L) {
          cur_end <- max(cur_end, x$end[i])
        } else {
          k <- k + 1L
          out[[k]] <- data.table(chromosome = chr, start = cur_start, end = cur_end)
          cur_start <- x$start[i]
          cur_end <- x$end[i]
        }
      }
    }
    k <- k + 1L
    out[[k]] <- data.table(chromosome = chr, start = cur_start, end = cur_end)
  }
  rbindlist(out[seq_len(k)])
}

read_bgz_regions <- function(path, ranges, batch_size = 200L) {
  ranges <- merge_query_ranges(ranges)
  if (!nrow(ranges)) return(data.table())

  col_names <- read_susie_snp_header(path)

  tabix <- find_tabix()
  if (is.na(tabix) || !file.exists(paste0(path, ".tbi"))) {
    warning("tabix or .tbi unavailable; falling back to full read for ", path)
    dt <- read_bgz(path)
    dt[, `:=`(chromosome = as.character(chromosome), position = as.integer(position))]
    return(dt[
      ranges,
      on = .(chromosome, position >= start, position <= end),
      nomatch = 0L
    ])
  }

  query <- paste0(ranges$chromosome, ":", ranges$start, "-", ranges$end)
  batches <- split(query, ceiling(seq_along(query) / batch_size))
  dts <- lapply(batches, function(q) {
    cmd <- paste(shQuote(tabix), shQuote(path), paste(shQuote(q), collapse = " "))
    tryCatch(
      fread(cmd = cmd, header = FALSE, col.names = col_names, showProgress = FALSE),
      error = function(e) data.table()
    )
  })
  dt <- rbindlist(dts, fill = TRUE)
  if (nrow(dt)) dt <- unique(dt)
  dt
}

parse_region <- function(x) {
  y <- sub("^chr", "", as.character(x))
  m <- tstrsplit(y, "[:-]")
  data.table(region_chr = m[[1]], region_start = as.integer(m[[2]]), region_end = as.integer(m[[3]]))
}

lift_cluster_windows_to_hg38 <- function(clusters) {
  clusters <- copy(clusters)
  clusters[, `:=`(
    cluster_row = .I,
    region_chr = sub("^chr", "", as.character(region_chr)),
    region_start_hg19 = as.integer(region_start),
    region_end_hg19 = as.integer(region_end)
  )]
  lifted_rows <- list()
  for (chr in sort(unique(clusters$region_chr))) {
    chr_clusters <- clusters[region_chr == chr]
    map_path <- sprintf(ONEK1K_QTL_LIFTOVER_POSMAP_PATTERN, chr)
    if (!file.exists(map_path)) {
      stop("Missing required OneK1K hg19-to-hg38 position map: ", map_path)
    }
    posmap <- fread(map_path, header = FALSE, col.names = c("old_id", "pos_hg38"))
    posmap[, `:=`(
      region_chr = sub("^chr", "", sub(":.*$", "", old_id)),
      pos_hg19 = as.integer(sub("^.*:", "", old_id)),
      pos_hg38 = as.integer(pos_hg38)
    )]
    posmap <- unique(posmap[
      region_chr == chr & is.finite(pos_hg19) & is.finite(pos_hg38),
      .(region_chr, pos_hg19, pos_hg38)
    ], by = c("region_chr", "pos_hg19"))
    if (!nrow(posmap)) {
      stop("Required OneK1K liftover map contains no valid chromosome ", chr, " records: ", map_path)
    }
    setkey(posmap, region_chr, pos_hg19)
    windows <- chr_clusters[, .(cluster_row, region_chr, start = region_start_hg19, end = region_end_hg19)]
    setkey(windows, region_chr, start, end)
    hits <- foverlaps(
      posmap[, .(region_chr, start = pos_hg19, end = pos_hg19, pos_hg38)],
      windows,
      by.x = c("region_chr", "start", "end"),
      by.y = c("region_chr", "start", "end"),
      nomatch = 0L
    )
    lifted <- hits[, .(
      n_lifted_variants_in_window = .N,
      region_start_hg38 = min(pos_hg38, na.rm = TRUE),
      region_end_hg38 = max(pos_hg38, na.rm = TRUE)
    ), by = cluster_row]
    missing_rows <- chr_clusters[!cluster_row %in% lifted$cluster_row, cluster_row]
    if (length(missing_rows)) {
      stop(
        "No lifted OneK1K variants were found in ", length(missing_rows),
        " chromosome ", chr, " cluster window(s); refusing an hg19-as-hg38 fallback."
      )
    }
    lifted[, gwas_region_source := "onek1k_variant_liftover_hg38_range"]
    lifted_rows[[chr]] <- lifted
  }
  lifted <- rbindlist(lifted_rows, fill = TRUE)
  clusters <- merge(clusters, lifted, by = "cluster_row", all.x = TRUE, sort = FALSE)
  if (clusters[!is.finite(region_start_hg38) | !is.finite(region_end_hg38), .N]) {
    stop("Incomplete hg19-to-hg38 cluster-window conversion; no coordinate fallback is permitted.")
  }
  clusters[, `:=`(
    region_start = pmax(0L, as.integer(region_start_hg38)),
    region_end = as.integer(region_end_hg38),
    region_source = paste0(region_source, "__", gwas_region_source)
  )]
  clusters[, c("cluster_row", "region_start_hg38", "region_end_hg38") := NULL]
  clusters[]
}

make_official_susie <- function(snp_dt, cred_dt) {
  snp_dt[, variant_key := standard_variant_id(sub("^chr", "", chromosome), position, allele2, allele1)]
  snp_dt[, `:=`(
    chr = sub("^chr", "", chromosome),
    pos = as.integer(position),
    effect_allele = as.character(allele2),
    other_allele = as.character(allele1),
    beta = as.numeric(beta),
    se = as.numeric(se),
    pvalue = as.numeric(p),
    af = as.numeric(maf),
    n = NA_integer_,
    pip = as.numeric(prob)
  )]
  setorder(snp_dt, chr, pos, effect_allele, other_allele)
  stats <- unique(snp_dt[, .(
    variant_key, chr, pos, marker_id = rsid, effect_allele, other_allele,
    beta, se, pvalue, af, n, pip
  )], by = "variant_key")

  signal_rows <- cred_dt[low_purity == FALSE]
  if (!nrow(signal_rows)) signal_rows <- cred_dt
  lbf <- matrix(NA_real_, nrow = nrow(signal_rows), ncol = nrow(stats))
  alpha <- matrix(0, nrow = nrow(signal_rows), ncol = nrow(stats))
  cs <- vector("list", nrow(signal_rows))

  for (i in seq_len(nrow(signal_rows))) {
    rg <- signal_rows$region[i]
    k <- as.integer(signal_rows$cs[i])
    lbf_col <- paste0("lbf_variable", k)
    alpha_col <- paste0("alpha", k)
    idx_region <- which(snp_dt$region == rg)
    idx_stats <- match(snp_dt$variant_key[idx_region], stats$variant_key)
    lbf[i, ] <- -Inf
    if (lbf_col %in% names(snp_dt)) lbf[i, idx_stats] <- as.numeric(snp_dt[[lbf_col]][idx_region])
    if (alpha_col %in% names(snp_dt)) alpha[i, idx_stats] <- as.numeric(snp_dt[[alpha_col]][idx_region])
    cs_idx <- idx_stats[which(snp_dt$region[idx_region] == rg & as.integer(snp_dt$cs[idx_region]) == k)]
    if (!length(cs_idx)) {
      a <- alpha[i, idx_stats]
      cs_idx <- idx_stats[order(a, decreasing = TRUE)][seq_len(min(length(idx_stats), max(1L, signal_rows$cs_size[i])))]
    }
    cs[[i]] <- sort(unique(cs_idx[!is.na(cs_idx)]))
  }

  keep <- lengths(cs) > 0L & rowSums(is.finite(lbf)) > 0L
  lbf <- lbf[keep, , drop = FALSE]
  alpha <- alpha[keep, , drop = FALSE]
  cs <- cs[keep]
  signal_rows <- signal_rows[keep]

  fit <- list(
    alpha = alpha,
    lbf_variable = lbf,
    pip = stats$pip,
    sets = list(
      cs = cs,
      cs_index = seq_along(cs),
      purity = as.matrix(signal_rows[, .(min_abs_corr = cs_min_r2, mean_abs_corr = cs_avg_r2)]),
      coverage = rep(SUSIE_COVERAGE, length(cs))
    )
  )
  class(fit) <- "susie"
  list(fit = fit, stats = stats, signals = signal_rows)
}

if (!file.exists(phenotype_manifest)) stop("Missing phenotype manifest: ", phenotype_manifest)
endpoints <- fread(phenotype_manifest)
required <- c(
  "phenocode", "phenotype", "endpoint_source", "trait_type", "category",
  "official_susie_in_bucket", "local_susie_cred", "local_susie_snp", "local_susie_tbi"
)
miss <- setdiff(required, names(endpoints))
if (length(miss)) stop("Missing phenotype columns: ", paste(miss, collapse = ", "))
endpoints <- endpoints[official_susie_in_bucket == TRUE]
if (!is.na(phenocode_filter)) endpoints <- endpoints[phenocode == phenocode_filter]
if (!is.na(endpoint_source_filter)) endpoints <- endpoints[endpoint_source == endpoint_source_filter]
setorder(endpoints, endpoint_source, phenocode)
if (!is.na(max_phenotypes)) endpoints <- head(endpoints, max_phenotypes)

if (!is.na(chunk_index) || !is.na(n_chunks)) {
  if (is.na(chunk_index) || is.na(n_chunks) || chunk_index < 1L || n_chunks < 1L || chunk_index > n_chunks) {
    stop("--chunk-index and --n-chunks must be valid positive integers")
  }
  chunk_size <- ceiling(nrow(endpoints) / n_chunks)
  start_i <- (chunk_index - 1L) * chunk_size + 1L
  end_i <- min(chunk_index * chunk_size, nrow(endpoints))
  endpoints <- if (start_i <= nrow(endpoints)) endpoints[start_i:end_i] else endpoints[0]
}

clusters <- lift_cluster_windows_to_hg38(fread(file.path(ROOT_DIR, "manifests/cluster_manifest.tsv")))
clusters[, `:=`(
  region_chr = sub("^chr", "", as.character(region_chr)),
  region_start = as.integer(region_start),
  region_end = as.integer(region_end)
)]
cluster_windows <- clusters[, .(
  celltype, cluster_id, chromosome = region_chr,
  region_start, region_end, region_source,
  start = region_start, end = region_end
)]
setkey(cluster_windows, chromosome, start, end)

rds_base <- file.path(out_dir, "rds")
chunk_dir <- file.path(out_dir, "official_import_chunks")
dir_create(rds_base)
dir_create(chunk_dir)

status_rows <- list()
cs_rows <- list()

for (i in seq_len(nrow(endpoints))) {
  ep <- endpoints[i]
  base_status <- data.table(
    phenocode = ep$phenocode,
    phenotype = ep$phenotype,
    endpoint_source = ep$endpoint_source,
    trait_type = ep$trait_type,
    category = ep$category
  )

  if (!file.exists(ep$local_susie_cred) || !file.exists(ep$local_susie_snp)) {
    status_rows[[length(status_rows) + 1L]] <- cbind(base_status, data.table(
      status = "missing_official_susie_files",
      n_overlapping_cluster_windows = 0L,
      n_imported_cluster_windows = 0L,
      message = "Missing local FinnGen official SUSIE.cred.bgz or SUSIE.snp.bgz"
    ))
    next
  }

  tryCatch({
    cred <- read_bgz(ep$local_susie_cred)
    if (!nrow(cred)) {
      status_rows[[length(status_rows) + 1L]] <- cbind(base_status, data.table(
        status = "empty_official_susie_cred",
        n_overlapping_cluster_windows = 0L,
        n_imported_cluster_windows = 0L,
        message = "Official FinnGen SuSiE credible set file has no rows"
      ))
      next
    }
    cred <- cbind(cred, parse_region(cred$region))
    cred_regions <- unique(cred[, .(
      region, region_chr,
      start = region_start,
      end = region_end
    )])
    setkey(cred_regions, region_chr, start, end)
    overlaps <- foverlaps(
      cluster_windows,
      cred_regions,
      by.x = c("chromosome", "start", "end"),
      by.y = c("region_chr", "start", "end"),
      nomatch = 0L
    )
    if (!nrow(overlaps)) {
      status_rows[[length(status_rows) + 1L]] <- cbind(base_status, data.table(
        status = "no_official_gwas_cs_overlap",
        n_overlapping_cluster_windows = 0L,
        n_imported_cluster_windows = 0L,
        message = "No official FinnGen SuSiE region overlaps any OneK1K cluster window"
      ))
      next
    }

    snp_ranges <- unique(overlaps[, .(
      chromosome,
      start = as.integer(region_start),
      end = as.integer(region_end)
    )])
    snp <- read_bgz_regions(ep$local_susie_snp, snp_ranges)
    if (!nrow(snp)) {
      status_rows[[length(status_rows) + 1L]] <- cbind(base_status, data.table(
        status = "empty_official_susie_snp",
        n_overlapping_cluster_windows = nrow(overlaps),
        n_imported_cluster_windows = 0L,
        message = "Official FinnGen SuSiE SNP file has no rows"
      ))
      next
    }
    snp[, `:=`(chr = sub("^chr", "", chromosome), pos = as.integer(position))]
    imported_n <- 0L

    for (j in seq_len(nrow(overlaps))) {
      task <- overlaps[j]
      overlap_cred <- cred[
        region == task$region &
          region_chr == task$chromosome &
          region_end >= task$region_start &
          region_start <= task$region_end
      ]
      snp_sub <- snp[
        region %in% overlap_cred$region &
          chr == task$chromosome &
          pos >= task$region_start &
          pos <= task$region_end
      ]
      if (!nrow(overlap_cred) || !nrow(snp_sub)) next
      obj <- make_official_susie(copy(snp_sub), copy(overlap_cred))
      if (!length(obj$fit$sets$cs)) next

      rds_dir <- file.path(rds_base, ep$phenocode, task$celltype, task$cluster_id)
      dir_create(rds_dir)
      rds <- file.path(rds_dir, "gwas.official_finngen_susie.rds")
      saveRDS(list(
        fit = obj$fit,
        stats = obj$stats,
        metadata = list(
          phenocode = ep$phenocode,
          phenotype = ep$phenotype,
          endpoint_source = ep$endpoint_source,
          trait_type = ep$trait_type,
          celltype = task$celltype,
          cluster_id = task$cluster_id,
          gwas_finemap_source = "FinnGen_R12_official_SuSiE_all_finemapped",
          ld_source = "FinnGen_R12_official_in_sample_LD",
          cred_file = ep$local_susie_cred,
          snp_file = ep$local_susie_snp
        )
      ), rds)

      cs_rows[[length(cs_rows) + 1L]] <- obj$signals[, .(
        phenocode = ep$phenocode,
        phenotype = ep$phenotype,
        endpoint_source = ep$endpoint_source,
        trait_type = ep$trait_type,
        category = ep$category,
        celltype = task$celltype,
        cluster_id = task$cluster_id,
        gwas_region = region,
        gwas_cs = cs,
        cs_log10bf,
        cs_avg_r2,
        cs_min_r2,
        low_purity,
        cs_size,
        gwas_susie_rds = rds,
        gwas_finemap_source = "FinnGen_R12_official_SuSiE_all_finemapped",
        ld_source = "FinnGen_R12_official_in_sample_LD"
      )]
      imported_n <- imported_n + 1L
    }

    status_rows[[length(status_rows) + 1L]] <- cbind(base_status, data.table(
      status = if (imported_n > 0L) "official_susie_imported" else "no_importable_official_cs",
      n_overlapping_cluster_windows = nrow(overlaps),
      n_imported_cluster_windows = imported_n,
      message = sprintf("Imported %d cluster windows from %d overlapping windows", imported_n, nrow(overlaps))
    ))
  }, error = function(e) {
    status_rows[[length(status_rows) + 1L]] <<- cbind(base_status, data.table(
      status = "failed",
      n_overlapping_cluster_windows = NA_integer_,
      n_imported_cluster_windows = NA_integer_,
      message = conditionMessage(e)
    ))
  })
  message(sprintf("[%d/%d] %s (%s)", i, nrow(endpoints), ep$phenocode, ep$endpoint_source))
}

empty_status_dt <- function() {
  data.table(
    phenocode = character(),
    phenotype = character(),
    endpoint_source = character(),
    trait_type = character(),
    category = character(),
    status = character(),
    n_overlapping_cluster_windows = integer(),
    n_imported_cluster_windows = integer(),
    message = character()
  )
}

empty_cs_dt <- function() {
  data.table(
    phenocode = character(),
    phenotype = character(),
    endpoint_source = character(),
    trait_type = character(),
    category = character(),
    celltype = character(),
    cluster_id = character(),
    gwas_region = character(),
    gwas_cs = integer(),
    cs_log10bf = numeric(),
    cs_avg_r2 = numeric(),
    cs_min_r2 = numeric(),
    low_purity = logical(),
    cs_size = integer(),
    gwas_susie_rds = character(),
    gwas_finemap_source = character(),
    ld_source = character()
  )
}

status_dt <- if (length(status_rows)) rbindlist(status_rows, fill = TRUE) else empty_status_dt()
cs_dt <- if (length(cs_rows)) rbindlist(cs_rows, fill = TRUE) else empty_cs_dt()
label_part <- if (!is.na(chunk_label) && nzchar(chunk_label)) paste0("_", chunk_label) else ""
suffix <- if (!is.na(chunk_index)) sprintf("%s_chunk_%03d", label_part, chunk_index) else label_part
fwrite(status_dt, file.path(chunk_dir, paste0("finngen_official_cs_import_status", suffix, ".tsv")), sep = "\t", quote = FALSE)
fwrite(cs_dt, file.path(chunk_dir, paste0("finngen_official_cs_by_cluster", suffix, ".tsv")), sep = "\t", quote = FALSE)
if (nrow(status_dt)) {
  print(status_dt[, .N, by = status][order(-N)])
} else {
  message("No phenotypes assigned to this chunk.")
}
