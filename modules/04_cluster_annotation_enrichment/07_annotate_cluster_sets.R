#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))
source("config/config.R")

cluster_dir <- file.path(ROOT_DIR, "results", "cluster_sets")
null_dir <- file.path(ROOT_DIR, "results", "null_sets")
ann_dir <- file.path(ROOT_DIR, "annotations", "processed")
out_dir <- file.path(ROOT_DIR, "results", "annotations")
dir_create(out_dir)

add_cov <- safe_fread(file.path(cluster_dir, "add_cov_correlated_clusters.tsv"))
null_add <- safe_fread(file.path(null_dir, "add_cov_sc_hurdle_null_neighbor_clusters.tsv"))
sets <- rbindlist(list(add_cov, null_add), fill = TRUE)
sets <- sets[num_genes %in% MAIN_CLUSTER_SIZES]
sets[, `:=`(
  set_id = paste(method, celltype, cluster_id, sep = "__"),
  log_cluster_length = log10(pmax(cluster_length, 1L))
)]

gene_model <- safe_fread(file.path(ann_dir, "gene_model.tsv"))
gene_model <- unique(gene_model[, .(gene_name, chr, start, end, strand)])
setkey(gene_model, gene_name)
go_bp <- safe_fread(file.path(ann_dir, "go_bp_gene_terms.tsv"))
paralog <- safe_fread(file.path(ann_dir, "paralog_pairs.tsv"))
abc <- safe_fread(file.path(ann_dir, "abc_links.tsv"))
ctcf <- safe_fread(file.path(ann_dir, "ctcf_peaks.tsv"))
tad <- safe_fread(file.path(ann_dir, "tad_boundaries.tsv"))

if (nrow(paralog)) {
  paralog[, key := apply(.SD, 1, function(z) paste(sort(z), collapse = "\t")), .SDcols = c("gene1", "gene2")]
  paralog_set <- unique(paralog$key)
} else {
  paralog_set <- character()
}

go_index <- if (nrow(go_bp)) split(go_bp$go_id, go_bp$gene_name) else list()

if (nrow(abc) && all(c("chr", "start", "end", "gene_name") %in% names(abc))) {
  if ("enhancer_name" %in% names(abc)) {
    abc[, enhancer_id := enhancer_name]
  } else {
    abc[, enhancer_id := paste(std_chr(chr), start, end, sep = ":")]
  }
  abc <- unique(abc[, .(gene_name, enhancer_id)])
  setkey(abc, gene_name)
}
abc_index <- if (nrow(abc)) split(abc$enhancer_id, abc$gene_name) else list()

make_position_index <- function(tab) {
  if (!nrow(tab) || !all(c("chr", "start", "end") %in% names(tab))) return(list())
  x <- copy(tab)
  x[, chr := std_chr(chr)]
  x <- x[chr %in% as.character(1:22)]
  pos <- unique(rbind(
    x[, .(chr, pos = as.integer(start))],
    x[, .(chr, pos = as.integer(end))]
  ))
  split_pos <- split(pos$pos, pos$chr)
  lapply(split_pos, function(v) sort(unique(v)))
}

ctcf_pos <- make_position_index(ctcf)
tad_pos <- make_position_index(tad)

has_shared_go <- function(gs) {
  terms <- unlist(go_index[intersect(gs, names(go_index))], use.names = FALSE)
  if (!length(terms)) return(FALSE)
  any(tabulate(match(terms, unique(terms))) >= 2L)
}

has_paralog <- function(gs) {
  if (length(gs) < 2 || !length(paralog_set)) return(FALSE)
  keys <- apply(combn(gs, 2), 2, function(z) paste(sort(z), collapse = "\t"))
  any(keys %in% paralog_set)
}

has_shared_abc <- function(gs) {
  vals <- abc_index[intersect(gs, names(abc_index))]
  if (!length(vals)) return(FALSE)
  enhancers <- unlist(vals, use.names = FALSE)
  any(duplicated(enhancers))
}

gene_pair_features <- function(gs) {
  gm <- gene_model[gs, nomatch = 0]
  if (nrow(gm) < 2) {
    return(list(
      shared_same_strand_promoter = FALSE,
      shared_opposite_strand_promoter = FALSE,
      same_strand_overlap_raw = FALSE,
      opposite_strand_overlap_raw = FALSE,
      same_strand_overlap = FALSE,
      opposite_strand_overlap = FALSE
    ))
  }
  same_prom <- FALSE
  opp_prom <- FALSE
  same_ov <- FALSE
  opp_ov <- FALSE
  for (i in seq_len(nrow(gm) - 1L)) {
    for (j in (i + 1L):nrow(gm)) {
      if (gm$chr[i] != gm$chr[j]) next
      strand_same <- !is.na(gm$strand[i]) && !is.na(gm$strand[j]) && gm$strand[i] == gm$strand[j]
      tss_i <- if (!is.na(gm$strand[i]) && gm$strand[i] == "-") gm$end[i] else gm$start[i]
      tss_j <- if (!is.na(gm$strand[j]) && gm$strand[j] == "-") gm$end[j] else gm$start[j]
      promoters_touch <- abs(tss_i - tss_j) <= PROMOTER_WINDOW_BP
      overlaps <- max(gm$start[i], gm$start[j]) <= min(gm$end[i], gm$end[j])
      if (promoters_touch && strand_same) same_prom <- TRUE
      if (promoters_touch && !strand_same) opp_prom <- TRUE
      if (overlaps && strand_same) same_ov <- TRUE
      if (overlaps && !strand_same) opp_ov <- TRUE
    }
  }
  shared_promoter <- same_prom || opp_prom
  list(
    shared_same_strand_promoter = same_prom,
    shared_opposite_strand_promoter = opp_prom,
    same_strand_overlap_raw = same_ov,
    opposite_strand_overlap_raw = opp_ov,
    same_strand_overlap = same_ov && !shared_promoter,
    opposite_strand_overlap = opp_ov && !shared_promoter
  )
}

cross_boundary <- function(chr, start, end, pos_index, expand = 0L) {
  v <- pos_index[[std_chr(chr)]]
  if (is.null(v) || !length(v)) return(FALSE)
  lower <- as.integer(start) + as.integer(expand)
  upper <- as.integer(end) - as.integer(expand)
  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) return(FALSE)
  idx <- findInterval(lower, v) + 1L
  idx <= length(v) && v[idx] < upper
}

precompute_spearman_types <- function(real_sets) {
  if (!nrow(real_sets)) {
    return(data.table(method = character(), celltype = character(), cluster_id = character(), correlation_type = character()))
  }
  out <- list()
  oi <- 0L
  for (ct in unique(real_sets$celltype)) {
    pair_file <- file.path(ROOT_DIR, "results", "pb_spearman", ct, "pb_spearman_local_pairs.tsv.gz")
    ct_sets <- real_sets[celltype == ct]
    if (!file.exists(pair_file)) {
      out[[oi <- oi + 1L]] <- ct_sets[, .(method, celltype, cluster_id, correlation_type = "unknown")]
      next
    }
    pairs <- fread_maybe_gz(pair_file, select = c("Gene1", "Gene2", "rho"))
    pairs[, key := paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "\t")]
    rho_map <- setNames(pairs$rho, pairs$key)
    ct_out <- vector("list", nrow(ct_sets))
    for (i in seq_len(nrow(ct_sets))) {
      gs <- unlist(strsplit(ct_sets$genes[i], ",", fixed = TRUE))
      keys <- apply(combn(gs, 2), 2, function(z) paste(sort(z), collapse = "\t"))
      vals <- unname(rho_map[keys])
      vals <- vals[is.finite(vals)]
      ctype <- if (!length(vals)) {
        "unknown"
      } else if (all(vals > 0)) {
        "positive_only"
      } else if (all(vals < 0)) {
        "negative_only"
      } else {
        "mixed"
      }
      ct_out[[i]] <- data.table(
        method = ct_sets$method[i],
        celltype = ct,
        cluster_id = ct_sets$cluster_id[i],
        correlation_type = ctype
      )
    }
    out[[oi <- oi + 1L]] <- rbindlist(ct_out)
    message("Precomputed donor-level Spearman cluster types for ", ct, ": ", nrow(ct_sets))
  }
  rbindlist(out, fill = TRUE)
}

spearman_type_map <- precompute_spearman_types(sets[is_correlated_cluster == TRUE & method == "add_cov_sc_hurdle"])
if (nrow(spearman_type_map)) setkey(spearman_type_map, method, celltype, cluster_id)

rows <- vector("list", nrow(sets))
for (i in seq_len(nrow(sets))) {
  s <- sets[i]
  gs <- unlist(strsplit(s$genes, ",", fixed = TRUE))
  pf <- gene_pair_features(gs)
  ctype <- if (!isTRUE(s$is_correlated_cluster)) {
    "null"
  } else if (s$method == "add_cov_sc_hurdle") {
    spearman_type_map[.(s$method, s$celltype, s$cluster_id), correlation_type]
  } else {
    "unknown"
  }
  if (!length(ctype) || is.na(ctype)) ctype <- "unknown"
  rows[[i]] <- data.table(
    method = s$method,
    celltype = s$celltype,
    cluster_id = s$cluster_id,
    set_id = s$set_id,
    chr = s$chr,
    num_genes = s$num_genes,
    cluster_start = s$cluster_start,
    cluster_end = s$cluster_end,
    cluster_length = s$cluster_length,
    log_cluster_length = s$log_cluster_length,
    genes = s$genes,
    is_correlated_cluster = as.logical(s$is_correlated_cluster),
    correlation_type = ctype,
    has_paralog = has_paralog(gs),
    has_shared_go_bp = has_shared_go(gs),
    has_shared_abc_enhancer = has_shared_abc(gs),
    has_shared_same_strand_promoter = pf$shared_same_strand_promoter,
    has_shared_opposite_strand_promoter = pf$shared_opposite_strand_promoter,
    has_same_strand_overlap_raw = pf$same_strand_overlap_raw,
    has_opposite_strand_overlap_raw = pf$opposite_strand_overlap_raw,
    has_same_strand_overlap = pf$same_strand_overlap,
    has_opposite_strand_overlap = pf$opposite_strand_overlap,
    cross_ctcf_peak = cross_boundary(s$chr, s$cluster_start, s$cluster_end, ctcf_pos, 0L),
    cross_tad_boundary = cross_boundary(s$chr, s$cluster_start, s$cluster_end, tad_pos, BOUNDARY_WINDOW_BP)
  )
  if (i %% 10000 == 0) message("Annotated sets: ", i, " / ", nrow(sets))
}

ann <- rbindlist(rows, fill = TRUE)
fwrite(ann, file.path(out_dir, "cluster_annotation_matrix.tsv"), sep = "\t")
fwrite(ann[, .N, by = .(method, is_correlated_cluster, correlation_type)], file.path(out_dir, "annotation_set_counts.tsv"), sep = "\t")
message("Wrote annotation matrix: ", nrow(ann), " rows")
