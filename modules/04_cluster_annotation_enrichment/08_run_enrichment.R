#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))
source("config/config.R")

ann_file <- file.path(ROOT_DIR, "results", "annotations", "cluster_annotation_matrix.tsv")
out_dir <- file.path(ROOT_DIR, "results", "enrichment")
dir_create(out_dir)
ann <- fread(ann_file)
ann <- ann[celltype %chin% CELLTYPES]
if (!nrow(ann)) stop("No annotation rows remain after applying the primary cell-type eligibility set.")

features <- c(
  "has_paralog", "has_shared_go_bp", "has_shared_abc_enhancer",
  "has_shared_same_strand_promoter", "has_shared_opposite_strand_promoter",
  "has_same_strand_overlap", "has_opposite_strand_overlap",
  "cross_ctcf_peak", "cross_tad_boundary"
)

expected_freq_qc <- function(y, group, min_expected = 1) {
  tab <- table(factor(y, levels = c(0, 1)), factor(group, levels = c(FALSE, TRUE)))
  row_tot <- rowSums(tab)
  col_tot <- colSums(tab)
  total <- sum(tab)
  if (!total || any(row_tot == 0) || any(col_tot == 0)) {
    return(list(pass = FALSE, min_expected = 0, tab = tab, reason = "complete_or_empty_table"))
  }
  expected <- outer(row_tot, col_tot) / total
  min_exp <- min(expected)
  list(pass = is.finite(min_exp) && min_exp >= min_expected, min_expected = min_exp, tab = tab, reason = "ok")
}

fit_one <- function(dt, feature, method, stratum) {
  dt <- copy(dt)
  dt[, y := as.integer(get(feature))]
  qc <- expected_freq_qc(dt$y, dt$is_correlated_cluster, min_expected = 1)
  if (!qc$pass) {
    return(data.table(
      method = method,
      stratum = stratum,
      annotation = feature,
      status = "skip_expected_lt_1",
      min_expected = qc$min_expected,
      n_sets = nrow(dt),
      n_correlated = sum(dt$is_correlated_cluster),
      n_null = sum(!dt$is_correlated_cluster),
      annotated_correlated = sum(dt$y == 1 & dt$is_correlated_cluster),
      annotated_null = sum(dt$y == 1 & !dt$is_correlated_cluster)
    ))
  }
  fit <- tryCatch(glm(y ~ is_correlated_cluster + num_genes + log_cluster_length, data = dt, family = binomial()), error = function(e) e)
  if (inherits(fit, "error")) {
    return(data.table(method = method, stratum = stratum, annotation = feature, status = paste0("glm_error:", fit$message)))
  }
  cf <- summary(fit)$coefficients
  term <- "is_correlated_clusterTRUE"
  if (!term %in% rownames(cf)) {
    return(data.table(method = method, stratum = stratum, annotation = feature, status = "missing_term"))
  }
  est <- cf[term, "Estimate"]
  se <- cf[term, "Std. Error"]
  data.table(
    method = method,
    stratum = stratum,
    annotation = feature,
    status = "ok",
    beta = est,
    se = se,
    OR = exp(est),
    CI_low = exp(est - 1.96 * se),
    CI_high = exp(est + 1.96 * se),
    pvalue = cf[term, "Pr(>|z|)"],
    min_expected = qc$min_expected,
    n_sets = nrow(dt),
    n_correlated = sum(dt$is_correlated_cluster),
    n_null = sum(!dt$is_correlated_cluster),
    annotated_correlated = sum(dt$y == 1 & dt$is_correlated_cluster),
    annotated_null = sum(dt$y == 1 & !dt$is_correlated_cluster)
  )
}

results <- list()
ri <- 0L
set.seed(20240701L)
for (method_name in "add_cov_sc_hurdle") {
  null_all <- ann[method == method_name & is_correlated_cluster == FALSE]
  corr_all <- ann[method == method_name & is_correlated_cluster == TRUE]
  # pcQTL paper (Lawrence & Montgomery) method: null clusters "sampled at a rate to
  # match the relative distributions of 2, 3, 4, or 5 gene clusters in the correlated
  # cluster set" -- downsample the null so its num_genes distribution matches the
  # correlated set's, maximizing total N subject to availability (seeded).
  sizes <- sort(unique(corr_all$num_genes))
  p_size <- prop.table(table(factor(corr_all$num_genes, levels = sizes)))
  avail <- table(factor(null_all$num_genes, levels = sizes))
  Ntot <- floor(min(as.numeric(avail[names(p_size)]) / as.numeric(p_size)))
  n_take <- pmin(as.integer(round(Ntot * as.numeric(p_size))), as.integer(avail[names(p_size)]))
  names(n_take) <- names(p_size)
  null_dt <- rbindlist(lapply(sizes, function(k) {
    pool <- null_all[num_genes == k]
    take <- n_take[as.character(k)]
    if (is.na(take) || take <= 0L || !nrow(pool)) return(pool[0])
    pool[sample(.N, min(.N, take))]
  }))
  message(method_name, ": size-matched null n=", nrow(null_dt), " (of ", nrow(null_all), " available)")
  for (stratum in c("all_correlated", "positive_only", "negative_only", "mixed")) {
    real_dt <- if (stratum == "all_correlated") {
      ann[method == method_name & is_correlated_cluster == TRUE]
    } else {
      ann[method == method_name & is_correlated_cluster == TRUE & correlation_type == stratum]
    }
    dt <- rbind(real_dt, null_dt, fill = TRUE)
    if (!nrow(real_dt) || !nrow(null_dt)) next
    for (feature in features) {
      ri <- ri + 1L
      results[[ri]] <- fit_one(dt, feature, method_name, stratum)
    }
  }
}
res <- rbindlist(results, fill = TRUE)
res[status == "ok", fdr := p.adjust(pvalue, method = "BH")]
fwrite(res, file.path(out_dir, "enrichment_by_method.tsv"), sep = "\t")

print(res[status == "ok"][order(fdr)][1:min(.N, 20)])
