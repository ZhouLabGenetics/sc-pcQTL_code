#!/usr/bin/env Rscript

# =============================================================================
# Method 2 - Step 2: Calculate SC hurdle associations (CHUNKED VERSION)
# =============================================================================
# Process a subset of genes on one chromosome to keep runtime under 6 hours
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(fasthurdle)
})

require_fasthurdle_version <- function(min_version = "1.1.1") {
  installed <- tryCatch(utils::packageVersion("fasthurdle"), error = function(e) NULL)
  if (is.null(installed)) {
    stop("Package 'fasthurdle' is not installed.")
  }
  if (installed < package_version(min_version)) {
    stop(sprintf(
      "fasthurdle >= %s is required (found %s). Please update fasthurdle.",
      min_version, as.character(installed)
    ))
  }
  as.character(installed)
}

get_this_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_idx <- grep("--file=", cmd_args)
  if (length(file_idx) > 0) {
    return(normalizePath(sub("--file=", "", cmd_args[file_idx])))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
  normalizePath(".")
}

script_dir <- dirname(get_this_path())
source(file.path(script_dir, "load_config.R"))
cfg <- load_config(script_dir)
fasthurdle_version <- require_fasthurdle_version("1.1.1")

option_list <- list(
  make_option("--chr", type = "integer", help = "Chromosome number"),
  make_option("--gene_start", type = "integer", help = "Start gene index (1-based)"),
  make_option("--gene_end", type = "integer", help = "End gene index (inclusive)"),
  make_option("--chunk_id", type = "integer", default = 1, help = "Chunk ID for naming"),
  make_option("--p_threshold", type = "double", default = 0.05),
  make_option("--nonzero_cutoff", type = "double", default = 0.01)
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$chr) || is.null(opt$gene_start) || is.null(opt$gene_end)) {
  stop("Please provide --chr, --gene_start, and --gene_end")
}

CHR <- opt$chr
GENE_START <- opt$gene_start
GENE_END <- opt$gene_end
CHUNK_ID <- opt$chunk_id
P_THRESHOLD <- opt$p_threshold
NONZERO_CUTOFF <- opt$nonzero_cutoff

cat(sprintf("\n=== Method 2 Step 2: SC Hurdle (Chunked) - Chr%d Chunk%d (%s) ===\n",
            CHR, CHUNK_ID, cfg$cell_type))
cat(sprintf("fasthurdle version: %s\n", fasthurdle_version))
cat(sprintf("Gene range: %d-%d\n", GENE_START, GENE_END))
cat(sprintf("Start time: %s\n\n", Sys.time()))

COUNT_FILE <- cfg$count_file
GENE_INFO_FILE <- cfg$gene_info_file
FILTER_FILE <- file.path(cfg$method2_results, "filtered_genes.tsv")

ASSOC_DIR <- cfg$chunk_dir
dir.create(ASSOC_DIR, recursive = TRUE, showWarnings = FALSE)

metadata_cols <- c(
  "CellID", "barcode", "IndividualID", "individual", "CellType",
  "sex", sprintf("pc%d", 1:6), "age", "pf1", "pf2",
  "total_read_counts", "log_total_read_counts"
)

load_library_size <- function(count_file, header_cols, id_col_name, cache_file) {
  available <- intersect(c(id_col_name, "total_read_counts", "log_total_read_counts"), header_cols)
  if (all(c("total_read_counts", "log_total_read_counts") %in% available)) {
    return(fread(count_file, select = available))
  }
  if (file.exists(cache_file)) {
    return(fread(cache_file))
  }

  cat("Computing library sizes from full count matrix...\n")
  full_dt <- fread(count_file)
  gene_cols <- setdiff(names(full_dt), metadata_cols)
  lib_dt <- full_dt[, ..id_col_name]
  lib_dt[, total_read_counts := rowSums(as.matrix(full_dt[, ..gene_cols]))]
  lib_dt[, log_total_read_counts := log(pmax(total_read_counts, 1))]
  fwrite(lib_dt, cache_file, sep = "\t", compress = "gzip")
  lib_dt
}

# Load gene metadata
cat("Loading gene metadata...\n")
gene_info <- fread(GENE_INFO_FILE)
filt <- fread(FILTER_FILE)
merged <- merge(gene_info, filt[keep == TRUE, .(gene_name, nonzero_prop)], by = "gene_name")
merged <- merged[!duplicated(gene_name)]
global_counts <- merged[, .N, by = chr_numeric]
global_counts <- global_counts[chr_numeric %in% seq_len(22)]
global_pairs <- global_counts[, sum(N * (N - 1) / 2)]
bonferroni_threshold <- P_THRESHOLD / global_pairs

chr_genes_all <- merged[chr_numeric == CHR][order(start), gene_name]

n_genes_chr <- length(chr_genes_all)
cat(sprintf("  Chromosome %d total filtered genes: %d\n", CHR, n_genes_chr))

# Select chunk
if (GENE_END > n_genes_chr) GENE_END <- n_genes_chr
chr_genes <- chr_genes_all[GENE_START:GENE_END]
n_genes <- length(chr_genes)

cat(sprintf("  Selected genes %d-%d: %d genes\n", GENE_START, GENE_END, n_genes))

if (n_genes < 2) {
  cat("  Not enough genes; skipping.\n")
  quit(status = 0)
}

# Load counts
cat("Loading single-cell counts...\n")
header_cols <- names(fread(COUNT_FILE, nrows = 0))
id_col <- intersect(c("CellID", "barcode"), header_cols)
if (!length(id_col)) stop("Count file must contain CellID or barcode column")
id_col_name <- id_col[1]
covariate_cols <- c("sex", sprintf("pc%d", 1:6), "age", "pf1", "pf2",
                    "total_read_counts", "log_total_read_counts")
select_cols <- unique(c(id_col_name, "individual", covariate_cols, chr_genes))
select_cols <- intersect(select_cols, header_cols)
count_dt <- fread(COUNT_FILE, select = select_cols)
count_mat <- as.matrix(count_dt[, ..chr_genes])
n_cells <- nrow(count_mat)

cov_dt <- copy(count_dt[, !..chr_genes])
if (!all(c("total_read_counts", "log_total_read_counts") %in% names(cov_dt))) {
  lib_cache <- file.path(cfg$method2_results, "library_size_cache.tsv.gz")
  lib_dt <- load_library_size(COUNT_FILE, header_cols, id_col_name, lib_cache)
  cov_dt <- merge(cov_dt, lib_dt, by = id_col_name, all.x = TRUE, sort = FALSE, suffixes = c("", ".cache"))
  if ("total_read_counts.cache" %in% names(cov_dt)) {
    cov_dt[is.na(total_read_counts), total_read_counts := total_read_counts.cache]
    cov_dt[, total_read_counts.cache := NULL]
  }
  if ("log_total_read_counts.cache" %in% names(cov_dt)) {
    cov_dt[is.na(log_total_read_counts), log_total_read_counts := log_total_read_counts.cache]
    cov_dt[, log_total_read_counts.cache := NULL]
  }
}
cov_dt[, age := as.numeric(age)]
for (pc in sprintf("pc%d", 1:6)) {
  if (pc %in% names(cov_dt)) cov_dt[[pc]] <- as.numeric(cov_dt[[pc]])
}
for (pf in c("pf1", "pf2")) {
  if (pf %in% names(cov_dt)) cov_dt[[pf]] <- as.numeric(cov_dt[[pf]])
}
if ("sex" %in% names(cov_dt)) {
  cov_dt[, sex := as.factor(sex)]
}
required_covars <- c("age", "sex", sprintf("pc%d", 1:6), "pf1", "pf2", "log_total_read_counts")
missing_covars <- setdiff(required_covars, names(cov_dt))
if (length(missing_covars)) {
  stop("Missing required covariate columns: ", paste(missing_covars, collapse = ", "))
}

# Calculate associations
pair_idx <- combn(n_genes, 2, simplify = FALSE)
n_pairs <- length(pair_idx)
cat(sprintf("  Cells: %d, Genes: %d, Pairs: %d\n", n_cells, n_genes, n_pairs))

run_pair <- function(i, j) {
  df <- data.frame(
    count_i = count_mat[, i],
    count_j = count_mat[, j],
    age = cov_dt$age,
    sex = cov_dt$sex,
    pc1 = cov_dt$pc1,
    pc2 = cov_dt$pc2,
    pc3 = cov_dt$pc3,
    pc4 = cov_dt$pc4,
    pc5 = cov_dt$pc5,
    pc6 = cov_dt$pc6,
    pf1 = cov_dt$pf1,
    pf2 = cov_dt$pf2,
    log_total_read_counts = cov_dt$log_total_read_counts
  )
  if (all(df$count_i == 0) || all(df$count_j == 0)) {
    return(list(p_count = NA_real_, p_zero = NA_real_))
  }

  get_active_covariates <- function(data, variables) {
    variables[vapply(variables, function(nm) {
      vals <- data[[nm]]
      length(unique(vals[!is.na(vals)])) > 1
    }, logical(1))]
  }

  fit_one_direction <- function(response, predictor) {
    base_covars <- get_active_covariates(
      df,
      c("age", "sex", sprintf("pc%d", 1:6), "pf1", "pf2")
    )
    zero_covars <- get_active_covariates(df, "log_total_read_counts")
    count_rhs <- paste(c(predictor, base_covars), collapse = " + ")
    zero_rhs <- paste(c(predictor, zero_covars, base_covars), collapse = " + ")
    model_formula <- as.formula(
      sprintf(
        "%s ~ %s + offset(log_total_read_counts) | %s",
        response,
        count_rhs,
        zero_rhs
      )
    )
    fit <- tryCatch(
      fasthurdle(
        model_formula,
        data = df,
        dist = "poisson",
        zero.dist = "binomial"
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(list(p_count = NA_real_, p_zero = NA_real_))
    }

    smry <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    if (is.null(smry)) {
      return(list(p_count = NA_real_, p_zero = NA_real_))
    }
    p_from_summary <- function(component_name, term) {
      component <- smry[[component_name]]
      if (is.null(component) || !term %in% rownames(component)) {
        return(NA_real_)
      }
      component[term, "Pr(>|z|)"]
    }
    list(
      p_count = p_from_summary("count", predictor),
      p_zero = p_from_summary("zero", predictor)
    )
  }

  fit_ij <- fit_one_direction("count_i", "count_j")
  fit_ji <- fit_one_direction("count_j", "count_i")
  p_count_vals <- c(fit_ij$p_count, fit_ji$p_count)
  p_zero_vals <- c(fit_ij$p_zero, fit_ji$p_zero)

  p_count <- min(p_count_vals, na.rm = TRUE)
  if (!is.finite(p_count)) p_count <- NA_real_
  p_zero <- min(p_zero_vals, na.rm = TRUE)
  if (!is.finite(p_zero)) p_zero <- NA_real_

  list(p_count = p_count, p_zero = p_zero)
}

pval_count <- matrix(NA_real_, n_genes, n_genes, dimnames = list(chr_genes, chr_genes))
pval_zero <- matrix(NA_real_, n_genes, n_genes, dimnames = list(chr_genes, chr_genes))

cat("Calculating associations...\n")
start_time <- Sys.time()

for (idx in seq_along(pair_idx)) {
  pair <- pair_idx[[idx]]
  i <- pair[1]; j <- pair[2]
  res <- run_pair(i, j)
  pval_count[i, j] <- pval_count[j, i] <- res$p_count
  pval_zero[i, j] <- pval_zero[j, i] <- res$p_zero

  if (idx %% 1000 == 0 || idx == n_pairs) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    rate <- idx / max(elapsed, 1e-6)
    eta <- (n_pairs - idx) / max(rate, 1e-6)
    cat(sprintf("  %d/%d pairs (%.1f%%) -- %.2f pairs/s, ETA %.1f min\n",
                idx, n_pairs, 100 * idx / n_pairs, rate, eta / 60))
  }
}

# Extract significant pairs
sig_pairs <- data.table()
for (i in 1:(n_genes - 1)) {
  for (j in (i + 1):n_genes) {
    p_count <- pval_count[i, j]
    p_zero <- pval_zero[i, j]
    # Check for NA before comparison to avoid "missing value where TRUE/FALSE needed"
    if ((!is.na(p_count) && p_count < bonferroni_threshold) ||
        (!is.na(p_zero) && p_zero < bonferroni_threshold)) {
      sig_pairs <- rbind(sig_pairs, data.table(
        Gene1 = chr_genes[i],
        Gene2 = chr_genes[j],
        Pvalue_count = p_count,
        Pvalue_zero = p_zero
      ))
    }
  }
}

# Save results
out_base <- file.path(ASSOC_DIR, sprintf("chr%d_chunk%03d", CHR, CHUNK_ID))
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

fwrite(as.data.table(pval_count, keep.rownames = "Gene"),
       file.path(out_base, "pvalues_count.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(as.data.table(pval_zero, keep.rownames = "Gene"),
       file.path(out_base, "pvalues_zero.tsv.gz"), sep = "\t", compress = "gzip")

if (nrow(sig_pairs) > 0) {
  fwrite(sig_pairs, file.path(out_base, "significant_pairs.tsv"), sep = "\t")
}

summary_dt <- data.table(
  Chromosome = CHR,
  ChunkID = CHUNK_ID,
  GeneStart = GENE_START,
  GeneEnd = GENE_END,
  NumGenes = n_genes,
  NumPairs = n_pairs,
  GlobalPairs = global_pairs,
  BonferroniThreshold = bonferroni_threshold,
  SignificantPairs = nrow(sig_pairs)
)
fwrite(summary_dt, file.path(out_base, "summary.tsv"), sep = "\t")

cat("\n=== Complete ===\n")
cat(sprintf("End time: %s\n", Sys.time()))
cat(sprintf("Output: %s\n", out_base))
