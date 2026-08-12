#!/usr/bin/env Rscript
# ============================================================================
# 01_collect_results.R
# Collect pcQTL results from primary-analysis cell types
# ============================================================================

.local_libs <- c(
  Sys.getenv("SC_PCQTL_R_LIBS", unset = ""),
  Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = "")
)
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# Setup
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  normalizePath(getwd())
}
MODULE_DIR <- get_script_dir()
RELEASE_ROOT <- normalizePath(file.path(MODULE_DIR, "..", "..", ".."), mustWork = TRUE)
source(file.path(RELEASE_ROOT, "config", "celltype_eligibility.R"))
WORKFLOW_ROOT <- Sys.getenv("SC_PCQTL_WORKFLOW_ROOT", unset = Sys.getenv("COQTL_WORKFLOW_ROOT", unset = ""))
if (!nzchar(WORKFLOW_ROOT) || !dir.exists(WORKFLOW_ROOT)) {
  stop("Set SC_PCQTL_WORKFLOW_ROOT or COQTL_WORKFLOW_ROOT to the completed workflow root.")
}
WORKFLOW_ROOT <- normalizePath(WORKFLOW_ROOT, mustWork = TRUE)
UPSTREAM_CELLTYPES_DIR <- normalizePath(
  Sys.getenv("COQTL_UPSTREAM_CELLTYPES_DIR",
             unset = file.path(WORKFLOW_ROOT, "03_analysis_celltypes", "01_upstream_main_pipeline_add_cov", "celltypes"))
)
SUMMARY_ROOT <- Sys.getenv(
  "SC_PCQTL_PCQTL_SUMMARY_ROOT",
  unset = file.path(
    WORKFLOW_ROOT, "03_analysis_celltypes",
    "02_downstream_analysis_modules_add_cov_fdr", "pcqtl_compare"
  )
)
OUT_DIR <- file.path(SUMMARY_ROOT, "data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CELLTYPES <- primary_celltypes(file.path(RELEASE_ROOT, "config", "celltype_eligibility.tsv"))

cat("=== Collecting pcQTL results ===\n")
cat("Celltypes:", length(CELLTYPES), "\n\n")

# Function to parse cluster_id
parse_cluster_id <- function(cluster_id) {
  # SC_chr10_cluster_001 -> chr=10, cluster_num=1
  parts <- str_match(cluster_id, "SC_chr(\\d+|X|Y)_cluster_(\\d+)")
  data.frame(
    chr = parts[,2],
    cluster_num = as.integer(parts[,3])
  )
}

# Function to parse MarkerID
parse_marker_id <- function(marker_id) {
  # "10:73582752" -> chr=10, pos=73582752
  parts <- str_split_fixed(marker_id, ":", 2)
  data.frame(
    snp_chr = parts[,1],
    snp_pos = as.numeric(parts[,2])
  )
}

# Helper: locate step2 SAIGE output for a given cluster/PC
locate_step2_file <- function(base_dir, cluster_id, pc) {
  candidates <- c(
    file.path(base_dir, cluster_id, pc),
    file.path(base_dir, cluster_id, paste0(pc, ".txt")),
    file.path(base_dir, cluster_id, paste0(pc, ".gz"))
  )
  exists <- file.exists(candidates)
  if (any(exists)) {
    return(candidates[which(exists)[1]])
  }
  NULL
}

# Helper: read step2 output and summarise SNP-level stats
summarise_step2 <- function(step2_file, fdr_threshold = 0.05) {
  if (is.null(step2_file) || !file.exists(step2_file)) {
    return(data.table(
      n_snps = NA_integer_,
      n_sig_snps = NA_integer_,
      min_snp_p = NA_real_,
      min_snp_fdr = NA_real_,
      lead_snp_marker = NA_character_,
      lead_snp_chr = NA_character_,
      lead_snp_pos = NA_real_,
      has_sig_snp = NA
    ))
  }

  dt <- tryCatch({
    fread(step2_file)
  }, error = function(e) {
    warning("Failed to read step2 file: ", step2_file, " (", e$message, ")")
    return(NULL)
  })

  if (is.null(dt) || nrow(dt) == 0) {
    return(data.table(
      n_snps = 0L,
      n_sig_snps = 0L,
      min_snp_p = NA_real_,
      min_snp_fdr = NA_real_,
      lead_snp_marker = NA_character_,
      lead_snp_chr = NA_character_,
      lead_snp_pos = NA_real_,
      has_sig_snp = FALSE
    ))
  }

  # Identify columns
  p_col <- intersect(c("p.value", "pvalue", "Pvalue", "P.value", "pval"), names(dt))
  marker_col <- intersect(c("MarkerID", "markerID", "SNPID", "rsID"), names(dt))
  chr_col <- intersect(c("CHR", "CHROM", "chromosome"), names(dt))
  pos_col <- intersect(c("POS", "Position", "BP"), names(dt))

  if (length(p_col) == 0) {
    warning("p-value column not found in ", step2_file)
    return(data.table(
      n_snps = nrow(dt),
      n_sig_snps = NA_integer_,
      min_snp_p = NA_real_,
      min_snp_fdr = NA_real_,
      lead_snp_marker = NA_character_,
      lead_snp_chr = NA_character_,
      lead_snp_pos = NA_real_,
      has_sig_snp = NA
    ))
  }

  pvals <- dt[[p_col[1]]]
  adj <- p.adjust(pvals, method = "BH")
  min_idx <- which.min(pvals)

  marker_val <- if (length(marker_col)) dt[[marker_col[1]]][min_idx] else NA_character_
  chr_val <- if (length(chr_col)) dt[[chr_col[1]]][min_idx] else NA_character_
  pos_val <- if (length(pos_col)) dt[[pos_col[1]]][min_idx] else NA_real_

  data.table(
    n_snps = nrow(dt),
    n_sig_snps = sum(adj < fdr_threshold),
    min_snp_p = pvals[min_idx],
    min_snp_fdr = adj[min_idx],
    lead_snp_marker = marker_val,
    lead_snp_chr = as.character(chr_val),
    lead_snp_pos = as.numeric(pos_val),
    has_sig_snp = sum(adj < fdr_threshold) > 0
  )
}

# Collect results from eligible cell types
all_results <- list()

for (ct in CELLTYPES) {
  cat("Processing:", ct, "...\n")

  # Paths
  pcqtl_dir <- file.path(UPSTREAM_CELLTYPES_DIR, ct, "pcQTL")
  step3_dir <- file.path(pcqtl_dir, "step3_saige/step3")
  step2_dir <- file.path(pcqtl_dir, "step3_saige/step2")
  cluster_file <- file.path(UPSTREAM_CELLTYPES_DIR, ct, "cluster_identification/results/method2_sc_hurdle/merged_clusters/cluster_summary.tsv")

  # Read cluster info
  if (!file.exists(cluster_file)) {
    cat("  WARNING: cluster_summary.tsv not found, skipping\n")
    next
  }

  cluster_info <- fread(cluster_file, header = TRUE)
  setnames(cluster_info,
           c("cluster_id", "chromosome", "cluster_size", "start_position",
             "end_position", "cluster_span_bp", "genes"))

  # Find all genePval files
  genepval_files <- list.files(step3_dir, pattern = "_genePval$",
                                recursive = TRUE, full.names = TRUE)

  if (length(genepval_files) == 0) {
    cat("  WARNING: no genePval files found, skipping\n")
    next
  }

  cat("  Found", length(genepval_files), "genePval files\n")

  # Read all genePval files
  results <- lapply(genepval_files, function(f) {
    # Parse file path: .../SC_chr10_cluster_001/PC1_genePval
    path_parts <- str_split(f, "/")[[1]]
    cluster_id <- path_parts[length(path_parts) - 1]
    pc_file <- path_parts[length(path_parts)]
    pc <- str_replace(pc_file, "_genePval$", "")

    # Read genePval
    d <- tryCatch({
      fread(f, header = TRUE)
    }, error = function(e) {
      cat("    ERROR reading", f, ":", e$message, "\n")
      NULL
    })
    if (is.null(d) || nrow(d) == 0) return(NULL)

    d$cluster_id <- cluster_id
    d$PC <- pc
    d$celltype <- ct

    # attach SNP-level summary from step2
    step2_file <- locate_step2_file(step2_dir, cluster_id, pc)
    snp_summary <- summarise_step2(step2_file)
    cbind(d, snp_summary)
  })

  # Combine
  results <- rbindlist(results[!sapply(results, is.null)], fill = TRUE)

  if (nrow(results) == 0) {
    cat("  WARNING: no valid results, skipping\n")
    next
  }

  # Parse cluster_id and marker_id
  cluster_parsed <- parse_cluster_id(results$cluster_id)
  marker_parsed <- parse_marker_id(results$top_MarkerID)

  results <- cbind(results, cluster_parsed, marker_parsed)

  # Merge with cluster info
  # Note: cluster_id is already in the correct format (SC_chr..._cluster_...)
  results <- merge(results,
                   cluster_info[, .(cluster_id, start_position, end_position,
                                   cluster_span_bp, genes)],
                   by = "cluster_id", all.x = TRUE)

  cat("  Collected", nrow(results), "results\n")
  all_results[[ct]] <- results
}

# Combine eligible cell types
cat("\n=== Combining primary-analysis cell types ===\n")
if (length(all_results) == 0) {
  stop(
    "No pcQTL genePval files were found under the configured upstream add_cov directories. ",
    "Check pcQTL config/output paths before rerunning pcqtl_compare."
  )
}
all_data <- rbindlist(all_results, fill = TRUE)

# Reorder columns
all_data <- all_data[, .(
  celltype, cluster_id, PC,
  chr, start_position, end_position, cluster_span_bp, genes,
  ACAT_p, top_MarkerID, top_pval,
  snp_chr, snp_pos,
  n_snps, n_sig_snps, min_snp_p, min_snp_fdr,
  lead_snp_marker, lead_snp_chr, lead_snp_pos, has_sig_snp
)]

# Sort
all_data <- all_data[order(celltype, chr, start_position, PC)]

# Summary
cat("\nTotal results:", nrow(all_data), "\n")
cat("By celltype:\n")
print(table(all_data$celltype))

# Save
out_file <- file.path(OUT_DIR, "all_pcqtl_results.tsv")
fwrite(all_data, out_file, sep = "\t")
cat("\n✓ Saved to:", out_file, "\n")

# Summary stats
cat("\n=== Summary Statistics ===\n")
cat("Celltypes:", uniqueN(all_data$celltype), "\n")
cat("Unique clusters:", uniqueN(all_data$cluster_id), "\n")
cat("Total tests:", nrow(all_data), "\n")
cat("Chromosomes:", paste(sort(unique(all_data$chr)), collapse = ", "), "\n")

cat("\nP-value distribution:\n")
print(summary(all_data$ACAT_p))

cat("\n✓ Step 1 complete!\n")
