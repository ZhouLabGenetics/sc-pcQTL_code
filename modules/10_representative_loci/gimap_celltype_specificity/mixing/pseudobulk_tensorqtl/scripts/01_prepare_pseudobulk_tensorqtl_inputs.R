#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
release_root <- normalizePath(file.path(dirname(script_arg), "..", "..", "..", "..", "..", ".."), mustWork = TRUE)
source(file.path(release_root, "config", "celltype_eligibility.R"))

require_path <- function(name, must_exist = TRUE) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value) || (must_exist && !dir.exists(value))) {
    stop("Set ", name, if (must_exist) " to an existing directory." else ".")
  }
  normalizePath(value, mustWork = must_exist)
}

base_dir <- require_path("SC_PCQTL_GIMAP_MIXING_ROOT", must_exist = FALSE)
cache_dir <- require_path("SC_PCQTL_GIMAP_SELECTED_READCOUNT_DIR")
for (d in c("phenotypes", "covariates", "data", "genotypes", "coloc_susie", "logs", "tmp")) {
  dir.create(file.path(base_dir, d), recursive = TRUE, showWarnings = FALSE)
}
set.seed(4437)

gimap_genes <- c("GIMAP8", "GIMAP7", "GIMAP4", "GIMAP6", "GIMAP2", "GIMAP1", "GIMAP5")
gene_id <- c(
  GIMAP2 = "ENSG00000106560.10",
  GIMAP6 = "ENSG00000133561.15",
  GIMAP4 = "ENSG00000133574.9",
  GIMAP8 = "ENSG00000171115.3",
  GIMAP7 = "ENSG00000179144.4",
  GIMAP5 = "ENSG00000196329.11",
  GIMAP1 = "ENSG00000213203.2"
)
# Use a common cluster-union interval for all GIMAP gene and PC phenotypes so TensorQTL tests the same cis variant set.
cluster_chr <- "7"
cluster_start_0based <- 150450628L
cluster_end <- 150737348L
cis_start <- 149450629L
cis_end <- 151737348L
celltypes <- primary_celltypes(file.path(release_root, "config", "celltype_eligibility.tsv"))
# Use only donor-level covariates for synthetic pseudobulk mixtures. PEER factors
# are cell-type-specific in the OneK1K single-cell workflow and are therefore
# not well defined after mixing cell types.
covars <- c("age", "sex", paste0("pc", 1:6))

read_tsv <- function(path) read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "", quote = "")
write_tsv <- function(x, path, row_names = FALSE, col_names = TRUE) write.table(x, path, sep = "\t", quote = FALSE, row.names = row_names, col.names = col_names, na = "")
num <- function(x) suppressWarnings(as.numeric(x))

aggregate_celltype <- function(ct) {
  path <- file.path(cache_dir, paste0(ct, "_gimap_covars.tsv"))
  if (!file.exists(path)) stop("Missing selected readcount cache: ", path)
  message("Reading and aggregating ", ct)
  dat <- read_tsv(path)
  dat$individual <- as.character(dat$individual)
  dat[gimap_genes] <- lapply(dat[gimap_genes], num)
  for (cc in covars) dat[[cc]] <- num(dat[[cc]])
  cell_count <- as.data.frame(table(dat$individual), stringsAsFactors = FALSE)
  names(cell_count) <- c("individual", "n_cells")
  expr <- aggregate(dat[gimap_genes], by = list(individual = dat$individual), FUN = mean, na.rm = TRUE)
  cv <- aggregate(dat[covars], by = list(individual = dat$individual), FUN = function(x) x[which(!is.na(x))[1]])
  out <- Reduce(function(a, b) merge(a, b, by = "individual", all = TRUE), list(expr, cv, cell_count))
  out$celltype <- ct
  out
}

by_ct <- setNames(lapply(celltypes, aggregate_celltype), celltypes)
agg <- do.call(rbind, by_ct)
write_tsv(agg, file.path(base_dir, "data", "gimap_donor_celltype_pseudobulk_expression.tsv"))
cell_counts <- aggregate(n_cells ~ celltype, data = agg, FUN = function(x) c(n_donors = length(x), median_cells = median(x), total_cells = sum(x)))
cell_counts <- do.call(data.frame, cell_counts)
names(cell_counts) <- c("celltype", "n_donors", "median_cells_per_donor", "total_cells")
write_tsv(cell_counts, file.path(base_dir, "data", "gimap_celltype_cell_count_summary.tsv"))

make_expression_matrix <- function(scenario, w_cd8 = NA_real_) {
  if (scenario == "cd8_nc_only") {
    sub <- agg[agg$celltype == "cd8_nc", , drop = FALSE]
    expr <- sub[, c("individual", gimap_genes, covars, "n_cells"), drop = FALSE]
    expr$scenario <- scenario; expr$cd8_nc_weight <- 1
    return(expr)
  }
  donors <- sort(unique(agg$individual))
  rows_out <- lapply(donors, function(id) {
    rows <- agg[agg$individual == id, , drop = FALSE]
    if (nrow(rows) == 0) return(NULL)
    if (scenario == "equal_celltype_weighted") {
      ww <- rep(1 / nrow(rows), nrow(rows))
    } else if (scenario == "observed_cell_count_weighted") {
      ww <- rows$n_cells / sum(rows$n_cells)
    } else if (grepl("^cd8_nc_weight_", scenario)) {
      cd8 <- rows[rows$celltype == "cd8_nc", , drop = FALSE]
      other <- rows[rows$celltype != "cd8_nc", , drop = FALSE]
      if (nrow(cd8) == 0 || nrow(other) == 0) return(NULL)
      other_w <- other$n_cells / sum(other$n_cells)
      expr_other <- colSums(sweep(as.matrix(other[gimap_genes]), 1, other_w, "*"), na.rm = TRUE)
      expr_cd8 <- as.numeric(cd8[1, gimap_genes])
      expr <- w_cd8 * expr_cd8 + (1 - w_cd8) * expr_other
      cv <- rows[1, covars, drop = FALSE]
      out <- data.frame(individual = id, as.list(expr), cv, n_cells = sum(rows$n_cells), stringsAsFactors = FALSE)
      names(out)[seq_len(length(gimap_genes)) + 1] <- gimap_genes
      return(out)
    } else {
      stop("Unknown scenario: ", scenario)
    }
    expr <- colSums(sweep(as.matrix(rows[gimap_genes]), 1, ww, "*"), na.rm = TRUE)
    cv <- rows[1, covars, drop = FALSE]
    out <- data.frame(individual = id, as.list(expr), cv, n_cells = sum(rows$n_cells), stringsAsFactors = FALSE)
    names(out)[seq_len(length(gimap_genes)) + 1] <- gimap_genes
    out
  })
  out <- do.call(rbind, rows_out)
  out$scenario <- scenario; out$cd8_nc_weight <- w_cd8
  out
}

scenarios <- c("cd8_nc_only", "observed_cell_count_weighted", "equal_celltype_weighted")
expr_list <- setNames(lapply(scenarios, make_expression_matrix), scenarios)
# cd8_nc_only is the 100% endpoint; do not create a duplicate weighted scenario.
for (w in c(0, 0.05, 0.10, 0.25, 0.50, 0.75)) {
  nm <- sprintf("cd8_nc_weight_%03d_of_100", as.integer(round(w * 100)))
  expr_list[[nm]] <- make_expression_matrix(nm, w_cd8 = w)
}

format_pheno_bed <- function(expr, scenario) {
  expr <- expr[complete.cases(expr[, c("individual", gimap_genes, covars), drop = FALSE]), , drop = FALSE]
  donors <- expr$individual
  gene_mat <- as.matrix(expr[, gimap_genes, drop = FALSE])
  storage.mode(gene_mat) <- "double"
  pca <- prcomp(gene_mat, center = TRUE, scale. = FALSE)
  pcs <- as.data.frame(pca$x[, seq_len(min(7, ncol(pca$x))), drop = FALSE])
  names(pcs) <- paste0("PC", seq_len(ncol(pcs)))
  if (ncol(pcs) < 7) for (j in (ncol(pcs) + 1):7) pcs[[paste0("PC", j)]] <- NA_real_

  gene_pheno <- data.frame(`#chr` = cluster_chr, start = cluster_start_0based, end = cluster_end,
                           phenotype_id = paste0(scenario, "__e__", gimap_genes), check.names = FALSE)
  gene_values <- as.data.frame(t(as.matrix(expr[, gimap_genes, drop = FALSE])))
  names(gene_values) <- donors
  gene_bed <- cbind(gene_pheno, gene_values)

  pc_pheno <- data.frame(`#chr` = cluster_chr, start = cluster_start_0based, end = cluster_end,
                         phenotype_id = paste0(scenario, "__", paste0("PC", 1:7)), check.names = FALSE)
  pc_values <- as.data.frame(t(as.matrix(pcs[, paste0("PC", 1:7), drop = FALSE])))
  names(pc_values) <- donors
  pc_bed <- cbind(pc_pheno, pc_values)
  cov <- expr[, covars, drop = FALSE]
  cov_t <- as.data.frame(t(as.matrix(cov)))
  names(cov_t) <- donors
  cov_t <- cbind(covariate_id = rownames(cov_t), cov_t)
  rownames(cov_t) <- NULL

  variance <- data.frame(scenario = scenario, PC = paste0("PC", seq_along(pca$sdev)), variance_explained = pca$sdev^2 / sum(pca$sdev^2))
  load <- as.data.frame(pca$rotation[, seq_len(min(7, ncol(pca$rotation))), drop = FALSE])
  names(load) <- paste0("PC", seq_len(ncol(load)))
  load$gene <- rownames(load); load$scenario <- scenario
  load <- load[, c("scenario", "gene", paste0("PC", seq_len(ncol(load) - 2))), drop = FALSE]
  summary <- data.frame(scenario = scenario, n_donors = length(donors), n_cells = sum(expr$n_cells), cd8_nc_weight = unique(expr$cd8_nc_weight)[1], stringsAsFactors = FALSE)
  list(gene_bed = gene_bed, pc_bed = pc_bed, cov = cov_t, variance = variance, loadings = load, summary = summary)
}

all_variance <- list(); all_loadings <- list(); all_summary <- list(); manifest <- list()
for (scenario in names(expr_list)) {
  message("Writing TensorQTL inputs for ", scenario)
  res <- format_pheno_bed(expr_list[[scenario]], scenario)
  scen_dir <- file.path(base_dir, "phenotypes", scenario); dir.create(scen_dir, recursive = TRUE, showWarnings = FALSE)
  cov_file <- file.path(base_dir, "covariates", paste0(scenario, ".covariates.tsv"))
  write_tsv(res$cov, cov_file)
  files <- c(
    single_gene_only = file.path(scen_dir, "single_gene_only.bed"),
    pc_only = file.path(scen_dir, "pc_only.bed")
  )
  write_tsv(res$gene_bed, files[["single_gene_only"]])
  write_tsv(res$pc_bed, files[["pc_only"]])
  for (nm in names(files)) {
    gz <- paste0(files[[nm]], ".gz")
    gzip_status <- system2("gzip", c("-f", "-c", files[[nm]]), stdout = gz)
    if (!identical(gzip_status, 0L)) stop("gzip failed for ", files[[nm]])
    manifest[[length(manifest) + 1]] <- data.frame(scenario = scenario, phenotype_set = nm, phenotype_bed = gz, covariates = cov_file, stringsAsFactors = FALSE)
  }
  all_variance[[scenario]] <- res$variance
  all_loadings[[scenario]] <- res$loadings
  all_summary[[scenario]] <- res$summary
}
manifest <- do.call(rbind, manifest); manifest$job_id <- seq_len(nrow(manifest)); manifest <- manifest[, c("job_id", "scenario", "phenotype_set", "phenotype_bed", "covariates")]
write_tsv(manifest, file.path(base_dir, "data", "tensorqtl_job_manifest.tsv"))
write_tsv(do.call(rbind, all_variance), file.path(base_dir, "data", "gimap_pseudobulk_pc_variance_explained.tsv"))
write_tsv(do.call(rbind, all_loadings), file.path(base_dir, "data", "gimap_pseudobulk_pc_loadings.tsv"))
write_tsv(do.call(rbind, all_summary), file.path(base_dir, "data", "gimap_pseudobulk_scenario_summary.tsv"))

# Reference metadata.
write_tsv(data.frame(chr = 7, start = cis_start, end = cis_end, cluster_start = cluster_start_0based + 1L, cluster_end = cluster_end), file.path(base_dir, "data", "gimap_tensorqtl_region.tsv"))

# Loading similarity to cd8_nc_only PC axis, useful for interpretation.
loadings <- do.call(rbind, all_loadings)
ref <- subset(loadings, scenario == "cd8_nc_only" & !is.na(PC1))
sims <- list()
for (sc in unique(loadings$scenario)) {
  cur <- subset(loadings, scenario == sc)
  for (pc_ref in paste0("PC", 1:7)) for (pc in paste0("PC", 1:7)) {
    if (!pc_ref %in% names(ref) || !pc %in% names(cur)) next
    m <- merge(ref[, c("gene", pc_ref)], cur[, c("gene", pc)], by = "gene")
    names(m) <- c("gene", "ref", "cur")
    sims[[length(sims)+1]] <- data.frame(reference = "cd8_nc_only", reference_PC = pc_ref, scenario = sc, PC = pc, loading_correlation = suppressWarnings(cor(m$ref, m$cur, use = "complete.obs")))
  }
}
write_tsv(do.call(rbind, sims), file.path(base_dir, "data", "gimap_pseudobulk_loading_similarity.tsv"))
message("Prepared TensorQTL pseudobulk inputs: ", base_dir)
