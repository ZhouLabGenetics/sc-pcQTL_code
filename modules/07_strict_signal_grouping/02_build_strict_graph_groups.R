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
source(file.path(module06_code_dir, "scripts", "common.R"))
root_dir <- ROOT_DIR

in_qtl_gwas <- file.path(
  root_dir,
  "results/coloc/qtl_gwas_susie_official_finngen_all_finemapped/qtl_gwas_susie_official_finngen_coloc_summary.tsv"
)
in_qtl_qtl <- file.path(
  root_dir,
  "06_strict_signal_grouping/results/qtl_qtl/eqtl_pcqtl_signal_coloc_edges.tsv"
)
out_dir <- file.path(root_dir, "06_strict_signal_grouping/results/strict_graph")
dir_create(out_dir)

if (!file.exists(in_qtl_gwas)) stop("Missing QTL-GWAS coloc summary: ", in_qtl_gwas)
if (!file.exists(in_qtl_qtl)) stop("Missing QTL-QTL coloc edges: ", in_qtl_qtl)

thresholds <- c(0.70, 0.75, 0.80)

make_qtl_node <- function(celltype, cluster_id, qtl_type, phenotype_id, signal_index) {
  paste("QTL", celltype, cluster_id, qtl_type, phenotype_id, paste0("SIG", signal_index), sep = "||")
}

make_gwas_node <- function(celltype, cluster_id, phenocode, signal_index) {
  paste("GWAS", celltype, cluster_id, phenocode, paste0("SIG", signal_index), sep = "||")
}

extract_gwas_lead <- function(gwas_susie_rds, signal_index) {
  out <- data.table(
    gwas_susie_rds = gwas_susie_rds,
    gwas_signal_index = as.integer(signal_index),
    gwas_lead_variant = NA_character_,
    gwas_lead_chr = NA_character_,
    gwas_lead_pos = NA_integer_,
    gwas_lead_effect_allele = NA_character_,
    gwas_lead_other_allele = NA_character_,
    gwas_lead_marker_id = NA_character_,
    gwas_lead_pip = NA_real_,
    lead_status = "failed"
  )
  tryCatch({
    fm <- readRDS(gwas_susie_rds)
    fit <- fm$fit
    stats <- as.data.table(fm$stats)
    sig <- as.integer(signal_index)
    cs <- fit$sets$cs %||% list()
    cs_index <- as.integer(fit$sets$cs_index %||% seq_along(cs))
    cs_pos <- match(sig, cs_index)
    if (is.na(cs_pos) && sig >= 1L && sig <= length(cs)) cs_pos <- sig
    if (is.na(cs_pos) || !length(cs[[cs_pos]])) {
      out[, lead_status := "missing_signal"]
      return(out)
    }
    idx <- as.integer(cs[[cs_pos]])
    pip <- as.numeric(fit$pip[idx])
    idx <- idx[is.finite(pip)]
    pip <- pip[is.finite(pip)]
    if (!length(idx)) {
      out[, lead_status := "missing_pip"]
      return(out)
    }
    lead_idx <- idx[which.max(pip)]
    lead <- stats[lead_idx]
    out[, `:=`(
      gwas_lead_variant = standard_variant_id(lead$chr, lead$pos, lead$effect_allele, lead$other_allele),
      gwas_lead_chr = as.character(lead$chr),
      gwas_lead_pos = as.integer(lead$pos),
      gwas_lead_effect_allele = as.character(lead$effect_allele),
      gwas_lead_other_allele = as.character(lead$other_allele),
      gwas_lead_marker_id = as.character(lead$marker_id %||% NA_character_),
      gwas_lead_pip = max(pip, na.rm = TRUE),
      lead_status = "ok"
    )]
    out
  }, error = function(e) {
    out[, lead_status := paste0("failed: ", conditionMessage(e))]
    out
  })
}

union_find_components <- function(node_ids, edge_from, edge_to) {
  node_ids <- unique(as.character(node_ids))
  idx <- setNames(seq_along(node_ids), node_ids)
  parent <- seq_along(node_ids)
  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  union_nodes <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra != rb) parent[rb] <<- ra
    invisible(NULL)
  }
  a <- unname(idx[as.character(edge_from)])
  b <- unname(idx[as.character(edge_to)])
  keep <- is.finite(a) & is.finite(b)
  for (i in which(keep)) union_nodes(a[i], b[i])
  roots <- vapply(seq_along(node_ids), find_root, integer(1))
  data.table(node_id = node_ids, component_root = paste0("root_", roots))
}

build_graph_for_threshold <- function(qg, qq, threshold, gwas_leads) {
  qg_pass <- qg[status == "ok" & is.finite(pph4) & pph4 > threshold &
                  !is.na(qtl_signal_index) & !is.na(gwas_signal_index)]
  qq_pass <- qq[status == "ok" & is.finite(pph4) & pph4 > threshold &
                  !is.na(pc_signal_index) & !is.na(eqtl_signal_index)]

  if (!nrow(qg_pass)) {
    return(list(
      edges = data.table(),
      nodes = data.table(),
      groups = data.table(),
      gwas_members = data.table(),
      summary = data.table(
        threshold = threshold,
        n_qtl_gwas_edges = 0L,
        n_qtl_qtl_edges = nrow(qq_pass),
        n_gwas_components = 0L,
        n_eQTL_events = 0L,
        n_pcQTL_specific_events = 0L,
        n_shared_events = 0L,
        n_eQTL_only_events = 0L,
        n_unique_eQTL_hits = 0L,
        n_unique_pcQTL_specific_hits = 0L,
        percent_event_increase = NA_real_,
        percent_unique_increase = NA_real_
      )
    ))
  }

  qg_pass[, qtl_node := make_qtl_node(celltype, cluster_id, qtl_type, qtl_phenotype_id, qtl_signal_index)]
  qg_pass[, gwas_node := make_gwas_node(celltype, cluster_id, phenocode, gwas_signal_index)]
  qq_pass[, pc_node := make_qtl_node(celltype, cluster_id, "pcQTL", pc_id, pc_signal_index)]
  qq_pass[, eqtl_node := make_qtl_node(celltype, cluster_id, "eQTL", gene_id, eqtl_signal_index)]

  qg_edges <- qg_pass[, .(
    threshold,
    edge_type = "QTL_GWAS",
    from_node = qtl_node,
    to_node = gwas_node,
    pph4,
    celltype,
    cluster_id,
    phenocode,
    phenotype,
    qtl_type,
    qtl_phenotype_id,
    qtl_signal_index = as.integer(qtl_signal_index),
    gwas_signal_index = as.integer(gwas_signal_index),
    n_coloc_variants
  )]
  qq_edges <- qq_pass[, .(
    threshold,
    edge_type = "QTL_QTL",
    from_node = pc_node,
    to_node = eqtl_node,
    pph4,
    celltype,
    cluster_id,
    phenocode = NA_character_,
    phenotype = NA_character_,
    qtl_type = NA_character_,
    qtl_phenotype_id = NA_character_,
    qtl_signal_index = NA_integer_,
    gwas_signal_index = NA_integer_,
    n_coloc_variants
  )]
  edges <- rbindlist(list(qg_edges, qq_edges), use.names = TRUE, fill = TRUE)

  qtl_nodes_qg <- qg_pass[, .(
    threshold,
    node_id = qtl_node,
    node_type = qtl_type,
    celltype,
    cluster_id,
    phenotype_id = qtl_phenotype_id,
    signal_index = as.integer(qtl_signal_index),
    phenocode = NA_character_,
    phenotype = NA_character_,
    susie_rds = qtl_susie_rds,
    gwas_lead_variant = NA_character_,
    gwas_lead_chr = NA_character_,
    gwas_lead_pos = NA_integer_,
    gwas_lead_pip = NA_real_
  )]
  gwas_nodes <- qg_pass[, .(
    threshold,
    node_id = gwas_node,
    node_type = "GWAS",
    celltype,
    cluster_id,
    phenotype_id = phenocode,
    signal_index = as.integer(gwas_signal_index),
    phenocode,
    phenotype,
    susie_rds = gwas_susie_rds
  )]
  gwas_nodes <- merge(
    gwas_nodes,
    gwas_leads[, .(
      gwas_susie_rds,
      gwas_signal_index,
      gwas_lead_variant,
      gwas_lead_chr,
      gwas_lead_pos,
      gwas_lead_pip,
      lead_status
    )],
    by.x = c("susie_rds", "signal_index"),
    by.y = c("gwas_susie_rds", "gwas_signal_index"),
    all.x = TRUE
  )
  gwas_nodes[is.na(gwas_lead_variant), gwas_lead_variant := paste0(phenocode, "::unknown_signal_", signal_index)]

  qtl_nodes_qq <- rbindlist(list(
    qq_pass[, .(
      threshold,
      node_id = pc_node,
      node_type = "pcQTL",
      celltype,
      cluster_id,
      phenotype_id = pc_id,
      signal_index = as.integer(pc_signal_index),
      phenocode = NA_character_,
      phenotype = NA_character_,
      susie_rds = pc_susie_rds,
      gwas_lead_variant = NA_character_,
      gwas_lead_chr = NA_character_,
      gwas_lead_pos = NA_integer_,
      gwas_lead_pip = NA_real_
    )],
    qq_pass[, .(
      threshold,
      node_id = eqtl_node,
      node_type = "eQTL",
      celltype,
      cluster_id,
      phenotype_id = gene_id,
      signal_index = as.integer(eqtl_signal_index),
      phenocode = NA_character_,
      phenotype = NA_character_,
      susie_rds = eqtl_susie_rds,
      gwas_lead_variant = NA_character_,
      gwas_lead_chr = NA_character_,
      gwas_lead_pos = NA_integer_,
      gwas_lead_pip = NA_real_
    )]
  ), use.names = TRUE, fill = TRUE)
  nodes <- rbindlist(list(qtl_nodes_qg, gwas_nodes, qtl_nodes_qq), use.names = TRUE, fill = TRUE)
  nodes <- unique(nodes, by = c("threshold", "node_id"))

  comp <- union_find_components(nodes$node_id, edges$from_node, edges$to_node)
  nodes <- merge(nodes, comp, by = "node_id", all.x = TRUE)
  nodes[, signal_group_id := sprintf(
    "strict_sg_h4_%s_%06d",
    gsub("\\.", "p", sprintf("%.2f", threshold)),
    frank(component_root, ties.method = "dense")
  )]

  groups <- nodes[, .(
    n_nodes = .N,
    n_gwas = sum(node_type == "GWAS", na.rm = TRUE),
    n_eqtl = sum(node_type == "eQTL", na.rm = TRUE),
    n_pcqtl = sum(node_type == "pcQTL", na.rm = TRUE),
    celltypes = paste(sort(unique(celltype)), collapse = ","),
    clusters = paste(sort(unique(cluster_id)), collapse = ","),
    phenocodes = paste(sort(unique(phenocode[!is.na(phenocode)])), collapse = ","),
    gwas_lead_variants = paste(sort(unique(gwas_lead_variant[!is.na(gwas_lead_variant) & gwas_lead_variant != ""])), collapse = ","),
    member_qtl_phenotypes = paste(sort(unique(phenotype_id[node_type != "GWAS"])), collapse = ","),
    member_nodes = paste(sort(unique(node_id)), collapse = ";")
  ), by = .(threshold, signal_group_id)]
  groups[, coloc_class := fcase(
    n_gwas > 0L & n_eqtl > 0L & n_pcqtl > 0L, "shared",
    n_gwas > 0L & n_eqtl > 0L & n_pcqtl == 0L, "eQTL_only",
    n_gwas > 0L & n_eqtl == 0L & n_pcqtl > 0L, "pcQTL_specific",
    n_gwas > 0L, "GWAS_only",
    default = "no_GWAS"
  )]
  groups <- groups[n_gwas > 0L]

  gwas_members <- nodes[node_type == "GWAS", .(
    threshold,
    signal_group_id,
    celltype,
    cluster_id,
    phenocode,
    phenotype,
    gwas_signal_index = as.integer(signal_index),
    gwas_lead_variant,
    gwas_lead_chr,
    gwas_lead_pos,
    gwas_lead_pip
  )]
  gwas_members <- merge(
    gwas_members,
    groups[, .(threshold, signal_group_id, coloc_class)],
    by = c("threshold", "signal_group_id"),
    all.x = TRUE
  )
  gwas_members[, gwas_hit_id := paste(phenocode, gwas_lead_variant, sep = "::")]

  n_eqtl_events <- groups[coloc_class %in% c("eQTL_only", "shared"), .N]
  n_pc_specific_events <- groups[coloc_class == "pcQTL_specific", .N]
  n_unique_eqtl_hits <- gwas_members[coloc_class %in% c("eQTL_only", "shared"), uniqueN(gwas_hit_id)]
  n_unique_pc_specific_hits <- gwas_members[coloc_class == "pcQTL_specific", uniqueN(gwas_hit_id)]

  summary <- data.table(
    threshold = threshold,
    n_qtl_gwas_edges = nrow(qg_pass),
    n_qtl_qtl_edges = nrow(qq_pass),
    n_gwas_components = nrow(groups),
    n_eQTL_events = n_eqtl_events,
    n_pcQTL_specific_events = n_pc_specific_events,
    n_shared_events = groups[coloc_class == "shared", .N],
    n_eQTL_only_events = groups[coloc_class == "eQTL_only", .N],
    n_unique_eQTL_hits = n_unique_eqtl_hits,
    n_unique_pcQTL_specific_hits = n_unique_pc_specific_hits,
    percent_event_increase = if (n_eqtl_events > 0L) n_pc_specific_events / n_eqtl_events else NA_real_,
    percent_unique_increase = if (n_unique_eqtl_hits > 0L) n_unique_pc_specific_hits / n_unique_eqtl_hits else NA_real_
  )

  list(edges = edges, nodes = nodes, groups = groups, gwas_members = gwas_members, summary = summary)
}

message("Reading QTL-GWAS coloc summary")
qg <- fread(in_qtl_gwas)
qg <- qg[celltype %chin% names(CELLTYPE_EQTL_MAP)]
qg[, `:=`(
  qtl_signal_index = as.integer(qtl_signal_index),
  gwas_signal_index = as.integer(gwas_signal_index),
  pph4 = as.numeric(pph4)
)]

message("Reading QTL-QTL signal-level coloc edges")
qq <- fread(in_qtl_qtl)
qq <- qq[celltype %chin% names(CELLTYPE_EQTL_MAP)]
qq[, `:=`(
  pc_signal_index = as.integer(pc_signal_index),
  eqtl_signal_index = as.integer(eqtl_signal_index),
  pph4 = as.numeric(pph4)
)]

message("Extracting GWAS signal lead variants for de-duplication")
gwas_signal_keys <- unique(qg[status == "ok" & is.finite(pph4) & pph4 > min(thresholds) &
                                !is.na(gwas_signal_index),
                              .(gwas_susie_rds, gwas_signal_index = as.integer(gwas_signal_index))])
setorder(gwas_signal_keys, gwas_susie_rds, gwas_signal_index)
gwas_leads <- rbindlist(lapply(seq_len(nrow(gwas_signal_keys)), function(i) {
  if (i %% 500L == 0L) message(sprintf("Extracted GWAS leads: %d / %d", i, nrow(gwas_signal_keys)))
  extract_gwas_lead(gwas_signal_keys$gwas_susie_rds[i], gwas_signal_keys$gwas_signal_index[i])
}), fill = TRUE)
fwrite(gwas_leads, file.path(out_dir, "strict_gwas_signal_leads.tsv"), sep = "\t")

message("Building strict connected-component signal groups")
graphs <- lapply(thresholds, function(thr) build_graph_for_threshold(qg, qq, thr, gwas_leads))
names(graphs) <- sprintf("h4_%s", gsub("\\.", "p", sprintf("%.2f", thresholds)))

edges <- rbindlist(lapply(graphs, `[[`, "edges"), fill = TRUE)
nodes <- rbindlist(lapply(graphs, `[[`, "nodes"), fill = TRUE)
groups <- rbindlist(lapply(graphs, `[[`, "groups"), fill = TRUE)
gwas_members <- rbindlist(lapply(graphs, `[[`, "gwas_members"), fill = TRUE)
summary <- rbindlist(lapply(graphs, `[[`, "summary"), fill = TRUE)

fwrite(edges, file.path(out_dir, "strict_graph_edges.tsv"), sep = "\t")
fwrite(nodes, file.path(out_dir, "strict_graph_nodes.tsv"), sep = "\t")
fwrite(groups, file.path(out_dir, "strict_signal_groups.tsv"), sep = "\t")
fwrite(gwas_members, file.path(out_dir, "strict_signal_group_gwas_members.tsv"), sep = "\t")
fwrite(summary, file.path(out_dir, "strict_additional_hit_summary.tsv"), sep = "\t")

qc <- rbindlist(list(
  data.table(metric = "qtl_gwas_rows", value = nrow(qg)),
  data.table(metric = "qtl_qtl_rows", value = nrow(qq)),
  data.table(metric = "qtl_qtl_ok_rows", value = qq[status == "ok", .N]),
  data.table(metric = "qtl_qtl_failed_rows", value = qq[status == "failed", .N]),
  data.table(metric = "gwas_signal_lead_rows", value = nrow(gwas_leads)),
  data.table(metric = "gwas_signal_lead_ok_rows", value = gwas_leads[lead_status == "ok", .N]),
  data.table(metric = "strict_group_rows", value = nrow(groups))
), fill = TRUE)
fwrite(qc, file.path(out_dir, "strict_graph_qc.tsv"), sep = "\t")

print(summary)
