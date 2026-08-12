joint_score_p <- function(stat_count, stat_detection, log.p = FALSE) {
  stat_joint <- stat_count + stat_detection
  stats::pchisq(stat_joint, df = 2, lower.tail = FALSE, log.p = log.p)
}

canonical_pair_key <- function(gene1, gene2) {
  left <- ifelse(gene1 <= gene2, gene1, gene2)
  right <- ifelse(gene1 <= gene2, gene2, gene1)
  paste(left, right, sep = "||")
}

canonical_gene_set <- function(genes) {
  values <- trimws(strsplit(genes, ",", fixed = TRUE)[[1]])
  paste(sort(unique(values[nzchar(values)])), collapse = ",")
}

jaccard_gene_sets <- function(left, right) {
  left <- unique(left)
  right <- unique(right)
  denom <- length(union(left, right))
  if (!denom) return(NA_real_)
  length(intersect(left, right)) / denom
}

best_jaccard_matches <- function(query_genes, query_chr, target_genes, target_chr) {
  if (!length(query_genes)) return(numeric())
  vapply(seq_along(query_genes), function(i) {
    candidates <- which(target_chr == query_chr[[i]])
    if (!length(candidates)) return(0)
    max(vapply(
      candidates,
      function(j) jaccard_gene_sets(query_genes[[i]], target_genes[[j]]),
      numeric(1)
    ))
  }, numeric(1))
}
