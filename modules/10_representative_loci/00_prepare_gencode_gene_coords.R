#!/usr/bin/env Rscript

.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
})

root <- Sys.getenv("SC_PCQTL_FORMAL_COLOC_ROOT", unset = normalizePath(file.path(getwd(), "..", "06_finngen_susie_coloc"), mustWork = FALSE))
work_root <- Sys.getenv(
  "SC_PCQTL_LOCUSZOOM_WORK_ROOT",
  unset = file.path(root, "10_publication_locuszoom_redesign")
)
gtf <- file.path(root, "resources/gencode/gencode.v38.annotation.gtf.gz")
out_dir <- file.path(work_root, "data/gencode")
out_file <- file.path(out_dir, "gencode.v38.gene_coordinates.tsv")
exon_file <- file.path(out_dir, "gencode.v38.protein_coding_exons.tsv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(gtf)) stop("Missing GENCODE GTF: ", gtf)

extract_attr <- function(x, key) {
  pat <- paste0(key, " \"([^\"]+)\"")
  m <- regexec(pat, x, perl = TRUE)
  v <- regmatches(x, m)
  out <- rep(NA_character_, length(x))
  hit <- lengths(v) >= 2
  out[hit] <- v[hit][[1]][2]
  if (sum(hit) > 1) {
    out[hit] <- vapply(v[hit], function(z) z[2], character(1))
  }
  out
}

read_gtf_feature <- function(feature) {
  fread(
    cmd = sprintf("zcat %s | awk '$3==\"%s\"'", shQuote(gtf), feature),
    header = FALSE,
    col.names = c("chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
    showProgress = FALSE
  )
}

dt <- read_gtf_feature("gene")

dt[, `:=`(
  gene_id = sub("\\..*$", "", extract_attr(attribute, "gene_id")),
  gene_name = extract_attr(attribute, "gene_name"),
  gene_type = extract_attr(attribute, "gene_type"),
  chrom = sub("^chr", "", chrom),
  start = as.integer(start),
  end = as.integer(end)
)]

coords <- unique(
  dt[!is.na(gene_name), .(gene_id, gene_name, chrom, start, end, strand, gene_type)],
  by = c("gene_id", "gene_name", "chrom", "start", "end", "strand")
)
setorder(coords, chrom, start, end, gene_name)
fwrite(coords, out_file, sep = "\t", quote = FALSE)

exons <- read_gtf_feature("exon")
exons[, `:=`(
  gene_name = extract_attr(attribute, "gene_name"),
  gene_type = extract_attr(attribute, "gene_type"),
  chrom = sub("^chr", "", chrom),
  start = as.integer(start),
  end = as.integer(end)
)]
exons <- unique(
  exons[
    gene_type == "protein_coding" & !is.na(gene_name) & is.finite(start) & is.finite(end),
    .(gene_name, chrom, start, end)
  ]
)
# GENCODE duplicates PAR1 records on chrX/chrY at identical coordinates. Keep
# the first copy while retaining PAR2 loci whose chrX/chrY coordinates differ.
exons <- unique(exons, by = c("gene_name", "start", "end"))
setorder(exons, gene_name, chrom, start, end)
exons[, previous_end := shift(cummax(end)), by = .(gene_name, chrom)]
exons[, exon_group := cumsum(is.na(previous_end) | start > previous_end), by = .(gene_name, chrom)]
collapsed_exons <- exons[, .(
  exon_start = min(start),
  exon_end = max(end)
), by = .(gene_name, chrom, exon_group)]
collapsed_exons[, exon_group := NULL]
setcolorder(collapsed_exons, c("gene_name", "exon_start", "exon_end", "chrom"))
setorder(collapsed_exons, gene_name, exon_start, exon_end, chrom)
fwrite(collapsed_exons, exon_file, sep = "\t", quote = FALSE)

message("Wrote ", out_file, " with ", nrow(coords), " genes")
message("Wrote ", exon_file, " with ", nrow(collapsed_exons), " collapsed exons")
