#!/usr/bin/env Rscript
.local_libs <- c(Sys.getenv("SC_PCQTL_R_LIBS", unset = ""), Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = ""))
.libPaths(c(.local_libs[nzchar(.local_libs)], .libPaths()))
suppressPackageStartupMessages(library(data.table))
source("config/config.R")

raw_dir <- file.path(ROOT_DIR, "annotations", "raw")
out_dir <- file.path(ROOT_DIR, "annotations", "processed")
dir_create(out_dir)

parse_gtf_attrs <- function(x, key) {
  sub(sprintf('.*%s "([^"]+)".*', key), "\\1", x)
}
gtf_file <- file.path(raw_dir, "gencode.v19.annotation.gtf.gz")
if (file.exists(gtf_file)) {
  gtf <- fread(cmd = paste("zcat", shQuote(gtf_file), "| awk '$3==\"gene\"'"), sep = "\t", header = FALSE)
  genes <- gtf[, .(
    chr = std_chr(V1),
    start = as.integer(V4),
    end = as.integer(V5),
    strand = V7,
    gene_name = parse_gtf_attrs(V9, "gene_name"),
    gene_id = parse_gtf_attrs(V9, "gene_id")
  )]
  genes <- genes[chr %in% as.character(1:22)]
} else {
  genes <- fread(file.path(ROOT_DIR, "results", "cluster_sets", "expressed_gene_universe.tsv"))
  genes[, strand := NA_character_]
}
fwrite(unique(genes), file.path(out_dir, "gene_model.tsv"), sep = "\t")

go_file <- file.path(raw_dir, "goa_human.gaf.gz")
if (file.exists(go_file)) {
  go <- fread(cmd = paste("zcat", shQuote(go_file), "| grep -v '^!'"), header = FALSE, sep = "\t", fill = TRUE, quote = "")
  go_bp <- go[V9 == "P", .(gene_name = V3, go_id = V5)]
} else {
  go_bp <- data.table(gene_name = character(), go_id = character())
}
fwrite(unique(go_bp), file.path(out_dir, "go_bp_gene_terms.tsv"), sep = "\t")

hom_file <- file.path(raw_dir, "Homo_sapiens.GRCh37.75.homologies.tsv.gz")
biomart_paralog_file <- file.path(raw_dir, "biomart_human_paralogs_grch37.tsv")
if (file.exists(hom_file)) {
  hom <- tryCatch(fread(hom_file), error = function(e) data.table())
  sym_cols <- grep("gene.*name|external.*name|homology.*gene", names(hom), ignore.case = TRUE, value = TRUE)
  if (length(sym_cols) >= 2) {
    paralog <- unique(hom[, .(gene1 = get(sym_cols[1]), gene2 = get(sym_cols[2]))])
    paralog <- paralog[!is.na(gene1) & !is.na(gene2) & gene1 != gene2]
  } else {
    paralog <- data.table(gene1 = character(), gene2 = character())
  }
} else if (file.exists(biomart_paralog_file)) {
  hom <- fread(biomart_paralog_file)
  req <- c("Gene name", "Human paralogue associated gene name")
  if (all(req %in% names(hom))) {
    paralog <- unique(hom[nzchar(get(req[1])) & nzchar(get(req[2])), .(
      gene1 = get(req[1]),
      gene2 = get(req[2])
    )])
    paralog <- paralog[!is.na(gene1) & !is.na(gene2) & gene1 != gene2]
    if (file.exists(file.path(ROOT_DIR, "results", "cluster_sets", "expressed_gene_universe.tsv"))) {
      univ <- unique(fread(file.path(ROOT_DIR, "results", "cluster_sets", "expressed_gene_universe.tsv"), select = "gene_name")$gene_name)
      paralog <- paralog[gene1 %in% univ & gene2 %in% univ]
    }
    paralog[, `:=`(a = pmin(gene1, gene2), b = pmax(gene1, gene2))]
    paralog <- unique(paralog[, .(gene1 = a, gene2 = b)])
  } else {
    paralog <- data.table(gene1 = character(), gene2 = character())
  }
} else {
  paralog <- data.table(gene1 = character(), gene2 = character())
}
fwrite(paralog, file.path(out_dir, "paralog_pairs.tsv"), sep = "\t")

empty_table <- function(cols) {
  as.data.table(setNames(replicate(length(cols), character(), simplify = FALSE), cols))
}

prepare_abc <- function() {
  nasser_file <- file.path(raw_dir, "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz")
  abc_file <- file.path(raw_dir, "ENCFF811DOR_GM12878_ABC_thresholded_GRCh38.bed.gz")
  lift_bin <- file.path(raw_dir, "liftOver")
  chain <- file.path(raw_dir, "hg38ToHg19.over.chain.gz")
  out_file <- file.path(out_dir, "abc_links.tsv")
  cols <- c("chr", "start", "end", "gene_name", "score")
  if (file.exists(nasser_file)) {
    abc <- fread_maybe_gz(nasser_file, showProgress = TRUE)
    req <- c("chr", "start", "end", "name", "class", "TargetGene", "ABC.Score")
    if (!all(req %in% names(abc))) {
      stop("Nasser2021 ABC file missing required columns: ", paste(setdiff(req, names(abc)), collapse = ", "))
    }
    abc <- abc[
      chr %chin% paste0("chr", 1:22) &
        !is.na(TargetGene) & TargetGene != "" &
        class != "promoter" &
        as.numeric(`ABC.Score`) > 0.1,
      .(
        chr = std_chr(chr),
        start = as.integer(start),
        end = as.integer(end),
        gene_name = TargetGene,
        score = as.numeric(`ABC.Score`),
        enhancer_name = name
      )
    ]
    abc <- unique(abc[chr %in% as.character(1:22), .(chr, start, end, gene_name, score, enhancer_name)])
    fwrite(abc, out_file, sep = "\t")
    message("Prepared Nasser2021 ABC links: ", nrow(abc), " links")
    return(invisible(TRUE))
  }
  if (!file.exists(abc_file) || !file.exists(lift_bin) || !file.exists(chain)) {
    fwrite(empty_table(cols), out_file, sep = "\t")
    return(invisible(FALSE))
  }
  abc <- fread_maybe_gz(abc_file)
  if ("#chr" %in% names(abc)) setnames(abc, "#chr", "chr")
  req <- c("chr", "start", "end", "class", "TargetGene", "Score")
  if (!all(req %in% names(abc))) {
    fwrite(empty_table(cols), out_file, sep = "\t")
    return(invisible(FALSE))
  }
  abc <- abc[chr %chin% paste0("chr", 1:22) & !is.na(TargetGene) & TargetGene != "" & class != "promoter"]
  if (!nrow(abc)) {
    fwrite(empty_table(cols), out_file, sep = "\t")
    return(invisible(FALSE))
  }
  abc[, abc_id := sprintf("abc_%d", .I)]
  bed38 <- file.path(raw_dir, "abc_gm12878_grch38_for_liftover.bed")
  bed19 <- file.path(raw_dir, "abc_gm12878_hg19_lifted.bed")
  unmapped <- file.path(raw_dir, "abc_gm12878_hg19_unmapped.bed")
  fwrite(abc[, .(chr, start = as.integer(start), end = as.integer(end), abc_id)], bed38, sep = "\t", col.names = FALSE)
  status <- system2(lift_bin, c(bed38, chain, bed19, unmapped), stdout = TRUE, stderr = TRUE)
  if (!file.exists(bed19) || file.info(bed19)$size == 0) {
    warning("ABC liftOver produced no mapped intervals: ", paste(status, collapse = "\n"))
    fwrite(empty_table(cols), out_file, sep = "\t")
    return(invisible(FALSE))
  }
  lifted <- fread(bed19, header = FALSE, col.names = c("chr", "start", "end", "abc_id"))
  lifted[, chr := std_chr(chr)]
  abc_lite <- abc[, .(abc_id, gene_name = TargetGene, score = as.numeric(Score))]
  x <- merge(lifted, abc_lite, by = "abc_id", allow.cartesian = TRUE)
  x <- unique(x[chr %in% as.character(1:22), .(chr, start = as.integer(start), end = as.integer(end), gene_name, score)])
  fwrite(x, out_file, sep = "\t")
  message("Prepared ABC links: ", nrow(x), " lifted links")
  invisible(TRUE)
}

prepare_ctcf <- function() {
  gm12878_file <- file.path(raw_dir, "wgEncodeSydhTfbsGm12878Ctcfsc15914c20StdPk.narrowPeak.gz")
  f <- file.path(raw_dir, "wgEncodeRegTfbsClusteredV3.bed.gz")
  out_file <- file.path(out_dir, "ctcf_peaks.tsv")
  cols <- c("chr", "start", "end")
  if (file.exists(gm12878_file)) {
    x <- fread_maybe_gz(gm12878_file, header = FALSE)
    if (ncol(x) < 3) stop("GM12878 CTCF narrowPeak has fewer than 3 columns")
    x <- unique(x[, .(chr = std_chr(V1), start = as.integer(V2), end = as.integer(V3))])
    x <- x[chr %in% as.character(1:22)]
    fwrite(x, out_file, sep = "\t")
    message("Prepared GM12878 CTCF peaks: ", nrow(x))
    return(invisible(TRUE))
  }
  if (!file.exists(f)) {
    fwrite(empty_table(cols), out_file, sep = "\t")
    return(invisible(FALSE))
  }
  x <- fread_maybe_gz(f, header = FALSE, select = 1:4, col.names = c("chr", "start", "end", "name"))
  x <- unique(x[name == "CTCF", .(chr = std_chr(chr), start = as.integer(start), end = as.integer(end))])
  x <- x[chr %in% as.character(1:22)]
  fwrite(x, out_file, sep = "\t")
  message("Prepared CTCF peaks: ", nrow(x))
  invisible(TRUE)
}

prepare_tad <- function() {
  f <- file.path(raw_dir, "GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz")
  out_file <- file.path(out_dir, "tad_boundaries.tsv")
  cols <- c("chr", "start", "end")
  if (!file.exists(f)) {
    fwrite(empty_table(cols), out_file, sep = "\t")
    return(invisible(FALSE))
  }
  x <- fread_maybe_gz(f)
  req <- c("chr1", "x1", "x2")
  if (!all(req %in% names(x))) {
    fwrite(empty_table(cols), out_file, sep = "\t")
    return(invisible(FALSE))
  }
  b1 <- x[, .(chr = std_chr(chr1), start = as.integer(x1), end = as.integer(x1) + 1L)]
  b2 <- x[, .(chr = std_chr(chr1), start = as.integer(x2), end = as.integer(x2) + 1L)]
  b <- unique(rbind(b1, b2)[chr %in% as.character(1:22)])
  fwrite(b, out_file, sep = "\t")
  message("Prepared TAD boundaries: ", nrow(b))
  invisible(TRUE)
}

prepare_abc()
prepare_ctcf()
prepare_tad()

message("Prepared annotations in ", out_dir)
