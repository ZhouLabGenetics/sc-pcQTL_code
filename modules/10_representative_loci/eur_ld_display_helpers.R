display_ld_key <- function(pos, ea, oa) {
  ea <- toupper(as.character(ea))
  oa <- toupper(as.character(oa))
  allele_pair <- ifelse(ea <= oa, paste0(ea, "/", oa), paste0(oa, "/", ea))
  paste(as.integer(pos), allele_pair, sep = ":")
}

read_1000g_eur_ref_window <- function(chr, start_hg19, end_hg19, formal_root) {
  eur_plink_dir <- Sys.getenv("EUR_PLINK_DIR", unset = "")
  if (!nzchar(eur_plink_dir)) {
    stop("Set EUR_PLINK_DIR to the 1000G Phase 3 EUR PLINK reference directory")
  }
  bfile <- file.path(eur_plink_dir, sprintf("1000G.EUR.QC.%s", sub("^chr", "", chr)))
  bim_file <- paste0(bfile, ".bim")
  if (!file.exists(bim_file)) stop("Missing 1000G EUR PLINK BIM: ", bim_file)
  bim <- fread(
    bim_file,
    header = FALSE,
    col.names = c("chr", "snp_id", "cm", "pos_hg19", "a1", "a2"),
    showProgress = FALSE
  )
  bim[, `:=`(
    chr = sub("^chr", "", as.character(chr)),
    pos_hg19 = as.integer(pos_hg19),
    a1 = toupper(as.character(a1)),
    a2 = toupper(as.character(a2))
  )]
  pad <- 20000L
  bim <- bim[pos_hg19 >= (as.integer(start_hg19) - pad) & pos_hg19 <= (as.integer(end_hg19) + pad)]
  if (!nrow(bim)) return(list(bfile = bfile, ref = data.table()))
  bim[, ref_row := .I]
  ref_lift <- liftover_qtl_positions(
    bim[, .(ref_row, chr, pos_hg19)],
    row_col = "ref_row"
  )
  bim <- merge(bim, ref_lift[, .(ref_row, pos_hg38)], by = "ref_row", all.x = FALSE, sort = FALSE)
  bim[, `:=`(
    pos_hg38 = as.integer(pos_hg38),
    eur_ld_key_hg38 = display_ld_key(pos_hg38, a1, a2)
  )]
  bim <- unique(bim[is.finite(pos_hg38)], by = "eur_ld_key_hg38")
  list(bfile = bfile, ref = bim)
}

run_plink_anchor_ld <- function(bfile, chr, start_hg19, end_hg19, anchor_snp, out_dir) {
  plink_bin <- Sys.getenv("PLINK", Sys.getenv("PLINK_BIN", "plink"))
  if (!file.exists(plink_bin)) plink_bin <- "plink"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_prefix <- file.path(out_dir, paste0("eur_ld_", chr, "_", gsub("[^A-Za-z0-9_.-]", "_", anchor_snp)))
  args <- c(
    "--bfile", bfile,
    "--chr", sub("^chr", "", chr),
    "--from-bp", as.character(as.integer(start_hg19)),
    "--to-bp", as.character(as.integer(end_hg19)),
    "--r2",
    "--ld-snp", anchor_snp,
    "--ld-window", "999999",
    "--ld-window-kb", "2000",
    "--ld-window-r2", "0",
    "--allow-no-sex",
    "--out", out_prefix
  )
  status <- system2(plink_bin, args, stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(status, "status")) && attr(status, "status") != 0) {
    warning("PLINK LD extraction failed for ", anchor_snp, ":\n", paste(status, collapse = "\n"))
    return(data.table(eur_ref_snp = anchor_snp, lead_r2 = 1))
  }
  ld_file <- paste0(out_prefix, ".ld")
  if (!file.exists(ld_file) || file.info(ld_file)$size == 0) {
    warning("PLINK did not produce LD file for ", anchor_snp)
    return(data.table(eur_ref_snp = anchor_snp, lead_r2 = 1))
  }
  ld <- fread(ld_file, showProgress = FALSE)
  if (!nrow(ld)) return(data.table(eur_ref_snp = anchor_snp, lead_r2 = 1))
  snp_cols <- intersect(c("SNP_A", "SNP_B"), names(ld))
  r2_col <- intersect(c("R2", "R"), names(ld))[1]
  if (!length(snp_cols) || is.na(r2_col)) return(data.table(eur_ref_snp = anchor_snp, lead_r2 = 1))
  out <- rbindlist(lapply(snp_cols, function(col) {
    data.table(eur_ref_snp = as.character(ld[[col]]), lead_r2 = as.numeric(ld[[r2_col]]))
  }))
  out <- out[is.finite(lead_r2)]
  out[, lead_r2 := pmax(pmin(lead_r2, 1), 0)]
  out <- unique(out[order(-lead_r2)], by = "eur_ref_snp")
  rbind(
    data.table(eur_ref_snp = anchor_snp, lead_r2 = 1),
    out[eur_ref_snp != anchor_snp],
    fill = TRUE
  )
}

add_1000g_eur_display_ld <- function(qstats, gstats, q_lead, g_lead, chr, formal_root, out_dir) {
  for (dt in list(qstats, gstats)) {
    drop_cols <- intersect(
      c("lead_r2", "ld_source_for_plot", "eur_ld_key_hg38", "eur_ref_snp", "eur_ref_pos_hg19", "eur_ref_pos_hg38"),
      names(dt)
    )
    if (length(drop_cols)) dt[, (drop_cols) := NULL]
  }
  start_hg19 <- min(qstats$pos_hg19, na.rm = TRUE)
  end_hg19 <- max(qstats$pos_hg19, na.rm = TRUE)
  ref_obj <- read_1000g_eur_ref_window(chr, start_hg19, end_hg19, formal_root)
  ref <- ref_obj$ref
  if (!nrow(ref)) {
    qstats[, `:=`(lead_r2 = NA_real_, ld_source_for_plot = "not in 1000G EUR Phase3 display reference")]
    gstats[, `:=`(lead_r2 = NA_real_, ld_source_for_plot = "not in 1000G EUR Phase3 display reference")]
    return(list(qstats = qstats, gstats = gstats, anchor_snp = NA_character_, ref_n = 0L))
  }
  ref_map <- ref[, .(
    eur_ld_key_hg38,
    eur_ref_snp = snp_id,
    eur_ref_pos_hg19 = pos_hg19,
    eur_ref_pos_hg38 = pos_hg38
  )]
  qstats[, eur_ld_key_hg38 := display_ld_key(pos, effect_allele, other_allele)]
  gstats[, eur_ld_key_hg38 := display_ld_key(pos, effect_allele, other_allele)]
  qstats <- merge(qstats, ref_map, by = "eur_ld_key_hg38", all.x = TRUE, sort = FALSE)
  gstats <- merge(gstats, ref_map, by = "eur_ld_key_hg38", all.x = TRUE, sort = FALSE)

  anchor_snp <- qstats[q_idx == q_lead$q_idx & !is.na(eur_ref_snp)][1, eur_ref_snp]
  if (is.na(anchor_snp) || !nzchar(anchor_snp)) {
    anchor_snp <- qstats[!is.na(eur_ref_snp)][order(-abs(alpha), pvalue)][1, eur_ref_snp]
  }
  if (is.na(anchor_snp) || !nzchar(anchor_snp)) {
    anchor_snp <- gstats[!is.na(eur_ref_snp)][order(-abs(alpha), pvalue)][1, eur_ref_snp]
  }
  if (is.na(anchor_snp) || !nzchar(anchor_snp)) {
    qstats[, `:=`(lead_r2 = NA_real_, ld_source_for_plot = "not in 1000G EUR Phase3 display reference")]
    gstats[, `:=`(lead_r2 = NA_real_, ld_source_for_plot = "not in 1000G EUR Phase3 display reference")]
    return(list(qstats = qstats, gstats = gstats, anchor_snp = NA_character_, ref_n = nrow(ref)))
  }

  ld_map <- run_plink_anchor_ld(
    ref_obj$bfile,
    chr,
    start_hg19,
    end_hg19,
    anchor_snp,
    out_dir
  )
  qstats <- merge(qstats, ld_map, by = "eur_ref_snp", all.x = TRUE, sort = FALSE)
  gstats <- merge(gstats, ld_map, by = "eur_ref_snp", all.x = TRUE, sort = FALSE)
  qstats[, `:=`(
    lead_r2 = fifelse(is.finite(lead_r2), lead_r2, NA_real_),
    ld_source_for_plot = fifelse(
      !is.na(eur_ref_snp),
      "1000G EUR Phase3 PLINK display reference",
      "not in 1000G EUR Phase3 display reference"
    )
  )]
  gstats[, `:=`(
    lead_r2 = fifelse(is.finite(lead_r2), lead_r2, NA_real_),
    ld_source_for_plot = fifelse(
      !is.na(eur_ref_snp),
      "1000G EUR Phase3 PLINK display reference",
      "not in 1000G EUR Phase3 display reference"
    )
  )]
  qstats[, eur_ld_key_hg38 := NULL]
  gstats[, eur_ld_key_hg38 := NULL]
  list(qstats = qstats, gstats = gstats, anchor_snp = anchor_snp, ref_n = nrow(ref))
}
