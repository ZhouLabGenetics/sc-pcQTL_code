#!/usr/bin/env Rscript
# Source plots for the prespecified-trait pcQTL-vs-eQTL S-LDSC analysis:
# two panels showing marginal enrichment and marginal/joint tau*.
.local_libs <- c(
  Sys.getenv("SC_PCQTL_R_LIBS", unset = ""),
  Sys.getenv("SC_PCQTL_SHARED_R_LIBS", unset = "")
)
.local_libs <- unlist(strsplit(.local_libs[nzchar(.local_libs)], .Platform$path.sep, fixed = TRUE), use.names = FALSE)
if (length(.local_libs)) .libPaths(c(.local_libs, .libPaths()))

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
module_dir <- normalizePath(
  Sys.getenv("SC_PCQTL_SLDSC_ROOT", unset = file.path(dirname(args_file), "..")),
  mustWork = TRUE
)
R <- file.path(module_dir, "results/heritability_enrichment")
analysis_prefix <- Sys.getenv("SLDSC_ANALYSIS_PREFIX", unset = "prespecified247_h2qc")
out_dir <- Sys.getenv(
  "SC_PCQTL_SLDSC_FIGURE_OUTPUT_DIR",
  unset = file.path(module_dir, "figures/heritability_enrichment/manuscript")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rd <- function(p) read.delim(p, stringsAsFactors = FALSE, check.names = FALSE)
num <- function(x) suppressWarnings(as.numeric(x))

# Final manuscript house style: theme_classic, bold left-aligned titles, and
# grey55 reference lines. Unified
# main-text QTL palette (Figures 3/4): pcQTL = orange (#F58518), eQTL = blue
# (#4C78A8) -- consistent with the "pcQTL (orange) vs eQTL (blue)" captions.
theme_ms <- function(base = 8.5) theme_classic(base_size = base) +
  theme(axis.title = element_text(color="black", face="bold"), axis.text = element_text(color="black"),
        axis.line = element_line(linewidth=0.3), axis.ticks = element_line(linewidth=0.28),
        plot.title = element_text(face="bold", size=base, hjust=0, margin=margin(b=2)),
        legend.title = element_text(face="bold"), plot.margin = margin(2,4,2,2))
qtl_cols <- c("pcQTL"="#F58518", "eQTL"="#4C78A8")
save_panel <- function(name, plot, width, height) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot,
         width=width, height=height, device=cairo_pdf)
  ggsave(file.path(out_dir, paste0(name, ".png")), plot,
         width=width, height=height, dpi=400, bg="white")
}

## ---- Supplementary Figure S8: genome-wide enrichment and tau* ----------------
mpc <- rd(file.path(R, paste0("eur_sldsc_qtlsig_", analysis_prefix, "_marg_pc_summary"), paste0(analysis_prefix, "_marg_pc_annotation_meta.tsv")))
meq <- rd(file.path(R, paste0("eur_sldsc_qtlsig_", analysis_prefix, "_marg_eq_summary"), paste0(analysis_prefix, "_marg_eq_annotation_meta.tsv")))
jt  <- rd(file.path(R, paste0("eur_sldsc_qtlsig_", analysis_prefix, "_joint_summary"), paste0(analysis_prefix, "_joint_annotation_meta.tsv")))
pick <- function(df, qtl, model) data.frame(
  qtl=qtl, model=model, tau=num(df$tau_star_re), lo=num(df$tau_star_ci_lo), hi=num(df$tau_star_ci_hi),
  enr=num(df$enrichment_re))
a_df <- rbind(
  transform(pick(mpc, "pcQTL", "Marginal"), row="pcQTL (marginal)"),
  transform(pick(meq, "eQTL", "Marginal"), row="eQTL (marginal)"),
  transform(pick(jt[grepl("pcQTL", jt$annotation),], "pcQTL", "Conditional"), row="pcQTL (joint)"),
  transform(pick(jt[grepl("eQTL", jt$annotation),], "eQTL", "Conditional"), row="eQTL (joint)"))
a_df$row <- factor(a_df$row, levels=rev(c("pcQTL (marginal)","eQTL (marginal)","pcQTL (joint)","eQTL (joint)")))

# Left: heritability enrichment bar chart (marginal model, %h2 / %SNP).
enr_df <- rbind(
  data.frame(qtl="pcQTL", enr=num(mpc$enrichment_re), lo=num(mpc$enrichment_ci_lo), hi=num(mpc$enrichment_ci_hi)),
  data.frame(qtl="eQTL",  enr=num(meq$enrichment_re), lo=num(meq$enrichment_ci_lo), hi=num(meq$enrichment_ci_hi)))
enr_df$qtl <- factor(enr_df$qtl, levels=c("pcQTL","eQTL"))
pa_enr <- ggplot(enr_df, aes(qtl, enr, fill=qtl)) +
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.3, color="grey55") +
  geom_col(width=0.62, color="grey30", linewidth=0.12) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.16, linewidth=0.45) +
  geom_text(aes(y=hi, label=sprintf("%.2f%s", enr, "×")), vjust=-0.6, size=2.5, color="black") +
  scale_fill_manual(values=qtl_cols, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(title="Heritability enrichment", x=NULL, y=expression("Enrichment ("*h^2*" % / SNP %)")) +
  theme_ms() + theme(axis.text.x=element_text(size=8))
save_panel("Fig_pcQTL_eQTL_herit_a_enrichment", pa_enr, width=1.75, height=1.9)

# Right: genome-wide standardized effect tau* forest (marginal + conditional).
pa <- ggplot(a_df, aes(tau, row, color=qtl)) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.3, color="grey55") +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0.2, linewidth=0.5, orientation="y") + geom_point(size=1.9) +
  scale_color_manual(values=qtl_cols, name=NULL) +
  labs(title="Genome-wide per-SNP heritability",
       x=expression(tau*"* (standardized)"), y=NULL) +
  theme_ms() + theme(legend.position="none", axis.text.y=element_text(size=7.2))
save_panel("Fig_pcQTL_eQTL_herit_a_genomewide", pa, width=2.9, height=1.9)

cat("wrote 2 S-LDSC manuscript panels, with PDF and PNG outputs, to",
    out_dir, "\n")
