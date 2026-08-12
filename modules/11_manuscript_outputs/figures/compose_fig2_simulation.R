#!/usr/bin/env Rscript
# Figure 2 (composite): pairwise hurdle co-expression calibration and power.
#   a  model-based null FPR              (fromstart simulation)
#   b  real-data permutation null FPR
#   c  model-based power, six settings   (fromstart simulation)
#   d  real-data permutation power, injected clusters
# Content is unchanged from the existing single panels; only the visual encoding,
# size/text ratio, and composition change. Run from the manuscript repository root.
#
# Colour groups the story into two families: our method and its two components all
# use ONE red family (so the components read as parts of the combined test), while
# the donor-aggregation baselines use a separate, more-visible BLUE family. Line
# style separates the combined test from its components (all same weight):
#   * Hurdle Union = our combined test -> SOLID, darkest red #A50F15.
#   * Hurdle Count = our count component -> DASHED, medium red #EF3B2C.
#   * Hurdle Zero  = our zero/detection component -> LONG-DASH, light red #FCAE91.
#   * Donor Pearson / Donor Spearman = donor-aggregation baselines -> SOLID, thin,
#     blue #3B6AA0 / #86ABD4 (a distinct family, visible but below the red method).
# So one colour family = "our method + its components", the other = "baselines";
# solid vs dashed/long-dash separates the combined test from its components.

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- normalizePath(dirname(args_file), mustWork = TRUE)
source(file.path(script_dir, "main_figure_common.R"))

out_dir <- file.path(OUTPUT_ROOT, "section2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- shared method styling -------------------------------------------------
LEG <- "Method / baseline"
# Legend / factor order tells the story: hero first, then our components, then baselines.
pal <- c(
  "Hurdle Union"   = "#A50F15",  # our method: darkest red (our family)
  "Hurdle Count"   = "#EF3B2C",  # count component: medium red (our family)
  "Hurdle Zero"    = "#FCAE91",  # zero component: light red (our family)
  "Donor Pearson"  = "#3B6AA0",  # baseline: blue family (distinct, more visible)
  "Donor Spearman" = "#86ABD4"   # baseline: light blue
)
lty <- c(
  "Hurdle Union"   = "solid",
  "Hurdle Count"   = "dashed",
  "Hurdle Zero"    = "longdash",
  "Donor Pearson"  = "solid",
  "Donor Spearman" = "solid"
)
lwd <- c(
  "Hurdle Union"   = 0.72,
  "Hurdle Count"   = 0.72,
  "Hurdle Zero"    = 0.72,
  "Donor Pearson"  = 0.62,
  "Donor Spearman" = 0.62
)
method_order <- names(pal)

threshold_levels <- c("0.05", "0.01", "0.005", "0.001")

label_threshold <- function(x) {
  out <- rep(NA_character_, length(x))
  out[abs(x - 0.05) < 1e-12]  <- "0.05"
  out[abs(x - 0.01) < 1e-12]  <- "0.01"
  out[abs(x - 0.005) < 1e-12] <- "0.005"
  out[abs(x - 0.001) < 1e-12] <- "0.001"
  out
}

# Manual scales for the three method aesthetics. Identical name + breaks merge them
# into a single legend showing each method's colour, line type, and line width.
method_scales <- function(name = LEG) list(
  scale_color_manual(values = pal, drop = FALSE, name = name),
  scale_linetype_manual(values = lty, drop = FALSE, name = name),
  scale_linewidth_manual(values = lwd, drop = FALSE, name = name)
)
# Redraw the hero (Hurdle Union) on top so it is never hidden under other curves.
union_on_top <- function(df) geom_line(
  data = df[method_component == "Hurdle Union"],
  aes(linetype = method_component, linewidth = method_component)
)
# Thin neutral diagonal = perfectly calibrated nominal FPR reference (panels a, b).
nominal_ref <- function(df) geom_line(
  data = unique(df[, .(threshold_label, p_threshold)]),
  aes(threshold_label, p_threshold, group = 1),
  color = "grey55", linetype = "dotted", linewidth = 0.45, inherit.aes = FALSE
)

## ===========================================================================
## fromstart (model-based) simulation: a (FPR) and c (power)
## ===========================================================================
fromstart_file <- Sys.getenv(
  "SC_PCQTL_MODEL_SIM_SUMMARY",
  unset = file.path(COQTL_WF, "02_simulations/00_from_start_simu/plots_target/alpha_curve_summary.tsv")
)
scenario_levels <- c("null", "zero_only_low", "count_only_low", "both_low",
                     "zero_only_high", "count_only_high", "both_high")

scenario_type_label <- function(s) {
  fcase(
    grepl("^zero_only",  s), "Zero-only signal",
    grepl("^count_only", s), "Count-only signal",
    grepl("^both",       s), "Zero+Count signal",
    default = NA_character_
  )
}
scenario_strength_label <- function(s) fifelse(grepl("_high$", s), "High", "Low")

curve_summary <- fread(fromstart_file)
curve_summary[, scenario := factor(as.character(scenario), levels = scenario_levels)]
curve_summary[, method_component := factor(as.character(method_component), levels = method_order)]
curve_summary[, threshold_label := factor(as.character(threshold_label), levels = threshold_levels)]

fpr_model <- curve_summary[scenario == "null"]
p_a <- ggplot(fpr_model,
              aes(threshold_label, mean_fpr, color = method_component, group = method_component)) +
  nominal_ref(fpr_model) +
  geom_line(aes(linetype = method_component, linewidth = method_component)) +
  union_on_top(fpr_model) +
  geom_point(size = 0.85) +
  scale_x_discrete(drop = FALSE) +
  method_scales() +
  labs(title = "Model-based null calibration",
       x = "Nominal p-value threshold", y = "Observed FPR") +
  theme_panel(9) + theme(legend.position = "none")

power_model <- curve_summary[scenario %in% c("zero_only_low", "count_only_low", "both_low",
                                             "zero_only_high", "count_only_high", "both_high")]
power_model[, signal_type := factor(scenario_type_label(as.character(scenario)),
                                    levels = c("Zero-only signal", "Count-only signal", "Zero+Count signal"))]
power_model[, signal_strength := factor(scenario_strength_label(as.character(scenario)),
                                        levels = c("Low", "High"))]
p_c <- ggplot(power_model,
              aes(threshold_label, mean_power, color = method_component, group = method_component)) +
  geom_line(aes(linetype = method_component, linewidth = method_component)) +
  union_on_top(power_model) +
  geom_point(size = 0.7) +
  scale_x_discrete(drop = FALSE) +
  method_scales() +
  facet_grid(signal_strength ~ signal_type, drop = FALSE) +
  labs(title = "Model-based power across six effect settings",
       x = "Nominal p-value threshold", y = "Observed power") +
  theme_panel(8.5) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey92", color = "grey70"),
        strip.text = element_text(face = "bold", size = 8))

## ===========================================================================
## real-data permutation-based simulation: b (FPR) and d (power)
## ===========================================================================
method_levels_raw <- c("SC_count", "SC_zero", "SC_combined", "PB")
method_labels <- c(SC_count = "Hurdle Count", SC_zero = "Hurdle Zero",
                   SC_combined = "Hurdle Union", PB = "Donor Spearman")
method_levels_sc <- c("Hurdle Union", "Hurdle Count", "Hurdle Zero", "Donor Spearman")

null_fpr_file <- Sys.getenv(
  "SC_PCQTL_REAL_DATA_PERMUTATION_NULL_FPR",
  unset = file.path(
    COQTL_WF,
    "02_simulations/01_main_text_results/01_null_calibration_final_sc/results/full_realcounts/fpr_comparison.tsv"
  )
)
fpr_dt <- fread(null_fpr_file)
fpr_dt <- fpr_dt[method %in% method_levels_raw]
fpr_dt[, method_component := factor(method_labels[method], levels = method_levels_sc)]
fpr_dt[, threshold_label := factor(label_threshold(p_threshold), levels = threshold_levels)]
fpr_dt <- fpr_dt[!is.na(threshold_label)][order(method_component, p_threshold)]
p_b <- ggplot(fpr_dt,
              aes(threshold_label, fpr, color = method_component, group = method_component)) +
  nominal_ref(fpr_dt) +
  geom_line(aes(linetype = method_component, linewidth = method_component)) +
  union_on_top(fpr_dt) +
  geom_point(size = 0.85) +
  scale_x_discrete(drop = FALSE) +
  method_scales() +
  labs(title = "Real-data permutation null calibration",
       x = "Nominal p-value threshold", y = "Observed FPR") +
  theme_panel(9) + theme(legend.position = "none")

power_compare_dir <- Sys.getenv(
  "SC_PCQTL_REAL_DATA_PERMUTATION_POWER_DIR",
  unset = file.path(
    COQTL_WF,
    "02_simulations/01_main_text_results/02_power_analysis_final_sc/results/full_realcounts/compare"
  )
)
selected_strengths <- c(0.02, 0.04, 0.06)
strength_dirs <- list.dirs(power_compare_dir, recursive = FALSE, full.names = TRUE)
power_dt <- rbindlist(lapply(strength_dirs, function(sd) {
  f <- file.path(sd, "power_comparison.tsv")
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  dt[, strength := as.numeric(sub("^strength_s", "", basename(sd))) / 100]
  dt
}), use.names = TRUE, fill = TRUE)
power_dt <- power_dt[method %in% method_levels_raw]
power_dt[, strength := round(strength, 2)]
power_dt <- power_dt[strength %in% selected_strengths]
power_dt[, method_component := factor(method_labels[method], levels = method_levels_sc)]
power_dt[, threshold_label := factor(label_threshold(p_threshold), levels = threshold_levels)]
power_dt <- power_dt[!is.na(threshold_label)]
power_dt[, signal_strength := factor(sprintf("Effect %.2f", strength),
                                     levels = sprintf("Effect %.2f", selected_strengths))]
power_dt <- power_dt[order(signal_strength, method_component, p_threshold)]
p_d <- ggplot(power_dt,
              aes(threshold_label, power, color = method_component, group = method_component)) +
  geom_line(aes(linetype = method_component, linewidth = method_component)) +
  union_on_top(power_dt) +
  geom_point(size = 0.7) +
  scale_x_discrete(drop = FALSE) +
  method_scales() +
  facet_wrap(~ signal_strength, nrow = 1) +
  labs(title = "Real-data permutation power across injected-cluster effect sizes",
       x = "Nominal p-value threshold", y = "Observed power") +
  theme_panel(8.5) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey92", color = "grey70"),
        strip.text = element_text(face = "bold", size = 8))

## ---- shared method legend (all five methods) -------------------------------
legend_src <- p_a +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(1.7, "lines"),
        legend.key.height = unit(0.55, "lines"),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.12, "lines"),
        legend.margin = margin(1, 2, 1, 2)) +
  guides(
    color    = guide_legend(nrow = 1, override.aes = list(size = 1.05)),
    linetype = guide_legend(nrow = 1),
    linewidth = guide_legend(nrow = 1)
  )
legend_grob <- extract_legend(legend_src)

## ---- assemble --------------------------------------------------------------
p_a <- add_tag(p_a, "a"); p_b <- add_tag(p_b, "b")
p_c <- add_tag(p_c, "c"); p_d <- add_tag(p_d, "d")

top_row <- arrangeGrob(p_a, p_b, ncol = 2)
body <- arrangeGrob(top_row, p_c, p_d, ncol = 1, heights = c(1.02, 1.45, 0.95))
figure <- arrangeGrob(body, legend_grob, ncol = 1, heights = c(1, 0.045))

save_composite(figure, file.path(out_dir, "main_simulation_composite"),
               width = 7.0, height = 8.1)
