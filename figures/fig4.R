#!/usr/bin/env Rscript
# =============================================================================
# figures/fig4.R  —  Figure 4. Accuracy of MIRAGE's cross-panel DAPI-anchored
#                    registration.
#
#   (a) Two-colour DAPI overlay, before vs after registration  [author-supplied]
#   (b) VALIS-internal target registration error, by arm       [computed]
#   (c) DAPI-overlap Dice, by arm                              [computed]
#
# WHAT THIS SCRIPT DOES NOT DO. It does not compute anything. (b) and (c) come out
# of build_arm_figs() / build_reg_figs(), the SAME builders analysis/paper_figures.Rmd
# calls, so the panel in the manuscript and the panel on the website are one object
# with one owner (CLAUDE.md: "one figure, one owner"). This script only strips the
# web titles, resizes the ink for print, lays the three panels out and exports.
#
# THE TWO SOURCES ARE NOT INTERCHANGEABLE, and the script says which it used:
#   data/registration_arms/  the study slides registered once per configuration.
#                            This is what the manuscript figure must use.
#   data/benchmark/          mirage's synthetic sweep, whose moving panels are
#                            shifted copies of channel 0. A residual there measures
#                            a known injected offset, not real cross-acquisition
#                            difficulty — usable for cost, NOT for an accuracy claim.
# Same precedence as the Rmd. Falling back is loud, because a figure built from the
# synthetic sweep and captioned as Figure 4 would overstate the result.
#
# Run:  Rscript figures/fig4.R
# =============================================================================

source(file.path(tryCatch(here::here(), error = function(e) normalizePath(".")),
                 "figures", "_common.R"))
suppressPackageStartupMessages(library(tibble))

root <- here_root

# --- (a) DAPI overlay --------------------------------------------------------
# Author-supplied. panel_slot() in analysis/paper_figures.Rmd looks for the same id
# in the same directory, so supplying the file once lights up both the website panel
# and this one. Any of the extensions below is fine; first match wins.
FIG4A_ID   <- "fig4a_dapi_overlay"
FIG4A_DIR  <- file.path(root, "output", "figures", "manual")
fig4a_path <- Sys.glob(file.path(FIG4A_DIR, paste0(FIG4A_ID, ".*")))
fig4a_path <- fig4a_path[tolower(tools::file_ext(fig4a_path)) %in%
                           c("png", "pdf", "jpg", "jpeg", "tif", "tiff", "svg")]

p4a <- image_panel(
  if (length(fig4a_path)) fig4a_path[1] else "",
  placeholder = paste0("(a) TO SUPPLY — save the two-colour DAPI overlay as\n",
                       "output/figures/manual/", FIG4A_ID, ".png and re-run"))

if (!length(fig4a_path))
  message("fig4: panel (a) not supplied — placeholder drawn at the right geometry. ",
          "Save it to ", file.path("output", "figures", "manual",
                                   paste0(FIG4A_ID, ".png")), " and re-run.")

# --- (b) and (c) TRE and Dice by arm ----------------------------------------
arm_figs <- list(); arm_source <- NULL

if (dir.exists(file.path(root, "data", "registration_arms"))) {
  suppressMessages(source(file.path(root, "code", "registration_arms.R")))
  man <- arm_manifest()
  seg <- if (nrow(man)) read_arms_seg_qc(man) else tibble::tibble()
  if (nrow(seg)) {
    arm_figs   <- build_arm_figs(seg, tibble::tibble(), man)[
      c("01_final_residual_um_by_arm", "02_final_dice_by_arm")]
    arm_source <- "real"
  }
}
if (!length(arm_figs) &&
    file.exists(file.path(root, "data", "benchmark", "param_matrix.csv"))) {
  suppressMessages(source(file.path(root, "code", "registration_accuracy_plots.R")))
  arm_figs   <- build_reg_figs(file.path(root, "data", "benchmark"))[
    c("06_tre_by_arm", "07_dice_by_arm")]
  arm_source <- "synthetic"
}
arm_figs <- arm_figs[!vapply(arm_figs, is.null, logical(1))]

if (identical(arm_source, "synthetic"))
  warning("fig4: built from the SYNTHETIC sweep (data/benchmark/). Usable for cost ",
          "and scaling; NOT a real-tissue accuracy claim. Put the per-arm runs under ",
          "data/registration_arms/ before submission.", call. = FALSE)
if (!length(arm_figs))
  stop("fig4: no arm data. Symlink the per-arm runs under data/registration_arms/<arm>/ ",
       "(preferred) or drop mirage's paper_data CSVs into data/benchmark/.")

# --- Assemble ----------------------------------------------------------------
# (a) spans the width because a before/after overlay is read as one image; (b) and
# (c) sit side by side because they are the SAME arms measured two ways and the
# reader compares them across, not down. Both are coord_flip()ed with the arm names
# on y, so each keeps its own labels — the duplication is the price of letting either
# panel be read on its own, which is what a reviewer does.
p4b <- for_panel(arm_figs[[1]])
# (c) DROPS the arm names and lets (b) carry them for the pair. They are the same
# arms in the same order — the manuscript legend says so in as many words ("plotted
# for the same arms as (b)") — so repeating a column of long labels between two
# adjacent panels spends ~30mm of a 190mm figure restating the caption. It was also
# what pushed (c)'s x-axis title past the right edge and clipped the closing bracket
# off "(unitless, 0-1)". Do NOT do this if the two panels are ever separated.
# The left margin is not cosmetic slack: with the arm names gone, (c)'s panel starts
# flush at its own left edge and theme_paper()'s topleft tag lands ON the axis rule.
p4c <- for_panel(arm_figs[[2]]) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 10))

save_panel(p4a, "p4a"); save_panel(p4b, "p4b"); save_panel(p4c, "p4c")

fig4 <- p4a / (p4b | p4c) +
  plot_layout(heights = c(0.95, 1), guides = "collect") +
  plot_annotation(tag_levels = TAG$tag_levels,
                  tag_prefix = TAG$tag_prefix, tag_suffix = TAG$tag_suffix) &
  theme(plot.tag = element_text(face = "bold"), legend.position = "bottom")

# Two columns: (b) and (c) each carry a full arm axis and will not read at 140mm.
export_figure(fig4, "Fig4", width_mm = MM[["two_col"]], height_mm = 185)

message("fig4: arm source = ", arm_source %||% "none")
