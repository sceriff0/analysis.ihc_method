#!/usr/bin/env Rscript
# =============================================================================
# figures/fig5.R  —  Figure 5. Proof-of-concept concordance of MIRAGE + FlowPath
#                    imaging immune quantification with two orthogonal
#                    transcriptomic proxies, six head-and-neck cases.
#
#   (a) Registered, phenotyped cases, cells coloured by population   [computed]
#   (b) mIF CD45+/all-cells, three hot vs three cold                 [computed]
#   (c) quanTIseq deconvolution vs imaging fraction, by population   [computed]
#
# As with fig4.R, nothing here is new science: (a)-(c) are paper_phenotype_map(),
# paper_immune_fraction_hotcold() and paper_deconv_scatter() out of
# code/paper_figures.R, the same objects analysis/paper_figures.Rmd prints. This
# script selects, strips, lays out and exports.
#
# ONE ARM. The cells come from one phenotyping arm (code/arms.R), the same one the
# website page names, because the arms' regions were annotated independently and
# cannot be pooled. Change ARM here and on the page together.
#
# ONE HOT/COLD AXIS. paper_hotcold_groups() ranks the Foy hotscore once and every
# panel reads that ranking: (b) groups on it, (c) colours by it, and (a) shows the
# most extreme case of each group that has imaging (paper_representative_cases()).
# Deriving (a)'s pair rather than pasting two patient IDs means the panels cannot
# drift: change the clinical table and all three move together. Override with
# FIG5A_CASES if a specific pair is wanted for image-quality reasons.
#
# NO TEST AND NO COEFFICIENT ON (b) OR (c). That is the legend's claim and it is
# enforced in the builders, not here — paper_immune_fraction_hotcold() deliberately
# never calls paired_spearman(), and paper_deconv_scatter() draws no fit line. Do
# not add either at assembly time: (c)'s axes use different denominators (imaging
# counts cells, deconvolution estimates a bulk-mixture fraction), so a regression
# line invites exactly the absolute-agreement reading the legend disclaims. Ranking
# is the claim.
#
# Run:  Rscript figures/fig5.R
# =============================================================================

source(file.path(tryCatch(here::here(), error = function(e) normalizePath(".")),
                 "figures", "_common.R"))
suppressPackageStartupMessages({
  library(dplyr); library(tibble)
})

root <- here_root
source(file.path(root, "code", "validation_helpers.R"))
source(file.path(root, "code", "paper_figures.R"))
source(file.path(root, "code", "arm_cells.R"))

# The phenotyping arm the manuscript panels are cut from. Same value as
# analysis/paper_figures.Rmd's ARM; the two must not diverge.
ARM            <- "massimo2"
# Set to a character vector of patient ids to pin panel (a); NULL derives them.
FIG5A_CASES    <- NULL
# "hot_score" ranks the continuous score and takes the top/bottom k. "immuno_phe"
# uses the clinical category instead. Same switch, same meaning, as the Rmd.
HOTCOLD_SOURCE <- "hot_score"

# --- Load: one arm, region cells + polygons, union metrics -------------------
# The same sequence as the page's `load-cells` chunk, so the objects are the
# page's objects. arm_annotations() errors without sf on purpose — a caller that
# asked for an outline should hear that it cannot be drawn. The map is still
# correct without one, so that degrades to NULL rather than aborting the figure.
as_spec <- arm_spec(ARM)
if (!dir.exists(as_spec$region_csv$path))
  stop("fig5: no cells for arm `", ARM, "` at ", as_spec$region_csv$path,
       ". Expected data/", as_spec$root, "/csv/<patient>/<patient>_<A|B|C>.csv")
as_cells  <- arm_cells(as_spec)
as_ucells <- arm_union_tier_cells(as_spec)
if (!nrow(as_cells)) stop("fig5: arm `", ARM, "` holds no cells.")

# An arm with no union tier hands back an EMPTY tibble with no columns.
pids     <- unique(c(as_cells$patient_id,
                     if ("patient_id" %in% names(as_ucells)) as_ucells$patient_id))
as_polys <- tryCatch(arm_annotations(as_spec, "region", patient_ids = pids),
  error = function(e) { warning("fig5: no annotation outlines (", conditionMessage(e),
                                ")", call. = FALSE); NULL })
as_upoly <- tryCatch(arm_annotations(as_spec, "union", patient_ids = pids),
                     error = function(e) NULL)
.prom    <- arm_promote_unregioned(as_spec, as_cells, as_polys, as_ucells, as_upoly)
as_cells <- .prom$cells; as_polys <- .prom$polys
as_union <- arm_metrics(as_spec, as_cells, as_polys, "union",
                        union_cells = as_ucells, union_polys = as_upoly)

# --- The hot/cold axis, shared by (a), (b) and (c) ---------------------------
groups <- paper_hotcold_groups(file.path(root, "data", "clinical_data.xlsx"),
                               patient_ids = as_union$patient_id,
                               source = HOTCOLD_SOURCE)
if (is.null(groups))
  stop("fig5: panel (b) needs data/clinical_data.xlsx (columns `ID PATIENT`, ",
       "`HOT score`) for the hot/cold axis.")

# --- (a) One representative case per group -----------------------------------
cases <- FIG5A_CASES %||%
  unname(paper_representative_cases(groups, available = as_cells$patient_id))
cases <- cases[!is.na(cases)]
if (!length(cases)) stop("fig5: no case is both grouped and present in the imaging.")

# THE MAPS ARE RASTERISED; EVERYTHING ELSE STAYS VECTOR. Two whole-slide maps are
# ~800k points, and as vector paths in a cairo PDF that is a 50 MB figure — five
# times the journal's cap, and one no reviewer's PDF viewer will scroll. Each map
# is drawn to a PNG at the panel's PLACED size and DPI, so nothing is lost that
# print could have shown, and placed as an image; the lineage key is drawn once as
# a real vector legend underneath, so its text is text.
MAP_W_MM <- MM[["two_col"]] / 2 - 4       # each map's placed width
MAP_H_MM <- 80
map_png <- function(p, id) {
  f <- file.path(here_root, "figures", "panels", paste0(id, ".png"))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  ggsave(f, p, width = MAP_W_MM, height = MAP_H_MM, units = "mm", dpi = DPI,
         device = ragg::agg_png, bg = "white")
  f
}
map_objs <- lapply(cases, function(pid) {
  # point_size is set for the PLACED size: ~0.1 mm is a cell's footprint at this
  # scale, so dense regions read as tissue and sparse ones as single cells. The
  # website builder's field-scaled default is tuned for a 9 in screen render.
  p <- paper_phenotype_map(as_cells, patient_id = pid, annots = as_polys,
                           title = NULL, point_size = 0.1)
  if (is.null(p)) return(NULL)
  for_panel(p, point_size = NULL) + theme(legend.position = "none")
})
names(map_objs) <- cases
map_objs <- map_objs[!vapply(map_objs, is.null, logical(1))]
if (!length(map_objs)) stop("fig5: paper_phenotype_map() returned nothing for ",
                            paste(cases, collapse = ", "))
maps <- Map(function(p, pid) image_panel(map_png(p, paste0("p5a_map_", pid))),
            map_objs, names(map_objs))

# The key: the union of the populations the maps actually show, drawn as an
# otherwise-empty ggplot whose only visible output is its legend. Same named
# scale as the maps, so the hues are the maps' hues by construction.
key_levels <- unique(unlist(lapply(map_objs, function(p) as.character(p$data$lineage))))
key_df     <- data.frame(x = 0, y = 0, lineage = lineage_legible(key_levels))
# Placed as a plot whose panel is void, so the legend is all that prints: no
# grob surgery, so it survives ggplot2's renaming of legend components.
p5a_key <- ggplot(key_df, aes(x, y, colour = lineage)) +
  geom_point(alpha = 0) +
  scale_colour_lineage(guide = guide_legend(override.aes = list(alpha = 1, size = 2.5),
                                            nrow = 2)) +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = BASE_PT),
        plot.margin = margin(0, 0, 0, 0))

# wrap_elements() so the two maps and their key read as ONE tagged panel. Without
# it, tag_levels tags each map separately and the real (b) becomes (c). The key sits
# UNDER the maps: it belongs to panel (a) alone — (b) and (c) do not use lineage
# colours, so the figure-wide legend strip must not carry a seven-population key
# for a figure two thirds of which never refers to it.
p5a <- wrap_elements(full = (Reduce(`|`, maps) / p5a_key) +
                       plot_layout(heights = c(1, 0.05)))

# --- (b) and (c) -------------------------------------------------------------
# Six points ARE the evidence on (b), and (c) has four per facet: drawn a little
# larger than the default print mark so they read as observations, not specks.
p5b <- for_panel(paper_immune_fraction_hotcold(as_union, groups), point_size = 1.4)
if (is.null(p5b)) stop("fig5: panel (b) empty after joining groups to metrics.")

paired_path <- file.path(root, "output", "paired_deconv.rds")
if (!file.exists(paired_path))
  stop("fig5: panel (c) needs output/paired_deconv.rds. Knit a molecular page ",
       "(analysis/molecular_", ARM, ".Rmd) first — its last chunk caches the paired ",
       "frame so this panel does not re-run immunedeconv.")
# `groups` is (b)'s axis; passing it is what makes hot in (b) mean hot in (c).
p5c <- for_panel(paper_deconv_scatter(readRDS(paired_path), method = "quantiseq",
                                      groups = groups), point_size = 1.3)
if (is.null(p5c)) stop("fig5: no quantiseq rows in ", paired_path)

# Placed sizes for the per-panel PDFs: 190 x 190 mm split 1 : 0.95 between the
# rows, bottom row split 1 : 2. (a)'s maps are also on disk as PNGs at the same
# placed size (figures/panels/p5a_map_<case>.png), for placing the image directly.
FIG5_W <- MM[["two_col"]]; FIG5_H <- 190
ROW_A  <- FIG5_H / 1.95; ROW_BC <- FIG5_H - ROW_A
save_panel(p5a, "p5a", FIG5_W, ROW_A)
save_panel(p5b, "p5b", FIG5_W / 3, ROW_BC)
save_panel(p5c, "p5c", FIG5_W * 2 / 3, ROW_BC)

# --- Assemble ----------------------------------------------------------------
# (b) is one axis with six points and (c) is four free-scaled facets, so the bottom
# row is split 1:2 rather than evenly — an equal split gives (b) whitespace it does
# not use and squeezes (c)'s facets below the width where the ranking is readable.
fig5 <- p5a / (p5b | p5c) +
  plot_layout(heights = c(1, 0.95), widths = c(1, 2), guides = "collect") +
  plot_annotation(tag_levels = TAG$tag_levels,
                  tag_prefix = TAG$tag_prefix, tag_suffix = TAG$tag_suffix) &
  theme(plot.tag = element_text(face = "bold"), legend.position = "bottom")

export_figure(fig5, "Fig5", width_mm = FIG5_W, height_mm = FIG5_H)

message("fig5: arm = ", ARM, " | panel (a) cases = ", paste(cases, collapse = ", "),
        " | groups = ", nrow(groups), " (", HOTCOLD_SOURCE, ")")
