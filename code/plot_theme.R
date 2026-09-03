# =============================================================================
# plot_theme.R  —  THE house figure style for this project.
#
# Every figure in analysis/*.Rmd and code/*.R goes through this file, so a plot
# rendered by clinical_flowpath.Rmd and one rendered by benchmarks.Rmd are visually
# indistinguishable apart from their content. It is modelled on the journal
# figure the PI supplied as the target (Nature-style multi-panel):
#
#   * white panel, NO gridlines — the axes carry the reading, not a grid
#   * thin black axis lines with short OUTWARD ticks
#   * black text (not grey) in a small humanist sans (Helvetica/Arial)
#   * no strip background — facet labels are plain small bold text
#   * compact legends, colourbars horizontal and short
#   * colourblind-safe categorical palette (Okabe-Ito), blue-white-red diverging
#     and single-hue sequential continuous ramps
#
# HOW IT IS APPLIED. Sourcing this file has three side effects, on purpose:
#   1. theme_set(theme_paper())          -> every plot inherits the theme
#   2. paper_geom_defaults()             -> geom_point/line/text/boxplot defaults
#   3. options(ggplot2.discrete.*)       -> unspecified discrete scales use `oi`
# (2) and (3) are what make a NEW plot publication-ready without the author
# remembering anything. (1) only works if the plot does not append its own
# `theme_*()`: an inline `theme_classic()` REPLACES the active theme wholesale.
# So the rule for this repo is: never call `theme_classic()`/`theme_bw()`/
# `theme_minimal()` in an analysis. Add bare `theme(...)` for per-plot tweaks,
# which layers on top of the house theme instead of discarding it.
#
# SIZING FOR PRINT. Do not hand-pick fig.width. What holds type at one size across
# the whole figure set is a constant SHRINK FACTOR from rendered width to printed
# column, not a constant fig.width — so ask fig_width() for the number:
#
#   ```{r fig-name, fig.width = fig_width("double"), fig.height = fig_height("double")}
#
# "single" (85mm) = 4.25in, "oneandhalf" (114mm) = 5.7in, "double" (180mm) = 9in.
# All three shrink by 0.787 and land the 10pt base type at 7.9pt on paper, inside
# the 5-8pt journal range. Rendering wider than the target column is deliberate: it
# also gives the workflowr site a readable raster. Picking 11in "because the panel
# is busy" does not make the panel bigger — it makes that figure's type 6.5pt while
# its neighbour's is 9.6pt. Split a busy panel or move to a facet instead.
#
# Dependencies: ggplot2 + grid only (no tidyverse/here), so benchmark_plots.R,
# validation_helpers.R and any bare Rscript can all source it.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

# --- unit helpers ------------------------------------------------------------
# ggplot2 measures line widths and geom text in mm, journals specify them in pt.
# A `linewidth` of 1 draws at ~2.13pt; `size` in geom_text() is pt/2.845. These
# two converters let the specs below (and call sites) be written in real points.
pt_line <- function(pt) pt / 2.13          # pt -> ggplot2 `linewidth`
pt_text <- function(pt) pt / 2.845276      # pt -> geom_text/geom_label `size`

# --- font --------------------------------------------------------------------
# The style calls for a humanist sans (Helvetica/Arial), and "" is how you ask for
# it SAFELY. Do not name the font here.
#
# WHY. A family name has to be registered in the DEVICE's font database, not
# merely installed on the machine. grDevices::pdf() ships the Type 1 base-14 set
# (Helvetica, Times, Courier, ...) and knows nothing about an OS face like "Arial",
# "Nimbus Sans" or "DejaVu Sans" — those are what systemfonts::system_fonts()
# reports, which is a DIFFERENT database. Naming an OS-only face makes every knit
# with `dev = c("png", "pdf")` die at the first text grob with
#   Error in grid.Call.graphics(C_text, ...) : invalid font type
# and the PNG (cairo) pass gives no warning of it, because cairo resolves OS fonts
# happily. That failure mode cost a build; hence this comment.
#
# "" means "the device's own default family" — which IS Helvetica on pdf() and a
# Helvetica-metric sans on cairo/ragg. Identical look, no font database to satisfy,
# works unchanged on a headless cluster node.
#
# To force a specific face (a journal insisting on Arial, say): register it with
# the device first (see ?pdfFonts / ?grDevices::Type1Font), then set
#   options(ihc.plot.family = "Arial")
# BEFORE sourcing this file. paper_family() validates the request against the
# device databases and degrades to "" with a warning rather than letting a knit
# fail three chunks in.
paper_family <- function(family = getOption("ihc.plot.family", "")) {
  if (!nzchar(family)) return("")
  known <- tryCatch(unique(c(names(grDevices::pdfFonts()),
                             names(grDevices::postscriptFonts()))),
                    error = function(e) character(0))
  if (family %in% known) return(family)
  warning("plot_theme: font family ", sQuote(family), " is not registered with the ",
          "pdf/postscript device, so it would abort any knit that renders PDFs. ",
          "Falling back to the device default. Register it with grDevices::pdfFonts() ",
          "first if you need it.", call. = FALSE)
  ""
}

# --- the theme ---------------------------------------------------------------
# Built on theme_classic() because that is the only built-in with the right
# grammar (axis lines, no grid, no panel border); everything below re-specifies
# the parts theme_classic gets wrong for print — grey-ish text, chunky lines,
# inward ticks, a grey strip background, and oversized legend keys.
#
# Args let a caller deviate deliberately without hand-rolling a theme:
#   base_size   type size in pt (10 = the house default; see SIZING FOR PRINT)
#   grid        "none" (default), "y", "x", or "both" — a faint reference grid,
#               for the rare panel (e.g. a wide dot plot) that is unreadable
#               without one
#   axis_lines  FALSE drops the L-shaped axis, for heatmaps/tile plots
theme_paper <- function(base_size = 10, base_family = paper_family(),
                        grid = c("none", "y", "x", "both"), axis_lines = TRUE) {
  grid <- match.arg(grid)
  line_col <- "black"
  grid_line <- element_line(linewidth = pt_line(0.25), colour = "grey92")

  th <- theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # -- text: black, tight, left-aligned titles over the WHOLE plot so a
      # title and its y-axis label do not fight for the same column.
      text             = element_text(colour = "black"),
      plot.title       = element_text(face = "bold", size = rel(1.0),
                                      margin = margin(b = 4)),
      plot.subtitle    = element_text(colour = "grey30", size = rel(0.85),
                                      margin = margin(b = 6)),
      plot.caption     = element_text(colour = "grey45", size = rel(0.65),
                                      hjust = 1, margin = margin(t = 6)),
      plot.title.position = "plot", plot.caption.position = "plot",
      # Panel letters (a, b, c...) for patchwork::plot_annotation(tag_levels = "a").
      # NOTE: a tag and a plot.title share the top-left slot and WILL overlap. That
      # is the correct trade-off, because a panel inside an assembled figure should
      # not carry its own title anyway — the caption does that work. When composing
      # the final figure, set labs(title = NULL, subtitle = NULL) on each panel and
      # let patchwork place the tags.
      plot.tag          = element_text(face = "bold", size = rel(1.3), hjust = 0),
      plot.tag.position = "topleft",

      # -- axes: hairline black rules, short OUTWARD ticks, black labels
      axis.line   = if (axis_lines) element_line(linewidth = pt_line(0.75),
                                                 colour = line_col,
                                                 lineend = "square")
                    else element_blank(),
      axis.ticks  = if (axis_lines) element_line(linewidth = pt_line(0.75),
                                                 colour = line_col)
                    else element_blank(),
      axis.ticks.length = unit(2, "pt"),
      axis.text   = element_text(colour = "black", size = rel(0.85)),
      axis.title  = element_text(colour = "black", size = rel(0.95)),

      # -- panel/grid
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = if (grid %in% c("y", "both")) grid_line else element_blank(),
      panel.grid.major.x = if (grid %in% c("x", "both")) grid_line else element_blank(),
      # 9pt, not less: with no panel border the only thing keeping one panel's
      # last x tick label off its neighbour's first one is this gap.
      panel.spacing    = unit(9, "pt"),

      # -- facets: no grey box, just small bold text sitting on the panel
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = rel(0.85),
                                      colour = "black",
                                      margin = margin(3, 3, 3, 3)),

      # -- legends: compact, top-left, no box. Colourbars are short and thin;
      # `legend.key.*` is also what sizes guide_colourbar() in ggplot2 >= 3.5.
      legend.position   = "top",
      legend.justification = "left",
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.title      = element_text(size = rel(0.85)),
      legend.text       = element_text(size = rel(0.8)),
      legend.key.height = unit(8, "pt"),
      legend.key.width  = unit(14, "pt"),
      legend.margin     = margin(0, 0, 2, 0),

      plot.margin = margin(6, 8, 4, 6)
    )
  th
}

# Heatmap / tile variant: cells ARE the panel, so the L-shaped axis and the ticks
# only add clutter. Keeps the type, legend and title styling identical.
theme_paper_tile <- function(base_size = 10, ...) {
  theme_paper(base_size = base_size, axis_lines = FALSE, ...) +
    theme(axis.ticks.length = unit(0, "pt"))
}

# Dense facet grid variant. theme_paper()'s 9pt gutter and borderless panels are
# tuned for a handful of facets; a grid like 5 methods x 4 lineages puts twenty
# small panels edge to edge, where one panel's point cloud reads as continuous with
# the next. Worse, such grids are almost always drawn with `scales = "free"`, so
# adjacent panels do NOT share an axis — running them together invites precisely the
# cross-panel comparison the free scales forbid. A hairline border and a wider gutter
# make each panel its own coordinate system, visibly. Type, colour and legend
# styling are untouched, so a gridded figure still matches every other one.
theme_paper_panels <- function(base_size = 10, spacing = 14, ...) {
  theme_paper(base_size = base_size, ...) +
    theme(panel.border  = element_rect(colour = "grey75", fill = NA, linewidth = 0.4),
          panel.spacing = unit(spacing, "pt"))
}

# --- palettes ----------------------------------------------------------------
# Okabe-Ito: the standard 8-colour colourblind-safe categorical palette. Ordered
# so the first two (blue, vermillion) are the maximally separable pair, which is
# what a 2-level scale gets.
oi <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
        "#E69F00", "#56B4E9", "#F0E442", "#000000")

# `oi` runs out at 8, and several figures here colour by `patient_id` with more
# patients than that. Extended set: the first 8 entries ARE `oi`, so a 3-patient
# and a 14-patient panel still open with the same blue; the tail adds darker,
# still-distinguishable hues from the Okabe-Ito barrier-free set. Past ~12
# categories colour stops separating anything — facet, or label points directly,
# rather than reaching for more colours.
oi_ext <- c(oi,
            "#004949", "#920000", "#490092", "#B66DFF",
            "#6DB6FF", "#924900", "#009292", "#FF6DB6")

# Diverging (ColorBrewer RdBu endpoints): for signed quantities around a real
# zero — log2 fold changes, correlations, differences. Blue = low/negative.
PAL_DIV <- c(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B")

# Sequential single-hue: for magnitudes with no meaningful midpoint (counts,
# densities, fractions). Perceptually ordered and safe in greyscale.
PAL_SEQ <- c(low = "#F7FBFF", high = "#08519C")

# Semantic colours for the two recurring ANNOTATION marks, so "the dashed x = y
# line" and "the fitted trend" look the same in every report. Reference lines are
# deliberately achromatic — a red identity line competes with the data for
# attention, and in these plots the data is what carries the finding.
REF_LINE <- "grey35"       # identity line, threshold, bias / limits of agreement
FIT_LINE <- "#0072B2"      # geom_smooth trend (= oi[1])

# Semantic colours for the immune "hot/cold" phenotype: HOT -> red, COLD -> light
# blue, intermediate/other -> orange/grey. Robust to case/spelling. Returns a
# named vector keyed by the given levels, for scale_colour_manual/scale_fill_manual.
hotcold_cols <- function(levels) {
  lv  <- as.character(levels)
  key <- toupper(trimws(lv))
  col <- ifelse(grepl("HOT|INFLAM", key),               "#D7191C",   # red
         ifelse(grepl("COLD|DESERT", key),              "#74ADD1",   # light blue
         ifelse(grepl("INTERMED|VARI|MIX|EXCLUD", key), "#FDAE61",   # orange
                                                        "grey65")))
  stats::setNames(col, lv)
}

# Order immune-phenotype levels cold -> intermediate -> hot (axis/legend order).
hotcold_order <- function(x) {
  lv   <- unique(as.character(x[!is.na(x)]))
  key  <- toupper(trimws(lv))
  rank <- ifelse(grepl("COLD|DESERT", key), 1L,
          ifelse(grepl("INTERMED|VARI|MIX|EXCLUD", key), 2L,
          ifelse(grepl("HOT|INFLAM", key), 3L, 4L)))
  factor(x, levels = lv[order(rank)])
}

# --- scale shorthands --------------------------------------------------------
# Thin wrappers so a call site names the SEMANTICS ("this is diverging") rather
# than repeating hex codes. Every `...` passes through to the underlying scale,
# so limits/name/labels/trans all still work.
scale_fill_div <- function(midpoint = 0, ...)
  scale_fill_gradient2(low = PAL_DIV[["low"]], mid = PAL_DIV[["mid"]],
                       high = PAL_DIV[["high"]], midpoint = midpoint,
                       na.value = "grey92", ...)

scale_colour_div <- function(midpoint = 0, ...)
  scale_colour_gradient2(low = PAL_DIV[["low"]], mid = PAL_DIV[["mid"]],
                         high = PAL_DIV[["high"]], midpoint = midpoint,
                         na.value = "grey92", ...)

scale_fill_seq <- function(...)
  scale_fill_gradient(low = PAL_SEQ[["low"]], high = PAL_SEQ[["high"]],
                      na.value = "grey92", ...)

scale_colour_seq <- function(...)
  scale_colour_gradient(low = PAL_SEQ[["low"]], high = PAL_SEQ[["high"]],
                        na.value = "grey92", ...)

scale_colour_oi <- function(...) scale_colour_manual(values = unname(oi), ...)
scale_fill_oi   <- function(...) scale_fill_manual(values = unname(oi), ...)

# ORDINAL discrete: levels that have an order (image size, tile count) or simply
# outnumber the 8 Okabe-Ito colours. `oi` is wrong for both — an unordered hue set
# hides the ordering, and scale_*_manual errors outright when it runs out of
# colours. One ramp (viridis D) for every such case, so ordered discrete looks the
# same everywhere; the repo previously mixed viridis options B, C and D.
scale_colour_ordinal <- function(...) scale_colour_viridis_d(option = "D", ...)
scale_fill_ordinal   <- function(...) scale_fill_viridis_d(option = "D", ...)
scale_color_ordinal  <- scale_colour_ordinal

# American/British spelling aliases, so a call site can use either.
scale_fill_diverging <- scale_fill_div
scale_color_div      <- scale_colour_div
scale_color_seq      <- scale_colour_seq
scale_color_oi       <- scale_colour_oi

# Short horizontal colourbar with the title above it, matching the reference
# figure's under-panel bars. ggplot2 3.5 moved bar sizing out of the guide and
# into the theme (and deprecated barwidth/barheight), so detect which API this
# installation has instead of pinning one.
# `title` defaults to waiver(), NOT NULL: in ggplot2 guides waiver() means
# "inherit the scale's name" while NULL means "draw no title at all", so a NULL
# default would silently strip the label off every colourbar it is used on.
guide_cbar <- function(title = waiver(), width = 72, height = 6, ...) {
  args <- list(title = title, title.position = "top", direction = "horizontal",
               ticks = FALSE, ...)
  if ("theme" %in% names(formals(guide_colourbar))) {
    args$theme <- theme(legend.key.width  = unit(width, "pt"),
                        legend.key.height = unit(height, "pt"))
    # 3.5+ renamed these to `theme` entries; drop the deprecated spellings.
    args$title.position <- NULL
    args$theme <- args$theme + theme(legend.title.position = "top")
  } else {
    args$barwidth  <- unit(width, "pt")
    args$barheight <- unit(height, "pt")
  }
  do.call(guide_colourbar, args)
}

# --- semantic categorical palettes -------------------------------------------
# These live HERE, not next to the figures that use them, and that placement is
# the whole point. A named palette defined inside code/paper_figures.R can only be
# obeyed by paper_figures.R; a second file drawing the same categories reaches for
# `oi` instead and assigns colours BY POSITION, so a lineage's hue silently depends
# on which levels happen to be present in that file's data frame. That is how CD8T
# came to be vermillion in Fig 3 and green in the composition panel. Anything that
# names a recurring category belongs in this file, and every call site uses the
# scale_*_() wrapper rather than passing the vector to scale_*_manual() by hand.
#
# The rule for adding one: it must be NAMED (level -> colour). Naming is what makes
# the colour independent of which levels a given panel happens to contain, and it
# also settles the `drop` question that an unnamed palette makes fraught. With an
# unnamed vector, dropping an unused level shifts every colour after it, so you are
# forced into drop = FALSE; with a named one, drop cannot move a single colour and
# only decides which legend KEYS are drawn. So these scales drop by default — an
# unused level otherwise leaves a labelled key with no swatch beside it, which reads
# as a rendering fault. Pass drop = FALSE deliberately when a figure SET needs one
# identical legend across panels (a whole-slide map and its inset, say).

# Cell lineages. A map with fourteen colours reads as noise at figure size, so the
# taxonomy is collapsed to the seven populations the paper argues about plus one
# catch-all; `lineage_legible()` does the collapsing.
#
# The five immune populations get saturated Okabe-Ito hues; the three structural
# classes are deliberately desaturated so they read as substrate rather than as
# findings. The hard part is that "desaturated" used to mean three greys, and at
# point size they blurred: Tumor vs other was dE 15.8 in CIE Lab and Tumor vs Stroma
# 17.3, both under the ~25 a small mark needs. They are now separated on BOTH axes a
# grey can vary in — lightness (L* 67 / 53 / 94) and hue (Tumor cool, Stroma warm) —
# which also keeps them apart under colour-vision deficiency, where hue collapses and
# only the lightness ladder survives. Worst pair in the whole palette is now dE 25.5.
# Re-check with convertColor(..., "Lab") before changing any of these three.
LEGIBLE_LINEAGES <- c("Tumor", "CD8T", "CD4T", "Treg", "NK", "Immune_other", "Stroma")

LINEAGE_COLS <- c(
  Tumor        = "#96A5B3",   # cool slate: the substrate the immune cells sit on
  CD8T         = "#D55E00",
  CD4T         = "#0072B2",
  Treg         = "#CC79A7",
  NK           = "#009E73",
  Immune_other = "#E69F00",
  Stroma       = "#8C7B6B",   # warm brown: same family as Tumor, opposite hue
  other        = "#EDEFF1"    # near-white: a catch-all should recede, not read
)

# What each lineage is CALLED in a legend. The keys are the analysis vocabulary and
# must stay as they are — cell_lineage() emits them and the tests assert them — but
# "Immune_other" is a column code, not a label a reader should meet in a figure.
LINEAGE_LABELS <- c(
  Tumor = "Tumour", CD8T = "CD8 T", CD4T = "CD4 T", Treg = "Treg", NK = "NK",
  Immune_other = "Immune (other)", Stroma = "Stroma", other = "Other",
  Unknown = "Unknown", Unclassified = "Unclassified")

# Marker-GATED populations (CD45+, CD3+ CD45+, GZMB+ NK) are not lineages — they
# are thresholded readouts that overlap the lineages — so they must not borrow a
# lineage hue and imply they are a disjoint population. They get achromatic greys
# on the same lightness ladder, which also groups them visually as "the gated set".
GATED_COLS <- c("CD45+" = "grey25", "CD3+ CD45+" = "grey50", "GZMB+ NK" = "grey72")

# Every category this project colours by lineage, in one lookup: the lineages,
# their gated companions, and the display spellings the two upstream tools use.
# `cell_lineage()` should normalise before we get here, but a scale that silently
# drops an unmapped level is worse than one that still draws it.
LINEAGE_PALETTE <- c(LINEAGE_COLS, GATED_COLS,
                     Unknown = "grey85", Unclassified = "grey85")

# Registration arms. Ordered so the valis micro-registration ladder (0 -> 1 -> 2)
# reads as a progression in the blue-green direction, with the two non-ladder arms
# in contrasting grey and vermillion.
ARM_KIND_COLS <- c("valis · micro 0" = "#0072B2", "valis · micro 1" = "#56B4E9",
                   "valis · micro 2" = "#009E73", "valis"           = "#7F8C8D",
                   "tiled (STARE)"   = "#D55E00")

# Collapse a raw lineage column to the legible subset, as an ORDERED factor with
# the full level set present (drop = FALSE then keeps colours stable across
# panels that happen to be missing a population).
lineage_legible <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !x %in% LEGIBLE_LINEAGES] <- "other"
  factor(x, levels = c(LEGIBLE_LINEAGES, "other"))
}

# `labels` defaults to the display names but stays overridable: a marker-gated panel
# passes its own, and the gated labels ("CD45+") are already reader-facing.
.lineage_labels <- function(breaks) {
  b <- as.character(breaks)
  ifelse(b %in% names(LINEAGE_LABELS), LINEAGE_LABELS[b], b)
}
scale_colour_lineage <- function(..., drop = TRUE, name = NULL, labels = .lineage_labels)
  scale_colour_manual(values = LINEAGE_PALETTE, drop = drop, name = name,
                      labels = labels, na.value = "grey85", ...)
scale_fill_lineage <- function(..., drop = TRUE, name = NULL, labels = .lineage_labels)
  scale_fill_manual(values = LINEAGE_PALETTE, drop = drop, name = name,
                    labels = labels, na.value = "grey85", ...)
scale_color_lineage <- scale_colour_lineage

# Cell COMPARTMENTS: the coarsest reading of a phenotype map, for the paired
# whole-slide / region views of Fig 5(a). On a whole slide the eye can hold one
# distinction — tumour against everything else — and on a region inset three:
# tumour, immune (the five immune lineages pooled) and stroma (the SMA+ call). One
# named palette serves both, so the overview and its inset agree on every hue and
# the grey that means "not tumour" on the slide is the grey that means "other" in
# the region. Red and green are the pair the manuscript asked for; they sit at
# L* 40 and 57, so the lightness ladder still separates them where hue collapses.
COMPARTMENTS <- c("Tumour", "Immune", "Stroma", "Other")

COMPARTMENT_COLS <- c(
  Tumour       = "#B2182B",   # red: PAL_DIV's high end, already the house red
  Immune       = "#009E73",   # Okabe-Ito bluish green
  Stroma       = "#B66DFF",   # violet, from the barrier-free extension in oi_ext
  Other        = "grey78",
  `Non-tumour` = "grey78"     # the overview's one non-tumour class, same grey
)

# Collapse a lineage vector (cell_lineage() output) to compartments. `binary = TRUE`
# is the overview's reading: Tumour against Non-tumour. Anything cell_lineage() did
# not resolve lands in Other / Non-tumour, never in a coloured class.
lineage_compartment <- function(x, binary = FALSE) {
  x <- as.character(x)
  comp <- ifelse(x %in% "Tumor", "Tumour",
          ifelse(x %in% c("CD8T", "CD4T", "Treg", "NK", "Immune_other"), "Immune",
          ifelse(x %in% "Stroma", "Stroma", "Other")))
  if (isTRUE(binary))
    return(factor(ifelse(comp == "Tumour", "Tumour", "Non-tumour"),
                  levels = c("Tumour", "Non-tumour")))
  factor(comp, levels = COMPARTMENTS)
}

scale_colour_compartment <- function(..., drop = TRUE, name = NULL)
  scale_colour_manual(values = COMPARTMENT_COLS, drop = drop, name = name,
                      na.value = "grey85", ...)
scale_fill_compartment <- function(..., drop = TRUE, name = NULL)
  scale_fill_manual(values = COMPARTMENT_COLS, drop = drop, name = name,
                    na.value = "grey85", ...)
scale_color_compartment <- scale_colour_compartment

# --- Measurement SCOPE -------------------------------------------------------
# The three nested answers to "which cells count?", which appear together on every
# arm's scope-comparison figure and separately elsewhere. A NAMED palette because
# the same three levels recur across panels that do not all contain all three: an
# unnamed scale would assign by factor POSITION, so `union` would be blue in a panel
# that happens to omit `whole_slide` and orange in one that does not.
#
#   whole_slide     every cell in the export, no polygon consulted at all
#   annotation_all  inside the dissolved union polygon (massimo1's annotation_all)
#   annotation_k    one pathologist region, ANNOTATION_1..k
SCOPE_COLS <- c(
  whole_slide    = "#000000",   # Okabe-Ito black — the widest set, the reference
  annotation_all = "#0072B2",   # blue  — one value per patient, like whole_slide
  annotation_k   = "#E69F00"    # orange — many values per patient
)

SCOPE_LABELS <- c(
  whole_slide    = "whole slide (no annotation)",
  annotation_all = "annotation_all (union polygon)",
  annotation_k   = "single annotation"
)

# Order widest -> narrowest, which is also nesting order, so a legend reads as a
# progression rather than alphabetically.
scope_factor <- function(x) factor(as.character(x), levels = names(SCOPE_COLS))

.scope_labels <- function(breaks) {
  b <- as.character(breaks)
  ifelse(b %in% names(SCOPE_LABELS), SCOPE_LABELS[b], b)
}
scale_colour_scope <- function(..., drop = TRUE, name = NULL, labels = .scope_labels)
  scale_colour_manual(values = SCOPE_COLS, drop = drop, name = name,
                      labels = labels, na.value = "grey85", ...)
scale_fill_scope <- function(..., drop = TRUE, name = NULL, labels = .scope_labels)
  scale_fill_manual(values = SCOPE_COLS, drop = drop, name = name,
                    labels = labels, na.value = "grey85", ...)
scale_color_scope <- scale_colour_scope

# --- Bulk-RNA DENOMINATORS ---------------------------------------------------
# The bulk-RNA comparison asks "does a marker's expression track the fraction of
# cells the gate calls positive", and `fraction of WHICH cells` has four answers
# because two independent restrictions cross:
#
#                      every cell        tumour-lineage cells only
#   whole slide        all_wholeslide    tumor_wholeslide
#   inside annotation  all_annotation    tumor_annotation
#
# The bulk RNA was extracted from TUMOUR material, so `all_wholeslide` — the
# comparison this page has always drawn — puts a whole-slide IHC fraction against
# an expression value from a tumour-restricted sample. The other three close that
# gap from the two directions available: spatially (restrict to the pathologist's
# polygon) and by cell identity (restrict to cells phenotyped Tumor). The bottom-
# right cell does both and is the closest match to what was sequenced.
#
# NAMED, not positional, for the usual reason: not every panel carries all four
# levels — a marker with no tumour-cell positives drops one — and an unnamed scale
# would then recolour the remaining three by position.
#
# HUE ENCODES THE CELL SET, LIGHTNESS THE SCOPE. Cool = every cell, warm = tumour
# cells; the lighter member of each pair is the whole slide and the darker is the
# annotation-restricted one, so the two restrictions read as two separate visual
# axes rather than four unrelated categories. All four are Okabe-Ito, so the warm/
# cool split survives the common colour-vision deficiencies where the light/dark
# ladder within each pair carries the scope on its own. Worst pair is the within-
# pair one, all_wholeslide vs all_annotation, at dE 26.4 in CIE Lab — just over the
# ~25 a small mark needs. Re-check with convertColor(..., "Lab") before changing
# either blue.
DENOM_COLS <- c(
  all_wholeslide   = "#56B4E9",   # sky blue    — every cell, no polygon (the reference)
  all_annotation   = "#0072B2",   # blue        — every cell inside the union polygon
  tumor_wholeslide = "#E69F00",   # orange      — tumour cells, whole slide
  tumor_annotation = "#D55E00"    # vermillion  — tumour cells inside the union polygon
)

# Spelled for a reader, not for the column. "Tumour" with the British u to match
# LINEAGE_LABELS, and the scope named by what it restricts to rather than by the
# directory it came from — `annotation_all` is a path, "inside annotation" is a fact.
DENOM_LABELS <- c(
  all_wholeslide   = "all cells, whole slide",
  all_annotation   = "all cells, inside annotation",
  tumor_wholeslide = "tumour cells, whole slide",
  tumor_annotation = "tumour cells, inside annotation"
)

# Order widest -> narrowest, which is also nesting order (tumour ∩ annotation is a
# subset of both of its parents), so a legend reads as a progressive restriction.
denom_factor <- function(x) factor(as.character(x), levels = names(DENOM_COLS))

.denom_labels <- function(breaks) {
  b <- as.character(breaks)
  ifelse(b %in% names(DENOM_LABELS), DENOM_LABELS[b], b)
}
scale_colour_denominator <- function(..., drop = TRUE, name = NULL, labels = .denom_labels)
  scale_colour_manual(values = DENOM_COLS, drop = drop, name = name,
                      labels = labels, na.value = "grey85", ...)
scale_fill_denominator <- function(..., drop = TRUE, name = NULL, labels = .denom_labels)
  scale_fill_manual(values = DENOM_COLS, drop = drop, name = name,
                    labels = labels, na.value = "grey85", ...)
scale_color_denominator <- scale_colour_denominator

scale_colour_arm <- function(..., drop = TRUE, name = NULL)
  scale_colour_manual(values = ARM_KIND_COLS, drop = drop, name = name, ...)
scale_fill_arm <- function(..., drop = TRUE, name = NULL)
  scale_fill_manual(values = ARM_KIND_COLS, drop = drop, name = name, ...)
scale_color_arm <- scale_colour_arm

# --- FlowPath panel titles ---------------------------------------------------
# Facet strips used to print the COLUMN name, which is three different spellings
# of one population across three pages: the clinical page counts CD45+ cells into
# `cd45_over_inside`, the molecular page counts the same cells into
# `CD45_posfrac`, and a reader comparing the two figures had to know that. Every
# FlowPath panel is titled by WHAT WAS COUNTED instead — "FlowPath CD45+ cells" —
# so the two strips read identically.
#
# The map only needs the names that are not already a bare marker; anything else
# falls through to "<name>+", which is what a `<MARKER>_posfrac` column means.
# An unmapped name is therefore still drawn, not turned into NA: a blank strip is
# a worse failure than a slightly wrong one, because it looks like a rendering bug
# rather than a missing entry.
FLOWPATH_PANEL_POPULATIONS <- c(
  cd45_over_inside    = "CD45+",
  cd3cd45_over_inside = "CD3+ CD45+",
  gzmb_nk_over_inside = "GZMB+ NK",
  frac_CD8T           = "CD8 T",
  frac_CD4T           = "CD4 T",
  frac_Treg           = "Treg",
  frac_NK             = "NK",
  Immune_other        = "immune (other)")

flowpath_panel_label <- function(x) {
  b   <- as.character(x)
  pop <- ifelse(b %in% names(FLOWPATH_PANEL_POPULATIONS),
                FLOWPATH_PANEL_POPULATIONS[b],
                paste0(sub("_(posfrac|z)$", "", b), "+"))
  stats::setNames(paste0("FlowPath ", pop, " cells"), b)
}

# The same thing as a ggplot2 labeller, so a call site reads
# `facet_wrap(~ marker, labeller = flowpath_labeller())` rather than repeating
# the as_labeller() wrapping at every one of the six call sites.
flowpath_labeller <- function() ggplot2::as_labeller(flowpath_panel_label)

# --- view clipping ------------------------------------------------------------
# "Remove the outliers" has two meanings and only one of them is honest here.
# Dropping the rows recomputes the boxes, the medians and the stated n, so the
# figure then disagrees with the summary table printed beside it — and a reader has
# no way to see that it does. These two crop the VIEW instead: every statistic is
# still computed over every observation, the axis just stops before the tail.
#
# A clipped axis that does not say how many points are outside it is a lie by
# omission, so the note is a separate function rather than optional: a call site
# that clips and forgets the caption is a diff you can spot.
#
# `by` is the grouping the figure already draws on its x axis, and passing it is
# not optional cosmetics. A registration ladder's stages differ by an order of
# magnitude BY DESIGN — native is tens of microns, micro is sub-micron — so a fence
# computed over the pooled vector calls the entire `native` stage an outlier and
# crops away the very comparison the figure exists to make. Computing the fence
# WITHIN each stage and taking the widest keeps every stage's bulk on screen while
# still cropping the one slide that blew up.
#
# The cut is Tukey's upper fence (Q3 + k*IQR), i.e. the same definition the boxplot
# beside it already uses to decide what counts as an outlier — so "clipped" and
# "drawn as an outlier point" agree instead of being two different thresholds. The
# view then ends on a real observation rather than on a quantile that may land in
# the middle of the gap.
clip_upper_ylim <- function(x, by = NULL, k = 1.5) {
  keep <- is.finite(x)
  x <- x[keep]
  if (!length(x)) return(NULL)
  g <- if (is.null(by)) rep("", length(x)) else as.character(by)[keep]
  fence <- function(v) {
    # Under four observations there is no distribution to call anything an outlier
    # against, so that group asks for its full range and the max below ignores it.
    if (length(v) < 4) return(max(v))
    q <- unname(stats::quantile(v, c(.25, .75)))
    inside <- v[v <= q[[2]] + k * (q[[2]] - q[[1]])]
    if (!length(inside)) max(v) else max(inside)
  }
  hi <- max(vapply(split(x, g), fence, numeric(1)))
  lo <- min(x)
  if (!is.finite(hi) || hi >= max(x)) return(NULL)   # nothing outside any fence
  c(lo, hi + 0.05 * (hi - lo))                       # pad so the top point clears the edge
}

# `id` is what one unit of `unit` IS. Without it the count is rows, and on any
# figure where a slide contributes one row per stage that reads "3 slides above the
# view" for a single slide seen three times — a number a reader can check against
# the stated n and find wrong.
clip_upper_note <- function(x, ylim, unit = "points", id = NULL) {
  if (is.null(ylim)) return(NULL)
  above <- is.finite(x) & x > ylim[[2]]
  n <- if (is.null(id)) sum(above) else length(unique(id[above]))
  if (n == 0) return(NULL)
  # "1 slides above the view" is the kind of thing a reviewer circles.
  unit <- if (n == 1) sub("s$", "", unit) else unit
  sprintf("y axis clipped: %d %s above the view; every box, median and n is computed over all of them",
          n, unit)
}

# Join a figure's standing caption to whatever notes this particular render earned.
# NULLs drop out, so a call site can pass clip_upper_note() unconditionally.
caption_with <- function(...) {
  parts <- Filter(nzchar, unlist(list(...)))
  if (!length(parts)) return(NULL)
  paste(parts, collapse = " · ")
}

# --- print sizing ------------------------------------------------------------
# What keeps type the SAME SIZE on paper across figures is not a constant
# fig.width — it is a constant SHRINK FACTOR between the rendered figure and the
# column it is placed in:
#
#     scale = target_column_mm / (fig.width_in * 25.4)
#
# The repo's established double-column convention is fig.width = 9in placed in a
# 180mm slot, i.e. scale = 0.787, which lands the 10pt base type at 7.9pt — inside
# the 5-8pt journal range. FIG_SCALE pins that factor, and the widths below are
# DERIVED from it, so a single-column figure is not simply "a narrower 9in figure"
# (which would print its type ~2x too large) but the width that shrinks by the
# same 0.787.
#
# Rendering wider than the target is deliberate: it also gives the workflowr site
# a readable raster, so one chunk serves both the website and the manuscript.
FIG_SCALE <- 0.787

# Journal column widths in mm (Nature-style: 85 single, 180 double).
FIG_COLUMN_MM <- c(single = 85, oneandhalf = 114, double = 180)

# Width in INCHES for a knitr chunk targeting the given column. Round to 0.05in so
# the numbers that end up in the Rmds stay readable.
fig_width <- function(column = c("double", "single", "oneandhalf")) {
  column <- match.arg(column)
  round(FIG_COLUMN_MM[[column]] / FIG_SCALE / 25.4 / 0.05) * 0.05
}

# Height for a target aspect ratio (height/width). Capped at the 230mm ceiling a
# figure may not exceed on a printed page, expressed at the same shrink factor.
fig_height <- function(column = c("double", "single", "oneandhalf"), aspect = 0.62) {
  w   <- fig_width(column)
  cap <- 230 / FIG_SCALE / 25.4
  round(min(w * aspect, cap) / 0.05) * 0.05
}

# A panel that will be PLACED BY HAND (paper_figures.R emits one PDF per panel and
# the figure is assembled in Affinity) has no column width to snap to, so the ladder
# above cannot help it. Declare the width the panel will occupy in the finished
# figure instead, in mm, and this returns the inches to render it at:
#
#   ```{r fig5b, fig.width = panel_width(68), fig.height = panel_height(72)}
#
# Two panels sharing a 180mm row declare 110 and 66; a full-width panel declares
# 180. Because FIG_SCALE * 25.4 is almost exactly 20, the conversion is "mm / 20",
# which makes a mismatched panel visible at a glance in the chunk header. This is
# the ONLY way to keep type uniform in a hand-assembled figure: an inch literal
# encodes a placement nobody wrote down, so nothing can check it.
panel_width <- function(mm) round(mm / FIG_SCALE / 25.4 / 0.05) * 0.05
panel_height <- panel_width

# --- reporting n -------------------------------------------------------------
# Every figure has to state its n somewhere. Two helpers, because there are two
# situations and they want different answers.
#
# `label_n(x)` is for a CATEGORICAL AXIS: it returns a labeller that appends the
# per-level count to each tick, so an unbalanced design is visible in the figure
# rather than buried in the caption.
#
#     scale_x_discrete(labels = label_n(df$lineage))
#
# It counts from the vector handed to it (the plotted data), not from the breaks,
# so it reports what was actually drawn. Levels present as a break but absent from
# the data are labelled without a count instead of "(n = NA)".
label_n <- function(x, sep = "\n") {
  counts <- table(as.character(x[!is.na(x)]))
  function(breaks) {
    n <- as.integer(counts[as.character(breaks)])
    ifelse(is.na(n) | !nzchar(as.character(breaks)),
           as.character(breaks),
           sprintf("%s%s(n = %d)", breaks, sep, n))
  }
}

# `n_note()` is the fallback for a figure with NO categorical axis (a scatter, a
# density, a heatmap): a phrase to paste into the subtitle. Naming the unit is
# required, because "n = 24" is ambiguous in this project — 24 patients, 24 slides
# and 24 runs are all plausible and mean different things.
n_note <- function(n, unit = "patients") {
  n <- if (length(n) > 1L) length(unique(n[!is.na(n)])) else as.integer(n)
  sprintf("n = %s %s", format(n, big.mark = ","), unit)
}

# Append `n_note()` to a subtitle, keeping the subtitle readable when it is NULL.
with_n <- function(subtitle, n, unit = "patients") {
  note <- n_note(n, unit)
  if (is.null(subtitle) || !nzchar(subtitle)) note else paste0(subtitle, " · ", note)
}

# --- geom defaults -----------------------------------------------------------
# The theme cannot reach inside a geom, so a default-sized geom_point (size 1.5,
# ~4pt on paper) and 0.5-linewidth lines stay chunky no matter how good the theme
# is. These bring the marks themselves down to print scale. Explicit sizes at a
# call site still win, so this only moves the plots nobody tuned by hand.
paper_geom_defaults <- function() {
  set <- function(geom, vals)
    try(update_geom_defaults(geom, vals), silent = TRUE)
  set("point",   list(size = 1.2, stroke = 0.3))
  set("line",    list(linewidth = pt_line(0.75)))
  set("path",    list(linewidth = pt_line(0.75)))
  set("segment", list(linewidth = pt_line(0.75)))
  set("hline",   list(linewidth = pt_line(0.5), colour = "grey40"))
  set("vline",   list(linewidth = pt_line(0.5), colour = "grey40"))
  set("abline",  list(linewidth = pt_line(0.5)))
  set("smooth",  list(linewidth = pt_line(1.0)))
  # Outline-only marks need more stroke than a data line: at 0.5pt a pale fill
  # colour (the light blue in the hot/cold scale, say) washes out to grey while
  # the thicker median segment inside the same box still reads as its true hue.
  set("boxplot", list(linewidth = pt_line(0.9)))
  set("bar",     list(linewidth = pt_line(0.9)))
  set("col",     list(linewidth = pt_line(0.9)))
  set("errorbar",list(linewidth = pt_line(0.75)))
  # In-panel annotation text at ~7pt. Deliberately does NOT set `family`: these
  # defaults feed annotate("text", ...) and geom_text() alike, and a family the pdf
  # device cannot resolve aborts the knit there (see the font section above).
  # Leaving it unset inherits "" = the device default, same as the theme.
  set("text",    list(size = pt_text(7), colour = "black"))
  set("label",   list(size = pt_text(7), colour = "black"))
  invisible(NULL)
}

# --- apply -------------------------------------------------------------------
# Side effects on source(): this is what makes the style automatic rather than
# something each Rmd has to remember. Idempotent, so sourcing twice is harmless.
theme_set(theme_paper())
paper_geom_defaults()

# Any discrete scale the author did NOT specify falls back to these instead of
# ggplot2's evenly-spaced hue rainbow (which is neither colourblind-safe nor
# printable in greyscale). Passed as a LIST of palettes: ggplot2 takes the first
# one long enough for the variable's levels, so <=8 categories get `oi` and the
# 9..16 range gets `oi_ext` rather than silently reverting to hue. (The list
# entries must be character vectors — ggplot2 rejects functions here.) Beyond 16
# levels ggplot2 does fall back to its own scale; that is a signal the plot needs
# faceting or direct labels, not a longer palette.
options(ggplot2.discrete.colour = list(unname(oi), unname(oi_ext)),
        ggplot2.discrete.fill   = list(unname(oi), unname(oi_ext)),
        ggplot2.continuous.colour = function(...) scale_colour_seq(...),
        ggplot2.continuous.fill   = function(...) scale_fill_seq(...))
