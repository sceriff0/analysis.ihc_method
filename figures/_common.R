# =============================================================================
# figures/_common.R  —  the journal spec, and the two things every figure script
#                       does to a panel before patchwork sees it.
#
# WHAT THIS IS FOR. `code/paper_figures.R` and `code/registration_arms.R` build
# panels for the WEBSITE, where workflowr renders them wide and the reader zooms.
# This directory builds the same panels for PAPER, where the figure is placed at a
# fixed column width and never resized again. Those are different jobs and the
# difference is not cosmetic — see SIZING below.
#
# WHY THIS FILE IS NOT `theme_set(theme_classic(base_size = 8, ...))`.
# The brief asked for exactly that line. It is the wrong line for this repo, twice:
#
#   1. `theme_classic()` REPLACES the active theme wholesale, so it silently
#      discards theme_paper()'s outward ticks, black (not grey) text, strip and
#      legend geometry. plot_theme.R's header documents this and the repo rule is
#      "never call theme_classic() in an analysis". theme_paper(base_size = 8) IS
#      theme_classic(base_size = 8) plus the corrections, so it is the same request,
#      honoured.
#   2. `base_family = "Arial"` aborts any knit that renders a PDF: "Arial" is an OS
#      face, not a face registered with grDevices::pdf(), whose base-14 set is
#      Type 1. paper_family() exists to catch exactly this and degrade to "" (the
#      device's own Helvetica-metric default, which is what Arial is standing in
#      for anyway). See the FONT block in plot_theme.R — that failure cost a build.
#
# WHY NOT scale_colour_manual(values = oi) EVERYWHERE, either.
# The brief asked for one Okabe-Ito palette applied by scale_*_manual(). The palette
# here IS Okabe-Ito, but applied through the NAMED scales in plot_theme.R
# (scale_*_lineage, scale_*_arm, hotcold_cols). `oi` is an UNNAMED vector, so ggplot
# assigns it by factor POSITION: a panel that happens to omit one level shifts every
# colour after it, which is how CD8T was vermillion in one panel and blue in another.
# The named scales are the same colours keyed by MEANING, so a population keeps its
# hue across panels that do not share a level set. Same palette, correct binding.
#
# SIZING — THE ONE REAL DEPARTURE FROM THE HOUSE STYLE.
# plot_theme.R's fig_width() renders WIDER than the target column and relies on the
# journal shrinking the page by 0.787, which lands its 10pt base type at 7.9pt on
# paper. This directory exports at the FINAL PRINTED SIZE and is never resized, so
# there is no shrink to absorb: the base size must be the size that reaches paper.
# Hence base_size = 8 with true-mm widths, and hence NOT fig_width() here. Both
# routes land 8pt type on paper; they just put the factor in different places.
# tests/testthat/test-plot-style.R bans literal fig.width in .Rmd CHUNK HEADERS —
# these are ggsave() calls in .R scripts, which is the sanctioned other route.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

here_root <- tryCatch(here::here(), error = function(e) normalizePath("."))
source(file.path(here_root, "code", "plot_theme.R"))

# --- The journal -------------------------------------------------------------
# Medical Image Analysis. Widths are the three column measures; the figure picks
# one and NEVER a number in between, because a 165mm figure is scaled to 190 or 140
# by the typesetter and the type size stops being what this file says it is.
JOURNAL <- "Medical Image Analysis"
MM      <- c(one_col = 90, one_half = 140, two_col = 190)
MAX_H   <- 230          # printed ceiling, matches plot_theme.R's FIG_MAX_H_MM
# 500, not 300: Fig 4(a) and Fig 5(a) put a microscopy image next to vector line
# art, and 300dpi visibly steps the diagonal edges of the line art when the raster
# half forces the whole page to be rasterised at the image's resolution.
DPI     <- 500
BASE_PT <- 8

# Tag style. MIA parenthesises; Genome Medicine would be bare bold lowercase.
TAG <- list(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

# --- The theme ---------------------------------------------------------------
# theme_paper() at the printed base size, plus the brief's two hairline specs. This
# is a bare theme() layered ON TOP, which is additive — the house theme survives.
theme_set(
  theme_paper(base_size = BASE_PT) +
    theme(axis.line  = element_line(linewidth = 0.3, colour = "black",
                                    lineend = "square"),
          axis.ticks = element_line(linewidth = 0.3, colour = "black"))
)

# --- Panel -> assembly-ready -------------------------------------------------
# Two jobs, both of which have to happen to the OBJECT rather than to the caller's
# memory, because these panels are built by functions shared with the website and
# those functions must keep their titles for the website.
#
# 1. DROP THE TITLES. plot_theme.R: "a tag and a plot.title share the top-left slot
#    and WILL overlap ... when composing the final figure, set labs(title = NULL,
#    subtitle = NULL) and let patchwork place the tags." The caption goes too — it
#    is the manuscript legend's job, and ARM_CAPTION repeated under two adjacent
#    panels is the same sentence twice.
# 2. SHRINK THE MARKS. The website panels use size 2.2-2.6 points, sized for a plot
#    rendered at 9in and read on screen. At 190mm placed, that is a blob. The brief
#    says ~0.8. This walks the built layers rather than rebuilding the plots,
#    so WHAT IS PLOTTED is untouched — only how big the ink is.
#
# A layer that MAPS size (aes(size = ...)) is left alone: overriding it would
# silently delete an encoded variable. Nothing in Fig 4 or 5 maps size today; this
# guard is here so that stays true if someone adds a bubble panel later.
# Everything that draws a STROKE rather than a mark or a glyph. Listed rather than
# inferred, because ggplot has no "is this geom a line" predicate and guessing from
# whether `linewidth` is a valid aesthetic catches GeomPoint too (its linewidth is
# the point's border).
LINE_GEOMS <- c("GeomBoxplot", "GeomErrorbar", "GeomCrossbar", "GeomLinerange",
                "GeomLine", "GeomPath", "GeomSegment", "GeomStep", "GeomSmooth",
                "GeomHline", "GeomVline", "GeomAbline", "GeomRect", "GeomPolygon")

for_panel <- function(p, point_size = 0.8, line_pt = 0.3, text_pt = BASE_PT) {
  if (is.null(p)) return(NULL)
  p <- p + labs(title = NULL, subtitle = NULL, caption = NULL)

  p$layers <- lapply(p$layers, function(ly) {
    cls    <- class(ly$geom)
    mapped <- names(ly$mapping %||% list())
    is_pt  <- any(c("GeomPoint", "GeomJitter") %in% cls)
    is_tx  <- any(c("GeomText", "GeomLabel") %in% cls)
    is_ln  <- any(LINE_GEOMS %in% cls)

    if (is_pt && !"size" %in% mapped)      ly$aes_params$size      <- point_size
    if (is_tx && !"size" %in% mapped)      ly$aes_params$size      <- pt_text(text_pt)
    # Unconditionally, not "only if already set": a boxplot the builder never gave a
    # linewidth to keeps ggplot's 0.5 default, which at 8pt base draws a box heavier
    # than the axis it sits on. That is what made the Fig 4 summary boxes read as
    # solid bars at print size.
    if (is_ln && !"linewidth" %in% mapped)  ly$aes_params$linewidth <- line_pt
    ly
  })
  p
}

# Every panel object is cached, so a co-author can rebuild one panel of a figure
# without re-running the pipeline that produced its data.
save_panel <- function(p, id) {
  if (is.null(p)) return(invisible(NULL))
  dir.create(file.path(here_root, "figures", "panels"), recursive = TRUE,
             showWarnings = FALSE)
  saveRDS(p, file.path(here_root, "figures", "panels", paste0(id, ".rds")))
  invisible(p)
}

# --- Image panels ------------------------------------------------------------
# An author-supplied micrograph, as a ggplot so patchwork can place it in the grid
# and tag it like any other panel. Returns a labelled placeholder (not an error and
# not NULL) when the file is absent, so the figure still assembles at the right
# geometry and the gap is visible rather than silent.
image_panel <- function(path, placeholder = "panel to supply") {
  if (!file.exists(path)) {
    return(ggplot() +
             annotate("text", 0, 0, label = placeholder, size = pt_text(BASE_PT),
                      colour = "grey45") +
             annotate("rect", xmin = -1, xmax = 1, ymin = -.5, ymax = .5,
                      fill = NA, colour = "grey80", linewidth = 0.3) +
             coord_cartesian(xlim = c(-1, 1), ylim = c(-.5, .5)) +
             theme_void())
  }
  magick::image_read(path) |> magick::image_ggplot()
}

# --- Export ------------------------------------------------------------------
# PDF is the submission format and TIFF is the fallback the production system asks
# for. Both are written at the SAME mm, so they are the same figure — the only
# difference is that one is resolution-independent.
#
# cairo_pdf, not pdf(): it embeds the font as a subset rather than referencing a
# base-14 name, so the file renders identically on a typesetter that has no
# Helvetica. LZW on the TIFF because the journal's 10MB cap is per-file and an
# uncompressed 190x230mm at 500dpi is ~50MB.
export_figure <- function(fig, name, width_mm, height_mm, dpi = DPI) {
  stopifnot(height_mm <= MAX_H)
  dir.create(file.path(here_root, "figures"), showWarnings = FALSE)
  pdf_p  <- file.path(here_root, "figures", paste0(name, ".pdf"))
  tiff_p <- file.path(here_root, "figures", paste0(name, ".tiff"))

  ggsave(pdf_p, fig, width = width_mm, height = height_mm, units = "mm",
         device = grDevices::cairo_pdf)
  ggsave(tiff_p, fig, width = width_mm, height = height_mm, units = "mm",
         dpi = dpi, device = grDevices::tiff, compression = "lzw", type = "cairo")

  sz <- function(f) round(file.info(f)$size / 1024^2, 2)
  # Reported, not silently accepted: over the cap the figure has to be re-exported
  # at a lower dpi and the author needs to know before submission, not after.
  if (sz(tiff_p) > 10)
    warning(name, ".tiff is ", sz(tiff_p), " MB — over the 10 MB cap. ",
            "Re-export with dpi = 300.", call. = FALSE)

  message(sprintf("%-6s %3.0f x %3.0f mm | PDF %5.2f MB | TIFF %5.2f MB @ %d dpi",
                  name, width_mm, height_mm, sz(pdf_p), sz(tiff_p), dpi))
  invisible(c(pdf = pdf_p, tiff = tiff_p))
}
