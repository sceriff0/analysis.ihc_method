# =============================================================================
# Standalone per-plot PDF export for the workflowr site.
#
# workflowr renders each chunk's plots to PNGs under the site figure path
# (docs/figure/<Rmd>/) for the website. When a chunk also sets
# `dev = c("png", "pdf")`, knitr emits a vector PDF of EVERY plot as well — one
# file per plot, so loops and multi-plot chunks each yield their own PDF, and
# base-R plots are captured just like ggplots.
#
# `export_pdf_figures(slug)` is meant to run in a final `include=FALSE` chunk of
# an analysis. It collects those generated PDFs and copies each into
# output/figures/<slug>/ so every figure is available as a single publication
# PDF, without disturbing the PNGs the HTML points at.
#
# It is deliberately fail-safe: any error is caught and downgraded to a message
# so a figure-export hiccup can never abort a knit. If no PDFs are found (e.g.
# `dev` was not set to include "pdf") it says so and does nothing.
#
# Dependencies: base R + knitr + here only (no tidyverse/sf), so it can be
# sourced from any analysis, including ones that do not load validation_helpers.
# It also sources code/placeholders.R, which attaches nothing on source (its dplyr /
# tibble calls are namespaced and run only in placeholder mode).
# =============================================================================

# Headless machines (e.g. HPC cluster nodes with no X server) cannot start the
# default X11-based png() device — knitr then fails with "unable to start device
# PNG". Route bitmap output through cairo, which needs no X11, whenever cairo is
# compiled in. Runs on source(), i.e. in the setup chunk before any plot, and
# only touches png()/jpeg() — pdf() is unaffected. No-op where cairo is absent.
if (isTRUE(capabilities("cairo"))) options(bitmapType = "cairo")

# PLACEHOLDER MODE (code/placeholders.R). Sourced here because every page sources
# this file: that is what puts placeholder_callout() and the export guard on every
# page without a second source line. placeholders.R attaches no package on source.
source(here::here("code", "placeholders.R"))

# With placeholder mode ON the export is DIVERTED, never merged: nothing is written
# to output/figures/<slug>/ (whose PDFs are what a manuscript is assembled from), the
# PDFs go to output/placeholders/figures/<slug>/<name>_PLACEHOLDER.pdf instead, and
# the page's sidecar list of synthetic values is written beside them. The divert is
# unconditional on the mode rather than on whether this page drew a synthetic point:
# a figure-level count cannot see every plot, and a wrong guess here is the one that
# puts a synthetic point into a paper figure.
export_pdf_figures <- function(slug, out_root = here::here("output", "figures"),
                               placeholder_root = placeholder_root()) {
  tryCatch({
    fp <- knitr::opts_chunk$get("fig.path")   # e.g. "figure/clinical_flowpath.Rmd/"

    # Candidate directories that may hold this knit's PDF figures. workflowr's
    # runtime fig.path is relative and the exact on-disk location varies by build
    # stage, so search the obvious roots and keep whatever exists.
    fig_roots <- character(0)
    if (!is.null(fp) && nzchar(fp)) {
      fig_roots <- c(fig_roots, fp, file.path(getwd(), fp), here::here(fp))
    }
    for (base in c(here::here("docs", "figure"),
                   here::here("analysis", "figure"),
                   here::here("figure"))) {
      if (dir.exists(base)) {
        subs <- list.dirs(base, recursive = FALSE)
        # keep only sub-directories whose name refers to this analysis
        fig_roots <- c(fig_roots, subs[grepl(slug, basename(subs), fixed = TRUE)])
      }
    }
    fig_roots <- unique(fig_roots[dir.exists(fig_roots)])

    pdfs <- unlist(lapply(fig_roots, function(d)
      list.files(d, pattern = "\\.pdf$", full.names = TRUE)), use.names = FALSE)
    pdfs <- unique(pdfs)

    if (!length(pdfs)) {
      message('export_pdf_figures("', slug, '"): no PDF figures found — ',
              'is `dev = c("png", "pdf")` set in this Rmd\'s setup chunk?')
      return(invisible(character(0)))
    }

    if (placeholder_mode()) {
      placeholder_write_sidecar(slug, root = placeholder_root)
      dest  <- file.path(placeholder_root, "figures", slug)
      names <- paste0(tools::file_path_sans_ext(basename(pdfs)), "_", PLACEHOLDER_TAG, ".pdf")
      dir.create(dest, recursive = TRUE, showWarnings = FALSE)
      ok <- file.copy(pdfs, file.path(dest, names), overwrite = TRUE)
      message(sprintf(paste0('export_pdf_figures("%s"): PLACEHOLDER MODE — refused %s; ',
                             'copied %d/%d PDF(s) to %s as *_%s.pdf'),
                      slug, file.path(out_root, slug), sum(ok), length(pdfs), dest,
                      PLACEHOLDER_TAG))
      return(invisible(file.path(dest, names)[ok]))
    }

    dest <- file.path(out_root, slug)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(pdfs, file.path(dest, basename(pdfs)), overwrite = TRUE)
    message(sprintf('export_pdf_figures("%s"): copied %d/%d PDF(s) to %s',
                    slug, sum(ok), length(pdfs), dest))
    invisible(file.path(dest, basename(pdfs))[ok])
  }, error = function(e) {
    message('export_pdf_figures("', slug, '") skipped: ', conditionMessage(e))
    invisible(character(0))
  })
}
