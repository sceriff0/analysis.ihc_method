# =============================================================================
# paper_figures.R  —  the panels the manuscript figures are cut from
#
# Every other figure script in code/ answers an ANALYSIS question and is read on
# the website. This one answers a LAYOUT question: it emits the specific panels
# named in the figure legends, at publication proportions, one panel per plot, so
# they can be exported as PDFs and assembled by hand in Affinity. Nothing here is
# new science — each function is a re-cut of a quantity some analysis page already
# computes, in the shape the legend asks for.
#
# ONE PANEL PER PLOT, NEVER A COMPOSITE. There is no patchwork/cowplot in the
# lockfile, and that is the right constraint here rather than a limitation: a
# whole-slide map and its inset are two PDFs, and the person assembling the figure
# places, scales and labels them. A composite baked in R would have to be taken
# apart again.
#
# WHY THESE FUNCTIONS AND NOT OTHERS. The legends ask for four things no analysis
# page currently draws:
#   paper_phenotype_map()          Fig 3(d) and Fig 5(a) — cells coloured by call
#   paper_immune_fraction_hotcold() Fig 5(b) — CD45+/all cells, hot vs cold
#   paper_deconv_scatter()          Fig 5(c) — one method, no fit, no coefficient
#   paper_lineage_table()           Additional file 4 — the mapping, as data
# The registration arm figures the legends call Fig 4(b)/(c) live in
# registration_accuracy_plots.R instead, because they belong to that page's data.
#
# Depends on validation_helpers.R (and through it cell_tables.R + plot_theme.R).
# sf is only needed when an annotation outline is drawn.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

# --- Spatial phenotype map (Fig 3d, Fig 5a) ----------------------------------
# A "legible subset of the tree" is the legend's phrase and it is load-bearing: a
# map with fourteen colours reads as noise at figure size. LEGIBLE_LINEAGES is that
# subset, ordered so the immune populations the paper argues about are the ones
# that get a distinct hue, and everything else collapses into one grey "other".
LEGIBLE_LINEAGES <- c("Tumor", "CD8T", "CD4T", "Treg", "NK", "Immune_other", "Stroma")

LINEAGE_COLS <- c(
  Tumor        = "#B0B7BE",   # grey: the background the immune cells sit on
  CD8T         = "#D55E00",
  CD4T         = "#0072B2",
  Treg         = "#CC79A7",
  NK           = "#009E73",
  Immune_other = "#E69F00",
  Stroma       = "#7F8C8D",
  other        = "#DEE2E6"
)

# A "nice" scale-bar length: the largest of 10/25/50/100/... µm that still fits in
# a sixth of the field. Picking a round number MATTERS — a 137 µm bar is unreadable
# as a bar, and hard-coding one length makes it invisible on a whole slide and
# wider than the field on an inset.
.nice_bar_um <- function(span_um) {
  target <- span_um / 6
  cand   <- as.vector(outer(c(1, 2.5, 5), 10^(0:5)))
  cand   <- cand[cand <= target]
  if (!length(cand)) return(signif(target, 1))
  max(cand)
}

# Cells coloured by phenotype call, in the geojson PIXEL frame.
#
# `zoom = c(xmin, xmax, ymin, ymax)` in pixels cuts the inset — the SAME function
# call, so the inset cannot drift from the map it is an inset of. Point size is
# left free by default and scaled to the field, because the size that reads as
# "one cell" on a whole slide is a smear at inset zoom.
#
# Y IS REVERSED. Image coordinates put the origin top-left and y increasing
# downward; a plot drawn with ggplot's default y-up is the slide upside down,
# which is not obvious from the plot but is obvious the moment it sits next to the
# image it came from.
paper_phenotype_map <- function(cells, patient_id = NULL, annots = NULL,
                                zoom = NULL, um_per_px = 0.325,
                                point_size = NULL, scale_bar = TRUE,
                                title = NULL) {
  if (!is.null(patient_id))
    cells <- dplyr::filter(cells, slide_key(.data$patient_id) == slide_key(!!patient_id))
  if (nrow(cells) == 0) {
    warning("paper_phenotype_map(): no cells for ", patient_id %||% "the given set")
    return(NULL)
  }

  xy  <- cell_centroids_px(cells, um_per_px)
  lin <- cell_lineage(cells$phenotype_clean)
  df  <- tibble::tibble(x = xy$x, y = xy$y, lineage = lin) |>
    dplyr::filter(is.finite(x), is.finite(y)) |>
    # An unmapped label is a vocabulary gap, not a cell type; it must not get a hue
    # and be read as a population. It joins "other" with everything else off-subset.
    dplyr::mutate(lineage = factor(ifelse(is.na(lineage) | !lineage %in% LEGIBLE_LINEAGES,
                                          "other", lineage),
                                   levels = c(LEGIBLE_LINEAGES, "other")))
  if (!is.null(zoom)) {
    stopifnot(length(zoom) == 4)
    df <- dplyr::filter(df, x >= zoom[1], x <= zoom[2], y >= zoom[3], y <= zoom[4])
    if (nrow(df) == 0) {
      warning("paper_phenotype_map(): the zoom window contains no cells")
      return(NULL)
    }
  }

  span_px <- max(diff(range(df$x)), diff(range(df$y)))
  if (is.null(point_size))
    point_size <- max(0.05, min(1.6, 900 / max(span_px, 1)))

  p <- ggplot(df, aes(x, y, colour = lineage)) +
    geom_point(size = point_size, shape = 16, alpha = .85) +
    scale_colour_manual(values = LINEAGE_COLS, drop = FALSE, name = NULL,
                        guide = guide_legend(override.aes = list(size = 2.5))) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = title %||% (if (!is.null(patient_id)) paste("Patient", patient_id) else NULL),
         x = NULL, y = NULL) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())

  # The pathologist's line, drawn as an outline only — filled, it would hide the
  # cells the panel exists to show.
  if (!is.null(annots) && requireNamespace("sf", quietly = TRUE)) {
    ap <- annots[slide_key(annots$patient_id) == slide_key(patient_id %||% ""), , drop = FALSE]
    if (nrow(ap) > 0) {
      coords <- do.call(rbind, lapply(seq_len(nrow(ap)), function(i) {
        m <- sf::st_coordinates(sf::st_geometry(ap)[i])
        data.frame(x = m[, "X"], y = m[, "Y"], grp = paste(i, m[, "L2"] %||% 1))
      }))
      p <- p + geom_path(data = coords, aes(x, y, group = grp),
                         inherit.aes = FALSE, colour = "grey20",
                         linewidth = pt_line(0.6))
    }
  }

  if (isTRUE(scale_bar)) {
    bar_um <- .nice_bar_um(span_px * um_per_px)
    bar_px <- bar_um / um_per_px
    x1 <- max(df$x) - bar_px; x0 <- max(df$x); y0 <- max(df$y) + span_px * 0.04
    p <- p +
      annotate("segment", x = x1, xend = x0, y = y0, yend = y0,
               linewidth = pt_line(2), colour = "grey10", lineend = "butt") +
      annotate("text", x = (x0 + x1) / 2, y = y0 + span_px * 0.035,
               label = paste0(bar_um, " µm"), size = pt_text(7), colour = "grey10")
  }
  p
}

# --- Immune fraction, hot vs cold (Fig 5b) -----------------------------------
# The legend's panel: CD45+/all-cells per case, two groups, EVERY point drawn.
# With three cases a side the points ARE the evidence, so they are never hidden
# behind the box — the box is scaffolding for the eye, not the result.
#
# NO TEST IS RUN AND NO COEFFICIENT IS PRINTED. That is the legend's claim, and it
# has to be true of the object, not just of the caption: `paired_spearman()` is
# deliberately not called here. The rest of the site reports rho freely; this panel
# is the one place that must not.
paper_immune_fraction_hotcold <- function(metrics, groups,
                                          value_col = "cd45_over_inside",
                                          group_col = "group",
                                          y_lab = "mIF CD45+ / all cells") {
  stopifnot(value_col %in% names(metrics))
  df <- metrics |>
    dplyr::mutate(.pid = slide_key(patient_id)) |>
    dplyr::inner_join(dplyr::mutate(groups, .pid = slide_key(patient_id)) |>
                        dplyr::select(.pid, group = dplyr::all_of(group_col)),
                      by = ".pid") |>
    dplyr::filter(is.finite(.data[[value_col]]), !is.na(group)) |>
    dplyr::mutate(group = hotcold_order(group))
  if (nrow(df) == 0) {
    warning("paper_immune_fraction_hotcold(): nothing to plot after joining groups")
    return(NULL)
  }

  ggplot(df, aes(group, .data[[value_col]], colour = group)) +
    geom_boxplot(outlier.shape = NA, width = .45, colour = "grey35",
                 linewidth = pt_line(0.6)) +
    geom_point(size = 2.6, alpha = .9,
               position = position_jitter(width = .07, height = 0, seed = 1)) +
    scale_colour_manual(values = hotcold_cols(levels(df$group)), guide = "none") +
    labs(x = NULL, y = y_lab,
         subtitle = sprintf("n = %s", paste(sprintf("%d %s", table(df$group),
                                                    names(table(df$group))),
                                            collapse = ", ")))
}

# --- Imaging vs deconvolution (Fig 5c) ---------------------------------------
# One method, faceted by population, raw paired points only.
#
# NO FIT LINE AND NO COEFFICIENT, for the same reason as above and one more: the
# two axes use different denominators (imaging counts cells, deconvolution
# estimates a mixture fraction), so a regression line would invite exactly the
# absolute-agreement reading the legend explicitly disclaims. Ranking is the claim;
# `free` scales per facet are what let ranking be read.
paper_deconv_scatter <- function(paired, method = "quantiseq",
                                 x_lab = "imaging fraction (of all cells)",
                                 y_lab = NULL, label_cases = FALSE) {
  stopifnot(all(c("method", "lineage", "ihc_frac", "score") %in% names(paired)))
  df <- dplyr::filter(paired, tolower(.data$method) == tolower(!!method))
  if (nrow(df) == 0) {
    warning("paper_deconv_scatter(): no rows for method '", method, "' — have: ",
            paste(sort(unique(paired$method)), collapse = ", "))
    return(NULL)
  }
  p <- ggplot(df, aes(ihc_frac, score)) +
    geom_point(size = 2.4, alpha = .85, colour = oi[1]) +
    facet_wrap(~ lineage, scales = "free") +
    labs(x = x_lab, y = y_lab %||% paste(method, "fraction"))
  if (isTRUE(label_cases) && "patient_id" %in% names(df))
    p <- p + geom_text(aes(label = patient_id), size = pt_text(6),
                       vjust = -0.9, colour = "grey35")
  p
}

# --- Additional file 4: the mapping, as data ---------------------------------
# The tier -> cell-type mapping the legends promise, generated from the SAME
# objects the analyses join on, so it cannot drift from what was actually computed.
# `deconv_to_lineage()` resolves four lineages and returns NA for everything else,
# so the deconvolution side of the table is exactly those four by construction.
paper_lineage_table <- function() {
  pheno <- phenotype_lineage_labels |>
    dplyr::transmute(side = "imaging (phenotype call)",
                     label = phenotype_clean, lineage)
  known <- c("T cell regulatory (Tregs)", "T cell CD8+", "T cell CD4+",
             "NK cell", "Macrophage M1", "Macrophage M2", "B cell",
             "Monocyte", "Neutrophil", "Dendritic cell")
  deconv <- tibble::tibble(side = "deconvolution (cell type)", label = known,
                           lineage = deconv_to_lineage(known)) |>
    dplyr::mutate(lineage = tidyr::replace_na(lineage, "(unmapped)"))
  dplyr::bind_rows(pheno, deconv) |> dplyr::arrange(side, lineage, label)
}
