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
# ONE PANEL PER PLOT, NEVER A COMPOSITE — still true of THIS file, for a reason that
# has changed. It used to be that there was no patchwork in the lockfile and figures
# were assembled by hand in Affinity. They are now assembled in code, by
# figures/fig4.R and figures/fig5.R, which is where patchwork lives.
#
# The rule survives the change because the split moved rather than disappeared: these
# functions own WHAT A PANEL SHOWS and the figures/ scripts own HOW PANELS SIT
# TOGETHER. Keeping them apart is what lets the same object be the website panel and
# the manuscript panel — analysis/paper_figures.Rmd prints these directly, titles and
# all, while figures/*.R strips the titles and lays them out. Bake a composite in here
# and the website gets a figure it cannot caption and the assembly scripts get a
# panel they cannot place.
#
# WHY THESE FUNCTIONS AND NOT OTHERS. The legends ask for four things no analysis
# page currently draws:
#   paper_phenotype_map()          Fig 3(d), Fig 5(a) — cells coloured by call
#   paper_phenotype_map_pair()     Fig 5(a) — a case as overview + region inset
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
# LEGIBLE_LINEAGES, LINEAGE_COLS and lineage_legible() now live in plot_theme.R,
# so the composition panels on the clinical page colour CD8T the same vermillion
# this map does. Reach for scale_colour_lineage(), never scale_colour_manual().

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
# `colour_by` picks the reading: "lineage" is the seven-population map (Fig 3d),
# "tumour" is tumour against everything else — the whole-slide overview of Fig
# 5(a) — and "compartment" is tumour / immune / stroma / other, the reading its
# region inset uses. All three come from the same cell_lineage() call collapsed
# by plot_theme.R (lineage_legible(), lineage_compartment()), so a cell that is
# tumour in one view cannot be anything else in another. Colours are the NAMED
# palettes behind scale_colour_lineage() / scale_colour_compartment().
#
# `zoom = c(xmin, xmax, ymin, ymax)` in pixels cuts the inset — the SAME function
# call, so the inset cannot drift from the map it is an inset of. The window is
# applied as coordinate limits, not only as a cell filter, so an annotation outline
# that continues past the edge is clipped rather than stretching the panel. Point
# size is left free by default and scaled to the field, because the size that reads
# as "one cell" on a whole slide is a smear at inset zoom.
#
# `highlight` names one annotation (its `annotation` value in `annots`) to draw
# with a heavier outline: the overview's locator for the region its inset shows.
# `legend_title` names the colour key — needed when two maps with different keys
# share a legend strip, where two untitled keys both starting "Tumour" read as one.
#
# A region-tier cell table lists a cell once per region file it appears in, so the
# cells are de-duplicated on cell_key_cols() first — otherwise a three-region slide
# draws every cell three times and reports three times its n.
#
# Y IS REVERSED. Image coordinates put the origin top-left and y increasing
# downward; a plot drawn with ggplot's default y-up is the slide upside down,
# which is not obvious from the plot but is obvious the moment it sits next to the
# image it came from.
paper_phenotype_map <- function(cells, patient_id = NULL, annots = NULL,
                                zoom = NULL, um_per_px = 0.325,
                                point_size = NULL, scale_bar = TRUE,
                                title = NULL,
                                colour_by = c("lineage", "tumour", "compartment"),
                                highlight = NULL, subtitle = NULL,
                                legend_title = NULL) {
  colour_by <- match.arg(colour_by)
  if (!is.null(patient_id))
    cells <- dplyr::filter(cells, slide_key(.data$patient_id) == slide_key(!!patient_id))
  keys <- intersect(cell_key_cols(cells), names(cells))
  if (length(keys))
    cells <- dplyr::distinct(cells, dplyr::across(dplyr::all_of(keys)), .keep_all = TRUE)
  if (nrow(cells) == 0) {
    warning("paper_phenotype_map(): no cells for ", patient_id %||% "the given set")
    return(NULL)
  }

  xy  <- cell_centroids_px(cells, um_per_px)
  lin <- cell_lineage(cells$phenotype_clean)
  df  <- tibble::tibble(x = xy$x, y = xy$y, lineage = lin) |>
    dplyr::filter(is.finite(x), is.finite(y)) |>
    # An unmapped label is a vocabulary gap, not a cell type; it must not get a hue
    # and be read as a population. lineage_legible() joins it to "other" with
    # everything else off-subset, keeping the full level set so drop = FALSE holds
    # the colours steady across a map and its inset. `call` is what the points are
    # coloured by; `lineage` stays so a caller can always recover the finer label.
    dplyr::mutate(lineage = lineage_legible(lineage),
                  call = switch(colour_by,
                                lineage     = lineage,
                                tumour      = lineage_compartment(lineage, binary = TRUE),
                                compartment = lineage_compartment(lineage)))
  if (!is.null(zoom)) {
    stopifnot(length(zoom) == 4)
    df <- dplyr::filter(df, x >= zoom[1], x <= zoom[2], y >= zoom[3], y <= zoom[4])
    if (nrow(df) == 0) {
      warning("paper_phenotype_map(): the zoom window contains no cells")
      return(NULL)
    }
  }

  # The field is the zoom window when there is one, else the extent of the cells:
  # point size and the scale bar are sized to what is drawn, not to the slide.
  field   <- zoom %||% c(range(df$x), range(df$y))
  span_px <- max(field[2] - field[1], field[4] - field[3])
  if (is.null(point_size))
    point_size <- max(0.05, min(1.6, 900 / max(span_px, 1)))

  scale_fn <- switch(colour_by,
                     lineage = scale_colour_lineage,
                     scale_colour_compartment)
  p <- ggplot(df, aes(x, y, colour = call)) +
    geom_point(size = point_size, shape = 16, alpha = .85) +
    # Unused levels are dropped. It is tempting to keep them so a map and its inset
    # carry identical legends, but the named palette ALREADY guarantees a population
    # is the same colour in both, and a kept-but-empty level draws a labelled key
    # with no swatch beside it — on a manuscript panel that reads as a broken figure.
    # A legend listing only the populations actually present is the accurate one.
    scale_fn(name = legend_title,
             guide = guide_legend(override.aes = list(size = 2.5))) +
    scale_y_reverse() +
    labs(title = title %||% (if (!is.null(patient_id)) paste("Patient", patient_id) else NULL),
         subtitle = subtitle, x = NULL, y = NULL) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
  # coord_fixed() alone for the whole slide; with a zoom the window becomes the
  # panel, edge to edge, so a highlighted outline running past it is cut, not drawn.
  p <- p + (if (is.null(zoom)) coord_fixed()
            else coord_fixed(xlim = zoom[1:2], ylim = zoom[3:4], expand = FALSE))

  # The pathologist's line, drawn as an outline only — filled, it would hide the
  # cells the panel exists to show.
  if (!is.null(annots) && requireNamespace("sf", quietly = TRUE)) {
    ap <- annots[slide_key(annots$patient_id) == slide_key(patient_id %||% ""), , drop = FALSE]
    if (nrow(ap) > 0) {
      ann_names <- if ("annotation" %in% names(ap)) as.character(ap$annotation)
                   else as.character(seq_len(nrow(ap)))
      coords <- do.call(rbind, lapply(seq_len(nrow(ap)), function(i) {
        m <- sf::st_coordinates(sf::st_geometry(ap)[i])
        data.frame(x = m[, "X"], y = m[, "Y"], grp = paste(i, m[, "L2"] %||% 1),
                   hi = !is.null(highlight) && ann_names[i] %in% highlight)
      }))
      p <- p +
        geom_path(data = coords[!coords$hi, , drop = FALSE], aes(x, y, group = grp),
                  inherit.aes = FALSE, colour = "grey20", linewidth = pt_line(0.6))
      if (any(coords$hi))
        p <- p +
          geom_path(data = coords[coords$hi, , drop = FALSE], aes(x, y, group = grp),
                    inherit.aes = FALSE, colour = "grey10", linewidth = pt_line(1.4))
    }
  }

  if (isTRUE(scale_bar)) {
    bar_um <- .nice_bar_um(span_px * um_per_px)
    bar_px <- bar_um / um_per_px
    if (is.null(zoom)) {
      # Below the tissue, where nothing is drawn.
      x0 <- field[2]; x1 <- x0 - bar_px; y0 <- field[4] + span_px * 0.04
      y_lab <- y0 + span_px * 0.035
    } else {
      # Inside the window, bottom-right: outside it would be clipped away.
      pad <- span_px * 0.04
      x0 <- field[2] - pad; x1 <- x0 - bar_px; y0 <- field[4] - pad - span_px * 0.05
      y_lab <- y0 + span_px * 0.035
    }
    p <- p +
      annotate("segment", x = x1, xend = x0, y = y0, yend = y0,
               linewidth = pt_line(2), colour = "grey10", lineend = "butt") +
      annotate("text", x = (x0 + x1) / 2, y = y_lab,
               label = paste0(bar_um, " µm"), size = pt_text(7), colour = "grey10")
  }
  p
}

# Which region the inset shows: the patient's annotated polygon holding the most
# tumour cells, read off the per-annotation metrics table (`n_tumor_inside`, the
# sf point-in-polygon count arm_metrics() already made) rather than recounted here,
# so the inset and the density tables agree on what "most tumour" means. Returns
# the annotation name and its bounding box in geojson pixels, padded by `pad` of
# the box on every side, as the `zoom` paper_phenotype_map() takes; NULL when the
# patient has no polygon with a count (a csv-flag fallback row has no geometry).
tumour_richest_region <- function(per, polys, patient_id, pad = 0.03) {
  stopifnot(all(c("patient_id", "annotation", "n_tumor_inside") %in% names(per)))
  if (is.null(polys) || !requireNamespace("sf", quietly = TRUE)) return(NULL)
  pid  <- slide_key(patient_id)
  pp   <- polys[slide_key(polys$patient_id) == pid, , drop = FALSE]
  cand <- per[slide_key(per$patient_id) == pid &
              per$annotation %in% as.character(pp$annotation) &
              is.finite(per$n_tumor_inside), , drop = FALSE]
  if (nrow(pp) == 0 || nrow(cand) == 0) return(NULL)
  ann  <- as.character(cand$annotation[which.max(cand$n_tumor_inside)])
  bb   <- sf::st_bbox(sf::st_geometry(pp[as.character(pp$annotation) == ann, , drop = FALSE]))
  w <- unname(bb["xmax"] - bb["xmin"]); h <- unname(bb["ymax"] - bb["ymin"])
  list(annotation = ann,
       n_tumor    = max(cand$n_tumor_inside),
       zoom       = unname(c(bb["xmin"] - pad * w, bb["xmax"] + pad * w,
                             bb["ymin"] - pad * h, bb["ymax"] + pad * h)))
}

# The two panels of one Fig 5(a) case, as a LIST — not a composite (see the header).
# `overview` is the whole slide with every annotation outlined, tumour in red and
# everything else grey, the inset's region drawn heavier; `inset` is that region,
# the one with the most tumour cells (tumour_richest_region()), coloured by
# compartment. The callers — analysis/paper_figures.Rmd and figures/fig5.R — put
# them side by side. When no region can be chosen (no polygons, no sf) the inset
# is NULL and the overview stands alone, so the page still knits.
paper_phenotype_map_pair <- function(cells, patient_id, annots = NULL, per = NULL,
                                     um_per_px = 0.325, title = NULL) {
  region <- if (!is.null(per)) tumour_richest_region(per, annots, patient_id) else NULL
  overview <- paper_phenotype_map(cells, patient_id, annots = annots,
                                  um_per_px = um_per_px, colour_by = "tumour",
                                  highlight = region$annotation,
                                  title = title %||% paste("Patient", patient_id),
                                  legend_title = "Whole slide")
  if (!is.null(overview))
    overview <- overview + labs(subtitle = with_n(NULL, nrow(overview$data), "cells"))
  inset <- NULL
  if (!is.null(region)) {
    inset <- paper_phenotype_map(cells, patient_id, annots = annots, zoom = region$zoom,
                                 um_per_px = um_per_px, colour_by = "compartment",
                                 highlight = region$annotation,
                                 title = region_label(region$annotation),
                                 legend_title = "Region")
    if (!is.null(inset))
      inset <- inset + labs(subtitle = with_n(NULL, nrow(inset$data), "cells"))
  }
  list(overview = overview, inset = inset, region = region)
}

# "ANNOTATION_2" is a filename token, not a panel title.
region_label <- function(annotation) {
  a <- as.character(annotation)
  k <- sub("^ANNOTATION_?", "", a, ignore.case = TRUE)
  ifelse(grepl("^ANNOTATION", a, ignore.case = TRUE), paste("Region", k),
         ifelse(tolower(a) == "whole_slide", "Whole slide", a))
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
                                          y_lab = "mIF CD45+ / all cells (unitless, 0-1)",
                                          summary = c("box", "median")) {
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

  # A BOXPLOT, AS THE LEGEND SAYS, WITH EVERY POINT OVER IT. Three cases a side is
  # too few for the hinges to mean much — they are computed from two intervals — so
  # the box is scaffolding for the eye and the points are the evidence. Outliers are
  # not drawn as a separate glyph: every case is already a point, and a second mark
  # for the same case would read as a seventh observation. `summary = "median"` gives
  # the earlier single-rule form if a caller wants it.
  summary_layer <- if (match.arg(summary) == "box") {
    geom_boxplot(outlier.shape = NA, width = .45, colour = "grey35",
                 fill = NA, linewidth = pt_line(0.6))
  } else {
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "errorbar", width = .38, colour = "grey35",
                 linewidth = pt_line(0.6))
  }

  ggplot(df, aes(group, .data[[value_col]], colour = group)) +
    summary_layer +
    geom_point(size = 2.6, alpha = .9,
               position = position_jitter(width = .07, height = 0, seed = 1)) +
    scale_colour_manual(values = hotcold_cols(levels(df$group)), guide = "none") +
    # n rides on the tick labels rather than a subtitle: the groups are unbalanced
    # by design and the reader needs the count attached to the group it describes.
    scale_x_discrete(labels = label_n(df$group)) +
    labs(x = NULL, y = y_lab, subtitle = "every case shown")
}

# --- Imaging vs deconvolution (Fig 5c) ---------------------------------------
# One method, one panel: the imaging fraction on x, the method's estimate on y, a
# point per (patient, population), coloured by population. Only the populations BOTH
# sides resolve are drawn — `paired` is already the inner join on (patient, lineage),
# and the filter below re-asserts that on the frame it is handed, so a stale cache
# with an unmapped row cannot put an unlabelled population on the panel.
#
# NO FIT LINE AND NO COEFFICIENT, for the reason 5(b) gives and one more: the two
# axes use different denominators (imaging counts cells, deconvolution estimates a
# mixture fraction), so a regression line would invite exactly the absolute-agreement
# reading the legend explicitly disclaims. Ranking is the claim.
#
# ONE PANEL, NOT FOUR FACETS. The facetted form put every population on its own free
# axis pair, which made ranking readable within a population but hid that the
# populations live at very different abundances. Sharing the axes shows the four as
# one cloud and lets the reader see, e.g., that CD8 T sits above Treg on both sides.
# Colour comes from scale_colour_lineage(), the named palette every other lineage
# panel uses, so CD8 T is the same vermillion here as in 5(a).
#
# SHAPE IS 5(b)'s AXIS. Pass the `groups` that panel (b) was built on and the points
# take the same hot / cold call on SHAPE, so a case that is hot in (b) is hot in (c)
# without stealing colour from the populations. Without `groups`, the clinical
# `Immuno-phenotype` the cache carries is used and labelled as such — a DIFFERENT
# variable from the hotscore rank, which is why it is the fallback and not the
# default. An older cache with neither plots one shape.
paper_deconv_scatter <- function(paired, method = "quantiseq", groups = NULL,
                                 x_lab = "Imaging fraction of all cells (mIF, 0-1)",
                                 y_lab = NULL, label_cases = FALSE) {
  stopifnot(all(c("method", "lineage", "ihc_frac", "score") %in% names(paired)))
  df <- dplyr::filter(paired, tolower(.data$method) == tolower(!!method))
  if (nrow(df) == 0) {
    warning("paper_deconv_scatter(): no rows for method '", method, "' — have: ",
            paste(sort(unique(paired$method)), collapse = ", "))
    return(NULL)
  }
  df <- df |>
    dplyr::filter(is.finite(.data$ihc_frac), is.finite(.data$score),
                  .data$lineage %in% LEGIBLE_LINEAGES) |>
    dplyr::mutate(lineage = lineage_legible(.data$lineage))
  if (nrow(df) == 0) {
    warning("paper_deconv_scatter(): no finite pairs on a shared population for '",
            method, "'")
    return(NULL)
  }

  # The hot/cold call rides on SHAPE. (b)'s hotscore groups win when given; the
  # clinical category the cached frame may carry is the fallback, and is named as
  # the different variable it is.
  shape_var <- NULL; shape_name <- NULL
  if (!is.null(groups) && nrow(groups) && "patient_id" %in% names(df)) {
    df <- df |>
      dplyr::mutate(.pid = slide_key(.data$patient_id)) |>
      dplyr::left_join(dplyr::transmute(groups, .pid = slide_key(patient_id), group),
                       by = ".pid") |>
      dplyr::select(-".pid")
    if (any(!is.na(df$group))) {
      df$group   <- hotcold_order(df$group)
      shape_var  <- "group"
      shape_name <- "Hotscore group"
    }
  }
  if (is.null(shape_var) && "immuno_phe" %in% names(df) && any(!is.na(df$immuno_phe))) {
    df$immuno_phe <- hotcold_order(df$immuno_phe)
    shape_var     <- "immuno_phe"
    shape_name    <- "Immuno-phenotype (clinical)"
  }

  n_pat <- if ("patient_id" %in% names(df)) df$patient_id else NULL
  p <- ggplot(df, aes(.data$ihc_frac, .data$score, colour = .data$lineage)) +
    (if (!is.null(shape_var))
       geom_point(aes(shape = .data[[shape_var]]), size = 2.4, alpha = .85)
     else geom_point(size = 2.4, alpha = .85)) +
    scale_colour_lineage(name = "Population",
                         guide = guide_legend(override.aes = list(size = 2.5))) +
    labs(x = x_lab, y = y_lab %||% paste(method_label(method), "fraction (0-1)"),
         subtitle = if (!is.null(n_pat)) with_n(NULL, n_pat, "patients") else NULL)
  if (!is.null(shape_var)) {
    # hotcold_order() keeps the spelling of the levels it was given (HOT / Hot / hot),
    # so shapes are keyed by the levels it returns, in its cold -> hot order.
    lv <- levels(df[[shape_var]])
    p <- p + scale_shape_manual(values = stats::setNames(c(16, 17, 15, 18)[seq_along(lv)], lv),
                                na.value = 1, name = shape_name)
  }
  if (isTRUE(label_cases) && "patient_id" %in% names(df))
    p <- p + geom_text(aes(label = .data$patient_id), size = pt_text(6),
                       vjust = -0.9, colour = "grey35", show.legend = FALSE)
  p
}

# The method name as a reader meets it. immunedeconv keys are lower-case slugs; the
# tools have house capitalisation that a legend should respect.
method_label <- function(method) {
  known <- c(quantiseq = "quanTIseq", epic = "EPIC", mcp_counter = "MCP-counter",
             xcell = "xCell", abis = "ABIS", timer = "TIMER",
             consensus_tme = "ConsensusTME", cibersort = "CIBERSORT",
             cibersort_abs = "CIBERSORT (abs.)")
  m <- tolower(method)
  ifelse(m %in% names(known), known[m], method)
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
  # immunedeconv's harmonised spellings, as quanTIseq / EPIC emit them. The CD4
  # helper label carries "(non-regulatory)" — see deconv_to_lineage() for why that
  # word order once cost the CD4T facet.
  known <- c("T cell regulatory (Tregs)", "T cell CD8+", "T cell CD4+ (non-regulatory)",
             "NK cell", "Macrophage M1", "Macrophage M2", "B cell",
             "Monocyte", "Neutrophil", "Myeloid dendritic cell", "uncharacterized cell")
  deconv <- tibble::tibble(side = "deconvolution (cell type)", label = known,
                           lineage = deconv_to_lineage(known)) |>
    dplyr::mutate(lineage = tidyr::replace_na(lineage, "(unmapped)"))
  dplyr::bind_rows(pheno, deconv) |> dplyr::arrange(side, lineage, label)
}

# --- The hot/cold axis (Fig 5a, 5b, 5c) --------------------------------------
# ONE derivation, shared by the website page and figures/fig5.R, because the
# legend's claim is that (a), (b) and (c) are the SAME six cases on the SAME axis.
# Two copies of this ranking — one per caller — is how (b) grouped by hotscore
# while (c) coloured by the clinical category, two different variables inside one
# figure.
#
# `clin` is the clinical CRF, either as the raw sheet (`ID PATIENT`, `HOT score`,
# `Immuno-phenotype`) or already normalised (`patient_id`, `hot_score`,
# `immuno_phe`), or a path to the xlsx. `source = "hot_score"` ranks the
# continuous Foy hotscore and takes the top / bottom k as hot / cold — the
# legend's selection axis. That IS a threshold, on rank rather than on value: the
# legend may say "no threshold on the score" but cannot say "no grouping".
# `source = "immuno_phe"` uses the clinical category instead; it is a different
# variable and a figure built on it must say so.
paper_hotcold_groups <- function(clin, patient_ids = NULL,
                                 source = c("hot_score", "immuno_phe"), k = 3) {
  source <- match.arg(source)
  if (is.character(clin) && length(clin) == 1) {
    if (!file.exists(clin)) return(NULL)
    clin <- readxl::read_excel(clin)
  }
  if (!"patient_id" %in% names(clin)) {
    if (!"ID PATIENT" %in% names(clin))
      stop("paper_hotcold_groups(): need `ID PATIENT` or `patient_id`")
    clin <- clin |>
      dplyr::filter(!is.na(.data[["ID PATIENT"]])) |>
      dplyr::mutate(patient_id = slide_key(.data[["ID PATIENT"]]))
  } else {
    clin <- dplyr::mutate(clin, patient_id = slide_key(patient_id))
  }
  # select(any_of()) so a sheet missing either column loses the column, not the chunk.
  clin <- clin |>
    dplyr::select(patient_id,
                  hot_score  = dplyr::any_of(c("hot_score", "HOT score")),
                  immuno_phe = dplyr::any_of(c("immuno_phe", "Immuno-phenotype"))) |>
    dplyr::distinct(patient_id, .keep_all = TRUE)
  if (!is.null(patient_ids))
    clin <- dplyr::filter(clin, patient_id %in% slide_key(patient_ids))
  if ("hot_score" %in% names(clin))
    clin$hot_score <- suppressWarnings(as.numeric(clin$hot_score))

  if (source == "hot_score") {
    if (!"hot_score" %in% names(clin)) return(NULL)
    sc <- dplyr::filter(clin, is.finite(hot_score))
    n  <- nrow(sc)
    if (n < 2) return(NULL)
    # k never exceeds half the cohort, so hot and cold cannot share a case.
    k <- min(k, floor(n / 2))
    r <- rank(-sc$hot_score, ties.method = "first")
    sc |>
      dplyr::mutate(group = dplyr::case_when(r <= k     ~ "hot",
                                             r > n - k  ~ "cold",
                                             TRUE       ~ NA_character_)) |>
      dplyr::filter(!is.na(group)) |>
      dplyr::select(patient_id, group, hot_score)
  } else {
    if (!"immuno_phe" %in% names(clin)) return(NULL)
    clin |>
      dplyr::filter(!is.na(immuno_phe)) |>
      dplyr::transmute(patient_id, group = as.character(immuno_phe),
                       hot_score = if ("hot_score" %in% names(clin)) hot_score else NA_real_)
  }
}

# One representative case per group for Fig 5(a): the most extreme hotscore in
# each group among the cases that actually have imaging. Derived, not pasted, so
# (a) cannot show a case (b) did not group. Returns c(hot = id, cold = id), with a
# missing group absent rather than NA.
paper_representative_cases <- function(groups, available) {
  if (is.null(groups) || nrow(groups) == 0) return(character(0))
  g <- dplyr::filter(groups, slide_key(patient_id) %in% slide_key(available))
  pick <- function(grp, decreasing) {
    x <- dplyr::filter(g, group == grp)
    if (nrow(x) == 0) return(NULL)
    if ("hot_score" %in% names(x) && any(is.finite(x$hot_score)))
      x <- x[order(x$hot_score, decreasing = decreasing), ]
    x$patient_id[1]
  }
  out <- c(hot = pick("hot", TRUE), cold = pick("cold", FALSE))
  out[!is.na(out)]
}
