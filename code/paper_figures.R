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
#   paper_phenotype_map()          Fig 5(a) — cells coloured by call
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
    # and be read as a population. lineage_legible() joins it to "other" with
    # everything else off-subset, keeping the full level set so drop = FALSE holds
    # the colours steady across a map and its inset.
    dplyr::mutate(lineage = lineage_legible(lineage))
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
    # Unused levels are dropped. It is tempting to keep them so a map and its inset
    # carry identical legends, but the named palette ALREADY guarantees a population
    # is the same colour in both, and a kept-but-empty level draws a labelled key
    # with no swatch beside it — on a manuscript panel that reads as a broken figure.
    # A legend listing only the populations actually present is the accurate one.
    scale_colour_lineage(guide = guide_legend(override.aes = list(size = 2.5))) +
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
# One method, faceted by population, raw paired points only.
#
# NO FIT LINE AND NO COEFFICIENT, for the same reason as above and one more: the
# two axes use different denominators (imaging counts cells, deconvolution
# estimates a mixture fraction), so a regression line would invite exactly the
# absolute-agreement reading the legend explicitly disclaims. Ranking is the claim;
# `free` scales per facet are what let ranking be read.
#
# COLOUR IS 5(b)'s AXIS. Pass the `groups` that panel (b) was built on and the
# points take the same hot / cold call, so a case that is hot in (b) is hot in
# (c). Without `groups`, the clinical `Immuno-phenotype` the cache carries is used
# and labelled as such — a DIFFERENT variable from the hotscore rank, which is why
# it is the fallback and not the default.
paper_deconv_scatter <- function(paired, method = "quantiseq", groups = NULL,
                                 x_lab = "Imaging fraction of all cells (unitless, 0-1)",
                                 y_lab = NULL, label_cases = FALSE,
                                 lineages = comparable_lineages) {
  stopifnot(all(c("method", "lineage", "ihc_frac", "score") %in% names(paired)))
  df <- dplyr::filter(paired, tolower(.data$method) == tolower(!!method))
  if (nrow(df) == 0) {
    warning("paper_deconv_scatter(): no rows for method '", method, "' — have: ",
            paste(sort(unique(paired$method)), collapse = ", "))
    return(NULL)
  }
  # A fixed facet order — the population list the legend prints — rather than
  # alphabetical, so the panel reads T cells first and the same way every time.
  df <- dplyr::mutate(df, lineage = factor(lineage, levels = union(lineages, unique(lineage))))

  colour_var <- NULL; colour_name <- NULL
  if (!is.null(groups) && nrow(groups) && "patient_id" %in% names(df)) {
    df <- df |>
      dplyr::mutate(.pid = slide_key(patient_id)) |>
      dplyr::left_join(dplyr::transmute(groups, .pid = slide_key(patient_id), group),
                       by = ".pid") |>
      dplyr::select(-.pid)
    if (any(!is.na(df$group))) {
      df$group    <- hotcold_order(df$group)
      colour_var  <- "group"
      colour_name <- "Hotscore group"
    }
  }
  if (is.null(colour_var) && "immuno_phe" %in% names(df) && any(!is.na(df$immuno_phe))) {
    df$immuno_phe <- hotcold_order(df$immuno_phe)
    colour_var    <- "immuno_phe"
    colour_name   <- "Immuno-phenotype (clinical)"
  }

  p <- ggplot(df, aes(ihc_frac, score)) +
    (if (!is.null(colour_var))
       geom_point(aes(colour = .data[[colour_var]]), size = 2.4, alpha = .85)
     else geom_point(size = 2.4, alpha = .85, colour = oi[1])) +
    facet_wrap(~ lineage, scales = "free", labeller = as_labeller(.lineage_labels)) +
    # Every panel has its own x and y range, because imaging counts cells while
    # deconvolution estimates a mixture fraction. theme_paper_panels() gives the
    # border and the wider gutter that say the axes are not shared — without them
    # four free-scaled panels sit edge to edge and invite a cross-panel comparison
    # this figure cannot support.
    theme_paper_panels() +
    # Four breaks, not ggplot's default five-to-seven: each facet has its own range,
    # and at print size the default overprints the x tick labels.
    scale_x_continuous(breaks = scales::breaks_extended(n = 4)) +
    scale_y_continuous(breaks = scales::breaks_extended(n = 4)) +
    labs(x = x_lab, y = y_lab %||% paste(method, "fraction"))
  if (!is.null(colour_var))
    p <- p + scale_colour_manual(values = hotcold_cols(levels(df[[colour_var]])),
                                 na.value = "grey70", name = colour_name)
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
