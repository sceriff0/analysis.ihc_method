# =============================================================================
# scope_compare.R  —  ONE quantity, THREE nested scopes, and how to compare a
# single number against several.
#
# The tumour fraction of a slide has three defensible answers, and they are not
# interchangeable:
#
#   whole_slide     every cell in the export. NO polygon is consulted at all, so
#                   this is the only one that exists for a slide nobody annotated
#                   and the only one that cannot be blamed on where a line was drawn.
#   annotation_all  the cells inside the DISSOLVED union polygon — massimo1's
#                   `annotation_all/`, or the dissolved per-region polygons for an
#                   arm without that tier. One value per patient, like whole_slide,
#                   but restricted to tissue the pathologist bounded.
#   annotation_k    one pathologist region. SEVERAL values per patient, and the
#                   only scope the per-region pathologist scores can be joined to.
#
# The three nest: annotation_k ⊆ annotation_all ⊆ whole_slide. So the difference
# between them is not noise — it is exactly the effect of the annotation, and the
# point of putting them on one figure is to size that effect rather than assume it
# is small.
#
# COMPARING ONE VALUE TO MANY is the actual difficulty, and there is no single
# right answer, so this file offers three and shows them together:
#
#   (a) DIRECTLY — whole_slide against annotation_all, one point per patient. Both
#       are single values, so this is an ordinary scatter against x = y. It answers
#       "does restricting to the annotation move the number, and which way?"
#   (b) VISUALLY — per patient, the spread of its ANNOTATION_k values drawn as a
#       range with the two single values marked on it. It answers "is the single
#       value inside the range its own regions span, or outside it?", which a
#       correlation cannot show and which is the question a reader actually has.
#   (c) QUANTITATIVELY — whole_slide against each AGGREGATE of that patient's
#       regions (cell-weighted mean, plain mean, median), reusing the aggregators
#       in aggregation_compare.R. It answers "which way of collapsing the regions
#       reproduces the unannotated number best?" and reports bias, not just
#       correlation: a method can rank patients perfectly and still read high.
#
# Cell-weighting is the aggregator to read first. A patient's regions differ in
# size by an order of magnitude, so a plain mean lets a 2,000-cell region and a
# 200,000-cell one vote equally — which is a statement about the pathologist's
# drawing habits, not about the tissue.
#
# Depends on validation_helpers.R (region_ratios_area, paired_cor3, slide_key) and
# aggregation_compare.R (IHC_AGG_LABELS, .wmean). No export column is named here:
# every one goes through cell_tables.R, as everywhere else.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

# --- The whole-slide scope ---------------------------------------------------
# One row per patient over ALL of that patient's cells, no polygon. Schema-identical
# to the metrics frames so it stacks with them without a rename.
#
# `cells` must be the COHORT set — one row per physical cell per patient
# (membership_data()$cells). Handing it the raw region tier would count every cell
# once per region file it appears in, which is the inflation arm_cohort_cells()
# exists to prevent, and the resulting fraction would still look plausible.
arm_wholeslide_ratios <- function(cells, um_per_px = 0.325) {
  if (is.null(cells) || nrow(cells) == 0) return(tibble::tibble())
  purrr::map_dfr(sort(unique(cells$patient_id)), function(pid) {
    cp <- dplyr::filter(cells, patient_id == pid)
    # area 0 -> NA densities. There is no polygon, so there is no area: reporting a
    # density here would require inventing a denominator.
    region_ratios_area(cp, 0, um_per_px) |>
      dplyr::mutate(patient_id = pid, annotation = "whole_slide",
                    source = "whole_slide", .before = 1)
  })
}

# --- The three scopes as one long frame --------------------------------------
# One row per (patient_id, scope, annotation). `scope` is the named-palette level;
# `annotation` keeps the region label so a per-region point can still be identified.
#
# An arm whose union rows are labelled "whole_slide" (massimo2's 24086, which has no
# annotation directory) is NOT relabelled: that patient genuinely has no annotated
# scope, so it contributes to whole_slide alone and is absent from the other two.
# Forcing it into annotation_all would claim a polygon nobody drew.
scope_table <- function(mem, metric = "tumor_over_inside", um_per_px = 0.325) {
  ws <- arm_wholeslide_ratios(mem$cells, um_per_px)

  pick <- function(d, scope) {
    if (is.null(d) || nrow(d) == 0) return(tibble::tibble())
    keep <- intersect(c("patient_id", "annotation", "source", "n_inside",
                        "n_tumor_inside", "area_mm2", metric), names(d))
    d |>
      dplyr::select(dplyr::all_of(keep)) |>
      dplyr::mutate(scope = scope, .after = patient_id)
  }

  un  <- if (!is.null(mem$union) && nrow(mem$union))
    dplyr::filter(mem$union, annotation == "union") else tibble::tibble()
  per <- if (!is.null(mem$per_annotation) && nrow(mem$per_annotation))
    dplyr::filter(mem$per_annotation, grepl("^ANNOTATION_[0-9]+$", annotation))
  else tibble::tibble()

  out <- dplyr::bind_rows(pick(ws, "whole_slide"),
                          pick(un, "annotation_all"),
                          pick(per, "annotation_k"))
  if (nrow(out) == 0) return(out)
  out |>
    dplyr::mutate(scope = scope_factor(scope),
                  value = .data[[metric]]) |>
    dplyr::arrange(patient_id, scope, annotation)
}

# --- (b) one value vs many: the per-patient spread ---------------------------
# Per patient: the two single values, and the min/median/max of its ANNOTATION_k
# values. `inside_range` is the readout — whether the unannotated number falls
# within the span its own regions cover.
#
# A patient with ONE region has a degenerate range (min == max), which is not a
# failure: 10338 and 15897 have a single annotation each, and the comparison is
# still meaningful — it is simply exact rather than a span.
scope_one_vs_many <- function(scope_tbl) {
  if (nrow(scope_tbl) == 0) return(tibble::tibble())
  single <- scope_tbl |>
    dplyr::filter(scope %in% c("whole_slide", "annotation_all")) |>
    dplyr::select(patient_id, scope, value) |>
    tidyr::pivot_wider(names_from = scope, values_from = value)

  many <- scope_tbl |>
    dplyr::filter(scope == "annotation_k", is.finite(value)) |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(n_ann      = dplyr::n(),
                     ann_min    = min(value),
                     ann_median = stats::median(value),
                     ann_max    = max(value),
                     n_cells_ann = sum(n_inside, na.rm = TRUE),
                     .groups = "drop")

  dplyr::full_join(single, many, by = "patient_id") |>
    dplyr::mutate(
      inside_range = dplyr::if_else(
        is.finite(whole_slide) & is.finite(ann_min) & is.finite(ann_max),
        whole_slide >= ann_min & whole_slide <= ann_max, NA),
      gap_to_union = whole_slide - annotation_all) |>
    dplyr::arrange(patient_id)
}

# --- (c) one value vs many: which aggregate reproduces it? -------------------
# Reuses aggregate_ihc()'s aggregators so a reader meets ONE vocabulary of ways to
# collapse regions across the whole site, not a second one invented here. Returns
# one row per (aggregator, patient) plus the agreement stats per aggregator.
scope_aggregator_pairs <- function(scope_tbl, mem, metric = "tumor_over_inside") {
  per <- if (!is.null(mem$per_annotation) && nrow(mem$per_annotation))
    dplyr::filter(mem$per_annotation, grepl("^ANNOTATION_[0-9]+$", annotation))
  else tibble::tibble()
  if (nrow(per) == 0) return(tibble::tibble())

  ws <- scope_tbl |>
    dplyr::filter(scope == "whole_slide") |>
    dplyr::select(patient_id, whole_slide = value)

  # aggregate_ihc() wants a `union` frame for its `pooled` aggregator; the union
  # scope IS that, so annotation_all enters the comparison as one aggregator among
  # the rest rather than as a separate figure.
  un <- if (!is.null(mem$union) && nrow(mem$union))
    dplyr::filter(mem$union, annotation == "union") else tibble::tibble()

  aggregate_ihc(per, un, metric) |>
    dplyr::inner_join(ws, by = "patient_id") |>
    dplyr::filter(is.finite(ihc_val), is.finite(whole_slide)) |>
    dplyr::mutate(ihc_lab = factor(IHC_AGG_LABELS[ihc_agg], levels = IHC_AGG_LABELS))
}

scope_aggregator_stats <- function(pairs) {
  if (nrow(pairs) == 0) return(tibble::tibble())
  pairs |>
    dplyr::group_by(ihc_agg, ihc_lab) |>
    dplyr::group_modify(~ {
      cc <- paired_cor3(.x$ihc_val, .x$whole_slide)
      # bias = mean(aggregate - whole slide): the systematic offset the correlation
      # coefficients are blind to. An aggregator can rank patients perfectly and
      # still read ten points high on every one of them.
      dplyr::mutate(cc, bias = mean(.x$ihc_val - .x$whole_slide))
    }) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(spearman))
}

# --- Figures -----------------------------------------------------------------
# (a) The two single values against each other. Both axes are the same fraction of
# the same quantity, so x = y is the target and vertical distance from it is the
# effect of restricting to the annotation.
plot_scope_scatter <- function(one_many, title = NULL) {
  d <- dplyr::filter(one_many, is.finite(whole_slide), is.finite(annotation_all))
  if (nrow(d) == 0) return(invisible())
  ggplot(d, aes(whole_slide, annotation_all)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = REF_LINE, linewidth = pt_line(0.5)) +
    geom_point(size = 2.2, colour = SCOPE_COLS[["annotation_all"]]) +
    # Plain geom_text, nudged, as paper_figures.R does. ggrepel would be neater but
    # is not in renv.lock, and that file is not hand-editable.
    geom_text(aes(label = patient_id), size = pt_text(6),
              vjust = -0.9, colour = "grey35") +
    coord_equal() +
    labs(x = "tumour fraction, whole slide (no annotation)",
         y = "tumour fraction, annotation_all",
         title = title,
         subtitle = with_n(NULL, d$patient_id, "patients")) 
}

# (b) The spread of a patient's regions, with the two single values marked on it.
# Patients on the y axis so the ranges are horizontal and directly comparable.
plot_scope_range <- function(one_many, scope_tbl, title = NULL) {
  d <- dplyr::filter(one_many, is.finite(ann_min) | is.finite(whole_slide))
  if (nrow(d) == 0) return(invisible())
  ord <- d$patient_id[order(dplyr::coalesce(d$whole_slide, d$ann_median))]
  d$patient_id <- factor(d$patient_id, levels = ord)
  pts <- scope_tbl |>
    dplyr::filter(scope == "annotation_k", is.finite(value)) |>
    dplyr::mutate(patient_id = factor(patient_id, levels = ord))

  ggplot(d, aes(y = patient_id)) +
    geom_segment(aes(x = ann_min, xend = ann_max, yend = patient_id),
                 colour = "grey75", linewidth = pt_line(2)) +
    # LAYER ORDER MATTERS for the single-annotation patients. Where a patient has one
    # region, annotation_all IS that region by construction, so the two markers sit
    # at the same x. Drawing the small region circle LAST puts it on top of the
    # larger triangle instead of under it — otherwise 10338 and 15897 look as though
    # they have no region at all, which is the opposite of what the figure says.
    geom_point(aes(x = annotation_all, colour = "annotation_all"),
               size = 2.8, shape = 17) +
    geom_point(aes(x = whole_slide, colour = "whole_slide"),
               size = 2.8, shape = 18) +
    geom_point(data = pts, aes(x = value, colour = "annotation_k"),
               size = 1.5, stroke = 0.7, shape = 21, fill = "white") +
    scale_colour_scope(name = NULL) +
    scale_y_discrete(labels = label_n(pts$patient_id)) +
    labs(x = "tumour fraction", y = NULL, title = title,
         subtitle = paste0("grey bar spans each patient's ANNOTATION_k values; ",
                           with_n(NULL, d$patient_id, "patients"))) +
    theme(legend.position = "top")
}

# (c) Whole slide against each aggregate of the regions, one facet per aggregator.
plot_scope_aggregators <- function(pairs, title = NULL) {
  if (nrow(pairs) == 0) return(invisible())
  ggplot(pairs, aes(whole_slide, ihc_val)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = REF_LINE, linewidth = pt_line(0.5)) +
    geom_point(size = 1.9, colour = SCOPE_COLS[["annotation_k"]]) +
    # 2 x 2, not 1 x 4. coord_equal() is not optional here — the dashed x = y line
    # only reads as "agreement" when it is actually at 45 degrees — but a single row
    # of four square panels leaves so much vertical slack on a double-column figure
    # that the title block collides with the y-axis label. A square grid fills the
    # canvas the aspect lock demands.
    facet_wrap(~ ihc_lab, nrow = 2) +
    coord_equal() +
    labs(x = "tumour fraction, whole slide (no annotation)",
         y = "tumour fraction, aggregated over annotations",
         title = title,
         subtitle = with_n(NULL, unique(pairs$patient_id), "patients")) +
    theme_paper_panels()
}
