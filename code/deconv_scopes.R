# =============================================================================
# deconv_scopes.R  —  WHICH CELLS the IHC side of a deconvolution comparison is
# counted over, when the answer is "more than one set".
#
# The molecular pages compare a deconvolution score against the IHC fraction over
# one arm's whole slide. The quanTIseq page asks the same question of every cell
# set the study has: each arm's whole slide, and each arm's pathologist polygons.
# That is four scopes over two arms, and this file is their registry — a scope
# names its arm, so no page can draw one without saying which arm it read.
#
#   scope                 cells                    kept
#   --------------------  -----------------------  ------------------------------
#   massimo1_wholeslide   FlowPath_csv_all/        every phenotyped cell
#   massimo1_annotation   FlowPath_csv_all/        inside the annotation_all/ polygon
#   massimo2_wholeslide   csv/, pooled + de-duped  every phenotyped cell
#   massimo2_annotation   csv/, pooled + de-duped  inside any annotation/ region
#                                                  (dissolved)
#
# THE TWO WHOLE-SLIDE SCOPES ARE NOT THE SAME NUMBER. The arms gate the same slides
# at different thresholds, so "every phenotyped cell" has one answer per arm.
#
# massimo1 HAS ONE ANNOTATION SCOPE, NOT TWO. Its `annotation_all` polygon is the
# exact union of that patient's `annotation_selected` regions (symmetric difference
# 0 on all four patients that have both, checked 2026-10-07), so a second scope
# would redraw the first. `annotation_all` is the one read because 10338 and 15897
# have a polygon there and no `_selected` region at all.
#
# A PATIENT WITH NO POLYGON IN A SCOPE'S SET IS CUT BY ITS EXPORT'S OWN FLAG.
# FlowPath writes an inside/outside call (Out_of_annotation) into every csv, so a
# slide the pathologist did not draw in that set still has one: 24086 in massimo2's
# `annotation`, and the csv-only patients in both arms. `.scope_source` records
# which rule decided, per cell.
#
# THIS DIFFERS FROM THE CLINICAL PAGES ON PURPOSE. There, 24086 in massimo2 is
# counted WHOLE under the arm's `bare_region_is` convention without its flag being
# read. Here a missing polygon always means "use the flag" —
# arm_cells_in_annotation(unannotated = "flag") — so an annotation-restricted
# fraction on this page need not equal that patient's union metrics row there.
#
# massimo1_inverted has no scope: it re-classifies massimo1's regions and borrows
# both polygon trees, so it is a third phenotyping, not a third annotation.
#
# Depends on validation_helpers.R (ihc_comparison_fraction, paired_spearman, the
# deconvolution mapping) and arm_cells.R (arm_cells_in_annotation). No export
# column is named here.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# `tier` is the polygon tier arm_cells_in_annotation() scores against; NA means no
# polygon is consulted at all. Row order is reading order: per arm, widest first.
DECONV_SCOPES <- tibble::tribble(
  ~scope,                ~arm,       ~tier,    ~label,                    ~kept,
  "massimo1_wholeslide", "massimo1", NA,       "Massimo1 — whole slide",  "every phenotyped cell",
  "massimo1_annotation", "massimo1", "union",  "Massimo1 — annotation",   "inside the annotation_all polygon",
  "massimo2_wholeslide", "massimo2", NA,       "Massimo2 — whole slide",  "every phenotyped cell",
  "massimo2_annotation", "massimo2", "region", "Massimo2 — annotation",   "inside any annotation region"
)

deconv_scope <- function(scope) {
  i <- match(scope, DECONV_SCOPES$scope)
  if (length(scope) != 1 || is.na(i))
    stop("unknown deconvolution scope '", paste(scope, collapse = "', '"), "' — have: ",
         paste(DECONV_SCOPES$scope, collapse = ", "))
  as.list(DECONV_SCOPES[i, ])
}

# `cells` (the scope's ARM's cohort cells, one row per physical cell) with two
# added columns:
#   in_scope       logical; NA where neither a polygon nor a flag could decide
#   .scope_source  "all_cells" (a whole-slide scope), "sf" (point-in-polygon) or
#                  "flag" (the export's Out_of_annotation column)
deconv_scope_cells <- function(scope, cells, spec = NULL, um_per_px = 0.325) {
  sc <- deconv_scope(scope)
  if (is.null(cells) || nrow(cells) == 0) return(cells)
  if (is.na(sc$tier))
    return(dplyr::mutate(cells, in_scope = TRUE, .scope_source = "all_cells"))
  if (is.null(spec)) spec <- arm_spec(sc$arm)
  arm_cells_in_annotation(spec, cells, um_per_px = um_per_px,
                          tier = sc$tier, unannotated = "flag") |>
    dplyr::rename(in_scope = in_annotation, .scope_source = .in_annotation_source)
}

# Provenance, one row per patient: how many cells the arm holds, how many the scope
# keeps, and which rule decided. A scope that kept 100 % under a polygon and one
# that kept 60 % under a flag draw the same-looking panel.
deconv_scope_inventory <- function(scoped) {
  if (is.null(scoped) || nrow(scoped) == 0 || !"in_scope" %in% names(scoped))
    return(tibble::tibble())
  scoped |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(
      membership   = paste(sort(unique(stats::na.omit(.scope_source))), collapse = "/"),
      n_cells      = dplyr::n(),
      n_in_scope   = sum(in_scope %in% TRUE),
      pct_in_scope = round(100 * n_in_scope / n_cells, 1),
      .groups = "drop") |>
    dplyr::arrange(patient_id)
}

# The IHC side for one scope: ihc_comparison_fraction() over the cells the scope
# keeps, so every population is counted and divided exactly as on the molecular
# pages and only the cell set changes. A patient with no cell in scope has no row.
deconv_scope_ihc <- function(scope, cells, spec = NULL, um_per_px = 0.325) {
  scoped <- deconv_scope_cells(scope, cells, spec = spec, um_per_px = um_per_px)
  if (is.null(scoped) || nrow(scoped) == 0)
    return(list(fractions = tibble::tibble(), inventory = tibble::tibble()))
  kept <- scoped[scoped$in_scope %in% TRUE, , drop = FALSE]
  list(
    fractions = if (nrow(kept)) dplyr::mutate(ihc_comparison_fraction(kept),
                                              scope = scope, .before = 1)
                else tibble::tibble(),
    inventory = dplyr::mutate(deconv_scope_inventory(scoped), scope = scope, .before = 1))
}

# Spearman per (scope, method, population) over a paired frame carrying `score` and
# `ihc_frac`; NA under 3 patients or a constant vector (paired_spearman).
deconv_scope_concordance <- function(paired) {
  paired |>
    dplyr::group_by(scope, method, lineage) |>
    dplyr::group_modify(~ paired_spearman(.x$score, .x$ihc_frac)) |>
    dplyr::ungroup()
}

# One population in one scope: IHC fraction (x) against the method's score (y), one
# point per patient, coloured by the clinical immuno-phenotype when `d` carries it.
# The dashed line is x = y, drawn only for a method whose score is a fraction of the
# same denominator as x — otherwise the two axes share no unit.
plot_deconv_scope_pair <- function(d, scope, lineage, method = "quantiseq") {
  sc   <- deconv_scope(scope)
  row  <- deconv_comparison_lineages[match(lineage, deconv_comparison_lineages$lineage), ]
  st   <- paired_spearman(d$score, d$ihc_frac)
  rho  <- if (is.finite(st$rho)) formatC(st$rho, format = "f", digits = 2) else "NA"
  cd45 <- identical(deconv_ihc_denominator(method), "CD45")
  if (!"immuno_phe" %in% names(d)) d$immuno_phe <- NA_character_
  d$immuno_phe <- hotcold_order(d$immuno_phe)

  ggplot(d, aes(ihc_frac, score)) +
    { if (deconv_is_additive(method))
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = REF_LINE) } +
    geom_smooth(method = "lm", se = FALSE, colour = FIT_LINE, formula = y ~ x) +
    geom_point(aes(colour = immuno_phe), alpha = 0.75, size = 1.8) +
    scale_colour_manual(values = hotcold_cols(levels(d$immuno_phe)),
                        na.value = "grey70", name = "Immuno-phenotype") +
    expand_limits(x = 0, y = 0) +
    labs(title = sprintf("%s: %s (%s)", method_label(method), lineage, row$ihc_gate),
         subtitle = with_n(sprintf("%s; Spearman rho = %s", sc$label, rho),
                           d$patient_id, "patients"),
         x = sprintf("IHC fraction of %s cells in scope (unitless, 0-1)",
                     if (cd45) "CD45+" else "all"),
         y = paste(method_label(method), "score (unitless)"))
}
