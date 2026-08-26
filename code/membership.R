# =============================================================================
# membership.R  —  which cells count as "inside the tumour annotation"
#
# ONE question, four answers, one interface. Every clinical report picks a mode
# and gets the same frames back, so the report body never branches on how
# membership was decided:
#
#   mode                 cells come from                    polygons
#   -------------------  ---------------------------------  ------------------------
#   "massimo1"           FlowPath_csv_selected/<pid>_a<k>    annotation_selected/ per
#                        .csv, plus the whole-slide          region; annotation_all/
#                        FlowPath_csv_all/<pid>/<pid>.csv    for the union
#   "massimo1_inverted"  csv_inverted-classification_        massimo1's polygons —
#                        modified-thrPANCK/<pid>_a<k>.csv    it re-classifies arm 1's
#                                                            regions, does not redraw
#   "massimo2"           csv/<pid>/<pid>_<A|B|C>.csv         annotation/<pid>/ per
#                                                            region; the union is
#                                                            those, dissolved
#   "mirage"             data/mirage/<patient>/              whatever `annots` the
#                        (mirage's own phenotyping)          caller supplies
#
# THE THREE MASSIMO ARMS ARE THE SAME SLIDES, THREE PHENOTYPINGS. Their diff
# isolates the METHOD: massimo1 vs massimo2 is two FlowPath exports with
# independently drawn regions, massimo1 vs massimo1_inverted is one region set
# classified two ways at different PANCK thresholds. "mirage" is the fourth
# phenotyper, on its own tree.
#
# REGIONS ARE ARM-LOCAL. massimo1's ANNOTATION_2 and massimo2's ANNOTATION_2 are
# different polygons over the same tissue — the pathologist drew each arm's set in
# a separate session. So each arm carries its own pathologist table
# (neoplastic_massimo1 / neoplastic_massimo2), and a cross-arm figure joins on
# PATIENT and never on (patient, annotation). code/arms.R states this in full.
#
# THREE MODES WERE REMOVED (2026-08-11) because their layouts no longer exist:
#   "geojson"   wanted a whole-slide data/flowpath/<patient>.csv plus a FLAT
#               data/annotation/<patient>_a<k>.geojson.
#   "flag"      wanted data/flowpath/per_annotation/<patient>_a<k>.csv and took
#               membership from the exporter's Out_of_annotation column rather
#               than from geometry.
#   "flag_old"  the same, with a per_annotation/old/ overlay on top.
# A FOURTH, "all_slide", was RENAMED to "massimo2" (2026-08-26) when the second
# and third arms arrived and one export stopped being "the" export. Same tree,
# same rules; only the name says which arm it is.
#
# mirage has no out-of-annotation flag of its own and never had a fallback: a
# patient with cells but no polygon contributes nothing there rather than being
# counted whole. massimo2 DOES count such a patient whole — but only because its
# layout states that a missing annotation directory means "everything is inside",
# which is a fact about the data rather than a guess (see code/arm_cells.R).
#
# Depends on validation_helpers.R (slide_key, region_ratios_area,
# ihc_annotation_metrics) and, through it, on cell_tables.R for every column read
# off an export — no export-specific column name appears in this file.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
  library(readr)
  library(stringr)
  library(purrr)
})

source(here::here("code", "mirage_cells.R"))   # the "mirage" mode's cell source
source(here::here("code", "arm_cells.R"))      # the three Massimo arms (pulls in arms.R)

# =============================================================================
# The one entry point the reports use
# =============================================================================
# `membership_data()` picks a mode and returns the SAME shape whichever it is, so a
# report body never asks how membership was decided. Every mode loads its own
# cells, so there is no plumbing above this point.
#
#   $mode            the mode string, echoed back
#   $cells           the per-patient cell set for the whole-cohort immune readout,
#                    one row per PHYSICAL cell
#   $per_annotation  metrics, one row per (patient_id, annotation)
#   $union           metrics, one row per patient
#   $inventory       provenance table — what was actually read, per patient/region
#   $has_area        TRUE when a polygon was read, i.e. when the dens_*_per_mm2
#                    columns carry values rather than NA
#   $overlap         arm_overlap_report() — whether the region files repeat cells
#   $reconciliation  arm_wholeslide_reconciliation() — only for an arm that ships a
#                    whole-slide tier, i.e. massimo1. Empty otherwise.
#
# `annots` is required for "mirage", whose cells come from a different tree than its
# polygons. The Massimo arms read their own polygons, matched to their own cell files
# region by region, so an externally supplied set would not line up — passing one is
# an error rather than a silently ignored argument. `neoplastic_data` is accepted and
# unused; it stays in the signature so the call sites need not change.
MEMBERSHIP_MODES <- c(ARM_MODES, "mirage")

membership_data <- function(mode = MEMBERSHIP_MODES, ihc_data = NULL,
                            neoplastic_data = NULL, annots = NULL,
                            um_per_px = 0.325) {
  mode <- match.arg(mode)

  # -- mirage: cells from the pipeline outdir, polygons from the annotation tree ---
  if (mode == "mirage") {
    if (is.null(annots)) stop("membership_data(\"mirage\") needs `annots`")
    cells <- load_mirage_cells()
    if (nrow(cells) == 0)
      warning("membership_data(\"", mode, "\"): no cells loaded — every panel will be empty")
    return(list(
      mode           = mode,
      cells          = cells,
      per_annotation = ihc_annotation_metrics(cells, annots, scope = "per_annotation",
                                              um_per_px = um_per_px),
      union          = ihc_annotation_metrics(cells, annots, scope = "union",
                                              um_per_px = um_per_px),
      inventory      = mirage_cells_inventory(cells),
      has_area       = TRUE,
      overlap        = tibble::tibble(),
      reconciliation = tibble::tibble()
    ))
  }

  # -- the region-matched Massimo arms ----------------------------------------
  # Cells and polygons come from the SAME arm and are paired by region, so this
  # branch loads both itself. has_area is per-patient, not per-cohort: a
  # whole-slide patient has no polygon and therefore no area, while its annotated
  # neighbours do, so the flag reports what at least one patient can support and
  # the metrics frame's own `source` column says which patient is which.
  if (!is.null(annots))
    stop("membership_data(\"", mode, "\") reads its own polygons; drop the `annots` argument")

  spec  <- arm_spec(mode)
  cells <- arm_cells(spec)
  ucell <- arm_union_tier_cells(spec)

  if (nrow(cells) == 0 && nrow(ucell) == 0) {
    warning("membership_data(\"", mode, "\"): no cells under ", spec$root_path,
            " — every panel will be empty")
    # TYPED empty metrics, not bare ones: an arm whose tree is not on disk yet is a
    # routine state, and its frames have to flow through the same selects and joins
    # a full arm's do. A zero-COLUMN tibble instead errors with "Column
    # `patient_id` doesn't exist" several chunks downstream, naming neither the arm
    # nor the directory that was missing.
    return(list(mode = mode, cells = cells,
                per_annotation = arm_empty_metrics(),
                union          = arm_empty_metrics(),
                inventory      = tibble::tibble(),
                has_area       = FALSE,
                overlap        = tibble::tibble(),
                reconciliation = tibble::tibble()))
  }

  # `$patient_id` on a 0-column tibble warns rather than returning NULL quietly, and
  # `ucell` IS 0-column for the two arms with no whole-slide tier — so a bare
  # c(cells$patient_id, ucell$patient_id) printed two spurious "unknown column"
  # warnings into every knit of a massimo2 or massimo1_inverted page. A page that
  # cries wolf on a normal load is a page whose real warnings stop being read.
  .pids  <- function(x) if (nrow(x) && "patient_id" %in% names(x)) x$patient_id else character(0)
  pids   <- unique(c(.pids(cells), .pids(ucell)))
  rpolys <- arm_annotations(spec, "region", patient_ids = pids)
  upolys <- arm_annotations(spec, "union",  patient_ids = pids)

  # A patient with a whole-slide export but no region files keeps a per-region row
  # by promoting its union polygon to ANNOTATION_1 — massimo1's 10338 and 15897.
  promoted <- arm_promote_unregioned(spec, cells, rpolys, ucell, upolys)
  cells    <- promoted$cells
  rpolys   <- promoted$polys

  list(
    mode           = mode,
    cells          = arm_cohort_cells(spec, cells, ucell),
    per_annotation = arm_metrics(spec, cells, rpolys, "per_annotation", um_per_px),
    union          = arm_metrics(spec, cells, rpolys, "union", um_per_px,
                                 union_cells = ucell, union_polys = upolys),
    inventory      = arm_inventory(cells),
    has_area       = (!is.null(rpolys) && nrow(rpolys) > 0) ||
                     (!is.null(upolys) && nrow(upolys) > 0),
    overlap        = arm_overlap_report(cells),
    reconciliation = arm_wholeslide_reconciliation(spec, cells, ucell)
  )
}
