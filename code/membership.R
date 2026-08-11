# =============================================================================
# membership.R  —  which cells count as "inside the tumour annotation"
#
# ONE question, three answers, one interface. Every clinical_data report picks a
# mode and gets the same three frames back, so the report body never branches on
# how membership was decided:
#
#   mode        inside a cell is ...                       cells come from
#   ----------  ----------------------------------------   -------------------------
#   "all_slide" its centroid falls in THAT REGION'S          data/all_slide/csv/<pid>/
#               pathologist polygon. The export and the        <pid>_<A|B|C>.csv
#               polygon are matched per region, so no      (+ .../annotation/<pid>/
#               cell is scored against a region it was       <pid>_<A|B|C>.geojson)
#               not exported for. A patient with NO
#               annotation directory is entirely inside
#               by the layout's own convention.
#   "mirage"    the same geometric rule, applied to        data/mirage/<patient>/
#               cells MIRAGE phenotyped itself rather        (+ the same polygons)
#               than FlowPath
#
# The diff between these two isolates the PHENOTYPING METHOD: same slides, same
# polygons, different tool deciding what each cell is.
#
# THREE MODES WERE REMOVED (2026-08-11) because their layouts no longer exist:
#   "geojson"   wanted a whole-slide data/flowpath/<patient>.csv plus a FLAT
#               data/annotation/<patient>_a<k>.geojson. The export is now per-region
#               and nested, so neither path is produced. "all_slide" is its successor
#               and strictly better: it scores each region against its own polygon
#               instead of one slide against all of them.
#   "flag"      wanted data/flowpath/per_annotation/<patient>_a<k>.csv and took
#               membership from the exporter's Out_of_annotation column rather than
#               from geometry.
#   "flag_old"  the same, with a per_annotation/old/ overlay on top. There is no
#               old/new split in the export any more.
# The flag modes are also why has_area used to be FALSE for some reports: only a
# polygon knows a region's AREA, so those modes returned NA for every density. Both
# surviving modes read polygons, so densities are always available now.
#
# mirage has no out-of-annotation flag of its own and never had a fallback: a patient
# with cells but no polygon contributes nothing there rather than being counted whole.
# "all_slide" DOES count such a patient whole — but only because its layout states
# that a missing annotation directory means "everything is inside", which is a fact
# about the data rather than a guess (see code/all_slide.R).
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
source(here::here("code", "all_slide.R"))      # the "all_slide" mode's cell + polygon source

# =============================================================================
# The one entry point the reports use
# =============================================================================
# `membership_data()` picks a mode and returns the SAME shape either way, so
# clinical_data's shared body never asks how membership was decided. Both modes load
# their own cells, so there is no plumbing above this point any more — the flag modes
# that needed it are gone.
#
#   $mode        the mode string, echoed back
#   $cells       the per-patient cell set for the whole-cohort immune readout, one
#                row per physical cell (de-duplicated across region files)
#   $per_annotation  metrics, one row per (patient_id, annotation)
#   $union           metrics, one row per patient
#   $inventory   provenance table — what was actually read, per patient and region
#   $has_area    TRUE when a polygon was read, i.e. when the dens_*_per_mm2 columns
#                carry values rather than NA
#
# `annots` is required for "mirage", whose cells come from a different tree than its
# polygons. "all_slide" reads its own polygons, matched to its own cell files region
# by region, so an externally supplied set would not line up — passing one is an
# error rather than a silently ignored argument. `neoplastic_data` is accepted and
# unused; it stays in the signature so the call sites need not change and so a future
# mode can label a fallback patient from its single scored annotation.
MEMBERSHIP_MODES <- c("all_slide", "mirage")

membership_data <- function(mode = MEMBERSHIP_MODES, ihc_data,
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
      has_area       = TRUE
    ))
  }

  # -- the region-matched polygon mode ----------------------------------------
  # Cells and polygons come from the SAME tree and are paired by region, so this
  # branch loads both itself. has_area is per-patient here, not per-cohort: a
  # whole-slide patient has no polygon and therefore no area, while its annotated
  # neighbours do, so the flag reports what at least one patient can support and
  # the metrics frame's own `source` column says which patient is which.
  if (mode == "all_slide") {
    if (!is.null(annots))
      stop("membership_data(\"all_slide\") reads its own polygons; drop the `annots` argument")
    cells <- all_slide_cells(ALL_SLIDE_DIR)
    if (nrow(cells) == 0) {
      warning("membership_data(\"all_slide\"): no cells under ", ALL_SLIDE_DIR,
              " — every panel will be empty")
      return(list(mode = mode, cells = cells, per_annotation = tibble::tibble(),
                  union = tibble::tibble(), inventory = tibble::tibble(),
                  has_area = FALSE))
    }
    polys <- all_slide_annotations(ALL_SLIDE_DIR, patient_ids = unique(cells$patient_id))
    return(list(
      mode           = mode,
      cells          = all_slide_union_cells(cells),
      per_annotation = all_slide_metrics(cells, polys, "per_annotation", um_per_px),
      union          = all_slide_metrics(cells, polys, "union", um_per_px),
      inventory      = all_slide_inventory(cells),
      has_area       = !is.null(polys) && nrow(polys) > 0
    ))
  }
}

