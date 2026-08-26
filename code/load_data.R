# =============================================================================
# load_data.R  —  the raw inputs, loaded once
#
# Sourced first by every analysis. Defines, from data/ (gitignored):
#   dds                     DESeq2 object from counts.RData, with DESeq() already run
#   clinical_data           the clinical CRF (clinical_data.xlsx)
#   neoplastic_massimo1     the pathologist's tumour-content scores for ARM 1, inline
#   neoplastic_massimo2     the same for ARM 2. Both LONG: one row per
#                           (SAMPLE, annotation) — see the note at their definition
#   counts_data             normalised counts, wide: one row per bulk-RNA Sample
#   ihc_massimo1            single cells, one row per cell per patient, ARM 1
#   ihc_massimo2            the same for ARM 2
#   ihc_massimo1_inverted   the same for ARM 3
#
# NOTHING HERE IS NAMED `ihc_data` OR `neoplastic_data` ANY MORE, and that is the
# point. Three arms phenotype the same slides three ways; a bare name would let a
# page use one of them without saying which, and the resulting figure would be
# indistinguishable from a figure about a different arm. Every page names its arm.
#
# It no longer prints colnames() on load. The callers still wrap the source() in
# capture.output(suppressMessages(...)), which now serves the LOADED / whole-slide
# -check messages load_arm_cells() emits per arm rather than a column dump.
#
# Derived quantities live in validation_helpers.R, membership rules in
# membership.R, the arm registry in arms.R and the cell-table schema in
# cell_tables.R. Nothing here reduces or reshapes: this file only loads.
# =============================================================================
library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)
library(stringr)
library(forcats)
library(DESeq2)
library(here)
library(readxl)
library(fs)
library("data.table")

# The cell-table schema adapter (read_cell_csv() and the accessors every analysis
# reads an export through) plus the arm reader, which needs slide_key() from
# validation_helpers. Sourcing validation_helpers here is not a cycle: it never reads
# anything this file defines.
source(here("code", "validation_helpers.R"))   # pulls in cell_tables.R + plot_theme.R
source(here("code", "arm_cells.R"))            # pulls in arms.R

dds <- get(load(here("data", "counts.RData")))
dds <- DESeq(dds)

clinical_data <- read_excel(here("data", "clinical_data.xlsx")) |>
  filter(!is.na(`ID PATIENT`)) |>
  mutate(`ID CRF PRESERVE` = gsub("-", ".", `ID CRF PRESERVE`))

# =============================================================================
# The pathologist's neoplastic-cellularity score — ONE TABLE PER ARM
# =============================================================================
# One row per SCORED REGION, and one table per arm, because THE ARMS' REGIONS ARE
# INDEPENDENTLY DRAWN. massimo1's ANNOTATION_2 and massimo2's ANNOTATION_2 are
# different polygons over the same tissue — the pathologist annotated each arm in a
# separate session — so a single shared table would silently score arm 1's regions
# with arm 2's percentages. code/arms.R states the rule; this is where it bites.
# massimo1_inverted re-classifies arm 1's own regions rather than redrawing them, so
# it reads neoplastic_massimo1 and has no table of its own.
#
# LONG, NOT WIDE (~ANNOTATION_1..3), for one reason: in arm 2, 24086 has no
# annotation drawn at all and its 75% refers to the WHOLE SLIDE. A wide frame can
# only put that in ANNOTATION_1, which would claim a region the pathologist never
# drew — and would then fail to join, because the metrics frame labels that
# patient's single row `whole_slide`. Long says what was actually scored.
#
# `annotation` must match the label the membership metrics emit, because the
# clinical pages join on (patient_id, annotation):
#   ANNOTATION_<k>  region k. In arm 2, k is the alphabet position of the export's
#                   letter suffix (A -> 1, B -> 2, C -> 3). In arm 1, k is the digit
#                   in the `_a<k>` suffix, taken literally.
#   whole_slide     no annotation directory, so every cell counts (arm_cells.R)

# --- ARM 2: the letter-suffixed export ---------------------------------------
# The region counts cross-check against the export exactly — 046 three csvs and three
# geojsons, 052 two, 5456 three, 10338 one, 15897 two, and 24086 a bare csv with no
# annotation directory. A mismatch between the two would silently drop a region from
# the correlation, so the clinical page prints the reconciliation.
#
# Values updated 2026-08-11 from the pathologist's re-read. They differ materially
# from the previous set, so the tumour-content correlation is NOT comparable to an
# earlier knit: 046 A 50->30, 052 70/50->50/75, 5456 A 70->80, 10338 75->80,
# 15897 gains B=75 (was single-annotation), and 24086 goes from three annotations
# (80/70/70) to one whole-slide score of 75.
neoplastic_massimo2 <- tribble(
  ~SAMPLE,  ~annotation,     ~path_pct,
  "046",    "ANNOTATION_1",         30,
  "046",    "ANNOTATION_2",         60,
  "046",    "ANNOTATION_3",         50,
  "052",    "ANNOTATION_1",         50,
  "052",    "ANNOTATION_2",         75,
  "5456",   "ANNOTATION_1",         80,
  "5456",   "ANNOTATION_2",         80,
  "5456",   "ANNOTATION_3",         70,
  "10338",  "ANNOTATION_1",         80,
  "15897",  "ANNOTATION_1",         60,
  "15897",  "ANNOTATION_2",         75,
  "24086",  "whole_slide",          75
)

# --- ARM 1: the `_a<k>` selected regions -------------------------------------
# Values read 2026-08-26 from data/Massimo1/thr_head&neck.xlsx, which carries one
# sheet per patient: marker thresholds, then a `neoplastic cellularity (%)` block
# keyed annotation_1..3. Its 11 scored regions match arm 1's 11 `_selected` exports
# exactly.
#
# THESE ARE THE PATHOLOGIST'S ORIGINAL READ; arm 2's are the RE-READ. That is the
# clearest evidence the two arms are separate annotation sessions rather than one
# relabelled: 046 A 50 vs 30, 052 70/50 vs 50/75, 5456 A 70 vs 80, and 24086 three
# regions at 80/70/70 here against a single whole-slide 75 there. Neither supersedes
# the other — each scores the polygons its own arm exported, which is why they are
# two tables and never one.
#
# The 13 rows are the regions arm 1 actually exports, and they are NOT arm 2's:
#   046, 5456        a1..a3        three selected regions each
#   052              a1, a2        two
#   24086            a1..a3        THREE — arm 2 has none for this patient at all,
#                                  where it is a bare whole-slide csv scored 75%.
#                                  Both are correct; the arms were annotated apart.
#   10338, 15897     ANNOTATION_1  no `_selected` files at all. Their single
#                                  `annotation_all` polygon is promoted to
#                                  ANNOTATION_1 by arm_promote_unregioned(), so the
#                                  score keyed here is the score for that polygon —
#                                  the whole slide as the pathologist bounded it,
#                                  not an unbounded slide. It is NOT `whole_slide`:
#                                  that label means "no polygon exists", and one does.
neoplastic_massimo1 <- tribble(
  ~SAMPLE,  ~annotation,     ~path_pct,
  "046",    "ANNOTATION_1",         50,
  "046",    "ANNOTATION_2",         60,
  "046",    "ANNOTATION_3",         50,
  "052",    "ANNOTATION_1",         70,
  "052",    "ANNOTATION_2",         50,
  "5456",   "ANNOTATION_1",         70,
  "5456",   "ANNOTATION_2",         80,
  "5456",   "ANNOTATION_3",         70,
  "24086",  "ANNOTATION_1",         80,
  "24086",  "ANNOTATION_2",         70,
  "24086",  "ANNOTATION_3",         70,
  # 10338 and 15897 have no `_selected` regions in arm 1 — one `annotation_all`
  # polygon each, promoted to ANNOTATION_1 by arm_promote_unregioned(). Their sheets
  # key the score as bare `annotation` (SINGULAR), not `annotation_<k>`, which is
  # exactly right for a patient with one region and is why a first pass reading only
  # `annotation_<k>` reported them unscored.
  #
  # 10338's sheet records a RANGE, "75-80"; 77.5 is its midpoint. It is the only
  # interpolated value in either table — every other number is transcribed. Arm 2
  # scored the same patient 80 against its own polygon.
  "10338",  "ANNOTATION_1",       77.5,
  "15897",  "ANNOTATION_1",         70
)

# The arm -> pathologist-table lookup, so a page says `neoplastic_for(ARM)` rather
# than branching. massimo1_inverted deliberately shares arm 1's scores: it is the
# same tissue, the same polygons, and only the CLASSIFICATION of the cells differs —
# and a pathologist's percentage is a property of the tissue, not of the classifier.
neoplastic_for <- function(arm = ARM_MODES) {
  arm <- match.arg(arm)
  switch(arm,
         massimo1          = neoplastic_massimo1,
         massimo1_inverted = neoplastic_massimo1,
         massimo2          = neoplastic_massimo2)
}

counts_data <-  counts(dds, normalized = TRUE) |>
  as_tibble(rownames = "GENE") |>
  pivot_longer(cols = -GENE, names_to = "Sample", values_to = "Expression") |>
  pivot_wider(names_from = GENE, values_from = Expression) |>
  filter(Sample %in% clinical_data$`ID CRF PRESERVE`)

# =============================================================================
# The whole-cohort cell set — ONE PER ARM
# =============================================================================
# ONE ROW PER PHYSICAL CELL PER PATIENT, which each arm reaches differently:
#
#   massimo1           reads FlowPath_csv_all/<pid>/<pid>.csv, a REAL whole-slide
#                      export. No de-duplication is involved at all.
#   massimo2           has no whole-slide tier, so it pools its region csvs and
#                      de-duplicates on cell_key_cols(). Each of a patient's region
#                      files covers the same slide, so a naive pool would count every
#                      cell two or three times and every cohort-level fraction on the
#                      site would be computed over an inflated denominator while
#                      still looking plausible.
#   massimo1_inverted  the same as massimo2, on arm 1's regions.
#
# THAT ASYMMETRY IS A FREE TEST, and load_arm_cells() takes it. massimo1 is the only
# arm publishing both a whole-slide export AND region files, so it is the only place
# the de-duplication the other two arms are FORCED to use can be checked rather than
# trusted: run the procedure on arm 1's region files and it should reproduce arm 1's
# own export. arm_wholeslide_reconciliation() prints the comparison, and flags the
# one direction that cannot be explained by the regions covering less than the slide.
#
# Which region a given cell belongs to is NOT decided here — only the polygon can say
# that, and that is membership.R's job via membership_data(<arm>). This file only loads.
load_arm_cells <- function(arm = ARM_MODES) {
  arm   <- match.arg(arm)
  spec  <- arm_spec(arm)
  cells <- arm_cells(spec)
  ucell <- arm_union_tier_cells(spec)

  if (nrow(cells) == 0 && nrow(ucell) == 0) {
    # Name the directory AND the likely cause. The commonest reason this fires is a
    # tree that predates the three-arm split: `data/all_slide/` is massimo2's export
    # under its old name, and nothing here reads it any more. Saying so beats a bare
    # "no cells", which reads like an empty dataset rather than an unmade symlink.
    hint <- ""
    legacy <- here::here("data", "all_slide")
    if (identical(arm, "massimo2") && dir.exists(legacy))
      hint <- paste0("\n  data/all_slide/ still exists — that IS this arm's export, ",
                     "under the name it had before the arms were split. Rename or ",
                     "symlink it:  ln -s all_slide data/massimo2")
    warning("load_arm_cells(\"", arm, "\"): no cells under ", spec$root_path,
            "\n  expected the region csvs at ", spec$region_csv$path,
            "\n  symlink the export:  ln -s <share>/", sub("massimo", "Massimo", spec$root),
            " data/", spec$root, hint)
    return(tibble::tibble())
  }

  out <- arm_cohort_cells(spec, cells, ucell)

  # Say what the cohort set actually is, rather than assuming it. Whether a patient's
  # region files repeat the same cells or partition them is a property of the
  # producer, not of the layout, and it sets every cohort-level denominator.
  rep <- arm_overlap_report(cells)
  message(sprintf("LOADED %s: %d cells across %d patients from %d region file(s); %s",
                  arm, nrow(out), dplyr::n_distinct(out$patient_id), nrow(rep),
                  attr(rep, "verdict") %||% "no region files"))
  attr(out, "overlap_report") <- rep

  recon <- arm_wholeslide_reconciliation(spec, cells, ucell)
  if (nrow(recon)) {
    message("  whole-slide check: ", attr(recon, "verdict"))
    attr(out, "reconciliation") <- recon
  }
  out
}

ihc_massimo1          <- load_arm_cells("massimo1")
ihc_massimo2          <- load_arm_cells("massimo2")
ihc_massimo1_inverted <- load_arm_cells("massimo1_inverted")

# The arm -> cell-set lookup, so a page that declares ARM once at the top gets both
# its cells and its pathologist table from that one constant.
ihc_for <- function(arm = ARM_MODES) {
  arm <- match.arg(arm)
  switch(arm,
         massimo1          = ihc_massimo1,
         massimo1_inverted = ihc_massimo1_inverted,
         massimo2          = ihc_massimo2)
}
