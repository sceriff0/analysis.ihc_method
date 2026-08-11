# =============================================================================
# load_data.R  —  the raw inputs, loaded once
#
# Sourced first by every analysis. Defines, from data/ (gitignored):
#   dds              DESeq2 object from counts.RData, with DESeq() already run
#   clinical_data    the clinical CRF (clinical_data.xlsx)
#   neoplastic_data  the pathologist tumour-content scores, inline, LONG:
#                    one row per (SAMPLE, annotation) — see the note at its definition
#   counts_data      normalised counts, wide: one row per bulk-RNA Sample
#   ihc_data         single cells, one row per cell, pooled and DE-DUPLICATED over
#                    the all-slide export's per-region csvs
#
# It prints colnames() for each on load — the callers wrap the source() in
# capture.output() to keep that out of the knitted page.
#
# Derived quantities live in validation_helpers.R, membership rules in
# membership.R, and the cell-table schema in cell_tables.R. Nothing here reduces
# or reshapes: this file only loads.
# =============================================================================
library(tidyverse)
library(DESeq2)
library(here)
library(readxl)
library(fs)
library("data.table")

# The cell-table schema adapter (read_cell_csv() and the accessors every analysis
# reads an export through) plus the all-slide reader, which needs slide_key() from
# validation_helpers. Sourcing validation_helpers here is not a cycle: it never reads
# anything this file defines — `ihc_data` appears in it only as a parameter name.
source(here("code", "validation_helpers.R"))   # pulls in cell_tables.R + plot_theme.R
source(here("code", "all_slide.R"))

dds <- get(load(here("data", "counts.RData")))
dds <- DESeq(dds)

clinical_data <- read_excel(here("data", "clinical_data.xlsx")) |>
  filter(!is.na(`ID PATIENT`)) |>
  mutate(`ID CRF PRESERVE` = gsub("-", ".", `ID CRF PRESERVE`))

# The pathologist's neoplastic-cellularity score, one row per SCORED REGION.
#
# LONG, NOT WIDE (~ANNOTATION_1..3), for one reason: 24086 has no annotation drawn at
# all, and its 75% refers to the WHOLE SLIDE. A wide frame can only put that in
# ANNOTATION_1, which would claim a region the pathologist never drew — and would then
# fail to join, because the metrics frame labels that patient's single row
# `whole_slide`. Long says what was actually scored.
#
# `annotation` here must match the label the membership metrics emit, because
# _clinical_data_body.Rmd joins on (patient_id, annotation):
#   ANNOTATION_<k>  region k, k being the alphabet position of the export's letter
#                   suffix (A -> 1, B -> 2, C -> 3)
#   whole_slide     no annotation directory, so every cell counts (see all_slide.R)
#
# The region counts cross-check against the export exactly — 046 three csvs and three
# geojsons, 052 two, 5456 three, 10338 one, 15897 two, and 24086 a bare csv with no
# annotation directory. A mismatch between the two would silently drop a region from
# the correlation, so clinical_data.Rmd prints the reconciliation.
#
# Values updated 2026-08-11 from the pathologist's re-read. They differ materially
# from the previous set, so the tumour-content correlation is NOT comparable to an
# earlier knit: 046 A 50->30, 052 70/50->50/75, 5456 A 70->80, 10338 75->80,
# 15897 gains B=75 (was single-annotation), and 24086 goes from three annotations
# (80/70/70) to one whole-slide score of 75.
neoplastic_data <- tribble(
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

counts_data <-  counts(dds, normalized = TRUE) |>
  as_tibble(rownames = "GENE") |>
  pivot_longer(cols = -GENE, names_to = "Sample", values_to = "Expression") |>
  pivot_wider(names_from = GENE, values_from = Expression) |>
  filter(Sample %in% clinical_data$`ID CRF PRESERVE`)

# The whole-cohort cell set, from the all-slide export (code/all_slide.R):
#
#   data/all_slide/csv/<patient>/<patient>_<A|B|C>.csv
#
# ONE ROW PER PHYSICAL CELL PER PATIENT. The export writes one csv per annotation
# REGION, and each of a patient's region files covers the same slide — so a naive
# pool counts every cell two or three times, and every cohort-level fraction on the
# site would be computed over an inflated denominator while still looking plausible.
# all_slide_union_cells() keys on cell_key_cols() (an id column plus the centroid
# pair) and keeps one copy.
#
# Which region a given cell belongs to is NOT decided here — only the polygon can say
# that, and that is code/membership.R's job via membership_data("all_slide"). This
# file only loads.
load_ihc_data <- function(root = ALL_SLIDE_DIR) {
  cells <- all_slide_cells(root)
  if (nrow(cells) == 0) {
    warning("load_ihc_data(): no cells under ", root,
            " — copy or symlink the all-slide export to data/all_slide/")
    return(cells)
  }
  out <- all_slide_union_cells(cells)
  # Say what the de-duplication actually did. Whether a patient's region files repeat
  # the same cells or partition them is a property of the producer, not of the layout,
  # and it sets every cohort-level denominator — so it is reported rather than assumed.
  rep <- all_slide_overlap_report(cells)
  message(sprintf("LOADED %d cells across %d patients from %d region file(s); %s",
                  nrow(out), dplyr::n_distinct(out$patient_id), nrow(rep),
                  attr(rep, "verdict")))
  attr(out, "overlap_report") <- rep
  out
}

ihc_data <- load_ihc_data()

colnames(clinical_data)
colnames(neoplastic_data)
colnames(counts_data)
colnames(ihc_data)

