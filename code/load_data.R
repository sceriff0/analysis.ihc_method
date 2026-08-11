# =============================================================================
# load_data.R  —  the raw inputs, loaded once
#
# Sourced first by every analysis. Defines, from data/ (gitignored):
#   dds              DESeq2 object from counts.RData, with DESeq() already run
#   clinical_data    the clinical CRF (clinical_data.xlsx)
#   neoplastic_data  the pathologist tumour-content scores, inline (six patients)
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

neoplastic_data <- tribble(
  ~SAMPLE, ~ANNOTATION_1, ~ANNOTATION_2, ~ANNOTATION_3,
  "24086",              80,            70,            70,
  "15897",              70,            NA,            NA,
  "052",                70,            50,            NA,
  "046",                50,            60,            50,
  "5456",               70,            80,            70,
  "10338",               75,            NA,            NA
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

