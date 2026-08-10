# =============================================================================
# load_data.R  —  the raw inputs, loaded once
#
# Sourced first by every analysis. Defines, from data/ (gitignored):
#   dds              DESeq2 object from counts.RData, with DESeq() already run
#   clinical_data    the clinical CRF (clinical_data.xlsx)
#   neoplastic_data  the pathologist tumour-content scores, inline (six patients)
#   counts_data      normalised counts, wide: one row per bulk-RNA Sample
#   ihc_data         single cells, one row per cell, pooled over data/flowpath/*.csv
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

# The cell-table schema adapter: read_cell_csv() and the accessors every analysis
# reads an export through. validation_helpers.R sources it too, but this file runs
# first and needs it here. See code/README.md.
source(here("code", "cell_tables.R"))

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

# ONE csv per patient: "<patient>.csv" (e.g. 046.csv), one row per cell; patient_id
# is the filename stem, never a column of the file. read_cell_csv() (cell_tables.R)
# does the rest: it derives phenotype_clean, coerces the boolean flags to logical so
# patients can be bound together, and accepts any of the three export schemas.
# No per-annotation splitting here — that is code/membership.R's job.
process_patient <- function(csv_path) {
  patient_id <- path_ext_remove(path_file(csv_path))     # "046.csv" -> "046"
  message(sprintf("LOADING PATIENT: %s", patient_id))
  read_cell_csv(csv_path, patient_id = patient_id)
}

load_ihc_data <- function() {
  csv_paths    <- dir_ls(here("data", "flowpath"), glob = "*.csv")  # top level only; skips old/
  safe_process <- possibly(process_patient, otherwise = NULL, quiet = FALSE)
  map(csv_paths, safe_process) |> list_rbind()
}

ihc_data <- load_ihc_data()

colnames(clinical_data)
colnames(neoplastic_data)
colnames(counts_data)
colnames(ihc_data)

