# =============================================================================
# mirage_cells.R  —  cells from mirage's own phenotyping
#
# The fourth cell source, alongside the FlowPath exports. mirage phenotypes every
# cell itself (bin/phenotype_cells.py, driven by panel.yaml) and publishes the
# result as plain per-patient CSVs, so no zarr reader is needed:
#
#   data/mirage/<patient>/phenotyping/phenotypes.csv        REQUIRED
#       label, phenotype, candidates, outcome, provenance,
#       pheno_score:<TYPE>, p_neg:<M>, p_pos:<M>, sign:<M>, state:<M>
#   data/mirage/<patient>/cell_properties/morphology.csv    REQUIRED
#       label, y, x (PIXELS), area, eccentricity, perimeter, convex_area,
#       axis_major_length, axis_minor_length, solidity
#   data/mirage/<patient>/quantification/merged_quant.csv   optional
#       label, "<MARKER>: <Compartment>: <Statistic>" intensities
#
# `<patient>` is the directory name — nothing inside the files identifies the
# patient. The three tables are joined on `label`, mirage's segmentation id, which
# is a real cell identity within a patient (unlike FlowPath's positional cell_id).
#
# WHY THE JOIN IS MANDATORY, NOT A CONVENIENCE
# phenotypes.csv carries no coordinates and no annotation flag. On its own it
# cannot answer "is this cell inside the tumour?" at all — which is why the mirage
# source only supports geojson membership (see membership.R). morphology.csv is
# what supplies the centroids, so a patient missing it is skipped rather than
# silently contributing every cell as "inside".
#
# COORDINATE FRAME
# morphology's x/y are regionprops centroids on the SEGMENT mask, i.e. pixels in
# the registered reference frame — the same frame the pathologist's QuPath geojson
# is drawn in. They are renamed x_px/y_px on load so cell_centroids_px() passes
# them through untouched; a micron export would be rescaled, and doing both would
# put every cell outside its annotation. If the mirage page shows implausibly few
# in-annotation cells, that assumption is the first thing to check — the
# annotation_membership_qc page exists for exactly this class of question.
#
# Depends on cell_tables.R (read_cell_csv, clean_phenotype) and on slide_key()
# from validation_helpers.R.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
  library(purrr)
})

MIRAGE_DIR <- here::here("data", "mirage")

# Where each table sits under a patient directory, relative to it.
MIRAGE_FILES <- c(phenotypes = "phenotyping/phenotypes.csv",
                  morphology = "cell_properties/morphology.csv",
                  quant      = "quantification/merged_quant.csv")

# Is this directory a phenotyped patient? Structure, not name: a mirage outdir has
# cohort-level directories sitting next to the patient ones — `phenotyping/` (the
# compiled panel from COMPILE_PANEL), `qc/`, `size_logs/`, `_UNROUTED_PUBLISH/` —
# and treating those as broken patients would warn four times on every knit. The
# presence of phenotypes.csv is what makes a directory a patient, so data/mirage can
# point straight at (or symlink) a raw outdir.
is_mirage_patient_dir <- function(dir) {
  file.exists(file.path(dir, MIRAGE_FILES[["phenotypes"]]))
}

# One patient. Returns a 0-row tibble and warns — rather than aborting the cohort
# load — when a directory IS a patient but a required table is missing, so a
# partially-processed results tree still yields the patients it does have.
read_mirage_patient <- function(dir, patient_id = fs::path_file(dir)) {
  paths <- file.path(dir, MIRAGE_FILES)
  names(paths) <- names(MIRAGE_FILES)

  if (!is_mirage_patient_dir(dir)) return(tibble::tibble())   # not a patient; silent
  if (!file.exists(paths[["morphology"]])) {
    warning("mirage: skipping ", patient_id, " — has ", MIRAGE_FILES[["phenotypes"]],
            " but no ", MIRAGE_FILES[["morphology"]], ", so its cells have no ",
            "coordinates and cannot be placed in an annotation")
    return(tibble::tibble())
  }

  pheno <- read_cell_csv(paths[["phenotypes"]])
  morph <- read_cell_csv(paths[["morphology"]], required = "label")
  if (nrow(pheno) == 0 || nrow(morph) == 0) return(tibble::tibble())

  # morphology carries no phenotype, so read_cell_csv leaves phenotype_clean NA on
  # it; drop that column before the join or it would overwrite the real one.
  morph <- morph |>
    dplyr::select(-dplyr::any_of("phenotype_clean")) |>
    dplyr::rename(x_px = x, y_px = y)

  cells <- dplyr::inner_join(pheno, morph, by = "label")
  if (nrow(cells) < nrow(pheno))
    warning("mirage: ", patient_id, " — ", nrow(pheno) - nrow(cells),
            " phenotyped cell(s) have no morphology row and were dropped")

  if (file.exists(paths[["quant"]])) {
    quant <- read_cell_csv(paths[["quant"]], required = "label")
    if (nrow(quant)) cells <- dplyr::left_join(
      cells, dplyr::select(quant, -dplyr::any_of("phenotype_clean")), by = "label")
  }

  dplyr::mutate(cells, patient_id = slide_key(patient_id), .before = 1)
}

# Every patient directory under `root`, as one long cell table. Empty (with a
# warning) when the directory does not exist, so a report that offers the mirage
# variant degrades to "no patients" instead of erroring on a clone with no data.
load_mirage_cells <- function(root = MIRAGE_DIR) {
  if (!fs::dir_exists(root)) {
    warning("mirage: directory not found: ", root)
    return(tibble::tibble())
  }
  dirs <- as.character(fs::dir_ls(root, type = "directory"))
  patients <- dirs[vapply(dirs, is_mirage_patient_dir, logical(1))]
  if (!length(patients)) {
    warning("mirage: no patient directories under ", root,
            " (a patient directory is one containing ", MIRAGE_FILES[["phenotypes"]], ")")
    return(tibble::tibble())
  }
  purrr::map_dfr(patients, read_mirage_patient)
}

# Provenance table for the report: which patients loaded, how many cells, and how
# many of them mirage actually resolved to a type. A large unresolved share is the
# headline caveat for any mirage panel — those cells stay in the denominators.
mirage_cells_inventory <- function(cells) {
  if (nrow(cells) == 0) return(tibble::tibble())
  cells |>
    dplyr::mutate(.unresolved = is_unresolved_phenotype(phenotype_clean)) |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(
      n_cells      = dplyr::n(),
      n_resolved   = sum(!.unresolved),
      pct_resolved = round(100 * mean(!.unresolved), 1),
      n_phenotypes = dplyr::n_distinct(phenotype_clean[!.unresolved]),
      .groups      = "drop")
}
