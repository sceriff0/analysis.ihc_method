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

# Is this directory a patient? Structure, not name: a mirage outdir has
# cohort-level directories sitting next to the patient ones — `phenotyping/` (the
# compiled panel from COMPILE_PANEL), `qc/`, `size_logs/`, `_UNROUTED_PUBLISH/` —
# and treating those as broken patients would warn four times on every knit.
#
# ANY of the three tables marks a patient, not phenotypes.csv alone. Phenotyping
# is an OPTIONAL pipeline stage: mirage ships it as a separate feature and a run
# built without it emits quantification/ and cell_properties/ but no
# phenotyping/. Gating on phenotypes.csv classified every such patient as "not a
# patient" — silently, since a non-patient directory is skipped without comment —
# and the whole cohort surfaced as one "no patient directories" warning that
# reads like an empty dataset rather than a pipeline built to a different spec.
is_mirage_patient_dir <- function(dir) {
  any(file.exists(file.path(dir, MIRAGE_FILES)))
}

# One patient. Returns a 0-row tibble and warns — rather than aborting the cohort
# load — when a directory IS a patient but a required table is missing, so a
# partially-processed results tree still yields the patients it does have.
read_mirage_patient <- function(dir, patient_id = fs::path_file(dir)) {
  paths <- file.path(dir, MIRAGE_FILES)
  names(paths) <- names(MIRAGE_FILES)

  if (!is_mirage_patient_dir(dir)) return(tibble::tibble())   # not a patient; silent
  if (!file.exists(paths[["morphology"]])) {
    warning("mirage: skipping ", patient_id, " — no ", MIRAGE_FILES[["morphology"]],
            ", so its cells have no coordinates and cannot be placed in an annotation")
    return(tibble::tibble())
  }

  morph <- read_cell_csv(paths[["morphology"]], required = "label")
  if (nrow(morph) == 0) return(tibble::tibble())
  morph <- dplyr::rename(morph, x_px = x, y_px = y)

  if (file.exists(paths[["phenotypes"]])) {
    pheno <- read_cell_csv(paths[["phenotypes"]])
    if (nrow(pheno) == 0) return(tibble::tibble())
    # morphology carries no phenotype, so read_cell_csv leaves phenotype_clean NA on
    # it; drop that column before the join or it would overwrite the real one.
    cells <- dplyr::inner_join(
      pheno, dplyr::select(morph, -dplyr::any_of("phenotype_clean")), by = "label")
    if (nrow(cells) < nrow(pheno))
      warning("mirage: ", patient_id, " — ", nrow(pheno) - nrow(cells),
              " phenotyped cell(s) have no morphology row and were dropped")
  } else {
    # No phenotyping stage in this run. Keep the cells: their coordinates are what
    # the membership metrics need, and phenotype is a STRATIFIER over those metrics,
    # not an input to them — counts, areas and densities are all still exact. What
    # is lost is the per-type breakdown, and phenotype_clean stays NA, which
    # is_unresolved_phenotype() already treats as unresolved so every downstream
    # denominator stays honest. load_mirage_cells() says so once for the cohort.
    cells <- morph
  }

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
            " (a patient directory is one containing any of: ",
            paste(MIRAGE_FILES, collapse = ", "), ")")
    return(tibble::tibble())
  }

  # Say it once, for the cohort, rather than once per patient. A run built without
  # the phenotyping stage is a legitimate configuration, not a broken copy — but a
  # panel that silently shows every cell as unresolved is indistinguishable from a
  # panel where phenotyping ran and failed, so the distinction is stated here.
  n_pheno <- sum(file.exists(file.path(patients, MIRAGE_FILES[["phenotypes"]])))
  if (n_pheno == 0)
    warning("mirage: no ", MIRAGE_FILES[["phenotypes"]], " under ", root,
            " — this pipeline run had no phenotyping stage. Cells, coordinates and ",
            "densities are exact; every cell is unresolved, so per-type panels will ",
            "be empty and phenotype shares are not interpretable.")
  else if (n_pheno < length(patients))
    warning("mirage: ", length(patients) - n_pheno, " of ", length(patients),
            " patient(s) have no ", MIRAGE_FILES[["phenotypes"]],
            " — their cells load unresolved")

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
