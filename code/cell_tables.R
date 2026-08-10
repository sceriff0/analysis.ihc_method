# =============================================================================
# cell_tables.R  —  one reader for every single-cell table this project consumes
#
# Three upstream tools export "one row per cell", with three different spellings
# of the same four facts (identity, phenotype, membership, positivity). Every
# analysis here goes through the accessors below instead of naming a column, so
# swapping the export format changes nothing downstream.
#
#   A. FlowPath  PhenotypeCsvExporter            <patient>.csv, <patient>_a<k>.csv
#      cell_id, phenotype, Out_of_annotation, Outlier, centroid_x, centroid_y (um),
#      area, perimeter, eccentricity, solidity, <col>_raw/_zscore/_sign
#      positive sign = "+"; blank = the gate tree never touched that column.
#
#   B. mirage  bin/phenotype_cells.py            <patient>_phenotypes.csv
#      label, phenotype, candidates, outcome, provenance,
#      pheno_score:<TYPE>, p_neg:<M>, p_pos:<M>, sign:<M>, state:<M>
#      sign  "1" pos / "0" neg / "·" never gated / "x" contradictory
#      state  1 pos / -1 neg / 0 free. No centroids and no membership flag: it is
#      keyed on the segmentation `label` and is meant to be joined, not used alone.
#
#   C. mirage  bin/join_flowpath.py --out-table  cohort table (csv/parquet)
#      B joined onto the SpatialData store: phenotype, label, x_px, y_px (PIXELS),
#      fp_out_of_annotation, fp_outlier, fp_matched, <M>, <M>_zscore, <M>_positive.
#
# The traps these accessors exist to close:
#   - `phenotype` is a free-text gate-branch name. This cohort's tree writes
#     "CD8 T cell (T_cytotoxic)", so the project reads the parenthetical — but B
#     and C write bare labels, and a plain regex silently yields NA for every cell.
#   - membership is `Out_of_annotation` in A and `fp_out_of_annotation` in C.
#   - positivity is `<M>_sign` / `sign:<M>` / `<M>_positive`.
#   - A's centroids are MICRONS, C's are already PIXELS. Dividing C by um_per_px
#     a second time puts every cell outside its annotation.
#
# Loaded by validation_helpers.R; needs nothing but base R + tibble.
# =============================================================================

# --- Dialect ----------------------------------------------------------------
# Which of the three exports a table came from. Only used for messages and for
# the centroid frame; every accessor probes columns directly, so an export that
# mixes dialects (a FlowPath csv joined to mirage scores) still resolves.
cell_dialect <- function(cells) {
  nm <- names(cells)
  if (any(startsWith(nm, "sign:")) || any(startsWith(nm, "pheno_score:"))) "mirage_phenotypes"
  else if ("fp_out_of_annotation" %in% nm || any(endsWith(nm, "_positive")))  "mirage_cohort"
  else "flowpath"
}

# --- Phenotype label --------------------------------------------------------
# The cell-type label to group by. This cohort's FlowPath gate tree encodes the
# type in parentheses ("CD8 T cell (T_cytotoxic)" -> "T_cytotoxic"); mirage and
# any other gate tree write the type bare. Take the parenthetical where there is
# one, the trimmed label otherwise, NA only for a genuinely empty cell.
clean_phenotype <- function(x) {
  x     <- as.character(x)
  inner <- sub("^[^(]*\\(([^)]*)\\).*$", "\\1", x)
  out   <- trimws(ifelse(grepl("\\(([^)]*)\\)", x), inner, x))
  ifelse(is.na(x) | out == "", NA_character_, out)
}

# `phenotype_clean` for a cell table, computed if the loader has not already.
cell_phenotype <- function(cells) {
  if ("phenotype_clean" %in% names(cells)) return(cells$phenotype_clean)
  if ("phenotype" %in% names(cells))       return(clean_phenotype(cells$phenotype))
  rep(NA_character_, nrow(cells))
}

# --- Membership -------------------------------------------------------------
# TRUE where the exporter flagged the cell as OUTSIDE the annotation. Covers
# FlowPath's "True"/"False" strings and the booleans mirage writes.
OUTSIDE_COL <- c("Out_of_annotation", "fp_out_of_annotation", "out_of_annotation")

has_outside_flag <- function(cells) any(OUTSIDE_COL %in% names(cells))

is_outside <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "t")
}

# Per-cell "outside the annotation" for a whole table; all-FALSE (everything
# inside) when no export carried a flag, so a flag-less table degrades to
# "one region = the whole slide" rather than erroring.
cell_outside <- function(cells) {
  col <- intersect(OUTSIDE_COL, names(cells))[1]
  if (is.na(col)) return(rep(FALSE, nrow(cells)))
  is_outside(cells[[col]])
}

# Every boolean-ish flag an export can carry. Coerced to real logical in place at
# read time, because readers guess the type from the values they happen to see —
# readr turns "True"/"False" into logical, fread and a hand-written csv leave it
# character, an export with no flagged cells at all can parse as all-NA — and
# binding two exports of the same cohort then fails on a type clash long before
# anyone looks at the numbers.
FLAG_COLS <- c(OUTSIDE_COL, "Outlier", "fp_outlier", "fp_matched")

normalise_cell_flags <- function(cells) {
  for (col in intersect(FLAG_COLS, names(cells)))
    cells[[col]] <- is_outside(cells[[col]])
  cells
}

# --- Positivity -------------------------------------------------------------
# TRUE where a sign column marks the cell positive for a marker. "1" is mirage's
# pos code and FlowPath's boolean-ish drift; "·"/"x"/"" are never-gated and
# contradictory, and read as not-positive.
is_pos <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  tolower(trimws(as.character(x))) %in% c("+", "pos", "positive", "yes", "true", "1")
}

# The column carrying marker positivity, in preference order across the three
# dialects. NA when the export never measured that marker.
marker_sign_col <- function(cells, marker) {
  cand <- c(paste0(marker, "_sign"), paste0("sign:", marker),
            paste0(marker, "_positive"), paste0("state:", marker))
  intersect(cand, names(cells))[1]
}

# Positivity for one marker. All-FALSE for a marker the export never gated, so a
# marker-defined subset degrades to n = 0 instead of aborting the report — check
# `marker_sign_col()` when a population is unexpectedly empty.
# `state:<M>` is three-valued (1 pos / -1 neg / 0 free), so only +1 counts.
marker_pos <- function(cells, marker) {
  col <- marker_sign_col(cells, marker)
  if (is.na(col)) return(rep(FALSE, nrow(cells)))
  v <- cells[[col]]
  if (startsWith(col, "state:")) return(suppressWarnings(as.numeric(v)) %in% 1)
  is_pos(v)
}

# The raw sign value as the export wrote it ("+"/"-"/"" for FlowPath, "1"/"0"/"·"/"x"
# for mirage, TRUE/FALSE for the cohort table). Mostly for display; prefer
# marker_pos() to decide positivity and marker_gated() to decide "was it measured".
marker_sign <- function(cells, marker) {
  col <- marker_sign_col(cells, marker)
  if (is.na(col)) return(rep(NA_character_, nrow(cells)))
  as.character(cells[[col]])
}

# TRUE where a gate actually evaluated this marker on this cell. All three exports
# distinguish "negative" from "never measured", and conflating them understates
# positivity in a way that reads as biology: FlowPath writes a blank sign, mirage
# writes "·" (free) or "x" (contradictory), and an absent column means the panel
# never carried the marker at all.
UNGATED_SIGNS <- c("", "·", "x", "na", "nan")

marker_gated <- function(cells, marker) {
  col <- marker_sign_col(cells, marker)
  if (is.na(col)) return(rep(FALSE, nrow(cells)))
  v <- cells[[col]]
  if (is.logical(v)) return(!is.na(v))
  if (startsWith(col, "state:")) return(suppressWarnings(as.numeric(v)) %in% c(-1, 1))
  !is.na(v) & !tolower(trimws(as.character(v))) %in% UNGATED_SIGNS
}

# Positivity for several markers at once: a tibble of logical columns named for the
# MARKERS, not for whichever column each export happened to use. This is what an
# analysis wants when it builds a co-expression matrix — `all_of(paste0(m, "_sign"))`
# both hard-codes the FlowPath spelling and aborts on a marker the panel lacks.
marker_matrix <- function(cells, markers) {
  tibble::as_tibble(stats::setNames(
    lapply(markers, function(m) marker_pos(cells, m)), markers))
}

# The intensity columns behind a marker, again in preference order per dialect:
#   z-score   <M>_zscore   (FlowPath and the mirage cohort table both use it)
#   value     <M>_raw      (FlowPath) or a bare <M> column (mirage cohort table)
# NA_real_ for a marker the export does not carry, so a partial panel yields empty
# panels rather than an error.
.marker_col <- function(cells, cand) {
  hit <- intersect(cand, names(cells))[1]
  if (is.na(hit)) NULL else as.numeric(cells[[hit]])
}

marker_zscore <- function(cells, marker) {
  .marker_col(cells, paste0(marker, c("_zscore", "_z"))) %||% rep(NA_real_, nrow(cells))
}

marker_value <- function(cells, marker) {
  .marker_col(cells, c(paste0(marker, "_raw"), marker)) %||% rep(NA_real_, nrow(cells))
}

# Default-on-empty. Same definition in aggregation_compare.R, cell_tables.R and
# registration_accuracy_plots.R: these files are sourced in different orders by
# different reports, so a variant that only tested is.null() would silently take
# over and let a zero-length result through where a default was intended. (base R's
# %||% is the is.null-only form, which this deliberately widens.)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# --- Coordinates ------------------------------------------------------------
# Cell centroids in the geojson PIXEL frame, which is what point-in-polygon
# needs. FlowPath exports microns and is rescaled; mirage's cohort table already
# exports pixels and is passed through untouched.
cell_centroids_px <- function(cells, um_per_px = 0.325) {
  if (all(c("x_px", "y_px") %in% names(cells)))
    return(data.frame(x = as.numeric(cells$x_px), y = as.numeric(cells$y_px)))
  if (!all(c("centroid_x", "centroid_y") %in% names(cells)))
    stop("cell table has no centroids: expected centroid_x/centroid_y (um) or x_px/y_px")
  data.frame(x = as.numeric(cells$centroid_x) / um_per_px,
             y = as.numeric(cells$centroid_y) / um_per_px)
}

# TRUE where a cell has usable coordinates in either frame.
has_centroids <- function(cells) {
  all(c("centroid_x", "centroid_y") %in% names(cells)) ||
    all(c("x_px", "y_px") %in% names(cells))
}

# --- Identity ---------------------------------------------------------------
# Columns that identify a cell for the union dedup (a cell inside two annotations
# must count once). Take an id column AND the centroid pair together, because the
# two readings of FlowPath's `cell_id` disagree and only the pair is safe under
# both: if a patient's per-annotation csvs share one whole-slide CellIndex, the id
# is already unique and the centroid adds nothing; if each csv is re-indexed from
# 0, the id alone would fuse two different cells that happen to share a row
# number. mirage's `label` is a segmentation id and is unique either way.
CELL_ID_COLS     <- c("label", "cell_id")
CELL_COORD_COLS  <- list(c("centroid_x", "centroid_y"), c("x_px", "y_px"))

cell_key_cols <- function(cells) {
  id    <- intersect(CELL_ID_COLS, names(cells))[1]
  coord <- Filter(function(k) all(k %in% names(cells)), CELL_COORD_COLS)
  c(if (!is.na(id)) id, if (length(coord)) coord[[1]])
}

# --- Loading ----------------------------------------------------------------
# Read one cell csv and attach the canonical derived columns. `patient_id` is the
# caller's (it comes from the filename, never from the file). Returns a 0-row
# tibble and warns — rather than aborting a whole cohort load — when the file is
# not a recognisable cell table.
read_cell_csv <- function(path, patient_id = NULL, required = "phenotype") {
  # fread where it is available (cell tables run to millions of rows), readr
  # otherwise. Either way normalise_cell_flags() settles the type differences.
  cells <- if (requireNamespace("data.table", quietly = TRUE))
    tibble::as_tibble(data.table::fread(path, showProgress = FALSE))
  else
    suppressWarnings(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  miss <- setdiff(required, names(cells))
  if (length(miss)) {
    warning("skipping ", path, " — not a cell table (missing ",
            paste(miss, collapse = ", "), ")")
    return(tibble::tibble())
  }
  cells <- normalise_cell_flags(cells)
  cells$phenotype_clean <- cell_phenotype(cells)
  if (!is.null(patient_id)) cells$patient_id <- patient_id
  cells
}
