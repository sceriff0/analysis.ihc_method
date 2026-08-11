# =============================================================================
# all_slide.R  —  the "all-slide" export layout: one CSV *and* one geojson per
# annotation region, nested under a per-patient directory.
#
# This is the fifth cell source, and the only one where the pathologist polygon
# and the cell export are matched 1:1 at the REGION level. The producer writes:
#
#   <root>/csv/<patient>/<patient>_<L>.csv          one export per region
#   <root>/csv/<patient>/<patient>.csv              whole-slide export (no regions)
#   <root>/annotation/<patient>/<patient>_<L>.geojson   the polygon for region <L>
#
# with <L> a LETTER (A, B, C ...), not the `_a<k>` digit the older per-annotation
# layout uses. On the cluster this is
# /hpcnfs/techunits/imaging/work/ATTEND/Mirage/all-slide_new; point ALL_SLIDE_DIR
# at a copy or a symlink of it.
#
# LETTERS ARE POSITIONS, NOT NAMES. A -> ANNOTATION_1, B -> ANNOTATION_2, ... so
# the regions line up with neoplastic_data's ANNOTATION_1..3 columns and with every
# other membership mode. A patient whose regions are named out of order, or that
# skips a letter, still maps by the letter's position in the alphabet — never by
# the order the files happen to be listed in.
#
# NO ANNOTATION DIRECTORY MEANS EVERYTHING IS INSIDE.
# This is the layout's own convention, stated by the people who produced it: a
# patient with no `annotation/<patient>/` contributes its whole slide as one
# region. That is a DIFFERENT statement from "this patient's polygon is missing",
# which every other mode treats as a reason to drop the patient, and it is why
# this file exists rather than the older loaders being widened — silently
# reinterpreting a missing polygon as "all in" everywhere else would turn a data
# problem into a 100%-inside patient without anyone noticing. Here it is the
# documented meaning, so it is applied deliberately and recorded in `source` as
# "whole_slide" rather than "sf".
#
# MEMBERSHIP INSIDE A REGION, in order of preference:
#   1. "sf"          the region's own geojson, point-in-polygon. Primary: it is the
#                    pathologist's line, and it is the only source that knows the
#                    region's AREA, so it is the only one that yields densities.
#   2. "flag"        the export's Out_of_annotation flag, when the region has a CSV
#                    but no readable polygon.
#   3. "whole_slide" every row counts, for a patient with no annotation directory.
#
# Depends on cell_tables.R (read_cell_csv and the accessors) and, at use time, on
# validation_helpers.R for slide_key() / region_ratios_area(). sf is a LAZY
# dependency, exactly as in validation_helpers.R: a cohort with no annotation
# directory at all must still load on a machine without sf.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
  library(purrr)
  library(tibble)
})

ALL_SLIDE_DIR <- here::here("data", "all_slide")

# Region letter -> the ANNOTATION_<k> label the rest of the project speaks.
# Position in the alphabet, so "C" is region 3 whether or not B exists.
all_slide_region <- function(letter) {
  idx <- match(toupper(letter), LETTERS)
  ifelse(is.na(idx), NA_character_, paste0("ANNOTATION_", idx))
}

# Split "<patient>_<L>" into its parts. A stem with no "_<L>" suffix is the
# whole-slide export, so it comes back with letter = NA.
.all_slide_parts <- function(stem) {
  m <- regmatches(stem, regexec("^(.*)_([A-Za-z])$", stem))[[1]]
  if (length(m) == 3) list(patient = m[2], letter = toupper(m[3]))
  else                list(patient = stem,  letter = NA_character_)
}

# --- Annotations -------------------------------------------------------------
# Every polygon under <root>/annotation/, as an sf frame with the same three
# columns load_annotations() returns, so anything that already consumes `annots`
# consumes these. Returns an EMPTY frame (not an error) when the directory is
# absent: a cohort can legitimately be all-whole-slide, and that must not stop the
# report at the loader.
all_slide_annotations <- function(root = ALL_SLIDE_DIR, patient_ids = NULL) {
  ann_root <- file.path(root, "annotation")
  if (!dir.exists(ann_root)) {
    message("all_slide: no annotation/ under ", root, " — every patient is whole-slide")
    return(NULL)
  }
  .require_sf("all_slide_annotations()")

  files <- fs::dir_ls(ann_root, glob = "*.geojson", recurse = TRUE, type = "file")
  if (length(files) == 0) {
    message("all_slide: annotation/ exists but holds no geojson — every patient is whole-slide")
    return(NULL)
  }

  rows <- purrr::map(as.character(files), function(path) {
    parts <- .all_slide_parts(fs::path_ext_remove(fs::path_file(path)))
    # The directory name is the patient of record; the filename stem only has to
    # agree with it. Trusting the stem alone would mis-key a file someone renamed.
    pid   <- fs::path_file(fs::path_dir(path))
    if (!is.na(parts$patient) && norm_id(parts$patient) != norm_id(pid))
      warning("all_slide: ", path, " sits under patient ", pid,
              " but its name says ", parts$patient, " — using the directory")
    # A lone bare <patient>.geojson is that patient's only region, i.e. region 1.
    ann <- if (is.na(parts$letter)) "ANNOTATION_1" else all_slide_region(parts$letter)
    geom <- tryCatch(read_polygon_geojson(path), error = function(e) {
      warning("all_slide: unreadable polygon, skipping ", path, " — ", conditionMessage(e))
      NULL
    })
    if (is.null(geom)) return(NULL)
    sf::st_sf(patient_id = slide_key(pid), annotation = ann, geometry = geom)
  })
  rows <- purrr::compact(rows)
  if (!length(rows)) return(NULL)

  out <- do.call(rbind, rows)
  if (!is.null(patient_ids))
    out <- out[norm_id(out$patient_id) %in% norm_id(patient_ids), , drop = FALSE]
  if (nrow(out) == 0) return(NULL)
  out
}

# --- Cells -------------------------------------------------------------------
# Every CSV under <root>/csv/, one long table keyed by (patient_id, annotation).
# `annotation` is the region the FILE belongs to — which region a given cell falls
# in is decided later, by all_slide_metrics(), because only the polygon can say so.
# A patient whose directory holds a bare <patient>.csv alongside region files is a
# contradiction in the layout, so the bare file is dropped with a warning rather
# than pooled in and double-counting every cell.
all_slide_cells <- function(root = ALL_SLIDE_DIR) {
  csv_root <- file.path(root, "csv")
  if (!dir.exists(csv_root)) {
    warning("all_slide: no csv/ under ", root, " — no cells loaded")
    return(tibble::tibble())
  }
  ann_root  <- file.path(root, "annotation")
  pat_dirs  <- as.character(fs::dir_ls(csv_root, type = "directory"))
  if (!length(pat_dirs)) {
    warning("all_slide: csv/ holds no patient directories — no cells loaded")
    return(tibble::tibble())
  }

  purrr::map_dfr(pat_dirs, function(dir) {
    pid   <- fs::path_file(dir)
    files <- as.character(fs::dir_ls(dir, glob = "*.csv", type = "file"))
    if (!length(files)) {
      warning("all_slide: patient ", pid, " has a csv/ directory but no csv — skipping")
      return(tibble::tibble())
    }
    # The layout's own rule: no annotation directory for this patient => every cell
    # is inside. Recorded per row so a report can tell these patients apart.
    has_ann <- dir.exists(file.path(ann_root, pid)) &&
               length(fs::dir_ls(file.path(ann_root, pid), glob = "*.geojson",
                                 type = "file", fail = FALSE)) > 0

    parsed  <- purrr::map(fs::path_ext_remove(fs::path_file(files)), .all_slide_parts)
    letters <- vapply(parsed, function(p) p$letter %||% NA_character_, character(1))
    regional <- !is.na(letters)
    if (any(regional) && any(!regional)) {
      warning("all_slide: patient ", pid, " has both region files (",
              paste(letters[regional], collapse = ", "),
              ") and a bare whole-slide csv — dropping the bare file, the regions win")
      files <- files[regional]; letters <- letters[regional]
    }

    purrr::map2_dfr(files, letters, function(path, letter) {
      ann   <- if (is.na(letter)) "whole_slide" else all_slide_region(letter)
      cells <- read_cell_csv(path, patient_id = slide_key(pid))
      if (nrow(cells) == 0) {
        warning("all_slide: ", path, " has no cells — skipping")
        return(tibble::tibble())
      }
      dplyr::mutate(cells,
                    annotation         = ann,
                    file_origin        = "all_slide",
                    has_annotation     = has_ann,
                    .membership_source = if (has_ann) "sf" else "whole_slide")
    })
  })
}

# Provenance: what was actually read, one row per file, for the report to print.
# `n_inside` is left to the metrics frame — this table answers "did my data arrive",
# which has to be answerable even when the geometry step later fails.
all_slide_inventory <- function(cells) {
  if (nrow(cells) == 0) return(tibble::tibble())
  cells |>
    dplyr::group_by(patient_id, annotation) |>
    dplyr::summarise(n_cells        = dplyr::n(),
                     has_annotation = any(has_annotation),
                     membership     = paste(sort(unique(.membership_source)), collapse = "/"),
                     .groups = "drop") |>
    dplyr::arrange(patient_id, annotation)
}

# --- Metrics -----------------------------------------------------------------
# One metrics row per (patient_id, annotation), schema-identical to
# ihc_annotation_metrics() so membership_data() can hand either to the same report.
#
# Per region: the cells of THAT region's file, cut down to the ones inside THAT
# region's polygon. Two files never contribute to the same region, so nothing is
# double-counted even where the polygons overlap.
#
# scope = "union" dissolves a patient's polygons and pools its files first, so a
# cell in the overlap of two regions is counted once — which is the whole reason
# the union is computed from geometry rather than by summing the per-region rows.
all_slide_metrics <- function(cells, annots, scope = c("per_annotation", "union"),
                              um_per_px = 0.325) {
  scope <- match.arg(scope)
  if (nrow(cells) == 0) return(tibble::tibble())

  has_poly <- !is.null(annots) && nrow(annots) > 0
  ann_key  <- if (has_poly) paste(slide_key(annots$patient_id), annots$annotation) else character(0)

  purrr::map_dfr(unique(cells$patient_id), function(pid) {
    cp <- dplyr::filter(cells, patient_id == pid)

    # -- whole-slide patient: no polygon by design, so every cell counts and there
    # -- is no area to divide by. NA densities, exactly as the flag modes return.
    if (!any(cp$has_annotation)) {
      return(region_ratios_area(cp, 0, um_per_px) |>
               dplyr::mutate(patient_id = pid, annotation = "whole_slide",
                             source = "whole_slide", .before = 1))
    }

    if (scope == "union") {
      polys <- if (has_poly) dplyr::filter(annots, slide_key(patient_id) == pid)
               else NULL
      # An annotated patient whose polygons could not be read (no sf, unreadable
      # geojson) must still produce a union row: dropping it here would shrink the
      # cohort in every union-scoped panel without saying so. Fall back to the
      # export's own flag on the pooled slide, and label the row `flag` so the
      # different membership rule is visible in the frame rather than inferred.
      if (is.null(polys) || nrow(polys) == 0) {
        pooled <- .all_slide_dedupe(cp)
        if (!has_outside_flag(pooled)) {
          warning("all_slide: ", pid, " has no readable polygon and no flag — skipping union")
          return(tibble::tibble())
        }
        return(region_ratios_area(pooled[!cell_outside(pooled), , drop = FALSE], 0, um_per_px) |>
                 dplyr::mutate(patient_id = pid, annotation = "union",
                               source = "flag", .before = 1))
      }
      # Pool the patient's region files, then de-duplicate: the same physical cell
      # appears in every region export of that slide when the producer exports the
      # whole slide per region, and cell_key_cols() is what identifies it.
      pooled <- .all_slide_dedupe(cp)
      poly_u <- sf::st_union(sf::st_geometry(polys))
      pts    <- sf::st_as_sf(cell_centroids_px(pooled, um_per_px),
                             coords = c("x", "y"), crs = sf::st_crs(polys))
      inside <- lengths(sf::st_within(pts, poly_u)) > 0
      return(region_ratios_area(pooled[inside, , drop = FALSE],
                                as.numeric(sum(sf::st_area(poly_u))), um_per_px) |>
               dplyr::mutate(patient_id = pid, annotation = "union",
                             source = "sf", .before = 1))
    }

    purrr::map_dfr(sort(unique(cp$annotation)), function(ann) {
      cells_r <- dplyr::filter(cp, annotation == ann)
      i <- match(paste(pid, ann), ann_key)
      if (is.na(i)) {
        # A region file with no polygon: fall back to the export's own flag if it
        # carries one, and say so in `source`. Never silently count it whole.
        if (!has_outside_flag(cells_r)) {
          warning("all_slide: ", pid, " ", ann, " has neither a polygon nor a flag — skipping")
          return(tibble::tibble())
        }
        return(region_ratios_area(cells_r[!cell_outside(cells_r), , drop = FALSE],
                                  0, um_per_px) |>
                 dplyr::mutate(patient_id = pid, annotation = ann,
                               source = "flag", .before = 1))
      }
      poly <- sf::st_geometry(annots)[i]
      pts  <- sf::st_as_sf(cell_centroids_px(cells_r, um_per_px),
                           coords = c("x", "y"), crs = sf::st_crs(annots))
      inside <- lengths(sf::st_within(pts, poly)) > 0
      region_ratios_area(cells_r[inside, , drop = FALSE],
                         as.numeric(sf::st_area(poly)), um_per_px) |>
        dplyr::mutate(patient_id = pid, annotation = ann, source = "sf", .before = 1)
    })
  })
}

# The same physical cell can appear in several of a patient's region files. Key on
# cell_key_cols() (cell_tables.R) when the export carries a real identity; fall back
# to the rounded centroid, which is stable across exports of the same segmentation.
.all_slide_dedupe <- function(cp) {
  keys <- intersect(cell_key_cols(cp), names(cp))
  if (length(keys)) return(dplyr::distinct(cp, dplyr::across(dplyr::all_of(keys)), .keep_all = TRUE))
  xy <- cell_centroids_px(cp, 1)
  cp[!duplicated(paste(round(xy$x, 2), round(xy$y, 2))), , drop = FALSE]
}

# The cell set for the whole-cohort readout: one row per physical cell per patient,
# pooled across that patient's region files. Regions overlap, files repeat cells, so
# this MUST dedupe — summing the region files would inflate every cohort-level count.
all_slide_union_cells <- function(cells) {
  if (nrow(cells) == 0) return(cells)
  purrr::map_dfr(unique(cells$patient_id), function(pid)
    .all_slide_dedupe(dplyr::filter(cells, patient_id == pid)))
}
