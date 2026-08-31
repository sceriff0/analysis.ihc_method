# =============================================================================
# arm_cells.R  —  ONE reader for all three phenotyping arms.
#
# Every arm publishes the same two things — cells and pathologist polygons —
# matched at the REGION level, so the pathologist's line and the cell export are
# paired 1:1 rather than one slide being scored against all of its polygons at
# once. What differs between arms is only WHERE the files sit and HOW a filename
# names its region, and both of those live in code/arms.R. This file knows about
# tiers and scopes; it never parses a path.
#
#   arm_cells(spec)              the region tier, one long table keyed
#                                (patient_id, annotation)
#   arm_union_tier_cells(spec)   the whole-slide tier, for the arms that ship one
#   arm_annotations(spec, tier)  the polygons of one tier
#   arm_metrics(...)             one metrics row per (patient, annotation), or per
#                                patient for scope = "union"
#   arm_cohort_cells(spec, ...)  one row per PHYSICAL cell per patient
#   arm_cells_in_annotation(...)  the same, plus a PER-CELL `in_annotation` flag,
#                                for a readout that re-cuts the denominator
#
# THE TWO TIERS ARE THE TWO SCOPES. An arm that ships an `_all` tier answers
# `union` from it directly — a real whole-slide export and a real union polygon,
# no inference. An arm without one derives the union the way this file always
# has: pool the region files, de-duplicate, dissolve the polygons.
#
# A PATIENT WITH NO REGION FILES FALLS BACK TO ITS UNION TIER, AS ANNOTATION_1.
# In massimo1, 10338 and 15897 have a whole-slide csv and a single `annotation_all`
# polygon but no `_selected` regions at all. Dropping them would shrink every
# union-scoped denominator without saying so; inventing regions for them would
# claim polygons nobody drew. So their one polygon serves as both the union row
# and ANNOTATION_1 — the same rule the flat legacy layout applied to a bare
# `<pid>.geojson`, applied here deliberately and recorded in `source`.
#
# NO ANNOTATION DIRECTORY MEANS EVERYTHING IS INSIDE — but only where the arm's
# spec says so (`bare_region_is = "whole_slide"`, which today is massimo2 alone).
# This is that layout's own convention, stated by the people who produced it: a
# patient with no `annotation/<pid>/` contributes its whole slide as one region.
# That is a DIFFERENT statement from "this patient's polygon is missing", which
# every other arm treats as a reason to drop the patient. Silently reinterpreting
# a missing polygon as "all in" everywhere would turn a data problem into a
# 100%-inside patient without anyone noticing, so it is per-arm and recorded in
# `source` as "whole_slide" rather than "sf".
#
# MEMBERSHIP INSIDE A REGION, in order of preference:
#   1. "sf"          the region's own geojson, point-in-polygon. Primary: it is the
#                    pathologist's line, and it is the only source that knows the
#                    region's AREA, so it is the only one that yields densities.
#   2. "flag"        the export's Out_of_annotation flag, when the region has a CSV
#                    but no readable polygon.
#   3. "whole_slide" every row counts, for a patient the arm's spec allows it for.
#
# Depends on arms.R (the registry), cell_tables.R (read_cell_csv and the
# accessors) and, at use time, on validation_helpers.R for slide_key() /
# region_ratios_area(). sf is a LAZY dependency, exactly as in
# validation_helpers.R: an arm with no readable polygon at all must still load on
# a machine without sf.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
  library(purrr)
  library(tibble)
})

source(here::here("code", "arms.R"))

# --- Reading one tier --------------------------------------------------------
# Every file of a tier, with the (patient, annotation) its name declares. Returns
# a 0-row frame rather than erroring when the directory is absent: an arm can
# legitimately lack a tier (massimo2 ships no `_all`), and a half-copied tree must
# fail at the report, naming the directory, not inside a loader.
.arm_tier_files <- function(tier, ext, what) {
  if (is.null(tier)) return(tibble::tibble())
  if (!dir.exists(tier$path)) {
    warning("arm: no ", what, " directory at ", tier$path)
    return(tibble::tibble())
  }
  files <- as.character(fs::dir_ls(tier$path, glob = paste0("*.", ext),
                                   recurse = TRUE, type = "file"))
  if (!length(files)) {
    warning("arm: no .", ext, " under ", tier$path)
    return(tibble::tibble())
  }
  keys <- lapply(files, arm_parse_name, pattern = tier$pattern)
  tibble::tibble(
    path       = files,
    patient_id = vapply(keys, `[[`, character(1), "patient"),
    annotation = vapply(keys, `[[`, character(1), "annotation"))
}

# --- Cells: the region tier --------------------------------------------------
# One long table, keyed by (patient_id, annotation). `annotation` is the region the
# FILE belongs to — which region a given cell falls in is decided later, by
# arm_metrics(), because only the polygon can say so.
#
# A file whose name carries no region is handled by the arm's own
# `bare_region_is` rule. Where that rule is NA the file is a contradiction in the
# layout, so it is dropped with a warning rather than pooled in as an extra region
# and double-counting every cell it holds.
# One tier's cells. Factored out of arm_cells() so the primary and fallback tiers
# are read by the SAME code — a fallback read a second way would be a second set of
# rules for the same question.
.arm_read_region_tier <- function(spec, tier, what) {
  meta <- .arm_tier_files(tier, "csv", what)
  if (nrow(meta) == 0) return(tibble::tibble())

  bare <- is.na(meta$annotation)
  if (any(bare)) {
    if (is.na(spec$bare_region_is)) {
      warning("arm ", spec$arm, ": ", sum(bare), " region csv(s) carry no region suffix (",
              paste(fs::path_file(meta$path[bare]), collapse = ", "),
              ") and this arm has no rule for one — dropping them")
      meta <- meta[!bare, , drop = FALSE]
    } else {
      # A patient with BOTH region files and a bare one is still a contradiction
      # even where bare files are legal: the regions are the finer statement, so
      # they win and the bare file is dropped rather than pooled on top of them.
      dup <- meta$patient_id[bare] %in% meta$patient_id[!bare]
      if (any(dup)) {
        warning("arm ", spec$arm, ": patient(s) ",
                paste(unique(meta$patient_id[bare][dup]), collapse = ", "),
                " have both region files and a bare whole-slide csv — dropping the ",
                "bare file, the regions win")
        meta <- meta[!(bare & meta$patient_id %in% meta$patient_id[!bare]), , drop = FALSE]
      }
      meta$annotation[is.na(meta$annotation)] <- spec$bare_region_is
    }
  }
  if (nrow(meta) == 0) return(tibble::tibble())

  # Which patients have a polygon at all, per the arm's region polygon tier. Held
  # per row so a report can tell an annotated patient from a whole-slide one
  # without re-reading the tree.
  poly_meta <- .arm_tier_files(spec$region_poly, "geojson", "region polygon")
  annotated <- if (nrow(poly_meta)) unique(slide_key(poly_meta$patient_id)) else character(0)

  purrr::pmap_dfr(meta, function(path, patient_id, annotation) {
    pid   <- slide_key(patient_id)
    cells <- read_cell_csv(path, patient_id = pid)
    if (nrow(cells) == 0) {
      warning("arm ", spec$arm, ": ", path, " has no cells — skipping")
      return(arm_empty_metrics())
    }
    has_ann <- pid %in% annotated
    # `ann` and `src` are settled BEFORE the mutate on purpose. mutate() evaluates its
    # arguments in order against the frame it is building, so a later argument that
    # names `annotation` sees the COLUMN just created — a length-nrow vector — not the
    # scalar this file meant. That silently made every whole-slide patient come out
    # `flag` instead of `whole_slide`, which is a different membership rule reported
    # under a plausible name.
    ann <- annotation
    src <- if (has_ann) "sf" else if (identical(ann, "whole_slide")) "whole_slide" else "flag"
    dplyr::mutate(cells,
                  arm                = spec$arm,
                  annotation         = ann,
                  classification     = tier$classification %||% "normal",
                  file_origin        = spec$arm,
                  has_annotation     = has_ann,
                  .membership_source = src)
  })
}

# The arm's region cells: its own tier, topped up from the fallback for the patients
# it does not cover.
#
# TOP-UP IS PER PATIENT, NOT PER FILE. A patient the arm re-ran is taken entirely
# from the arm's own tier; a patient it did not is taken entirely from the fallback.
# Mixing the two within one patient would put re-classified and ordinary cells in the
# same denominator, and no column could then say what that patient's number means.
arm_cells <- function(spec) {
  primary <- .arm_read_region_tier(spec, spec$region_csv, "region csv")
  if (!arm_has_fallback(spec)) return(primary)

  fb <- .arm_read_region_tier(spec, spec$region_csv_fallback, "fallback region csv")
  if (nrow(fb) == 0) return(primary)
  have <- if (nrow(primary)) unique(primary$patient_id) else character(0)
  fb   <- dplyr::filter(fb, !patient_id %in% have)
  if (nrow(fb) == 0) return(primary)

  message(sprintf("arm %s: %d patient(s) not re-run (%s) load their %s cells from %s",
                  spec$arm, dplyr::n_distinct(fb$patient_id),
                  paste(sort(unique(fb$patient_id)), collapse = ", "),
                  spec$region_csv_fallback$classification, spec$region_csv_fallback$dir))
  if (nrow(primary) == 0) fb else dplyr::bind_rows(primary, fb)
}

# --- Cells: the whole-slide tier ---------------------------------------------
# The arms that publish a real whole-slide export read it here rather than
# reconstructing it. One row per cell per patient by construction, so nothing is
# de-duplicated — that is the whole point of having the tier. Empty for an arm
# with no `_all`, whose caller then derives the set instead.
arm_union_tier_cells <- function(spec) {
  if (!arm_has_union_tier(spec)) return(tibble::tibble())
  meta <- .arm_tier_files(spec$union_csv, "csv", "whole-slide csv")
  if (nrow(meta) == 0) return(tibble::tibble())

  # An arm with a fallback re-ran only some patients, and its whole-slide tier is the
  # FALLBACK's export. Serving a re-run patient from it would put ordinary cells in
  # that patient's union row while its per-region rows hold re-classified ones — the
  # same patient reading two different ways inside one arm, with no column to say so.
  if (arm_has_fallback(spec)) {
    rerun <- .arm_read_region_tier(spec, spec$region_csv, "region csv")
    if (nrow(rerun)) {
      drop <- slide_key(meta$patient_id) %in% unique(rerun$patient_id)
      if (any(drop))
        message(sprintf("arm %s: %s re-run by this arm, so their union comes from the ",
                        spec$arm, paste(sort(unique(slide_key(meta$patient_id)[drop])),
                                        collapse = ", ")),
                "pooled region files, not the whole-slide export")
      meta <- meta[!drop, , drop = FALSE]
    }
  }
  if (nrow(meta) == 0) return(tibble::tibble())

  purrr::pmap_dfr(meta, function(path, patient_id, annotation) {
    pid   <- slide_key(patient_id)
    cells <- read_cell_csv(path, patient_id = pid)
    if (nrow(cells) == 0) {
      warning("arm ", spec$arm, ": ", path, " has no cells — skipping")
      return(arm_empty_metrics())
    }
    dplyr::mutate(cells, arm = spec$arm, annotation = "union",
                  classification = spec$union_csv$classification %||% "normal",
                  file_origin = spec$arm, has_annotation = TRUE,
                  .membership_source = "sf")
  })
}

# --- Polygons ----------------------------------------------------------------
# One tier's polygons, as an sf frame of (patient_id, annotation, geometry).
#
# `tier = "region"` labels each polygon ANNOTATION_<k>. `tier = "union"` labels
# every polygon "union" — a bare `<pid>.geojson` in an `_all` tree is the
# patient's dissolved boundary, NOT its region 1, and reading it as ANNOTATION_1
# would score the whole slide against a region the pathologist never drew.
arm_annotations <- function(spec, tier = c("region", "union"), patient_ids = NULL) {
  tier <- match.arg(tier)
  key  <- if (tier == "region") "region_poly" else "union_poly"
  if (is.null(spec[[key]])) return(NULL)
  .require_sf("arm_annotations()")

  meta <- .arm_tier_files(spec[[key]], "geojson", paste(tier, "polygon"))
  if (nrow(meta) == 0) return(NULL)
  meta$patient_id <- slide_key(meta$patient_id)
  if (tier == "union") meta$annotation <- "union"

  drop <- is.na(meta$annotation)
  if (any(drop)) {
    warning("arm ", spec$arm, ": ", sum(drop), " ", tier,
            " polygon(s) carry no region suffix — skipping ",
            paste(fs::path_file(meta$path[drop]), collapse = ", "))
    meta <- meta[!drop, , drop = FALSE]
  }
  if (!is.null(patient_ids))
    meta <- meta[meta$patient_id %in% slide_key(patient_ids), , drop = FALSE]
  if (nrow(meta) == 0) return(NULL)

  polys <- purrr::pmap(meta, function(path, patient_id, annotation) {
    geom <- tryCatch(read_polygon_geojson(path), error = function(e) {
      warning("arm ", spec$arm, ": unreadable polygon, skipping ", path,
              " — ", conditionMessage(e))
      NULL
    })
    if (is.null(geom)) return(NULL)
    sf::st_sf(patient_id = patient_id, annotation = annotation, geometry = geom)
  })
  polys <- purrr::compact(polys)
  if (!length(polys)) return(NULL)
  do.call(rbind, polys)
}

# --- The union-tier fallback for a patient with no regions -------------------
# 10338 and 15897 in massimo1: a whole-slide csv and one `annotation_all` polygon,
# no `_selected` files at all. Promote both to ANNOTATION_1 so the patient appears
# in the per-annotation frame instead of vanishing from it. Returns its arguments
# unchanged for an arm with no union tier, so the rule costs nothing where it does
# not apply.
#
# The promotion is deliberately NOT a general "a missing region means the whole
# slide": it fires only for a patient the arm has NO region file for, and only
# when that patient does have a union-tier polygon to promote. A patient with two
# of three regions exported keeps two regions and is reported as such.
arm_promote_unregioned <- function(spec, cells, polys, union_cells, union_polys) {
  if (!arm_has_union_tier(spec) || nrow(union_cells) == 0)
    return(list(cells = cells, polys = polys))

  regioned <- if (nrow(cells)) unique(cells$patient_id) else character(0)
  orphans  <- setdiff(unique(union_cells$patient_id), regioned)
  if (!length(orphans)) return(list(cells = cells, polys = polys))

  message(sprintf("arm %s: %s ha%s a whole-slide export but no region files — ",
                  spec$arm, paste(orphans, collapse = ", "),
                  if (length(orphans) == 1) "s" else "ve"),
          "promoting the union polygon to ANNOTATION_1 so they keep a per-region row")

  extra_cells <- union_cells |>
    dplyr::filter(patient_id %in% orphans) |>
    dplyr::mutate(annotation = "ANNOTATION_1")
  cells <- if (nrow(cells)) dplyr::bind_rows(cells, extra_cells) else extra_cells

  if (!is.null(union_polys) && nrow(union_polys)) {
    extra_polys <- union_polys[union_polys$patient_id %in% orphans, , drop = FALSE]
    if (nrow(extra_polys)) {
      extra_polys$annotation <- "ANNOTATION_1"
      polys <- if (is.null(polys) || nrow(polys) == 0) extra_polys
               else rbind(polys, extra_polys)
    }
  }
  list(cells = cells, polys = polys)
}

# --- Provenance --------------------------------------------------------------
# What was actually read, one row per (patient, region), for the report to print.
# `n_inside` is left to the metrics frame — this table answers "did my data
# arrive", which has to be answerable even when the geometry step later fails.
arm_inventory <- function(cells) {
  if (nrow(cells) == 0) return(tibble::tibble())
  cells |>
    dplyr::group_by(arm, patient_id, annotation, classification) |>
    dplyr::summarise(n_cells        = dplyr::n(),
                     has_annotation = any(has_annotation),
                     membership     = paste(sort(unique(.membership_source)), collapse = "/"),
                     .groups = "drop") |>
    dplyr::arrange(arm, patient_id, annotation)
}

# THE EMPTY METRICS FRAME IS TYPED, NOT BARE.
#
# `tibble::tibble()` has zero COLUMNS, so a downstream
# `select(patient_id, annotation, source, n_inside)` does not return nothing — it
# ERRORS with "Column `patient_id` doesn't exist", three chunks after the actual
# problem and naming none of it. An arm with no data on disk is a routine state
# (a fresh clone, a tree not symlinked yet), so it has to flow through the same
# joins and selects as a full one and simply produce empty output.
#
# The schema is taken from region_ratios_area() itself rather than written out, so
# it cannot drift from the real thing: region_ratios() already handles a zero-row
# cell table, so calling it on one yields the exact column set and types, and the
# [0, ] drops the placeholder row.
arm_empty_metrics <- function() {
  rr <- region_ratios_area(tibble::tibble(), 0, 1)
  dplyr::mutate(rr,
                patient_id = NA_character_,
                annotation = NA_character_,
                source     = NA_character_,
                .before    = 1)[0, , drop = FALSE]
}

# --- Metrics -----------------------------------------------------------------
# One metrics row per (patient_id, annotation), schema-identical to
# ihc_annotation_metrics() so membership_data() can hand any arm to the same report.
#
# Per region: the cells of THAT region's file, cut down to the ones inside THAT
# region's polygon. Two files never contribute to the same region, so nothing is
# double-counted even where the polygons overlap.
#
# scope = "union" prefers the arm's own whole-slide tier — a real export against a
# real union polygon, with no de-duplication in the path at all. Only an arm
# without that tier dissolves its region polygons and pools its region files,
# which is why the union is computed from geometry rather than by summing the
# per-region rows: a cell in the overlap of two regions must be counted once.
arm_metrics <- function(spec, cells, annots, scope = c("per_annotation", "union"),
                        um_per_px = 0.325, union_cells = NULL, union_polys = NULL) {
  scope <- match.arg(scope)

  # -- the dedicated whole-slide tier, where the arm has one -------------------
  # It may cover only SOME patients: an arm with a fallback serves its re-run
  # patients from their pooled region files instead (see arm_union_tier_cells()).
  # So this collects the rows it can and lets the rest fall through to the pooled
  # path below, rather than returning and dropping them.
  tier_rows <- arm_empty_metrics()
  if (scope == "union" && arm_has_union_tier(spec) &&
      !is.null(union_cells) && nrow(union_cells) > 0) {
    tier_rows <- purrr::map_dfr(unique(union_cells$patient_id), function(pid) {
      cp <- dplyr::filter(union_cells, patient_id == pid)
      pu <- if (!is.null(union_polys))
              union_polys[union_polys$patient_id == pid, , drop = FALSE] else NULL
      if (is.null(pu) || nrow(pu) == 0) {
        # A whole-slide csv with no union polygon. Fall back to the export's own
        # flag and SAY so — never count it whole, because unlike massimo2's
        # convention nothing here states that a missing polygon means "all in".
        if (!has_outside_flag(cp)) {
          warning("arm ", spec$arm, ": ", pid,
                  " has a whole-slide csv but no union polygon and no flag — skipping union")
          return(arm_empty_metrics())
        }
        return(region_ratios_area(cp[!cell_outside(cp), , drop = FALSE], 0, um_per_px) |>
                 dplyr::mutate(patient_id = pid, annotation = "union",
                               source = "flag", .before = 1))
      }
      poly_u <- sf::st_union(sf::st_geometry(pu))
      pts    <- sf::st_as_sf(cell_centroids_px(cp, um_per_px),
                             coords = c("x", "y"), crs = sf::st_crs(pu))
      inside <- lengths(sf::st_within(pts, poly_u)) > 0
      region_ratios_area(cp[inside, , drop = FALSE],
                         as.numeric(sum(sf::st_area(poly_u))), um_per_px) |>
        dplyr::mutate(patient_id = pid, annotation = "union", source = "sf", .before = 1)
    })
    # Whatever the tier covered must not also be pooled, or the patient gets two
    # union rows and every cohort mean counts it twice.
    if (nrow(cells))
      cells <- dplyr::filter(cells, !patient_id %in% unique(union_cells$patient_id))
  }

  if (nrow(cells) == 0)
    return(if (nrow(tier_rows)) tier_rows else arm_empty_metrics())
  has_poly <- !is.null(annots) && nrow(annots) > 0
  ann_key  <- if (has_poly) paste(slide_key(annots$patient_id), annots$annotation) else character(0)

  rest <- purrr::map_dfr(unique(cells$patient_id), function(pid) {
    cp <- dplyr::filter(cells, patient_id == pid)

    # -- whole-slide patient: no polygon by the arm's own convention, so every cell
    # -- counts and there is no area to divide by. NA densities, as the flag modes.
    if (!any(cp$has_annotation)) {
      return(region_ratios_area(cp, 0, um_per_px) |>
               dplyr::mutate(patient_id = pid, annotation = "whole_slide",
                             source = "whole_slide", .before = 1))
    }

    if (scope == "union") {
      polys <- if (has_poly) dplyr::filter(annots, slide_key(patient_id) == pid) else NULL
      # An annotated patient whose polygons could not be read (no sf, unreadable
      # geojson) must still produce a union row: dropping it here would shrink the
      # cohort in every union-scoped panel without saying so. Fall back to the
      # export's own flag on the pooled slide, and label the row `flag` so the
      # different membership rule is visible in the frame rather than inferred.
      pooled <- .arm_dedupe(cp)
      if (is.null(polys) || nrow(polys) == 0) {
        if (!has_outside_flag(pooled)) {
          warning("arm ", spec$arm, ": ", pid,
                  " has no readable polygon and no flag — skipping union")
          return(arm_empty_metrics())
        }
        return(region_ratios_area(pooled[!cell_outside(pooled), , drop = FALSE], 0, um_per_px) |>
                 dplyr::mutate(patient_id = pid, annotation = "union",
                               source = "flag", .before = 1))
      }
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
          warning("arm ", spec$arm, ": ", pid, " ", ann,
                  " has neither a polygon nor a flag — skipping")
          return(arm_empty_metrics())
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

  if (nrow(tier_rows) == 0) return(rest)
  dplyr::arrange(dplyr::bind_rows(tier_rows, rest), patient_id, annotation)
}

# The same physical cell can appear in several of a patient's region files. Key on
# cell_key_cols() (cell_tables.R) when the export carries a real identity; fall back
# to the rounded centroid, which is stable across exports of the same segmentation.
.arm_dedupe <- function(cp) {
  keys <- intersect(cell_key_cols(cp), names(cp))
  if (length(keys)) return(dplyr::distinct(cp, dplyr::across(dplyr::all_of(keys)), .keep_all = TRUE))
  xy <- cell_centroids_px(cp, 1)
  cp[!duplicated(paste(round(xy$x, 2), round(xy$y, 2))), , drop = FALSE]
}

# WHAT THE REGION FILES ACTUALLY CONTAIN — report it, do not assume it.
#
# There are two possibilities and the layout alone cannot tell them apart:
#
#   (a) each region csv holds the WHOLE SLIDE, with Out_of_annotation computed for
#       that region. Then a patient's region files repeat the same cells, and pooling
#       them triples that patient's cell count.
#   (b) each holds ONLY the cells inside its region. Then the files are disjoint and
#       pooling them is exactly right.
#
# The presence of an Out_of_annotation column hints at (a) — you only need a flag if
# the file contains cells the region excludes — but a producer can write it either way,
# and guessing wrong silently changes every cohort-level denominator on the site while
# leaving every number plausible.
#
# So the code is correct under BOTH: the per-region metrics intersect each file with
# its own polygon (a no-op under (b)), and the cohort cell set de-duplicates (a no-op
# under (b)). This function reports which one the data on disk actually is, so the
# answer arrives as a printed table on the first knit rather than as an assumption
# nobody revisits.
arm_overlap_report <- function(cells) {
  if (nrow(cells) == 0) return(tibble::tibble())
  out <- purrr::map_dfr(unique(cells$patient_id), function(pid) {
    cp <- dplyr::filter(cells, patient_id == pid)
    n_files <- dplyr::n_distinct(cp$annotation)
    n_rows  <- nrow(cp)
    n_uniq  <- nrow(.arm_dedupe(cp))
    tibble::tibble(patient_id = pid, n_region_files = n_files,
                   n_rows = n_rows, n_unique_cells = n_uniq,
                   rows_per_cell = round(n_rows / max(n_uniq, 1), 2))
  })
  # One region file cannot repeat itself, so a single-file patient is uninformative.
  multi <- dplyr::filter(out, n_region_files > 1)
  verdict <- if (nrow(multi) == 0) "no multi-region patient — cannot tell"
    else if (all(multi$rows_per_cell > 1.5)) "whole-slide exports (files repeat cells)"
    else if (all(multi$rows_per_cell < 1.1)) "region-restricted exports (files are disjoint)"
    else "MIXED — check the producer; some patients repeat cells and others do not"
  attr(out, "verdict") <- verdict
  out
}

# --- The cohort cell set -----------------------------------------------------
# One row per physical cell per patient — the denominator behind every whole-slide
# readout (marker_qc, molecular_hot_cold).
#
# An arm with an `_all` tier reads it straight off disk. An arm without one pools
# its region files and de-duplicates, which is the only thing it can do and is
# correct under both export shapes above.
arm_cohort_cells <- function(spec, cells = NULL, union_cells = NULL) {
  if (is.null(cells)) cells <- arm_cells(spec)
  if (arm_has_union_tier(spec) && is.null(union_cells))
    union_cells <- arm_union_tier_cells(spec)
  if (is.null(union_cells)) union_cells <- tibble::tibble()

  # PER PATIENT, the best source it has. The whole-slide tier where one exists for
  # that patient — a real export, no de-duplication in the path at all — and the
  # pooled, de-duplicated region files otherwise. An arm with a fallback has both
  # kinds at once: its re-run patients pool, the rest read their whole-slide csv.
  covered <- if (nrow(union_cells)) unique(union_cells$patient_id) else character(0)
  if (nrow(cells) == 0) return(union_cells)
  rest <- dplyr::filter(cells, !patient_id %in% covered)
  pooled <- if (nrow(rest))
    purrr::map_dfr(unique(rest$patient_id), function(pid)
      .arm_dedupe(dplyr::filter(rest, patient_id == pid)))
  else tibble::tibble()

  if (nrow(union_cells) == 0) return(pooled)
  if (nrow(pooled) == 0) return(union_cells)
  dplyr::bind_rows(union_cells, pooled)
}

# --- Per-cell membership, where a metrics ROW is not enough -------------------
# arm_metrics() already asks the polygon which cells are inside — and then throws
# the per-cell verdict away, because a metrics row is an aggregate. A readout that
# RE-CUTS the denominator cannot be built from that aggregate after the fact: the
# bulk-RNA comparison needs a per-marker positive fraction over "tumour cells inside
# the annotation", and no combination of tumor_over_inside and cd45_over_inside
# yields it. So the verdict is exposed here rather than recomputed at the call site,
# which is how validation_helpers.R ended up with a second copy of the parser.
#
# Returns `cells` with two added columns:
#   in_annotation          logical, NA where no rule could decide
#   .in_annotation_source  which rule did decide, one of:
#     "sf"           point-in-polygon against the arm's UNION polygon — massimo1's
#                    `annotation_all`, massimo2's dissolved region polygons. The
#                    SAME polygon arm_metrics(scope = "union") scores against, so a
#                    fraction computed here and a union metrics row agree by
#                    construction rather than by coincidence.
#     "flag"         the export's own Out_of_annotation column, for a patient with
#                    cells but no readable polygon.
#     "whole_slide"  every cell counts — ONLY for an arm whose `bare_region_is`
#                    states that an unannotated patient is annotated whole
#                    (massimo2's 24086). No other arm may take this branch, and it
#                    OUTRANKS the flag, because arm_metrics() never consults the flag
#                    for such a patient either.
#     NA             no polygon, no flag, no convention. in_annotation is NA and the
#                    patient DROPS OUT of the restricted readout, rather than being
#                    counted whole under a rule nobody wrote down — a silently
#                    whole-counted patient is indistinguishable from a correct one.
#
# THE UNION POLYGON IS PREFERRED OVER THE DISSOLVED REGIONS, and for massimo1 that
# matters: annotation_all is a boundary the pathologist actually drew, while
# dissolving the three `_selected` regions produces a shape nobody drew. The regions
# nest inside the union (annotation_k ⊆ annotation_all), so preferring the union is
# the WIDER and less presumptuous of the two. An arm with no `_all` tier has only the
# dissolved regions and uses them.
#
# sf is a LAZY dependency here as everywhere: a machine without it falls through to
# the flag rather than failing to load the page.
arm_cells_in_annotation <- function(spec, cells, um_per_px = 0.325) {
  if (is.null(cells) || nrow(cells) == 0) return(cells)
  pids <- unique(cells$patient_id)

  polys <- NULL
  if (requireNamespace("sf", quietly = TRUE)) {
    polys <- arm_annotations(spec, "union", patient_ids = pids)
    if (is.null(polys) || nrow(polys) == 0)
      polys <- arm_annotations(spec, "region", patient_ids = pids)
  } else {
    warning("arm_cells_in_annotation(): sf is not installed, so no polygon can be ",
            "read — falling back to the export's Out_of_annotation flag")
  }
  poly_ids <- if (!is.null(polys) && nrow(polys)) unique(polys$patient_id) else character(0)

  purrr::map_dfr(pids, function(pid) {
    cp <- dplyr::filter(cells, patient_id == pid)

    if (pid %in% poly_ids) {
      pp     <- polys[polys$patient_id == pid, , drop = FALSE]
      poly_u <- sf::st_union(sf::st_geometry(pp))
      xy     <- cell_centroids_px(cp, um_per_px)
      ok     <- is.finite(xy$x) & is.finite(xy$y)
      # A cell with no usable centroid cannot be placed, so it is NA rather than
      # outside: "outside" is a geometric claim and there is no geometry to make it.
      if (any(ok)) {
        inside     <- rep(NA, nrow(cp))
        pts        <- sf::st_as_sf(xy[ok, , drop = FALSE], coords = c("x", "y"),
                                   crs = sf::st_crs(pp))
        inside[ok] <- lengths(sf::st_within(pts, poly_u)) > 0
        return(dplyr::mutate(cp, in_annotation = inside, .in_annotation_source = "sf"))
      }
    }

    # THE CONVENTION IS CHECKED BEFORE THE FLAG, and the order matters.
    #
    # arm_metrics() tests `!any(cp$has_annotation)` FIRST and returns every cell for
    # such a patient, never consulting the flag. So a patient with no annotation
    # directory at all in an arm whose spec allows it (massimo2's 24086) must be
    # counted whole here too — taking its flag instead would keep ~a quarter fewer
    # cells than its own union metrics row, and the two numbers would disagree on
    # the same page with nothing to say which was which.
    #
    # `has_annotation` comes off the LISTING of the polygon tier, not off a
    # successful parse, so a patient whose geojson exists but is unreadable is still
    # "annotated" and correctly falls through to the flag below — which is again what
    # arm_metrics() does.
    unannotated <- !("has_annotation" %in% names(cp)) || !any(cp$has_annotation %in% TRUE)
    if (unannotated && identical(spec$bare_region_is, "whole_slide"))
      return(dplyr::mutate(cp, in_annotation = TRUE,
                           .in_annotation_source = "whole_slide"))

    if (has_outside_flag(cp))
      return(dplyr::mutate(cp, in_annotation = !cell_outside(cp),
                           .in_annotation_source = "flag"))

    warning("arm ", spec$arm, ": ", pid, " has no readable union polygon and no ",
            "Out_of_annotation flag — its in-annotation membership is NA and it is ",
            "excluded from annotation-restricted readouts")
    dplyr::mutate(cp, in_annotation = NA, .in_annotation_source = NA_character_)
  })
}

# Provenance for the above, one row per patient: how many cells the arm holds, how
# many of them the annotation keeps, and which rule decided. Printed beside every
# annotation-restricted figure, because a restriction that quietly kept 100% (a
# `whole_slide` row) and one that kept 40% (an `sf` row) produce the same-looking
# panel and mean entirely different things.
arm_in_annotation_inventory <- function(cells) {
  if (is.null(cells) || nrow(cells) == 0 || !"in_annotation" %in% names(cells))
    return(tibble::tibble())
  # Settled BEFORE the mutate, as everywhere in this file: an expression inside
  # mutate() that names the frame it is building is a trap waiting for the day a
  # column called `phenotype_clean` is added upstream.
  is_tumor <- cell_lineage(cell_phenotype(cells)) %in% "Tumor"
  cells |>
    dplyr::mutate(.tumor = is_tumor) |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(
      membership       = paste(sort(unique(stats::na.omit(.in_annotation_source))),
                               collapse = "/"),
      n_cells          = dplyr::n(),
      n_in_annotation  = sum(in_annotation %in% TRUE),
      n_tumor          = sum(.tumor),
      n_tumor_in_annotation = sum(.tumor & in_annotation %in% TRUE),
      pct_in_annotation = round(100 * sum(in_annotation %in% TRUE) / dplyr::n(), 1),
      .groups = "drop") |>
    dplyr::arrange(patient_id)
}

# CHECK THE PROCEDURE AGAINST THE GROUND TRUTH, where an arm supplies both.
#
# massimo1 is the only arm that publishes a real whole-slide export AND a set of
# region files, so it is the only one where the de-duplication the other two arms
# are FORCED to use can be verified rather than trusted. Arms 2 and 3 have no `_all`
# tier and must reconstruct their cohort set by pooling and de-duplicating; if that
# procedure is sound, running it on arm 1's region files should reproduce arm 1's
# own whole-slide export.
#
# A per-patient disagreement is not automatically an error — the `_selected` regions
# need not cover the whole slide, so the deduped union is expected to be a SUBSET —
# but the direction and size of the gap is the readout. Deduped > whole-slide is
# the one that cannot be explained by coverage and means the key is wrong.
arm_wholeslide_reconciliation <- function(spec, cells = NULL, union_cells = NULL) {
  if (!arm_has_union_tier(spec)) return(tibble::tibble())
  if (is.null(cells))       cells       <- arm_cells(spec)
  if (is.null(union_cells)) union_cells <- arm_union_tier_cells(spec)
  if (nrow(cells) == 0 || nrow(union_cells) == 0) return(tibble::tibble())

  pids <- sort(union(unique(cells$patient_id), unique(union_cells$patient_id)))
  out <- purrr::map_dfr(pids, function(pid) {
    cp <- dplyr::filter(cells, patient_id == pid)
    up <- dplyr::filter(union_cells, patient_id == pid)
    n_dedup <- if (nrow(cp)) nrow(.arm_dedupe(cp)) else 0L
    tibble::tibble(patient_id = pid,
                   n_region_files = dplyr::n_distinct(cp$annotation),
                   n_wholeslide   = nrow(up),
                   n_dedup_union  = n_dedup,
                   pct_of_slide   = if (nrow(up)) round(100 * n_dedup / nrow(up), 1) else NA_real_)
  })
  over <- dplyr::filter(out, n_wholeslide > 0, n_dedup_union > n_wholeslide)
  attr(out, "verdict") <- if (nrow(over))
      paste0("DE-DUPLICATION IS UNDER-MERGING for ", paste(over$patient_id, collapse = ", "),
             " — the pooled regions yield MORE cells than the whole-slide export, so ",
             "cell_key_cols() is not identifying cells. Every arm without an `_all` ",
             "tier inflates its denominators by the same factor.")
    else "de-duplicated regions are a subset of the whole-slide export, as expected"
  out
}
