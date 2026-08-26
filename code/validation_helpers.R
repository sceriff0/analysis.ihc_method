# =============================================================================
# validation_helpers.R  —  shared logic for the IHC-validation reports
#   (clinical_flowpath.Rmd, clinical_mirage.Rmd, clinical_membership_qc.Rmd,
#    molecular_hot_cold.Rmd, marker_qc.Rmd).
#
# The single-cell IHC table is the "measurement under test". These helpers derive
# the quantities each report validates against an independent reference
# (pathologist / bulk-RNA deconvolution / bulk-RNA marker genes).
#
# This file owns the DERIVED QUANTITIES (region ratios, composition, marker and
# lineage fractions, invasive-margin metrics, agreement statistics). It does not
# own the cell-table schema: every column read off an export goes through
# code/cell_tables.R, so a FlowPath csv, a mirage `*_phenotypes.csv` and a mirage
# join_flowpath cohort table all work here unchanged.
#
# Geometric membership: sf point-in-polygon on the raw pathologist geojson, cells
# mapped into the geojson pixel frame by cell_centroids_px() (fixed 0.325 um/px for
# a micron export). The exporter's own out-of-annotation flag is the fallback for
# patients with no geojson. code/membership.R builds on both.
#
# renv: needs sf, jsonlite (+ the dplyr/tidyr/readr/stringr/purrr/tibble/ggplot2
#   set, here, fs — see .need below).
#   renv::install(c("sf", "jsonlite")); renv::snapshot()
# =============================================================================

# Declared package by package rather than via the tidyverse metapackage, so this
# file can be sourced (and unit-tested) without attaching the whole of it. The Rmds
# now declare the same packages one by one in their own setup chunks, so the
# metapackage is no longer attached anywhere in the project.
.need <- c("dplyr", "tidyr", "readr", "stringr", "purrr", "tibble", "ggplot2",
           "here", "fs")
.missing <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing))
  stop("validation_helpers.R needs: ", paste(.missing, collapse = ", "),
       "\n  renv::install(c(", paste(sprintf('"%s"', .missing), collapse = ", "), "))",
       call. = FALSE)
suppressPackageStartupMessages(lapply(.need, library, character.only = TRUE))

# sf is a LAZY requirement, checked by the geometry functions that use it (every
# call below is sf::-namespaced, so the namespace suffices — no attach). Only the
# polygon paths need it: the flag-membership reports read a precomputed in/out
# column and no geojson at all, yet could not previously source on a machine
# without sf installed, for a dependency they never call.
.require_sf <- function(what) {
  if (!requireNamespace("sf", quietly = TRUE))
    stop(what, " needs the sf package (geojson polygons).",
         "\n  renv::install(\"sf\")", call. = FALSE)
}

# --- House plot style -------------------------------------------------------
# One definition for the whole project, in code/plot_theme.R: it supplies the
# `oi` palette, `hotcold_cols()`/`hotcold_order()`, `theme_paper()` and the
# scale shorthands, AND applies the theme + geom defaults on source(). Sourcing
# it here means every analysis that loads these helpers is styled automatically;
# nothing below needs to call theme_set().
source(here::here("code", "plot_theme.R"))

# --- Cell-table schema adapter ----------------------------------------------
# Every column the analyses read off a single-cell export goes through
# code/cell_tables.R: clean_phenotype(), cell_outside()/is_outside(), is_pos()/
# marker_pos(), cell_centroids_px(), cell_key_cols(). That file, not this one,
# is where the FlowPath / mirage-phenotypes / mirage-cohort spellings are
# reconciled — see its header for the three schemas.
source(here::here("code", "cell_tables.R"))
source(here::here("code", "arms.R"))    # ARM_MODES / arm_spec() / arm_parse_name()

# --- ID handling ------------------------------------------------------------
# All datasets share one patient/slide ID but differ in punctuation (clinical
# CRF uses dots, neoplastic uses bare digits, IHC comes from a filename). Strip
# non-alphanumerics and upper-case so the shared key joins; leading zeros are
# kept (they are significant, e.g. "052").
norm_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

# Normalise a SLIDE-space id (IHC patient_id, neoplastic SAMPLE, clinical
# `ID PATIENT`). These are numeric slide codes that sometimes carry an alpha
# prefix (e.g. "EPM - 052"); keep only the digits so "EPM - 052" -> "052" matches
# the IHC "052". Falls back to norm_id() for any purely non-numeric id.
norm_slide_id <- function(x) {
  x <- as.character(x)
  d <- gsub("[^0-9]", "", x)
  ifelse(d == "" | is.na(x), norm_id(x), d)
}

# Known clinical `ID PATIENT` data-entry errors: the value on the LEFT is what the
# clinical sheet reads, the RIGHT is the true slide/neoplastic id it should match.
# Default handles the verified transposition (clinical "15879" == slide "15897").
SLIDE_ID_FIXES <- c("15879" = "15897")

# THE canonical slide-space join key: normalise + apply the typo fixes. Every
# join on a slide id (IHC, neoplastic, clinical `ID PATIENT`) must use this so no
# call site silently omits the correction and drops a patient.
slide_key <- function(x, fixes = SLIDE_ID_FIXES) {
  s   <- norm_slide_id(x)
  hit <- s %in% names(fixes)
  s[hit] <- unname(fixes[s[hit]])
  s
}

# Crosswalk between the two id systems, read from clinical_data:
#   slide_id (norm_slide_id of `ID PATIENT`)  <->  crf_id (norm_id of `ID CRF PRESERVE`)
# `ID CRF PRESERVE` also keys counts_data$Sample, so this bridges IHC/neoplastic
# (slide space) to bulk RNA (CRF space). `slide_id_fixes` corrects known clinical
# data-entry errors: the default handles the verified transposition where clinical
# reads "15879" but the slide/neoplastic id is "15897". Extend or clear as needed.
id_crosswalk <- function(clinical_data,
                         patient_col = "ID PATIENT",
                         crf_col     = "ID CRF PRESERVE",
                         slide_id_fixes = SLIDE_ID_FIXES) {
  slide <- slide_key(clinical_data[[patient_col]], slide_id_fixes)
  crf   <- norm_id(clinical_data[[crf_col]])
  keep  <- slide != "" & crf != "" & !is.na(slide) & !is.na(crf)
  unique(tibble::tibble(slide_id = slide[keep], crf_id = crf[keep]))
}

# Named vector  crf_id -> slide_id  built DIRECTLY from clinical_data (base R: no
# tibble, no join). Map bulk-RNA (CRF space) sample ids to IHC slide ids with
# `slide[crf_id]`. Prefer this over id_crosswalk() + join in the RNA reports —
# dplyr join column-propagation is unreliable on this project's dplyr build.
crf_to_slide_map <- function(clinical_data,
                             patient_col = "ID PATIENT",
                             crf_col     = "ID CRF PRESERVE",
                             slide_id_fixes = SLIDE_ID_FIXES) {
  slide <- slide_key(clinical_data[[patient_col]], slide_id_fixes)
  crf   <- norm_id(clinical_data[[crf_col]])
  keep  <- !is.na(slide) & !is.na(crf) & slide != "" & crf != ""
  stats::setNames(slide[keep], crf[keep])
}

# Arcsine square-root ("angular") transform for proportions, the classic variance-
# stabiliser for data bounded in [0, 1]: it de-compresses the crowded 0/1 tails so
# a fraction's variance stops shrinking near the bounds. Values are clamped to
# [0, 1] first (guards tiny FP overshoot); anything outside — e.g. the unbounded
# tumor_over_cd45 ratio — becomes NA, since the transform is only meaningful for
# true proportions. Range of the result is [0, pi/2].
asin_sqrt <- function(p) {
  p <- ifelse(is.finite(p) & p >= 0 & p <= 1, p, NA_real_)
  asin(sqrt(p))
}

# --- Phenotype / cell-type vocabulary ---------------------------------------
# The label -> lineage table and cell_lineage() live in code/cell_tables.R, with the
# rest of the cross-tool vocabulary: reconciling FlowPath's and mirage's names is a
# reading concern, not an analysis one, and keeping it there lets it be tested
# without sf. What stays here is the analysis choice made on top of it.

# Lineages with a clean deconvolution counterpart, used for the method comparison.
comparable_lineages <- c("CD8T", "CD4T", "Treg", "NK")

# IHC protein marker -> canonical gene symbol(s). Multi-gene markers are averaged
# downstream. `wrongL1CAM` and DAPI are intentionally excluded.
marker_gene_map <- tibble::tribble(
  ~marker,     ~genes,
  "CD45",      "PTPRC",
  "CD3",       "CD3D,CD3E,CD3G",
  "CD8",       "CD8A,CD8B",
  "CD4",       "CD4",
  "GZMB",      "GZMB",
  "FOXP3",     "FOXP3",
  "CD14",      "CD14",
  "CD163",     "CD163",
  "CD56",      "NCAM1",
  "SMA",       "ACTA2",
  "PANCK",     "EPCAM,KRT5,KRT8,KRT14,KRT18,KRT19",
  "VIMENTIN",  "VIM",
  "CD74",      "CD74",
  "L1CAM",     "L1CAM",
  "P53",       "TP53",
  "PDL1",      "CD274",
  "PD1",       "PDCD1",
  "ARID1A",    "ARID1A",
  "FSP1",      "S100A4"
)

# Markers carrying a `<marker>_sign` column in ihc_data (wrongL1CAM excluded).
ihc_markers <- marker_gene_map$marker

# --- GeoJSON reading (sf) ----------------------------------------------------
# Parse a QuPath-style annotation geojson into an sf polygon. Handles a bare
# Feature, a FeatureCollection, or a raw geometry; closes open rings; dissolves
# multiple features in one file into a single (multi)polygon. jsonlite is used
# (simplifyVector = FALSE) so nested coordinate arrays parse predictably across
# GDAL versions.
.ring_matrix <- function(ring) {
  m <- do.call(rbind, lapply(ring, function(pt) c(as.numeric(pt[[1]]), as.numeric(pt[[2]]))))
  if (!isTRUE(all.equal(m[1, ], m[nrow(m), ]))) m <- rbind(m, m[1, ])  # close ring
  m
}

read_polygon_geojson <- function(path) {
  .require_sf("read_polygon_geojson()")
  j <- jsonlite::fromJSON(path, simplifyVector = FALSE)

  # FOUR SHAPES, because QuPath writes more than one and the arms differ.
  #   {"type":"FeatureCollection","features":[...]}   the standard container
  #   {"type":"Feature","geometry":{...}}             one annotation
  #   [{"type":"Feature",...}, ...]                   a BARE ARRAY of features —
  #                                                   no envelope, so `$type` is NULL
  #   {"type":"Polygon", ...}                         a bare geometry
  #
  # The bare array is what `annotation_all/` holds. Without a branch for it, `$type`
  # is NULL, the geometry-only fallback wraps the whole ARRAY as one geometry, and
  # `g$type` is NULL too — so it stopped with `unsupported geometry type: ` and an
  # EMPTY name. load_annotations()/arm_annotations() catch that and warn-and-skip, so
  # the symptom was every union polygon in arm 1 silently missing: union rows falling
  # back to the export flag, and 10338/15897 losing the polygon their promoted
  # ANNOTATION_1 depends on. A whole scope quietly degraded on a parser branch.
  .is_feature <- function(x) is.list(x) && !is.null(x$geometry)
  feats <- if (identical(j$type, "FeatureCollection")) j$features
           else if (identical(j$type, "Feature"))       list(j)
           else if (is.null(names(j)) && length(j) > 0 &&
                    all(vapply(j, .is_feature, logical(1)))) j
           else                                          list(list(geometry = j))

  geoms <- lapply(feats, function(f) {
    g <- f$geometry
    if (identical(g$type, "Polygon")) {
      sf::st_polygon(lapply(g$coordinates, .ring_matrix))
    } else if (identical(g$type, "MultiPolygon")) {
      sf::st_multipolygon(lapply(g$coordinates,
        function(poly) lapply(poly, .ring_matrix)))
    } else {
      stop(sprintf("unsupported geometry type '%s' in %s",
                   g$type %||% "<none>", path))
    }
  })
  # planar image pixel coordinates -> no CRS; dissolve to one geometry.
  # st_make_valid: QuPath annotation rings are often self-intersecting, and
  # st_within() silently returns NO matches on an invalid polygon (every cell
  # reads as outside). Repairing the geometry is what makes point-in-polygon work.
  sf::st_make_valid(sf::st_union(sf::st_sfc(geoms)))
}

# Load every annotation polygon under `dir`, keyed by patient and ANNOTATION_<k>.
#
# TWO LAYOUTS, ONE READER. The producer changed shape, and both shapes have to load:
#
#   nested (current)  <dir>/<patient>/<patient>_<A|B|C>.geojson
#   flat   (legacy)   <dir>/<patient>_a<k>.geojson   (or a bare <patient>.geojson)
#
# The nested tree is what the cluster share writes, and its region suffix is a LETTER
# whose alphabet position is the annotation index — so `C` is ANNOTATION_3 whether or
# not `B` exists. The flat tree used a digit. Rather than ask every caller which it
# has, this recurses and reads either, because the four pages that consume `annots`
# (clinical_data, clinical_data_mirage, periphery, annotation_membership_qc) do not
# care and should not have to.
#
# THE DIRECTORY NAME WINS over the filename stem when they disagree. A file that was
# renamed by hand would otherwise be keyed to the wrong patient, silently, and the
# only symptom would be an implausible in-annotation count.
#
# Returns NULL rather than erroring when there is nothing to read: a cohort can
# legitimately be all-whole-slide (a patient with no annotation directory is entirely
# inside — see code/arm_cells.R), and that must not stop a report at the loader.
#
# THERE IS NO DEFAULT DIRECTORY ANY MORE. `ANNOTATION_DIR()` used to answer "the"
# annotation tree, which stopped being a question with one answer when the second and
# third phenotyping arms arrived: three trees now hold polygons and they are drawn
# INDEPENDENTLY of each other, so a page that took the default would silently compare
# one arm's cells against another arm's regions. `arm_annotation_dir(arm)` makes the
# arm part of the call, and `arm_annotations(spec, tier)` in arm_cells.R is what the
# arm-aware pages actually use. This generic reader survives for a geojson tree that
# belongs to no arm.
arm_annotation_dir <- function(arm = ARM_MODES, tier = c("region", "union")) {
  spec <- arm_spec(arm)
  tier <- match.arg(tier)
  key  <- if (tier == "region") "region_poly" else "union_poly"
  if (is.null(spec[[key]]))
    stop("arm \"", spec$arm, "\" has no ", tier, " polygon tier")
  spec[[key]]$path
}

# "A" -> ANNOTATION_1, by alphabet position, not by file order.
annotation_from_letter <- function(letter) {
  idx <- match(toupper(letter), LETTERS)
  ifelse(is.na(idx), NA_character_, paste0("ANNOTATION_", idx))
}

# Patient + annotation for one per-region FILE, given the directory it sat in. Used
# for geojsons here and for the cell csvs in annotation_membership_qc(): the two trees
# name their regions identically, and having one parser is what keeps a geojson and its
# csv agreeing on which region they are.
# DELEGATES TO arm_parse_name(); it does not re-implement it. There is exactly one
# place in the project that knows how a filename names its region, and it is
# code/arms.R — a second copy here is how a geojson and the csv it is compared
# against would drift into disagreeing about which region they are.
#
# The one thing this adds on top is the LEGACY DEFAULT: a bare `<patient>.geojson`
# in a tree that belongs to no arm is that patient's sole annotation, so it is
# ANNOTATION_1. An arm's tiers deliberately do NOT apply that rule — a bare geojson
# in massimo1's `_all` tree is a union polygon, not region 1 — which is why
# arm_parse_name() returns NA there and leaves the reading to its caller.
.annotation_key <- function(path, root) {
  nested <- !identical(fs::path_norm(fs::path_dir(path)), fs::path_norm(root))
  for (pat in if (nested) c("nested_letter", "nested_digit")
              else        c("flat_letter", "flat_digit")) {
    k <- arm_parse_name(path, pat)
    if (!is.na(k$annotation)) return(k)
  }
  list(patient = arm_parse_name(path, if (nested) "nested_bare" else "flat_digit")$patient,
       annotation = "ANNOTATION_1")
}

load_annotations <- function(dir, patient_ids = NULL) {
  .require_sf("load_annotations()")
  if (!dir.exists(dir)) {
    message("load_annotations(): no annotation directory at ", dir,
            " — every patient will be treated as whole-slide")
    return(NULL)
  }
  files <- fs::dir_ls(dir, glob = "*.geojson", recurse = TRUE, type = "file")
  if (length(files) == 0) {
    message("load_annotations(): no geojson under ", dir)
    return(NULL)
  }

  keys <- lapply(as.character(files), .annotation_key, root = dir)
  meta <- tibble::tibble(
    path       = as.character(files),
    patient_id = vapply(keys, `[[`, character(1), "patient"),
    annotation = vapply(keys, `[[`, character(1), "annotation"))

  if (!is.null(patient_ids))
    meta <- meta[norm_id(meta$patient_id) %in% norm_id(patient_ids), , drop = FALSE]
  if (nrow(meta) == 0) {
    message("load_annotations(): no geojson matched the requested patients")
    return(NULL)
  }

  polys <- purrr::pmap(meta, function(path, patient_id, annotation) {
    geom <- tryCatch(read_polygon_geojson(path), error = function(e) {
      warning("load_annotations(): unreadable polygon, skipping ", path,
              " — ", conditionMessage(e))
      NULL
    })
    if (is.null(geom)) return(NULL)
    sf::st_sf(patient_id = patient_id, annotation = annotation, geometry = geom)
  })
  polys <- purrr::compact(polys)
  if (!length(polys)) return(NULL)
  do.call(rbind, polys)  # sf provides rbind(); preserves geometry across versions
}

# --- Cell-in-annotation metrics ---------------------------------------------
# Lineages tracked in every region (immune subsets + stroma).
region_lineages <- c("CD8T", "CD4T", "Treg", "NK", "Immune_other", "Stroma")

# Per-region counts + ratios for one set of cells inside a region. Emits the three
# ATTEND denominators as raw counts (n_inside / n_tumor_inside / n_cd45_inside) and
# per-lineage counts (n_<lineage>), so composition can be normalised any of ATTEND's
# three ways downstream (see region_composition()). `frac_<lineage>` is the default
# normalisation (per all cells inside). Returns one row.
region_ratios <- function(cells) {
  n_inside <- nrow(cells)
  is_tumor <- if (n_inside) stringr::str_detect(tidyr::replace_na(cells$phenotype_clean, ""), "Tumor") else logical(0)
  is_cd45  <- if (n_inside) marker_pos(cells, "CD45") else logical(0)
  # Marker-gated (not phenotype-gated) subsets: CD3+CD45+ = marker-defined T cells;
  # GZMB+ NK = cytotoxic/activated NK (lineage NK AND granzyme-B positive).
  is_cd3   <- if (n_inside) marker_pos(cells, "CD3")  else logical(0)
  is_gzmb  <- if (n_inside) marker_pos(cells, "GZMB") else logical(0)
  lin      <- if (n_inside) cell_lineage(cells$phenotype_clean) else character(0)
  is_nk    <- if (n_inside) lin %in% "NK" else logical(0)
  # "Clean" cell set = cells inside the region the exporter did NOT flag as an
  # Outlier and that carry a resolved phenotype. is_unresolved_phenotype() covers
  # FlowPath's "Unknown" AND mirage's four reserved outcomes (Unclassified /
  # Ambiguous / Conflict / Artefact); matching only "Unknown" would keep mirage's
  # unresolved cells in the clean denominator and drop FlowPath's, which is exactly
  # the asymmetry that makes the two tools look different when they are not.
  # An export with no `Outlier` column has no outliers (same convention as
  # marker_pos(): absent means "not measured", which cannot exclude a cell).
  is_out     <- if (n_inside && "Outlier" %in% names(cells)) is_pos(cells$Outlier)
                else rep(FALSE, n_inside)
  is_unknown <- if (n_inside) is_unresolved_phenotype(cells$phenotype_clean)
                else logical(0)
  keep_clean <- if (n_inside) !is_out & !is_unknown else logical(0)
  n_tumor  <- sum(is_tumor, na.rm = TRUE)
  n_cd45   <- sum(is_cd45,  na.rm = TRUE)
  n_nk     <- sum(is_nk,    na.rm = TRUE)
  n_cd3cd45 <- sum(is_cd45 & is_cd3, na.rm = TRUE)   # CD45+ AND CD3+ (marker T cells)
  n_gzmb_nk <- sum(is_nk  & is_gzmb, na.rm = TRUE)   # NK lineage AND GZMB+
  n_inside_clean <- sum(keep_clean, na.rm = TRUE)              # cells minus outliers/Unknown
  n_tumor_clean  <- sum(is_tumor & keep_clean, na.rm = TRUE)   # tumour cells in that clean set
  ncount   <- function(l) sum(lin == l, na.rm = TRUE)
  safe     <- function(num, den) if (den > 0) num / den else NA_real_

  out <- tibble::tibble(
    n_inside          = n_inside,
    n_inside_clean    = n_inside_clean,
    n_tumor_inside    = n_tumor,
    n_tumor_clean     = n_tumor_clean,
    n_cd45_inside     = n_cd45,
    n_cd3cd45_inside  = n_cd3cd45,
    n_gzmb_nk_inside  = n_gzmb_nk,
    tumor_over_inside = safe(n_tumor, n_inside),
    # Same tumour fraction, but over the "clean" denominator (Outlier-flagged and
    # Unknown-phenotype cells removed from both numerator and denominator), so it
    # stays a proper 0..1 fraction comparable to the pathologist score.
    tumor_over_inside_clean = safe(n_tumor_clean, n_inside_clean),
    cd45_over_inside  = safe(n_cd45,  n_inside),
    tumor_over_cd45   = safe(n_tumor, n_cd45),  # tumour cells per CD45+ cell inside
    # CD3+CD45+ T cells: as a share of all cells, and of the CD45+ compartment.
    cd3cd45_over_inside = safe(n_cd3cd45, n_inside),
    cd3cd45_over_cd45   = safe(n_cd3cd45, n_cd45),
    # GZMB+ NK: as a share of all cells, and the GZMB+ fraction WITHIN NK cells
    # (an NK-activation readout).
    gzmb_nk_over_inside = safe(n_gzmb_nk, n_inside),
    gzmb_nk_over_nk     = safe(n_gzmb_nk, n_nk)
  )
  for (l in region_lineages) out[[paste0("n_", l)]]    <- ncount(l)
  for (l in region_lineages) out[[paste0("frac_", l)]] <- safe(ncount(l), n_inside)
  out
}

# ATTEND-style multi-normalisation composition (mirrors code/attend_ihc.R
# ihc_celltype_metrics). From a region-metrics table (rows from region_ratios,
# carrying patient_id/n_* counts), returns LONG per lineage with the same three
# denominators ATTEND uses:
#   frac_inside = lineage / all cells inside   (overall composition)
#   frac_tumor  = lineage / tumour cells inside (immune-to-tumour density)
#   frac_cd45   = lineage / CD45+ cells inside  (composition of the immune compartment)
# `markers` optionally adds marker-gated populations (label -> count column from
# region_ratios, e.g. c("CD3+ CD45+" = "n_cd3cd45_inside")) alongside the phenotype
# lineages; default NULL keeps the lineage-only behaviour every existing caller
# relies on. Each population is materialised under a collision-free "pop::" name so
# a count that is also a denominator (n_cd45_inside) can appear as a population too.
region_composition <- function(region_metrics,
                               lineages = c("CD8T", "CD4T", "Treg", "NK", "Immune_other"),
                               markers  = NULL) {
  id_cols <- intersect(c("patient_id", "annotation", "source"), names(region_metrics))
  pop_map <- c(stats::setNames(paste0("n_", lineages), lineages),
               markers[markers %in% names(region_metrics)])

  wide <- region_metrics |>
    dplyr::select(dplyr::all_of(c(id_cols, "n_inside", "n_tumor_inside", "n_cd45_inside")))
  for (lab in names(pop_map)) wide[[paste0("pop::", lab)]] <- region_metrics[[pop_map[[lab]]]]

  wide |>
    tidyr::pivot_longer(dplyr::starts_with("pop::"),
                        names_to = "lineage", names_prefix = "pop::", values_to = "n_cell") |>
    dplyr::mutate(
      lineage     = factor(lineage, levels = names(pop_map)),  # keep call-site order
      frac_inside = dplyr::if_else(n_inside       > 0, n_cell / n_inside,       NA_real_),
      frac_tumor  = dplyr::if_else(n_tumor_inside > 0, n_cell / n_tumor_inside, NA_real_),
      frac_cd45   = dplyr::if_else(n_cd45_inside  > 0, n_cell / n_cd45_inside,  NA_real_)
    )
}

# sf point-in-polygon membership for ONE patient's cells against its polygons.
# The geojson is in PIXELS; cell_centroids_px() puts the cells in that frame by the
# FIXED conversion centroid / um_per_px (0.325) for a micron export, or passes a
# mirage x_px/y_px export straight through. (No per-patient scale search — that
# heuristic is what produced the wrong pixel/micron scale before.) Returns one
# metrics row per polygon (per_annotation) or one over the dissolved union.
.annotation_metrics_sf <- function(cells_p, polys_p, scope, um_per_px) {
  .require_sf("geojson membership")
  pts <- sf::st_as_sf(cell_centroids_px(cells_p, um_per_px),
                      coords = c("x", "y"), crs = sf::st_crs(polys_p))
  # The geojson polygon gives the region AREA, so region_ratios_area() adds area_mm2
  # and cell DENSITIES (cells / mm^2, incl. tumour density) inside the annotation.
  if (scope == "union") {
    poly_u  <- sf::st_union(sf::st_geometry(polys_p))
    inside  <- lengths(sf::st_within(pts, poly_u)) > 0
    area_px <- as.numeric(sum(sf::st_area(poly_u)))
    region_ratios_area(cells_p[inside, , drop = FALSE], area_px, um_per_px) |>
      dplyr::mutate(annotation = "union", source = "sf", .before = 1)
  } else {
    within <- sf::st_within(pts, polys_p)          # per-cell list of polygon indices
    purrr::map_dfr(seq_len(nrow(polys_p)), function(i) {
      inside_i <- vapply(within, function(idx) i %in% idx, logical(1))
      area_px  <- as.numeric(sf::st_area(sf::st_geometry(polys_p)[i]))
      region_ratios_area(cells_p[inside_i, , drop = FALSE], area_px, um_per_px) |>
        dplyr::mutate(annotation = polys_p$annotation[i], source = "sf", .before = 1)
    })
  }
}

# Region metrics per patient. PRIMARY membership = sf point-in-polygon on the raw
# pathologist geojson (`annots`, fixed um_per_px scale); FALLBACK = the FlowPath
# `Out_of_annotation` flag (a cell is inside when the flag is FALSE) for patients
# with NO geojson. `ihc_data` supplies the cells (centroids + optional flag).
# `scope` = "per_annotation" (one row per polygon) or "union" (dissolved). The
# `source` column records sf vs csv provenance per patient; the CSV fallback has no
# per-polygon breakdown, so it always emits a single row labelled annotation "csv".
ihc_annotation_metrics <- function(ihc_data, annots,
                                   scope = c("per_annotation", "union"),
                                   use_csv_fallback = TRUE, um_per_px = 0.325) {
  scope    <- match.arg(scope)
  ihc_data <- dplyr::mutate(ihc_data, .pid = slide_key(patient_id))
  annots   <- dplyr::mutate(annots,   .pid = slide_key(patient_id))
  ann_ids  <- unique(annots$.pid)

  purrr::map_dfr(unique(ihc_data$.pid), function(pid) {
    cells_p <- dplyr::filter(ihc_data, .pid == pid)
    xy      <- cell_centroids_px(cells_p, um_per_px)
    cells_p <- cells_p[is.finite(xy$x) & is.finite(xy$y), , drop = FALSE]
    if (nrow(cells_p) == 0) return(tibble::tibble())

    if (pid %in% ann_ids) {
      polys_p <- dplyr::filter(annots, .pid == pid)
      .annotation_metrics_sf(cells_p, polys_p, scope, um_per_px) |>
        dplyr::mutate(patient_id = pid, .before = 1)
    } else if (use_csv_fallback && has_outside_flag(cells_p)) {
      # no polygon -> no area, so area_mm2 / densities come back NA (schema matches sf).
      inside <- !cell_outside(cells_p)
      region_ratios_area(cells_p[inside, , drop = FALSE], 0, um_per_px) |>
        dplyr::mutate(patient_id = pid, annotation = "csv", source = "csv", .before = 1)
    } else {
      tibble::tibble()  # no geojson and no CSV flag -> nothing to report
    }
  })
}

# QC: does sf point-in-polygon agree with the precomputed `Out_of_annotation` flag?
# For every per-annotation FlowPath csv `<patient>_a<k>.csv` in `dir` (each carries a
# flag computed for ANNOTATION_k), compute BOTH memberships on the SAME cells —
#   flag inside = !Out_of_annotation
#   sf inside   = st_within(centroid / um_per_px, the ANNOTATION_k geojson polygon)
# — so agreement is a direct per-cell comparison (no cross-file cell matching). One
# row per (patient, annotation): set sizes, overlap (Jaccard, % of cells agreeing),
# and the tumour fraction inside under each method. A large flag-vs-sf gap on a slide
# means the fixed um_per_px is wrong FOR THAT slide (sf catching the wrong cells).
# TAKES AN ARM, NOT A DIRECTORY. The check is only meaningful when the csvs and the
# polygons come from the SAME arm — an arm's regions are drawn independently of every
# other arm's, so comparing one arm's flag against another's polygon would report a
# disagreement that is real geometry rather than a coordinate-scale problem. Passing
# `spec` makes that pairing structural: the csvs come from spec$region_csv and are
# keyed with THAT tier's own naming pattern, so a flat `<pid>_a<k>.csv` tree and a
# nested `<pid>/<pid>_<L>.csv` one both read without the caller saying which it has.
annotation_membership_qc <- function(spec, annots, um_per_px = 0.325) {
  .require_sf("annotation_membership_qc()")
  if (is.character(spec)) spec <- arm_spec(spec)
  tier <- spec$region_csv
  dir  <- tier$path
  if (!dir.exists(dir)) {
    warning("annotation_membership_qc(): no directory at ", dir)
    return(tibble::tibble())
  }
  files  <- fs::dir_ls(dir, glob = "*.csv", recurse = TRUE, type = "file")
  if (length(files) == 0) {
    warning("annotation_membership_qc(): no csv under ", dir)
    return(tibble::tibble())
  }
  annots <- dplyr::mutate(annots, .pid = slide_key(patient_id))

  purrr::map_dfr(as.character(files), function(path) {
    key  <- arm_parse_name(path, tier$pattern)
    pid  <- slide_key(key$patient)
    ann  <- if (is.na(key$annotation)) spec$bare_region_is else key$annotation
    if (is.na(ann))
      return(tibble::tibble(patient_id = pid, annotation = NA_character_,
                            note = "csv carries no region suffix and this arm has no rule for one"))

    cells <- tibble::as_tibble(data.table::fread(path))
    if (!has_centroids(cells) || !has_outside_flag(cells))
      return(tibble::tibble(patient_id = pid, annotation = ann,
                            note = "csv missing centroid / out-of-annotation flag"))
    poly <- dplyr::filter(annots, .pid == pid, annotation == ann)
    if (nrow(poly) == 0)
      return(tibble::tibble(patient_id = pid, annotation = ann,
                            n_cells = nrow(cells), note = "no matching geojson"))

    pts     <- sf::st_as_sf(cell_centroids_px(cells, um_per_px),
                            coords = c("x", "y"), crs = sf::st_crs(poly))
    sf_in   <- lengths(sf::st_within(pts, sf::st_geometry(poly))) > 0
    flag_in <- !cell_outside(cells)
    is_tum  <- stringr::str_detect(tidyr::replace_na(cell_phenotype(cells), ""), "Tumor")

    n_union <- sum(sf_in | flag_in)
    frac    <- function(sel) if (sum(sel) > 0) sum(is_tum & sel) / sum(sel) else NA_real_
    tibble::tibble(
      patient_id      = pid, annotation = ann,
      n_cells         = nrow(cells),
      n_flag_inside   = sum(flag_in),
      n_sf_inside     = sum(sf_in),
      n_both          = sum(sf_in & flag_in),
      jaccard         = if (n_union > 0) sum(sf_in & flag_in) / n_union else NA_real_,
      pct_cells_agree = mean(sf_in == flag_in),
      tumor_frac_flag = frac(flag_in),
      tumor_frac_sf   = frac(sf_in),
      note            = NA_character_
    )
  })
}

# --- Whole-slide IHC quantities (for the RNA comparisons) -------------------
# Per-patient fraction of cells positive for each marker, plus the mean z-score as a
# sensitivity readout. Whole slide (bulk RNA has no annotation boundary). Returns
# wide: patient_id + <marker>_posfrac + <marker>_z. A marker the export never gated
# comes back 0 (no positive cells) with an NA z-score, so a partial panel narrows the
# comparison instead of aborting it.
ihc_marker_fraction <- function(ihc_data, markers = ihc_markers) {
  wide <- tibble::tibble(patient_id = slide_key(ihc_data$patient_id))
  for (m in markers) {
    wide[[paste0(m, "_posfrac")]] <- marker_pos(ihc_data, m)
    wide[[paste0(m, "_z")]]       <- marker_zscore(ihc_data, m)
  }
  # An all-NA column (a marker with no z-score in this export) must come back NA,
  # not the NaN that mean(na.rm = TRUE) returns on an empty vector.
  wide |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(dplyr::across(dplyr::everything(),
                                   ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)),
                     .groups = "drop")
}

# Long form of the per-cell marker table: one row per (cell, marker) carrying the
# value / z-score / positivity triple. Used by the internal-QC reports (phenotype-
# marker concordance, per-channel usability). wrongL1CAM/DAPI are excluded via the
# default marker list. Built marker by marker through the cell_tables.R accessors,
# so the triple is found under any of the three export spellings.
ihc_marker_long <- function(ihc_data, markers = ihc_markers) {
  keep <- intersect(c("cell_id", "label", "patient_id", "phenotype_clean"),
                    names(ihc_data))
  base <- dplyr::mutate(ihc_data[keep], patient_id = slide_key(patient_id))
  purrr::map_dfr(markers, function(m) dplyr::mutate(
    base,
    marker = m,
    raw    = marker_value(ihc_data, m),
    zscore = marker_zscore(ihc_data, m),
    sign   = marker_sign(ihc_data, m),    # verbatim, for display
    gated  = marker_gated(ihc_data, m),   # FALSE = never measured, NOT negative
    pos    = marker_pos(ihc_data, m)
  ))
}

# Per-patient lineage composition over ALL cells (whole slide), for deconvolution
# comparison. Long: patient_id, lineage, n, frac.
ihc_lineage_fraction <- function(ihc_data) {
  ihc_data |>
    dplyr::mutate(patient_id = slide_key(patient_id)) |>
    dplyr::mutate(lineage = cell_lineage(phenotype_clean)) |>
    dplyr::count(patient_id, lineage, name = "n") |>
    dplyr::group_by(patient_id) |>
    dplyr::mutate(frac = n / sum(n)) |>
    dplyr::ungroup()
}

# Map an immunedeconv cell_type string to one of `comparable_lineages` (or NA).
# Pattern-based because the exact strings differ across methods.
deconv_to_lineage <- function(cell_type) {
  ct <- tolower(cell_type)
  dplyr::case_when(
    stringr::str_detect(ct, "regulatory|treg")        ~ "Treg",
    stringr::str_detect(ct, "cd8")                    ~ "CD8T",
    stringr::str_detect(ct, "cd4")                    ~ "CD4T",
    stringr::str_detect(ct, "^nk|nk cell|natural killer") ~ "NK",
    TRUE                                              ~ NA_character_
  )
}

# =============================================================================
# Tumour border / periphery (invasive-margin) metrics
# -----------------------------------------------------------------------------
# The tumour annotation gives an "inside vs outside" partition; the biology at
# the tumour-host interface needs a third region — the invasive margin (IM), a
# band centred on the annotation boundary. Literature grounding for the band
# half-width `d` (microns each side of the border):
#   * Consensus Immunoscore (Galon 2014 J Pathol; Pages 2018 Lancet) scores CD3/
#     CD8 densities in the tumour core (CT) and a 500 um invasive margin (IM) —
#     the de-facto standard, and the default in QuPath IM workflows.
#   * Reviews put the IM width in a 200-500 um range, with some invasive-margin
#     detection algorithms using up to 1000 um.
#   * Immune spatial phenotypes (Chen & Mellman; Hegde 2016) are defined by WHERE
#     CD8 sits relative to this border: inflamed/hot = throughout the core,
#     excluded/cold = trapped at the margin, desert/cold = sparse everywhere.
# So we report several thresholds (default 100 / 250 / 500 um) and, per threshold,
# both the margin band and the eroded "deep core" so a margin-vs-core contrast can
# recover the inflamed/excluded/desert axis.
#
# All geometry is done in the polygon's PIXEL coordinate frame. Cell centroids are
# in microns and the QuPath geojson is in pixels, so centroids are mapped in by
# dividing by `um_per_px` (= 0.325): x_px = centroid_x / 0.325. One pixel is then
# `um_per_px` microns, so a micron threshold `d` becomes `d / um_per_px` pixel
# units, and a pixel area is scaled by `um_per_px^2` to reach microns^2.

# From a dissolved core polygon (sfc) and a signed buffer distance `d` (polygon
# units), return the invasive-margin band (+/- d around the boundary) and the
# eroded interior "deep core" (further than d inside the boundary). If d exceeds
# the core's inradius the eroded core is empty (sfc of length 0) — the caller
# treats that as "no deep core at this threshold".
margin_regions <- function(core, d) {
  outer <- suppressWarnings(sf::st_buffer(core, d))
  inner <- suppressWarnings(sf::st_buffer(core, -d))
  band  <- if (length(inner) == 0 || all(sf::st_is_empty(inner))) outer
           else suppressWarnings(sf::st_difference(outer, inner))
  list(margin = band, core = inner)
}

# region_ratios() augmented with region AREA and area-normalised DENSITIES
# (cells / mm^2), the native Immunoscore unit. `area_units2` is the region area in
# squared polygon units; `um_per_unit` converts it to mm^2. Densities are NA when
# the region is empty/degenerate so they never masquerade as a real zero.
region_ratios_area <- function(cells, area_units2, um_per_unit) {
  rr       <- region_ratios(cells)
  area_mm2 <- area_units2 * (um_per_unit^2) / 1e6         # units^2 -> um^2 -> mm^2
  dens     <- function(n) if (is.finite(area_mm2) && area_mm2 > 0) n / area_mm2 else NA_real_
  rr$area_mm2           <- if (is.finite(area_mm2) && area_mm2 > 0) area_mm2 else NA_real_
  rr$dens_all_per_mm2   <- dens(rr$n_inside)
  rr$dens_tumor_per_mm2 <- dens(rr$n_tumor_inside)   # tumour cells / mm^2 (neoplastic density)
  rr$dens_cd45_per_mm2  <- dens(rr$n_cd45_inside)
  for (l in c("CD8T", "CD4T", "Treg", "NK"))
    rr[[paste0("dens_", l, "_per_mm2")]] <- dens(rr[[paste0("n_", l)]])
  rr
}

# Invasive-margin metrics per patient. For every requested `thresholds_um` and,
# per `scope`, either the dissolved union (one core) or each single annotation
# (one core each), returns TWO rows — region = "margin" (the +/- d band) and
# region = "core" (interior beyond d) — carrying every region_ratios_area column.
# Only patients WITH a geojson polygon are covered: a margin band cannot be built
# from the binary CSV Out_of_annotation flag, so there is no CSV fallback here.
# `um_per_unit` and `source` are echoed so the physical band width is auditable.
ihc_periphery_metrics <- function(ihc_data, annots,
                                  scope = c("union", "per_annotation"),
                                  thresholds_um = c(100, 250, 500),
                                  um_per_px = 0.325) {
  .require_sf("ihc_periphery_metrics()")
  scope    <- match.arg(scope)
  ihc_data <- dplyr::mutate(ihc_data, .pid = slide_key(patient_id))
  annots   <- dplyr::mutate(annots,   .pid = slide_key(patient_id))
  ann_ids  <- unique(annots$.pid)

  purrr::map_dfr(intersect(unique(ihc_data$.pid), ann_ids), function(pid) {
    cells_p <- dplyr::filter(ihc_data, .pid == pid)
    xy      <- cell_centroids_px(cells_p, um_per_px)
    cells_p <- cells_p[is.finite(xy$x) & is.finite(xy$y), , drop = FALSE]
    if (nrow(cells_p) == 0) return(tibble::tibble())
    polys_p <- dplyr::filter(annots, .pid == pid)
    upu     <- um_per_px    # microns per polygon (pixel) unit
    pts     <- sf::st_as_sf(cell_centroids_px(cells_p, um_per_px),
                            coords = c("x", "y"), crs = sf::st_crs(polys_p))

    cores <- if (scope == "union") {
      list(union = sf::st_union(sf::st_geometry(polys_p)))
    } else {
      stats::setNames(lapply(seq_len(nrow(polys_p)),
                             function(i) sf::st_geometry(polys_p)[i]),
                      polys_p$annotation)
    }

    purrr::imap_dfr(cores, function(core, ann_label) {
      purrr::map_dfr(thresholds_um, function(d_um) {
        rg <- margin_regions(core, d_um / upu)
        row_for <- function(region_geom, region_label) {
          empty  <- length(region_geom) == 0 || all(sf::st_is_empty(region_geom))
          inside <- if (empty) rep(FALSE, nrow(cells_p))
                    else lengths(sf::st_within(pts, region_geom)) > 0
          area_u2 <- if (empty) 0 else as.numeric(sum(sf::st_area(region_geom)))
          region_ratios_area(cells_p[inside, , drop = FALSE], area_u2, upu) |>
            dplyr::mutate(region = region_label, .before = 1)
        }
        dplyr::bind_rows(row_for(rg$margin, "margin"),
                         row_for(rg$core,   "core")) |>
          dplyr::mutate(patient_id = pid, annotation = ann_label,
                        threshold_um = d_um, um_per_unit = upu,
                        source = "sf:um2px(/0.325)", .before = 1)
      })
    })
  })
}

# Spearman rho + p + n for two paired numeric vectors, as a one-row tibble.
# Guards the small-n / zero-variance cases that abort cor.test().
paired_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3 || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(tibble::tibble(n = length(x), rho = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman"))
  tibble::tibble(n = length(x), rho = unname(ct$estimate), p = ct$p.value)
}

# Pearson r, Spearman rho and Kendall tau for two paired numeric vectors, as a
# one-row tibble with each estimate + its p-value + the shared n. Same small-n /
# zero-variance guard as paired_spearman(). Use where all three are wanted side by
# side (Pearson = linear, Spearman/Kendall = rank/monotonic and robust to outliers).
paired_cor3 <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]; n <- length(x)
  na_row <- tibble::tibble(n = n,
                           pearson = NA_real_,  p_pearson  = NA_real_,
                           spearman = NA_real_, p_spearman = NA_real_,
                           kendall = NA_real_,  p_kendall  = NA_real_)
  if (n < 3 || stats::sd(x) == 0 || stats::sd(y) == 0) return(na_row)
  ct <- function(m) suppressWarnings(stats::cor.test(x, y, method = m))
  pe <- ct("pearson"); sp <- ct("spearman"); ke <- ct("kendall")
  tibble::tibble(n = n,
                 pearson  = unname(pe$estimate), p_pearson  = pe$p.value,
                 spearman = unname(sp$estimate), p_spearman = sp$p.value,
                 kendall  = unname(ke$estimate), p_kendall  = ke$p.value)
}

# Pairwise cross-source agreement, for comparing several deconvolution methods (and
# IHC) against each other rather than only against one reference. `df` is long with
# columns source / patient_id / lineage / value (one value per source x lineage x
# patient). For every ordered source pair it correlates each lineage across the
# shared patients (Spearman, rank-based so differing score scales are fine) and
# averages the per-lineage rho, returning a symmetric long table (a == b -> rho 1).
# Different methods put scores on different scales; only rank agreement is meaningful.
pairwise_agreement <- function(df, sources = sort(unique(df$source))) {
  grid <- expand.grid(a = sources, b = sources, stringsAsFactors = FALSE)
  purrr::pmap_dfr(grid, function(a, b) {
    if (a == b)
      return(tibble::tibble(a = a, b = b, mean_rho = 1, n_lineages = NA_integer_))
    j <- dplyr::inner_join(
      dplyr::filter(df, source == a),
      dplyr::filter(df, source == b),
      by = c("patient_id", "lineage"), suffix = c("_a", "_b"))
    if (nrow(j) == 0)
      return(tibble::tibble(a = a, b = b, mean_rho = NA_real_, n_lineages = 0L))
    per <- j |>
      dplyr::group_by(lineage) |>
      dplyr::group_modify(~ paired_spearman(.x$value_a, .x$value_b)) |>
      dplyr::ungroup()
    tibble::tibble(a = a, b = b,
                   mean_rho   = mean(per$rho, na.rm = TRUE),
                   n_lineages = sum(is.finite(per$rho)))
  })
}
