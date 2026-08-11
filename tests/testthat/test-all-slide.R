# The "all-slide" export layout: one CSV *and* one geojson per annotation region,
# nested per patient. These tests build a synthetic tree rather than needing the
# cluster share, so they run anywhere — the geometry tests skip without sf, but the
# layout, letter-mapping and no-annotation rules are checked unconditionally,
# because those are the parts that decide whether a patient is counted at all.
source(here::here("code", "cell_tables.R"))
source(here::here("code", "validation_helpers.R"))
source(here::here("code", "all_slide.R"))

# A cell table spread over a known box, so a polygon covering half of it is
# expected to keep roughly half the cells.
.cells_csv <- function(path, n, seed = 1) {
  set.seed(seed)
  readr::write_csv(tibble::tibble(
    cell_id           = seq_len(n),
    x_px              = seq(0, 3999, length.out = n),
    y_px              = seq(0, 2999, length.out = n),
    phenotype         = rep_len(c("PANCK+Tumor", "T cytotoxic", "Immune", "Unknown"), n),
    Out_of_annotation = rep_len(c("False", "False", "False", "True"), n)), path)
}

.poly_geojson <- function(path, x0, y0, x1, y1) {
  jsonlite::write_json(list(
    type = "Feature",
    geometry = list(type = "Polygon", coordinates = list(list(
      c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0)))),
    properties = structure(list(), names = character(0))), path, auto_unbox = TRUE)
}

# regions = named list of patient -> region letters; NULL letters means one bare
# whole-slide csv. `annotate` names the patients that get an annotation directory.
all_slide_tree <- function(regions = list(`046` = c("A", "B"), `24086` = NULL),
                           annotate = "046", n = 200) {
  root <- file.path(tempdir(), paste0("allslide-", as.integer(Sys.time()), "-", sample(1e6, 1)))
  for (pid in names(regions)) {
    ls <- regions[[pid]]
    dir.create(file.path(root, "csv", pid), recursive = TRUE, showWarnings = FALSE)
    if (is.null(ls)) {
      .cells_csv(file.path(root, "csv", pid, paste0(pid, ".csv")), n)
    } else {
      for (L in ls) .cells_csv(file.path(root, "csv", pid, sprintf("%s_%s.csv", pid, L)), n)
    }
    if (pid %in% annotate && !is.null(ls)) {
      dir.create(file.path(root, "annotation", pid), recursive = TRUE, showWarnings = FALSE)
      # Region A covers the lower-left half, B the upper-right — disjoint, so the
      # union is a strictly larger area than either.
      boxes <- list(A = c(0, 0, 2000, 1500), B = c(2000, 1500, 4000, 3000),
                    C = c(0, 1500, 2000, 3000))
      for (L in ls) do.call(.poly_geojson,
        c(list(file.path(root, "annotation", pid, sprintf("%s_%s.geojson", pid, L))),
          as.list(boxes[[L]])))
    }
  }
  root
}

test_that("region letters map to ANNOTATION_<k> by alphabet position", {
  expect_equal(all_slide_region("A"), "ANNOTATION_1")
  expect_equal(all_slide_region("B"), "ANNOTATION_2")
  # C stays region 3 whether or not B was exported — position, not file order.
  expect_equal(all_slide_region("c"), "ANNOTATION_3")
  expect_true(is.na(all_slide_region("?")))
})

test_that("a patient with no annotation directory is entirely inside", {
  root  <- all_slide_tree()
  cells <- all_slide_cells(root)
  ws    <- dplyr::filter(cells, patient_id == "24086")

  expect_equal(nrow(ws), 200)
  expect_equal(unique(ws$annotation), "whole_slide")
  expect_false(any(ws$has_annotation))
  expect_equal(unique(ws$.membership_source), "whole_slide")

  # The rule that matters: every cell counts, INCLUDING the ones the export's own
  # flag calls out-of-annotation. The layout says there is no annotation to be out
  # of, so the flag must not be allowed to quietly drop a quarter of the slide.
  m <- all_slide_metrics(cells, NULL, "per_annotation")
  expect_equal(m$n_inside[m$patient_id == "24086"], 200)
  expect_equal(m$source[m$patient_id == "24086"], "whole_slide")
})

test_that("region files are keyed by letter and never pooled across regions", {
  cells <- all_slide_cells(all_slide_tree())
  a <- dplyr::filter(cells, patient_id == "046")
  expect_setequal(unique(a$annotation), c("ANNOTATION_1", "ANNOTATION_2"))
  expect_equal(nrow(a), 400)          # two files x 200, kept apart
  expect_true(all(a$has_annotation))
})

test_that("the union cell set de-duplicates cells repeated across region files", {
  cells <- all_slide_cells(all_slide_tree())
  # Both of 046's region files export the SAME slide, so pooling them naively
  # double-counts every cell. The union must be one slide's worth.
  u <- all_slide_union_cells(cells)
  expect_equal(sum(u$patient_id == "046"), 200)
  expect_equal(sum(u$patient_id == "24086"), 200)
})

test_that("metrics fall back to the export flag when a region has no polygon", {
  # Region files present, annotation directory absent for that patient => the
  # patient is whole-slide; but a region file whose polygon is unreadable while
  # OTHER regions have one must fall back to the flag, not be counted whole.
  cells <- all_slide_cells(all_slide_tree())
  m <- all_slide_metrics(cells, NULL, "per_annotation")
  a <- dplyr::filter(m, patient_id == "046")
  expect_equal(unique(a$source), "flag")
  expect_equal(sum(a$n_inside), 300)     # 150 flagged-inside per region file
  expect_true(all(is.na(a$area_mm2)))    # no polygon -> no area -> no densities
})

test_that("metrics keep the ihc_annotation_metrics column schema", {
  cells <- all_slide_cells(all_slide_tree())
  m     <- all_slide_metrics(cells, NULL, "per_annotation")
  ref   <- names(region_ratios_area(cells[seq_len(10), , drop = FALSE], 0, 0.325))
  expect_true(all(ref %in% names(m)))
  expect_true(all(c("patient_id", "annotation", "source") %in% names(m)))
})

test_that("a missing csv/ directory warns and yields no cells rather than erroring", {
  expect_warning(out <- all_slide_cells(file.path(tempdir(), "definitely-not-there")),
                 "no csv/")
  expect_equal(nrow(out), 0)
})

test_that("the inventory reports what was read, per patient and region", {
  inv <- all_slide_inventory(all_slide_cells(all_slide_tree()))
  expect_equal(nrow(inv), 3)                       # 046 x2 regions + 24086 whole slide
  expect_true(all(c("n_cells", "has_annotation", "membership") %in% names(inv)))
  expect_equal(inv$n_cells, c(200L, 200L, 200L))
})

# ---- geometry: needs sf ------------------------------------------------------

test_that("polygons load with patient and region from the tree, not the filename", {
  skip_if_not_installed("sf")
  root  <- all_slide_tree()
  polys <- all_slide_annotations(root)
  expect_equal(nrow(polys), 2)
  expect_setequal(polys$annotation, c("ANNOTATION_1", "ANNOTATION_2"))
  expect_equal(unique(polys$patient_id), "046")
})

test_that("each region is scored against its OWN polygon", {
  skip_if_not_installed("sf")
  root  <- all_slide_tree()
  cells <- all_slide_cells(root)
  polys <- all_slide_annotations(root)
  m     <- all_slide_metrics(cells, polys, "per_annotation")
  a     <- dplyr::filter(m, patient_id == "046")

  expect_equal(unique(a$source), "sf")
  # The polygon gives the area, so densities exist here where the flag fallback
  # returned NA — that difference is the whole reason to prefer the geometry.
  expect_true(all(is.finite(a$area_mm2)))
  expect_true(all(a$n_inside > 0))
  # The two boxes are disjoint and the cells run along the diagonal, so region 1
  # (lower-left) must hold cells that region 2 does not.
  expect_true(all(a$n_inside < 200))
})

test_that("the union dissolves the polygons and counts each cell once", {
  skip_if_not_installed("sf")
  root  <- all_slide_tree()
  cells <- all_slide_cells(root)
  polys <- all_slide_annotations(root)
  per   <- all_slide_metrics(cells, polys, "per_annotation")
  uni   <- all_slide_metrics(cells, polys, "union")

  u <- dplyr::filter(uni, patient_id == "046")
  expect_equal(nrow(u), 1)
  expect_equal(u$annotation, "union")
  # Area is additive over disjoint regions; cell counts are NOT, because the union
  # is computed on the de-duplicated slide rather than by adding the region rows.
  expect_equal(u$area_mm2, sum(per$area_mm2[per$patient_id == "046"]), tolerance = 1e-6)
  expect_lte(u$n_inside, 200)
})

test_that("a whole-slide patient still appears in the union frame", {
  skip_if_not_installed("sf")
  root  <- all_slide_tree()
  cells <- all_slide_cells(root)
  polys <- all_slide_annotations(root)
  uni   <- all_slide_metrics(cells, polys, "union")
  # 24086 has no polygon to dissolve; dropping it here would silently shrink the
  # cohort from six patients to five in every union-scoped panel.
  expect_true("24086" %in% uni$patient_id)
  expect_equal(uni$source[uni$patient_id == "24086"], "whole_slide")
})

# ---- what the region files actually contain ----------------------------------
# Whether a patient's region files repeat the same cells or partition them is a
# property of the PRODUCER, not of the layout, and it sets every cohort-level
# denominator. The code is correct either way; these tests pin that the report says
# which one the data on disk is, so nobody has to assume.

.region_csv <- function(path, rows) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(rows, path)
}

test_that("repeated cells across region files are reported as whole-slide exports", {
  root <- file.path(tempdir(), paste0("ws-", sample(1e6, 1)))
  rows <- tibble::tibble(cell_id = 1:50, x_px = seq(0, 490, 10), y_px = 5,
                         phenotype = "Immune")
  for (L in c("A", "B", "C"))
    .region_csv(file.path(root, "csv", "046", sprintf("046_%s.csv", L)), rows)
  dir.create(file.path(root, "annotation", "046"), recursive = TRUE)

  cells <- all_slide_cells(root)
  rep   <- all_slide_overlap_report(cells)
  expect_equal(rep$n_rows, 150L)
  expect_equal(rep$n_unique_cells, 50L)
  expect_equal(rep$rows_per_cell, 3)
  expect_match(attr(rep, "verdict"), "whole-slide exports")
  # And the cohort cell set must be one slide's worth, not three.
  expect_equal(nrow(all_slide_union_cells(cells)), 50)
})

test_that("disjoint region files are reported as region-restricted exports", {
  root <- file.path(tempdir(), paste0("rr-", sample(1e6, 1)))
  for (k in seq_along(c("A", "B"))) {
    L <- c("A", "B")[k]
    .region_csv(file.path(root, "csv", "046", sprintf("046_%s.csv", L)),
                tibble::tibble(cell_id = seq_len(50) + (k - 1) * 50,
                               x_px = seq(0, 490, 10) + (k - 1) * 1000,
                               y_px = 5, phenotype = "Immune"))
  }
  dir.create(file.path(root, "annotation", "046"), recursive = TRUE)

  rep <- all_slide_overlap_report(all_slide_cells(root))
  expect_equal(rep$rows_per_cell, 1)
  expect_match(attr(rep, "verdict"), "region-restricted")
})

test_that("a single-region patient cannot answer the question, and says so", {
  root <- file.path(tempdir(), paste0("one-", sample(1e6, 1)))
  .region_csv(file.path(root, "csv", "10338", "10338_A.csv"),
              tibble::tibble(cell_id = 1:10, x_px = 1:10, y_px = 1, phenotype = "Immune"))
  dir.create(file.path(root, "annotation", "10338"), recursive = TRUE)
  rep <- all_slide_overlap_report(all_slide_cells(root))
  expect_match(attr(rep, "verdict"), "cannot tell")
})

test_that("a mixed cohort is flagged rather than averaged into one verdict", {
  root <- file.path(tempdir(), paste0("mix-", sample(1e6, 1)))
  rows <- tibble::tibble(cell_id = 1:40, x_px = seq(0, 390, 10), y_px = 5,
                         phenotype = "Immune")
  for (L in c("A", "B"))                                   # 046: repeats
    .region_csv(file.path(root, "csv", "046", sprintf("046_%s.csv", L)), rows)
  for (k in 1:2)                                           # 052: disjoint
    .region_csv(file.path(root, "csv", "052", sprintf("052_%s.csv", c("A","B")[k])),
                dplyr::mutate(rows, cell_id = cell_id + (k - 1) * 40,
                              x_px = x_px + (k - 1) * 1000))
  for (p in c("046", "052")) dir.create(file.path(root, "annotation", p), recursive = TRUE)

  expect_match(attr(all_slide_overlap_report(all_slide_cells(root)), "verdict"), "MIXED")
})
