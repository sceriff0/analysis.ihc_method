# The arm reader: three phenotyping arms, two file layouts, one set of rules.
# These tests build synthetic trees rather than needing the cluster share, so they
# run anywhere — the geometry tests skip without sf, but the layout, region-keying
# and promotion rules are checked unconditionally, because those decide whether a
# patient is counted at all.
source(here::here("code", "cell_tables.R"))
source(here::here("code", "validation_helpers.R"))
source(here::here("code", "arm_cells.R"))

# A cell table spread over a known box, so a polygon covering half of it keeps
# roughly half the cells. `offset` shifts the ids so two files can be made to hold
# genuinely different cells rather than the same ones.
.cells_csv <- function(path, n, offset = 0L) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(tibble::tibble(
    cell_id           = seq_len(n) + offset,
    x_px              = seq(0, 3999, length.out = n),
    y_px              = seq(0, 2999, length.out = n),
    phenotype         = rep_len(c("PANCK+Tumor", "T cytotoxic", "Immune", "Unknown"), n),
    Out_of_annotation = rep_len(c("False", "False", "False", "True"), n)), path)
}

.poly_geojson <- function(path, x0, y0, x1, y1) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(list(
    type = "Feature",
    geometry = list(type = "Polygon", coordinates = list(list(
      c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0)))),
    properties = structure(list(), names = character(0))), path, auto_unbox = TRUE)
}

.boxes <- list(c(0, 0, 2000, 1500), c(2000, 1500, 4000, 3000), c(0, 1500, 2000, 3000))

.tmp_data <- function() {
  d <- file.path(tempdir(), paste0("armdata-", as.integer(Sys.time()), "-", sample(1e6, 1)))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# --- massimo2: nested, letter-suffixed, no `_all` tier -----------------------
# `regions` maps patient -> region letters; NULL letters means one bare whole-slide
# csv. `annotate` names the patients that get an annotation directory at all — the
# ones that do NOT are this layout's "everything is inside" case.
.massimo2_tree <- function(data_dir,
                           regions  = list(`046` = c("A", "B"), `24086` = NULL),
                           annotate = "046", n = 200) {
  root <- file.path(data_dir, "massimo2")
  for (pid in names(regions)) {
    ls <- regions[[pid]]
    if (is.null(ls)) {
      .cells_csv(file.path(root, "csv", pid, paste0(pid, ".csv")), n)
      next
    }
    # Each region file holds the SAME cells — this layout's region exports repeat
    # the whole slide, so the cohort set must de-duplicate.
    for (L in ls) .cells_csv(file.path(root, "csv", pid, sprintf("%s_%s.csv", pid, L)), n)
    if (pid %in% annotate)
      for (L in ls) do.call(.poly_geojson,
        c(list(file.path(root, "annotation", pid, sprintf("%s_%s.geojson", pid, L))),
          as.list(.boxes[[match(toupper(L), LETTERS)]])))
  }
  root
}

# --- massimo1: flat digit regions + a real whole-slide `_all` tier -----------
# `regions` maps patient -> region indices (integer(0) for a patient with a
# whole-slide export and no selected regions at all — massimo1's 10338 and 15897).
.massimo1_tree <- function(data_dir, regions = list(`046` = 1:2, `10338` = integer(0)),
                           n = 200, n_all = 300, inverted = character(0)) {
  root <- file.path(data_dir, "massimo1")
  for (pid in names(regions)) {
    ks <- regions[[pid]]
    for (k in ks) {
      .cells_csv(file.path(root, "FlowPath_csv_selected", sprintf("%s_a%d.csv", pid, k)), n)
      do.call(.poly_geojson,
        c(list(file.path(root, "annotation_selected", pid, sprintf("%s_a%d.geojson", pid, k))),
          as.list(.boxes[[k]])))
    }
    if (pid %in% inverted)
      for (k in ks)
        .cells_csv(file.path(root, "csv_inverted-classification_modified-thrPANCK",
                             sprintf("%s_a%d.csv", pid, k)), n)
    # The whole-slide tier: a strict superset of the region cells, as a real export is.
    .cells_csv(file.path(root, "FlowPath_csv_all", pid, paste0(pid, ".csv")), n_all)
    do.call(.poly_geojson,
      c(list(file.path(root, "annotation_all", pid, paste0(pid, ".geojson"))),
        list(-1, -1, 4001, 3001)))   # strictly outside the cell extent: st_within is strict
  }
  root
}

.spec2 <- function(d) arm_spec("massimo2", data_dir = d)
.spec1 <- function(d) arm_spec("massimo1", data_dir = d)

# =============================================================================
# massimo2 — the letter layout
# =============================================================================
test_that("massimo2 region letters key to ANNOTATION_<k> by alphabet position", {
  d <- .tmp_data(); .massimo2_tree(d, regions = list(`5456` = c("A", "C")), annotate = "5456")
  cells <- arm_cells(.spec2(d))
  # C is region 3 whether or not B was exported — never by file order.
  expect_setequal(unique(cells$annotation), c("ANNOTATION_1", "ANNOTATION_3"))
})

test_that("a massimo2 patient with no annotation directory is entirely inside", {
  # This layout's OWN stated convention, and the reason it is per-arm rather than
  # global: everywhere else a missing polygon is a reason to drop the patient.
  d <- .tmp_data(); .massimo2_tree(d)
  cells <- arm_cells(.spec2(d))
  bare  <- dplyr::filter(cells, patient_id == "24086")
  expect_gt(nrow(bare), 0)
  expect_equal(unique(bare$annotation), "whole_slide")
  expect_false(any(bare$has_annotation))
  expect_equal(unique(bare$.membership_source), "whole_slide")

  m <- arm_metrics(.spec2(d), cells, NULL, "per_annotation")
  expect_equal(dplyr::filter(m, patient_id == "24086")$source, "whole_slide")
})

test_that("massimo2 region files are keyed per region and never pooled across regions", {
  d <- .tmp_data(); .massimo2_tree(d)
  cells <- dplyr::filter(arm_cells(.spec2(d)), patient_id == "046")
  expect_setequal(unique(cells$annotation), c("ANNOTATION_1", "ANNOTATION_2"))
  # Two files, same cells in each: the raw table is twice the unique count.
  expect_equal(nrow(cells), 400)
})

test_that("massimo2's cohort set de-duplicates cells repeated across region files", {
  # It has no whole-slide tier, so pooling and de-duplicating is the only route —
  # and a naive pool would double every cohort-level denominator while still
  # looking plausible.
  d <- .tmp_data(); .massimo2_tree(d)
  cells <- arm_cells(.spec2(d))
  u     <- arm_cohort_cells(.spec2(d), cells, tibble::tibble())
  expect_equal(nrow(dplyr::filter(u, patient_id == "046")), 200)
  expect_lt(nrow(u), nrow(cells))
})

test_that("massimo2 has no whole-slide tier, so its reconciliation is empty", {
  d <- .tmp_data(); .massimo2_tree(d)
  expect_equal(nrow(arm_wholeslide_reconciliation(.spec2(d))), 0)
})

# =============================================================================
# massimo1 — flat digit regions plus a real whole-slide tier
# =============================================================================
test_that("massimo1 keys its FLAT region csvs from the filename stem", {
  # No directory to trust in this tier, so the stem is the only key. A wrong parse
  # here attributes one patient's cells to another silently.
  d <- .tmp_data(); .massimo1_tree(d, regions = list(`046` = 1:3, `052` = 1:2))
  cells <- arm_cells(.spec1(d))
  expect_setequal(unique(cells$patient_id), c("046", "052"))
  expect_setequal(unique(dplyr::filter(cells, patient_id == "046")$annotation),
                  paste0("ANNOTATION_", 1:3))
  expect_setequal(unique(dplyr::filter(cells, patient_id == "052")$annotation),
                  paste0("ANNOTATION_", 1:2))
})

test_that("massimo1's cohort set is the whole-slide export, NOT the de-duplicated regions", {
  # It is the only arm publishing a real whole-slide csv, so it should never have
  # to reconstruct one. Reconstructing would silently discard the cells that lie
  # outside every selected region.
  d <- .tmp_data(); .massimo1_tree(d, regions = list(`046` = 1:2), n = 200, n_all = 300)
  spec  <- .spec1(d)
  u     <- arm_cohort_cells(spec, arm_cells(spec), arm_union_tier_cells(spec))
  expect_equal(nrow(u), 300)                 # the `_all` file, not the 200 deduped
  expect_equal(unique(u$annotation), "union")
})

test_that("the whole-slide reconciliation reports the deduped union as a subset", {
  # The check that makes arm 1 the ground truth for arms 2 and 3: run their forced
  # procedure on arm 1's regions and it must not exceed arm 1's own export.
  d <- .tmp_data(); .massimo1_tree(d, regions = list(`046` = 1:2), n = 200, n_all = 300)
  spec <- .spec1(d)
  rec  <- arm_wholeslide_reconciliation(spec)
  r046 <- dplyr::filter(rec, patient_id == "046")
  expect_equal(r046$n_wholeslide, 300)
  expect_equal(r046$n_dedup_union, 200)      # two files, same cells, merged to one copy
  expect_match(attr(rec, "verdict"), "as expected")
})

test_that("the reconciliation FLAGS de-duplication that under-merges", {
  # If cell_key_cols() ever stopped identifying cells, the pooled regions would
  # yield MORE cells than the whole-slide export — the one direction that coverage
  # cannot explain, and the one that inflates every arm without an `_all` tier.
  d    <- .tmp_data()
  root <- .massimo1_tree(d, regions = list(`046` = 1:2), n = 200, n_all = 300)
  # Rewrite region 2 with disjoint ids so the two files no longer merge.
  .cells_csv(file.path(root, "FlowPath_csv_selected", "046_a2.csv"), 200, offset = 10000L)
  rec <- arm_wholeslide_reconciliation(.spec1(d))
  expect_equal(dplyr::filter(rec, patient_id == "046")$n_dedup_union, 400)
  expect_match(attr(rec, "verdict"), "UNDER-MERGING")
})

test_that("a patient with a whole-slide export but no regions is promoted to ANNOTATION_1", {
  # massimo1's 10338 and 15897. Dropping them would shrink every union-scoped
  # denominator without saying so; inventing regions would claim polygons nobody drew.
  d <- .tmp_data(); .massimo1_tree(d, regions = list(`046` = 1:2, `10338` = integer(0)))
  spec  <- .spec1(d)
  cells <- arm_cells(spec)
  expect_false("10338" %in% cells$patient_id)      # no region files at all

  ucell <- arm_union_tier_cells(spec)
  polys <- arm_annotations(spec, "region")
  upoly <- arm_annotations(spec, "union")
  out   <- expect_message(arm_promote_unregioned(spec, cells, polys, ucell, upoly),
                          "promoting the union polygon")
  p     <- dplyr::filter(out$cells, patient_id == "10338")
  expect_gt(nrow(p), 0)
  expect_equal(unique(p$annotation), "ANNOTATION_1")
  # 046 has its own regions and must NOT be promoted on top of them.
  expect_setequal(unique(dplyr::filter(out$cells, patient_id == "046")$annotation),
                  c("ANNOTATION_1", "ANNOTATION_2"))
})

test_that("the promotion does not fire for a patient that has SOME of its regions", {
  # Two of three exported is a coverage fact to report, not a reason to substitute
  # the whole slide for the missing region.
  d <- .tmp_data(); .massimo1_tree(d, regions = list(`046` = 1:2))
  spec <- .spec1(d)
  out  <- arm_promote_unregioned(spec, arm_cells(spec), NULL,
                                 arm_union_tier_cells(spec), NULL)
  expect_setequal(unique(out$cells$annotation), c("ANNOTATION_1", "ANNOTATION_2"))
})

# =============================================================================
# massimo1_inverted — arm 1's regions, re-classified
# =============================================================================
test_that("the inverted arm reads its own cells against massimo1's polygons", {
  d <- .tmp_data()
  .massimo1_tree(d, regions = list(`052` = 1:2), inverted = "052")
  spec  <- arm_spec("massimo1_inverted", data_dir = d)
  cells <- arm_cells(spec)
  expect_equal(unique(cells$arm), "massimo1_inverted")
  expect_setequal(unique(cells$annotation), c("ANNOTATION_1", "ANNOTATION_2"))

  skip_if_not_installed("sf")
  polys <- arm_annotations(spec, "region")
  expect_equal(nrow(polys), 2)                      # borrowed from annotation_selected
  expect_setequal(polys$annotation, c("ANNOTATION_1", "ANNOTATION_2"))
})

test_that("the inverted arm has no whole-slide tier and pools its regions instead", {
  d <- .tmp_data()
  .massimo1_tree(d, regions = list(`052` = 1:2), inverted = "052")
  spec <- arm_spec("massimo1_inverted", data_dir = d)
  expect_equal(nrow(arm_union_tier_cells(spec)), 0)
  expect_equal(nrow(arm_cohort_cells(spec)), 200)   # deduped, not 400
})

# =============================================================================
# Shared rules
# =============================================================================
test_that("metrics keep the ihc_annotation_metrics column schema", {
  d <- .tmp_data(); .massimo2_tree(d)
  spec  <- .spec2(d)
  m     <- arm_metrics(spec, arm_cells(spec), NULL, "per_annotation")
  expect_true(all(c("patient_id", "annotation", "source", "n_inside") %in% names(m)))
})

test_that("metrics fall back to the export flag when a region has no polygon", {
  # Never silently count it whole: `source` says which rule was used, always.
  d <- .tmp_data(); .massimo2_tree(d, regions = list(`046` = c("A", "B")), annotate = "046")
  spec <- .spec2(d)
  m    <- arm_metrics(spec, arm_cells(spec), NULL, "per_annotation")
  expect_setequal(unique(dplyr::filter(m, patient_id == "046")$source), "flag")
})

test_that("each region is scored against its OWN polygon", {
  skip_if_not_installed("sf")
  d <- .tmp_data(); .massimo2_tree(d)
  spec  <- .spec2(d)
  cells <- arm_cells(spec)
  m     <- arm_metrics(spec, cells, arm_annotations(spec, "region"), "per_annotation")
  m046  <- dplyr::filter(m, patient_id == "046")
  expect_setequal(m046$source, "sf")
  # The two boxes are disjoint, so neither region can keep all of the slide's cells.
  expect_true(all(m046$n_inside < 200))
})

test_that("massimo1's union comes from its own whole-slide tier, not from dissolving", {
  skip_if_not_installed("sf")
  d <- .tmp_data(); .massimo1_tree(d, regions = list(`046` = 1:2), n = 200, n_all = 300)
  spec <- .spec1(d)
  u <- arm_metrics(spec, arm_cells(spec), arm_annotations(spec, "region"), "union",
                   union_cells = arm_union_tier_cells(spec),
                   union_polys = arm_annotations(spec, "union"))
  # annotation_all covers the whole box, so every cell of the whole-slide export is in.
  expect_equal(dplyr::filter(u, patient_id == "046")$n_inside, 300)
  expect_equal(dplyr::filter(u, patient_id == "046")$source, "sf")
})

test_that("an arm with NO data still yields metrics frames the pages can select from", {
  # THE REGRESSION. A bare tibble::tibble() has zero COLUMNS, so the clinical page's
  # `select(patient_id, annotation, source, n_inside)` did not return nothing — it
  # errored with "Column `patient_id` doesn't exist", three chunks after the real
  # problem and naming neither the arm nor the missing directory. An arm whose tree
  # is not on disk yet is routine, so its frames must flow through the same selects
  # and joins a full arm's do.
  spec <- arm_spec("massimo2", data_dir = file.path(tempdir(), "definitely-not-there"))
  m    <- suppressWarnings(arm_metrics(spec, arm_cells(spec), NULL, "per_annotation"))
  expect_equal(nrow(m), 0)
  expect_gt(ncol(m), 0)
  expect_true(all(c("patient_id", "annotation", "source", "n_inside") %in% names(m)))
  # The exact expression that failed on the cluster.
  expect_silent(dplyr::select(m, patient_id, annotation, source, n_inside))
  # ... and it must still join, so the reconciliation table reports
  # "scored, NOT exported" rather than aborting the knit.
  scored <- tibble::tibble(patient_id = "046", annotation = "ANNOTATION_1", path_pct = 30)
  expect_equal(nrow(dplyr::full_join(scored, m, by = c("patient_id", "annotation"))), 1)
})

test_that("the empty schema is taken from region_ratios_area, so it cannot drift", {
  # Written-out column lists rot. Deriving the empty frame from the real producer
  # means a new metric column appears in both at once.
  full  <- region_ratios_area(tibble::tibble(), 0, 1)
  empty <- arm_empty_metrics()
  expect_setequal(names(empty), c("patient_id", "annotation", "source", names(full)))
})

test_that("a missing region directory warns and yields no cells rather than erroring", {
  spec <- arm_spec("massimo2", data_dir = file.path(tempdir(), "definitely-not-there"))
  expect_warning(out <- arm_cells(spec), "no region csv directory")
  expect_equal(nrow(out), 0)
})

test_that("the inventory reports what was read, per arm, patient and region", {
  d <- .tmp_data(); .massimo2_tree(d)
  inv <- arm_inventory(arm_cells(.spec2(d)))
  expect_true(all(c("arm", "patient_id", "annotation", "n_cells", "membership") %in% names(inv)))
  expect_equal(unique(inv$arm), "massimo2")
  expect_equal(sum(inv$n_cells), 600)   # 046 two region files + 24086 bare
})

test_that("polygons take patient and region from the tree, not the filename", {
  skip_if_not_installed("sf")
  d    <- .tmp_data()
  root <- .massimo2_tree(d)
  file.rename(file.path(root, "annotation", "046", "046_A.geojson"),
              file.path(root, "annotation", "046", "999_A.geojson"))
  polys <- arm_annotations(.spec2(d), "region")
  expect_equal(unique(polys$patient_id), "046")
})

test_that("a union-tier polygon is labelled union, never ANNOTATION_1", {
  # A bare `<pid>.geojson` in an `_all` tree is the patient's dissolved boundary.
  # Reading it as region 1 would score the whole slide against a region nobody drew.
  skip_if_not_installed("sf")
  d <- .tmp_data(); .massimo1_tree(d, regions = list(`046` = 1:2))
  expect_equal(unique(arm_annotations(.spec1(d), "union")$annotation), "union")
})

test_that("the overlap report says whether region files repeat cells", {
  d <- .tmp_data(); .massimo2_tree(d)
  rep <- arm_overlap_report(arm_cells(.spec2(d)))
  expect_match(attr(rep, "verdict"), "repeat cells")
  expect_equal(dplyr::filter(rep, patient_id == "046")$rows_per_cell, 2)
})

# =============================================================================
# Per-cell membership — the flag a re-cut denominator needs
# =============================================================================
# arm_metrics() asks the polygon which cells are inside and then aggregates the
# answer away. The bulk-RNA comparison re-cuts the denominator (tumour cells inside
# the annotation), which no aggregate can supply after the fact, so the per-cell
# verdict is exposed. These tests pin WHICH polygon it scores against and what it
# does when there is none, because both silently change every restricted fraction.

test_that("in-annotation membership scores against the arm's UNION polygon", {
  skip_if_not_installed("sf")
  d <- .tmp_data(); .massimo1_tree(d, regions = list(`046` = 1:2), n = 200, n_all = 300)
  spec <- .spec1(d)
  sc   <- arm_cells_in_annotation(spec, arm_cohort_cells(spec))
  # annotation_all covers the whole cell extent, so every whole-slide cell is inside —
  # and it must be the UNION polygon that says so, not the two `_selected` boxes,
  # which between them cover only part of the slide.
  expect_equal(nrow(sc), 300)
  expect_true(all(sc$in_annotation))
  expect_equal(unique(sc$.in_annotation_source), "sf")
})

test_that("a smaller union polygon keeps only the cells inside it", {
  skip_if_not_installed("sf")
  d <- .tmp_data(); .massimo1_tree(d, regions = list(`046` = 1L), n = 100, n_all = 300)
  # Shrink annotation_all to the lower-left quarter of the cell extent. The cells run
  # along the box diagonal, so a quarter-box keeps roughly a quarter of them — the
  # point is that the count is strictly between 0 and all, i.e. geometry ran.
  .poly_geojson(file.path(d, "massimo1", "annotation_all", "046", "046.geojson"),
                -1, -1, 1000, 750)
  spec <- .spec1(d)
  sc   <- arm_cells_in_annotation(spec, arm_cohort_cells(spec))
  expect_equal(unique(sc$.in_annotation_source), "sf")
  expect_gt(sum(sc$in_annotation), 0)
  expect_lt(sum(sc$in_annotation), nrow(sc))
})

test_that("massimo2 dissolves its region polygons, having no union tier", {
  skip_if_not_installed("sf")
  d <- .tmp_data()
  # 046 gets regions A and B, whose boxes are the lower-left and upper-right
  # quarters — so the DISSOLVED pair keeps more cells than either alone and still
  # fewer than the whole slide.
  .massimo2_tree(d, regions = list(`046` = c("A", "B")), annotate = "046", n = 200)
  spec <- .spec2(d)
  sc   <- arm_cells_in_annotation(spec, arm_cohort_cells(spec))
  expect_equal(unique(sc$.in_annotation_source), "sf")
  expect_gt(sum(sc$in_annotation), 0)
  expect_lt(sum(sc$in_annotation), nrow(sc))
})

test_that("a massimo2 patient with no annotation directory is entirely in-annotation", {
  # The arm's OWN stated convention (`bare_region_is = "whole_slide"`), reached only
  # because 24086 has no polygon AND its export carries no usable flag path here.
  # The source column has to say `whole_slide`, not `sf`: a 100%-inside patient under
  # a convention and one under a polygon look identical in the fraction.
  d <- .tmp_data(); .massimo2_tree(d, regions = list(`24086` = NULL), annotate = character(0))
  spec <- .spec2(d)
  sc   <- suppressWarnings(arm_cells_in_annotation(spec, arm_cohort_cells(spec)))
  expect_true(all(sc$in_annotation %in% TRUE))
  expect_equal(unique(sc$.in_annotation_source), "whole_slide")
  # THE REGRESSION. The synthetic export also carries Out_of_annotation, which flags
  # one cell in four as outside — so an implementation that reached for the flag
  # first would keep 75% of this patient's cells while arm_metrics(), which never
  # consults the flag for an unannotated patient, keeps 100%. Two numbers for one
  # patient on one page, both plausible.
  expect_equal(sum(sc$in_annotation), nrow(sc))
  expect_lt(sum(!cell_outside(sc)), nrow(sc))
})

test_that("with no polygon the export's own flag decides, and says so", {
  # A massimo1 patient whose annotation_all geojson is unreadable. The arm has no
  # "count it whole" convention, so the only remaining source is the flag.
  d <- .tmp_data(); .massimo1_tree(d, regions = list(`046` = 1L), n = 100, n_all = 200)
  writeLines("not json", file.path(d, "massimo1", "annotation_all", "046", "046.geojson"))
  writeLines("not json", file.path(d, "massimo1", "annotation_selected", "046", "046_a1.geojson"))
  spec <- .spec1(d)
  sc   <- suppressWarnings(arm_cells_in_annotation(spec, arm_cohort_cells(spec)))
  expect_equal(unique(sc$.in_annotation_source), "flag")
  # .cells_csv flags one cell in every four as outside.
  expect_equal(sum(sc$in_annotation), sum(!cell_outside(sc)))
})

test_that("the inventory separates the tumour subset from the annotation subset", {
  skip_if_not_installed("sf")
  d <- .tmp_data(); .massimo1_tree(d, regions = list(`046` = 1L), n = 100, n_all = 200)
  spec <- .spec1(d)
  inv  <- arm_in_annotation_inventory(arm_cells_in_annotation(spec, arm_cohort_cells(spec)))
  expect_equal(inv$n_cells, 200)
  expect_equal(inv$n_in_annotation, 200)          # annotation_all covers the extent
  expect_equal(inv$n_tumor, 50)                   # one cell in four is PANCK+Tumor
  expect_equal(inv$n_tumor_in_annotation, 50)
  expect_equal(inv$pct_in_annotation, 100)
})

test_that("the tumour subset is taken from cell_lineage, not from the label text", {
  # The two phenotypers spell the same population differently, so a grepl("Tumor")
  # here would agree with cell_lineage() by accident on FlowPath and disagree the
  # moment a panel renames anything. Pin the accessor.
  d <- .tmp_data(); .massimo1_tree(d, regions = list(`046` = 1L), n = 100, n_all = 200)
  spec  <- .spec1(d)
  cells <- arm_cohort_cells(spec)
  expect_equal(sum(cell_lineage(cell_phenotype(cells)) %in% "Tumor"), 50)
})
