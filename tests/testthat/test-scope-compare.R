# ONE quantity, THREE nested scopes. These tests pin the two things a reader relies
# on and cannot check by eye: that the scopes are built from the right cell sets, and
# that a patient missing a scope is ABSENT from it rather than silently substituted.
source(here::here("code", "cell_tables.R"))
source(here::here("code", "validation_helpers.R"))
source(here::here("code", "aggregation_compare.R"))
source(here::here("code", "scope_compare.R"))

# A minimal cell table: `frac` of the cells are tumour.
.cells <- function(pid, n, frac) {
  tibble::tibble(
    patient_id      = pid,
    cell_id         = seq_len(n),
    phenotype       = c(rep("PANCK+Tumor", round(n * frac)),
                        rep("T cytotoxic", n - round(n * frac))),
    phenotype_clean = c(rep("PANCK+Tumor", round(n * frac)),
                        rep("T cytotoxic", n - round(n * frac))))
}

# A metrics-shaped row, as arm_metrics() emits.
.metrics <- function(pid, ann, frac, n = 100, source = "sf") {
  region_ratios_area(.cells(pid, n, frac), 0, 1) |>
    dplyr::mutate(patient_id = pid, annotation = ann, source = source, .before = 1)
}

test_that("the whole-slide scope counts every cell and claims no area", {
  # No polygon is consulted, so there IS no area: reporting a density would require
  # inventing a denominator.
  cells <- dplyr::bind_rows(.cells("046", 100, 0.4), .cells("052", 200, 0.6))
  ws <- arm_wholeslide_ratios(cells)
  expect_equal(nrow(ws), 2)
  expect_equal(unique(ws$annotation), "whole_slide")
  expect_equal(unique(ws$source), "whole_slide")
  expect_true(all(is.na(ws$area_mm2)))
  expect_equal(ws$n_inside, c(100, 200))
  expect_equal(round(ws$tumor_over_inside, 3), c(0.4, 0.6))
})

test_that("scope_table stacks the three scopes with the named palette's levels", {
  mem <- list(
    cells = dplyr::bind_rows(.cells("046", 100, 0.30), .cells("052", 100, 0.50)),
    union = dplyr::bind_rows(.metrics("046", "union", 0.40),
                             .metrics("052", "union", 0.55)),
    per_annotation = dplyr::bind_rows(.metrics("046", "ANNOTATION_1", 0.45),
                                      .metrics("046", "ANNOTATION_2", 0.35),
                                      .metrics("052", "ANNOTATION_1", 0.55)))
  st <- scope_table(mem)
  expect_setequal(levels(st$scope), names(SCOPE_COLS))
  expect_equal(as.integer(table(st$scope)), c(2L, 2L, 3L))
  # The scopes NEST, so the whole-slide row must exist for every patient that has
  # any cells at all — it is the one scope no annotation can remove.
  expect_setequal(dplyr::filter(st, scope == "whole_slide")$patient_id, c("046", "052"))
})

# A patient with no annotation directory (massimo2's 24086) gets a union row
# labelled `whole_slide`. Whether that counts as an ANNOTATED scope is the arm's
# own call, stated in the registry as `bare_region_is` — so the next two tests are
# the same fixture read under two arms, and they must disagree.
.mem_bare <- function(mode = NULL) {
  m <- list(
    cells = dplyr::bind_rows(.cells("046", 100, 0.3), .cells("24086", 100, 0.7)),
    union = dplyr::bind_rows(.metrics("046", "union", 0.4),
                             .metrics("24086", "whole_slide", 0.7, source = "whole_slide")),
    per_annotation = .metrics("046", "ANNOTATION_1", 0.45))
  if (!is.null(mode)) m$mode <- mode
  m
}

test_that("massimo2's bare-region patient IS its own annotation_all", {
  # arm 2 declares bare_region_is = "whole_slide": an unsuffixed region file means
  # the whole slide IS the annotated region. So 24086's annotation_all and its
  # whole_slide are the same number BY CONSTRUCTION, and it belongs in both scopes.
  st  <- scope_table(.mem_bare("massimo2"))
  s24 <- dplyr::filter(st, patient_id == "24086")
  expect_setequal(as.character(s24$scope), c("whole_slide", "annotation_all"))
  # Same number in both — not a second measurement, the same union row carried
  # through. If these ever differ, annotation_all stopped reading the arm's row.
  expect_equal(dplyr::n_distinct(s24$value), 1L)
  # The provenance survives the promotion: the scope says annotation_all, `source`
  # still says membership was decided by the arm's convention, not by a polygon.
  expect_equal(unique(dplyr::filter(s24, scope == "annotation_all")$source),
               "whole_slide")
  # It has no ANNOTATION_k row, and gets none invented for it.
  expect_false("annotation_k" %in% as.character(s24$scope))
})

test_that("an arm with no bare-region convention leaves such a patient out", {
  # massimo1 does NOT declare bare_region_is, so there a `whole_slide` union row is
  # a missing polygon -- a geojson that failed to parse, a tree linked to the wrong
  # root -- not a stated convention. Promoting it would turn that data problem into
  # a silently 100 %-inside patient sitting on the x = y line looking like a result.
  for (mode in list("massimo1", NULL)) {   # NULL = a hand-built mem / the mirage mode
    st  <- scope_table(.mem_bare(mode))
    s24 <- dplyr::filter(st, patient_id == "24086")
    expect_equal(nrow(s24), 1)
    expect_equal(as.character(s24$scope), "whole_slide")
  }
})

test_that("one-vs-many reports whether the single value sits inside the region span", {
  # The comparison a correlation cannot show, and the reason the range figure exists.
  mem <- list(
    cells = dplyr::bind_rows(.cells("in", 100, 0.45), .cells("out", 100, 0.10)),
    union = dplyr::bind_rows(.metrics("in", "union", 0.45), .metrics("out", "union", 0.50)),
    per_annotation = dplyr::bind_rows(
      .metrics("in",  "ANNOTATION_1", 0.40), .metrics("in",  "ANNOTATION_2", 0.50),
      .metrics("out", "ANNOTATION_1", 0.40), .metrics("out", "ANNOTATION_2", 0.60)))
  om <- scope_one_vs_many(scope_table(mem))
  expect_true(dplyr::filter(om, patient_id == "in")$inside_range)
  expect_false(dplyr::filter(om, patient_id == "out")$inside_range)
  # gap is signed: whole slide MINUS union, so a negative gap means the annotation is
  # tumour-richer than the rest of the slide — the direction a pathologist's drawing
  # normally produces.
  expect_lt(dplyr::filter(om, patient_id == "out")$gap_to_union, 0)
})

test_that("a single-annotation patient gets a degenerate range, not a dropped row", {
  # 10338 and 15897 have one region each. min == max is exact, not missing — and the
  # patient must still appear, which is the whole point of including them.
  mem <- list(
    cells = .cells("10338", 100, 0.30),
    union = .metrics("10338", "union", 0.52),
    per_annotation = .metrics("10338", "ANNOTATION_1", 0.52))
  om <- scope_one_vs_many(scope_table(mem))
  expect_equal(nrow(om), 1)
  expect_equal(om$n_ann, 1)
  expect_equal(om$ann_min, om$ann_max)
  expect_false(om$inside_range)          # 0.30 is outside [0.52, 0.52]
})

test_that("the aggregators are the project's existing vocabulary, not a second one", {
  # A reader should meet ONE set of ways to collapse annotations across the site.
  mem <- list(
    cells = dplyr::bind_rows(.cells("046", 100, 0.3), .cells("052", 100, 0.5)),
    union = dplyr::bind_rows(.metrics("046", "union", 0.4), .metrics("052", "union", 0.55)),
    per_annotation = dplyr::bind_rows(
      .metrics("046", "ANNOTATION_1", 0.45), .metrics("046", "ANNOTATION_2", 0.35),
      .metrics("052", "ANNOTATION_1", 0.55), .metrics("052", "ANNOTATION_2", 0.50)))
  pairs <- scope_aggregator_pairs(scope_table(mem), mem)
  expect_setequal(unique(pairs$ihc_agg), names(IHC_AGG_LABELS))
  expect_true(all(c("whole_slide", "ihc_val") %in% names(pairs)))

  st <- scope_aggregator_stats(pairs)
  expect_true(all(c("spearman", "bias") %in% names(st)))
  # bias is signed and is the point: an aggregator can rank patients perfectly and
  # still read high on every one of them.
  expect_equal(nrow(st), length(IHC_AGG_LABELS))
})

test_that("every figure degrades to invisible() on empty input rather than erroring", {
  # These run on pages whose arm may have no data at all; a missing tree must not
  # abort the knit three chunks later.
  empty <- tibble::tibble()
  expect_null(plot_scope_scatter(tibble::tibble(whole_slide = numeric(0),
                                                annotation_all = numeric(0),
                                                patient_id = character(0))))
  expect_null(plot_scope_range(tibble::tibble(whole_slide = numeric(0),
                                              ann_min = numeric(0), ann_max = numeric(0),
                                              ann_median = numeric(0),
                                              annotation_all = numeric(0),
                                              patient_id = character(0)), empty))
  expect_null(plot_scope_aggregators(empty))
  expect_equal(nrow(scope_table(list(cells = empty, union = empty, per_annotation = empty))), 0)
  expect_equal(nrow(arm_wholeslide_ratios(empty)), 0)
})

test_that("the scope palette is NAMED, so a level keeps its colour when one is absent", {
  # An unnamed scale assigns by factor POSITION: annotation_all would be blue in a
  # panel containing all three scopes and black in one that omits whole_slide.
  expect_setequal(names(SCOPE_COLS), c("whole_slide", "annotation_all", "annotation_k"))
  expect_equal(levels(scope_factor(c("annotation_k", "whole_slide"))), names(SCOPE_COLS))
  expect_setequal(names(SCOPE_LABELS), names(SCOPE_COLS))
})
