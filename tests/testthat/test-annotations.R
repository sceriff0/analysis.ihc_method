# The annotation loader reads TWO layouts, because the producer changed shape and both
# have to load: the current nested tree (<patient>/<patient>_<LETTER>.geojson) and the
# flat legacy one (<patient>_a<k>.geojson). Four pages consume `annots` and none of
# them should have to know which is on disk.
source(here::here("code", "cell_tables.R"))
source(here::here("code", "validation_helpers.R"))

.poly <- function(path, x0, y0, x1, y1) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(list(
    type = "Feature",
    geometry = list(type = "Polygon", coordinates = list(list(
      c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0)))),
    properties = structure(list(), names = character(0))), path, auto_unbox = TRUE)
}

nested_tree <- function() {
  root <- file.path(tempdir(), paste0("ann-nested-", sample(1e6, 1)))
  .poly(file.path(root, "046", "046_A.geojson"), 0, 0, 2000, 1500)
  .poly(file.path(root, "046", "046_C.geojson"), 0, 1500, 2000, 3000)   # B skipped
  .poly(file.path(root, "052", "052_A.geojson"), 0, 0, 1000, 1000)
  root
}

flat_tree <- function() {
  root <- file.path(tempdir(), paste0("ann-flat-", sample(1e6, 1)))
  .poly(file.path(root, "046_a1.geojson"), 0, 0, 2000, 1500)
  .poly(file.path(root, "046_a2.geojson"), 2000, 1500, 4000, 3000)
  .poly(file.path(root, "052.geojson"), 0, 0, 1000, 1000)   # bare = sole ANNOTATION_1
  root
}

test_that("region letters map by alphabet position, not by file order", {
  # C stays region 3 even though B was never exported — otherwise a patient's regions
  # would silently renumber and stop lining up with neoplastic_data's ANNOTATION_1..3.
  expect_equal(annotation_from_letter("A"), "ANNOTATION_1")
  expect_equal(annotation_from_letter("c"), "ANNOTATION_3")
  expect_true(is.na(annotation_from_letter("1")))
})

test_that(".annotation_key reads both layouts without needing sf", {
  # The key parser is where the layout logic lives, and it is pure string work — so it
  # is tested directly rather than only through the sf-gated loader above. This is also
  # the parser annotation_membership_qc() uses on the cell csvs, which is what makes a
  # csv and the polygon it is compared against agree on which region they are.
  root <- "/tmp/tree"

  k <- .annotation_key(file.path(root, "046", "046_A.geojson"), root)
  expect_equal(k$patient, "046"); expect_equal(k$annotation, "ANNOTATION_1")

  k <- .annotation_key(file.path(root, "046", "046_C.csv"), root)   # csvs too
  expect_equal(k$patient, "046"); expect_equal(k$annotation, "ANNOTATION_3")

  # Nested: the DIRECTORY is the patient of record, whatever the stem says.
  k <- .annotation_key(file.path(root, "052", "renamed_B.geojson"), root)
  expect_equal(k$patient, "052"); expect_equal(k$annotation, "ANNOTATION_2")

  # Flat legacy: the stem carries the patient, and the suffix is a DIGIT.
  k <- .annotation_key(file.path(root, "046_a2.geojson"), root)
  expect_equal(k$patient, "046"); expect_equal(k$annotation, "ANNOTATION_2")

  # Bare name, flat: that patient's sole annotation.
  k <- .annotation_key(file.path(root, "24086.csv"), root)
  expect_equal(k$patient, "24086"); expect_equal(k$annotation, "ANNOTATION_1")

  # Bare name, nested: the directory names the patient; still region 1.
  k <- .annotation_key(file.path(root, "24086", "24086.csv"), root)
  expect_equal(k$patient, "24086"); expect_equal(k$annotation, "ANNOTATION_1")
})

test_that("the nested layout keys patient from the directory and region from the letter", {
  skip_if_not_installed("sf")
  a <- load_annotations(nested_tree())
  expect_equal(nrow(a), 3)
  expect_setequal(a$annotation[a$patient_id == "046"], c("ANNOTATION_1", "ANNOTATION_3"))
  expect_equal(a$annotation[a$patient_id == "052"], "ANNOTATION_1")
})

test_that("the flat legacy layout still loads, digits and bare names both", {
  skip_if_not_installed("sf")
  a <- load_annotations(flat_tree())
  expect_equal(nrow(a), 3)
  expect_setequal(a$annotation[a$patient_id == "046"], c("ANNOTATION_1", "ANNOTATION_2"))
  # A bare <patient>.geojson is that patient's sole annotation, and in neoplastic_data
  # such patients carry tumour content only in ANNOTATION_1.
  expect_equal(a$annotation[a$patient_id == "052"], "ANNOTATION_1")
})

test_that("the directory name wins when the filename stem disagrees", {
  skip_if_not_installed("sf")
  # A hand-renamed file would otherwise be keyed to the wrong patient, silently, and
  # the only symptom would be an implausible in-annotation count.
  root <- file.path(tempdir(), paste0("ann-mismatch-", sample(1e6, 1)))
  .poly(file.path(root, "046", "WRONGNAME_A.geojson"), 0, 0, 100, 100)
  a <- load_annotations(root)
  expect_equal(a$patient_id, "046")
  expect_equal(a$annotation, "ANNOTATION_1")
})

test_that("patient_ids filters, and an empty result is NULL rather than an error", {
  skip_if_not_installed("sf")
  root <- nested_tree()
  expect_equal(unique(load_annotations(root, patient_ids = "046")$patient_id), "046")
  expect_null(load_annotations(root, patient_ids = "99999"))
})

test_that("a missing or empty annotation directory returns NULL, not an error", {
  skip_if_not_installed("sf")
  # A cohort can legitimately be all-whole-slide, so this must not stop a report at
  # the loader.
  expect_null(load_annotations(file.path(tempdir(), "definitely-absent")))
  empty <- file.path(tempdir(), paste0("ann-empty-", sample(1e6, 1)))
  dir.create(empty, recursive = TRUE)
  expect_null(load_annotations(empty))
})

test_that("only two membership modes survive", {
  source(here::here("code", "membership.R"))
  expect_setequal(MEMBERSHIP_MODES, c("all_slide", "mirage"))
  # The removed modes must not be silently accepted by match.arg's partial matching.
  expect_error(membership_data("flag", tibble::tibble()), "should be one of")
  expect_error(membership_data("geojson", tibble::tibble()), "should be one of")
})

test_that("all_slide refuses an externally supplied annotation set", {
  source(here::here("code", "membership.R"))
  # Its polygons are paired to its cell files region by region, so an outside set
  # would not line up. Erroring beats ignoring the argument.
  expect_error(membership_data("all_slide", tibble::tibble(), annots = "anything"),
               "reads its own polygons")
})
