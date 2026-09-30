# The arm registry: which files belong to which arm, and how a filename names its
# region. Pure path logic — no data on disk, no sf — so it runs anywhere and is the
# first thing to fail if a producer renames a directory.
source(here::here("code", "arms.R"))

test_that("every arm declares a region tier and resolves to a root", {
  for (a in ARM_MODES) {
    spec <- arm_spec(a, data_dir = "/data")
    expect_equal(spec$arm, a)
    expect_false(is.null(spec$region_csv), info = a)
    expect_false(is.null(spec$region_poly), info = a)
    expect_equal(spec$region_csv$path, file.path("/data", spec$root, spec$region_csv$dir))
  }
})

test_that("the inverted arm reads massimo1's polygons, not its own", {
  # It re-classifies arm 1's regions; it does not redraw them. If this ever
  # diverged, the inverted arm's cells would be scored against polygons nobody
  # drew for them and the disagreement would read as a threshold effect.
  m1  <- arm_spec("massimo1", data_dir = "/data")
  inv <- arm_spec("massimo1_inverted", data_dir = "/data")
  expect_equal(inv$root, m1$root)
  expect_equal(inv$region_poly$path, m1$region_poly$path)
  expect_equal(inv$union_poly$path,  m1$union_poly$path)
  # ... but its CELLS are its own.
  expect_false(identical(inv$region_csv$path, m1$region_csv$path))
})

test_that("massimo2 has no whole-slide tier; both massimo1 arms do", {
  # arm 1's `_all` export is what makes it the ground truth for the de-duplication
  # massimo2 is forced to use. The inverted arm shares that export, but only for the
  # patients it did not re-run — see the fallback test below.
  expect_true(arm_has_union_tier(arm_spec("massimo1", data_dir = "/data")))
  expect_true(arm_has_union_tier(arm_spec("massimo1_inverted", data_dir = "/data")))
  expect_false(arm_has_union_tier(arm_spec("massimo2", data_dir = "/data")))
})

test_that("only the inverted arm falls back to another tier", {
  # The modified-PANCK run covers 052 and 5456 only. Reading just those would make
  # every cohort number a statement about two patients while looking like a statement
  # about the study, so the other four load their ordinary massimo1 cells.
  inv <- arm_spec("massimo1_inverted", data_dir = "/data")
  expect_true(arm_has_fallback(inv))
  expect_equal(inv$region_csv_fallback$path,
               arm_spec("massimo1", data_dir = "/data")$region_csv$path)
  # The two tiers must be distinguishable in the output, or a reader cannot tell
  # which patients were actually re-classified.
  expect_equal(inv$region_csv$classification, "inverted")
  expect_equal(inv$region_csv_fallback$classification, "normal")
  for (a in c("massimo1", "massimo2"))
    expect_false(arm_has_fallback(arm_spec(a, data_dir = "/data")))
})

test_that("letters map to region index by alphabet position, not file order", {
  expect_equal(arm_letter_index("A"), 1L)
  expect_equal(arm_letter_index("c"), 3L)   # C is region 3 whether or not B exists
  expect_true(is.na(arm_letter_index("?")))
})

test_that("each naming pattern keys the patient and the region it declares", {
  p <- function(path, pat) arm_parse_name(path, pat)

  expect_equal(p("csv/046/046_A.csv", "nested_letter"),
               list(patient = "046", annotation = "ANNOTATION_1"))
  expect_equal(p("csv/5456/5456_C.csv", "nested_letter")$annotation, "ANNOTATION_3")
  expect_equal(p("FlowPath_csv_selected/24086_a3.csv", "flat_digit"),
               list(patient = "24086", annotation = "ANNOTATION_3"))
  expect_equal(p("annotation_selected/052/052_a2.geojson", "nested_digit"),
               list(patient = "052", annotation = "ANNOTATION_2"))
})

test_that("the directory name wins over the filename stem in every nested pattern", {
  # A file renamed by hand would otherwise be keyed to the wrong patient silently,
  # and the only symptom is an implausible in-annotation count.
  expect_equal(arm_parse_name("csv/046/999_A.csv", "nested_letter")$patient, "046")
  expect_equal(arm_parse_name("annotation_selected/052/xx_a2.geojson", "nested_digit")$patient, "052")
  expect_equal(arm_parse_name("annotation_all/10338/whatever.geojson", "nested_bare")$patient, "10338")
  # flat_digit has no directory to trust, so the stem is the only key there.
  expect_equal(arm_parse_name("FlowPath_csv_selected/046_a1.csv", "flat_digit")$patient, "046")
})

test_that("a file with no region suffix comes back with annotation NA, not a guessed region", {
  # The CALLER decides what a bare file means, because it differs by tier: a bare
  # csv in massimo2's region tier is a whole-slide export, while a bare geojson in
  # massimo1's `_all` tier is that patient's union polygon. Deciding here would
  # force one reading onto the other.
  expect_true(is.na(arm_parse_name("csv/24086/24086.csv", "nested_letter")$annotation))
  expect_true(is.na(arm_parse_name("annotation_all/10338/10338.geojson", "nested_bare")$annotation))
})

test_that("one arm's parser never silently claims another arm's filenames", {
  # `046_a1.csv` ends in a digit, so the LETTER parser must fall through to the bare
  # reading rather than inventing a region. Pointing a tree at the wrong spec has to
  # produce a visible "no region suffix" path, not a plausible wrong answer.
  expect_true(is.na(arm_parse_name("csv/046/046_a1.csv", "nested_letter")$annotation))
  # ... and symmetrically, a letter-suffixed file handed to the digit parser.
  expect_true(is.na(arm_parse_name("annotation_selected/046/046_B.geojson", "nested_digit")$annotation))
})

test_that("flat_letter exists for the legacy tree even though no arm uses it", {
  # The axis is nested/flat x letter/digit/bare. A registry describing five of six
  # cases invites the sixth to be re-implemented elsewhere — which is how
  # validation_helpers.R ended up with a second copy of this parser.
  expect_equal(arm_parse_name("ann/046_B.geojson", "flat_letter"),
               list(patient = "046", annotation = "ANNOTATION_2"))
  expect_true(is.na(arm_parse_name("ann/046_a2.geojson", "flat_letter")$annotation))
})

test_that(".annotation_key delegates here rather than keeping its own parser", {
  source(here::here("code", "cell_tables.R"))
  source(here::here("code", "validation_helpers.R"))
  root <- "/tmp/tree"
  # A bare file is the one thing the legacy reader adds on top: it means ANNOTATION_1
  # in a tree that belongs to no arm, whereas an ARM's bare geojson is a union polygon.
  expect_equal(.annotation_key(file.path(root, "24086", "24086.csv"), root)$annotation,
               "ANNOTATION_1")
  expect_true(is.na(arm_parse_name(file.path(root, "24086", "24086.csv"),
                                   "nested_bare")$annotation))
})

test_that("arm_tier_status reports which tiers are actually on disk", {
  spec <- arm_spec("massimo2", data_dir = file.path(tempdir(), "definitely-not-there"))
  st   <- arm_tier_status(spec)
  expect_equal(nrow(st), 6)   # + region_csv_fallback, flag_csv
  expect_false(any(st$exists))
  expect_true(all(is.na(st$dir[st$tier %in% c("union_csv", "union_poly")])))
})

# =============================================================================
# csv-only patients: membership from the export's own flag
# =============================================================================
test_that("the csv-only patients are recognised under every spelling slide_key() would give them", {
  expect_true(all(arm_flag_patient(FLAG_MEMBERSHIP_PATIENTS)))
  expect_equal(FLAG_MEMBERSHIP_PATIENTS, c("15370", "36672_2"))
  # The same normalisation as norm_slide_id(): digits only, so a file named
  # "36672_2", a clinical "36672-2" and the key "366722" are one patient.
  expect_true(arm_flag_patient("36672-2"))
  expect_true(arm_flag_patient("366722"))
  expect_true(arm_flag_patient("EPM - 15370"))
  # And nobody else — 24086 keeps massimo2's whole-slide convention. 153070 is
  # 15370's FILE name, not an id: the registry maps it, nothing parses it.
  expect_false(any(arm_flag_patient(c("046", "24086", "36672", "1537", "153070", "36766_2"))))
  expect_length(arm_flag_patient(character(0)), 0)
})

test_that("every arm points the csv-only patients at the same data/flowpath files", {
  d <- file.path(tempdir(), paste0("flagdata-", sample(1e6, 1)))
  dir.create(file.path(d, "flowpath"), recursive = TRUE)
  file.create(file.path(d, "flowpath", "153070.csv"))
  for (arm in ARM_MODES) {
    ff <- arm_flag_files(arm_spec(arm, data_dir = d))
    expect_equal(ff$patient_id, c("15370", "36672_2"))
    expect_equal(basename(ff$path), c("153070.csv", "36672_2.csv"))
    expect_equal(ff$exists, c(TRUE, FALSE))
  }
})
