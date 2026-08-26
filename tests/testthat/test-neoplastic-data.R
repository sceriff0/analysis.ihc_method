# The pathologist's neoplastic-cellularity reference — ONE TABLE PER ARM. These tests
# pin the SHAPE and the join key, not the values: a value typo is caught by the
# reconciliation table the clinical page prints, but a shape or vocabulary change
# breaks the concordance join silently, dropping pairs while still producing a
# plausible correlation.
source(here::here("code", "cell_tables.R"))
source(here::here("code", "validation_helpers.R"))

# Read the tribbles out of load_data.R without executing the file: load_data.R eagerly
# reads counts.RData and clinical_data.xlsx, which are off-repo.
.neo <- function(name) {
  src  <- readLines(here::here("code", "load_data.R"))
  from <- grep(paste0("^", name, " <- tribble\\("), src)
  expect_length(from, 1)
  to   <- from - 1 + grep("^\\)", src[from:length(src)])[1]
  e <- new.env()
  e$tribble <- tibble::tribble
  eval(parse(text = paste(src[from:to], collapse = "\n")), envir = e)
  e[[name]]
}

m1 <- .neo("neoplastic_massimo1")
m2 <- .neo("neoplastic_massimo2")

test_that("both arms are LONG: one row per scored region, not one per patient", {
  # Wide (~ANNOTATION_1..3) cannot express a whole-slide score without claiming a
  # region the pathologist never drew.
  for (nm in c("m1", "m2")) {
    neo <- get(nm)
    expect_setequal(names(neo), c("SAMPLE", "annotation", "path_pct"))
    expect_gt(nrow(neo), dplyr::n_distinct(neo$SAMPLE))
  }
})

test_that("annotation labels are the vocabulary the metrics frames emit", {
  # The clinical pages inner-join on (patient_id, annotation). A label outside this
  # vocabulary silently drops that region from every correlation.
  ok <- c(paste0("ANNOTATION_", 1:3), "whole_slide")
  for (nm in c("m1", "m2")) {
    neo <- get(nm)
    expect_true(all(neo$annotation %in% ok),
                info = paste(nm, "unexpected:",
                             paste(setdiff(neo$annotation, ok), collapse = ", ")))
  }
})

test_that("region numbering is contiguous from 1 per patient, in both arms", {
  # A patient with an ANNOTATION_3 and no ANNOTATION_2 means a region went missing
  # on one side of the join.
  for (nm in c("m1", "m2")) {
    by_pt <- get(nm) |>
      dplyr::filter(annotation != "whole_slide") |>
      dplyr::mutate(k = as.integer(sub("^ANNOTATION_", "", annotation))) |>
      dplyr::group_by(SAMPLE) |>
      dplyr::summarise(ks = list(sort(k)), .groups = "drop")
    for (i in seq_len(nrow(by_pt)))
      expect_equal(by_pt$ks[[i]], seq_along(by_pt$ks[[i]]),
                   info = paste(nm, "patient", by_pt$SAMPLE[i]))
  }
})

test_that("a whole-slide patient has exactly one score and no regions", {
  for (nm in c("m1", "m2")) {
    neo <- get(nm)
    ws  <- dplyr::filter(neo, annotation == "whole_slide")
    for (pid in unique(ws$SAMPLE))
      expect_equal(nrow(dplyr::filter(neo, SAMPLE == pid)), 1,
                   info = paste(nm, pid, "should have one whole-slide score only"))
  }
})

test_that("slide_key normalises every SAMPLE to the export's patient id", {
  # The join is on slide_key(SAMPLE); if that does not equal the csv directory name
  # the patient never matches.
  for (nm in c("m1", "m2")) expect_equal(slide_key(get(nm)$SAMPLE), get(nm)$SAMPLE)
})

# --- Arm 2: the scored table -------------------------------------------------
test_that("arm 2's percentages are percentages, not fractions", {
  # path_frac = path_pct / 100 downstream; a value already in 0..1 would come out
  # 100x too small and read as near-zero tumour content.
  expect_true(all(m2$path_pct >= 1 & m2$path_pct <= 100))
})

test_that("arm 2 keeps 24086 as a whole-slide score", {
  # It has a bare csv and no annotation directory in that export, so `whole_slide`
  # is the only label its metrics row can carry.
  expect_equal(dplyr::filter(m2, SAMPLE == "24086")$annotation, "whole_slide")
})

# --- Arm 1: read from thr_head&neck.xlsx -------------------------------------
test_that("arm 1's scored regions are percentages, and its unscored ones stay NA", {
  # Filled 2026-08-26 from data/Massimo1/thr_head&neck.xlsx (one sheet per patient,
  # each with a `neoplastic cellularity (%)` block keyed annotation_1..3).
  scored <- m1$path_pct[!is.na(m1$path_pct)]
  expect_equal(length(scored), 11)                    # the 11 `_selected` regions
  expect_true(all(scored >= 1 & scored <= 100))       # percentages, not fractions
  # 10338 and 15897 have no `_selected` regions and no cellularity block, so their
  # promoted ANNOTATION_1 is exported-but-unscored. NA keeps the pair out of the
  # correlation instead of inventing a plausible number that would plot.
  expect_setequal(m1$SAMPLE[is.na(m1$path_pct)], c("10338", "15897"))
  expect_type(m1$path_pct, "double")
})

test_that("the two arms record DIFFERENT reads of the same tissue", {
  # arm 1's xlsx is the pathologist's ORIGINAL read, arm 2's the re-read. If these
  # ever became equal, the two tables would be one table and the whole per-arm split
  # would be ceremony — so the difference is asserted, not assumed.
  both <- dplyr::inner_join(m1, m2, by = c("SAMPLE", "annotation"),
                            suffix = c("_1", "_2"))
  differing <- dplyr::filter(both, !is.na(path_pct_1), !is.na(path_pct_2),
                             path_pct_1 != path_pct_2)
  expect_gt(nrow(differing), 0)
  # 046 ANNOTATION_1 is the clearest case: 50 in arm 1, 30 in arm 2.
  expect_equal(dplyr::filter(m1, SAMPLE == "046", annotation == "ANNOTATION_1")$path_pct, 50)
  expect_equal(dplyr::filter(m2, SAMPLE == "046", annotation == "ANNOTATION_1")$path_pct, 30)
})

test_that("the arm 1 stub covers exactly the regions arm 1 exports", {
  # 046/5456/24086 three selected regions each, 052 two, and 10338/15897 one apiece
  # from the union polygon promoted to ANNOTATION_1 by arm_promote_unregioned().
  # A row here with no export is a score for a region nobody drew; an export with no
  # row silently drops that region from the correlation.
  want <- c(`046` = 3, `052` = 2, `5456` = 3, `24086` = 3, `10338` = 1, `15897` = 1)
  got  <- table(m1$SAMPLE)
  expect_equal(as.integer(got[names(want)]), unname(want))
  expect_equal(nrow(m1), 13)
})

test_that("arm 1 has no whole_slide row", {
  # 10338 and 15897 DO have a polygon — their single annotation_all boundary — so
  # they are ANNOTATION_1, not `whole_slide`. That label means "no polygon exists",
  # which is arm 2's 24086 and nothing in arm 1.
  expect_false("whole_slide" %in% m1$annotation)
})

# --- The two arms are NOT interchangeable ------------------------------------
test_that("the arms disagree about 24086, and both tables are entitled to", {
  # THE case that proves regions are arm-local: three annotated regions in arm 1,
  # none at all in arm 2. A single shared table could not represent both, and using
  # one arm's percentages for the other would score polygons that do not overlap.
  expect_equal(nrow(dplyr::filter(m1, SAMPLE == "24086")), 3)
  expect_equal(dplyr::filter(m2, SAMPLE == "24086")$annotation, "whole_slide")
})

test_that("neoplastic_for maps the inverted arm to arm 1's scores", {
  # massimo1_inverted re-classifies arm 1's own regions rather than redrawing them,
  # so it shares arm 1's table — a pathologist's percentage is a property of the
  # tissue, not of the classifier. Checked on the source text, since load_data.R
  # cannot be executed here.
  src <- paste(readLines(here::here("code", "load_data.R")), collapse = "\n")
  expect_match(src, "massimo1_inverted\\s*=\\s*neoplastic_massimo1")
})
