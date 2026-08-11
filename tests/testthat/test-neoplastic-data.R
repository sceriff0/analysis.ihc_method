# The pathologist's neoplastic-cellularity reference. These tests pin the SHAPE and the
# join key, not the values — a value typo is caught by the reconciliation table the
# clinical_data page prints, but a shape or vocabulary change breaks the concordance
# join silently, dropping pairs while still producing a plausible correlation.
source(here::here("code", "cell_tables.R"))
source(here::here("code", "validation_helpers.R"))

# Read the tribble out of load_data.R without executing the file: load_data.R eagerly
# reads counts.RData and clinical_data.xlsx, which are off-repo.
neo <- local({
  src  <- readLines(here::here("code", "load_data.R"))
  from <- grep("^neoplastic_data <- tribble\\(", src)
  to   <- from - 1 + grep("^\\)", src[from:length(src)])[1]
  e <- new.env()
  e$tribble <- tibble::tribble
  eval(parse(text = paste(src[from:to], collapse = "\n")), envir = e)
  e$neoplastic_data
})

test_that("it is LONG: one row per scored region, not one per patient", {
  # Wide (~ANNOTATION_1..3) cannot express a whole-slide score without claiming a
  # region the pathologist never drew.
  expect_setequal(names(neo), c("SAMPLE", "annotation", "path_pct"))
  expect_gt(nrow(neo), dplyr::n_distinct(neo$SAMPLE))
})

test_that("its annotation labels are the vocabulary the metrics frames emit", {
  # _clinical_data_body.Rmd inner-joins on (patient_id, annotation). A label outside
  # this vocabulary silently drops that region from every correlation.
  ok <- c(paste0("ANNOTATION_", 1:3), "whole_slide")
  expect_true(all(neo$annotation %in% ok),
              info = paste("unexpected:", paste(setdiff(neo$annotation, ok), collapse = ", ")))
})

test_that("region numbering is contiguous from 1 per patient", {
  # ANNOTATION_k is the alphabet position of the export's letter suffix, so a patient
  # with an ANNOTATION_3 and no ANNOTATION_2 means a region went missing on one side.
  by_pt <- neo |>
    dplyr::filter(annotation != "whole_slide") |>
    dplyr::mutate(k = as.integer(sub("^ANNOTATION_", "", annotation))) |>
    dplyr::group_by(SAMPLE) |>
    dplyr::summarise(ks = list(sort(k)), .groups = "drop")
  for (i in seq_len(nrow(by_pt)))
    expect_equal(by_pt$ks[[i]], seq_along(by_pt$ks[[i]]),
                 info = paste("patient", by_pt$SAMPLE[i]))
})

test_that("a whole-slide patient has exactly one score and no regions", {
  ws <- dplyr::filter(neo, annotation == "whole_slide")
  for (pid in unique(ws$SAMPLE)) {
    rows <- dplyr::filter(neo, SAMPLE == pid)
    expect_equal(nrow(rows), 1, info = paste(pid, "should have one whole-slide score only"))
  }
})

test_that("percentages are percentages, not fractions", {
  # path_frac = path_pct / 100 downstream; a value already in 0..1 would come out 100x
  # too small and read as near-zero tumour content.
  expect_true(all(neo$path_pct >= 1 & neo$path_pct <= 100))
})

test_that("slide_key normalises every SAMPLE to the export's patient id", {
  # The join is on slide_key(SAMPLE); if that does not equal the csv directory name the
  # patient never matches.
  expect_equal(slide_key(neo$SAMPLE), neo$SAMPLE)
})
