# The mirage cell source: phenotypes.csv + morphology.csv joined on `label`, one
# directory per patient. These tests build a synthetic results tree rather than
# needing a pipeline run — no off-repo data.
source(here::here("code", "cell_tables.R"))
source(here::here("code", "validation_helpers.R"))
source(here::here("code", "mirage_cells.R"))

mirage_tree <- function(patients = list(`046` = 6), quant = TRUE, morphology = TRUE,
                        phenotypes = TRUE) {
  root <- file.path(tempdir(), paste0("mirage-", as.integer(Sys.time()), "-", sample(1e6, 1)))
  for (pid in names(patients)) {
    n <- patients[[pid]]
    d <- file.path(root, pid)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    if (phenotypes) {
      dir.create(file.path(d, "phenotyping"), recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(tibble::tibble(
        label       = seq_len(n),
        phenotype   = rep_len(c("PANCK_Tumor", "T_cytotoxic", "Unclassified"), n),
        `sign:CD3`  = rep_len(c("1", "0", "·"), n),
        `state:PD1` = rep_len(c(1L, -1L, 0L), n)),
        file.path(d, "phenotyping", "phenotypes.csv"))
    }
    if (morphology) {
      dir.create(file.path(d, "cell_properties"), recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(tibble::tibble(
        label = seq_len(n), y = seq_len(n) * 10, x = seq_len(n) * 20,
        area = 100, eccentricity = .5, perimeter = 40, convex_area = 110,
        axis_major_length = 12, axis_minor_length = 8, solidity = .9),
        file.path(d, "cell_properties", "morphology.csv"))
    }
    if (quant) {
      dir.create(file.path(d, "quantification"), recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(tibble::tibble(
        label = seq_len(n), `CD3: Cytoplasm: Median` = seq_len(n) * 1.5),
        file.path(d, "quantification", "merged_quant.csv"))
    }
  }
  root
}

test_that("phenotypes and morphology are joined on label, per patient directory", {
  cells <- load_mirage_cells(mirage_tree(list(`046` = 6, `052` = 3)))
  expect_equal(nrow(cells), 9)
  expect_setequal(unique(cells$patient_id), c("046", "052"))
  expect_true(all(c("phenotype", "phenotype_clean", "x_px", "y_px", "area") %in% names(cells)))
})

test_that("the patient id comes from the directory and is slide_key normalised", {
  # nothing inside mirage's files identifies the patient
  cells <- load_mirage_cells(mirage_tree(list(`EPM - 052` = 4)))
  expect_equal(unique(cells$patient_id), "052")
})

test_that("morphology's x/y become x_px/y_px and are NOT rescaled", {
  cells <- load_mirage_cells(mirage_tree(list(`046` = 3)))
  # mirage centroids are already pixels; applying the micron conversion as well
  # would put every cell outside its annotation
  expect_equal(cell_centroids_px(cells, um_per_px = 0.325)$x, c(20, 40, 60))
  expect_equal(cell_centroids_px(cells, um_per_px = 0.325)$y, c(10, 20, 30))
  expect_false(any(c("x", "y") %in% names(cells)))
})

test_that("a patient with no morphology is skipped, not silently placed", {
  # without centroids there is no way to decide membership, and mirage has no
  # out-of-annotation flag to fall back on
  expect_warning(cells <- load_mirage_cells(
    mirage_tree(list(`046` = 3), morphology = FALSE)), "no cell_properties/morphology.csv")
  expect_equal(nrow(cells), 0)
})

test_that("merged_quant is optional and supplies intensities when present", {
  with_q <- load_mirage_cells(mirage_tree(list(`046` = 3)))
  expect_equal(marker_value(with_q, "CD3"), c(1.5, 3.0, 4.5))
  no_q <- load_mirage_cells(mirage_tree(list(`046` = 3), quant = FALSE))
  expect_equal(nrow(no_q), 3)
  expect_true(all(is.na(marker_value(no_q, "CD3"))))
  expect_equal(marker_pos(no_q, "CD3"), c(TRUE, FALSE, FALSE))   # signs still work
})

test_that("a missing mirage directory warns and yields no cells", {
  expect_warning(cells <- load_mirage_cells(file.path(tempdir(), "no-such-mirage-dir")),
                 "directory not found")
  expect_equal(nrow(cells), 0)
})

test_that("the inventory reports the unresolved share", {
  inv <- mirage_cells_inventory(load_mirage_cells(mirage_tree(list(`046` = 6))))
  expect_equal(inv$n_cells, 6)
  expect_equal(inv$n_resolved, 4)          # 2 of the 6 are "Unclassified"
  expect_equal(inv$n_phenotypes, 2)        # PANCK_Tumor, T_cytotoxic
})

test_that("mirage cells flow through the shared metric helpers", {
  cells <- load_mirage_cells(mirage_tree(list(`046` = 6)))
  r <- region_ratios(cells)
  expect_equal(r$n_inside, 6)
  expect_equal(r$n_inside_clean, 4)        # the unresolved cells leave the clean set
  expect_equal(r$n_tumor_inside, 2)        # PANCK_Tumor matched despite the _ spelling
  expect_false(any(is.na(cell_lineage(cells$phenotype_clean))))
})

test_that("cohort-level directories in a raw outdir are skipped silently", {
  # A mirage outdir carries phenotyping/ (COMPILE_PANEL), qc/, size_logs/ and
  # _UNROUTED_PUBLISH/ alongside the patient dirs. Warning about each of them on
  # every knit would train the reader to ignore the warnings that matter, so a
  # patient is recognised by structure — it has phenotyping/phenotypes.csv.
  root <- mirage_tree(list(`046` = 3))
  for (d in c("phenotyping", "qc", "size_logs", "_UNROUTED_PUBLISH"))
    dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)
  writeLines("model", file.path(root, "phenotyping", "model_config.json"))

  expect_no_warning(cells <- load_mirage_cells(root))
  expect_equal(unique(cells$patient_id), "046")
  expect_false(is_mirage_patient_dir(file.path(root, "phenotyping")))
})


# ---------------------------------------------------------------------------
# Runs built WITHOUT the phenotyping stage
# ---------------------------------------------------------------------------
# mirage ships phenotyping as a separate, optional feature. A pipeline built
# without it emits quantification/ and cell_properties/ and no phenotyping/.
# is_mirage_patient_dir() used to gate on phenotypes.csv alone, so every such
# patient was classified "not a patient" and skipped without comment — the whole
# cohort then surfaced as a single "no patient directories" warning that reads
# like an empty dataset rather than a differently-built pipeline.

test_that("a patient with no phenotyping stage is still a patient", {
  root <- mirage_tree(list(`046` = 6), phenotypes = FALSE)
  expect_true(is_mirage_patient_dir(file.path(root, "046")))
})

test_that("cells load without phenotypes, keeping coordinates and intensities", {
  root <- mirage_tree(list(`046` = 6, `052` = 3), phenotypes = FALSE)
  cells <- suppressWarnings(load_mirage_cells(root))
  expect_equal(nrow(cells), 9)
  expect_setequal(unique(cells$patient_id), c("046", "052"))
  # coordinates survive under the names the membership metrics use, unrescaled
  expect_true(all(c("x_px", "y_px", "area") %in% names(cells)))
  expect_equal(sort(cells$x_px[cells$patient_id == "052"]), c(20, 40, 60))
  # quantification still joins
  expect_true("CD3: Cytoplasm: Median" %in% names(cells))
  # and every cell is unresolved rather than silently typed, so the inventory's
  # denominators stay honest instead of reporting a 100%-resolved cohort
  expect_true(all(is_unresolved_phenotype(cells$phenotype_clean)))
  inv <- mirage_cells_inventory(cells)
  expect_equal(sum(inv$n_resolved), 0)
  expect_equal(sum(inv$n_cells), 9)
})

test_that("a cohort with no phenotyping says so once, and says what is still exact", {
  root <- mirage_tree(list(`046` = 6, `052` = 3), phenotypes = FALSE)
  w <- testthat::capture_warnings(load_mirage_cells(root))
  hit <- grep("no phenotyping stage", w, value = TRUE)
  expect_length(hit, 1)                       # once for the cohort, not per patient
  expect_match(hit, "densities are exact")
})

test_that("a partially phenotyped cohort names the count", {
  root <- mirage_tree(list(`046` = 6))        # phenotyped
  d <- file.path(root, "052")                 # not
  dir.create(file.path(d, "cell_properties"), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(tibble::tibble(
    label = 1:3, y = c(10, 20, 30), x = c(20, 40, 60),
    area = 100, eccentricity = .5, perimeter = 40, convex_area = 110,
    axis_major_length = 12, axis_minor_length = 8, solidity = .9),
    file.path(d, "cell_properties", "morphology.csv"))
  w <- testthat::capture_warnings(load_mirage_cells(root))
  expect_match(paste(w, collapse = " "), "1 of 2 patient\\(s\\) have no")
})

test_that("a patient with only quantification is a patient, but is skipped loudly", {
  # It has no coordinates, so its cells cannot be placed in an annotation. That is
  # a broken patient rather than a non-patient, and the difference must be audible.
  root <- mirage_tree(list(`046` = 6), phenotypes = FALSE, morphology = FALSE)
  expect_true(is_mirage_patient_dir(file.path(root, "046")))
  expect_warning(cells <- load_mirage_cells(root), "no cell_properties/morphology.csv")
  expect_equal(nrow(cells), 0)
})
