# The compatibility contract: the same four facts (phenotype, membership,
# positivity, position) must come back identically from all three upstream export
# schemas — FlowPath's PhenotypeCsvExporter, mirage's phenotype_cells.py, and
# mirage's join_flowpath cohort table. See code/cell_tables.R for the schemas.
source(here::here("code", "cell_tables.R"))

flowpath <- tibble::tibble(
  cell_id           = 0:1,
  phenotype         = c("CD8 T cell (T_cytotoxic)", "Stromal (Stroma)"),
  Out_of_annotation = c("True", "False"),          # FlowPath writes the strings
  Outlier           = c("False", "False"),
  centroid_x        = c(325, 650),                 # microns
  centroid_y        = c(0, 325),
  CD3_sign          = c("+", ""),                  # "" = never gated
  CD45_sign         = c("+", "-"))

mirage_phenotypes <- tibble::tibble(
  label       = 1:2,
  phenotype   = c("T_cytotoxic", "Stroma"),        # bare label, no parentheses
  `sign:CD3`  = c("1", "·"),                  # · = never gated
  `sign:CD45` = c("1", "0"),
  `state:PD1` = c(1L, -1L))                        # 1 pos / -1 neg / 0 free

mirage_cohort <- tibble::tibble(
  label                = 1:2,
  phenotype            = c("T_cytotoxic", "Stroma"),
  fp_out_of_annotation = c(TRUE, FALSE),
  fp_outlier           = c(FALSE, FALSE),
  x_px                 = c(1000, 2000),            # already pixels
  y_px                 = c(0, 1000),
  CD3_positive         = c(TRUE, FALSE),
  CD45_positive        = c(TRUE, FALSE))

test_that("each export is recognised as its own dialect", {
  expect_equal(cell_dialect(flowpath), "flowpath")
  expect_equal(cell_dialect(mirage_phenotypes), "mirage_phenotypes")
  expect_equal(cell_dialect(mirage_cohort), "mirage_cohort")
})

test_that("clean_phenotype takes the parenthetical only when there is one", {
  # The regression this guards: a bare label under the old `(?<=\\().*?(?=\\))`
  # extraction became NA, so every mirage-phenotyped cell silently dropped out of
  # the tumour and lineage counts instead of failing loudly.
  expect_equal(clean_phenotype("CD8 T cell (T_cytotoxic)"), "T_cytotoxic")
  expect_equal(clean_phenotype("PANCK_Tumor"), "PANCK_Tumor")
  expect_equal(clean_phenotype(c("  ", NA)), c(NA_character_, NA_character_))
  expect_equal(cell_phenotype(mirage_phenotypes), c("T_cytotoxic", "Stroma"))
})

test_that("membership resolves across the flag spellings", {
  expect_equal(cell_outside(flowpath),          c(TRUE, FALSE))
  expect_equal(cell_outside(mirage_cohort),     c(TRUE, FALSE))
  # phenotype_cells.py carries no flag: no cell is excluded rather than all of them
  expect_equal(cell_outside(mirage_phenotypes), c(FALSE, FALSE))
  expect_true(has_outside_flag(flowpath))
  expect_false(has_outside_flag(mirage_phenotypes))
})

test_that("positivity resolves across _sign / sign: / _positive", {
  for (cells in list(flowpath, mirage_phenotypes, mirage_cohort)) {
    expect_equal(marker_pos(cells, "CD3"),  c(TRUE, FALSE))
    expect_equal(marker_pos(cells, "CD45"), c(TRUE, FALSE))
  }
  # three-valued state columns: only +1 is positive
  expect_equal(marker_pos(mirage_phenotypes, "PD1"), c(TRUE, FALSE))
  # a marker the export never gated yields no cells, it does not error
  expect_equal(marker_pos(flowpath, "FOXP3"), c(FALSE, FALSE))
  expect_true(is.na(marker_sign_col(flowpath, "FOXP3")))
})

test_that("centroids land in the same pixel frame from either unit", {
  expect_equal(cell_centroids_px(flowpath, um_per_px = 0.325),
               data.frame(x = c(1000, 2000), y = c(0, 1000)))
  # mirage's cohort table is already in pixels — dividing again would put every
  # cell outside its annotation
  expect_equal(cell_centroids_px(mirage_cohort, um_per_px = 0.325),
               data.frame(x = c(1000, 2000), y = c(0, 1000)))
  expect_error(cell_centroids_px(mirage_phenotypes), "no centroids")
})

test_that("the union dedup key pairs an id with the centroid", {
  # cell_id alone is unsafe: FlowPath re-indexes from 0 per export, so two cells in
  # different per-annotation files can share a row number.
  expect_equal(cell_key_cols(flowpath), c("cell_id", "centroid_x", "centroid_y"))
  expect_equal(cell_key_cols(mirage_cohort), c("label", "x_px", "y_px"))
  expect_equal(cell_key_cols(mirage_phenotypes), "label")
})

test_that("flag columns are coerced to logical so two exports can be bound", {
  # readr reads "True"/"False" as logical and a hand-written csv as character;
  # binding the two then fails on a type clash rather than on the numbers.
  a <- normalise_cell_flags(flowpath)
  b <- normalise_cell_flags(dplyr::mutate(flowpath, Out_of_annotation = c(TRUE, FALSE),
                                                    Outlier = c(FALSE, FALSE)))
  expect_type(a$Out_of_annotation, "logical")
  expect_type(a$Outlier, "logical")
  expect_no_error(dplyr::bind_rows(a, b))
})

test_that("read_cell_csv attaches phenotype_clean and skips a non-cell table", {
  f <- tempfile(fileext = ".csv")
  readr::write_csv(flowpath, f)
  cells <- read_cell_csv(f, patient_id = "046")
  expect_equal(cells$phenotype_clean, c("T_cytotoxic", "Stroma"))
  expect_equal(unique(cells$patient_id), "046")
  expect_type(cells$Out_of_annotation, "logical")

  g <- tempfile(fileext = ".csv")
  readr::write_csv(tibble::tibble(a = 1:2), g)
  expect_warning(out <- read_cell_csv(g), "not a cell table")
  expect_equal(nrow(out), 0)
})

test_that("marker_gated separates 'never measured' from 'negative'", {
  # The distinction all three exports make and a plain sign read loses: FlowPath
  # leaves an ungated column blank, mirage writes "·" (free) or "x" (contradictory).
  ungated <- tibble::tibble(CD3_sign = c("+", "-", ""), `sign:CD8` = c("1", "0", "·"))
  expect_equal(marker_gated(ungated, "CD3"), c(TRUE, TRUE, FALSE))
  expect_equal(marker_gated(ungated, "CD8"), c(TRUE, TRUE, FALSE))
  expect_equal(marker_pos(ungated, "CD3"),   c(TRUE, FALSE, FALSE))
  # a marker the panel never carried is ungated everywhere, not negative everywhere
  expect_equal(marker_gated(ungated, "FOXP3"), c(FALSE, FALSE, FALSE))
  expect_equal(marker_sign(flowpath, "CD3"), c("+", ""))
  expect_true(all(is.na(marker_sign(flowpath, "FOXP3"))))
})

test_that("marker_matrix names its columns for the markers, not the export", {
  for (cells in list(flowpath, mirage_phenotypes, mirage_cohort)) {
    m <- marker_matrix(cells, c("CD3", "CD45", "FOXP3"))
    expect_named(m, c("CD3", "CD45", "FOXP3"))
    expect_equal(m$CD3, c(TRUE, FALSE))
    expect_equal(m$FOXP3, c(FALSE, FALSE))   # absent marker: no cells, no error
  }
})

test_that("marker value and z-score resolve per dialect", {
  fp <- dplyr::mutate(flowpath, CD3_raw = c(10, 2), CD3_zscore = c(1.5, -1))
  mc <- dplyr::mutate(mirage_cohort, CD3 = c(10, 2), CD3_zscore = c(1.5, -1))
  expect_equal(marker_value(fp, "CD3"),  c(10, 2))   # <M>_raw
  expect_equal(marker_value(mc, "CD3"),  c(10, 2))   # bare <M>
  expect_equal(marker_zscore(fp, "CD3"), c(1.5, -1))
  expect_equal(marker_zscore(mc, "CD3"), c(1.5, -1))
  # phenotype_cells.py carries no intensities at all
  expect_true(all(is.na(marker_value(mirage_phenotypes, "CD3"))))
  expect_true(all(is.na(marker_zscore(mirage_phenotypes, "CD3"))))
})

# ── the mirage / FlowPath phenotype vocabulary ───────────────────────────────
# The two tools name the SAME taxonomy differently. An unmapped label does not
# error, it joins to lineage NA and quietly empties the composition panels, so
# this is the table that has to stay in sync with mirage's panel.yaml.

test_that("both vocabularies collapse to the same lineages", {
  flowpath_labels <- c("PANCK+Tumor", "VIM+Tumor", "T helper", "T cytotoxic",
                       "Activated T cytotoxic", "CD8+ T reg", "CD4+ Treg",
                       "Natural Killer", "Activated Natural Killer", "Immune", "Stroma")
  mirage_labels   <- c("PANCK_Tumor", "VIM_Tumor", "T_helper", "T_cytotoxic",
                       "Activated_T_cytotoxic", "CD8_Treg", "CD4_Treg",
                       "NK_cell", "Activated_NK", "Immune", "Stroma")
  expect_equal(cell_lineage(flowpath_labels), cell_lineage(mirage_labels))
  expect_false(any(is.na(cell_lineage(flowpath_labels))))
})

test_that("every phenotype mirage's panel can emit is mapped", {
  # mirage/panel.yaml `phenotypes:` plus palette.py RESERVED. If mirage adds a leaf,
  # this fails and the lineage table needs the new name.
  panel <- c("Immune", "Activated_T_cytotoxic", "CD8_Treg", "T_cytotoxic", "CD4_Treg",
             "T_helper", "Activated_NK", "NK_cell", "Myeloid", "Macrophage_M2",
             "Stroma", "PANCK_Tumor", "VIM_Tumor")
  reserved <- c("Ambiguous", "Conflict", "Artefact", "Unclassified")
  expect_false(any(is.na(cell_lineage(panel))))
  expect_equal(unname(cell_lineage(reserved)), rep("Unknown", 4))
  # mirage-only leaves fold onto the branch FlowPath dead-ends at
  expect_equal(unname(cell_lineage(c("Myeloid", "Macrophage_M2"))),
               c("Immune_other", "Immune_other"))
})

test_that("unresolved covers both tools' spellings", {
  expect_true(all(is_unresolved_phenotype(
    c("Unknown", "Unclassified", "Ambiguous", "Conflict", "Artefact", NA, ""))))
  expect_false(any(is_unresolved_phenotype(c("T_helper", "PANCK+Tumor", "Stroma"))))
})

test_that("mirage's measurement-named intensity columns resolve", {
  q <- tibble::tibble(label = 1:2, `CD3: Cytoplasm: Median` = c(10, 2),
                      `PANCK: Cytoplasm: Median` = c(1, 9))
  expect_equal(marker_measurement_col(q, "CD3"), "CD3: Cytoplasm: Median")
  expect_equal(marker_value(q, "CD3"), c(10, 2))
  expect_true(all(is.na(marker_value(q, "FOXP3"))))
  # must not match a marker that merely shares a prefix
  expect_true(is.na(marker_measurement_col(q, "CD")))
})
