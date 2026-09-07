# PD-L1 is counted three ways in region_ratios() — every PD-L1+ cell, and the PD-L1+
# cells of the tumour and CD45+ compartments. The clinical page's hot/cold figure
# draws PD-L1+ (and CD45+) cells PER TUMOUR CELL; the within-compartment
# positivities are kept as columns. These tests pin which numerator goes over which
# denominator, that the counts are read through marker_pos() (so the mirage
# spelling works too), and that an export with no PDL1 gate degrades to zero
# rather than erroring.
source(here::here("code", "cell_tables.R"))
source(here::here("code", "validation_helpers.R"))

# Eight cells: 4 tumour, 4 immune (all CD45+). PD-L1+ on 2 tumour and 1 immune cell.
.pdl1_cells <- function(pdl1_col = "PDL1_sign", cd45_col = "CD45_sign") {
  d <- tibble::tibble(
    patient_id      = "046",
    cell_id         = 1:8,
    phenotype_clean = c(rep("PANCK+Tumor", 4), rep("T cytotoxic", 4)))
  d[[cd45_col]] <- c(rep("-", 4), rep("+", 4))
  d[[pdl1_col]] <- c("+", "+", "-", "-", "+", "-", "-", "-")
  d
}

test_that("PD-L1 counts are split by compartment and the proportions use their own denominators", {
  r <- region_ratios(.pdl1_cells())
  expect_equal(r$n_pdl1_inside,       3)
  expect_equal(r$n_pdl1_tumor_inside, 2)
  expect_equal(r$n_pdl1_cd45_inside,  1)
  expect_equal(r$pdl1_over_inside, 3 / 8)
  expect_equal(r$pdl1_over_tumor,  3 / 4)   # ALL PD-L1+ cells per tumour cell
  expect_equal(r$cd45_over_tumor,  4 / 4)   # CD45+ cells per tumour cell
  expect_equal(r$pdl1_pos_in_tumor, 2 / 4)  # PD-L1+ fraction OF tumour cells
  expect_equal(r$pdl1_pos_in_cd45,  1 / 4)  # PD-L1+ fraction OF CD45+ cells
})

test_that("the mirage spelling of the PDL1 sign column is read too", {
  r <- region_ratios(.pdl1_cells(pdl1_col = "sign:PDL1", cd45_col = "sign:CD45"))
  expect_equal(r$n_pdl1_inside, 3)
  expect_equal(r$pdl1_over_tumor, 0.75)
})

test_that("an export that never gated PDL1 counts zero PD-L1+ cells, not an error", {
  d <- .pdl1_cells(); d$PDL1_sign <- NULL
  r <- region_ratios(d)
  expect_equal(r$n_pdl1_inside, 0)
  expect_equal(r$pdl1_over_tumor, 0)      # 0 / 4: the denominator still exists
})

test_that("a zero denominator makes the PD-L1 proportion NA", {
  d <- .pdl1_cells(); d$CD45_sign <- "-"
  r <- region_ratios(d)
  expect_true(is.na(r$pdl1_pos_in_cd45))
  expect_equal(r$n_pdl1_cd45_inside, 0)
})

test_that("a region with no tumour cells makes the per-tumour-cell ratios NA", {
  d <- .pdl1_cells(); d$phenotype_clean <- "T cytotoxic"
  r <- region_ratios(d)
  expect_true(is.na(r$cd45_over_tumor))
  expect_true(is.na(r$pdl1_over_tumor))
  expect_equal(r$n_pdl1_inside, 3)        # the count itself is unaffected
})

test_that("the hot/cold compartment panels carry their denominator in the strip", {
  source(here::here("code", "plot_theme.R"))
  lab <- flowpath_panel_label(c("tumor_over_inside", "cd45_over_tumor", "pdl1_over_tumor"))
  expect_match(lab[["tumor_over_inside"]], "tumour / all")
  expect_match(lab[["cd45_over_tumor"]],   "CD45\\+ / tumour")
  expect_match(lab[["pdl1_over_tumor"]],   "PD-L1\\+ / tumour")
})
