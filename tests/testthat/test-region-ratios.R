# PD-L1 is counted three ways in region_ratios() — every PD-L1+ cell, and the PD-L1+
# cells of the tumour and CD45+ compartments — and the within-compartment
# positivities are what the clinical page's standalone PD-L1 figure draws. These
# tests pin that the numerators are restricted to the right compartment, that the
# counts are read through marker_pos() (so the mirage spelling works too), and that
# an export with no PDL1 gate degrades to zero rather than erroring.
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
  expect_equal(r$pdl1_over_tumor,  2 / 4)
  expect_equal(r$pdl1_over_cd45,   1 / 4)
})

test_that("the mirage spelling of the PDL1 sign column is read too", {
  r <- region_ratios(.pdl1_cells(pdl1_col = "sign:PDL1", cd45_col = "sign:CD45"))
  expect_equal(r$n_pdl1_inside, 3)
  expect_equal(r$pdl1_over_tumor, 0.5)
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
  expect_true(is.na(r$pdl1_over_cd45))
  expect_equal(r$n_pdl1_cd45_inside, 0)
})

test_that("region_composition() accepts PD-L1+ as a marker population under the panel convention", {
  m <- region_ratios(.pdl1_cells()) |> dplyr::mutate(patient_id = "046", .before = 1)
  comp <- region_composition(m, markers = c("PD-L1+" = "n_pdl1_inside"))
  pd <- comp[comp$lineage == "PD-L1+", ]
  expect_equal(nrow(pd), 1)
  expect_equal(pd$frac_inside, 3 / 8)
  expect_equal(pd$frac_tumor,  3 / 4)   # ALL PD-L1+ cells over the tumour count
  expect_equal(pd$frac_cd45,   3 / 4)
})
