# The marker -> gene(s) -> slide bridge behind every IHC-vs-bulk-RNA figure.
source(here::here("code", "cell_tables.R"))
source(here::here("code", "validation_helpers.R"))

.counts <- function() tibble::tibble(
  Sample = c("A.001.004", "A.001.005", "A.009.999"),
  PTPRC  = c(10, 20, 30), CD3D = c(1, 2, 3), CD3E = c(3, 4, 5), CD8A = c(7, 8, 9))

test_that("a multi-gene marker is the mean of its genes, keyed to the slide id", {
  crf2slide <- c(A001004 = "366722", A001005 = "15370")
  r <- rna_marker_expression(.counts(), crf2slide)
  cd3 <- dplyr::filter(r$expr, marker == "CD3")
  expect_equal(cd3$rna_expr[cd3$patient_id == "366722"], mean(c(1, 3)))
  # The sample with no slide is dropped from the slide-keyed frame ...
  expect_false(any(is.na(r$expr$patient_id)))
  expect_equal(sort(unique(r$expr$patient_id)), c("15370", "366722"))
  # ... but kept in CRF space, which the immuno-phenotype boxplots read.
  expect_equal(dplyr::n_distinct(r$expr_crf$crf_id), 3)
})

test_that("genes absent from the matrix are reported and the marker uses the rest", {
  r <- rna_marker_expression(.counts(), c(A001004 = "046"))
  expect_true("CD8B" %in% r$missing_genes)
  expect_equal(r$expr$rna_expr[r$expr$marker == "CD8"], 7)
  expect_equal(r$genes[["CD8"]], "CD8A")
})
