# The manuscript panels. These tests are about the CLAIMS the figure legends make
# — that a panel reports no coefficient, that the population list is fixed, that an
# unmapped phenotype cannot acquire a colour — because those are the things a
# careless edit would quietly break while the plot still renders.
source(here::here("code", "cell_tables.R"))
source(here::here("code", "validation_helpers.R"))
source(here::here("code", "paper_figures.R"))

.cells <- function(n = 400, phen = NULL) {
  set.seed(3)
  cells <- tibble::tibble(
    patient_id = "046",
    cell_id    = seq_len(n),
    x_px       = seq(0, 3999, length.out = n),
    y_px       = seq(0, 2999, length.out = n),
    phenotype  = rep_len(phen %||% c("PANCK+Tumor", "T cytotoxic", "T helper",
                                     "Natural Killer", "Immune", "Stroma"), n))
  cells$phenotype_clean <- cell_phenotype(cells)
  cells
}

test_that("the phenotype map builds and keeps only the legible lineage subset", {
  p <- paper_phenotype_map(.cells(), patient_id = "046")
  expect_s3_class(p, "ggplot")
  expect_true(all(levels(p$data$lineage) %in% c(LEGIBLE_LINEAGES, "other")))
})

test_that("a phenotype neither vocabulary knows becomes 'other', never a hue", {
  # The failure this guards against is silent: an unmapped label given its own
  # colour reads as a discovered population. It must land in the grey bucket.
  p <- paper_phenotype_map(.cells(phen = c("PANCK+Tumor", "SomeNewLabel")),
                           patient_id = "046")
  expect_true("other" %in% as.character(p$data$lineage))
  expect_false("SomeNewLabel" %in% as.character(p$data$lineage))
})

test_that("the y axis is reversed so the map is not the slide upside down", {
  p <- paper_phenotype_map(.cells(), patient_id = "046")
  expect_true(any(vapply(p$scales$scales,
                         function(s) inherits(s, "ScaleContinuousPosition") &&
                                     isTRUE(s$trans$name == "reverse"), logical(1))))
})

test_that("the zoom window cuts an inset from the same call, not a different one", {
  full <- paper_phenotype_map(.cells(), patient_id = "046")
  ins  <- paper_phenotype_map(.cells(), patient_id = "046", zoom = c(0, 800, 0, 600))
  expect_lt(nrow(ins$data), nrow(full$data))
  expect_true(all(ins$data$x <= 800))
  expect_true(all(ins$data$y <= 600))
})

test_that("an empty zoom window warns and returns NULL rather than an empty panel", {
  expect_warning(p <- paper_phenotype_map(.cells(), patient_id = "046",
                                          zoom = c(1e6, 1e6 + 1, 1e6, 1e6 + 1)),
                 "no cells")
  expect_null(p)
})

test_that("the scale bar picks a round length that fits the field", {
  expect_equal(.nice_bar_um(6000), 1000)   # 6000/6 = 1000 exactly
  expect_equal(.nice_bar_um(1300), 100)    # target ~217 -> largest nice <= it
  expect_true(.nice_bar_um(30) <= 5)       # a tiny inset still gets a bar
})

test_that("the hot/cold panel reports no correlation coefficient", {
  # The legend's claim, enforced on the object: no layer may carry a statistic, and
  # no label may contain rho/r/p. Every other page reports Spearman freely; this
  # panel is the one that must not.
  m <- tibble::tibble(patient_id = c("046", "052", "10338", "15897", "5456", "24086"),
                      cd45_over_inside = c(.31, .28, .35, .09, .11, .07))
  g <- tibble::tibble(patient_id = m$patient_id,
                      group = c("hot", "hot", "hot", "cold", "cold", "cold"))
  p <- paper_immune_fraction_hotcold(m, g)
  expect_s3_class(p, "ggplot")
  labs <- unlist(p$labels)
  expect_false(any(grepl("rho|\\br\\b|p *=|p-value|CI", labs, ignore.case = TRUE)))
  expect_false(any(vapply(p$layers, function(l) inherits(l$stat, "StatSmooth"), logical(1))))
  # Every case is drawn — with three a side the points are the evidence.
  expect_equal(nrow(p$data), 6)
})

test_that("the hot/cold panel orders cold before hot", {
  m <- tibble::tibble(patient_id = c("a", "b"), cd45_over_inside = c(.1, .3))
  g <- tibble::tibble(patient_id = c("a", "b"), group = c("cold", "hot"))
  p <- paper_immune_fraction_hotcold(m, g)
  expect_equal(levels(p$data$group)[1], "cold")
})

test_that("the deconvolution panel keeps one method and draws no fit line", {
  paired <- tidyr::expand_grid(method = c("quantiseq", "epic"),
                               lineage = c("CD8T", "CD4T", "Treg", "NK"),
                               patient_id = as.character(1:6))
  paired$ihc_frac <- runif(nrow(paired)); paired$score <- runif(nrow(paired))
  p <- paper_deconv_scatter(paired, "quantiseq")
  expect_equal(unique(p$data$method), "quantiseq")
  expect_false(any(vapply(p$layers, function(l) inherits(l$stat, "StatSmooth"), logical(1))))
})

test_that("asking for a method that was not run warns and names what is available", {
  paired <- tibble::tibble(method = "epic", lineage = "CD8T",
                           ihc_frac = .1, score = .2, patient_id = "1")
  expect_warning(p <- paper_deconv_scatter(paired, "quantiseq"), "epic")
  expect_null(p)
})

test_that("the lineage table maps exactly four deconvolution types and no more", {
  # Additional file 4's content, and the answer to the legend's "full list of
  # populations plotted". deconv_to_lineage() resolves four; everything else is the
  # catch-all the legend calls Tier 3. If this number changes, Fig 5(c) changes.
  t <- paper_lineage_table()
  d <- dplyr::filter(t, side == "deconvolution (cell type)")
  expect_setequal(setdiff(unique(d$lineage), "(unmapped)"),
                  c("CD4T", "CD8T", "Treg", "NK"))
  expect_gt(sum(d$lineage == "(unmapped)"), 0)
  # The imaging side is generated from the join table itself, so it cannot drift.
  expect_equal(sum(t$side == "imaging (phenotype call)"), nrow(phenotype_lineage_labels))
})
