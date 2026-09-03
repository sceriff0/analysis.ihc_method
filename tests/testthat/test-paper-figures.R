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

test_that("the tumour and compartment readings colour from the named compartment palette", {
  ov <- paper_phenotype_map(.cells(), patient_id = "046", colour_by = "tumour")
  expect_equal(levels(ov$data$call), c("Tumour", "Non-tumour"))
  # the finer lineage stays on the frame, so a caller can always recover it
  expect_true(all(levels(ov$data$lineage) %in% c(LEGIBLE_LINEAGES, "other")))
  cols <- unique(ggplot2::ggplot_build(ov)$data[[1]]$colour)
  expect_setequal(cols, unname(COMPARTMENT_COLS[c("Tumour", "Non-tumour")]))

  ins <- paper_phenotype_map(.cells(), patient_id = "046", colour_by = "compartment")
  expect_equal(levels(ins$data$call), COMPARTMENTS)
  cols <- unique(ggplot2::ggplot_build(ins)$data[[1]]$colour)
  expect_true(all(cols %in% unname(COMPARTMENT_COLS)))
})

test_that("lineage_compartment() pools the immune lineages and never colours an unmapped call", {
  x <- c("Tumor", "CD8T", "CD4T", "Treg", "NK", "Immune_other", "Stroma", "other", NA)
  expect_equal(as.character(lineage_compartment(x)),
               c("Tumour", "Immune", "Immune", "Immune", "Immune", "Immune",
                 "Stroma", "Other", "Other"))
  expect_equal(as.character(lineage_compartment(x, binary = TRUE)),
               c("Tumour", rep("Non-tumour", 8)))
  # the overview's grey and the inset's grey are the SAME grey
  expect_equal(unname(COMPARTMENT_COLS["Non-tumour"]), unname(COMPARTMENT_COLS["Other"]))
})

test_that("a region-tier table that lists a cell per region file draws each cell once", {
  one   <- .cells()
  three <- dplyr::bind_rows(one, one, one)
  expect_equal(nrow(paper_phenotype_map(three, patient_id = "046")$data),
               nrow(paper_phenotype_map(one,   patient_id = "046")$data))
})

test_that("a zoom becomes the coordinate window, so outlines past it are clipped", {
  p <- paper_phenotype_map(.cells(), patient_id = "046", zoom = c(0, 800, 0, 600))
  expect_equal(p$coordinates$limits$x, c(0, 800))
  expect_equal(p$coordinates$limits$y, c(0, 600))
  expect_null(paper_phenotype_map(.cells(), patient_id = "046")$coordinates$limits$x)
})

# The fixture cells sit on the diagonal of a 4000 x 3000 field, so ANNOTATION_2 is
# placed where that diagonal crosses it — an inset window must hold some cells.
.polys <- function() {
  sq <- function(x0, y0, w) sf::st_polygon(list(rbind(c(x0, y0), c(x0 + w, y0),
                                                       c(x0 + w, y0 + w), c(x0, y0 + w),
                                                       c(x0, y0))))
  sf::st_sf(patient_id = c("046", "046", "052"),
            annotation = c("ANNOTATION_1", "ANNOTATION_2", "ANNOTATION_1"),
            geometry = sf::st_sfc(sq(0, 0, 1000), sq(2000, 1500, 500), sq(0, 0, 100)))
}

test_that("the inset region is the polygon with the most tumour cells, boxed and padded", {
  skip_if_not_installed("sf")
  per <- tibble::tibble(patient_id = c("046", "046", "052"),
                        annotation = c("ANNOTATION_1", "ANNOTATION_2", "ANNOTATION_1"),
                        n_tumor_inside = c(10L, 500L, 3L))
  r <- tumour_richest_region(per, .polys(), "046", pad = 0)
  expect_equal(r$annotation, "ANNOTATION_2")
  expect_equal(r$n_tumor, 500L)
  expect_equal(r$zoom, c(2000, 2500, 1500, 2000))
  # padding widens the box on every side by the given fraction of it
  expect_equal(tumour_richest_region(per, .polys(), "046", pad = 0.1)$zoom,
               c(1950, 2550, 1450, 2050))
  # a per row with no polygon behind it (csv fallback) cannot be chosen
  per2 <- tibble::tibble(patient_id = "046", annotation = "csv", n_tumor_inside = 9L)
  expect_null(tumour_richest_region(per2, .polys(), "046"))
  expect_null(tumour_richest_region(per, .polys(), "99999"))
})

test_that("the pair is two panels, the overview tumour-only and the inset by compartment", {
  skip_if_not_installed("sf")
  per <- tibble::tibble(patient_id = c("046", "046"),
                        annotation = c("ANNOTATION_1", "ANNOTATION_2"),
                        n_tumor_inside = c(10L, 500L))
  pr <- paper_phenotype_map_pair(.cells(), "046", annots = .polys(), per = per)
  expect_s3_class(pr$overview, "ggplot")
  expect_s3_class(pr$inset, "ggplot")
  expect_equal(pr$region$annotation, "ANNOTATION_2")
  expect_equal(levels(pr$overview$data$call), c("Tumour", "Non-tumour"))
  expect_equal(levels(pr$inset$data$call), COMPARTMENTS)
  expect_equal(pr$inset$labels$title, "Region 2")
  expect_match(pr$overview$labels$subtitle, "cells$")
  # without a per table there is no inset, and the overview still comes back
  alone <- paper_phenotype_map_pair(.cells(), "046", annots = .polys())
  expect_s3_class(alone$overview, "ggplot")
  expect_null(alone$inset)
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
  # One panel, populations distinguished by colour from the named lineage palette,
  # never by facet — the four share axes so their abundances can be compared.
  expect_s3_class(p$facet, "FacetNull")
  built <- ggplot2::ggplot_build(p)
  cols  <- unique(built$data[[1]]$colour)
  expect_true(all(cols %in% unname(LINEAGE_COLS[c("CD8T", "CD4T", "Treg", "NK")])))
  expect_match(p$labels$y, "quanTIseq")
})

test_that("the deconvolution panel keeps only populations both sides resolve", {
  paired <- tibble::tibble(method = "quantiseq",
                           lineage = c("CD8T", "CD8T", "B cell", NA, "NK"),
                           patient_id = c("1", "2", "1", "2", "1"),
                           ihc_frac = c(.1, .2, .3, .4, NA), score = c(.1, .2, .3, .4, .5))
  p <- paper_deconv_scatter(paired, "quantiseq")
  expect_equal(nrow(p$data), 2L)
  expect_equal(as.character(unique(p$data$lineage)), "CD8T")
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
