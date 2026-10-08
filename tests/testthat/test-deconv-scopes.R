# The cell scopes behind the quanTIseq page. A scope that silently reads the wrong
# cell set still knits and still reports a rho, so each test pins one rule about
# WHICH cells a scope keeps. No data on disk: cells are injected, and the arm spec
# points at an empty directory so every patient is "one with no polygon".
source(here::here("code", "cell_tables.R"))
source(here::here("code", "validation_helpers.R"))
source(here::here("code", "arm_cells.R"))
source(here::here("code", "deconv_scopes.R"))

# Two patients, eight cells each; the LAST four of each are flagged outside. Inside,
# patient A is half CD8 T cells and patient B has none.
.scope_cells <- function() tibble::tibble(
  patient_id        = rep(c("046", "052"), each = 8),
  cell_id           = rep(1:8, 2),
  centroid_x        = rep(1:8, 2), centroid_y = rep(1:8, 2),
  phenotype         = c("T cytotoxic", "T cytotoxic", "PANCK+Tumor", "Stroma",
                        rep("T helper", 4),
                        "PANCK+Tumor", "PANCK+Tumor", "Stroma", "Natural Killer",
                        rep("T cytotoxic", 4)),
  phenotype_clean   = phenotype,
  Out_of_annotation = rep(c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE), 2))

.empty_spec <- function(arm) {
  d <- file.path(tempdir(), paste0("scopes-", sample(1e6, 1)))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  arm_spec(arm, data_dir = d)
}

test_that("every scope names a real arm and an existing polygon tier", {
  expect_true(all(DECONV_SCOPES$arm %in% ARM_MODES))
  expect_true(all(DECONV_SCOPES$tier %in% c(NA, "union", "region")))
  expect_equal(anyDuplicated(DECONV_SCOPES$scope), 0)
  # A "union" scope on an arm with no union tier would have no polygon to read and
  # fall through to the flag for EVERY patient, looking like a result.
  for (i in which(DECONV_SCOPES$tier %in% "union"))
    expect_true(arm_has_union_tier(ARM_SPECS[[DECONV_SCOPES$arm[i]]]))
})

test_that("an unknown scope is an error, not an empty panel", {
  expect_error(deconv_scope("massimo3_annotation"), "unknown deconvolution scope")
})

test_that("a whole-slide scope keeps every cell and consults nothing", {
  sc <- deconv_scope_cells("massimo2_wholeslide", .scope_cells())
  expect_true(all(sc$in_scope))
  expect_equal(unique(sc$.scope_source), "all_cells")
})

test_that("an annotation scope cuts a patient with no polygon by its export's flag", {
  sc <- suppressWarnings(deconv_scope_cells("massimo2_annotation", .scope_cells(),
                                            spec = .empty_spec("massimo2")))
  expect_equal(unique(sc$.scope_source), "flag")
  expect_equal(sc$in_scope, !sc$Out_of_annotation)
})

test_that("the IHC fraction is taken over the cells the scope keeps, and only those", {
  cells <- .scope_cells()
  whole <- deconv_scope_ihc("massimo2_wholeslide", cells)
  ann   <- suppressWarnings(deconv_scope_ihc("massimo2_annotation", cells,
                                             spec = .empty_spec("massimo2")))
  cd8 <- function(r, pid) r$fractions$frac_all[r$fractions$lineage == "CD8T" &
                                                r$fractions$patient_id == pid]
  expect_equal(cd8(whole, "046"), 2 / 8)
  expect_equal(cd8(ann,   "046"), 2 / 4)
  expect_equal(cd8(whole, "052"), 4 / 8)
  expect_equal(cd8(ann,   "052"), 0)        # completed to 0, not a missing row
  expect_equal(unique(ann$fractions$scope), "massimo2_annotation")
  expect_equal(ann$inventory$n_in_scope, c(4, 4))
  expect_equal(ann$inventory$pct_in_scope, c(50, 50))
  expect_equal(ann$inventory$membership, c("flag", "flag"))
})

test_that("concordance is one Spearman per scope and population", {
  paired <- tibble::tibble(
    scope = rep(c("massimo1_wholeslide", "massimo2_wholeslide"), each = 4),
    method = "quantiseq", lineage = "CD8T", patient_id = rep(letters[1:4], 2),
    score = c(1, 2, 3, 4, 1, 2, 3, 4), ihc_frac = c(1, 2, 3, 4, 4, 3, 2, 1))
  out <- deconv_scope_concordance(paired)
  expect_equal(out$rho, c(1, -1))
  expect_equal(out$n, c(4L, 4L))
})

test_that("the scatter names its scope and states its n in patients", {
  d <- tibble::tibble(patient_id = letters[1:4], score = c(.1, .2, .3, .4),
                      ihc_frac = c(.1, .3, .2, .4), immuno_phe = c("hot", "cold", NA, "hot"))
  p <- plot_deconv_scope_pair(d, "massimo2_annotation", "CD8T")
  expect_match(p$labels$subtitle, "Massimo2 — annotation")
  expect_match(p$labels$subtitle, "n = 4 patients")
  expect_match(p$labels$title, "quanTIseq: CD8T")
  expect_silent(ggplot2::ggplot_build(p))
})

# --- Single phenotype labels ---------------------------------------------------
test_that("every leaf faces a population whose IHC side actually contains it", {
  # A leaf set against a population it is not part of would correlate two unrelated
  # things and report a rho.
  for (i in seq_len(nrow(DECONV_LEAVES))) {
    row  <- deconv_comparison_lineages[deconv_comparison_lineages$lineage == DECONV_LEAVES$lineage[i], ]
    toks <- sub("\\[.*$", "", trimws(strsplit(row$ihc_phenotypes, ";")[[1]]))
    expect_true(pheno_join_key(DECONV_LEAVES$leaf[i]) %in% pheno_join_key(toks),
                info = paste(DECONV_LEAVES$leaf[i], "is not a leaf of", DECONV_LEAVES$lineage[i]))
  }
  expect_equal(anyDuplicated(pheno_join_key(DECONV_LEAVES$leaf)), 0)
})

test_that("a leaf fraction counts that label alone, over the cells the scope keeps", {
  ann <- suppressWarnings(deconv_scope_ihc("massimo2_annotation", .scope_cells(),
                                           spec = .empty_spec("massimo2")))
  f <- function(pid, leaf) ann$leaves$ihc_frac[ann$leaves$patient_id == pid & ann$leaves$leaf == leaf]
  expect_equal(f("046", "T cytotoxic"), 2 / 4)
  expect_equal(f("046", "Activated T cytotoxic"), 0)      # completed to 0
  expect_equal(f("052", "Natural Killer"), 1 / 4)
  expect_equal(unique(ann$leaves$lineage[ann$leaves$leaf == "CD8+ T reg"]), "CD8T")
  expect_equal(ann$unmapped, character(0))
})

test_that("a label the leaf table does not list is reported, not dropped in silence", {
  cells <- .scope_cells()
  cells$phenotype[1] <- cells$phenotype_clean[1] <- "Plasma cell"
  expect_equal(deconv_scope_ihc("massimo2_wholeslide", cells)$unmapped, "Plasma cell")
})

test_that("the leaf scatter names the label on x and the population on y", {
  d <- tibble::tibble(patient_id = letters[1:4], score = c(.1, .2, .3, .4),
                      ihc_frac = c(.1, .3, .2, .4))
  p <- plot_deconv_scope_pair(d, "massimo1_wholeslide", "CD8T", leaf = "CD8+ T reg")
  expect_match(p$labels$title, "CD8T vs phenotype 'CD8\\+ T reg'")
  expect_match(p$labels$x, "'CD8\\+ T reg' fraction")
  expect_match(p$labels$y, "quanTIseq CD8T score")
  expect_false(any(vapply(p$layers, function(l) inherits(l$geom, "GeomAbline"), logical(1))))
})
