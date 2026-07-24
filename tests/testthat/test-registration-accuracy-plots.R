# Smoke tests for the registration-accuracy figure builder. Each writes a tiny
# schema-shaped CSV to a tempdir and checks the right ggplot is produced — no
# off-repo sweep data required.
source(here::here("code", "registration_accuracy_plots.R"))

tmp_data <- function() {
  d <- file.path(tempdir(), paste0("regfig-", as.integer(Sys.time()), "-", sample(1e6, 1)))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

test_that("empty dir yields an empty figure list without error", {
  expect_length(build_reg_figs(tmp_data()), 0)
})

test_that("registration_valis_rtre.csv builds the rTRE slopegraph and n_matches panel", {
  d <- tmp_data()
  # columns as mirage make_tables.py emits them: run_id + summary_csv + VALIS's own columns
  readr::write_csv(tibble::tibble(
    run_id         = "run0000",
    summary_csv    = "P001_summary.csv",
    img_name       = c("P001_mov1", "P001_mov2"),
    original_rTRE  = c(50.2, 40.0),
    rigid_rTRE     = c(10.1, 9.0),
    non_rigid_rTRE = c(5.1, 4.0),
    n_matches      = c(100, 80)
  ), file.path(d, "registration_valis_rtre.csv"))
  figs <- build_reg_figs(d)
  expect_true("01_valis_rtre_by_stage" %in% names(figs))
  expect_true("01b_valis_n_matches"   %in% names(figs))
  expect_s3_class(figs[["01_valis_rtre_by_stage"]], "ggplot")
})

test_that("registration_accuracy.csv builds the overlap Dice and displacement figures", {
  d <- tmp_data()
  readr::write_csv(tibble::tibble(
    run_id             = "r1",
    moving             = rep(c("P001_mov1", "P001_mov2"), each = 4),
    stage              = rep(c("native", "rigid", "non_rigid", "micro"), 2),
    dice_matched       = c(.10, .55, .72, .74, .12, .50, .70, .71),
    displacement_um_p50 = c(9.0, 3.1, 1.3, 1.2, 8.5, 3.4, 1.5, 1.4),
    displacement_um_p90 = c(18.0, 6.2, 2.6, 2.4, 17.0, 6.8, 3.0, 2.8)
  ), file.path(d, "registration_accuracy.csv"))
  figs <- build_reg_figs(d)
  expect_true(all(c("02_overlap_dice_by_stage", "02b_displacement_um_by_stage") %in% names(figs)))
  expect_s3_class(figs[["02b_displacement_um_by_stage"]], "ggplot")
})

test_that("param_matrix.csv builds the accuracy-vs-cost Pareto figure", {
  d <- tmp_data()
  readr::write_csv(tibble::tibble(
    run_id                = c("r1", "r2", "r3"),
    cpu_hours             = c(1.2, 2.5, 4.1),
    reg_displacement_um_p50 = c(2.1, 1.4, 1.3)
  ), file.path(d, "param_matrix.csv"))
  figs <- build_reg_figs(d)
  expect_true("04_accuracy_vs_cost" %in% names(figs))
  expect_s3_class(figs[["04_accuracy_vs_cost"]], "ggplot")
})

test_that("param_matrix.csv builds the VALIS-vs-overlap agreement figure", {
  d <- tmp_data()
  readr::write_csv(tibble::tibble(
    run_id            = c("r1", "r2", "r3"),
    reg_dice_matched  = c(0.72, 0.80, 0.83),
    valis_non_rigid_D = c(5.1, 3.8, 3.2)
  ), file.path(d, "param_matrix.csv"))
  figs <- build_reg_figs(d)
  expect_true("05_valis_vs_overlap_agreement" %in% names(figs))
})

test_that("feature_dist/*.json builds the (legacy) distance-reduction figure when present", {
  skip_if_not_installed("jsonlite")
  d <- tmp_data()
  dir.create(file.path(d, "feature_dist"))
  jsonlite::write_json(list(
    moving_image = "P001_mov1",
    improvement  = list(distance_reduction_percent = 62.5)
  ), file.path(d, "feature_dist", "P001_mov1.json"), auto_unbox = TRUE)
  figs <- build_reg_figs(d)
  expect_true("03_feature_distance_reduction" %in% names(figs))
})
