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
  # This fixture carries no preset/micro column, so the arm figures warn that they
  # cannot group runs. That is tested on its own below; here it is incidental.
  figs <- suppressWarnings(build_reg_figs(d))
  expect_true("05_valis_vs_overlap_agreement" %in% names(figs))
})

# --- the paper's arm comparison (Fig 4b/4c) ----------------------------------

test_that("arms are built from whichever preset/micro columns param_matrix carries", {
  d <- tmp_data()
  readr::write_csv(tibble::tibble(
    run_id               = c("r1", "r2", "r3", "r4"),
    reg_max_image_dim    = c(850, 850, 2000, 2000),
    reg_micro_reg        = c(TRUE, FALSE, TRUE, FALSE),
    valis_non_rigid_rTRE = c(0.0021, 0.0034, 0.0018, 0.0029),
    reg_dice_matched     = c(0.91, 0.86, 0.93, 0.88)
  ), file.path(d, "param_matrix.csv"))
  figs <- build_reg_figs(d)

  expect_true(all(c("06_tre_by_arm", "07_dice_by_arm") %in% names(figs)))
  arms <- unique(as.character(figs[["06_tre_by_arm"]]$data$arm))
  # Two presets x micro on/off = four configurations, and each is one arm.
  expect_equal(sum(figs[["06_tre_by_arm"]]$data$kind == "configuration"), 4)
  expect_true(any(grepl("max_image_dim", arms)))
  expect_true(any(grepl("micro", arms)))
})

test_that("the baselines come from the stage columns, not from extra runs", {
  d <- tmp_data()
  readr::write_csv(tibble::tibble(
    run_id               = c("r1", "r2"),
    reg_micro_reg        = c(TRUE, FALSE),
    valis_non_rigid_rTRE = c(0.0021, 0.0034),
    reg_dice_matched     = c(0.91, 0.86)
  ), file.path(d, "param_matrix.csv"))
  # The same runs' before-and-after columns ARE the baselines: no extra sweep.
  readr::write_csv(tibble::tibble(
    img_name        = c("s1", "s2"),
    original_rTRE   = c(0.021, 0.018),
    rigid_rTRE      = c(0.0051, 0.0043),
    non_rigid_rTRE  = c(0.0021, 0.0019)
  ), file.path(d, "registration_valis_rtre.csv"))
  readr::write_csv(tibble::tibble(
    moving       = rep(c("s1", "s2"), each = 3),
    stage        = rep(c("native", "rigid", "non_rigid"), 2),
    dice_matched = c(0.41, 0.79, 0.91, 0.38, 0.81, 0.92)
  ), file.path(d, "registration_accuracy.csv"))
  figs <- build_reg_figs(d)

  for (nm in c("06_tre_by_arm", "07_dice_by_arm")) {
    arms <- unique(as.character(figs[[nm]]$data$arm[figs[[nm]]$data$kind == "baseline"]))
    expect_setequal(arms, c("no registration", "rigid only"))
  }
})

test_that("a param_matrix with no arm columns warns rather than pretending one arm", {
  d <- tmp_data()
  readr::write_csv(tibble::tibble(
    run_id               = c("r1", "r2"),
    valis_non_rigid_rTRE = c(0.0021, 0.0034),
    reg_dice_matched     = c(0.91, 0.86)
  ), file.path(d, "param_matrix.csv"))
  expect_warning(figs <- build_reg_figs(d), "no preset/micro column")
  # Still the right SHAPE of figure — one point per run, labelled by run.
  expect_equal(sum(figs[["06_tre_by_arm"]]$data$kind == "configuration"), 2)
})

test_that("the arm figures are skipped, not broken, without param_matrix.csv", {
  d <- tmp_data()
  readr::write_csv(tibble::tibble(
    img_name = "s1", original_rTRE = 0.02, rigid_rTRE = 0.005, non_rigid_rTRE = 0.002
  ), file.path(d, "registration_valis_rtre.csv"))
  figs <- build_reg_figs(d)
  expect_false(any(c("06_tre_by_arm", "07_dice_by_arm") %in% names(figs)))
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
