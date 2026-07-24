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
