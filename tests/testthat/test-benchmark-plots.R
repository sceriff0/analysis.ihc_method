# benchmark_plots.R must source end-to-end from a minimal sweep, and must NOT build
# the two registration-accuracy figures: they read the same param_matrix.csv columns
# as registration_accuracy_plots.R (§4 accuracy-vs-cost, §5 VALIS-vs-overlap) from the
# same data/benchmark/ directory, and are owned by that page. Re-adding them here —
# the easy mistake when merging upstream mirage plots.R — puts the identical figure on
# two pages under two numbers, which this test exists to catch.

# measurements.csv is read unconditionally at source time and every figure block runs
# top-to-bottom over it, so the fixture has to be richer than the figure under test:
#   - powerlaw_plot() (figs 01/02) needs >= 2 distinct input_gb per process, or it
#     returns NULL and save_fig(NULL, ...) errors on `NULL + labs(...)`
#   - fig 07 is unconditional and eagerly group_bys n_channels / n_register_images /
#     target_px — a missing column errors outright rather than skipping the figure
minimal_sweep <- function() {
  d <- file.path(tempdir(), paste0("bench-", as.integer(Sys.time()), "-", sample(1e6, 1)))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(tibble::tibble(
    run_id = c("r1", "r2"), process = "MIRAGE:REGISTER",
    peak_rss_gb = c(1, 2), realtime_s = c(1, 2),
    input_gb = c(1, 2), varied_axis = "baseline",
    n_channels = 2, n_register_images = 2, target_px = 512), file.path(d, "measurements.csv"))
  d
}

build_bench_figs <- function(dir) {
  bench_figs <- list()                    # the sourced script fills this via save_fig
  commandArgs <- function(...) dir        # shadow so `adir` resolves to the fixture
  source(here::here("code", "benchmark_plots.R"), local = TRUE)
  bench_figs
}

test_that("a minimal sweep builds the resource figures", {
  figs <- build_bench_figs(minimal_sweep())
  expect_true(any(grepl("^0[12]_", names(figs))))
  expect_true("07_stage_memory_heatmap" %in% names(figs))
})

test_that("registration accuracy is left to registration_accuracy_plots.R", {
  d <- minimal_sweep()
  # A param_matrix carrying BOTH accuracy signals — everything the retired figs 11 and
  # 17 needed. They must still not appear.
  readr::write_csv(tibble::tibble(
    run_id                  = c("r1", "r2", "r3"),
    cpu_hours               = c(1.2, 2.5, 4.1),
    reg_displacement_um_p50 = c(2.1, 1.4, 1.3),
    reg_dice_matched        = c(0.72, 0.80, 0.83),
    valis_non_rigid_D       = c(5.1, 3.8, 3.2)), file.path(d, "param_matrix.csv"))
  figs <- build_bench_figs(d)
  expect_false(any(grepl("accuracy_vs_cost|valis", names(figs))))
})
