# Fig 11 must read the CURRENT mirage accuracy schema (param_matrix.csv,
# reg_displacement_um_p50) — not the retired quality.csv / reg_tre_median_px.
test_that("benchmark fig 11 reads param_matrix reg_displacement_um_p50", {
  d <- file.path(tempdir(), paste0("bench-", as.integer(Sys.time()), "-", sample(1e6, 1)))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  # measurements.csv is read unconditionally at source time and every figure block runs
  # top-to-bottom over it. Two problems with the brief's 1-row fixture surfaced when sourcing:
  #   1. powerlaw_plot() (figs 01/02) needs >=2 distinct input_gb per process to fit a power law
  #      and return a real ggplot; with 1 row it returns NULL and save_fig(NULL, ...) errors on
  #      `NULL + labs(...)`.
  #   2. fig 7 (STAGE-COST HEATMAP) is unconditional and eagerly filters/group_bys on
  #      n_channels, n_register_images, target_px — missing columns error immediately
  #      ("object not found"), not just an empty/skipped figure.
  # So the fixture carries two distinct input_gb rows (power-law fit succeeds) and the
  # n_channels/n_register_images/target_px columns fig 7 needs, letting the whole script source
  # end-to-end without touching what fig 11 actually asserts (param_matrix.csv-driven).
  readr::write_csv(tibble::tibble(
    run_id = c("r1", "r2"), process = "MIRAGE:REGISTER",
    peak_rss_gb = c(1, 2), realtime_s = c(1, 2),
    input_gb = c(1, 2), varied_axis = "baseline",
    n_channels = 2, n_register_images = 2, target_px = 512), file.path(d, "measurements.csv"))
  readr::write_csv(tibble::tibble(
    run_id = c("r1", "r2"), cpu_hours = c(1.2, 2.5),
    reg_displacement_um_p50 = c(2.1, 1.4)), file.path(d, "param_matrix.csv"))
  bench_figs <- list()                    # sourced script fills this via save_fig
  commandArgs <- function(...) d          # shadow so `adir` resolves to the fixture
  source(here::here("code", "benchmark_plots.R"), local = TRUE)
  expect_true("11_accuracy_vs_cost" %in% names(bench_figs))
})
