# Smoke tests for the ANHIR figure module. Each builds small schema-shaped
# cases/aggregates frames, writes them to a tempdir as the two hand-off CSVs,
# loads them with anhir_load() and checks each plot function returns a ggplot —
# no off-repo benchmark data required.
source(here::here("code", "anhir_plots.R"))

tmp_data <- function() {
  d <- file.path(tempdir(), paste0("anhir-", as.integer(Sys.time()), "-", sample(1e6, 1)))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

METHODS <- c("initial", "bunwarpj", "valis", "tiled")

# Two training cases (scored) and one evaluation case (unscored, metrics NA),
# each registered by all four methods — the shape pull_to_ihc_method.sh hands over.
synth_cases <- function() {
  scored_rtre <- c(0.080, 0.020, 0.004, 0.005,    # case 0
                   0.060, 0.015, 0.003, 0.0035)   # case 1
  tibble::tibble(
    case_id      = rep(c(0L, 1L, 2L), each = 4),
    tissue       = rep(c("COAD", "lung-lesion", "kidney"), each = 4),
    scale        = "scale-25pc",
    status       = rep(c("training", "training", "evaluation"), each = 4),
    source_image = "dataset/s.jpg",
    target_image = "dataset/t.jpg",
    method       = rep(METHODS, 3),
    n_landmarks  = 80L,
    scored       = rep(c(TRUE, TRUE, FALSE), each = 4),
    rtre_median  = c(scored_rtre, rep(NA_real_, 4)),
    rtre_mean    = c(scored_rtre * 1.2, rep(NA_real_, 4)),
    rtre_max     = c(scored_rtre * 3, rep(NA_real_, 4)),
    tre_median_px = c(scored_rtre * 4000, rep(NA_real_, 4)),
    robustness   = c(0, .70, .95, .90, 0, .60, .90, .92, rep(NA_real_, 4)),
    rank_median_rtre = c(4, 3, 1, 2, 4, 3, 1, 2, rep(NA_real_, 4)),
    time_min     = c(0, 3, 12, 8, 0, 3.5, 11, 9, rep(NA_real_, 4))
  )
}

synth_aggregates <- function() {
  one <- function(subset, n) tibble::tibble(
    method = METHODS, subset = subset, n_cases = n,
    avg_median_rtre = c(0.07, 0.0175, 0.0035, 0.00425),
    med_median_rtre = c(0.07, 0.0175, 0.0035, 0.00425),
    avg_mean_rtre   = c(0.084, 0.021, 0.0042, 0.0051),
    avg_max_rtre    = c(0.21, 0.0525, 0.0105, 0.01275),
    avg_robustness  = c(0, .65, .925, .91),
    med_robustness  = c(0, .65, .925, .91),
    avg_rank_median_rtre = c(4, 3, 1, 2),
    avg_time_min    = c(0, 3.25, 11.5, 8.5))
  dplyr::bind_rows(one("all", 2L), one("training", 2L),
                   dplyr::mutate(one("evaluation", 0L),
                                 dplyr::across(-c(method, subset, n_cases), ~ NA_real_)))
}

write_fixture <- function(d, cases = synth_cases(), aggregates = synth_aggregates()) {
  readr::write_csv(cases, file.path(d, "anhir_cases.csv"))
  readr::write_csv(aggregates, file.path(d, "anhir_aggregates.csv"))
  d
}

# ggplot2 4 hides user labels behind get_labs(); 3.x keeps them in $labels.
plot_title <- function(p) {
  if (exists("get_labs", asNamespace("ggplot2"))) ggplot2::get_labs(p)$title else p$labels$title
}

PLOT_FUNS <- list(
  plot_anhir_rtre_by_method = function(a) plot_anhir_rtre_by_method(a$cases),
  plot_anhir_rtre_by_tissue = function(a) plot_anhir_rtre_by_tissue(a$cases),
  plot_anhir_robustness     = function(a) plot_anhir_robustness(a$cases),
  plot_anhir_rank           = function(a) plot_anhir_rank(a$aggregates))

test_that("anhir_load reads both tables with typed columns", {
  a <- anhir_load(write_fixture(tmp_data()))
  expect_named(a, c("cases", "aggregates"))
  expect_equal(nrow(a$cases), 12)
  expect_true(is.logical(a$cases$scored))
  expect_true(is.numeric(a$cases$rtre_median))
  expect_true(is.numeric(a$aggregates$avg_rank_median_rtre))
  expect_equal(sum(a$cases$scored), 8)
})

test_that("anhir_load stops naming a missing column", {
  d <- write_fixture(tmp_data(), cases = dplyr::select(synth_cases(), -robustness))
  expect_error(anhir_load(d), "anhir_cases[.]csv.*robustness")
  d <- write_fixture(tmp_data(), aggregates = dplyr::select(synth_aggregates(), -avg_rank_median_rtre))
  expect_error(anhir_load(d), "anhir_aggregates[.]csv.*avg_rank_median_rtre")
})

test_that("anhir_load stops on an absent file, naming the hand-off script", {
  expect_error(anhir_load(tmp_data()), "pull_to_ihc_method")
})

test_that("an all-NA metric column still loads as numeric", {
  # Every case unscored: readr guesses the NA-only columns as logical.
  unscored <- dplyr::mutate(synth_cases(), scored = FALSE,
                            dplyr::across(c(rtre_median, rtre_mean, rtre_max, tre_median_px,
                                            robustness, rank_median_rtre, time_min), ~ NA_real_))
  a <- anhir_load(write_fixture(tmp_data(), cases = unscored))
  expect_true(is.numeric(a$cases$rtre_median))
  expect_true(is.numeric(a$cases$robustness))
})

test_that("every plot function returns a ggplot on a scored fixture", {
  a <- anhir_load(write_fixture(tmp_data()))
  for (nm in names(PLOT_FUNS)) {
    p <- PLOT_FUNS[[nm]](a)
    expect_true(inherits(p, "ggplot"), info = nm)
    expect_match(plot_title(p), "ANHIR", info = nm)
    # inherits() is satisfied by a plot that dies at draw time (a scale on a
    # missing aesthetic, a labeller on the wrong breaks); building it is not.
    expect_no_error(ggplot2::ggplot_build(p), message = nm)
  }
})

test_that("the case figures draw scored cases only", {
  a <- anhir_load(write_fixture(tmp_data()))
  p <- plot_anhir_rtre_by_method(a$cases)
  expect_equal(nrow(p$data), 8)
  expect_true(all(p$data$scored))
  expect_setequal(unique(p$data$status), "training")
})

test_that("robustness leaves out `initial`, which it is defined against", {
  a <- anhir_load(write_fixture(tmp_data()))
  p <- plot_anhir_robustness(a$cases)
  expect_false("initial" %in% as.character(p$data$method))
  expect_equal(nrow(p$data), 6)
})

test_that("the rank figure ranks subset == 'all' only, best method on top", {
  a <- anhir_load(write_fixture(tmp_data()))
  p <- plot_anhir_rank(a$aggregates)
  expect_equal(nrow(p$data), 4)
  expect_setequal(unique(p$data$subset), "all")
  # coord_flip() draws the last level at the top; VALIS has rank 1 in the fixture.
  expect_equal(tail(levels(p$data$method_lab), 1), "VALIS")
})

test_that("the summary table covers `all` and `training`, ranked, for kable", {
  a <- anhir_load(write_fixture(tmp_data()))
  tbl <- anhir_summary_table(a$aggregates)
  expect_s3_class(tbl, "data.frame")
  expect_setequal(unique(tbl$subset), c("all", "training"))
  expect_true(all(c("method", "n_cases", "mean_rank", "avg_median_rtre") %in% names(tbl)))
  expect_equal(tbl$method[tbl$subset == "all"][1], "VALIS")
  expect_equal(nrow(tbl), 8)
})

test_that("a cases table with zero scored rows draws an informative empty figure", {
  unscored <- dplyr::mutate(synth_cases(), scored = FALSE,
                            dplyr::across(c(rtre_median, rtre_mean, rtre_max, tre_median_px,
                                            robustness, rank_median_rtre, time_min), ~ NA_real_))
  no_all <- dplyr::filter(synth_aggregates(), subset != "all")
  a <- anhir_load(write_fixture(tmp_data(), cases = unscored, aggregates = no_all))
  for (nm in names(PLOT_FUNS)) {
    expect_no_error(p <- PLOT_FUNS[[nm]](a), message = nm)
    expect_true(inherits(p, "ggplot"), info = nm)
    expect_match(plot_title(p), "ANHIR", info = nm)
    expect_no_error(ggplot2::ggplot_build(p), message = nm)
  }
  # Zero rows, not an error, from the table too.
  expect_equal(nrow(anhir_summary_table(dplyr::filter(no_all, subset == "evaluation"))), 0)
})

test_that("a method the palette does not name still draws", {
  extra <- dplyr::mutate(synth_cases(), method = ifelse(method == "tiled", "ashlar", method))
  a <- anhir_load(write_fixture(tmp_data(), cases = extra))
  p <- plot_anhir_rtre_by_method(a$cases)
  expect_true(inherits(p, "ggplot"))
  expect_true("ashlar" %in% levels(p$data$method))
})
