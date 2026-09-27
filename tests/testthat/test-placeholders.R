# The opt-in placeholder mode (code/placeholders.R). What these tests pin is the
# INTEGRITY GUARD, not the arithmetic: off means untouched, on means every synthetic
# row is flagged and no real row is edited, the draw is reproducible, a figure that
# draws a synthetic point says so, and nothing synthetic reaches output/figures/.
source(here::here("code", "placeholders.R"))
source(here::here("code", "pdf_export.R"))

# Two arms x two patients x two stages; arm "b" never ran for patient "p2", and one
# real row has an NA metric.
ph_real <- function() {
  tibble::tibble(
    arm     = c("a", "a", "a", "a", "b", "b"),
    backend = "valis",
    patient = c("p1", "p1", "p2", "p2", "p1", "p1"),
    stage   = factor(c("rigid", "micro", "rigid", "micro", "rigid", "micro"),
                     levels = c("rigid", "micro")),
    disp    = c(4, 1.0, 5, 1.2, 3, NA),
    dice    = c(.60, .80, .55, .78, .62, .81))
}
ph_expected <- function() {
  tidyr::expand_grid(arm = c("a", "b", "c"), patient = c("p1", "p2"),
                     stage = c("rigid", "micro")) |>
    dplyr::mutate(backend = "valis")
}
ph_call <- function(real = ph_real(), ...) {
  placeholder_fill(real, ph_expected(), keys = c("arm", "patient", "stage"),
                   metrics = c("disp", "dice"),
                   levels = list(arm = c("arm", "stage"), backend = c("backend", "stage")),
                   ranges = list(dice = c(0, 1), disp = c(1e-6, Inf)),
                   log_scale = "disp", what = "test", ...)
}

test_that("the mode is OFF by default and then adds nothing — not even a column", {
  withr::local_options(ihc.placeholder_missing = NULL)
  withr::local_envvar(IHC_PLACEHOLDER_MISSING = NA)
  expect_false(placeholder_mode())
  real <- ph_real()
  out  <- suppressMessages(ph_call(real))
  expect_identical(out, real)
  expect_equal(nrow(out), nrow(real))
})

test_that("the env var and the option both switch it on; an explicit option wins", {
  withr::local_options(ihc.placeholder_missing = NULL)
  withr::local_envvar(IHC_PLACEHOLDER_MISSING = "1")
  expect_true(placeholder_mode())
  withr::local_options(ihc.placeholder_missing = FALSE)
  expect_false(placeholder_mode())
  withr::local_options(ihc.placeholder_missing = TRUE)
  withr::local_envvar(IHC_PLACEHOLDER_MISSING = NA)
  expect_true(placeholder_mode())
})

test_that("placeholder mode flags every synthetic row and never alters a real one", {
  withr::local_options(ihc.placeholder_missing = TRUE)
  real <- ph_real()
  out  <- suppressMessages(ph_call(real))
  expect_true(all(c("is_placeholder", "placeholder_rule") %in% names(out)))
  kept <- dplyr::select(out[!out$is_placeholder, ], -is_placeholder, -placeholder_rule)
  expect_identical(kept, real)                       # real rows byte-identical
  syn <- out[out$is_placeholder, ]
  # 12 expected keys - 6 real = 6 missing rows, plus ONE companion for the NA disp.
  expect_equal(nrow(syn), 7)
  expect_true(all(nzchar(syn$placeholder_rule)))
  # The real row with the NA keeps its NA; the value is on a companion row.
  expect_true(is.na(out$disp[!out$is_placeholder & out$arm == "b" & out$stage == "micro"]))
  comp <- syn[grepl("^NA metric", syn$placeholder_rule), ]
  expect_equal(nrow(comp), 1)
  expect_true(is.finite(comp$disp)); expect_true(is.na(comp$dice))
  # Ranges hold and the factor stays a factor with its own levels.
  expect_true(all(syn$dice[is.finite(syn$dice)] >= 0 & syn$dice[is.finite(syn$dice)] <= 1))
  expect_true(all(syn$disp[is.finite(syn$disp)] > 0))
  expect_s3_class(out$stage, "factor")
  expect_identical(levels(out$stage), c("rigid", "micro"))
})

test_that("the rule names the fallback level: same arm first, then the neighbours", {
  withr::local_options(ihc.placeholder_missing = TRUE)
  out <- suppressMessages(ph_call())
  syn <- out[out$is_placeholder & grepl("^no row", out$placeholder_rule), ]
  # Arm "b" has a real rigid disp, so its missing rigid rows use it; it has NO real
  # micro disp (the one row is NA), so its micro rows fall through to the neighbours.
  expect_true(all(grepl("disp arm mean", syn$placeholder_rule[syn$arm == "b" & syn$stage == "rigid"])))
  expect_true(all(grepl("disp backend mean", syn$placeholder_rule[syn$arm == "b" & syn$stage == "micro"])))
  # Arm "c" has no data at all: it falls through to the backend level.
  expect_true(all(grepl("disp backend mean", syn$placeholder_rule[syn$arm == "c"])))
  # No real value anywhere -> the fixed prior, and it says so.
  real <- dplyr::mutate(ph_real(), dice = NA_real_)
  expect_false(any(grepl("fixed prior", suppressMessages(ph_call(real))$placeholder_rule)))
  pr <- suppressMessages(placeholder_fill(
    ph_real(), ph_expected(), keys = c("arm", "patient", "stage"), metrics = "disp",
    levels = list(), what = "prior-test", na_structure = NULL,
    enabled = TRUE, priors = list(disp = list(mean = 1, sd = .1))))
  expect_true(any(grepl("global", pr$placeholder_rule)))
})

test_that("a metric absent BY CONSTRUCTION is never invented", {
  withr::local_options(ihc.placeholder_missing = TRUE)
  # `delta` exists only at the micro stage (like d_disp_um_vs_rigid at rigid).
  real <- dplyr::mutate(ph_real(), delta = ifelse(stage == "micro", -2, NA_real_))
  out  <- suppressMessages(placeholder_fill(
    real, ph_expected(), keys = c("arm", "patient", "stage"), metrics = c("disp", "delta"),
    levels = list(arm = c("arm", "stage")), na_structure = "stage", what = "structure"))
  syn <- out[out$is_placeholder, ]
  expect_true(all(is.na(syn$delta[syn$stage == "rigid"])))
  expect_true(all(is.finite(syn$delta[syn$stage == "micro" & grepl("^no row", syn$placeholder_rule)])))
})

test_that("the draw is deterministic under the seed and leaves the caller's RNG alone", {
  withr::local_options(ihc.placeholder_missing = TRUE, ihc.placeholder_seed = 11L)
  set.seed(99); before <- .Random.seed
  a <- suppressMessages(ph_call())
  expect_identical(.Random.seed, before)
  b <- suppressMessages(ph_call())
  expect_identical(a, b)
  withr::local_options(ihc.placeholder_seed = 12L)
  c <- suppressMessages(ph_call())
  expect_false(identical(a$disp, c$disp))
  expect_identical(a[!a$is_placeholder, ], c[!c$is_placeholder, ])
})

test_that("a second pass is idempotent", {
  withr::local_options(ihc.placeholder_missing = TRUE)
  a <- suppressMessages(ph_call())
  b <- suppressMessages(ph_call(a))
  expect_identical(a, b)
})

test_that("the sidecar registry lists every synthetic value", {
  withr::local_options(ihc.placeholder_missing = TRUE)
  placeholder_reset()
  out <- suppressMessages(ph_call())
  reg <- placeholder_registry()
  n_vals <- sum(is.finite(out$disp[out$is_placeholder])) +
            sum(is.finite(out$dice[out$is_placeholder]))
  expect_equal(nrow(reg), n_vals)
  expect_true(all(reg$table == "test"))
  root <- withr::local_tempdir()
  path <- suppressMessages(placeholder_write_sidecar("slugx", root = root))
  expect_true(file.exists(path))
  expect_equal(nrow(utils::read.csv(path)), n_vals)
  placeholder_reset()
})

# --- the figure helper ----------------------------------------------------------
ph_plot <- function(d) {
  ggplot2::ggplot(d, ggplot2::aes(arm, disp)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(width = .1, height = 0) +
    ggplot2::geom_line(ggplot2::aes(group = patient)) +
    ggplot2::labs(subtitle = "the real subtitle", caption = "cap")
}

test_that("the watermark is added only when the plot draws a placeholder", {
  real_p <- ph_plot(ph_real())
  expect_identical(placeholder_style(real_p), real_p)          # no flag column
  withr::local_options(ihc.placeholder_missing = TRUE)
  flagged_real <- suppressMessages(ph_call(dplyr::filter(ph_real(), !is.na(disp)),
                                           enabled = TRUE))
  none <- ph_plot(flagged_real[!flagged_real$is_placeholder, ])
  expect_identical(placeholder_style(none), none)               # flag column, zero TRUE

  d <- suppressMessages(ph_call())
  d <- d[is.finite(d$disp), ]
  p <- placeholder_style(ph_plot(d))
  n <- sum(d$is_placeholder)
  expect_match(p$labels$subtitle, sprintf("^PLACEHOLDER — %d synthetic points", n))
  expect_match(p$labels$subtitle, "the real subtitle")
  labs <- unlist(lapply(p$layers, function(l) l$aes_params$label))
  expect_true("PLACEHOLDER" %in% labs)
  # point and line layers are split in two; the boxplot is not.
  expect_equal(length(p$layers), 3 + 2 + 1)
  b <- ggplot2::ggplot_build(p)
  pts <- vapply(b$data, nrow, integer(1))
  expect_equal(pts[[2]], sum(!d$is_placeholder))   # real points, filled
  expect_equal(pts[[3]], sum(d$is_placeholder))    # synthetic points, hollow
  expect_equal(b$data[[3]]$shape[1], 1)
  # Styling twice does not stack.
  expect_identical(placeholder_style(p)$labels$subtitle, p$labels$subtitle)
})

# --- the export guard -----------------------------------------------------------
fake_knit_pdf <- function() {
  fig <- withr::local_tempdir(.local_envir = parent.frame())
  grDevices::pdf(file.path(fig, "fig-1.pdf")); plot(1); grDevices::dev.off()
  knitr::opts_chunk$set(fig.path = paste0(fig, "/"))
  withr::defer(knitr::opts_chunk$set(fig.path = "figure/"), envir = parent.frame())
  fig
}

test_that("export_pdf_figures writes output/figures/ normally", {
  withr::local_options(ihc.placeholder_missing = FALSE)
  fake_knit_pdf()
  out <- withr::local_tempdir(); ph <- withr::local_tempdir()
  suppressMessages(export_pdf_figures("pagex", out_root = out, placeholder_root = ph))
  expect_true(file.exists(file.path(out, "pagex", "fig-1.pdf")))
  expect_length(list.files(ph, recursive = TRUE), 0)
})

test_that("in placeholder mode export_pdf_figures NEVER writes output/figures/", {
  withr::local_options(ihc.placeholder_missing = TRUE)
  fake_knit_pdf()
  out <- withr::local_tempdir(); ph <- withr::local_tempdir()
  suppressMessages(export_pdf_figures("pagex", out_root = out, placeholder_root = ph))
  expect_length(list.files(out, recursive = TRUE), 0)
  expect_true(file.exists(file.path(ph, "figures", "pagex", "fig-1_PLACEHOLDER.pdf")))
  expect_true(file.exists(file.path(ph, "pagex_placeholders.csv")))
})

test_that("manuscript exporters refuse, and derived tables move aside", {
  withr::local_options(ihc.placeholder_missing = FALSE)
  expect_invisible(placeholder_refuse("figures/fig4.pdf"))
  expect_identical(placeholder_output_path("/x/output/t.csv"), "/x/output/t.csv")
  withr::local_options(ihc.placeholder_missing = TRUE)
  expect_error(placeholder_refuse("figures/fig4.pdf"), "refusing")
  root <- withr::local_tempdir()
  expect_identical(placeholder_output_path("/x/output/t.csv", root = root),
                   file.path(root, "t_PLACEHOLDER.csv"))
})

test_that("the page callout is empty when off and loud when on", {
  withr::local_options(ihc.placeholder_missing = FALSE)
  expect_identical(placeholder_callout(), "")
  withr::local_options(ihc.placeholder_missing = TRUE)
  expect_match(placeholder_callout(), "PLACEHOLDER MODE IS ON")
})
