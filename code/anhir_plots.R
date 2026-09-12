# =============================================================================
# anhir_plots.R — figures for analysis/benchmark_anhir.Rmd: how mirage's two
# registration backends score on the public ANHIR challenge, against the
# challenge's own baseline. Twin of registration_accuracy_plots.R in style, but a
# different contract: every plot function here is PURE — it takes a data frame
# and returns a ggplot, never touching disk — and `anhir_load(dir)` is the one
# reader. mirage's benchmarks/pull_to_ihc_method.sh drops the two CSVs into
# data/benchmark/; re-knit.
#
# Vocabulary (https://anhir.grand-challenge.org/Performance_Metrics/):
#   TRE          target registration error, px, per landmark after warping
#   rTRE         TRE / image diagonal — unitless, comparable across image sizes
#   robustness   fraction of a case's landmarks whose rTRE improved over `initial`
#   rank         the challenge's primary metric: each method's rank on per-case
#                median rTRE, averaged over cases (1 = best, lower = better)
# Only TRAINING cases carry target landmarks and are scored locally; evaluation
# cases are scored server-side, so their metric columns are NA (`scored = FALSE`).
# =============================================================================
.need <- c("ggplot2", "dplyr", "readr", "tidyr", "tibble")
.missing <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing))
  stop("Missing R packages: ", paste(.missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(lapply(.need, library, character.only = TRUE))

source(here::here("code", "plot_theme.R"))   # house theme + ANHIR_METHOD_COLS + label_n()

ANHIR_CAPTION <- "ANHIR challenge · training cases, scored locally on the released landmarks"

# The hand-off contract with mirage's benchmarks/anhir/. A renamed column is a
# stop() naming it, not a figure that silently draws the wrong thing.
ANHIR_CASE_COLS <- c("case_id", "tissue", "scale", "status", "source_image", "target_image",
                     "method", "n_landmarks", "scored", "rtre_median", "rtre_mean",
                     "rtre_max", "tre_median_px", "robustness", "rank_median_rtre",
                     "time_min")
ANHIR_AGG_COLS  <- c("method", "subset", "n_cases", "avg_median_rtre", "med_median_rtre",
                     "avg_mean_rtre", "avg_max_rtre", "avg_robustness", "med_robustness",
                     "avg_rank_median_rtre", "avg_time_min")
.ANHIR_CASE_NUM <- c("n_landmarks", "rtre_median", "rtre_mean", "rtre_max", "tre_median_px",
                     "robustness", "rank_median_rtre", "time_min")
.ANHIR_AGG_NUM  <- setdiff(ANHIR_AGG_COLS, c("method", "subset"))

# Reading order for the method axis: the two baselines first, then mirage's two
# backends. The palette (plot_theme.R) is keyed on what the CSV spells.
ANHIR_METHOD_LEVELS <- names(ANHIR_METHOD_COLS)

# --- loading ------------------------------------------------------------------
.anhir_read <- function(path, required, numeric) {
  if (!file.exists(path))
    stop("anhir_load: ", path, " does not exist — run mirage's ",
         "benchmarks/pull_to_ihc_method.sh first.", call. = FALSE)
  d <- readr::read_csv(path, show_col_types = FALSE)
  missing <- setdiff(required, names(d))
  if (length(missing))
    stop("anhir_load: ", basename(path), " is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  # An all-NA metric column (every case unscored) is guessed as logical by readr;
  # coerce so a downstream is.finite() sees a double either way.
  d %>% dplyr::mutate(dplyr::across(dplyr::all_of(numeric),
                                    ~ suppressWarnings(as.numeric(.x))))
}

# Read anhir_cases.csv + anhir_aggregates.csv from `dir`. Returns
# list(cases = <tibble>, aggregates = <tibble>); stops, naming the file and the
# columns, when either is absent or malformed.
anhir_load <- function(dir = here::here("data", "benchmark")) {
  cases <- .anhir_read(file.path(dir, "anhir_cases.csv"), ANHIR_CASE_COLS, .ANHIR_CASE_NUM) %>%
    # pandas writes True/False, R writes TRUE/FALSE; readr accepts both, but a
    # quoted or oddly-cased value would arrive as character.
    dplyr::mutate(scored = as.logical(scored))
  aggregates <- .anhir_read(file.path(dir, "anhir_aggregates.csv"), ANHIR_AGG_COLS, .ANHIR_AGG_NUM)
  list(cases = cases, aggregates = aggregates)
}

# --- shared helpers -----------------------------------------------------------
# Methods in reading order; an unlisted method (a future backend) is appended
# rather than dropped, so it still draws — in an unnamed palette colour.
.anhir_method_factor <- function(x) {
  x <- as.character(x)
  extra <- setdiff(unique(x[!is.na(x)]), ANHIR_METHOD_LEVELS)
  factor(x, levels = c(ANHIR_METHOD_LEVELS, sort(extra)))
}

# Add `method_lab`, the display-name factor the x axis is drawn on, keeping
# `method` (the key) for the named colour scale. Both are needed: the palette is
# keyed on the CSV spelling, the tick on what a reader should see.
.anhir_label_methods <- function(d) {
  d$method <- .anhir_method_factor(d$method)
  d$method_lab <- factor(.anhir_method_labels(d$method),
                         levels = .anhir_method_labels(levels(d$method)))
  d
}

# Scored cases with a finite value of `metric`. `positive` states the log-axis
# requirement (rTRE of exactly 0 does not occur, but filtering here says so rather
# than letting scale_y_log10() drop the row with a warning nobody reads).
.anhir_scored <- function(cases, metric, positive = FALSE) {
  d <- cases %>%
    dplyr::filter(scored %in% TRUE, is.finite(.data[[metric]]))
  if (positive) d <- dplyr::filter(d, .data[[metric]] > 0)
  .anhir_label_methods(d)
}

# The zero-scored-rows figure. Every plot function returns this instead of
# erroring, so a page knitted from the evaluation split alone — or before the
# harness has run — still renders and says why the panel is blank.
.anhir_empty <- function(title, what = "scored cases") {
  ggplot() +
    labs(title = title,
         subtitle = sprintf("No %s in the table — nothing to draw.", what),
         caption = ANHIR_CAPTION) +
    theme_void()
}

# --- figures ------------------------------------------------------------------

# Per-method distribution of the per-case MEDIAN rTRE, one facet per status.
# Only training cases are scored, so on a real hand-off this is one facet; the
# facet is kept so a harness that does score evaluation cases (a server-side
# leaderboard dump) lands in the same figure without a code change.
plot_anhir_rtre_by_method <- function(cases) {
  title <- "ANHIR: registration error by method"
  d <- .anhir_scored(cases, "rtre_median", positive = TRUE)
  if (!nrow(d)) return(.anhir_empty(title))
  ggplot(d, aes(method_lab, rtre_median, colour = method)) +
    geom_boxplot(outlier.shape = NA, width = .55, colour = "grey35") +
    geom_jitter(width = .12, height = 0, alpha = .7, size = 1.6) +
    # rTRE spans two orders of magnitude between `initial` and a good non-rigid
    # result; on a linear axis every registered method collapses onto zero.
    scale_y_log10() +
    scale_colour_anhir_method(guide = "none") +
    # label_n counts across facets; with one scored status that is the per-method n.
    scale_x_discrete(labels = label_n(d$method_lab)) +
    facet_wrap(~ status) +
    labs(title = title,
         subtitle = paste("One point per case: the median rTRE over that case's landmarks.",
                          "Lower = better.", n_note(d$case_id, "cases")),
         x = NULL, y = "median rTRE per case (fraction of image diagonal, log10)",
         caption = ANHIR_CAPTION)
}

# Method x tissue: the median over cases of the per-case median rTRE, tissues
# ordered hardest-first by their pooled median so the line profile reads as a
# difficulty ladder. n per tissue counts CASES, not (case, method) rows.
plot_anhir_rtre_by_tissue <- function(cases) {
  title <- "ANHIR: registration error by tissue"
  d <- .anhir_scored(cases, "rtre_median", positive = TRUE)
  if (!nrow(d)) return(.anhir_empty(title))
  d$tissue <- stats::reorder(factor(d$tissue), -d$rtre_median, FUN = stats::median)
  med <- d %>%
    dplyr::group_by(method, tissue) %>%
    dplyr::summarise(rtre = stats::median(rtre_median), .groups = "drop")
  per_tissue <- dplyr::distinct(d, tissue, case_id)
  ggplot(med, aes(tissue, rtre, colour = method, group = method)) +
    geom_line(alpha = .8) +
    geom_point(size = 2) +
    scale_y_log10() +
    scale_colour_anhir_method(name = NULL) +
    scale_x_discrete(labels = label_n(per_tissue$tissue)) +
    labs(title = title,
         subtitle = paste("Median over cases of the per-case median rTRE, per tissue;",
                          "hardest tissue on the left. Lower = better.",
                          n_note(d$case_id, "cases")),
         x = NULL, y = "median rTRE (fraction of image diagonal, log10)",
         caption = ANHIR_CAPTION) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# Robustness per method. `initial` is left out: robustness is DEFINED as the
# fraction of landmarks improved over `initial`, so its own value is identically
# zero and would only add a floor the axis already shows. The dashed line is the
# challenge's "robust" cut (> 0.5), the same one the `robust` aggregate subset uses.
plot_anhir_robustness <- function(cases) {
  title <- "ANHIR: robustness by method"
  d <- .anhir_scored(cases, "robustness") %>%
    dplyr::filter(as.character(method) != "initial") %>%
    .anhir_label_methods()
  if (!nrow(d)) return(.anhir_empty(title))
  ggplot(d, aes(method_lab, robustness, colour = method)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = REF_LINE) +
    geom_boxplot(outlier.shape = NA, width = .55, colour = "grey35") +
    geom_jitter(width = .12, height = 0, alpha = .7, size = 1.6) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_colour_anhir_method(guide = "none") +
    scale_x_discrete(labels = label_n(d$method_lab)) +
    labs(title = title,
         subtitle = paste("Fraction of a case's landmarks whose rTRE improved over no",
                          "registration; dashed = the challenge's robust cut (0.5). Higher = better.",
                          n_note(d$case_id, "cases")),
         x = NULL, y = "robustness (fraction of landmarks improved)",
         caption = ANHIR_CAPTION)
}

# The challenge's primary ranking metric, subset == "all": each method's rank on
# per-case median rTRE, averaged over the cases it was scored on. Best on top.
plot_anhir_rank <- function(aggregates) {
  title <- "ANHIR: mean rank of per-case median rTRE"
  d <- aggregates %>%
    dplyr::filter(subset == "all", is.finite(avg_rank_median_rtre)) %>%
    .anhir_label_methods()
  if (!nrow(d)) return(.anhir_empty(title, "ranked methods (subset == \"all\")"))
  # coord_flip() draws the LAST level at the top, so reverse the order by rank.
  d$method_lab <- stats::reorder(d$method_lab, -d$avg_rank_median_rtre)
  n_by <- stats::setNames(d$n_cases, as.character(d$method_lab))
  ggplot(d, aes(method_lab, avg_rank_median_rtre, fill = method)) +
    geom_col(width = .6) +
    geom_text(aes(label = sprintf("%.2f", avg_rank_median_rtre)),
              hjust = -0.15, size = pt_text(7)) +
    coord_flip() +
    scale_fill_anhir_method(guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, .15))) +
    # n here is the aggregate's own n_cases, not a row count: one row per method.
    scale_x_discrete(labels = function(b) sprintf("%s (n = %d)", b, as.integer(n_by[b]))) +
    labs(title = title,
         subtitle = paste("The challenge's primary metric: the method's rank among methods",
                          "on each case, averaged over cases. 1 = best; lower = better."),
         x = NULL, y = "mean rank of per-case median rTRE (lower = better)",
         caption = ANHIR_CAPTION)
}

# --- table --------------------------------------------------------------------
# One row per (subset, method) for the `all` and `training` subsets, methods in
# rank order within a subset, for knitr::kable(). Column names are plain so the
# page can pass `digits` without renaming.
anhir_summary_table <- function(aggregates) {
  aggregates %>%
    dplyr::filter(subset %in% c("all", "training")) %>%
    dplyr::mutate(subset = factor(subset, levels = c("all", "training")),
                  method = .anhir_method_factor(method)) %>%
    dplyr::arrange(subset, avg_rank_median_rtre, method) %>%
    dplyr::transmute(subset         = as.character(subset),
                     method         = .anhir_method_labels(method),
                     n_cases,
                     mean_rank      = avg_rank_median_rtre,
                     avg_median_rtre,
                     med_median_rtre,
                     avg_max_rtre,
                     avg_robustness,
                     avg_time_min)
}
