# =============================================================================
# placeholders.R  —  OPT-IN synthetic stand-ins for benchmark points not yet run
#
# WHY THIS EXISTS. The registration arms, the synthetic sweep and the ANHIR legs run
# on the cluster for days, and every page renders only what has arrived: an arm
# with no QC yet is simply absent from every panel. To SEE the final shape of a
# figure while the runs are in flight, this file can fill each expected-but-missing
# point with a synthetic value drawn from the points that DID arrive.
#
# THE INTEGRITY GUARD — every rule below is load-bearing, none is cosmetic.
#
#   1. OFF BY DEFAULT. Nothing here does anything unless
#        options(ihc.placeholder_missing = TRUE)       or
#        IHC_PLACEHOLDER_MISSING=1                    (env var, so `wflow_build`
#                                                       can be driven from a shell)
#      With the mode off, placeholder_fill() returns its input UNCHANGED — not even a
#      flag column is added — so a default knit is byte-identical to one made before
#      this file existed. An explicitly set option wins over the env var, so
#      options(ihc.placeholder_missing = FALSE) switches it off for one session.
#   2. REAL ROWS ARE NEVER EDITED. A missing value is filled on a NEW row
#      (`is_placeholder = TRUE`); a real row whose metric is NA keeps its NA, and the
#      synthetic value rides on a companion row carrying only the filled metric(s).
#   3. EVERY SYNTHETIC ROW SAYS HOW IT WAS MADE: `placeholder_rule` names, per
#      metric, the fallback level the value came from and the n behind it.
#   4. EVERY FIGURE THAT DRAWS ONE SAYS SO: placeholder_style() draws synthetic
#      points hollow, synthetic lines dashed, synthetic bars/tiles faded, stamps a
#      PLACEHOLDER watermark in each panel and prefixes the subtitle with the count.
#   5. NOTHING SYNTHETIC CAN BECOME A PAPER FIGURE. export_pdf_figures() (pdf_export.R)
#      refuses to write output/figures/<slug>/ while the mode is on and diverts to
#      output/placeholders/figures/<slug>/*_PLACEHOLDER.pdf; figures/_common.R's
#      exporters stop() outright; and every synthetic point is listed in the sidecar
#      output/placeholders/<slug>_placeholders.csv.
#
# THE SYNTHESIS MODEL, in order — the first level with any real value wins:
#   caller-supplied `levels`  e.g. same arm (same stage) -> same micro depth -> same
#                             tier -> same backend -> same stage. mean ± sd of the
#                             real values at that level, one normal draw.
#   global                    every real value of that metric.
#   prior                     a FIXED, documented guess (PLACEHOLDER_PRIORS below),
#                             used only when no real value exists anywhere.
# Metrics named in `log_scale` are modelled on log() so a positive, skewed quantity
# (a residual in µm, an rTRE) stays positive and keeps its spread; everything is then
# clipped to `ranges`. A level with n = 1 has no sd; its spread falls back to the
# prior's sd (else 10% of the mean), so one arrived point is not copied verbatim.
#
# DETERMINISTIC. The draw uses its own seed (option `ihc.placeholder_seed`, else env
# IHC_PLACEHOLDER_SEED, else PLACEHOLDER_DEFAULT_SEED) offset per table, and restores
# the caller's RNG state afterwards, so turning the mode on never perturbs a jitter
# or a bootstrap elsewhere on the page.
#
# Dependencies: base R + tibble + dplyr (namespaced — nothing is attached), plus
# ggplot2 for the plot helper only. Sourcing this file loads no package, which is
# what lets pdf_export.R source it and stay dependency-light. Testable with no data
# on disk: tests/testthat/test-placeholders.R.
# =============================================================================

PLACEHOLDER_DEFAULT_SEED <- 20260927L
PLACEHOLDER_TAG          <- "PLACEHOLDER"

# The fixed prior, used ONLY when a metric has no real value anywhere in its table.
# These are round guesses of the right order of magnitude, not estimates; a figure
# drawn from them says nothing except where the points will land on the axis. For a
# metric listed in `log_scale`, `mean` and `sd` are on the natural-log scale.
PLACEHOLDER_PRIORS <- list(
  # segmentation-overlap QC (registration_arms.R, registration_accuracy_plots.R)
  disp_um_p50             = list(mean = log(2),     sd = 0.7),
  displacement_um_p50     = list(mean = log(2),     sd = 0.7),
  displacement_um_p90     = list(mean = log(5),     sd = 0.7),
  reg_displacement_um_p50 = list(mean = log(2),     sd = 0.7),
  dice_matched            = list(mean = 0.6,        sd = 0.1),
  reg_dice_matched        = list(mean = 0.6,        sd = 0.1),
  pair_fraction           = list(mean = 0.7,        sd = 0.1),
  d_disp_um_vs_rigid      = list(mean = -0.5,       sd = 1),
  # VALIS's own error (fraction of the diagonal, or a distance) and STARE's (pixels)
  rTRE                    = list(mean = log(0.01),  sd = 0.8),
  D                       = list(mean = log(20),    sd = 0.8),
  tre_px                  = list(mean = log(5),     sd = 0.8),
  # ANHIR
  rtre_median             = list(mean = log(0.01),  sd = 1),
  rtre_mean               = list(mean = log(0.012), sd = 1),
  rtre_max                = list(mean = log(0.04),  sd = 1),
  tre_median_px           = list(mean = log(40),    sd = 1),
  robustness              = list(mean = 0.7,        sd = 0.15),
  # cost
  cpu_hours               = list(mean = log(5),     sd = 0.8))

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# --- the switch ----------------------------------------------------------------
.ph_truthy <- function(x) tolower(trimws(as.character(x))) %in% c("1", "true", "yes", "on")

placeholder_mode <- function() {
  opt <- getOption("ihc.placeholder_missing", NULL)
  if (!is.null(opt)) return(isTRUE(as.logical(opt)) || .ph_truthy(opt))
  .ph_truthy(Sys.getenv("IHC_PLACEHOLDER_MISSING", ""))
}

placeholder_seed <- function() {
  s <- getOption("ihc.placeholder_seed", NULL)
  if (is.null(s)) s <- Sys.getenv("IHC_PLACEHOLDER_SEED", "")
  s <- suppressWarnings(as.integer(s))
  if (length(s) != 1 || is.na(s)) PLACEHOLDER_DEFAULT_SEED else s
}

# --- the registry: what was synthesised this session, for the sidecar -----------
# Re-sourcing this file must not wipe what an earlier source already recorded.
if (!exists(".ph_registry") || !is.environment(.ph_registry)) {
  .ph_registry <- new.env(parent = emptyenv())
  .ph_registry$rows <- list()
}

placeholder_reset <- function() { .ph_registry$rows <- list(); invisible(NULL) }

placeholder_registry <- function() {
  if (!length(.ph_registry$rows))
    return(tibble::tibble(table = character(), key = character(), metric = character(),
                          value = numeric(), rule = character()))
  dplyr::distinct(dplyr::bind_rows(.ph_registry$rows))
}

# --- helpers --------------------------------------------------------------------
# Run `expr` under `seed` and put the caller's RNG state back afterwards.
.ph_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
            rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  force(expr)
}

# A stable per-table offset, so two tables filled in one knit do not draw the same
# noise sequence.
.ph_table_seed <- function(seed, what) {
  u <- utf8ToInt(what)
  as.integer((as.numeric(seed) + sum(u * seq_along(u))) %% .Machine$integer.max)
}

# Rows of `pool` whose `cols` equal `row`'s. NA matches NA: a tiled arm's NA micro
# depth is a real value ("no depth"), not an unknown one.
.ph_match <- function(pool, row, cols) {
  keep <- rep(TRUE, nrow(pool))
  for (c in cols) {
    a <- as.character(pool[[c]])
    b <- as.character(row[[c]][1])
    keep <- keep & if (is.na(b)) is.na(a) else (!is.na(a) & a == b)
  }
  keep
}

.ph_prior_for <- function(m, priors) {
  if (!is.null(priors[[m]])) return(priors[[m]])
  # VALIS columns are named per stage (rigid_rTRE, non_rigid_D, ...) and STARE's per
  # percentile (rigid_p50, after_p90): fall back to the unit's prior.
  if (grepl("_rTRE$", m)) return(priors[["rTRE"]])
  if (grepl("_D$", m))    return(priors[["D"]])
  if (grepl("^(rigid|after|coarse)", m)) return(priors[["tre_px"]])
  NULL
}

# One synthetic value for metric `m` of `row`, walking the fallback levels.
.ph_draw <- function(row, m, real, levels, ranges, log_scale, priors, min_n) {
  lg  <- m %in% log_scale
  val <- suppressWarnings(as.numeric(real[[m]]))
  ok  <- is.finite(val) & (!lg | val > 0)
  mu <- NA_real_; s <- NA_real_; rule <- NULL
  for (nm in names(levels)) {
    cols <- levels[[nm]]
    if (!all(cols %in% names(real)) || !all(cols %in% names(row))) next
    v <- val[ok & .ph_match(real, row, cols)]
    if (!length(v)) next
    v  <- if (lg) log(v) else v
    mu <- mean(v); s <- if (length(v) >= min_n) stats::sd(v) else NA_real_
    rule <- sprintf("%s %s mean±sd (n=%d)", m, nm, length(v))
    break
  }
  if (is.null(rule) && any(ok)) {
    v  <- if (lg) log(val[ok]) else val[ok]
    mu <- mean(v); s <- if (length(v) >= min_n) stats::sd(v) else NA_real_
    rule <- sprintf("%s global mean±sd (n=%d)", m, length(v))
  }
  pr <- .ph_prior_for(m, priors)
  if (is.null(rule)) {
    if (is.null(pr))
      return(list(value = NA_real_, rule = sprintf("%s none (no data, no prior)", m)))
    mu <- pr$mean; s <- pr$sd
    rule <- sprintf("%s fixed prior", m)
  }
  if (!is.finite(s) || s == 0)
    s <- if (!is.null(pr)) pr$sd else if (lg) 0.1 else 0.1 * abs(mu) + 1e-9
  x <- stats::rnorm(1, mu, s)
  if (lg) x <- exp(x)
  r <- ranges[[m]]
  if (!is.null(r)) x <- min(max(x, r[1]), r[2])
  list(value = x, rule = rule)
}

# Is metric `m` one that rows shaped like `row` actually carry? A metric absent by
# construction (the delta-vs-rigid of the rigid stage itself) must stay NA, or the
# placeholder would invent a quantity the producer never defines.
.ph_expected_metric <- function(row, m, real, na_structure) {
  ok <- is.finite(suppressWarnings(as.numeric(real[[m]])))
  if (!any(ok)) return(FALSE)
  if (is.null(na_structure) || !length(na_structure)) return(TRUE)
  cols <- intersect(na_structure, intersect(names(real), names(row)))
  if (!length(cols)) return(TRUE)
  any(ok & .ph_match(real, row, cols))
}

.ph_key_string <- function(d, keys) {
  if (!nrow(d)) return(character(0))
  do.call(paste, c(lapply(keys, function(k) paste0(k, "=", as.character(d[[k]]))),
                   list(sep = "; ")))
}

# --- the core -------------------------------------------------------------------
# Fill expected-but-missing rows (and NA metrics) of `real` with synthetic values.
#
#   real          the loaded table, exactly as the loader built it.
#   expected      one row per EXPECTED key, carrying the key columns plus any
#                 descriptor columns the synthetic rows should have (arm, backend,
#                 stage_index, ...). Rows whose key is already in `real` are ignored.
#   keys          the columns that identify one point.
#   metrics       numeric columns to synthesise; only those present in `real` are.
#   levels        NAMED list of grouping-column vectors, most specific first. The
#                 name is what `placeholder_rule` reports.
#   ranges        named list of c(lo, hi) clips.
#   log_scale     metrics modelled on log().
#   priors        named list of list(mean, sd); defaults to PLACEHOLDER_PRIORS.
#   na_structure  columns within which "this metric exists" is decided (see
#                 .ph_expected_metric); NULL = anywhere in the table.
#   what          the table's name, for the sidecar and the seed offset.
#
# Returns `real` UNCHANGED when the mode is off; otherwise `real` plus two columns
# (`is_placeholder` FALSE, `placeholder_rule` NA), bound to the synthetic rows.
placeholder_fill <- function(real, expected, keys, metrics, levels = list(),
                             ranges = list(), log_scale = character(0),
                             priors = PLACEHOLDER_PRIORS, na_structure = NULL,
                             what = "table", seed = placeholder_seed(),
                             enabled = placeholder_mode(), min_n = 2L) {
  if (!isTRUE(enabled)) return(real)
  if (is.null(real) || !nrow(real)) return(real)
  real <- tibble::as_tibble(real)
  # Idempotent: a second pass starts again from the real rows only.
  if ("is_placeholder" %in% names(real))
    real <- dplyr::select(real[!real$is_placeholder %in% TRUE, , drop = FALSE],
                          -dplyr::any_of(c("is_placeholder", "placeholder_rule")))
  if (!all(keys %in% names(real))) {
    warning("placeholder_fill(", what, "): key column(s) missing from the table: ",
            paste(setdiff(keys, names(real)), collapse = ", "), " — nothing synthesised")
    return(real)
  }
  metrics <- intersect(metrics, names(real))
  marked  <- dplyr::mutate(real, is_placeholder = FALSE, placeholder_rule = NA_character_)
  if (!length(metrics)) return(marked)

  # 1. Keys with no real row at all.
  miss <- if (!is.null(expected) && nrow(expected) && all(keys %in% names(expected)))
    dplyr::anti_join(dplyr::distinct(tibble::as_tibble(expected)),
                     dplyr::distinct(real[keys]), by = keys)
  else tibble::tibble()
  if (nrow(miss)) {
    miss <- dplyr::select(miss, -dplyr::any_of(c(metrics, "is_placeholder", "placeholder_rule")))
    miss <- dplyr::arrange(miss, dplyr::across(dplyr::all_of(keys)))
    miss$.why  <- rep("no row", nrow(miss))
    miss$.need <- lapply(seq_len(nrow(miss)), function(i)
      Filter(function(m) .ph_expected_metric(miss[i, ], m, real, na_structure), metrics))
  }

  # 2. Real rows whose metric is NA where rows of the same structure carry it. The
  #    real row stays as it is; the value goes on a companion row.
  need_na <- lapply(seq_len(nrow(real)), function(i) {
    Filter(function(m) !is.finite(suppressWarnings(as.numeric(real[[m]][i]))) &&
             .ph_expected_metric(real[i, ], m, real, na_structure), metrics)
  })
  has_na <- vapply(need_na, length, integer(1)) > 0
  comp <- real[has_na, , drop = FALSE]
  if (nrow(comp)) {
    comp$.why  <- rep("NA metric", nrow(comp))
    comp$.need <- need_na[has_na]
    for (m in metrics) comp[[m]] <- NA_real_
  }

  syn <- dplyr::bind_rows(miss, comp)
  if (!nrow(syn)) return(marked)
  for (m in metrics) syn[[m]] <- NA_real_

  reg <- list()
  syn <- .ph_with_seed(.ph_table_seed(seed, what), {
    rules <- character(nrow(syn))
    for (i in seq_len(nrow(syn))) {
      rr <- character(0)
      for (m in syn$.need[[i]]) {
        dr <- .ph_draw(syn[i, ], m, real, levels, ranges, log_scale, priors, min_n)
        syn[[m]][i] <- dr$value
        rr <- c(rr, dr$rule)
        reg[[length(reg) + 1]] <- tibble::tibble(
          table = what, key = .ph_key_string(syn[i, ], keys), metric = m,
          value = dr$value, rule = paste0(syn$.why[i], ": ", dr$rule))
      }
      rules[i] <- paste0(syn$.why[i], " | ",
                         if (length(rr)) paste(rr, collapse = "; ") else "no metric expected")
    }
    syn$placeholder_rule <- rules
    syn
  })
  syn$is_placeholder <- TRUE
  syn <- dplyr::select(syn, -".why", -".need")
  # Keep the real frame's column types (a stage axis ordered by its factor levels
  # must not turn into an alphabetical character axis because one row is synthetic).
  for (c in intersect(names(real), names(syn))) {
    if (is.factor(real[[c]]) && !is.factor(syn[[c]]))
      syn[[c]] <- factor(as.character(syn[[c]]), levels = levels(real[[c]]))
    else if (is.integer(real[[c]]) && !is.integer(syn[[c]]) &&
             (is.numeric(syn[[c]]) || is.logical(syn[[c]])))
      syn[[c]] <- as.integer(syn[[c]])
    else if (is.double(real[[c]]) && is.logical(syn[[c]]))
      syn[[c]] <- as.numeric(syn[[c]])
    else if (is.character(real[[c]]) && is.logical(syn[[c]]))
      syn[[c]] <- as.character(syn[[c]])
    else if (is.logical(real[[c]]) && !is.logical(syn[[c]]) && all(is.na(syn[[c]])))
      syn[[c]] <- as.logical(syn[[c]])
  }
  if (length(reg)) .ph_registry$rows <- c(.ph_registry$rows, list(dplyr::bind_rows(reg)))
  message(sprintf(
    "placeholders: %s — %d synthetic row(s) (%d with no real row, %d NA-metric companions)",
    what, nrow(syn), nrow(miss), nrow(comp)))
  dplyr::bind_rows(marked, syn)
}

placeholder_count <- function(d) {
  if (is.null(d) || !is.data.frame(d) || !"is_placeholder" %in% names(d)) return(0L)
  sum(d$is_placeholder %in% TRUE)
}

# --- output guards --------------------------------------------------------------
placeholder_root <- function() here::here("output", "placeholders")

# Every synthetic value this session made, one row per (table, key, metric).
placeholder_write_sidecar <- function(slug, root = placeholder_root()) {
  if (!placeholder_mode()) return(invisible(NULL))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(root, paste0(slug, "_placeholders.csv"))
  reg  <- placeholder_registry()
  utils::write.csv(reg, path, row.names = FALSE)
  message(sprintf("placeholders: %d synthetic value(s) listed in %s", nrow(reg), path))
  invisible(path)
}

# Where a derived table (Additional file 2, ...) may be written. In placeholder mode
# the file moves under output/placeholders/ with a _PLACEHOLDER suffix, so a table
# carrying synthetic medians can never overwrite the real one.
placeholder_output_path <- function(path, root = placeholder_root()) {
  if (!placeholder_mode()) return(path)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  ext  <- tools::file_ext(path)
  stem <- tools::file_path_sans_ext(basename(path))
  file.path(root, paste0(stem, "_", PLACEHOLDER_TAG, if (nzchar(ext)) paste0(".", ext)))
}

# For exporters that write manuscript figures (figures/_common.R): there is no
# placeholder version of a paper figure, so they refuse.
placeholder_refuse <- function(what) {
  if (placeholder_mode())
    stop("placeholder mode is ON (ihc.placeholder_missing / IHC_PLACEHOLDER_MISSING): ",
         "refusing to write ", what, ". A paper figure must not carry synthetic points — ",
         "unset the switch and re-run.", call. = FALSE)
  invisible(TRUE)
}

# The page-level warning. Empty string when the mode is off, so an unconditional
# `cat(placeholder_callout())` in a results="asis" chunk costs a default page nothing.
placeholder_callout <- function() {
  if (!placeholder_mode()) return("")
  paste0(
    '\n\n<div style="border:3px solid #B2182B;background:#FDECEA;padding:0.8em 1em;',
    'margin:1em 0;">\n\n',
    "**", PLACEHOLDER_TAG, " MODE IS ON — this page is NOT a result.**\n\n",
    "Benchmark points that have not arrived yet are filled with SYNTHETIC values drawn ",
    "from the points that did (same arm, then neighbouring arms, then the metric's global ",
    "distribution, then a fixed prior — `code/placeholders.R`). Synthetic points are drawn ",
    "hollow / dashed / faded and every affected figure carries a PLACEHOLDER watermark; ",
    "aggregated tables carry an `n_placeholder` column. PDFs go to ",
    "`output/placeholders/figures/` with a `_PLACEHOLDER` suffix, never to ",
    "`output/figures/`, and every synthetic value is listed in ",
    "`output/placeholders/<page>_placeholders.csv`. Seed: ", placeholder_seed(), ". ",
    "Turn it off by unsetting `IHC_PLACEHOLDER_MISSING` / ",
    "`options(ihc.placeholder_missing = FALSE)`.\n\n</div>\n\n")
}

# --- the figure helper ------------------------------------------------------------
# How many synthetic rows a plot draws: the plot data's, else the largest layer's.
placeholder_plot_count <- function(p) {
  n <- placeholder_count(p$data)
  if (n > 0) return(n)
  for (ly in p$layers) if (is.data.frame(ly$data)) n <- max(n, placeholder_count(ly$data))
  n
}

.ph_keep <- function(d, want) {
  if (!is.data.frame(d)) return(d)
  if (!"is_placeholder" %in% names(d)) return(if (want) d[0, , drop = FALSE] else d)
  d[(d$is_placeholder %in% TRUE) == want, , drop = FALSE]
}

# A child of `layer` drawing only the real (want = FALSE) or synthetic rows. The
# original layer object is never modified: ggproto objects are environments, and a
# plot printed twice must not be filtered twice.
.ph_layer_subset <- function(layer, want, params = list()) {
  ch   <- ggplot2::ggproto(NULL, layer)
  base <- layer$data
  ch$data <- if (is.data.frame(base)) .ph_keep(base, want)
             else if (is.function(base)) function(d) .ph_keep(base(d), want)
             else function(d) .ph_keep(d, want)
  if (length(params)) ch$aes_params <- utils::modifyList(layer$aes_params %||% list(), params)
  ch
}

# Mark a plot's synthetic points. A no-op (returns `p` itself) when the plot draws
# none, so builders call it unconditionally.
#   points (incl. jitter)  -> hollow circle, faded
#   lines / paths          -> dashed, faded
#   bars / tiles / rects   -> faded
#   boxplots, text         -> unchanged: they summarise real + synthetic together,
#                             which the watermark and the subtitle state
placeholder_style <- function(p) {
  if (is.null(p) || !inherits(p, "ggplot")) return(p)
  n <- placeholder_plot_count(p)
  if (n == 0) return(p)
  old_sub <- p$labels$subtitle
  if (is.character(old_sub) && startsWith(old_sub[1], PLACEHOLDER_TAG)) return(p)  # styled

  out <- list()
  for (ly in p$layers) {
    g <- ly$geom
    params <- if (inherits(g, "GeomPoint")) list(shape = 1, alpha = 0.6, stroke = 0.7)
              else if (inherits(g, "GeomPath")) list(linetype = "22", alpha = 0.5)
              else if (inherits(g, "GeomRect")) list(alpha = 0.35)
              else NULL
    if (is.null(params)) { out <- c(out, list(ly)); next }
    out <- c(out, list(.ph_layer_subset(ly, FALSE), .ph_layer_subset(ly, TRUE, params)))
  }
  p$layers <- out
  lab <- sprintf("%s — %d synthetic point%s (hollow / dashed / faded) — NOT a result",
                 PLACEHOLDER_TAG, n, if (n == 1) "" else "s")
  p +
    ggplot2::annotate("text", x = Inf, y = Inf, label = PLACEHOLDER_TAG,
                      hjust = 1.05, vjust = 1.4, colour = "#B2182B", alpha = 0.55,
                      fontface = "bold", size = 3.2) +
    ggplot2::labs(
      subtitle = if (is.null(old_sub)) lab else paste0(lab, "\n", old_sub),
      caption  = paste(c(p$labels$caption,
                         "Synthetic points: code/placeholders.R; listed in output/placeholders/."),
                       collapse = "\n"))
}
