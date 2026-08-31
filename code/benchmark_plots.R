# ============================================================================
# benchmark_plots.R  —  Mirage benchmark figures (vendored fork)
# ============================================================================
# Upstream: ../mirage `benchmarks/analysis/plots.R` on branch bench/reconcile-main
# (the branch that carries benchmarks/ — main does not). This is a FORK, not a
# mirror: re-vendoring means merging, not copying. Deliberate differences:
#   1. default `adir` reads from data/benchmark/  (drop the sweep CSVs there)
#   2. save_fig collects each ggplot into the in-memory `bench_figs` list (no
#      files on disk); benchmarks.Rmd sources this file and renders them inline
#   3. the house style comes from code/plot_theme.R, so these figures cannot
#      drift away from the validation figures
#   4. upstream figs 03/03b/10/15/16 are dropped — they compare the classic and
#      distributed registration paths, and mirage now has a single path
#   5. upstream figs 11 and 17 are dropped here and rendered by
#      code/registration_accuracy_plots.R (§4, §5) instead: same CSV, same
#      columns, one owner. Do not re-add them when merging upstream.
#
# INPUT CSVs — data/benchmark/ is a merge of TWO mirage output sets:
#   from `python -m benchmarks.analysis.make_figures` (writes benchmarks/analysis/)
#     measurements.csv    REQUIRED. One row per (run x PROCESS): peak_rss_gb,
#                         peak_vmem_gb, realtime_s, duration_s, cpus, input_gb,
#                         read_gb, write_gb + every swept param.
#     resource_stats.csv  per (process, config): n_reps + mean/std/cv
#     run_cost.csv        per run: cpu_hours, gpu_hours, wall_clock_s, bottleneck_stage
#   from `python benchmarks/analysis/make_tables.py` (writes benchmarks/paper_data/)
#     runs_master.csv, param_matrix.csv, segmentation_agreement.csv
#     (+ registration_accuracy.csv / registration_valis_rtre.csv / scaling_fits.csv /
#      segmentation_eval.csv, which the registration_accuracy page reads)
# Only measurements.csv is required; a figure whose CSV or column is absent is skipped.
# ============================================================================
.need <- c("ggplot2", "dplyr", "readr", "tidyr", "stringr", "forcats", "purrr", "scales")
.missing <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing))
  stop("Missing R packages: ", paste(.missing, collapse = ", "),
       "\n  install.packages(c(", paste(sprintf('"%s"', .missing), collapse = ", "), "))",
       call. = FALSE)
suppressPackageStartupMessages(lapply(.need, library, character.only = TRUE))

adir  <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else here::here("data", "benchmark")
CAPTION <- "Mirage benchmark sweep · mean over replicate runs · SLURM-isolated per-process resources"
# Collect each figure into a named list of ggplot objects. benchmarks.Rmd sources
# this file and renders them INLINE from the data (no PNG files on disk).
bench_figs <- list()
# NOTE — no width/height arguments, deliberately. This used to be
# `save_fig(p, name, w = 8, h = 5)` and every call site passed a size (8x5, 11x8,
# 12x8...), but the body never read them: the figures are collected into a list
# and drawn by the Rmd chunk, so the size that actually applies is the chunk's
# fig.width/fig.height. Twelve call sites therefore documented sizes that were not
# the sizes, and "make this figure wider" edits here did nothing. The real knob is
# fig_width() in analysis/benchmark_pipeline.Rmd.
save_fig <- function(p, name) {
  bench_figs[[name]] <<- p + labs(caption = CAPTION)
  invisible(NULL)
}

# House style. This file used to carry its own copy of the theme and the `oi`
# palette (vendored from mirage); both now live in code/plot_theme.R so the
# benchmark figures and the validation figures cannot drift apart. Sourcing it
# applies theme_set() + the geom defaults and defines `oi`.
source(here::here("code", "plot_theme.R"))

m <- read_csv(file.path(adir, "measurements.csv"), show_col_types = FALSE) %>%
  mutate(proc = str_replace(process, ".*:", ""),                 # leaf process name
         input_gb = as.numeric(input_gb))
size_axes <- c("baseline", "scaling_grid", "registration_grid", "distributed_grid",
               "target_px", "n_channels")
# I/O volume per process (read+write GiB) — present only if the trace carried rchar/wchar
# (load.py parses them into read_gb/write_gb; older CSVs won't have the columns).
has_io <- all(c("read_gb", "write_gb") %in% names(m))
if (has_io) m <- m %>% mutate(total_io_gb = read_gb + write_gb)

# POWER-LAW fit per process: lm(log10(y) ~ log10(x)). The slope beta is the scaling exponent (beta=1 linear,
# >1 super-linear, <1 sub-linear) — the paper number. beta/R2 are surfaced in each facet strip and the
# fit is drawn as a curve on LINEAR axes (no log-log), so `powerlaw` returns a fine x-grid, not just
# the two endpoints (a power law is a straight line only in log-log space; on linear axes it curves).
powerlaw <- function(df, xcol, ycol) {
  parts <- lapply(split(df, df$proc), function(d) {
    if (length(unique(d[[xcol]])) < 2) return(NULL)
    f <- lm(log10(d[[ycol]]) ~ log10(d[[xcol]]))
    b <- unname(coef(f)[2]); a <- unname(coef(f)[1]); r2 <- summary(f)$r.squared
    xr <- range(d[[xcol]])
    xs <- 10 ^ seq(log10(xr[1]), log10(xr[2]), length.out = 80)   # smooth curve for linear axes
    data.frame(proc = d$proc[1], x = xs, y = 10 ^ (a + b * log10(xs)), exponent = b, r2 = r2)
  })
  do.call(rbind, parts)
}
powerlaw_plot <- function(df, ycol, point_col, title, ylab) {
  d <- df %>% filter(varied_axis %in% size_axes, is.finite(input_gb), input_gb > 0, .data[[ycol]] > 0)
  pl <- powerlaw(d, "input_gb", ycol)
  if (is.null(pl) || !nrow(pl)) return(NULL)
  # beta/R2 go into the facet strip label (declutters the panel: no overlapping in-panel text box).
  strip <- pl %>% group_by(proc) %>%
    summarise(l = sprintf("%s  (beta=%.2f, R2=%.2f)", first(proc), first(exponent), first(r2)),
              .groups = "drop")
  lookup <- setNames(strip$l, strip$proc)
  relabel <- function(v) ifelse(is.na(lookup[v]), v, lookup[v])   # procs with no fit keep bare name
  ggplot(d, aes(input_gb, .data[[ycol]])) +
    geom_point(alpha = .6, colour = point_col) +
    geom_line(data = pl, aes(x, y), colour = oi[2], linewidth = .6) +
    # free_y, not free: every process was run on the SAME set of input sizes, so a
    # free x gave each panel a different span of the same axis and made the sweeps
    # look like they covered different experiments. y stays free because peak RSS
    # and runtime differ by orders of magnitude between processes, and the panel is
    # read for the SHAPE of the curve (beta, in the strip), which a shared y flattens.
    facet_wrap(~ proc, scales = "free_y", labeller = labeller(proc = relabel)) +
    labs(title = title,
         subtitle = "Linear axes; curve = fitted power law. beta = log-log slope (1 = linear, >1 super-linear, <1 sub-linear).",
         x = "input (GiB)", y = ylab)
}

# ── 1. MEMORY SCALING per process (the headline) — peak RSS vs input, power law (linear axes) ──
save_fig(powerlaw_plot(m, "peak_rss_gb", oi[1], "Peak memory scaling per process (power law)",
                       "peak RSS (GiB)"), "01_memory_scaling_per_process")

# ── 2. TIME SCALING per process — realtime vs input, power law (linear axes) ──
save_fig(powerlaw_plot(m, "realtime_s", oi[3], "Runtime scaling per process (power law)",
                       "realtime (s)"), "02_time_scaling_per_process")

# ── 2b. I/O VOLUME SCALING per process — bytes moved (read+write) vs input, power law ──
if (has_io && any(is.finite(m$total_io_gb) & m$total_io_gb > 0)) {
  io_fig <- powerlaw_plot(m, "total_io_gb", oi[6],
                          "I/O volume scaling per process (power law)", "read + write (GiB)")
  if (!is.null(io_fig)) save_fig(io_fig, "02b_io_volume_scaling")
}

# ── 4. N-IMAGE REGISTRATION — REGISTER cost vs number of slides ──
reg <- m %>% filter(proc == "REGISTER", varied_axis %in% c("registration_grid", "baseline", "scaling_grid"))
if (nrow(reg) > 0) {
  p4 <- reg %>% group_by(target_px, n_channels, n_register_images) %>%
    summarise(peak_rss_gb = mean(peak_rss_gb), realtime_s = mean(realtime_s), .groups = "drop") %>%
    ggplot(aes(n_register_images, peak_rss_gb, colour = factor(target_px))) +
    geom_line() + geom_point(size = 2) +
    facet_wrap(~ n_channels, labeller = label_both) +
    scale_colour_ordinal(name = "size (px)") +
    labs(title = "N-image registration: peak RAM vs slide count",
         subtitle = "Co-registering more slides to one reference; coloured by image size.",
         x = "Slides co-registered (n; 1 reference + N-1 moving)",
         y = "Peak RSS (GiB)")
  save_fig(p4, "04_nimage_registration_ram")
}

# ── 5. OFAT KNOB EFFECTS — one panel per single-knob axis ──
# For each OFAT axis, plot the most-affected process's realtime vs the knob value.
# Only the true single-knob OFAT axes belong here; memory_mode / skip_micro_registration have no
# dedicated figure (the classic/distributed-path comparison that read them was retired) and the
# segmentation tile knobs go to plots 9/9b (per method).
knob_targets <- tribble(
  ~axis,                       ~proc,          ~metric,
  "preproc_n_iter",            "PREPROCESS",   "realtime_s",
  "preproc_overlap",           "PREPROCESS",   "realtime_s",
  "preproc_pool_workers",      "PREPROCESS",   "realtime_s",
  "seg_gpu",                   "SEGMENT",      "realtime_s",
  "quantify_compartments",     "QUANTIFY",     "realtime_s",
  "expanded_quantification",   "QUANTIFY",     "realtime_s"
)
knob_df <- pmap_dfr(knob_targets, function(axis, proc, metric) {
  if (!axis %in% names(m)) return(NULL)
  m %>% filter(varied_axis %in% c(axis, "baseline"), proc == !!proc) %>%
    transmute(axis = axis, proc = proc,
              value = as.character(.data[[axis]]), y = .data[[metric]], metric = metric)
})
if (nrow(knob_df) > 0) {
  p5 <- knob_df %>% group_by(axis, proc, metric, value) %>%
    summarise(y = mean(y), .groups = "drop") %>%
    ggplot(aes(fct_inseq(value), y)) +
    geom_col(fill = oi[1], width = .6) +
    # free_x, the mirror image of the scaling plots above: every panel measures the
    # same quantity (mean realtime in seconds), so y is SHARED and the panels can be
    # ranked against each other — which knob actually costs is the question this
    # figure exists to answer, and a free y made every knob look equally expensive.
    # x must stay free: the knob VALUES differ per panel (iteration counts vs a
    # TRUE/FALSE), and a shared discrete x would draw the union of all of them in
    # every panel as empty slots.
    facet_wrap(~ paste0(axis, "  (", proc, ")"), scales = "free_x", ncol = 3) +
    labs(title = "OFAT knob effects (single param varied off baseline)",
         subtitle = "Mean realtime (s) per knob value; each panel is one knob. Shared y, so knobs are compared on one runtime scale; x is per-knob because the knob values differ.",
         x = NULL, y = "realtime (s)") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(p5, "05_ofat_knob_effects")
}

# ── 6. REPLICATE VARIANCE — mean +/- sd from the repeats ──
stats_path <- file.path(adir, "resource_stats.csv")
if (file.exists(stats_path)) {
  st <- read_csv(stats_path, show_col_types = FALSE)
  if (nrow(st) > 0 && "peak_rss_gb_mean" %in% names(st)) {
    p6 <- st %>% mutate(proc = str_replace(process, ".*:", "")) %>%
      filter(n_reps > 1) %>%
      ggplot(aes(reorder(proc, peak_rss_gb_mean), peak_rss_gb_mean)) +
      geom_col(fill = oi[6], width = .6) +
      geom_errorbar(aes(ymin = peak_rss_gb_mean - peak_rss_gb_std,
                        ymax = peak_rss_gb_mean + peak_rss_gb_std), width = .3) +
      coord_flip() +
      labs(title = "Peak RSS by process (mean +/- sd across repeats)", x = NULL, y = "peak RSS (GiB)")
    save_fig(p6, "06_replicate_variance")
  }
}

# ── 7. STAGE-COST HEATMAP — process x size, fill = peak RSS ──
p7 <- m %>% filter(varied_axis %in% size_axes, n_channels == 2, n_register_images == 2) %>%
  group_by(proc, target_px) %>% summarise(peak_rss_gb = mean(peak_rss_gb), .groups = "drop") %>%
  ggplot(aes(factor(target_px), fct_reorder(proc, peak_rss_gb), fill = peak_rss_gb)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_seq(transform = "log10", name = "peak RSS (GiB)", guide = guide_cbar()) +
  labs(title = "Where the memory goes",
       subtitle = "Peak RSS by stage x image size (log colour). Darker = the memory bottleneck at that size.",
       x = "image size (px)", y = NULL)
save_fig(p7, "07_stage_memory_heatmap")

# ── 7b. I/O VOLUME by stage — mean bytes read vs written per process (which stage is I/O-heavy) ──
if (has_io) {
  io_stage <- m %>% filter(varied_axis %in% size_axes) %>%
    group_by(proc) %>%
    summarise(read = mean(read_gb, na.rm = TRUE), write = mean(write_gb, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_longer(c(read, write), names_to = "direction", values_to = "gb") %>%
    filter(is.finite(gb))
  if (nrow(io_stage) > 0 && any(io_stage$gb > 0)) {
    p7b <- io_stage %>%
      ggplot(aes(fct_reorder(proc, gb, .fun = sum), gb, fill = direction)) +
      geom_col(width = .6) + coord_flip() +
      scale_fill_manual(values = oi[c(1, 2)], name = NULL) +
      labs(title = "I/O volume by stage",
           subtitle = "Mean bytes read/written per process (trace rchar/wchar) — the I/O bottleneck, stacked read + write.",
           x = NULL, y = "I/O volume (GiB)")
    save_fig(p7b, "07b_stage_io_split")
  }
}

# ── 8. CHANNEL EFFECT — does 2 vs 4 channels shift the memory scaling? ──
p8 <- m %>% filter(varied_axis %in% size_axes, is.finite(input_gb), input_gb > 0,
                   proc %in% c("REGISTER","PREPROCESS","SEGMENT","QUANTIFY")) %>%
  ggplot(aes(input_gb, peak_rss_gb, colour = factor(n_channels))) +
  geom_point(alpha = .6) + geom_smooth(method = "lm", se = FALSE, linewidth = .6, formula = y ~ x) +
  # Shared x (the same input sizes were run for every process), free y (peak RSS
  # differs by orders of magnitude between processes and the reading is the slope).
  facet_wrap(~ proc, scales = "free_y") +
  scale_colour_manual(values = oi[c(1,2)], name = "channels") +
  labs(title = "Channel-count effect on memory scaling", x = "input (GiB)", y = "peak RSS (GiB)")
save_fig(p8, "08_channel_effect")

# ── 9. SEGMENTATION METHODS — each backend with its own parameter sweep ──
seg <- m %>% filter(str_starts(varied_axis, "segmentation_grid"), proc == "SEGMENT")
if (nrow(seg) > 0) {
  p9 <- seg %>%
    ggplot(aes(seg_method, realtime_s, colour = seg_method)) +
    geom_boxplot(outlier.shape = NA, width = .5) +
    scale_x_discrete(labels = label_n(seg$seg_method)) +
    geom_jitter(width = .12, alpha = .5, size = 1) +
    scale_colour_manual(values = oi, guide = "none") +
    labs(title = "Segmentation methods compared",
         subtitle = "Box = IQR across each method's own parameter sweep; points = individual configs.",
         x = NULL, y = "SEGMENT realtime (s)")
  save_fig(p9, "09_segmentation_methods")

  # StarDist tile grid effect (its own params)
  sd <- seg %>% filter(seg_method == "stardist")
  if (nrow(sd) > 0) {
    p9b <- sd %>% group_by(seg_n_tiles_x, seg_n_tiles_y) %>%
      summarise(peak_rss_gb = mean(peak_rss_gb), .groups = "drop") %>%
      ggplot(aes(factor(seg_n_tiles_x), factor(seg_n_tiles_y), fill = peak_rss_gb)) +
      geom_tile(colour = "white", linewidth = 0.3) +
      scale_fill_seq(name = "peak RSS (GiB)", guide = guide_cbar()) +
      labs(title = "StarDist tiling: peak RSS vs tile grid",
           x = "Tiles across image width (n)", y = "Tiles across image height (n)")
    save_fig(p9b, "09b_stardist_tile_grid")
  }
}

# Helper: read an optional analysis CSV, returning NULL if absent/empty (keeps plots robust to
# failed runs / signals the sweep didn't produce).
read_opt <- function(name) {
  p <- file.path(adir, name)
  if (!file.exists(p)) return(NULL)
  d <- suppressWarnings(read_csv(p, show_col_types = FALSE))
  if (nrow(d) == 0) NULL else d
}

# ── 11 / 17. REGISTRATION ACCURACY — deliberately NOT built here ──
# Both accuracy views (residual-vs-cost, and VALIS feature error vs DAPI-overlap Dice)
# read the same param_matrix.csv columns as registration_accuracy_plots.R §4 and §5,
# from the same data/benchmark/ directory — rendering them on both pages would put the
# identical figure under two numbers. This page owns resource, cost and segmentation;
# registration accuracy is the registration_accuracy page's subject.
pm   <- read_opt("param_matrix.csv")
cost <- read_opt("run_cost.csv")

# ── 12. SEGMENTATION cell counts by method ──
# n_cells and seg_method both live in param_matrix.csv now (make_tables carries seg_method through);
# fall back to runs_master.csv only if a stripped param_matrix lacks seg_method.
rm_tbl <- read_opt("runs_master.csv")
sc <- NULL
if (!is.null(pm) && "n_cells" %in% names(pm)) {
  sc <- pm
  if (!("seg_method" %in% names(sc)) && !is.null(rm_tbl) && "seg_method" %in% names(rm_tbl))
    sc <- sc %>% left_join(rm_tbl %>% select(any_of(c("run_id", "seg_method"))), by = "run_id")
  sc <- sc %>% filter(is.finite(n_cells))
}
if (!is.null(sc) && nrow(sc) > 0) {
  if ("seg_method" %in% names(sc)) {
    p12 <- ggplot(sc, aes(seg_method, n_cells, colour = seg_method)) +
      geom_boxplot(outlier.shape = NA, width = .5) + geom_jitter(width = .12, alpha = .5) +
      scale_x_discrete(labels = label_n(sc$seg_method)) +
      scale_colour_manual(values = oi, guide = "none") +
      labs(title = "Segmentation: cells detected per method",
           subtitle = "Spread = each method's own parameter sweep. Large gaps = methods disagree on cell count.",
           x = NULL, y = "Cells detected (n, max mask label)")
  } else {
    p12 <- ggplot(sc, aes("all runs", n_cells)) +
      geom_boxplot(outlier.shape = NA, width = .4) + geom_jitter(width = .1, alpha = .5) +
      scale_x_discrete(labels = label_n(rep("all runs", nrow(sc)))) +
      labs(title = "Segmentation: cells detected",
           subtitle = "seg_method unavailable — pooled distribution.",
           x = NULL, y = "Cells detected (n, max mask label)")
  }
  save_fig(p12, "12_segmentation_cell_counts")
}
agree <- read_opt("segmentation_agreement.csv")
if (!is.null(agree) && "instance_f1" %in% names(agree)) {
  p12b <- agree %>% mutate(pair = paste(method_a, "vs", method_b)) %>%
    ggplot(aes(pair, instance_f1)) +
    geom_col(width = .6, fill = oi[1]) +
    geom_text(aes(label = sprintf("count ratio %.2f", cell_count_ratio)), vjust = -.4, size = 3) +
    ylim(0, 1) +
    labs(title = "Segmentation cross-method agreement (instance F1)",
         subtitle = "IoU-matched per-cell F1 between methods (1 = agree on every cell); label = cell-count ratio.",
         x = NULL, y = "Instance F1, IoU-matched (unitless, 0-1)")
  save_fig(p12b, "12b_segmentation_agreement")
}

# ── 13. END-TO-END COST — CPU-hours (and wall-clock) vs image size ──
if (!is.null(cost) && "target_px" %in% names(cost)) {
  size_cost <- cost %>% filter(varied_axis %in% size_axes) %>%
    group_by(target_px) %>%
    summarise(cpu_hours = mean(cpu_hours),
              wall_clock_h = mean(wall_clock_s, na.rm = TRUE) / 3600, .groups = "drop")
  if (nrow(size_cost) > 0) {
    p13 <- size_cost %>% pivot_longer(c(cpu_hours, wall_clock_h), names_to = "metric", values_to = "hours") %>%
      filter(is.finite(hours)) %>%
      ggplot(aes(target_px, hours, colour = metric)) +
      geom_line(linewidth = .8) + geom_point(size = 2) +
      scale_colour_manual(values = oi[c(1, 2)],
        labels = c(cpu_hours = "CPU-hours", wall_clock_h = "wall-clock (h)"), name = NULL) +
      labs(title = "End-to-end pipeline cost vs image size",
           subtitle = "Total compute (CPU-hours) and wall-clock per slide.",
           x = "Image size (px)", y = "Time (hours)")
    save_fig(p13, "13_end_to_end_cost")
  }
}

# ── 14. BOTTLENECK STAGE by image size — which stage dominates wall-clock where ──
if (!is.null(cost) && all(c("bottleneck_stage", "target_px") %in% names(cost))) {
  bn <- cost %>% filter(varied_axis %in% size_axes, !is.na(bottleneck_stage))
  if (nrow(bn) > 0) {
    p14 <- bn %>% count(target_px, bottleneck_stage) %>%
      ggplot(aes(factor(target_px), n, fill = bottleneck_stage)) +
      geom_col(position = "fill") +
      scale_fill_manual(values = oi, name = "bottleneck") +
      scale_y_continuous(labels = percent_format()) +
      labs(title = "Pipeline bottleneck by image size",
           subtitle = "Share of runs whose slowest single process is each stage — the bottleneck shifts with size.",
           x = "Image size (px)", y = "Share of runs (unitless, 0-1)")
    save_fig(p14, "14_bottleneck_by_size")
  }
}

message("Built ", length(bench_figs), " figure(s) from ", normalizePath(adir))
