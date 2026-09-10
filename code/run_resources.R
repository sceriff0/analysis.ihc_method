# =============================================================================
# run_resources.R — readers + figures for ONE real mirage run's resource profile
# =============================================================================
# Reads the four tables `benchmarks/pull_run_resources.sh` (mirage, `benchmarking`
# branch) copies into data/run_resources/:
#
#   run_resources_tasks.csv       one row per Nextflow task
#   run_resources_processes.csv   one row per process
#   run_resources_fits.csv        one row per process x target: target ~ input_gb
#   run_resources_summary.csv     one row: the whole run
#   run_resources.dict.md         every column, defined
#
# They come from Nextflow's trace.txt (realtime, %cpu, cpus, memory requested,
# peak_rss) joined with mirage's size_logs/input_sizes.csv on (process, tag). Column
# meanings live in the .dict.md next to the data and are not repeated here.
#
# Same contract as run_qc.R: every reader returns an EMPTY tibble for a missing file,
# and build_run_resources_figs() returns a named list of house-styled ggplots with
# any figure whose input is absent simply missing -- the page renders what it gets.

.need <- c("dplyr", "tidyr", "readr", "ggplot2", "forcats", "here", "tibble")
.missing <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing))
  stop("run_resources.R needs: ", paste(.missing, collapse = ", "), call. = FALSE)

RUN_RESOURCES_ROOT <- here::here("data", "run_resources")
RUN_RESOURCES_CAPTION <- paste(
  "Nextflow trace.txt joined with mirage size_logs/input_sizes.csv;",
  "one run, real slides. benchmarks/analysis/run_resources.py")

.read_rr <- function(name, root) {
  path <- file.path(root, paste0(name, ".csv"))
  if (!file.exists(path)) return(tibble::tibble())
  tryCatch(readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
           error = function(e) tibble::tibble())
}

run_resources_tables <- function(root = RUN_RESOURCES_ROOT) {
  list(tasks     = .read_rr("run_resources_tasks", root),
       processes = .read_rr("run_resources_processes", root),
       fits      = .read_rr("run_resources_fits", root),
       summary   = .read_rr("run_resources_summary", root))
}

# Short process labels for axes: the trace already strips the workflow prefix,
# so this only shortens the two very long names.
.short_proc <- function(x) {
  x <- gsub("GENERATE_", "GEN_", x)
  gsub("EXTRACT_", "EXTR_", x)
}

# =============================================================================
# Figures
# =============================================================================
build_run_resources_figs <- function(root = RUN_RESOURCES_ROOT,
                                     tables = run_resources_tables(root)) {
  figs <- list()
  lab  <- function(...) ggplot2::labs(..., caption = RUN_RESOURCES_CAPTION)
  pr   <- tables$processes
  tk   <- tables$tasks
  fits <- tables$fits

  # -- §1 where the wall-time went: per-process task hours, longest task marked -----
  if (nrow(pr)) {
    d <- pr |>
      dplyr::mutate(process = forcats::fct_reorder(.short_proc(process), realtime_total_h))
    figs[["01_wall_by_process"]] <-
      ggplot2::ggplot(d, ggplot2::aes(realtime_total_h, process)) +
      ggplot2::geom_col(fill = "grey70") +
      ggplot2::geom_point(ggplot2::aes(x = realtime_max_h), shape = 124, size = 4) +
      lab(x = "task wall-time, hours (bar = sum over tasks; tick = longest single task)",
          y = NULL,
          title = "Where the run's time went",
          subtitle = "The tick is what a process costs on the critical path when its tasks run in parallel")
  }

  # -- §2 CPU-hours: used vs reserved ---------------------------------------------
  if (nrow(pr) && any(!is.na(pr$cpu_h_reserved_total))) {
    d <- pr |>
      dplyr::select(process, used = cpu_h_used_total, reserved = cpu_h_reserved_total) |>
      tidyr::pivot_longer(-process, names_to = "kind", values_to = "cpu_h") |>
      dplyr::mutate(process = forcats::fct_reorder(.short_proc(process), cpu_h, .fun = max),
                    kind = factor(kind, levels = c("reserved", "used")))
    figs[["02_cpu_hours_used_vs_reserved"]] <-
      ggplot2::ggplot(d, ggplot2::aes(cpu_h, process, fill = kind)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.65) +
      ggplot2::scale_fill_manual(values = c(reserved = "grey80", used = "black"), name = NULL) +
      lab(x = "CPU-hours", y = NULL,
          title = "CPU reserved against CPU actually used",
          subtitle = "The gap is idle allocation: cores requested for a task that ran on fewer")
  }

  # -- §3 memory: peak RSS against the request, per task ---------------------------
  if (nrow(tk) && any(!is.na(tk$memory_req_gb) & !is.na(tk$peak_rss_gb))) {
    d <- tk |>
      dplyr::filter(!is.na(memory_req_gb), !is.na(peak_rss_gb)) |>
      dplyr::mutate(process = .short_proc(process), failed = as.logical(failed))
    lim <- max(c(d$memory_req_gb, d$peak_rss_gb), na.rm = TRUE)
    figs[["03_peak_rss_vs_request"]] <-
      ggplot2::ggplot(d, ggplot2::aes(memory_req_gb, peak_rss_gb, colour = process)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
      ggplot2::geom_point(ggplot2::aes(shape = failed), size = 2, alpha = 0.8) +
      ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                                  labels = c(`FALSE` = "completed", `TRUE` = "failed"), name = NULL) +
      ggplot2::coord_equal(xlim = c(0, lim), ylim = c(0, lim)) +
      lab(x = "memory requested, GiB", y = "peak RSS, GiB", colour = NULL,
          title = "Peak memory against the reservation, one point per task",
          subtitle = "On the dashed line a task used exactly what it asked for; below it the difference sat idle")
  }

  # -- §4 scaling inside this run: peak RSS and wall-time against input size --------
  if (nrow(tk) && any(!is.na(tk$input_gb))) {
    lv <- c("peak RSS, GiB", "wall-time, hours")
    d <- tk |>
      dplyr::filter(!is.na(input_gb)) |>
      dplyr::transmute(process = .short_proc(process), input_gb,
                       `peak RSS, GiB` = peak_rss_gb, `wall-time, hours` = realtime_s / 3600) |>
      tidyr::pivot_longer(c(`peak RSS, GiB`, `wall-time, hours`),
                          names_to = "target", values_to = "value") |>
      dplyr::filter(!is.na(value)) |>
      dplyr::mutate(target = factor(target, levels = lv))
    lines <- NULL
    if (nrow(fits)) {
      lines <- fits |>
        dplyr::filter(as.logical(fit_ok)) |>
        dplyr::mutate(process = .short_proc(process),
                      # realtime fits are in seconds; the panel plots hours
                      slope = ifelse(target == "realtime_s", slope / 3600, slope),
                      intercept = ifelse(target == "realtime_s", intercept / 3600, intercept),
                      target = factor(ifelse(target == "realtime_s", lv[2], lv[1]), levels = lv))
    }
    p <- ggplot2::ggplot(d, ggplot2::aes(input_gb, value, colour = process)) +
      ggplot2::geom_point(size = 2, alpha = 0.8) +
      ggplot2::facet_wrap(~ target, scales = "free_y") +
      lab(x = "task input, GiB (size_logs)", y = NULL, colour = NULL,
          title = "Cost against input size, inside this run",
          subtitle = paste("A line is drawn only where a process ran on >= 3 inputs of different size;",
                           "a single point per process is a point, not a law"))
    if (!is.null(lines) && nrow(lines))
      p <- p + ggplot2::geom_abline(data = lines,
                                    ggplot2::aes(slope = slope, intercept = intercept, colour = process),
                                    linewidth = 0.4)
    figs[["04_cost_vs_input"]] <- p
  }

  figs
}

# The per-process table the page prints above the figures: the columns a reader
# wants first, rounded, with the RSS-per-input slope beside them when there is one.
run_resources_process_table <- function(tables) {
  pr <- tables$processes
  if (!nrow(pr)) return(tibble::tibble())
  rss_fit <- if (nrow(tables$fits))
    tables$fits |>
      dplyr::filter(target == "peak_rss_gb") |>
      dplyr::transmute(process,
                       `RSS per GiB input` = ifelse(as.logical(fit_ok), round(slope, 2), NA_real_),
                       `fit n` = n)
  else tibble::tibble(process = character(), `RSS per GiB input` = numeric(), `fit n` = integer())
  pr |>
    dplyr::transmute(process,
                     tasks = n_tasks, failed = n_failed,
                     `wall h (sum)` = round(realtime_total_h, 2),
                     `longest task h` = round(realtime_max_h, 2),
                     `share of task time` = round(wall_share, 3),
                     `CPU-h used` = round(cpu_h_used_total, 2),
                     `CPU-h reserved` = round(cpu_h_reserved_total, 2),
                     `CPU efficiency` = round(cpu_efficiency, 2),
                     `peak RSS GiB` = round(peak_rss_max_gb, 1),
                     `mem requested GiB` = round(memory_req_max_gb, 1),
                     `RSS / request (max)` = round(rss_utilisation_max, 2),
                     `input GiB (mean)` = round(input_gb_mean, 2)) |>
    dplyr::left_join(rss_fit, by = "process") |>
    dplyr::arrange(dplyr::desc(`wall h (sum)`))
}
