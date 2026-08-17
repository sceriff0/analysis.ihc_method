# =============================================================================
# registration_accuracy_plots.R — landmark-free, in-pipeline registration
# accuracy figures for analysis/benchmark_registration.Rmd. Twin of
# benchmark_plots.R: reads mirage paper_data CSVs from a directory (default
# data/benchmark/) and returns a named list of house-styled ggplots, skipping
# any figure whose CSV/columns are absent. Drop the sweep outputs in and re-knit.
# =============================================================================
.need <- c("ggplot2", "dplyr", "readr", "tidyr", "stringr", "tibble")
.missing <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing))
  stop("Missing R packages: ", paste(.missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(lapply(.need, library, character.only = TRUE))

source(here::here("code", "plot_theme.R"))   # house theme + oi/oi_ext palettes

REG_CAPTION          <- "Mirage registration QC · landmark-free · per moving slide"
STAGE_LEVELS_RTRE    <- c("original", "rigid", "non_rigid")            # VALIS rTRE stages
STAGE_LEVELS_OVERLAP <- c("native", "rigid", "non_rigid", "micro")     # warp_seg_qc stages

# Return NULL for a missing/empty CSV so a figure guarded on it is skipped.
.reg_read_opt <- function(dir, name) {
  p <- file.path(dir, name)
  if (!file.exists(p)) return(NULL)
  d <- suppressWarnings(readr::read_csv(p, show_col_types = FALSE))
  if (nrow(d) == 0) NULL else d
}

# Default-on-empty. Same definition in aggregation_compare.R, cell_tables.R and
# registration_accuracy_plots.R: these files are sourced in different orders by
# different reports, so a variant that only tested is.null() would silently take
# over and let a zero-length result through where a default was intended. (base R's
# %||% is the is.null-only form, which this deliberately widens.)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

build_reg_figs <- function(dir = here::here("data", "benchmark")) {
  figs <- list()

  # -- §1 VALIS self-reported feature rTRE, per moving slide, across stages -----
  # mirage auto-emits registration_valis_rtre.csv (make_tables.py). Columns are VALIS's
  # verbatim + run_id/summary_csv, so detect the id column and the stage columns rather than
  # hard-coding them: prefer the relative rTRE columns, fall back to the raw distance (_D) ones.
  vs <- .reg_read_opt(dir, "registration_valis_rtre.csv")
  if (!is.null(vs)) {
    id_col <- intersect(c("img_name", "name", "filename", "summary_csv"), names(vs))[1]
    rtre_cols <- grep("_rTRE$", names(vs), value = TRUE)
    metric_lab <- "VALIS rTRE (relative)"
    if (!length(rtre_cols)) {
      rtre_cols <- grep("_D$", names(vs), value = TRUE)
      metric_lab <- "VALIS matched-feature distance (D)"
    }
    rtre_cols <- rtre_cols[sub("_(rTRE|D)$", "", rtre_cols) %in% STAGE_LEVELS_RTRE]
    if (!is.na(id_col) && length(rtre_cols) >= 2) {
      long <- vs %>%
        dplyr::select(dplyr::all_of(c(id_col, rtre_cols))) %>%
        tidyr::pivot_longer(dplyr::all_of(rtre_cols),
                            names_to = "stage", values_to = "rTRE") %>%
        dplyr::mutate(stage = factor(sub("_(rTRE|D)$", "", stage), levels = STAGE_LEVELS_RTRE)) %>%
        dplyr::filter(is.finite(rTRE))
      if (nrow(long)) {
        n_slide <- dplyr::n_distinct(long[[id_col]])
        figs[["01_valis_rtre_by_stage"]] <-
          ggplot(long, aes(stage, rTRE, group = .data[[id_col]], colour = .data[[id_col]])) +
          geom_line(alpha = .8) + geom_point() +
          scale_colour_manual(values = rep_len(oi_ext, n_slide), guide = "none") +
          labs(title = "VALIS registration error by stage",
               subtitle = "Self-reported feature error per moving slide; lower = better.",
               x = NULL, y = metric_lab, caption = REG_CAPTION)
      }
    }
    if ("n_matches" %in% names(vs) && any(is.finite(vs$n_matches)) && !is.na(id_col)) {
      figs[["01b_valis_n_matches"]] <-
        ggplot(vs, aes(stats::reorder(.data[[id_col]], n_matches), n_matches)) +
        geom_col(fill = oi[1], width = .7) + coord_flip() +
        labs(title = "Feature matches behind each rTRE estimate",
             subtitle = "Correspondences VALIS used per slide; few matches = low-confidence estimate.",
             x = NULL, y = "feature matches (n)", caption = REG_CAPTION)
    }
  }

  # -- §2 independent overlap accuracy (DAPI-nucleus Dice + centroid residual) --
  ra <- .reg_read_opt(dir, "registration_accuracy.csv")
  if (!is.null(ra) && all(c("stage", "dice_matched") %in% names(ra))) {
    ra <- ra %>%
      dplyr::mutate(stage = factor(stage, levels = STAGE_LEVELS_OVERLAP)) %>%
      dplyr::filter(!is.na(stage))
    has_slide <- "moving" %in% names(ra)
    dice_p <- ggplot(ra, aes(stage, dice_matched)) +
      geom_boxplot(outlier.shape = NA, width = .5)
    if (has_slide)
      dice_p <- dice_p + geom_line(aes(group = moving), alpha = .25)
    figs[["02_overlap_dice_by_stage"]] <- dice_p +
      geom_jitter(width = .10, alpha = .5) +
      scale_x_discrete(labels = label_n(ra$stage)) +
      labs(title = "Nucleus-overlap Dice by registration stage",
           subtitle = "Independent check via DAPI segmentation overlap (not VALIS features); higher = better.",
           x = NULL, y = "Matched-nucleus Dice (unitless, 0-1)", caption = REG_CAPTION)

    if (all(c("displacement_um_p50", "displacement_um_p90") %in% names(ra))) {
      disp_long <- ra %>%
        tidyr::pivot_longer(c(displacement_um_p50, displacement_um_p90),
                            names_to = "pct", values_to = "um") %>%
        dplyr::mutate(pct = dplyr::recode(pct,
          displacement_um_p50 = "median", displacement_um_p90 = "90th pct")) %>%
        dplyr::filter(is.finite(um))
      figs[["02b_displacement_um_by_stage"]] <-
        ggplot(disp_long, aes(stage, um, colour = pct)) +
        geom_boxplot(outlier.shape = NA, width = .5, position = position_dodge(.6)) +
        # Count one percentile only. disp_long is pivoted long over median/90th,
        # so counting every row would report twice the number of runs.
        scale_x_discrete(labels = label_n(disp_long$stage[disp_long$pct == "median"])) +
        scale_colour_manual(values = oi[c(1, 2)], name = NULL) +
        labs(title = "Centroid residual displacement by stage",
             subtitle = "Matched-nucleus centroid distance in physical units; lower = tighter alignment.",
             x = NULL, y = "displacement (µm)", caption = REG_CAPTION)
    }
  }

  # -- §3 feature-distance improvement (OPTIONAL / LEGACY) ----------------------
  # mirage no longer emits feature_dist/*.json by default (the sweep uses reg_qc=2 + VALIS rTRE);
  # this renders only if a run set enable_feature_error. Needs jsonlite.
  fd_dir <- file.path(dir, "feature_dist")
  if (dir.exists(fd_dir) && requireNamespace("jsonlite", quietly = TRUE)) {
    jf <- list.files(fd_dir, pattern = "\\.json$", full.names = TRUE)
    if (length(jf)) {
      rows <- lapply(jf, function(j) {
        x <- tryCatch(jsonlite::fromJSON(j), error = function(e) NULL)
        if (is.null(x) || is.null(x$improvement$distance_reduction_percent)) return(NULL)
        data.frame(moving = x$moving_image %||% basename(j),
                   reduction_pct = as.numeric(x$improvement$distance_reduction_percent))
      })
      fd <- do.call(rbind, rows)
      if (!is.null(fd) && nrow(fd)) {
        figs[["03_feature_distance_reduction"]] <-
          ggplot(fd, aes(stats::reorder(moving, reduction_pct), reduction_pct)) +
          geom_col(fill = oi[3], width = .7) + coord_flip() +
          labs(title = "Feature-distance reduction after registration (legacy)",
               subtitle = "Per moving slide: percent drop in mean matched-feature distance (before → after).",
               x = NULL, y = "distance reduction (%)", caption = REG_CAPTION)
      }
    }
  }

  # -- §4 accuracy vs cost (Pareto) --------------------------------------------
  pm <- .reg_read_opt(dir, "param_matrix.csv")
  if (!is.null(pm) && all(c("reg_displacement_um_p50", "cpu_hours") %in% names(pm))) {
    pmf <- pm %>% dplyr::filter(is.finite(reg_displacement_um_p50), is.finite(cpu_hours))
    if (nrow(pmf)) {
      figs[["04_accuracy_vs_cost"]] <-
        ggplot(pmf, aes(cpu_hours, reg_displacement_um_p50)) +
        geom_point(size = 2, alpha = .8, colour = oi[1]) +
        labs(title = "Registration accuracy vs cost",
             subtitle = "Lower-left is better: less residual for fewer CPU-hours. One point per config.",
             x = "registration CPU-hours", y = "residual displacement, median (µm)",
             caption = REG_CAPTION)
    }
  }

  # Human labels for the two interchangeable VALIS error columns. They are NOT the
  # same quantity: rTRE is normalised by the image diagonal (unitless), D is a raw
  # pixel distance, so a figure that prints one label for both would misstate the axis.
  VALIS_ERR_LABS <- c(
    valis_non_rigid_rTRE = "VALIS relative TRE (unitless, fraction of image diagonal)",
    valis_non_rigid_D    = "VALIS mean feature distance (px)")

  # -- §5 agreement of the two independent estimates ---------------------------
  # The paper's thesis: VALIS's own feature error and the segmentation-overlap Dice are computed by
  # DIFFERENT methods yet should track per run. Both are pre-joined in param_matrix.csv, so this is a
  # single scatter — no extra plumbing. Prefer the relative rTRE median, fall back to the distance.
  if (!is.null(pm) && "reg_dice_matched" %in% names(pm)) {
    valis_col <- intersect(c("valis_non_rigid_rTRE", "valis_non_rigid_D"), names(pm))[1]
    if (!is.na(valis_col)) {
      ag <- pm %>% dplyr::filter(is.finite(.data[[valis_col]]), is.finite(reg_dice_matched))
      if (nrow(ag) > 1) {
        figs[["05_valis_vs_overlap_agreement"]] <-
          ggplot(ag, aes(.data[[valis_col]], reg_dice_matched)) +
          geom_point(size = 3, alpha = .8, colour = oi[1]) +
          # x used to be `valis_col` itself, i.e. the literal column name printed on
          # the axis. Which of the two columns was picked also changes the UNIT — rTRE
          # is a fraction of the image diagonal, D is a pixel distance — so the label
          # has to be looked up, not reused.
          labs(title = "Registration accuracy: VALIS vs segmentation-overlap",
               subtitle = paste("Independent estimates per run — VALIS feature error (x) vs",
                                "matched-nucleus Dice (y). They should track.",
                                n_note(nrow(ag), "runs")),
               x = VALIS_ERR_LABS[[valis_col]],
               y = "Matched-nucleus Dice (unitless, 0-1)", caption = REG_CAPTION)
      }
    }
  }

  # -- §6 THE PAPER'S ARM COMPARISON (Fig 4b/4c) -------------------------------
  # §1 and §2 answer "does registration help?" by walking the STAGES of one run.
  # This answers the different question the paper actually asks: "which CONFIG do we
  # ship?" — one point per arm, arms being the sweep's registration configurations
  # (accuracy preset x micro-registration on/off).
  #
  # THE TWO BASELINES COST NOTHING EXTRA. "No registration" and "rigid only" are not
  # separate runs to be launched: they are the `original_*` and `rigid_*` COLUMNS
  # every run already reports, i.e. that run's own before-and-after. They are folded
  # in here as two pseudo-arms so the baseline and the arms are read off one axis.
  # They are drawn from the same runs as the arms, so they are paired with them, not
  # independent of them — which is why they are labelled and ordered apart.
  arm_tbl <- .reg_arm_table(pm)
  if (!is.null(arm_tbl)) {
    base_tre <- .reg_stage_baselines(vs, c(original = "no registration", rigid = "rigid only"))

    if ("tre" %in% names(arm_tbl) && any(is.finite(arm_tbl$tre))) {
      d <- dplyr::bind_rows(
        dplyr::transmute(dplyr::filter(arm_tbl, is.finite(tre)),
                         arm, value = tre, kind = "configuration"),
        if (!is.null(base_tre)) dplyr::transmute(base_tre, arm = label, value, kind = "baseline"))
      d$arm <- stats::reorder(factor(d$arm), d$value)
      figs[["06_tre_by_arm"]] <-
        ggplot(d, aes(arm, value, colour = kind)) +
        geom_boxplot(outlier.shape = NA, width = .55, colour = "grey35") +
        geom_jitter(width = .12, height = 0, alpha = .75, size = 2) +
        scale_colour_manual(values = c(configuration = oi[1], baseline = oi[2]), name = NULL) +
        # sep = " " because coord_flip() puts these ticks beside a horizontal box,
        # where the default newline would double the left margin for no gain.
        scale_x_discrete(labels = label_n(d$arm, sep = " ")) +
        coord_flip() +
        labs(title = "Registration error by configuration",
             subtitle = paste("VALIS target registration error per arm, with the",
                              "no-registration and rigid-only baselines. Lower = better."),
             x = NULL, y = attr(arm_tbl, "tre_label") %||% "VALIS TRE", caption = REG_CAPTION)
    }

    if ("dice" %in% names(arm_tbl) && any(is.finite(arm_tbl$dice))) {
      base_dice <- .reg_long_baselines(ra, "dice_matched",
                                       c(native = "no registration", rigid = "rigid only"))
      d <- dplyr::bind_rows(
        dplyr::transmute(dplyr::filter(arm_tbl, is.finite(dice)),
                         arm, value = dice, kind = "configuration"),
        if (!is.null(base_dice)) dplyr::transmute(base_dice, arm = label, value, kind = "baseline"))
      d$arm <- stats::reorder(factor(d$arm), d$value)
      figs[["07_dice_by_arm"]] <-
        ggplot(d, aes(arm, value, colour = kind)) +
        geom_boxplot(outlier.shape = NA, width = .55, colour = "grey35") +
        geom_jitter(width = .12, height = 0, alpha = .75, size = 2) +
        scale_colour_manual(values = c(configuration = oi[1], baseline = oi[2]), name = NULL) +
        # sep = " " because coord_flip() puts these ticks beside a horizontal box,
        # where the default newline would double the left margin for no gain.
        scale_x_discrete(labels = label_n(d$arm, sep = " ")) +
        coord_flip() +
        labs(title = "DAPI-nucleus overlap Dice by configuration",
             subtitle = paste("Segmentation-overlap Dice per arm, same arms as the TRE",
                              "figure. Independent of VALIS's features. Higher = better."),
             x = NULL, y = "Matched-nucleus Dice (unitless, 0-1)", caption = REG_CAPTION)
    }
  }

  figs
}

# --- §6 helpers ---------------------------------------------------------------
# What makes two runs DIFFERENT arms. The sweep names its knobs differently across
# mirage versions, so the arm label is assembled from whichever of these columns
# param_matrix.csv actually carries rather than from a hard-coded pair — a renamed
# knob then costs a generic label, not a missing figure.
REG_ARM_PATTERNS <- c(preset = "preset|accuracy|max_?dim|max_image_dim",
                      micro  = "micro")

# One row per arm: its label plus the two headline metrics. Returns NULL when
# param_matrix.csv is absent or carries neither metric, so §6 skips as a unit.
# `tre_label` rides along as an attribute because which TRE column exists (relative
# rTRE vs raw distance D) determines what the axis is allowed to claim.
.reg_arm_table <- function(pm) {
  if (is.null(pm)) return(NULL)
  tre_col  <- intersect(c("valis_non_rigid_rTRE", "valis_non_rigid_D"), names(pm))[1]
  dice_col <- intersect("reg_dice_matched", names(pm))[1]
  if (is.na(tre_col) && is.na(dice_col)) return(NULL)

  knobs <- unlist(lapply(REG_ARM_PATTERNS, function(pat)
    grep(pat, names(pm), value = TRUE, ignore.case = TRUE)[1]))
  knobs <- unique(knobs[!is.na(knobs)])

  if (!length(knobs)) {
    # No knob column survived. Every run is then its own arm, labelled by run id —
    # honest, and still the right SHAPE of figure, rather than pretending one arm.
    id <- intersect(c("run_id", "run", "config", "name"), names(pm))[1]
    warning("registration arms: no preset/micro column in param_matrix.csv — ",
            "labelling each run separately. Arms will not be grouped.")
    arm <- if (is.na(id)) paste("run", seq_len(nrow(pm))) else as.character(pm[[id]])
  } else {
    arm <- do.call(paste, c(lapply(knobs, function(k)
      paste0(sub("^reg_", "", k), "=", pm[[k]])), list(sep = " · ")))
  }

  out <- tibble::tibble(
    arm  = arm,
    tre  = if (is.na(tre_col))  NA_real_ else suppressWarnings(as.numeric(pm[[tre_col]])),
    dice = if (is.na(dice_col)) NA_real_ else suppressWarnings(as.numeric(pm[[dice_col]])))
  attr(out, "tre_label") <- if (is.na(tre_col)) NULL
    else if (grepl("rTRE$", tre_col)) "VALIS rTRE (relative)" else "VALIS matched-feature distance (D)"
  out
}

# Turn named stage columns into labelled baseline rows. `stages` maps the column
# stem to the label the figure should show. Returns NULL when none are present, so
# a figure without baselines still draws its arms rather than failing.
.reg_stage_baselines <- function(df, stages) {
  if (is.null(df)) return(NULL)
  rows <- lapply(names(stages), function(st) {
    col <- grep(paste0("^", st, "_(rTRE|D)$"), names(df), value = TRUE)[1]
    if (is.na(col)) col <- if (st %in% names(df)) st else NA_character_
    if (is.na(col)) return(NULL)
    v <- suppressWarnings(as.numeric(df[[col]]))
    v <- v[is.finite(v)]
    if (!length(v)) return(NULL)
    tibble::tibble(label = unname(stages[[st]]), value = v)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}

# The Dice table is LONG (one row per run x moving slide x stage) where the VALIS
# table is WIDE (one column per stage), so its baselines are rows to filter, not
# columns to find. Same contract as .reg_stage_baselines(): labelled values or NULL.
.reg_long_baselines <- function(df, value_col, stages) {
  if (is.null(df) || !all(c("stage", value_col) %in% names(df))) return(NULL)
  rows <- lapply(names(stages), function(st) {
    v <- suppressWarnings(as.numeric(df[[value_col]][as.character(df$stage) == st]))
    v <- v[is.finite(v)]
    if (!length(v)) return(NULL)
    tibble::tibble(label = unname(stages[[st]]), value = v)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}

# Module-level list for the Rmd (harmless on an empty data dir: returns list()).
reg_figs <- build_reg_figs()
message("registration_accuracy_plots.R: built ", length(reg_figs), " figure(s)")
