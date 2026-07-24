# =============================================================================
# registration_accuracy_plots.R — landmark-free, in-pipeline registration
# accuracy figures for analysis/registration_accuracy.Rmd. Twin of
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
      labs(title = "Nucleus-overlap Dice by registration stage",
           subtitle = "Independent check via DAPI segmentation overlap (not VALIS features); higher = better.",
           x = NULL, y = "matched-nucleus Dice", caption = REG_CAPTION)

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
        scale_colour_manual(values = oi[c(1, 2)], name = NULL) +
        labs(title = "Centroid residual displacement by stage",
             subtitle = "Matched-nucleus centroid distance in physical units; lower = tighter alignment.",
             x = NULL, y = "displacement (µm)", caption = REG_CAPTION)
    }
  }

  figs
}

# Module-level list for the Rmd (harmless on an empty data dir: returns list()).
reg_figs <- build_reg_figs()
message("registration_accuracy_plots.R: built ", length(reg_figs), " figure(s)")
