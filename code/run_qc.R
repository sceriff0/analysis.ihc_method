# =============================================================================
# run_qc.R  —  mirage's QC of the run on THE STUDY SAMPLES
#
# The benchmark pages (benchmarks.Rmd, registration_accuracy.Rmd) read the tables
# mirage's `benchmarks/` sweep produces on SYNTHETIC images: rescaled copies with a
# known injected offset, run to measure cost and scaling. This file reads the QC the
# pipeline emits for a NORMAL run, on the real patient slides — the numbers that say
# whether *this cohort* was registered and phenotyped well enough to analyse.
#
# Everything lives under the same data/mirage/<patient>/ tree the cell loader reads,
# so a run already copied (or symlinked) in for clinical_data_mirage needs nothing
# further:
#
#   qc/registration/*_seg_qc.json    staged segmentation-overlap QC (reg_qc = 2).
#                                    Per (moving slide, stage): matched-nucleus Dice,
#                                    per-pair IoU, centroid residual in px and um,
#                                    and the delta against the rigid anchor.
#   registered/summary/*.csv         VALIS's OWN error, written during register().
#                                    Columns are whatever VALIS declared (rTRE / D /
#                                    n_matches variants) — read generically.
#   qc/registration/*_tre.json       STARE's own TRE, when the tiled path ran.
#                                    coarse/rigid/post-refinement percentiles PLUS a
#                                    per-tile spatial map VALIS cannot give.
#   phenotyping/phenotype_qc.json    conformal calibration: the alpha actually used,
#                                    whether CRC ran, markers that degraded.
#   phenotyping/constraint_audit.csv observed vs nominal rate per panel constraint.
#
# THREE INDEPENDENT REGISTRATION-ACCURACY SIGNALS, and that is the point
#   VALIS rTRE and STARE TRE are *intrinsic*: each method grading itself from the
#   features it registered on. They cannot be compared to each other directly — only
#   one of them ran — but each can be compared to the segmentation-overlap Dice,
#   which is computed by a different method entirely (DAPI nuclei, not features) and
#   so is the independent check on whichever intrinsic number you have. An intrinsic
#   score can look excellent while the overlap says otherwise: that is the failure
#   mode of grading yourself on your own correspondences.
#
# ONE TRAP IN THE VALIS SUMMARIES
#   register_micro() re-runs measure_error() and OVERWRITES <name>_summary.csv, so
#   that file is POST-micro. <name>_summary_premicro.csv (written only at
#   micro_reg >= 2) is the one that isolates the non-rigid stage. They are tagged
#   `stage_scope` here rather than pooled — averaging them would mix a pre- and a
#   post-micro estimate of the same slide.
#
# Every reader returns an empty tibble for a missing/unparseable artifact rather
# than erroring, so a partial run still reports what it has.
# =============================================================================
.need <- c("dplyr", "tidyr", "readr", "purrr", "tibble", "ggplot2", "jsonlite", "fs")
.missing <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing))
  stop("run_qc.R needs: ", paste(.missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(lapply(.need, library, character.only = TRUE))

# slide_key() normalises the patient id read off each directory name. plot_theme.R
# comes with it, so the figures below are house-styled with nothing further to call.
# This file does NOT depend on mirage_cells.R: the QC page reads the run's QC
# artifacts, not its cells, and must render even when no cell table is present.
source(here::here("code", "validation_helpers.R"))

RUN_QC_ROOT <- here::here("data", "mirage")   # the same tree the cell loader reads

RUN_QC_CAPTION      <- "mirage run QC · study samples · landmark-free"
QC_STAGE_LEVELS     <- c("native", "rigid", "non_rigid", "micro")   # warp_seg_qc stages
VALIS_STAGE_LEVELS  <- c("original", "rigid", "non_rigid")          # VALIS rTRE stages

# --- helpers ----------------------------------------------------------------
# Patient directories under `root`, keyed by the directory name. A mirage outdir
# also carries cohort-level dirs (qc/, phenotyping/, size_logs/), so a patient is one
# with a qc/ or phenotyping/ subdirectory of its own.
.qc_patient_dirs <- function(root) {
  if (!fs::dir_exists(root)) return(character(0))
  dirs <- as.character(fs::dir_ls(root, type = "directory"))
  dirs[fs::dir_exists(file.path(dirs, "qc")) |
       fs::file_exists(file.path(dirs, "phenotyping", "phenotype_qc.json"))]
}

.read_json <- function(path) tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE),
                                      error = function(e) NULL)

# Every file matching `glob` under one patient's subdirectory.
.qc_files <- function(dir, sub, glob) {
  d <- file.path(dir, sub)
  if (!fs::dir_exists(d)) return(character(0))
  as.character(fs::dir_ls(d, glob = glob))
}

# `%||%` is defined in cell_tables.R / plot_theme.R's siblings; keep run_qc.R
# self-contained so it can be sourced alone.
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# --- 1. staged segmentation-overlap QC (reg_qc = 2) --------------------------
# One row per (patient, moving slide, stage). Mirrors mirage's own
# benchmarks/analysis/lib/quality.py harvest_registration_qc(), which reads the same
# JSONs out of a sweep — the schema is shared, only the tree differs.
# The column set read_seg_qc() promises, so an empty run still satisfies every
# `nrow()` and column check downstream.
SEG_QC_EMPTY <- tibble::tibble(
  patient_id = character(), moving = character(),
  stage = factor(levels = QC_STAGE_LEVELS),
  n_pairs = numeric(), pair_fraction = numeric(), iou_mean = numeric(),
  iou_p50 = numeric(), dice_matched = numeric(), disp_um_p50 = numeric(),
  disp_um_p90 = numeric(), disp_px_p50 = numeric(),
  d_dice_vs_rigid = numeric(), d_disp_um_vs_rigid = numeric())

read_seg_qc <- function(root = RUN_QC_ROOT) {
  out <- purrr::map_dfr(.qc_patient_dirs(root), function(dir) {
    purrr::map_dfr(.qc_files(dir, "qc/registration", "*_seg_qc.json"), function(f) {
      d <- .read_json(f)
      if (is.null(d) || is.null(d$stages)) return(tibble::tibble())
      stages <- d$stage_order %||% names(d$stages)
      purrr::map_dfr(stages, function(st) {
        s  <- d$stages[[st]]
        dv <- (d$delta_vs_anchor %||% list())[[st]] %||% list()
        tibble::tibble(
          patient_id       = slide_key(d$patient_id %||% fs::path_file(dir)),
          moving           = d$moving %||% fs::path_ext_remove(fs::path_file(f)),
          stage            = st,
          n_pairs          = as.numeric(s$n_pairs %||% NA),
          pair_fraction    = as.numeric((d$matching %||% list())$pair_fraction %||% NA),
          iou_mean         = as.numeric(s$iou_mean %||% NA),
          iou_p50          = as.numeric(s$iou_p50 %||% NA),
          dice_matched     = as.numeric(s$dice_matched %||% NA),
          disp_um_p50      = as.numeric(s$displacement_um_p50 %||% NA),
          disp_um_p90      = as.numeric(s$displacement_um_p90 %||% NA),
          disp_px_p50      = as.numeric(s$displacement_px_p50 %||% NA),
          d_dice_vs_rigid  = as.numeric(dv$dice_matched %||% NA),
          d_disp_um_vs_rigid = as.numeric(dv$displacement_um_p50 %||% NA))
      })
    })
  })
  # An empty result has no columns at all, and `stage` would then resolve to base R's
  # stats::stage() rather than to a column — the error a reader hits on their first
  # knit, before any QC has been copied in. Return the typed empty frame instead.
  if (nrow(out) == 0) return(SEG_QC_EMPTY)
  out |>
    dplyr::mutate(stage = factor(stage, levels = QC_STAGE_LEVELS)) |>
    dplyr::filter(!is.na(stage))
}

# --- 2. VALIS's own error ----------------------------------------------------
# Column-agnostic: VALIS's schema varies by build, so keep what it wrote and tag the
# provenance. `stage_scope` separates the post-micro summary from the pre-micro one.
read_valis_summary <- function(root = RUN_QC_ROOT) {
  purrr::map_dfr(.qc_patient_dirs(root), function(dir) {
    purrr::map_dfr(.qc_files(dir, "registered/summary", "*.csv"), function(f) {
      d <- tryCatch(readr::read_csv(f, show_col_types = FALSE, progress = FALSE),
                    error = function(e) NULL)
      if (is.null(d) || nrow(d) == 0) return(tibble::tibble())
      nm <- fs::path_file(f)
      dplyr::mutate(d,
        patient_id  = slide_key(fs::path_file(dir)),
        summary_csv = nm,
        stage_scope = if (grepl("premicro", nm)) "pre-micro" else "post-micro",
        .before = 1)
    })
  })
}

# VALIS's per-stage error in long form: detect the stage columns rather than naming
# them, preferring the RELATIVE rTRE (a fraction of the image diagonal, comparable
# across slide sizes) over the raw distance D.
valis_error_long <- function(valis) {
  if (nrow(valis) == 0) return(tibble::tibble())
  cols <- grep("_rTRE$", names(valis), value = TRUE)
  metric <- "VALIS rTRE (relative)"
  if (!length(cols)) {
    cols <- grep("_D$", names(valis), value = TRUE)
    metric <- "VALIS matched-feature distance (D)"
  }
  cols <- cols[sub("_(rTRE|D)$", "", cols) %in% VALIS_STAGE_LEVELS]
  if (!length(cols)) return(tibble::tibble())
  id <- intersect(c("img_name", "name", "filename"), names(valis))[1]
  valis |>
    dplyr::mutate(slide = if (is.na(id)) summary_csv else .data[[id]]) |>
    dplyr::select(patient_id, slide, stage_scope, dplyr::all_of(cols)) |>
    tidyr::pivot_longer(dplyr::all_of(cols), names_to = "stage", values_to = "error") |>
    dplyr::mutate(stage  = factor(sub("_(rTRE|D)$", "", stage), levels = VALIS_STAGE_LEVELS),
                  metric = metric) |>
    dplyr::filter(is.finite(error))
}

# --- 3. STARE's own TRE (tiled path only) ------------------------------------
# Per patient/moving summary. Absent when the run used the VALIS path.
read_stare_tre <- function(root = RUN_QC_ROOT) {
  purrr::map_dfr(.qc_patient_dirs(root), function(dir) {
    purrr::map_dfr(.qc_files(dir, "qc/registration", "*_tre.json"), function(f) {
      d <- .read_json(f)
      if (is.null(d) || is.null(d$coarse_tre_px)) return(tibble::tibble())
      pct <- function(x, k) as.numeric((d[[x]] %||% list())[[k]] %||% NA)
      tibble::tibble(
        patient_id     = slide_key(fs::path_file(dir)),
        moving         = sub("_tre$", "", fs::path_ext_remove(fs::path_file(f))),
        coarse_tre_px  = as.numeric(d$coarse_tre_px),
        n_inliers      = as.numeric(d$n_inliers %||% NA),
        n_tiles        = as.numeric(d$n_tiles %||% NA),
        mesh_refined   = isTRUE(d$mesh_refined),
        rigid_p50      = pct("rigid_tre_px", "p50"),
        rigid_p90      = pct("rigid_tre_px", "p90"),
        rigid_max      = pct("rigid_tre_px", "max"),
        after_p50      = pct("residual_after_px", "p50"),
        after_p90      = pct("residual_after_px", "p90"))
    })
  })
}

# The per-tile records behind those percentiles: a SPATIAL map of registration error
# across the slide, which the feature-summary numbers average away. This is what the
# tiled path gives that VALIS does not.
read_stare_tiles <- function(root = RUN_QC_ROOT) {
  purrr::map_dfr(.qc_patient_dirs(root), function(dir) {
    purrr::map_dfr(.qc_files(dir, "qc/registration", "*_tre.json"), function(f) {
      d <- .read_json(f)
      if (is.null(d$tiles) || !NROW(d$tiles)) return(tibble::tibble())
      tibble::as_tibble(d$tiles) |>
        dplyr::mutate(patient_id = slide_key(fs::path_file(dir)),
                      moving = sub("_tre$", "", fs::path_ext_remove(fs::path_file(f))),
                      .before = 1)
    })
  })
}

# --- 4. phenotyping calibration ---------------------------------------------
# `chosen_alpha` below `alpha_target` means CRC tightened the risk budget; a marker
# in `degraded_markers` had too little calibration data and fell back to a reporting-
# only sign, so populations gated on it are the ones to distrust.
read_phenotype_qc <- function(root = RUN_QC_ROOT) {
  purrr::map_dfr(.qc_patient_dirs(root), function(dir) {
    f <- file.path(dir, "phenotyping", "phenotype_qc.json")
    if (!file.exists(f)) return(tibble::tibble())
    d <- .read_json(f)
    if (is.null(d)) return(tibble::tibble())
    tibble::tibble(
      patient_id     = slide_key(fs::path_file(dir)),
      n_cells        = as.numeric(d$n_cells %||% NA),
      chosen_alpha   = as.numeric(d$chosen_alpha %||% NA),
      alpha_target   = as.numeric(d$alpha_target %||% NA),
      crc_ran        = isTRUE(d$crc_ran),
      reporting_mode = isTRUE(d$reporting_mode),
      n_degraded     = length(d$degraded_markers %||% character(0)),
      degraded       = paste(d$degraded_markers %||% character(0), collapse = ", "),
      density_radius = as.numeric(d$density_radius %||% NA),
      n_bins         = as.numeric(d$n_bins %||% NA))
  })
}

# Observed vs nominal co-expression rate for every constraint in the panel. A
# `verdict` other than pass means the data contradicts a biological exclusivity the
# panel asserts — usually spillover or a mis-set gate, and a reason to distrust the
# phenotypes that depend on those markers.
read_constraint_audit <- function(root = RUN_QC_ROOT) {
  purrr::map_dfr(.qc_patient_dirs(root), function(dir) {
    f <- file.path(dir, "phenotyping", "constraint_audit.csv")
    if (!file.exists(f)) return(tibble::tibble())
    d <- tryCatch(readr::read_csv(f, show_col_types = FALSE, progress = FALSE),
                  error = function(e) NULL)
    if (is.null(d) || nrow(d) == 0) return(tibble::tibble())
    dplyr::mutate(d, patient_id = slide_key(fs::path_file(dir)), .before = 1)
  })
}

# Everything at once, for the report.
run_qc_tables <- function(root = RUN_QC_ROOT) {
  list(seg_qc     = read_seg_qc(root),
       valis      = read_valis_summary(root),
       stare      = read_stare_tre(root),
       stare_tiles= read_stare_tiles(root),
       pheno_qc   = read_phenotype_qc(root),
       audit      = read_constraint_audit(root))
}

# =============================================================================
# Figures
# =============================================================================
# Same contract as registration_accuracy_plots.R / benchmark_plots.R: a named list
# of house-styled ggplots, with any figure whose input is absent simply missing. The
# report renders whatever comes back, so a VALIS run and a tiled run each show their
# own intrinsic section without either page needing to know which ran.
build_run_qc_figs <- function(root = RUN_QC_ROOT, tables = run_qc_tables(root)) {
  figs <- list()
  lab  <- function(...) labs(..., caption = RUN_QC_CAPTION)

  # -- §1 VALIS intrinsic error, per slide across stages -----------------------
  vl <- valis_error_long(tables$valis)
  if (nrow(vl)) {
    figs[["01_valis_error_by_stage"]] <-
      ggplot(vl, aes(stage, error, group = interaction(patient_id, slide),
                     colour = patient_id)) +
      geom_line(alpha = .7) + geom_point(size = 1.8) +
      scale_colour_oi(name = "patient") +
      facet_wrap(~ stage_scope) +
      lab(title = "VALIS registration error by stage",
          subtitle = paste("VALIS grading itself from its own feature matches, one line",
                           "per moving slide. Lower = better."),
          x = NULL, y = vl$metric[1])
    if ("n_matches" %in% names(tables$valis)) {
      # Label each bar with the slide VALIS named, falling back to the summary file
      # plus a row index when the build wrote no name column.
      id <- intersect(c("img_name", "name", "filename"), names(tables$valis))[1]
      nm <- tables$valis |>
        dplyr::filter(is.finite(n_matches)) |>
        dplyr::mutate(slide = if (is.na(id)) paste(summary_csv, dplyr::row_number())
                              else as.character(.data[[id]]),
                      slide = paste(patient_id, slide))
      if (nrow(nm))
        figs[["01b_valis_n_matches"]] <-
          ggplot(nm, aes(stats::reorder(slide, n_matches), n_matches, fill = patient_id)) +
          geom_col(width = .7) + coord_flip() + scale_fill_oi(name = "patient") +
          lab(title = "Feature matches behind each VALIS estimate",
              subtitle = "Few correspondences means a low-confidence error estimate.",
              x = NULL, y = "feature matches (n)")
    }
  }

  # -- §2 STARE intrinsic TRE ---------------------------------------------------
  st <- tables$stare
  if (nrow(st)) {
    long <- st |>
      dplyr::select(patient_id, moving, coarse_tre_px, rigid_p50, rigid_p90, after_p50) |>
      tidyr::pivot_longer(-c(patient_id, moving), names_to = "metric", values_to = "px") |>
      dplyr::filter(is.finite(px)) |>
      dplyr::mutate(metric = factor(metric,
        levels = c("coarse_tre_px", "rigid_p50", "rigid_p90", "after_p50"),
        labels = c("coarse (global rigid)", "per-tile rigid, median",
                   "per-tile rigid, p90", "after mesh, median")))
    figs[["02_stare_tre"]] <-
      ggplot(long, aes(metric, px, colour = patient_id)) +
      geom_point(size = 2.6, alpha = .85,
                 position = position_jitter(width = .12, height = 0)) +
      scale_colour_oi(name = "patient") +
      lab(title = "STARE intrinsic TRE (tiled registration path)",
          subtitle = paste("The method's own target-registration error, from the features it",
                           "registered on. 'after mesh' is the post-refinement residual."),
          x = NULL, y = "TRE (pixels)") +
      theme(axis.text.x = element_text(angle = 20, hjust = 1))
  }

  # -- §2b STARE spatial TRE map — what a summary number averages away ----------
  ti <- tables$stare_tiles
  if (nrow(ti) && all(c("ix", "iy", "tre_rigid") %in% names(ti))) {
    figs[["02b_stare_tre_map"]] <-
      ggplot(ti, aes(ix, iy, fill = tre_rigid)) +
      geom_tile() + coord_equal() +
      scale_fill_seq(name = "TRE (px)") +
      facet_wrap(~ patient_id + moving) +
      guides(fill = guide_cbar()) +
      lab(title = "Where the registration error sits on the slide",
          subtitle = paste("Per-tile rigid-stage TRE. A hot corner or edge is a local failure",
                           "a single summary number hides — tissue folds and slide edges."),
          x = "tile column", y = "tile row") +
      theme_paper_tile()
  }

  # -- §3 the independent cross-check: DAPI-nucleus overlap ---------------------
  sq <- tables$seg_qc
  if (nrow(sq)) {
    figs[["03_overlap_dice_by_stage"]] <-
      ggplot(sq, aes(stage, dice_matched)) +
      geom_boxplot(outlier.shape = NA, width = .5) +
      geom_line(aes(group = interaction(patient_id, moving), colour = patient_id), alpha = .4) +
      geom_point(aes(colour = patient_id), alpha = .8) +
      scale_colour_oi(name = "patient") +
      lab(title = "Matched-nucleus Dice by registration stage",
          subtitle = paste("Computed from DAPI segmentation overlap, NOT from the features the",
                           "registration used — the independent check. Higher = better."),
          x = NULL, y = "matched-nucleus Dice")

    if (any(is.finite(sq$disp_um_p50)))
      figs[["03b_displacement_um_by_stage"]] <-
        sq |>
        tidyr::pivot_longer(c(disp_um_p50, disp_um_p90), names_to = "pct", values_to = "um") |>
        dplyr::mutate(pct = dplyr::recode(pct, disp_um_p50 = "median", disp_um_p90 = "90th pct")) |>
        dplyr::filter(is.finite(um)) |>
        ggplot(aes(stage, um, colour = pct)) +
        geom_boxplot(outlier.shape = NA, width = .5, position = position_dodge(.6)) +
        scale_colour_manual(values = oi[c(1, 2)], name = NULL) +
        lab(title = "Centroid residual displacement by stage",
            subtitle = "Physical units, so this is the number to quote as spatial resolution.",
            x = NULL, y = "displacement (µm)")

    # Pair fraction is the sanity gate on every number above it: a low fraction means
    # the Dice was computed over a thin, and probably biased, subset of cells.
    if (any(is.finite(sq$pair_fraction)))
      figs[["03c_pair_fraction"]] <-
        sq |>
        dplyr::distinct(patient_id, moving, pair_fraction) |>
        dplyr::filter(is.finite(pair_fraction)) |>
        ggplot(aes(stats::reorder(paste(patient_id, moving), pair_fraction),
                   pair_fraction, fill = patient_id)) +
        geom_col(width = .7) + coord_flip() +
        geom_hline(yintercept = .5, linetype = "dashed", colour = REF_LINE) +
        scale_fill_oi(name = "patient") +
        lab(title = "Fraction of cells that paired between slides",
            subtitle = paste("The evidence base for every accuracy number. Below the dashed 0.5",
                             "line the pairing is thin and the Dice is not trustworthy."),
            x = NULL, y = "pair fraction")
  }

  # -- §4 do the intrinsic and the independent estimate agree? ------------------
  # Per (patient, moving) at the final stage: whichever intrinsic number the run has,
  # against the overlap Dice. Disagreement is the interesting outcome — it means the
  # method's self-assessment is not measuring what a second method sees.
  if (nrow(sq)) {
    final <- sq |>
      dplyr::filter(!is.na(dice_matched)) |>
      dplyr::group_by(patient_id, moving) |>
      dplyr::slice_max(as.integer(stage), n = 1, with_ties = FALSE) |>
      dplyr::ungroup()
    intrinsic <- if (nrow(st))
      dplyr::transmute(st, patient_id, moving,
                       intrinsic = dplyr::coalesce(after_p50, rigid_p50),
                       what = "STARE TRE (px)")
    else if (nrow(vl))
      vl |>
        dplyr::filter(stage == dplyr::last(levels(droplevels(stage)))) |>
        dplyr::transmute(patient_id, moving = slide, intrinsic = error,
                         what = metric)
    else NULL
    if (!is.null(intrinsic) && nrow(intrinsic)) {
      ag <- dplyr::inner_join(final, intrinsic, by = c("patient_id", "moving")) |>
        dplyr::filter(is.finite(intrinsic), is.finite(dice_matched))
      if (nrow(ag) > 1)
        figs[["04_intrinsic_vs_overlap"]] <-
          ggplot(ag, aes(intrinsic, dice_matched, colour = patient_id)) +
          geom_point(size = 3, alpha = .85) +
          scale_colour_oi(name = "patient") +
          lab(title = "Self-reported error vs the independent overlap check",
              subtitle = paste("One point per moving slide at the final stage. A method grading",
                               "itself can look good while nucleus overlap disagrees; that gap",
                               "is what this plots."),
              x = ag$what[1], y = "matched-nucleus Dice")
    }
  }

  # -- §5 phenotyping calibration ----------------------------------------------
  pq <- tables$pheno_qc
  if (nrow(pq) && any(is.finite(pq$chosen_alpha))) {
    figs[["05_phenotype_alpha"]] <-
      ggplot(pq, aes(stats::reorder(patient_id, chosen_alpha), chosen_alpha,
                     fill = crc_ran)) +
      geom_col(width = .7) + coord_flip() +
      geom_hline(aes(yintercept = alpha_target), linetype = "dashed", colour = REF_LINE) +
      scale_fill_manual(values = c(`TRUE` = oi[1], `FALSE` = oi[2]),
                        name = "CRC ran", labels = c(`TRUE` = "yes", `FALSE` = "reporting only")) +
      lab(title = "Conformal risk budget actually used, per patient",
          subtitle = paste("Dashed line is the requested alpha_target. Below it, CRC tightened",
                           "the budget; 'reporting only' means calibration was too thin to",
                           "certify a risk level at all."),
          x = NULL, y = "chosen alpha")
  }

  # -- §5b constraint audit: does the data contradict the panel? ----------------
  au <- tables$audit
  if (nrow(au) && all(c("markers", "observed", "nominal") %in% names(au)))
    figs[["05b_constraint_audit"]] <-
      au |>
      dplyr::filter(is.finite(observed)) |>
      ggplot(aes(observed, stats::reorder(markers, observed), colour = patient_id)) +
      geom_point(aes(x = nominal), shape = 124, size = 4, colour = REF_LINE) +
      geom_point(size = 2.5, alpha = .85) +
      scale_colour_oi(name = "patient") +
      lab(title = "Constraint audit — observed vs asserted co-expression rate",
          subtitle = paste("Bars mark the rate the panel asserts; points are what the slide shows.",
                           "A point far right of its bar means spillover or a mis-set gate on",
                           "those markers."),
          x = "co-expression rate", y = NULL)

  figs
}
