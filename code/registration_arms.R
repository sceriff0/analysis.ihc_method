# =============================================================================
# registration_arms.R  —  the registration ARM sweep, run on the REAL slides
#
# Three pages now read registration accuracy, and they are not interchangeable:
#
#   benchmark_registration.Rmd  mirage's benchmarks/ sweep on SYNTHETIC images with a
#                              known injected offset. Measures cost and scaling.
#   run_qc.Rmd                 ONE run's QC on the study slides. "Was this cohort
#                              registered well enough to analyse?"
#   registration_arms.Rmd      THIS file. The same study slides registered N times,
#                              once per configuration. "Which configuration do we
#                              ship, and what does it buy?"
#
# The arm sweep is the one the manuscript's Fig 4(b)/(c) actually needs: an arm
# ranking measured on real tissue rather than on a synthetic offset.
#
# THE AXES
#   registration_method  the BACKEND: `valis` (default) or `tiled` (STARE, JVM-free,
#                  internally tiled). This is a different axis from the two below, not
#                  a third level of them: memory_mode and reg_micro_reg are VALIS-only
#                  params, so a tiled arm has NEITHER and carries NA for both. A tiled
#                  arm is "the other backend at its defaults", one point of comparison,
#                  not a cell of the preset x depth grid.
#   memory_mode    VALIS accuracy preset. `low` = BRISK/RANSAC at small dims,
#                  `high` = SuperPoint/SuperGlue at larger dims. NOTE: these are
#                  DIFFERENT FEATURE MATCHERS, not one matcher at two resolutions.
#   reg_micro_reg  micro-registration DEPTH, nested, not a boolean:
#                    0 = none
#                    1 = micro-rigid only (refines slide.M)
#                    2 = + micro non-rigid (register_micro)
#
# THE BACKENDS DO NOT SHARE A STAGE VOCABULARY. lib/WarpBackends.groovy:
#   valis -> native, rigid, non_rigid, micro
#   tiled -> native, rigid, refined
# The segmentation-overlap SCORER is method-agnostic (bin/warp_seg_qc.py takes
# `--method` and builds its warper from either a VALIS registrar pickle or a STARE
# transform manifest), so the metric itself IS comparable across backends — which is
# the whole reason a tiled arm can join this page at all. What is not comparable is
# the stage axis: only `native` is a shared spelling with a shared meaning. `rigid`
# is shared as a WORD and not as an operation (VALIS: affine, composed with
# micro-rigid at depth >= 1; STARE: the coarse global anchor before mesh refinement).
#
# The backends also report their OWN error in different units and different files:
# VALIS writes rTRE (a fraction of the image diagonal) or a raw distance into
# registered/summary/*.csv; STARE writes TRE in PIXELS into qc/registration/*_tre.json.
# Those never share an axis here. They are separate figures on purpose.
#
# WHY `rigid` CANNOT BE COMPARED ACROSS ARMS — the trap this file exists to close.
# mirage's staged QC defines the `rigid` stage as the rigid transform AFTER
# MicroRigidRegistrar refined it (docs/parameters.md: "At >=1 the QC `rigid` stage
# means affine o micro-rigid"). So `rigid` means:
#     micro_reg = 0  ->  affine alone
#     micro_reg >= 1 ->  affine o micro-rigid
# Plotting "rigid" across a micro_reg-crossed sweep therefore plots two different
# transforms on one axis and reads as a micro-registration effect that is really a
# definition change. arm_comparable_stages() refuses to return `rigid` for a mixed
# sweep, and comparable_stage_note() says why, so the guard is visible in the report
# rather than buried here.
#
# AND WHY `native` IS ONLY NEARLY COMPARABLE. `native` is the untransformed
# segmentation, so it looks like a clean no-registration baseline — but the cell
# PAIRING is established once at the `rigid` anchor, and arms with different rigid
# transforms pair different cells. `native` is then scored over a different pair set
# per arm. Close enough to read as a baseline, not close enough to difference against.
#
# THE ARM'S ANSWER IS ITS LAST STAGE, and which stage that is depends on the arm:
# micro_reg = 0 and 1 emit no `micro` stage at all (at 1 the micro-rigid refinement
# has already been folded into `rigid`), so ranking arms on "the micro stage" would
# silently drop four of six arms. arm_final_stage() takes the last stage each run
# actually reported, which is that run's shipped output by construction.
#
# LAYOUT
#   data/registration_arms/<arm>/<patient>/qc/registration/*_seg_qc.json
#   data/registration_arms/<arm>/<patient>/registered/summary/*.csv
#   data/registration_arms/arms.csv        (optional, and strongly preferred)
#
# Every reader here is run_qc.R's, called once per arm directory. The artifacts are
# the same artifacts; only the number of trees changed. Adding a parallel parser
# would have been a second place for the schema to drift.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(fs)
  library(purrr)
  library(tibble)
})

source(here::here("code", "run_qc.R"))    # the readers, and QC_STAGE_LEVELS

ARMS_DIR <- here::here("data", "registration_arms")

ARM_CAPTION <- "mirage staged registration QC (reg_qc = 2), study slides, one run per arm."

# --- Which arm is which ------------------------------------------------------
# MANIFEST FIRST, name-parsing second. A mislabelled arm does not fail: it produces
# a clean figure with the conclusion inverted, which is the worst failure mode
# available here. `arms.csv` lets the person who launched the runs state the mapping
# instead of encoding it in a directory name and hoping the regex agrees.
#
#   arm_dir,memory_mode,micro_reg
#   valis_high_micro2,high,2
#
# Without one, the directory name is parsed for `high`/`low` and for a micro depth
# written any of the usual ways (micro2, micro_2, micro-reg-2, mr2). Anything that
# does not parse keeps the directory name as its label and gets NA knobs, so it
# still appears in the figures — unlabelled, rather than dropped.
.parse_arm_dir <- function(nm) {
  low <- tolower(nm)
  # The backend first: a tiled arm has no preset and no micro depth, so reading those
  # off its name would invent knob values it was never run with.
  backend <- if (grepl("tiled|stare", low)) "tiled" else "valis"
  if (backend == "tiled")
    return(tibble::tibble(backend = backend,
                          memory_mode = NA_character_, micro_reg = NA_integer_))
  mode <- if (grepl("high", low) && !grepl("low", low)) "high"
          else if (grepl("low", low)) "low"
          else NA_character_
  m <- regmatches(low, regexec("(?:micro|mr)[^0-9]{0,6}([0-2])", low))[[1]]
  micro <- if (length(m) == 2) as.integer(m[2]) else NA_integer_
  tibble::tibble(backend = backend, memory_mode = mode, micro_reg = micro)
}

arm_manifest <- function(root = ARMS_DIR) {
  if (!fs::dir_exists(root)) {
    warning("registration arms: no directory at ", root)
    return(tibble::tibble())
  }
  dirs <- as.character(fs::dir_ls(root))
  dirs <- dirs[fs::dir_exists(dirs)]
  # Only directories that actually hold a run: a stray `figures/` or `logs/` beside
  # the arms is not an arm, and counting it as one puts an empty box in every panel.
  dirs <- dirs[vapply(dirs, function(d) length(.qc_patient_dirs(d)) > 0, logical(1))]
  if (!length(dirs)) {
    warning("registration arms: ", root, " holds no directory with mirage QC artifacts")
    return(tibble::tibble())
  }

  out <- tibble::tibble(arm_dir = fs::path_file(dirs), path = dirs) |>
    dplyr::bind_cols(purrr::map_dfr(fs::path_file(dirs), .parse_arm_dir))

  man_path <- file.path(root, "arms.csv")
  if (file.exists(man_path)) {
    man <- tryCatch(readr::read_csv(man_path, show_col_types = FALSE, progress = FALSE),
                    error = function(e) NULL)
    if (!is.null(man) && "arm_dir" %in% names(man)) {
      keep <- intersect(c("arm_dir", "backend", "memory_mode", "micro_reg", "label"), names(man))
      out  <- out |>
        dplyr::select(-dplyr::any_of(setdiff(keep, "arm_dir"))) |>
        dplyr::left_join(dplyr::select(man, dplyr::all_of(keep)), by = "arm_dir")
      message("registration arms: labels taken from ", man_path)
    } else {
      warning("registration arms: ", man_path, " has no `arm_dir` column — ignoring it")
    }
  }

  if (!"label" %in% names(out)) out$label <- NA_character_
  out |>
    dplyr::mutate(
      micro_reg = suppressWarnings(as.integer(micro_reg)),
      backend   = dplyr::coalesce(backend, "valis"),
      # A tiled arm is named for its backend, not for knobs it does not have.
      arm = dplyr::coalesce(label, dplyr::case_when(
        backend == "tiled" ~ "tiled (STARE, defaults)",
        !is.na(memory_mode) & !is.na(micro_reg) ~
          sprintf("%s / micro %d", memory_mode, micro_reg),
        TRUE ~ arm_dir))) |>
    # VALIS arms first, grouped by preset then depth; the tiled comparator last, since
    # it is a different backend rather than another cell of the grid.
    dplyr::arrange(backend != "valis", dplyr::desc(memory_mode), micro_reg, arm_dir)
}

# --- Reading every arm -------------------------------------------------------
# run_qc.R's reader, once per arm tree. `arm` and the two knobs ride along so every
# downstream figure can facet on them without re-deriving anything.
read_arms_seg_qc <- function(manifest = arm_manifest()) {
  if (nrow(manifest) == 0) return(tibble::tibble())
  purrr::pmap_dfr(manifest, function(arm_dir, path, backend, memory_mode, micro_reg,
                                     label, arm) {
    d <- read_seg_qc(path)
    if (nrow(d) == 0) {
      warning("registration arms: no seg_qc under ", path)
      return(tibble::tibble())
    }
    dplyr::mutate(d, arm = arm, arm_dir = arm_dir, backend = backend,
                  memory_mode = memory_mode, micro_reg = micro_reg, .before = 1)
  })
}

read_arms_valis <- function(manifest = arm_manifest()) {
  if (nrow(manifest) == 0 || !"backend" %in% names(manifest)) return(tibble::tibble())
  # VALIS arms only, by construction: a tiled run writes no registered/summary/*.csv.
  # Its backend-native error is read by read_arms_stare_tre() instead, in pixels.
  purrr::pmap_dfr(dplyr::filter(manifest, backend == "valis"),
                  function(arm_dir, path, backend, memory_mode, micro_reg, label, arm) {
    d <- read_valis_summary(path)
    if (nrow(d) == 0) return(tibble::tibble())
    dplyr::mutate(d, arm = arm, arm_dir = arm_dir, backend = backend,
                  memory_mode = memory_mode, micro_reg = micro_reg, .before = 1)
  })
}

# --- The comparability guard -------------------------------------------------
# Which stages may be put on one axis across THIS set of arms.
#
# `native` and each run's final stage always may. `rigid` and `non_rigid` may only
# when every arm shares one micro_reg, because at micro_reg >= 1 `rigid` silently
# absorbs the micro-rigid refinement (see the header). Returning a vector rather
# than a boolean lets the caller filter instead of remembering the rule.
arm_comparable_stages <- function(seg) {
  if (nrow(seg) == 0) return(character(0))
  # Backends first. Crossing them is the stronger restriction: the two stage lists
  # share only `native` as both a spelling AND a meaning, so no amount of matching
  # micro depth makes `rigid` comparable between a VALIS and a STARE run.
  if (dplyr::n_distinct(seg$backend %||% "valis") > 1) return("native")
  depths <- unique(stats::na.omit(seg$micro_reg))
  present <- intersect(QC_STAGE_LEVELS, unique(as.character(seg$stage)))
  if (length(depths) <= 1) return(present)
  c("native")
}

comparable_stage_note <- function(seg) {
  if (dplyr::n_distinct(seg$backend %||% "valis") > 1)
    return(paste0(
      "This sweep crosses REGISTRATION BACKENDS (",
      paste(sort(unique(seg$backend)), collapse = " and "),
      "). They do not share a stage vocabulary — VALIS reports `native → rigid →",
      " non_rigid → micro`, STARE reports `native → rigid → refined` — and `rigid`",
      " is a shared word rather than a shared operation (VALIS: affine, composed",
      " with micro-rigid at depth ≥ 1; STARE: the coarse global anchor before mesh",
      " refinement). Only `native` and each run's FINAL stage may be put on one",
      " axis. The segmentation-overlap METRIC is backend-agnostic — mirage's scorer",
      " builds its warper from either a VALIS registrar or a STARE manifest — which",
      " is what makes the final-stage comparison legitimate at all."))
  depths <- sort(unique(stats::na.omit(seg$micro_reg)))
  if (length(depths) <= 1)
    return(paste0("All arms share `reg_micro_reg = ", depths[1] %||% "?",
                  "`, so every stage means the same thing across arms and all",
                  " stages are directly comparable."))
  paste0(
    "This sweep crosses `reg_micro_reg` = ", paste(depths, collapse = ", "),
    ". mirage defines the QC `rigid` stage as the rigid transform *after*",
    " `MicroRigidRegistrar` refined it, so `rigid` means affine alone at depth 0 and",
    " affine ∘ micro-rigid at depth ≥ 1 — two different transforms. Per-stage",
    " comparisons across arms are therefore restricted to `native` (untransformed)",
    " and to each run's FINAL stage, which is that run's shipped output. The stage",
    " ladder is still read *within* an arm, where the definition is fixed.")
}

# Each run's last reported stage: the transform that arm actually ships. Not a
# constant across arms — depth 0 and 1 emit no `micro` stage at all — which is
# exactly why it is derived per (arm, slide) instead of being named once.
arm_final_stage <- function(seg) {
  if (nrow(seg) == 0) return(seg)
  # Ordered by the run's OWN stage_index, never by the shared factor levels. The two
  # backends' vocabularies interleave in any single global ordering, so `as.integer(stage)`
  # would be comparing a STARE stage's position against a VALIS one. stage_index comes
  # from each run's `stage_order`, so "the last stage" is that run's own last stage.
  ord <- if ("stage_index" %in% names(seg)) seg$stage_index else as.integer(seg$stage)
  seg |>
    dplyr::mutate(.ord = ord) |>
    dplyr::filter(!is.na(stage)) |>
    dplyr::group_by(arm, arm_dir, backend, memory_mode, micro_reg, patient_id, moving) |>
    dplyr::slice_max(order_by = .ord, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(-.ord) |>
    dplyr::rename(final_stage = stage)
}

# --- STARE's own error -------------------------------------------------------
# The tiled backend reports its intrinsic TRE in qc/registration/*_tre.json, in
# PIXELS, with a per-tile breakdown VALIS has no equivalent of. Kept in its own
# reader and its own figure: putting a pixel TRE on the same axis as VALIS's
# diagonal-relative rTRE would be a unit error dressed as a comparison.
read_arms_stare_tre <- function(manifest = arm_manifest()) {
  if (nrow(manifest) == 0 || !"backend" %in% names(manifest)) return(tibble::tibble())
  tiled <- dplyr::filter(manifest, backend == "tiled")
  if (nrow(tiled) == 0) return(tibble::tibble())
  purrr::pmap_dfr(tiled, function(arm_dir, path, backend, memory_mode, micro_reg,
                                  label, arm) {
    d <- read_stare_tre(path)
    if (nrow(d) == 0) return(tibble::tibble())
    dplyr::mutate(d, arm = arm, arm_dir = arm_dir, backend = backend, .before = 1)
  })
}

# Recover a minimal manifest from an already-read seg frame, for a caller that passed
# no manifest. Only the tiled arms matter, and only their paths.
read_arms_stare_tre_from <- function(seg) {
  if (nrow(seg) == 0 || !"backend" %in% names(seg)) return(tibble::tibble())
  tl <- dplyr::distinct(dplyr::filter(seg, backend == "tiled"), arm, arm_dir, backend)
  if (nrow(tl) == 0) return(tibble::tibble())
  read_arms_stare_tre(dplyr::mutate(tl, path = file.path(ARMS_DIR, arm_dir),
                                    memory_mode = NA_character_,
                                    micro_reg = NA_integer_, label = NA_character_))
}

read_arms_stare_tiles <- function(manifest = arm_manifest()) {
  if (nrow(manifest) == 0 || !"backend" %in% names(manifest)) return(tibble::tibble())
  tiled <- dplyr::filter(manifest, backend == "tiled")
  if (nrow(tiled) == 0) return(tibble::tibble())
  purrr::pmap_dfr(tiled, function(arm_dir, path, backend, memory_mode, micro_reg,
                                  label, arm) {
    d <- read_stare_tiles(path)
    if (nrow(d) == 0) return(tibble::tibble())
    dplyr::mutate(d, arm = arm, arm_dir = arm_dir, .before = 1)
  })
}

# --- Figures -----------------------------------------------------------------
# Same contract as the other figure builders: a named list, skipping any figure
# whose input is absent, so the page renders against a partial sweep.
build_arm_figs <- function(seg = read_arms_seg_qc(), valis = read_arms_valis(),
                           manifest = NULL) {
  figs <- list()
  if (nrow(seg) == 0) return(figs)

  # Arm display order comes from the manifest when the caller has one; otherwise from
  # the data itself. Never from a re-read of the global root — a builder that reaches
  # past its arguments for a global is right only when the global happens to match.
  fin      <- arm_final_stage(seg)
  arm_lvls <- if (!is.null(manifest) && nrow(manifest)) unique(manifest$arm)
              else unique(seg$arm)
  .arm_f   <- function(x) factor(x, levels = intersect(arm_lvls, unique(x)))

  # -- 1. THE ARM RANKING. Residual displacement in microns at each arm's final
  # transform. Physical units, so it is the number to quote; one point per slide,
  # because with a handful of slides the points are the evidence.
  if (any(is.finite(fin$disp_um_p50))) {
    d <- dplyr::filter(fin, is.finite(disp_um_p50)) |> dplyr::mutate(arm = .arm_f(arm))
    figs[["01_final_residual_um_by_arm"]] <-
      ggplot(d, aes(arm, disp_um_p50)) +
      geom_boxplot(outlier.shape = NA, width = .5, colour = "grey35") +
      scale_x_discrete(labels = label_n(d$arm, sep = " ")) +
      geom_jitter(aes(colour = .arm_kind(backend, micro_reg)), width = .12, height = 0,
                  alpha = .85, size = 2.2) +
      scale_colour_arm() +
      coord_flip() +
      labs(title = "Residual alignment error at each arm's final transform",
           subtitle = paste("Matched-nucleus centroid residual, median per slide.",
                            "Physical units; lower = tighter. One point per slide."),
           x = NULL, y = "residual displacement, median (µm)", caption = ARM_CAPTION)
  }

  # -- 2. The same ranking by overlap Dice — a different measurement of the same
  # alignment, so agreement between figures 1 and 2 is evidence, not restatement.
  if (any(is.finite(fin$dice_matched))) {
    d <- dplyr::filter(fin, is.finite(dice_matched)) |> dplyr::mutate(arm = .arm_f(arm))
    figs[["02_final_dice_by_arm"]] <-
      ggplot(d, aes(arm, dice_matched)) +
      geom_boxplot(outlier.shape = NA, width = .5, colour = "grey35") +
      scale_x_discrete(labels = label_n(d$arm, sep = " ")) +
      geom_jitter(aes(colour = .arm_kind(backend, micro_reg)), width = .12, height = 0,
                  alpha = .85, size = 2.2) +
      scale_colour_arm() +
      coord_flip() +
      labs(title = "Matched-nucleus Dice at each arm's final transform",
           subtitle = "Higher = better. Same arms and slides as the residual figure.",
           x = NULL, y = "Matched-nucleus Dice (unitless, 0-1)", caption = ARM_CAPTION)
  }

  # -- 3. The stage ladder, WITHIN each arm. This is where the per-stage story is
  # legible, because inside one arm `rigid` has a fixed meaning. Faceting by arm is
  # the guard made visual: the stages are never put on one shared axis.
  if (any(is.finite(seg$disp_um_p50))) {
    d <- dplyr::filter(seg, is.finite(disp_um_p50)) |> dplyr::mutate(arm = .arm_f(arm))
    figs[["03_stage_ladder_within_arm"]] <-
      ggplot(d, aes(stage, disp_um_p50, group = interaction(patient_id, moving))) +
      geom_line(alpha = .35) + geom_point(size = 1.6, alpha = .8, colour = oi[1]) +
      facet_wrap(~ arm) +
      scale_y_log10() +
      labs(title = "What each stage bought, within each arm",
           subtitle = paste("One line per moving slide. Read DOWN the ladder inside a",
                            "panel; do not read `rigid` across panels — at",
                            "micro_reg ≥ 1 it already contains micro-rigid."),
           x = NULL, y = "residual displacement, median (µm, log)", caption = ARM_CAPTION) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  }

  # -- 4. Is the pairing thick enough to trust the arm at all? mirage's own rule:
  # below ~0.5 the later stages are scored on a biased subset of cells that happened
  # to land close. An arm that wins on residual while pairing thinly has not won.
  if (any(is.finite(seg$pair_fraction))) {
    d <- seg |>
      dplyr::filter(is.finite(pair_fraction)) |>
      dplyr::distinct(arm, micro_reg, patient_id, moving, pair_fraction) |>
      dplyr::mutate(arm = .arm_f(arm))
    figs[["04_pair_fraction_by_arm"]] <-
      ggplot(d, aes(arm, pair_fraction)) +
      geom_hline(yintercept = 0.5, linetype = "dashed", colour = REF_LINE) +
      geom_boxplot(outlier.shape = NA, width = .5, colour = "grey35") +
      scale_x_discrete(labels = label_n(d$arm, sep = " ")) +
      geom_jitter(width = .12, height = 0, alpha = .8, size = 2, colour = oi[3]) +
      coord_flip() + ylim(0, 1) +
      labs(title = "How much of the slide each arm could actually pair",
           subtitle = paste("Fraction of nuclei matched at the rigid anchor. Below the",
                            "dashed 0.5 line the later stages are measured on a biased",
                            "subset — check this before believing a ranking."),
           x = NULL, y = "Pair fraction (unitless, 0-1)", caption = ARM_CAPTION)
  }

  # -- 5. The two knobs, separated. The ranking figures order arms by outcome; this
  # one asks which KNOB moved the outcome, which is the question a reader has next.
  if (any(is.finite(fin$disp_um_p50)) &&
      dplyr::n_distinct(stats::na.omit(fin$memory_mode)) > 1) {
    # VALIS arms only: the two knobs do not exist on the tiled backend, so including
    # it here would put a point on a grid it was never run on.
    d <- dplyr::filter(fin, backend == "valis",
                       is.finite(disp_um_p50), !is.na(memory_mode), !is.na(micro_reg))
    if (nrow(d)) {
      figs[["05_knob_effects"]] <-
        ggplot(d, aes(factor(micro_reg), disp_um_p50, colour = memory_mode)) +
        geom_boxplot(outlier.shape = NA, width = .55,
                     position = position_dodge(.7)) +
        geom_point(position = position_jitterdodge(jitter.width = .12, dodge.width = .7),
                   alpha = .8, size = 1.9) +
        scale_colour_manual(values = unname(oi[c(1, 2)]), name = "memory_mode") +
        labs(title = "Which knob moved the result",
             subtitle = paste("Final-transform residual by micro-registration depth,",
                              "split by VALIS accuracy preset. Note the presets use",
                              "different feature MATCHERS (BRISK/RANSAC vs",
                              "SuperPoint/SuperGlue), not one matcher at two scales."),
             x = "reg_micro_reg (0 = none, 1 = micro-rigid, 2 = + micro non-rigid)",
             y = "residual displacement, median (µm)", caption = ARM_CAPTION)
    }
  }

  # -- 6. Did micro-registration ever make it WORSE? mirage caught-and-continues a
  # failed micro-registration, so a regression is silent by design; delta_vs_anchor
  # is where it surfaces. A POSITIVE displacement delta is worse.
  dv <- dplyr::filter(seg, is.finite(d_disp_um_vs_rigid), stage != "rigid")
  if (nrow(dv)) {
    figs[["06_delta_vs_rigid_anchor"]] <-
      ggplot(dplyr::mutate(dv, arm = .arm_f(arm)),
             aes(stage, d_disp_um_vs_rigid, colour = arm)) +
      geom_hline(yintercept = 0, colour = REF_LINE) +
      geom_jitter(width = .15, height = 0, alpha = .8, size = 2) +
      scale_colour_manual(values = rep_len(oi_ext, dplyr::n_distinct(dv$arm)), name = NULL) +
      labs(title = "Change against the rigid anchor, per stage",
           subtitle = paste("Residual minus the rigid stage's.",
                            "ABOVE the zero line means that stage made alignment WORSE",
                            "— the failure mode micro-registration hides."),
           x = NULL, y = "Δ residual vs rigid (µm)", caption = ARM_CAPTION)
  }

  # -- 7. VALIS grading itself, ONE figure: the stage axis, faceted by arm.
  #
  # This replaces a per-arm summary at each arm's final stage. That collapsed an arm to
  # one number and threw away the ladder, which is the interesting part — and the ladder
  # is what makes the arms comparable at all, because the STAGE MEANINGS differ by depth
  # and faceting is what fixes them (same reason figure 3 facets).
  #
  # Three stages, two columns. VALIS's error_df is `from`/`filename`, `rigid_D`,
  # `non_rigid_D`; micro-registration has no column of its own because it UPDATES the
  # non-rigid field, so the micro value lives in the difference between the pre-micro and
  # final files rather than in a column. valis_error_long() does that reconstruction.
  #
  # The depth-0 and depth-1 arms therefore show NO micro box, and that blank is the
  # finding: micro-registration did not run. A duplicated non_rigid box would instead
  # read as "micro bought nothing".
  vl <- if (nrow(valis)) valis_error_long(valis) else tibble::tibble()
  if (nrow(vl) && "arm" %in% names(vl)) {
    vf <- dplyr::filter(vl, is.finite(error)) |> dplyr::mutate(arm = .arm_f(arm))
    if (nrow(vf)) {
      med <- vf |>
        dplyr::group_by(arm, stage) |>
        dplyr::summarise(error = stats::median(error, na.rm = TRUE), .groups = "drop")
      figs[["07_valis_intrinsic_by_arm"]] <-
        ggplot(vf, aes(stage, error)) +
        geom_boxplot(outlier.shape = NA, width = .5, colour = "grey35") +
        scale_x_discrete(labels = label_n(vf$stage)) +
        geom_jitter(width = .12, height = 0, alpha = .8, size = 1.8, colour = oi[4]) +
        geom_text(data = med, aes(label = signif(error, 3)), vjust = -1.0,
                  size = pt_text(6.5), colour = "grey15") +
        facet_wrap(~ arm) +
        labs(title = "VALIS's own reported error, by stage, within each arm",
             subtitle = paste("Independent of every segmentation-overlap metric above:",
                              "VALIS grading itself from feature correspondences.",
                              "Label = median. A missing `micro` box means the arm wrote",
                              "no pre-micro summary (reg_micro_reg < 2), so",
                              "micro-registration never ran — not that it gained nothing.",
                              "Read DOWN a panel; the stage meanings differ across depths,",
                              "which is why this facets rather than sharing an axis."),
             x = NULL, y = vf$metric[1] %||% "VALIS rTRE / distance",
             caption = ARM_CAPTION) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
    }
  }

  # -- 8. STARE's own reported error, in its own units and its own figure. Kept
  # apart from figure 7 (VALIS's rTRE) deliberately: pixels versus a fraction of the
  # image diagonal is not one axis, and a reader who sees them side by side will
  # compare the numbers anyway.
  st <- if (!is.null(manifest)) read_arms_stare_tre(manifest)
        else read_arms_stare_tre_from(seg)
  if (nrow(st)) {
    long <- st |>
      dplyr::select(arm, patient_id, moving, rigid_p50, after_p50) |>
      tidyr::pivot_longer(c(rigid_p50, after_p50), names_to = "stage", values_to = "tre_px") |>
      dplyr::mutate(stage = factor(dplyr::recode(stage,
                      rigid_p50 = "rigid anchor", after_p50 = "after refinement"),
                      levels = c("rigid anchor", "after refinement"))) |>
      dplyr::filter(is.finite(tre_px))
    if (nrow(long))
      figs[["08_stare_intrinsic_tre_px"]] <-
        ggplot(long, aes(stage, tre_px, group = interaction(patient_id, moving))) +
        geom_line(alpha = .35) +
        geom_point(size = 2.2, alpha = .85, colour = oi[5 %% length(oi) + 1]) +
        labs(title = "STARE's own reported error, before and after mesh refinement",
             subtitle = paste("Tiled backend only, in PIXELS — not comparable with",
                              "VALIS's diagonal-relative rTRE. One line per moving slide."),
             x = NULL, y = "TRE (px)", caption = ARM_CAPTION)
  }

  # -- 9. What only the tiled backend can show: WHERE on the slide the error is.
  # VALIS reports one number per slide; STARE reports per tile, so a systematic
  # regional failure is visible instead of averaged away.
  tiles <- if (!is.null(manifest)) read_arms_stare_tiles(manifest) else tibble::tibble()
  xy    <- intersect(c("x", "y"), names(tiles))
  tre_c <- intersect(c("tre_px", "residual_px", "tre"), names(tiles))[1]
  if (nrow(tiles) && length(xy) == 2 && !is.na(tre_c)) {
    d <- dplyr::filter(tiles, is.finite(.data[[tre_c]]))
    if (nrow(d))
      figs[["09_stare_tile_error_map"]] <-
        ggplot(d, aes(x, y, fill = .data[[tre_c]])) +
        geom_tile() +
        scale_fill_seq(name = "TRE (px)") +
        scale_y_reverse() + coord_fixed() +
        facet_wrap(~ patient_id) +
        labs(title = "Where the residual error sits, per tile",
             subtitle = paste("Tiled backend only — the spatial breakdown VALIS has no",
                              "equivalent of. A bright region is a local failure a",
                              "slide-level median hides."),
             x = NULL, y = NULL, caption = ARM_CAPTION) +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid = element_blank())
  }

  figs
}

# How to colour an arm: by BACKEND, with the VALIS micro depth as the sub-level. A
# tiled arm has no depth, so colouring by depth alone would draw it as NA and read as
# missing data rather than as the other backend. ARM_KIND_COLS itself lives in
# plot_theme.R with the other recurring category palettes; use scale_colour_arm().

.arm_kind <- function(backend, micro_reg) {
  lab <- ifelse(backend == "tiled", "tiled (STARE)",
                ifelse(is.na(micro_reg), "valis",
                       paste0("valis · micro ", micro_reg)))
  factor(lab, levels = names(ARM_KIND_COLS))
}

# The ranking, as a table — Additional-file material, and the numbers to quote.
arm_ranking_table <- function(seg = read_arms_seg_qc()) {
  fin <- arm_final_stage(seg)
  if (nrow(fin) == 0) return(tibble::tibble())
  fin |>
    dplyr::group_by(arm, backend, memory_mode, micro_reg) |>
    dplyr::summarise(
      n_slides      = dplyr::n(),
      final_stage   = paste(sort(unique(as.character(final_stage))), collapse = "/"),
      disp_um_p50   = stats::median(disp_um_p50, na.rm = TRUE),
      dice_matched  = stats::median(dice_matched, na.rm = TRUE),
      pair_fraction = stats::median(pair_fraction, na.rm = TRUE),
      .groups = "drop") |>
    dplyr::arrange(disp_um_p50)
}
