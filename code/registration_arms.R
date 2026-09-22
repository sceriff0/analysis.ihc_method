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
#   registration_method  the BACKEND: `valis` (default), `tiled` (STARE, JVM-free,
#                  internally tiled) or `ashlar` (labsyspharm ASHLAR, the external
#                  baseline). This is a different axis from the two below, not a third
#                  level of them: memory_mode and reg_micro_reg are VALIS-only params, so
#                  a tiled or ashlar arm has NEITHER and carries NA for both. Each is
#                  "another backend", a point of comparison against the whole preset x
#                  depth grid, not a cell of it.
#
#                  The ashlar arms fan out over GRID GRANULARITY (`ashlar_tile1024`,
#                  `ashlar_tile4096`) because ASHLAR takes one independent shift per
#                  tile: a finer grid buys it more local freedom, so tile size is a
#                  FAIRNESS knob against STARE's reg_tiled_tile, not a cost knob.
#                  ASHLAR attempts NO non-rigid warp at all, so read it against VALIS's
#                  `rigid` stage for the like-for-like number and against `micro` to see
#                  what non-rigid buys; reporting only the second overstates VALIS.
#   memory_mode    VALIS accuracy preset. ONE feature matcher at two resolutions:
#                  every tier uses SuperPoint + SuperGlue (mirage
#                  bin/utils/valis_config.py, MEMORY_PRESETS); `low` detects on a
#                  512 px processed image and fits the non-rigid field at 512 px,
#                  `medium` at 1024 / 1024 px, `high` at 2048 / 2048 px. An earlier
#                  note here — and a comment in mirage's benchmarks/configs/arms.yaml —
#                  said `low` was BRISK/RANSAC. It was never true of the shipped code;
#                  the manuscript legend's "differing in feature-matching resolution,
#                  not in the feature matcher" is the correct statement.
#   reg_micro_reg  micro-registration DEPTH, nested, not a boolean:
#                    0 = none
#                    1 = micro-rigid only (refines slide.M)
#                    2 = + micro non-rigid (register_micro)
#
# THE BACKENDS DO NOT SHARE A STAGE VOCABULARY. lib/WarpBackends.groovy:
#   valis  -> native, rigid, non_rigid, micro
#   tiled  -> native, rigid, refined
#   ashlar -> native, rigid, refined
# The segmentation-overlap SCORER is method-agnostic (bin/warp_seg_qc.py takes
# `--method` and builds its warper from either a VALIS registrar pickle or a STARE
# transform manifest), so the metric itself IS comparable across backends — which is
# the whole reason a tiled arm can join this page at all. ASHLAR shares the tiled
# vocabulary because it shares the ARTIFACT: bin/ashlar_solve.py rewrites ASHLAR's
# per-tile placements into the same M0 + mesh manifest STARE emits, which is the only
# reason an ashlar arm can be ranked here rather than in a separate table of
# residual-TRE numbers that share no column with this one. What is not comparable is
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
  # `scales` for pseudo_log_trans() on the signed delta figure. It is a hard
  # dependency of ggplot2, so this only makes the use declared, not conditional.
  library(scales)
})

source(here::here("code", "run_qc.R"))    # the readers, and QC_STAGE_LEVELS

ARMS_DIR <- here::here("data", "registration_arms")

ARM_CAPTION <- "mirage staged registration QC (reg_qc = 2), study slides, one run per arm."

# The unit of every final-transform panel, said the same way in each. A slide is one
# staining round — one set of channels — and mirage scores each moving slide against
# its patient's reference, so one QC record is one PAIR of channel sets.
PAIR_UNIT_NOTE <- paste("One point per channel pair (reference vs moving slide),",
                        "pooled over all pairs and all patients.")

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
  # The backend first: a non-VALIS arm has no preset and no micro depth, so reading those
  # off its name would invent knob values it was never run with.
  #
  # ASHLAR is tested BEFORE the tiled pattern, and the `else` is VALIS -- so a backend
  # added upstream without a rule here is silently read as VALIS, complete with a preset
  # and a depth parsed out of a name that never had them. That is not hypothetical: the
  # ashlar arm dirs are `ashlar_tile1024` / `ashlar_tile4096`, and the old two-way test
  # would have called both of them VALIS with memory_mode = NA.
  backend <- if (grepl("ashlar", low)) "ashlar"
             else if (grepl("tiled|stare", low)) "tiled"
             else "valis"
  if (backend %in% c("tiled", "ashlar"))
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
      # A non-VALIS arm is named for its backend, not for knobs it does not have. The
      # ashlar arms differ only by grid granularity, which IS in the directory name, so
      # the fallback recovers it rather than collapsing both arms onto one label -- two
      # identically-labelled boxes would read as a duplicated bar, not as two arms.
      arm = dplyr::coalesce(label, dplyr::case_when(
        backend == "tiled" ~ "tiled (STARE, defaults)",
        backend == "ashlar" ~ {
          t <- regmatches(tolower(arm_dir),
                          regexec("tile[^0-9]{0,2}([0-9]+)", tolower(arm_dir)))
          t <- vapply(t, function(m) if (length(m) == 2) m[2] else NA_character_,
                      character(1))
          ifelse(is.na(t), "ashlar", sprintf("ashlar (tile %s)", t))
        },
        !is.na(memory_mode) & !is.na(micro_reg) ~
          sprintf("%s / micro %d", memory_mode, micro_reg),
        TRUE ~ arm_dir))) |>
    # VALIS arms first, grouped by preset then depth; the other BACKENDS after them,
    # alphabetically by directory (ashlar_tile1024, ashlar_tile4096, tiled_defaults),
    # since each is a different backend rather than another cell of the VALIS grid.
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
  # transform. Physical units, so it is the number to quote. One point per CHANNEL
  # PAIR (one reference-vs-moving QC record), pooled over every pair and every
  # patient: that is the main-text panel. The same points split by patient and by
  # pair are the S-figures at the end of this builder.
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
           # The unit note takes its own line in all three main panels: run on, the
           # subtitle overruns a double-column panel and the clipped part is the unit.
           subtitle = paste0("Matched-nucleus centroid residual, median per channel pair. ",
                             "Physical units; lower = tighter.\n", PAIR_UNIT_NOTE),
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
           subtitle = paste0("Higher = better. Same arms and channel pairs as the ",
                             "residual figure.\n", PAIR_UNIT_NOTE),
           x = NULL, y = "Matched-nucleus Dice (unitless, 0-1)", caption = ARM_CAPTION)
  }

  # -- 2b. VALIS grading itself, at each arm's final transform — the third main-text
  # panel, and the only one of the three that owes nothing to the segmentation.
  # Same unit as 1 and 2 (a channel pair), but NOT the same pairs: VALIS reports a
  # slide against the neighbour it was aligned toward (`from -> to`), the overlap QC
  # reports every moving slide against the reference. VALIS arms only, because the
  # other backends write no registered/summary/*.csv; STARE's pixel TRE stays in
  # figure 8, on its own axis.
  #
  # "Final" is the last VALIS stage the slide reported, `original` excluded. A single
  # ordering is safe here — unlike arm_final_stage() — because only one backend's
  # vocabulary is ever in this frame.
  vl <- if (nrow(valis)) valis_error_long(valis) else tibble::tibble()
  vf <- if (nrow(vl) && "arm" %in% names(vl))
    dplyr::filter(vl, is.finite(error)) |> dplyr::mutate(arm = .arm_f(arm))
  else tibble::tibble()
  vfin <- tibble::tibble()
  if (nrow(vf)) {
    vfin <- vf |>
      dplyr::filter(stage != "original") |>
      dplyr::group_by(arm, patient_id, slide) |>
      dplyr::slice_max(order_by = as.integer(stage), n = 1, with_ties = FALSE) |>
      dplyr::ungroup()
    for (k in setdiff(c("backend", "micro_reg"), names(vfin))) vfin[[k]] <- NA
    vfin$backend <- dplyr::coalesce(as.character(vfin$backend), "valis")
  }
  if (nrow(vfin)) {
    figs[["02b_final_valis_error_by_arm"]] <-
      ggplot(vfin, aes(arm, error)) +
      geom_boxplot(outlier.shape = NA, width = .5, colour = "grey35") +
      scale_x_discrete(labels = label_n(vfin$arm, sep = " ")) +
      geom_jitter(aes(colour = .arm_kind(backend, micro_reg)), width = .12, height = 0,
                  alpha = .85, size = 2.2) +
      scale_colour_arm() +
      coord_flip() +
      labs(title = "VALIS's own reported error at each arm's final transform",
           subtitle = paste0("Lower = better. VALIS arms only; independent of the two ",
                             "segmentation-overlap panels.\n", PAIR_UNIT_NOTE),
           x = NULL, y = vfin$metric[1] %||% "VALIS rTRE / distance",
           caption = ARM_CAPTION)
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
    # SYMMETRIC log, not log10. The delta is signed on purpose: negative means the
    # stage improved on the rigid anchor, positive means it regressed, and the zero
    # line between them is the whole reading. A plain scale_y_log10() would drop
    # every negative value — i.e. exactly the runs where micro-registration worked —
    # and leave a figure showing only failures. pseudo_log_trans is linear within
    # ±sigma of zero and logarithmic in both tails, so the small deltas near zero
    # stay separable and a 100 µm blow-up no longer sets the range for all of them.
    #
    # sigma is in µm: 0.1 µm is well below the ~0.5 µm pixel size these residuals
    # are measured at, so the linear window covers "indistinguishable from no
    # change" and nothing wider.
    delta_trans  <- scales::pseudo_log_trans(sigma = 0.1)
    delta_breaks <- c(-100, -10, -1, -0.1, 0, 0.1, 1, 10, 100)
    figs[["06_delta_vs_rigid_anchor"]] <-
      ggplot(dplyr::mutate(dv, arm = .arm_f(arm)),
             aes(stage, d_disp_um_vs_rigid, colour = arm)) +
      geom_hline(yintercept = 0, colour = REF_LINE) +
      geom_jitter(width = .15, height = 0, alpha = .8, size = 2) +
      scale_colour_manual(values = rep_len(oi_ext, dplyr::n_distinct(dv$arm)), name = NULL) +
      scale_y_continuous(trans = delta_trans, breaks = delta_breaks) +
      scale_x_discrete(labels = label_n(dv$stage)) +
      labs(title = "Change against the rigid anchor, per stage",
           # Explicit line break: the two sentences together overrun a double-column
           # panel, and a clipped subtitle is the one figure defect a reader cannot
           # work around.
           subtitle = paste0(
             paste("Residual minus the rigid stage's.",
                   "ABOVE the zero line means that stage made alignment WORSE",
                   "— the failure mode micro-registration hides."),
             "\n",
             paste("Axis is symmetric-log: linear within ±0.1 µm of zero (below the",
                   "pixel size), logarithmic in both tails.")),
           x = NULL, y = "Δ residual vs rigid (µm, symlog)", caption = ARM_CAPTION)
  }

  # -- 7. VALIS grading itself, ONE figure: the stage axis, faceted by arm.
  #
  # The companion of figure 2b, which summarises each arm at its final stage for the
  # main text. On its own that collapses an arm to one box and throws away the ladder,
  # which is the interesting part — and the ladder is what makes the arms comparable at
  # all, because the STAGE MEANINGS differ by depth and faceting is what fixes them
  # (same reason figure 3 facets).
  #
  # Three stages, two columns. VALIS's error_df is `from`/`filename`, `rigid_D`,
  # `non_rigid_D`; micro-registration has no column of its own because it UPDATES the
  # non-rigid field, so the micro value lives in the difference between the pre-micro and
  # final files rather than in a column. valis_error_long() does that reconstruction.
  #
  # The depth-0 and depth-1 arms therefore show NO micro box, and that blank is the
  # finding: micro-registration did not run. A duplicated non_rigid box would instead
  # read as "micro bought nothing".
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

  # -- S1-S6. THE SUPPLEMENTARY SPLITS of the three main panels. The main panels pool
  # every channel pair of every patient into one box per arm, which answers "which
  # arm" and hides WHY an arm is wide. The two splits separate the two possible
  # reasons: a hard patient (tissue, section quality) or a hard channel pair (two
  # rounds that share little structure). Same points as the main panels, regrouped —
  # nothing is re-summarised, so a point can be followed from one figure to the next.
  splits <- list(
    residual_um = list(d = dplyr::filter(fin, is.finite(disp_um_p50)), y = "disp_um_p50",
                       what = "Residual alignment error",
                       ylab = "residual displacement, median (µm)"),
    dice        = list(d = dplyr::filter(fin, is.finite(dice_matched)), y = "dice_matched",
                       what = "Matched-nucleus Dice",
                       ylab = "Matched-nucleus Dice (unitless, 0-1)"),
    valis_error = list(d = vfin, y = "error", what = "VALIS's own reported error",
                       ylab = if (nrow(vfin)) vfin$metric[1] else NA_character_))
  # The by-pair split only works if a pair is spelled the same in every patient. When
  # NO label is shared between patients the figure still draws — one point per box —
  # and looks like a finding. It is a naming problem (see .channel_pair()), so say so.
  if (dplyr::n_distinct(fin$patient_id) > 1 &&
      all(tapply(fin$patient_id, fin$pair, dplyr::n_distinct) == 1))
    warning("registration arms: no channel-pair label is shared between patients, so ",
            "the by-channel-pair figures hold one patient per box. The slide names ",
            "probably embed the patient id somewhere other than a leading `<patient>_`.")
  i <- 0
  for (nm in names(splits)) {
    sp <- splits[[nm]]
    for (by in c("patient", "channel_pair")) {
      i <- i + 1
      if (nrow(sp$d) == 0) next
      figs[[sprintf("S%d_%s_by_%s", i, nm, by)]] <-
        .arm_split_fig(dplyr::mutate(sp$d, arm = .arm_f(arm)), sp$y, by, sp$what, sp$ylab)
    }
  }

  figs
}

# One supplementary split: the final-transform metric `y`, one panel per arm, grouped
# by patient (each box = that patient's channel pairs) or by channel pair (each box =
# that pair across patients).
#
# n goes in the subtitle, not on the ticks: label_n() counts from the whole vector, so
# on a faceted axis it would print a patient's n pooled over every arm — a number that
# describes no box on the page.
.arm_split_fig <- function(d, y, by = c("patient", "channel_pair"), what, ylab) {
  by  <- match.arg(by)
  col <- if (by == "patient") "patient_id" else "pair"
  sub <- if (by == "patient")
    "Each box is ONE PATIENT, pooled across its channel pairs; one point per pair."
  else
    "Each box is ONE CHANNEL PAIR, pooled across patients; one point per patient."
  ggplot(d, aes(.data[[col]], .data[[y]])) +
    geom_boxplot(outlier.shape = NA, width = .5, colour = "grey35") +
    geom_jitter(aes(colour = .arm_kind(backend, micro_reg)), width = .12, height = 0,
                alpha = .85, size = 1.8) +
    scale_colour_arm() +
    coord_flip() +
    facet_wrap(~ arm) +
    labs(title = paste0(what, " at each arm's final transform, by ",
                        if (by == "patient") "patient" else "channel pair"),
         subtitle = paste0(sub, " ", n_note(d$patient_id, "patients"), ", ",
                           n_note(d$pair, "channel pairs"), "."),
         x = NULL, y = ylab, caption = ARM_CAPTION)
}

# How to colour an arm: by BACKEND, with the VALIS micro depth as the sub-level. A
# tiled arm has no depth, so colouring by depth alone would draw it as NA and read as
# missing data rather than as the other backend. ARM_KIND_COLS itself lives in
# plot_theme.R with the other recurring category palettes; use scale_colour_arm().

.arm_kind <- function(backend, micro_reg, baseline = FALSE) {
  # ASHLAR named explicitly: it has no micro depth either, so without its own branch
  # it fell through to the plain "valis" level and was drawn, and keyed, as VALIS.
  # A baseline row (Fig 4's no-registration / rigid-only) is not an arm of any
  # backend, so it is tested first and takes the one non-arm hue.
  baseline <- rep_len(baseline, length(backend))   # the default is a scalar FALSE
  lab <- ifelse(baseline, "baseline",
         ifelse(backend == "tiled", "tiled (STARE)",
         ifelse(backend == "ashlar", "ashlar",
                ifelse(is.na(micro_reg), "valis",
                       paste0("valis · micro ", micro_reg)))))
  factor(lab, levels = names(ARM_KIND_COLS))
}

# --- The manuscript's Fig 4(b)/(c) and Additional file 2 ---------------------
# The website's build_arm_figs() keeps every arm's ladder and every diagnostic. The
# manuscript asks a narrower question — which configuration, against what — and
# these two panels answer only that: one point per slide at each arm's FINAL
# transform, plus the two baselines the legend names.
#
# THE BASELINES ARE STAGES OF THE DEPTH-0 RUNS, NOT EXTRA ARMS. Nothing was run
# "without registration": `native` is the untransformed segmentation scored by the
# same QC, and `rigid` is the affine transform alone. Both are read from the
# `reg_micro_reg = 0` arms ONLY, because at depth >= 1 mirage's `rigid` stage is
# affine ∘ micro-rigid (see the header) and calling that "rigid only" would be the
# inverted-conclusion figure this file exists to prevent. The tiled arm is a
# different backend with a different `rigid` and contributes to neither baseline.
# Each baseline therefore holds one point per slide per depth-0 preset, and its n
# says so.
#
# TRE IS VALIS'S OWN NUMBER; DICE IS NOT. Panel (b) plots the error VALIS reports
# from its feature correspondences (`*_rTRE`, a fraction of the image diagonal, or
# `*_D` in physical units when that is what the run wrote) — the metric the legend
# ranks arms on. Panel (c) plots the matched-nucleus Dice from mirage's
# segmentation-overlap QC, which never sees those features, which is what makes it
# the independent corroboration. The tiled backend writes no VALIS summary, so it
# appears in (c) and not in (b).
BASELINE_NONE  <- "no registration"
BASELINE_RIGID <- "rigid only"

.arm_levels <- function(seg, manifest) {
  arms <- if (!is.null(manifest) && nrow(manifest)) unique(manifest$arm) else unique(seg$arm)
  # Baselines FIRST: coord_flip() draws the first level at the bottom, and the
  # reader expects the floor of the axis to be the floor of the ranking.
  c(BASELINE_NONE, BASELINE_RIGID, setdiff(arms, c(BASELINE_NONE, BASELINE_RIGID)))
}

# The two baseline rows from a staged frame. `stage_col` is `stage` for the
# segmentation QC and for VALIS's long form alike; `native_stage` names the
# no-registration stage in that vocabulary (`native` in the QC, `original` in VALIS).
.arm_baselines <- function(staged, native_stage) {
  if (nrow(staged) == 0) return(staged)
  d0 <- dplyr::filter(staged, backend == "valis", !is.na(micro_reg), micro_reg == 0L)
  dplyr::bind_rows(
    dplyr::filter(d0, as.character(stage) == native_stage) |>
      dplyr::mutate(arm = BASELINE_NONE),
    dplyr::filter(d0, as.character(stage) == "rigid") |>
      dplyr::mutate(arm = BASELINE_RIGID)) |>
    dplyr::mutate(is_baseline = TRUE) |>
    dplyr::rename(final_stage = stage)
}

# VALIS's error at each run's final stage, one row per (arm, slide). The final stage
# is `micro` when a pre-micro summary proves micro ran, else `non_rigid` — the same
# reconstruction valis_error_long() does for the website ladder.
.valis_final <- function(vl) {
  if (nrow(vl) == 0) return(vl)
  vl |>
    dplyr::filter(as.character(stage) != "original") |>
    dplyr::group_by(arm, arm_dir, backend, memory_mode, micro_reg, patient_id, slide) |>
    dplyr::slice_max(order_by = as.integer(stage), n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::rename(final_stage = stage) |>
    dplyr::mutate(is_baseline = FALSE)
}

.arm_paper_panel <- function(d, y, y_lab, title, subtitle, lvls, log_y = FALSE) {
  d <- d |>
    dplyr::filter(is.finite(.data[[y]])) |>
    dplyr::mutate(arm  = factor(arm, levels = intersect(lvls, unique(arm))),
                  kind = .arm_kind(backend, micro_reg, is_baseline))
  if (nrow(d) == 0) return(NULL)
  p <- ggplot(d, aes(arm, .data[[y]])) +
    geom_boxplot(outlier.shape = NA, width = .5, colour = "grey35") +
    scale_x_discrete(labels = label_n(d$arm, sep = " ")) +
    geom_jitter(aes(colour = kind), width = .12, height = 0, alpha = .85, size = 2.2) +
    scale_colour_arm() +
    coord_flip() +
    labs(title = title, subtitle = subtitle, x = NULL, y = y_lab, caption = ARM_CAPTION)
  # The baselines sit one to two orders of magnitude above the arms. On a linear
  # axis every arm collapses onto zero and the ranking — the point of the panel —
  # is unreadable; log keeps the ratios legible without hiding the baselines.
  if (log_y) p <- p + scale_y_log10()
  p
}

#
# THE SAME ARMS ON BOTH PANELS. The legend says (c) is "plotted for the same arms
# as (b)", and (b) can only hold VALIS arms (the tiled backend reports no VALIS
# error), so by default both panels hold the VALIS preset x depth grid and nothing
# else. That also lets figures/fig4.R drop (c)'s arm labels and read the rows off
# (b). The STARE comparator stays on the website page and in Additional file 2;
# pass backends = c("valis", "tiled") to draw it on (c) — and then keep (c)'s labels.
build_arm_paper_figs <- function(seg = read_arms_seg_qc(), valis = read_arms_valis(),
                                 manifest = NULL, backends = "valis") {
  figs <- list()
  if (nrow(seg) == 0) return(figs)
  seg <- dplyr::filter(seg, backend %in% backends)
  if (nrow(seg) == 0) return(figs)
  if (!is.null(manifest) && nrow(manifest))
    manifest <- dplyr::filter(manifest, backend %in% backends)
  lvls <- .arm_levels(seg, manifest)

  # -- (b) VALIS-internal TRE, by arm, with baselines.
  vl <- if (nrow(valis)) valis_error_long(valis) else tibble::tibble()
  if (nrow(vl) && "arm" %in% names(vl)) {
    d <- dplyr::bind_rows(.valis_final(vl), .arm_baselines(vl, "original"))
    figs[["tre_by_arm"]] <- .arm_paper_panel(
      d, "error", vl$metric[1] %||% "VALIS rTRE",
      title    = "Registration error reported by VALIS, at each arm's final transform",
      subtitle = paste("One point per moving slide. Baselines are the untransformed",
                       "and affine-only stages of the depth-0 runs."),
      lvls = lvls, log_y = TRUE)
  }

  # -- (c) matched-nucleus Dice, by arm, with baselines.
  fin <- arm_final_stage(seg) |> dplyr::mutate(is_baseline = FALSE)
  d   <- dplyr::bind_rows(fin, .arm_baselines(seg, "native"))
  figs[["dice_by_arm"]] <- .arm_paper_panel(
    d, "dice_matched", "Matched-nucleus Dice (unitless, 0-1)",
    title    = "DAPI-nucleus overlap at each arm's final transform",
    subtitle = paste("Independent of the features VALIS registered on.",
                     "Same arms and slides as the TRE panel; higher = better."),
    lvls = lvls)

  figs[!vapply(figs, is.null, logical(1))]
}

# Additional file 2: one row per arm plus the two baselines, both metrics side by
# side. Medians over slides; `n_slides` counts what each median is over, which for
# a baseline is slides x depth-0 presets.
arm_paper_table <- function(seg = read_arms_seg_qc(), valis = read_arms_valis(),
                            manifest = NULL) {
  if (nrow(seg) == 0) return(tibble::tibble())
  lvls <- .arm_levels(seg, manifest)
  keys <- c("arm", "backend", "memory_mode", "micro_reg")

  qc <- dplyr::bind_rows(arm_final_stage(seg) |> dplyr::mutate(is_baseline = FALSE),
                         .arm_baselines(seg, "native")) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(
      n_slides      = dplyr::n(),
      stage         = paste(sort(unique(as.character(final_stage))), collapse = "/"),
      disp_um_p50   = stats::median(disp_um_p50, na.rm = TRUE),
      dice_matched  = stats::median(dice_matched, na.rm = TRUE),
      pair_fraction = stats::median(pair_fraction, na.rm = TRUE),
      .groups = "drop")

  vl <- if (nrow(valis)) valis_error_long(valis) else tibble::tibble()
  tre <- if (nrow(vl) && "arm" %in% names(vl)) {
    dplyr::bind_rows(.valis_final(vl), .arm_baselines(vl, "original")) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
      dplyr::summarise(valis_tre = stats::median(error, na.rm = TRUE),
                       valis_tre_metric = metric[1], .groups = "drop")
  } else {
    tibble::tibble(arm = character(), backend = character(), memory_mode = character(),
                   micro_reg = integer(), valis_tre = numeric(), valis_tre_metric = character())
  }

  # Baseline rows carry no preset or depth of their own: they pool the depth-0 runs.
  collapse <- function(x) dplyr::mutate(x,
    memory_mode = ifelse(arm %in% c(BASELINE_NONE, BASELINE_RIGID), NA_character_, memory_mode),
    micro_reg   = ifelse(arm %in% c(BASELINE_NONE, BASELINE_RIGID), NA_integer_, micro_reg))
  qc  <- collapse(qc)  |> dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(n_slides = sum(n_slides), stage = paste(unique(stage), collapse = "/"),
                     dplyr::across(c(disp_um_p50, dice_matched, pair_fraction),
                                   ~ stats::median(.x, na.rm = TRUE)), .groups = "drop")
  tre <- collapse(tre) |> dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(valis_tre = stats::median(valis_tre, na.rm = TRUE),
                     valis_tre_metric = valis_tre_metric[1], .groups = "drop")

  dplyr::left_join(qc, tre, by = keys) |>
    dplyr::select(arm, backend, memory_mode, micro_reg, n_slides, stage,
                  valis_tre, valis_tre_metric, disp_um_p50, dice_matched, pair_fraction) |>
    dplyr::mutate(arm = factor(arm, levels = intersect(lvls, unique(arm)))) |>
    dplyr::arrange(arm) |>
    dplyr::mutate(arm = as.character(arm))
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
