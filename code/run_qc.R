# =============================================================================
# run_qc.R  —  mirage's QC of the run on THE STUDY SAMPLES
#
# The benchmark pages (benchmarks.Rmd, benchmark_registration.Rmd) read the tables
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
#   registered/summary/**/*.csv      VALIS's OWN error, written during register().
#                                    Columns are whatever VALIS declared (rTRE / D /
#                                    n_matches variants) — read generically. Note the
#                                    RECURSION: Nextflow preserves REGISTER's
#                                    `preprocessed/data/` producer subdirectory, so
#                                    these actually sit two levels down. See
#                                    .qc_files_in() for why that is not cosmetic.
#   qc/registration/*_tre.json       STARE's own TRE, when the tiled path ran.
#                                    coarse/rigid/post-refinement percentiles PLUS a
#                                    per-tile spatial map VALIS cannot give.
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
# warp_seg_qc stages, BOTH backends. lib/WarpBackends.groovy declares one stage list
# per registration method and they are not the same vocabulary:
#   valis -> native, rigid, non_rigid, micro
#   tiled -> native, rigid, refined          (STARE: mesh refinement of the anchor)
# Only `native` and `rigid` are shared spellings, and `rigid` does NOT mean the same
# operation in both (VALIS: affine, composed with micro-rigid at reg_micro_reg >= 1;
# STARE: the coarse global anchor before mesh refinement). This vector exists to stop
# a stage being DROPPED, which is what a VALIS-only list did to every tiled run's
# final stage — see stage_index below for how ordering is actually decided.
QC_STAGE_LEVELS     <- c("native", "rigid", "refined", "non_rigid", "micro")
# VALIS's OWN error, as a stage axis. THREE stages from TWO columns.
#
# VALIS's error_df schema is `from`/`filename`, `rigid_D`, `non_rigid_D` — that is the
# whole of it. Micro-registration is not modelled as a fourth stage with a column of
# its own; it is an UPDATE TO THE NON-RIGID FIELD. So the micro number does not live
# in a column at all: it lives in the difference between the two files.
# mirage's own report recovers the axis the same way, from the FILE rather than the
# column (bin/generate_qc_report.py:_RECONCILE_TRE_SOURCE).
#
# `original` is kept in the vocabulary but is normally ABSENT: VALIS's per-patient
# summary has no such column. The benchmark sweep's registration_valis_rtre.csv does
# (see registration_accuracy_plots.R) — a different artifact, from make_tables.py — and
# a build that emits one here should plot rather than be silently dropped, which is why
# the level survives and the column is detected rather than named.
VALIS_STAGE_LEVELS  <- c("original", "rigid", "non_rigid", "micro")

# --- helpers ----------------------------------------------------------------
# Where each artifact sits inside a patient directory. Also the definition of "is a
# patient": a directory carrying ANY of them. Detecting on one of them only (say a
# qc/ subdirectory) hides a run that produced a different subset — a VALIS run with
# reg_qc < 2 and no phenotyping has registered/summary/ and nothing else, and would
# have rendered an empty page rather than its rTRE.
QC_ARTIFACTS <- c(seg_qc = "qc/registration",
                  stare  = "qc/registration",
                  valis  = "registered/summary")

# SEARCH THE ARTIFACT DIRECTORIES RECURSIVELY, NOT ONE LEVEL DEEP.
# Nextflow's publishDir PRESERVES the producer subdirectory: a process that writes into
# a named subdirectory of its task dir and publishes with a `pattern:` naming that
# subdirectory gets the subdirectory carried into the published path. mirage's
# lib/Layout.groovy states this and names the deepest case — REGISTER declares its VALIS
# error tables as `path("preprocessed/data/*.csv")` (modules/local/register.nf) and
# publishes that pattern into `<pid>/registered/summary`, so they land TWO levels down at
#
#     <pid>/registered/summary/preprocessed/data/<name>_summary.csv
#
# not in `registered/summary/` itself. A one-level `dir_ls` therefore found no VALIS
# summaries on any real run, and — because .qc_patient_dirs() detects a patient by
# finding FILES — a VALIS run with reg_qc < 2 and no phenotyping was not detected as a
# patient at all: `registered/summary/` holds only a `preprocessed/` DIRECTORY, so the
# page reported "No mirage QC found" for a run whose rTRE was sitting right there.
# Recursing is also what keeps this robust to the next producer subdirectory, since
# nothing stops a future emit from being deeper still.
.qc_files_in <- function(d, glob = NULL) {
  if (!fs::dir_exists(d)) return(character(0))
  as.character(fs::dir_ls(d, recurse = TRUE, type = "file", glob = glob))
}

# Patient directories under `root`. A mirage outdir also carries cohort-level
# directories (qc/, phenotyping/, size_logs/, _UNROUTED_PUBLISH/) beside the patient
# ones; those hold no per-patient artifact and so never match.
#
# A SYMLINKED PATIENT DIRECTORY IS STILL A PATIENT DIRECTORY. `dir_ls(type =
# "directory")` filters on the entry's OWN type, and a symlink's type is "symlink" —
# so pulling one run's patients together under data/mirage/ by symlinking each of them
# (rather than symlinking data/mirage itself) yielded zero patients and an empty page.
# List everything and keep what resolves to a directory, which follows the link.
.qc_patient_dirs <- function(root) {
  if (!fs::dir_exists(root)) return(character(0))
  dirs <- as.character(fs::dir_ls(root))
  dirs <- dirs[fs::dir_exists(dirs)]
  has_any <- vapply(dirs, function(d)
    any(vapply(unique(QC_ARTIFACTS), function(sub)
      length(.qc_files_in(file.path(d, sub))) > 0, logical(1))), logical(1))
  dirs[has_any]
}

.read_json <- function(path) tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE),
                                      error = function(e) NULL)

# Every file matching `glob` anywhere under one patient's subdirectory.
.qc_files <- function(dir, sub, glob) .qc_files_in(file.path(dir, sub), glob)

# The join key for "the same moving slide", across artifacts that spell it differently.
#
# `moving` is written by three producers and they do not agree. On the VALIS path they
# do, by a documented invariant: seg_qc names the slide from the native image
# (mirage's seg_qc.nf resolves `meta.qc_slide` by stripping .ome.tif/.ome.tiff, and its
# comment requires it to equal the registrar's slide_dict key), and VALIS's own summary
# is keyed by that same name. On the TILED path they differ by exactly the patient
# prefix: TILED_SOLVE's `_tre.json` is named `<patient_id>_<channels>_tre.json` while
# seg_qc carries the native stem, which may or may not carry the patient id.
#
# So normalise rather than guess: drop any image extension and a leading `<patient>_`.
# This is a strict WIDENING of the old exact match — anything that joined before still
# joins — and it is what stops §4 from silently plotting nothing when the two spellings
# differ only by that prefix.
.slide_token <- function(x, patient_id = NULL) {
  s <- basename(as.character(x))
  s <- sub("\\.(ome\\.tiff?|tiff?|qptiff|svs|ndpi)$", "", s, ignore.case = TRUE)
  if (is.null(patient_id) || !length(patient_id)) return(s)
  pre <- paste0(as.character(patient_id), "_")   # recycled against s, one pid per row
  ifelse(startsWith(s, pre), substring(s, nchar(pre) + 1L), s)
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
  patient_id = character(), moving = character(), slide_token = character(),
  stage = factor(levels = QC_STAGE_LEVELS),
  stage_index = integer(),
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
      mv     <- d$moving %||% fs::path_ext_remove(fs::path_file(f))
      purrr::map_dfr(seq_along(stages), function(i) {
        st <- stages[[i]]
        s  <- d$stages[[st]]
        dv <- (d$delta_vs_anchor %||% list())[[st]] %||% list()
        tibble::tibble(
          patient_id       = slide_key(d$patient_id %||% fs::path_file(dir)),
          moving           = mv,
          slide_token      = .slide_token(mv, fs::path_file(dir)),
          stage            = st,
          # This run's OWN position for the stage, straight from its `stage_order`.
          # "Which stage did this run end on" must not be answered by a global factor
          # ordering: the two warp backends have different vocabularies, so any single
          # ordering would silently rank one backend's stages against the other's.
          # Read off the producer's own list, it is right for both by construction.
          stage_index      = i,
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
        # The RAW directory name, before slide_key() normalises it: it is the
        # `meta.patient_id` mirage prefixes slide names with, so it is what
        # .slide_token() has to strip. slide_key("EPM - 052") is "052", which is not.
        patient_dir = fs::path_file(dir),
        summary_csv = nm,
        # A `_summary_premicro.csv` is written ONLY at reg_micro_reg = 2
        # (bin/register.py: "below it register_micro never runs, so the summary already
        # IS the pre-micro one"). So the plain summary is post-micro only when a
        # premicro sibling proves micro ran; on its own it is simply the final — and
        # calling it "post-micro" told the reader a stage happened that did not.
        stage_scope = if (grepl("premicro", nm)) "pre-micro" else "final",
        .before = 1)
    })
  })
}

# VALIS's per-stage error in long form, on ONE semantic stage axis.
#
# WHY THERE IS NO MICRO COLUMN. VALIS's error_df carries `rigid_D` and `non_rigid_D`
# and nothing else. Micro-registration is not a fourth stage with its own column — it
# is an update to the non-rigid displacement field: `register_micro()` re-runs
# `measure_error()` and OVERWRITES `<name>_summary.csv`, composing the micro residual
# into that same field. So the micro number does not live in a column. **It lives in
# the difference between the two files.**
#
# The stage axis is therefore recovered from WHICH FILE a value came from, not from
# which column — exactly as mirage's own report does
# (bin/generate_qc_report.py:_RECONCILE_TRE_SOURCE):
#
#   rigid      rigid_D      final       (unchanged by micro; read from either file)
#   non_rigid  non_rigid_D  PRE-MICRO   (the only file that still isolates it)
#   micro      non_rigid_D  FINAL       (same column, other file)
#
# and `_read_intrinsic_tre()` sorts the two files purely on the `_summary_premicro.csv`
# filename suffix.
#
# This is why the old pre-/post-micro FACET looked like a duplicate: `rigid_D` is the
# same number in both files, so one of two stages was drawn twice, and the one that did
# differ was labelled as the same stage in both panels when the two values are two
# different stages.
#
# A BLANK MICRO STAGE IS A CLAIM, NOT A GAP. With no pre-micro file (reg_micro_reg < 2,
# the shipped default) the final summary's non_rigid_D IS the pre-micro value, so it
# becomes `non_rigid` and NO micro stage is emitted. Emitting one anyway — a
# byte-for-byte duplicate of non_rigid — would read as "micro bought nothing", which is
# a different statement from "micro did not run". mirage makes the same choice on the
# segmentation-overlap side (docs/registration_qc.md: when micro-registration is skipped
# or raises and is caught, the QC omits the stage rather than duplicating non_rigid).
#
# STARE never appears here at all: the tiled backend writes no registered/summary/*.csv,
# so it has no micro stage by construction. Its intrinsic TRE is read from *_tre.json
# and plotted separately, in pixels (see registration_arms.R). mirage's report folds both
# into one slide dict by reusing the rigid_D / non_rigid_D keys; this project keeps them
# apart because the units differ and a shared axis invites a comparison that is not one.
valis_error_long <- function(valis) {
  if (nrow(valis) == 0) return(tibble::tibble())
  cols <- grep("_rTRE$", names(valis), value = TRUE)
  metric <- "VALIS rTRE (relative)"
  if (!length(cols)) {
    cols <- grep("_D$", names(valis), value = TRUE)
    metric <- "VALIS matched-feature distance (D)"
  }
  cols <- cols[sub("_(rTRE|D)$", "", cols) %in%
                 c("original", "rigid", "non_rigid")]   # the columns that exist
  if (!length(cols)) return(tibble::tibble())
  id <- intersect(c("img_name", "name", "filename"), names(valis))[1]

  # Columns a caller may have attached before handing the frame over (the arm sweep
  # adds five). Carried through by name rather than dropped by a fixed select: this
  # function returning a narrower frame than it was given is how an arm-aware caller
  # loses its grouping and its figure quietly stops being built.
  extra <- intersect(c("arm", "arm_dir", "backend", "memory_mode", "micro_reg"),
                     names(valis))

  long <- valis |>
    dplyr::mutate(slide = if (is.na(id)) summary_csv else .data[[id]],
                  slide_token = .slide_token(slide, patient_dir)) |>
    dplyr::select(patient_id, slide, slide_token, stage_scope,
                  dplyr::all_of(extra), dplyr::all_of(cols)) |>
    tidyr::pivot_longer(dplyr::all_of(cols), names_to = ".col", values_to = "error") |>
    dplyr::mutate(.col = sub("_(rTRE|D)$", "", .col)) |>
    dplyr::filter(is.finite(error))
  if (nrow(long) == 0) return(tibble::tibble())

  # Did micro actually run for this slide? Only a pre-micro sibling proves it.
  # Grouped by the caller's extra keys too: the same slide appears once per arm in an
  # arm sweep, and pooling them would let one arm's pre-micro summary decide another
  # arm's stage labels.
  grp <- c("patient_id", "slide", extra)
  micro_ran <- long |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) |>
    dplyr::summarise(has_premicro = any(stage_scope == "pre-micro"), .groups = "drop")

  long |>
    dplyr::left_join(micro_ran, by = grp) |>
    dplyr::mutate(stage = dplyr::case_when(
      # The invariant stages are read from the final summary only; the pre-micro copy
      # of them is the same number and would double every box.
      .col %in% c("original", "rigid") & stage_scope == "final" ~ .col,
      .col %in% c("original", "rigid")                         ~ NA_character_,
      .col == "non_rigid" & stage_scope == "pre-micro"          ~ "non_rigid",
      # Final non_rigid: the micro stage when micro ran, otherwise the non-rigid stage.
      .col == "non_rigid" & has_premicro                        ~ "micro",
      .col == "non_rigid"                                       ~ "non_rigid",
      TRUE                                                      ~ NA_character_)) |>
    dplyr::filter(!is.na(stage)) |>
    dplyr::mutate(stage  = factor(stage, levels = VALIS_STAGE_LEVELS),
                  metric = metric) |>
    dplyr::select(patient_id, slide, slide_token, dplyr::all_of(extra),
                  stage, error, metric, source_file = stage_scope)
}

# --- 3. STARE's own TRE (tiled path only) ------------------------------------
# Per patient/moving summary. Absent when the run used the VALIS path.
read_stare_tre <- function(root = RUN_QC_ROOT) {
  purrr::map_dfr(.qc_patient_dirs(root), function(dir) {
    purrr::map_dfr(.qc_files(dir, "qc/registration", "*_tre.json"), function(f) {
      d <- .read_json(f)
      if (is.null(d) || is.null(d$coarse_tre_px)) return(tibble::tibble())
      pct <- function(x, k) as.numeric((d[[x]] %||% list())[[k]] %||% NA)
      # From the FILENAME (`<patient_id>_<channels>_tre.json`), not from the JSON's own
      # `moving` field: TILED_SOLVE sets that to the channel set alone, which carries
      # less than the filename does. .slide_token() then reconciles the two spellings.
      mv <- sub("_tre$", "", fs::path_ext_remove(fs::path_file(f)))
      tibble::tibble(
        patient_id     = slide_key(fs::path_file(dir)),
        moving         = mv,
        slide_token    = .slide_token(mv, fs::path_file(dir)),
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

# Everything at once, for the report.
run_qc_tables <- function(root = RUN_QC_ROOT) {
  list(seg_qc     = read_seg_qc(root),
       valis      = read_valis_summary(root),
       stare      = read_stare_tre(root),
       stare_tiles= read_stare_tiles(root))
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
    # Boxplot per stage with the median PRINTED. Boxes only: the per-slide points and
    # their patient legend were dropped because they were what the figure was read
    # for, and the box already summarises them — and reading a value off a log axis
    # is guesswork, so the number is written on the box. Outliers stay drawn, since
    # without the overlay they would otherwise vanish from the figure entirely.
    med <- vl |>
      dplyr::group_by(stage) |>
      dplyr::summarise(error = stats::median(error, na.rm = TRUE),
                       n = dplyr::n(), .groups = "drop")
    micro_ran <- "micro" %in% as.character(vl$stage)
    # Log y axis. The error falls by an order of magnitude or more from `original`
    # to `non_rigid`, and a slide or two register badly enough to set a linear range
    # for every stage — on a linear axis the registered boxes sat on zero and the
    # figure used to clip its top to stay readable. A log axis shows every stage
    # and every slide at once, so nothing is cropped and the caption has nothing
    # to declare. The medians were computed above, over every slide, BEFORE the
    # positivity filter below; that filter exists only because log10(0) is not a
    # coordinate, and an error of exactly 0 does not occur in practice.
    vl <- dplyr::filter(vl, error > 0)
    figs[["01_valis_error_by_stage"]] <-
      ggplot(vl, aes(stage, error)) +
      geom_boxplot(width = .5, colour = "grey35", outlier.size = 1.2) +
      scale_x_discrete(labels = label_n(vl$stage)) +
      geom_text(data = med, aes(label = signif(error, 3)),
                vjust = -1.1, size = pt_text(7), colour = "grey15") +
      scale_y_log10() +
      lab(title = "VALIS registration error by stage",
           subtitle = paste0(
             "VALIS grading itself from its own feature matches; lower = better. ",
             "Label = median. VALIS has no micro column — micro-registration updates the ",
             "non-rigid field, so the stage axis comes from which FILE a value is in. ",
             if (micro_ran)
               paste("`non_rigid` is the pre-micro summary; `micro` is the same",
                     "`non_rigid_D` column in the final one.")
             else
               paste("No `micro` box: this run wrote no pre-micro summary, so",
                     "reg_micro_reg < 2 and micro-registration never ran. A blank stage",
                     "means it did not run, not that it bought nothing.")),
           x = NULL, y = paste(vl$metric[1], "(log10)"))
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
    # Boxes only. The per-(patient, moving) connecting lines and their points made
    # this panel unreadable as soon as a run carried more than a couple of moving
    # slides: the spaghetti crossed every box and the dots sat on the median line
    # that is the actual reading. The per-slide ladder is not lost — that is what
    # 01_valis_error_by_stage and §3c are for.
    figs[["03_overlap_dice_by_stage"]] <-
      ggplot(sq, aes(stage, dice_matched)) +
      geom_boxplot(outlier.shape = NA, width = .5) +
      scale_x_discrete(labels = label_n(sq$stage)) +
      lab(title = "Matched-nucleus Dice by registration stage",
          subtitle = paste("Computed from DAPI segmentation overlap, NOT from the features the",
                           "registration used — the independent check. Higher = better."),
          x = NULL, y = "Matched-nucleus Dice (unitless, 0-1)")

    if (any(is.finite(sq$disp_um_p50)))
      figs[["03b_displacement_um_by_stage"]] <-
        sq |>
        tidyr::pivot_longer(c(disp_um_p50, disp_um_p90), names_to = "pct", values_to = "um") |>
        dplyr::mutate(pct = dplyr::recode(pct, disp_um_p50 = "median", disp_um_p90 = "90th pct")) |>
        # log10 below needs a strictly positive residual; filtering here states that
        # rather than letting scale_y_log10() drop rows with a knit-log warning.
        dplyr::filter(is.finite(um), um > 0) |>
        ggplot(aes(stage, um, colour = pct)) +
        geom_boxplot(outlier.shape = NA, width = .5, position = position_dodge(.6)) +
        # Residual spans orders of magnitude across the stage ladder — native is
        # tens of um, micro is sub-um — so on a linear axis every registered stage
        # collapses onto zero and the improvement that matters is invisible.
        scale_y_log10() +
        # Counted from `sq` (one row per run), NOT from the piped long frame: that
        # frame is pivoted over median/90th and would state twice the runs per box.
        # ... and over the rows the log axis actually draws, so the stated n does not
        # include a run the `um > 0` filter above removed.
        scale_x_discrete(labels = label_n(
          sq$stage[is.finite(sq$disp_um_p50) & sq$disp_um_p50 > 0])) +
        scale_colour_manual(values = oi[c(1, 2)], name = NULL) +
        lab(title = "Centroid residual displacement by stage",
            subtitle = "Physical units, so this is the number to quote as spatial resolution.",
            x = NULL, y = "displacement (µm, log10)")

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
            x = NULL, y = "Pair fraction (unitless, 0-1)")
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
    # Joined on `slide_token`, NOT on `moving`: the two sides name the same slide
    # differently on the tiled path (see .slide_token()), and an exact join on the raw
    # spelling drops every row — silently, since an empty join just removes the figure.
    intrinsic <- if (nrow(st))
      dplyr::transmute(st, patient_id, slide_token,
                       intrinsic = dplyr::coalesce(after_p50, rigid_p50),
                       what = "STARE TRE (px)")
    else if (nrow(vl))
      vl |>
        dplyr::filter(stage == dplyr::last(levels(droplevels(stage)))) |>
        dplyr::transmute(patient_id, slide_token, intrinsic = error, what = metric)
    else NULL
    if (!is.null(intrinsic) && nrow(intrinsic)) {
      ag <- dplyr::inner_join(final, intrinsic, by = c("patient_id", "slide_token")) |>
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
              x = ag$what[1], y = "Matched-nucleus Dice (unitless, 0-1)")
    }
  }

  figs
}
