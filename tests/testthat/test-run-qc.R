# mirage's QC of a run on the study samples. These build a synthetic outdir rather
# than needing a pipeline run — the schemas are pinned in code/run_qc.R's header.
source(here::here("code", "run_qc.R"))

qc_fixture <- function(patients = list(`046` = c("cycle2", "cycle3")),
                       tiled = FALSE) {
  root <- file.path(tempdir(), paste0("runqc-", as.integer(Sys.time()), "-", sample(1e6, 1)))
  for (pid in names(patients)) {
    d <- file.path(root, pid)
    dir.create(file.path(d, "qc", "registration"), recursive = TRUE, showWarnings = FALSE)
    for (mv in patients[[pid]]) {
      stages <- c("native", "rigid", "non_rigid", "micro")
      dice   <- c(.08, .55, .71, .73)
      disp   <- c(12, 3.2, 1.4, 1.3)
      jsonlite::write_json(list(
        patient_id = pid, moving = mv, stage_order = stages,
        matching = list(pair_fraction = 0.62),
        stages = stats::setNames(lapply(seq_along(stages), function(i) list(
          n_pairs = 4000L, iou_mean = dice[i] * .9, iou_p50 = dice[i] * .92,
          dice_matched = dice[i], displacement_px_p50 = disp[i] / 0.325,
          displacement_um_p50 = disp[i], displacement_um_p90 = disp[i] * 2)), stages),
        delta_vs_anchor = list(micro = list(dice_matched = .18,
                                            displacement_um_p50 = -1.9))),
        file.path(d, "qc", "registration", paste0(mv, "_seg_qc.json")), auto_unbox = TRUE)
      if (tiled) {
        tiles <- expand.grid(ix = 0:2, iy = 0:1)
        tiles$cx <- tiles$ix * 512; tiles$cy <- tiles$iy * 512
        tiles$tre_rigid <- 1 + tiles$ix; tiles$tre_after <- tiles$tre_rigid * .3
        jsonlite::write_json(list(
          coarse_tre_px = 4.2, n_inliers = 800L, n_tiles = nrow(tiles), mesh_refined = TRUE,
          rigid_tre_px = list(mean = 2, p50 = 2, p90 = 3, max = 3),
          residual_after_px = list(mean = .6, p50 = .6, p90 = .9, max = .9),
          tiles = tiles),
          file.path(d, "qc", "registration", paste0(mv, "_tre.json")), auto_unbox = TRUE)
      }
    }
    if (!tiled) {
      # NOT `registered/summary/` itself. REGISTER declares these as
      # path("preprocessed/data/*.csv") and Nextflow's publishDir preserves the
      # producer subdirectory, so a real run puts them two levels down. Writing them
      # flat here is what let a one-level dir_ls() pass its tests while finding
      # nothing on every actual pipeline output.
      sdir <- file.path(d, "registered", "summary", "preprocessed", "data")
      dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
      mv <- patients[[pid]]
      # post-micro AND pre-micro: register_micro() overwrites the plain summary, so
      # the two files are different measurements of the same slide
      for (suffix in c("_summary.csv", "_summary_premicro.csv"))
        readr::write_csv(tibble::tibble(
          img_name = mv, original_rTRE = 50, rigid_rTRE = 10,
          non_rigid_rTRE = if (suffix == "_summary.csv") 3 else 5,
          n_matches = 120L),
          file.path(sdir, paste0(pid, suffix)))
    }
  }
  # cohort-level directories a real outdir also carries
  dir.create(file.path(root, "qc"), showWarnings = FALSE)
  dir.create(file.path(root, "size_logs"), showWarnings = FALSE)
  root
}

test_that("staged overlap QC is read one row per (patient, moving, stage)", {
  sq <- read_seg_qc(qc_fixture(list(`046` = c("cycle2", "cycle3"), `052` = "cycle2")))
  expect_equal(nrow(sq), 3 * 4)
  expect_setequal(levels(sq$stage), QC_STAGE_LEVELS)
  final <- dplyr::filter(sq, stage == "micro")
  expect_equal(unique(final$dice_matched), .73)
  expect_equal(unique(final$disp_um_p50), 1.3)
  expect_equal(unique(final$d_disp_um_vs_rigid), -1.9)   # negative = tightened
})

test_that("an empty or absent tree yields typed empties, not an error", {
  # the first-knit path: `stage` would otherwise resolve to stats::stage()
  empty <- file.path(tempdir(), paste0("norunqc-", sample(1e6, 1)))
  dir.create(empty, showWarnings = FALSE)
  sq <- read_seg_qc(empty)
  expect_equal(nrow(sq), 0)
  expect_true(all(c("patient_id", "stage", "dice_matched") %in% names(sq)))
  expect_equal(sum(vapply(run_qc_tables(empty), nrow, integer(1))), 0)
  expect_length(build_run_qc_figs(empty, run_qc_tables(empty)), 0)
})

test_that("cohort-level directories are not read as patients", {
  sq <- read_seg_qc(qc_fixture())
  expect_setequal(unique(sq$patient_id), "046")
})

test_that("the two VALIS summaries resolve to ONE stage axis, not two panels", {
  # register_micro() re-runs measure_error() and overwrites <name>_summary.csv,
  # composing the micro residual into the same field. Across the two files:
  #   original / rigid  are IDENTICAL (micro cannot change them)
  #   non_rigid         differs, and the two values are two DIFFERENT STAGES —
  #                     pre-micro is the non-rigid stage, final is the micro stage.
  # Faceting by source file therefore duplicated two stages and mislabelled the third,
  # which is why the panels looked the same. mirage's own report makes exactly these
  # assignments (bin/generate_qc_report.py:_RECONCILE_TRE_SOURCE).
  v <- read_valis_summary(qc_fixture())
  expect_setequal(v$stage_scope, c("final", "pre-micro"))
  long <- valis_error_long(v)

  expect_setequal(unique(as.character(long$stage)),
                  c("original", "rigid", "non_rigid", "micro"))
  # The invariant stages appear ONCE, from the final summary — not twice.
  expect_equal(sum(long$stage == "rigid"), 2)          # two moving slides, one row each
  expect_equal(unique(long$source_file[long$stage == "rigid"]), "final")
  # The pre-micro non_rigid column is the non-rigid stage...
  expect_equal(unique(long$error[long$stage == "non_rigid"]), 5)
  expect_equal(unique(long$source_file[long$stage == "non_rigid"]), "pre-micro")
  # ...and the final one is the micro stage.
  expect_equal(unique(long$error[long$stage == "micro"]), 3)
  expect_equal(unique(long$source_file[long$stage == "micro"]), "final")
  expect_true(grepl("rTRE", long$metric[1]))
})

test_that("a lone summary is not called post-micro, and yields no micro stage", {
  # A `_summary_premicro.csv` is written ONLY at reg_micro_reg = 2. Below that the plain
  # summary already IS the pre-micro one, so labelling it "post-micro" claimed a stage
  # that never ran, and inventing an empty `micro` level would imply the same.
  root <- file.path(tempdir(), paste0("nomicro-", sample(1e6, 1)))
  d <- file.path(root, "046", "registered", "summary")
  dir.create(d, recursive = TRUE)
  readr::write_csv(tibble::tibble(img_name = "c2", original_rTRE = 50, rigid_rTRE = 10,
                                  non_rigid_rTRE = 3, n_matches = 99L),
                   file.path(d, "046_summary.csv"))
  v <- read_valis_summary(root)
  expect_equal(unique(v$stage_scope), "final")
  expect_false("post-micro" %in% v$stage_scope)

  long <- valis_error_long(v)
  expect_setequal(unique(as.character(long$stage)), c("original", "rigid", "non_rigid"))
  expect_false("micro" %in% as.character(long$stage))
  # The final non_rigid becomes the NON-RIGID stage here, not the micro one.
  expect_equal(long$error[long$stage == "non_rigid"], 3)
})

test_that("the stage figure is a bare boxplot on a log axis", {
  figs <- build_run_qc_figs(qc_fixture())
  p <- figs[["01_valis_error_by_stage"]]
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_equal(unname(geoms), "GeomBoxplot")   # no points, no printed medians
  # No facet: the pre/post split was the thing that made two near-identical panels.
  expect_s3_class(p$facet, "FacetNull")
  # Log y: on a linear axis the registered stages sat on zero under one bad slide.
  # ggplot2 < 3.5 keeps the transform in $trans and $transform is a method; 3.5+
  # renamed the field to $transform. Read whichever one is the object.
  y  <- p$scales$get_scales("y")
  tr <- if (is.function(y$transform)) y$trans else y$transform
  expect_equal(tr$name, "log-10")
})

test_that("valis_error_long carries a caller's extra grouping columns through", {
  # The arm sweep attaches `arm` and friends before handing the frame over. A fixed
  # select() dropped them, and the arm figure then silently stopped being built.
  v <- read_valis_summary(qc_fixture()) |> dplyr::mutate(arm = "high / micro 2")
  long <- valis_error_long(v)
  expect_true("arm" %in% names(long))
  expect_equal(unique(long$arm), "high / micro 2")
})

test_that("VALIS summaries are found under their producer subdirectory, and flat", {
  # The regression: Nextflow preserves REGISTER's `preprocessed/data/` on publish, so a
  # one-level dir_ls() found nothing on a real run — and, because a patient is detected
  # by finding FILES, `registered/summary/` holding only a directory meant a VALIS run
  # with reg_qc < 2 was not detected at all. Both depths must work: a hand-copied or
  # flattened tree is still a legitimate input.
  for (sub in list(c("preprocessed", "data"), character(0))) {
    root <- file.path(tempdir(), paste0("depth-", sample(1e6, 1)))
    d <- do.call(file.path, as.list(c(root, "046", "registered", "summary", sub)))
    dir.create(d, recursive = TRUE)
    readr::write_csv(tibble::tibble(img_name = "cycle2", original_rTRE = 50,
                                    rigid_rTRE = 10, non_rigid_rTRE = 3, n_matches = 99L),
                     file.path(d, "046_summary.csv"))
    expect_equal(basename(.qc_patient_dirs(root)), "046")
    expect_equal(nrow(read_valis_summary(root)), 1)
    # original + rigid + non_rigid; no premicro sibling, so no micro stage.
    expect_equal(nrow(valis_error_long(read_valis_summary(root))), 3)
  }
})

test_that("a symlinked patient directory is still a patient directory", {
  # Assembling data/mirage/ by symlinking individual patient directories out of one or
  # more run outdirs is a normal way to work off-repo. dir_ls(type = "directory") sees a
  # symlink's OWN type ("symlink") and skipped every one of them.
  src <- qc_fixture(list(`046` = "cycle2"))
  root <- file.path(tempdir(), paste0("symlinked-", sample(1e6, 1)))
  dir.create(root, recursive = TRUE)
  skip_if_not(file.symlink(file.path(src, "046"), file.path(root, "046")),
              "filesystem does not support symlinks")
  expect_equal(basename(.qc_patient_dirs(root)), "046")
  expect_equal(nrow(read_seg_qc(root)), 4)                 # four stages, one slide
  expect_true("01_valis_error_by_stage" %in% names(build_run_qc_figs(root)))
})

test_that("the intrinsic-vs-overlap view joins across the patient-prefix spelling", {
  # TILED_SOLVE names its artifact <patient_id>_<channels>_tre.json while seg_qc carries
  # the native image stem, which need not repeat the patient id. Joining on the raw
  # `moving` string dropped every row, and an empty join removes §4 silently.
  root <- file.path(tempdir(), paste0("token-", sample(1e6, 1)))
  d <- file.path(root, "046", "qc", "registration")
  dir.create(d, recursive = TRUE)
  stages <- c("rigid", "micro")
  tiles  <- data.frame(ix = 0:1, iy = 0L, cx = c(0, 512), cy = 0, tre_rigid = c(1, 2))
  for (mv in c("cycle2", "cycle3")) {
    jsonlite::write_json(list(
      patient_id = "046", moving = paste0(mv, ".ome.tiff"), stage_order = stages,
      matching = list(pair_fraction = 0.7),
      stages = stats::setNames(lapply(c(.5, .8), function(x) list(
        n_pairs = 100L, dice_matched = x, displacement_um_p50 = 2,
        displacement_um_p90 = 4)), stages)),
      file.path(d, paste0(mv, "_seg_qc.json")), auto_unbox = TRUE)
    jsonlite::write_json(list(
      coarse_tre_px = 4.2, n_inliers = 800L, n_tiles = 2L, mesh_refined = TRUE,
      moving = mv, rigid_tre_px = list(mean = 2, p50 = 2, p90 = 3, max = 3),
      residual_after_px = list(mean = .6, p50 = .6, p90 = .9, max = .9), tiles = tiles),
      file.path(d, paste0("046_", mv, "_tre.json")), auto_unbox = TRUE)  # prefixed
  }

  sq <- read_seg_qc(root); st <- read_stare_tre(root)
  expect_equal(nrow(dplyr::inner_join(sq, st, by = c("patient_id", "moving"))), 0)
  expect_setequal(sq$slide_token, st$slide_token)
  expect_true("04_intrinsic_vs_overlap" %in% names(build_run_qc_figs(root)))
})

test_that("STARE TRE and its per-tile spatial map are read when the tiled path ran", {
  root <- qc_fixture(tiled = TRUE)
  st <- read_stare_tre(root)
  expect_equal(nrow(st), 2)                       # one per moving slide
  expect_equal(unique(st$coarse_tre_px), 4.2)
  expect_equal(unique(st$after_p50), .6)
  tiles <- read_stare_tiles(root)
  expect_equal(nrow(tiles), 2 * 6)
  expect_true(all(c("ix", "iy", "tre_rigid", "patient_id", "moving") %in% names(tiles)))
  # a VALIS run has no STARE artifacts at all
  expect_equal(nrow(read_stare_tre(qc_fixture(tiled = FALSE))), 0)
})

test_that("each registration path builds its own intrinsic section", {
  valis_figs <- build_run_qc_figs(qc_fixture(tiled = FALSE))
  stare_figs <- build_run_qc_figs(qc_fixture(tiled = TRUE))
  expect_true("01_valis_error_by_stage" %in% names(valis_figs))
  expect_false(any(grepl("stare", names(valis_figs))))
  expect_true(all(c("02_stare_tre", "02b_stare_tre_map") %in% names(stare_figs)))
  expect_false(any(grepl("valis", names(stare_figs))))
  # the independent overlap check and the agreement view exist either way
  for (f in list(valis_figs, stare_figs))
    expect_true(all(c("03_overlap_dice_by_stage", "03c_pair_fraction") %in% names(f)))
})

test_that("every figure actually renders", {
  # a ggplot object can be malformed and only fail at draw time (aes() referring to
  # something that is not a column, say)
  figs <- c(build_run_qc_figs(qc_fixture(tiled = FALSE)),
            build_run_qc_figs(qc_fixture(tiled = TRUE)))
  expect_gt(length(figs), 8)
  for (nm in names(figs))
    expect_no_error(ggplot2::ggsave(file.path(tempdir(), paste0(nm, ".png")),
                                    figs[[nm]], width = 6, height = 4, dpi = 50))
})

test_that("a patient is detected from ANY artifact subset it produced", {
  # Runs differ in what they emit: reg_qc < 2 produces no *_seg_qc.json, the VALIS
  # path produces no *_tre.json. Keying
  # detection on one location would render an empty page for the others.
  root <- file.path(tempdir(), paste0("subset-", sample(1e6, 1)))

  # (a) VALIS summaries only — no qc/ directory at all. Published at their real depth,
  # so `registered/summary/` holds a DIRECTORY and no files: detection has to recurse
  # or this patient disappears entirely and the page claims the run produced no QC.
  vdir <- file.path(root, "046", "registered", "summary", "preprocessed", "data")
  dir.create(vdir, recursive = TRUE)
  readr::write_csv(tibble::tibble(img_name = "cycle2", original_rTRE = 50,
                                  rigid_rTRE = 10, non_rigid_rTRE = 3, n_matches = 99L),
                   file.path(vdir, "046_summary.csv"))
  # (b) staged seg_qc only — no VALIS summary copied in. The mirror image of (a), and
  # a real shape: the summaries live under a producer subdirectory that is easy to miss
  # when hand-copying a run.
  qdir <- file.path(root, "052", "qc", "registration")
  dir.create(qdir, recursive = TRUE)
  jsonlite::write_json(list(
    patient_id = "052", moving = "052_cycle2", stages_separable = TRUE,
    stage_order = c("native", "rigid", "non_rigid"),
    matching = list(anchor_stage = "rigid", n_pairs = 100, pair_fraction = 0.8),
    stages = list(
      native    = list(n_pairs = 100, dice_matched = .4, displacement_um_p50 = 20),
      rigid     = list(n_pairs = 100, dice_matched = .8, displacement_um_p50 = 4),
      non_rigid = list(n_pairs = 100, dice_matched = .9, displacement_um_p50 = 1.2))),
    file.path(qdir, "052_cycle2_seg_qc.json"), auto_unbox = TRUE)
  # (c) cohort-level directories, one of them non-empty — must never be a patient
  dir.create(file.path(root, "size_logs"), recursive = TRUE)
  dir.create(file.path(root, "qc"), recursive = TRUE)
  writeLines("<html>", file.path(root, "qc", "report.html"))

  expect_setequal(basename(.qc_patient_dirs(root)), c("046", "052"))
  expect_equal(nrow(read_valis_summary(root)), 1)
  expect_gt(nrow(read_seg_qc(root)), 0)
  # and the page still builds the section each subset supports
  figs <- build_run_qc_figs(root)
  expect_true("01_valis_error_by_stage" %in% names(figs))
  expect_true("03_overlap_dice_by_stage" %in% names(figs))
})

test_that("a phenotyping-only directory is no longer a patient", {
  # Conformal phenotyping QC was removed, so `phenotyping/` is not an artifact this
  # page reads. A directory holding only it must not be detected as a patient — doing
  # so would put an empty row in every table for a run with no registration QC at all.
  root <- file.path(tempdir(), paste0("phenoonly-", sample(1e6, 1)))
  dir.create(file.path(root, "052", "phenotyping"), recursive = TRUE)
  writeLines("label,phenotype", file.path(root, "052", "phenotyping", "phenotypes.csv"))
  expect_equal(length(.qc_patient_dirs(root)), 0)
  expect_equal(nrow(run_qc_tables(root)$seg_qc), 0)
})

# ---- VALIS's real two-column schema -----------------------------------------
# The fixtures above write an `original_rTRE` column to exercise the tolerance path,
# but VALIS's actual error_df is `from`/`filename`, `rigid_D`, `non_rigid_D` and nothing
# else (bin/generate_qc_report.py iterates exactly those two metric columns). These
# tests pin the reconstruction against the real schema.

.valis_real <- function(root, pid = "046", premicro = TRUE, nr_final = 3, nr_pre = 5) {
  d <- file.path(root, pid, "registered", "summary", "preprocessed", "data")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  for (sfx in c("_summary.csv", if (premicro) "_summary_premicro.csv"))
    readr::write_csv(tibble::tibble(
      filename       = c("cycle2", "cycle3"),
      rigid_D        = c(10, 11),
      non_rigid_D    = if (sfx == "_summary.csv") nr_final else nr_pre),
      file.path(d, paste0(pid, sfx)))
  root
}

test_that("two columns and two files reconstruct three stages", {
  # Micro-registration has no column: it updates the non-rigid field, so its value
  # lives in the DIFFERENCE BETWEEN THE FILES. rigid_D is unchanged by it.
  long <- valis_error_long(read_valis_summary(
    .valis_real(file.path(tempdir(), paste0("vreal-", sample(1e6, 1))))))

  expect_setequal(unique(as.character(long$stage)), c("rigid", "non_rigid", "micro"))
  expect_false("original" %in% as.character(long$stage))
  # rigid comes from the final file and appears ONCE, not once per file.
  expect_equal(sum(long$stage == "rigid"), 2)                    # two moving slides
  expect_setequal(long$error[long$stage == "rigid"], c(10, 11))
  # Same column, two files, two stages.
  expect_equal(unique(long$error[long$stage == "non_rigid"]), 5)
  expect_equal(unique(long$source_file[long$stage == "non_rigid"]), "pre-micro")
  expect_equal(unique(long$error[long$stage == "micro"]), 3)
  expect_equal(unique(long$source_file[long$stage == "micro"]), "final")
  # `_D` columns, so the metric label must say distance, not relative rTRE.
  expect_match(long$metric[1], "distance")
})

test_that("a blank micro stage means micro did not run, not that it gained nothing", {
  # No premicro file => reg_micro_reg < 2. The final non_rigid_D IS the pre-micro value,
  # so it becomes `non_rigid` and NO micro stage is emitted. A duplicated non_rigid box
  # would read as "micro bought nothing", a different claim. mirage's seg-QC side omits
  # the stage for the same reason.
  long <- valis_error_long(read_valis_summary(
    .valis_real(file.path(tempdir(), paste0("vnomicro-", sample(1e6, 1))),
                premicro = FALSE)))
  expect_setequal(unique(as.character(long$stage)), c("rigid", "non_rigid"))
  expect_false("micro" %in% as.character(long$stage))
  expect_equal(unique(long$error[long$stage == "non_rigid"]), 3)   # final, not duplicated
  expect_equal(unique(long$source_file), "final")
})
