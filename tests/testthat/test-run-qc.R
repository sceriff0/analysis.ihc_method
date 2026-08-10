# mirage's QC of a run on the study samples. These build a synthetic outdir rather
# than needing a pipeline run — the schemas are pinned in code/run_qc.R's header.
source(here::here("code", "run_qc.R"))

qc_fixture <- function(patients = list(`046` = c("cycle2", "cycle3")),
                       tiled = FALSE, phenotyping = TRUE) {
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
      dir.create(file.path(d, "registered", "summary"), recursive = TRUE, showWarnings = FALSE)
      mv <- patients[[pid]]
      # post-micro AND pre-micro: register_micro() overwrites the plain summary, so
      # the two files are different measurements of the same slide
      for (suffix in c("_summary.csv", "_summary_premicro.csv"))
        readr::write_csv(tibble::tibble(
          img_name = mv, original_rTRE = 50, rigid_rTRE = 10,
          non_rigid_rTRE = if (suffix == "_summary.csv") 3 else 5,
          n_matches = 120L),
          file.path(d, "registered", "summary", paste0(pid, suffix)))
    }
    if (phenotyping) {
      dir.create(file.path(d, "phenotyping"), recursive = TRUE, showWarnings = FALSE)
      jsonlite::write_json(list(
        chosen_alpha = .038, alpha_target = .05, crc_ran = TRUE, reporting_mode = FALSE,
        degraded_markers = c("FOXP3"), n_cells = 100000L, density_radius = 30, n_bins = 4L),
        file.path(d, "phenotyping", "phenotype_qc.json"), auto_unbox = TRUE)
      readr::write_csv(tibble::tibble(
        id = 1:2, markers = c("PANCK|CD45", "CD4|CD8"), observed = c(.002, .031),
        nominal = c(0, .01), density_corr = 0, neighbour_contact_corr = 0,
        verdict = c("pass", "warn")),
        file.path(d, "phenotyping", "constraint_audit.csv"))
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

test_that("VALIS pre- and post-micro summaries are tagged, never pooled", {
  # register_micro() overwrites <name>_summary.csv, so the plain file is post-micro
  # and only the premicro one isolates the non-rigid stage
  v <- read_valis_summary(qc_fixture())
  expect_setequal(v$stage_scope, c("post-micro", "pre-micro"))
  long <- valis_error_long(v)
  # the fixture has two moving slides, so one row each per scope
  expect_equal(unique(dplyr::filter(long, stage == "non_rigid",
                                    stage_scope == "post-micro")$error), 3)
  expect_equal(unique(dplyr::filter(long, stage == "non_rigid",
                                    stage_scope == "pre-micro")$error), 5)
  expect_true(grepl("rTRE", long$metric[1]))
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

test_that("phenotyping calibration and the constraint audit are read per patient", {
  root <- qc_fixture(list(`046` = "cycle2", `EPM - 052` = "cycle2"))
  pq <- read_phenotype_qc(root)
  expect_setequal(pq$patient_id, c("046", "052"))  # slide_key normalises the dir name
  expect_true(all(pq$crc_ran))
  expect_equal(unique(pq$n_degraded), 1)
  au <- read_constraint_audit(root)
  expect_equal(nrow(au), 4)
  expect_equal(sum(au$verdict == "warn"), 2)
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
  # path produces no *_tre.json, a registration-only run has no phenotyping. Keying
  # detection on one location would render an empty page for the others.
  root <- file.path(tempdir(), paste0("subset-", sample(1e6, 1)))

  # (a) VALIS summaries only — no qc/ directory at all
  dir.create(file.path(root, "046", "registered", "summary"), recursive = TRUE)
  readr::write_csv(tibble::tibble(img_name = "cycle2", original_rTRE = 50,
                                  rigid_rTRE = 10, non_rigid_rTRE = 3, n_matches = 99L),
                   file.path(root, "046", "registered", "summary", "046_summary.csv"))
  # (b) phenotyping only — no registration QC at all
  dir.create(file.path(root, "052", "phenotyping"), recursive = TRUE)
  jsonlite::write_json(list(chosen_alpha = .04, alpha_target = .05, crc_ran = TRUE,
                            reporting_mode = FALSE, degraded_markers = character(0),
                            n_cells = 10L, density_radius = 30, n_bins = 4L),
                       file.path(root, "052", "phenotyping", "phenotype_qc.json"),
                       auto_unbox = TRUE)
  # (c) cohort-level directories, one of them non-empty — must never be a patient
  dir.create(file.path(root, "size_logs"), recursive = TRUE)
  dir.create(file.path(root, "qc"), recursive = TRUE)
  writeLines("<html>", file.path(root, "qc", "report.html"))

  expect_setequal(basename(.qc_patient_dirs(root)), c("046", "052"))
  expect_equal(nrow(read_valis_summary(root)), 1)
  expect_equal(nrow(read_phenotype_qc(root)), 1)
  # and the page still builds the section each subset supports
  figs <- build_run_qc_figs(root)
  expect_true("01_valis_error_by_stage" %in% names(figs))
  expect_true("05_phenotype_alpha" %in% names(figs))
})
