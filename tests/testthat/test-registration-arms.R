# The arm sweep: the study slides registered once per configuration. These tests are
# mostly about the COMPARABILITY RULES, because that is where this analysis can be
# wrong while still producing a clean, plausible figure — mirage's `rigid` QC stage
# changes meaning with micro-registration depth, and depths 0 and 1 emit no `micro`
# stage at all.
source(here::here("code", "registration_arms.R"))

# A synthetic sweep shaped like a real mirage outdir set. `micro` is emitted only at
# depth 2, exactly as the pipeline does — that asymmetry is the point of the fixture.
arms_tree <- function(modes = c("high", "low"), depths = 0:2,
                      patients = c("046", "052"), manifest = FALSE,
                      tiled = FALSE) {
  root <- file.path(tempdir(), paste0("arms-", as.integer(Sys.time()), "-", sample(1e6, 1)))
  if (tiled) {
    # The tiled backend's OWN stage vocabulary: native -> rigid -> refined. Writing
    # `refined` is the point of this fixture — a VALIS-only stage list silently
    # dropped it, and the arm's final stage came out as `rigid`.
    for (p in patients) {
      d <- file.path(root, "tiled_defaults", p, "qc", "registration")
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      vals <- c(native = 21, rigid = 3.4, refined = 1.15)
      st <- lapply(names(vals), function(s) list(
        n_pairs = 1000, iou_mean = 0.9 - vals[[s]] / 40,
        dice_matched = 0.95 - vals[[s]] / 60,
        displacement_px_p50 = vals[[s]] / 0.325,
        displacement_um_p50 = vals[[s]]))
      names(st) <- names(vals)
      jsonlite::write_json(list(
        patient_id = p, moving = paste0(p, "_cycle2"), stages_separable = TRUE,
        stage_order = c("native", "rigid", "refined"),
        matching = list(anchor_stage = "rigid", n_pairs = 1000, pair_fraction = 0.84),
        stages = st,
        delta_vs_anchor = list(
          native  = list(displacement_um_p50 = vals[["native"]]  - vals[["rigid"]]),
          refined = list(displacement_um_p50 = vals[["refined"]] - vals[["rigid"]]))),
        file.path(d, paste0(p, "_cycle2_seg_qc.json")), auto_unbox = TRUE)
      # STARE's intrinsic TRE lives in a different file, in PIXELS, with tiles.
      jsonlite::write_json(list(
        moving = "DAPI", coarse_tre_px = 8.4, n_inliers = 1420, n_tiles = 4,
        mesh_refined = TRUE,
        rigid_tre_px = list(p50 = 3.4, p90 = 7.1, max = 12),
        residual_after_px = list(p50 = 1.15, p90 = 2.6),
        tiles = data.frame(x = c(0, 1, 0, 1), y = c(0, 0, 1, 1),
                           tre_px = c(0.9, 1.4, 2.2, 1.1))),
        file.path(d, paste0(p, "_DAPI_tre.json")), auto_unbox = TRUE)
    }
  }
  for (mode in modes) {
    base <- if (mode == "high") 1.0 else 1.6
    for (dp in depths) {
      arm <- sprintf("valis_%s_micro%d", mode, dp)
      stages <- c("native", "rigid", "non_rigid", if (dp == 2) "micro")
      for (p in patients) {
        d <- file.path(root, arm, p, "qc", "registration")
        dir.create(d, recursive = TRUE, showWarnings = FALSE)
        # micro-rigid folds into `rigid`, so rigid is tighter at depth >= 1 — the
        # definition change these tests exist to protect against.
        rigid <- base * if (dp == 0) 4.0 else 3.0
        vals  <- c(native = base * 20, rigid = rigid,
                   non_rigid = rigid * .35, micro = rigid * .28)
        st <- lapply(stages, function(s) list(
          n_pairs = 1000, iou_mean = 0.9 - vals[[s]] / 40,
          dice_matched = 0.95 - vals[[s]] / 60,
          displacement_px_p50 = vals[[s]] / 0.325,
          displacement_um_p50 = vals[[s]],
          displacement_um_p90 = vals[[s]] * 2.2))
        names(st) <- stages
        dv <- lapply(setdiff(stages, "rigid"), function(s) list(
          dice_matched = st[[s]]$dice_matched - st$rigid$dice_matched,
          displacement_um_p50 = st[[s]]$displacement_um_p50 - st$rigid$displacement_um_p50))
        names(dv) <- setdiff(stages, "rigid")
        jsonlite::write_json(list(
          patient_id = p, moving = paste0(p, "_cycle2"), reference = paste0(p, "_cycle1"),
          stages_separable = TRUE, stage_order = stages,
          matching = list(anchor_stage = "rigid", n_pairs = 1000,
                          pair_fraction = if (mode == "high") 0.9 else 0.62),
          stages = st, delta_vs_anchor = dv),
          file.path(d, paste0(p, "_cycle2_seg_qc.json")), auto_unbox = TRUE)
      }
    }
  }
  if (manifest) {
    readr::write_csv(tibble::tibble(
      arm_dir     = as.vector(outer(modes, depths, function(m, d) sprintf("valis_%s_micro%d", m, d))),
      memory_mode = rep(modes, times = length(depths)),
      micro_reg   = rep(depths, each = length(modes)),
      label       = as.vector(outer(modes, depths, function(m, d) sprintf("PRESET-%s D%d", m, d)))),
      file.path(root, "arms.csv"))
  }
  root
}

test_that("directory names parse into the two knobs", {
  expect_equal(.parse_arm_dir("valis_high_micro2")$memory_mode, "high")
  expect_equal(.parse_arm_dir("valis_high_micro2")$micro_reg, 2L)
  expect_equal(.parse_arm_dir("low-mr0")$memory_mode, "low")
  expect_equal(.parse_arm_dir("low-mr0")$micro_reg, 0L)
  # Unparseable names must not vanish; they keep NA knobs and are labelled by dir.
  expect_true(is.na(.parse_arm_dir("run_alpha")$memory_mode))
  expect_true(is.na(.parse_arm_dir("run_alpha")$micro_reg))
})

test_that("a preset x depth sweep is six arms, not four", {
  man <- arm_manifest(arms_tree())
  expect_equal(nrow(man), 6)
  expect_setequal(man$micro_reg, c(0L, 0L, 1L, 1L, 2L, 2L))
  expect_setequal(unique(man$memory_mode), c("high", "low"))
})

test_that("arms.csv overrides the name parse", {
  man <- arm_manifest(arms_tree(manifest = TRUE))
  expect_true(all(grepl("^PRESET-", man$arm)))
  expect_equal(nrow(man), 6)
})

test_that("a directory with no mirage QC is not counted as an arm", {
  root <- arms_tree()
  dir.create(file.path(root, "figures"), showWarnings = FALSE)   # a sibling, not an arm
  expect_equal(nrow(arm_manifest(root)), 6)
})

test_that("only depth 2 reports a micro stage", {
  seg <- read_arms_seg_qc(arm_manifest(arms_tree()))
  by_depth <- seg |> dplyr::filter(stage == "micro") |> dplyr::distinct(micro_reg)
  expect_equal(sort(by_depth$micro_reg), 2L)
  # Every arm still reports the other three, so nothing was dropped wholesale.
  expect_setequal(unique(as.character(seg$stage[seg$micro_reg == 0])),
                  c("native", "rigid", "non_rigid"))
})

test_that("the final stage is per-arm, not a fixed stage name", {
  fin <- arm_final_stage(read_arms_seg_qc(arm_manifest(arms_tree())))
  # Ranking on "the micro stage" would have kept only the two depth-2 arms.
  expect_equal(nrow(dplyr::distinct(fin, arm)), 6)
  expect_equal(unique(as.character(fin$final_stage[fin$micro_reg == 2])), "micro")
  expect_setequal(unique(as.character(fin$final_stage[fin$micro_reg < 2])), "non_rigid")
})

test_that("a depth-crossed sweep refuses to compare `rigid` across arms", {
  # THE central guard. At depth >= 1 mirage's `rigid` stage already contains the
  # micro-rigid refinement, so `rigid` on a shared axis is two different transforms.
  seg <- read_arms_seg_qc(arm_manifest(arms_tree()))
  expect_equal(arm_comparable_stages(seg), "native")
  expect_false("rigid" %in% arm_comparable_stages(seg))
  expect_match(comparable_stage_note(seg), "micro-rigid")
})

test_that("a single-depth, single-backend sweep allows every stage it reported", {
  seg <- read_arms_seg_qc(arm_manifest(arms_tree(depths = 2)))
  # Every stage PRESENT, not every stage the vocabulary knows: a VALIS-only sweep must
  # not advertise the tiled backend's `refined` as comparable when nothing emitted it.
  expect_setequal(arm_comparable_stages(seg), c("native", "rigid", "non_rigid", "micro"))
  expect_false("refined" %in% arm_comparable_stages(seg))
  expect_match(comparable_stage_note(seg), "same thing across arms")
})

test_that("the ranking table has one row per arm and is ordered by residual", {
  tbl <- arm_ranking_table(read_arms_seg_qc(arm_manifest(arms_tree())))
  expect_equal(nrow(tbl), 6)
  expect_equal(tbl$disp_um_p50, sort(tbl$disp_um_p50))
  # pair_fraction rides along because it gates every other number in the table.
  expect_true("pair_fraction" %in% names(tbl))
})

test_that("the figures build and the stage ladder is faceted, never pooled", {
  man  <- arm_manifest(arms_tree())
  figs <- build_arm_figs(read_arms_seg_qc(man), tibble::tibble(), man)
  expect_true(all(c("01_final_residual_um_by_arm", "02_final_dice_by_arm",
                    "03_stage_ladder_within_arm", "04_pair_fraction_by_arm") %in% names(figs)))
  # Faceting by arm IS the comparability guard made visual: the stages must never
  # share one axis across arms.
  expect_false(inherits(figs[["03_stage_ladder_within_arm"]]$facet, "FacetNull"))
})

test_that("an empty arms directory warns rather than erroring", {
  expect_warning(man <- arm_manifest(file.path(tempdir(), "no-arms-here")), "no directory")
  expect_equal(nrow(man), 0)
  expect_equal(length(build_arm_figs(tibble::tibble(), tibble::tibble())), 0)
})


# ---- the tiled (STARE) backend ----------------------------------------------

test_that("a tiled directory is read as the other backend, with no invented knobs", {
  # reg_micro_reg and memory_mode are VALIS-only params. Parsing them off a tiled
  # run's name would attribute knob values to a run that never had them.
  d <- .parse_arm_dir("tiled_defaults")
  expect_equal(d$backend, "tiled")
  expect_true(is.na(d$memory_mode))
  expect_true(is.na(d$micro_reg))
  expect_equal(.parse_arm_dir("stare_run")$backend, "tiled")
  expect_equal(.parse_arm_dir("valis_high_micro2")$backend, "valis")
})

test_that("a tiled arm is not scanned for a preset even when its name contains one", {
  # "tiled_high_throughput" must not become memory_mode = high.
  d <- .parse_arm_dir("tiled_high_throughput")
  expect_equal(d$backend, "tiled")
  expect_true(is.na(d$memory_mode))
})

test_that("the tiled `refined` stage survives the read", {
  # THE REGRESSION THIS GUARDS. read_seg_qc() factors `stage` against a known list and
  # drops non-matches; with a VALIS-only list every `refined` row vanished, no error
  # and no warning, and the arm's final stage came out as `rigid`.
  seg <- read_arms_seg_qc(arm_manifest(arms_tree(modes = "high", depths = 2, tiled = TRUE)))
  tl  <- dplyr::filter(seg, backend == "tiled")
  expect_gt(nrow(tl), 0)
  expect_true("refined" %in% as.character(tl$stage))
  expect_setequal(unique(as.character(tl$stage)), c("native", "rigid", "refined"))
})

test_that("each backend's final stage is its own last stage", {
  fin <- arm_final_stage(read_arms_seg_qc(
    arm_manifest(arms_tree(modes = "high", depths = 2, tiled = TRUE))))
  expect_equal(unique(as.character(fin$final_stage[fin$backend == "tiled"])), "refined")
  expect_equal(unique(as.character(fin$final_stage[fin$backend == "valis"])), "micro")
})

test_that("stage_index comes from the run's own stage_order, not a global ordering", {
  # The two vocabularies interleave under any single factor ordering, so ranking by
  # factor level would compare a STARE stage's position against a VALIS one.
  seg <- read_arms_seg_qc(arm_manifest(arms_tree(modes = "high", depths = 2, tiled = TRUE)))
  expect_true("stage_index" %in% names(seg))
  tl <- dplyr::filter(seg, backend == "tiled") |> dplyr::distinct(stage, stage_index)
  expect_equal(tl$stage_index[tl$stage == "native"], 1L)
  expect_equal(tl$stage_index[tl$stage == "refined"], 3L)
})

test_that("crossing backends restricts comparison to native, whatever the depth", {
  # Stronger than the depth restriction: matching micro depth cannot make `rigid`
  # comparable between backends, because it is a different operation in each.
  seg <- read_arms_seg_qc(arm_manifest(arms_tree(modes = "high", depths = 2, tiled = TRUE)))
  expect_equal(arm_comparable_stages(seg), "native")
  expect_match(comparable_stage_note(seg), "BACKENDS")
  expect_match(comparable_stage_note(seg), "shared word rather than a shared operation")
})

test_that("the tiled arm joins the ranking with NA knobs, not a dropped row", {
  tbl <- arm_ranking_table(read_arms_seg_qc(arm_manifest(arms_tree(tiled = TRUE))))
  expect_equal(nrow(tbl), 7)                       # six VALIS cells + one comparator
  tl <- dplyr::filter(tbl, backend == "tiled")
  expect_equal(nrow(tl), 1)
  expect_true(is.na(tl$micro_reg))
  expect_equal(tl$final_stage, "refined")
})

test_that("the knob-effects figure excludes the tiled arm", {
  # A tiled run has no position on a preset x depth grid; plotting it there would
  # place a point on axes it was never run against.
  man  <- arm_manifest(arms_tree(tiled = TRUE))
  figs <- build_arm_figs(read_arms_seg_qc(man), tibble::tibble(), man)
  expect_false("tiled" %in% figs[["05_knob_effects"]]$data$backend)
})

test_that("STARE's own error and tile map are read into their own figures", {
  root <- arms_tree(modes = "high", depths = 2, tiled = TRUE)
  man  <- arm_manifest(root)
  st   <- read_arms_stare_tre(man)
  expect_gt(nrow(st), 0)
  expect_true(all(st$backend == "tiled"))          # VALIS arms contribute nothing
  expect_true(all(is.finite(st$after_p50)))
  figs <- build_arm_figs(read_arms_seg_qc(man), tibble::tibble(), man)
  expect_true("08_stare_intrinsic_tre_px" %in% names(figs))
  expect_true("09_stare_tile_error_map" %in% names(figs))
  # In PIXELS, and labelled as such, so it is never read against VALIS's rTRE.
  expect_match(figs[["08_stare_intrinsic_tre_px"]]$labels$y, "px")
})

test_that("a VALIS-only sweep produces no STARE figures", {
  man  <- arm_manifest(arms_tree(tiled = FALSE))
  expect_equal(nrow(read_arms_stare_tre(man)), 0)
  figs <- build_arm_figs(read_arms_seg_qc(man), tibble::tibble(), man)
  expect_false(any(grepl("stare", names(figs))))
})
