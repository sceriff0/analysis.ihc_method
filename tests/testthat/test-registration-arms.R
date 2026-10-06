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
                      tiled = FALSE, movings = "cycle2", valis = FALSE) {
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
      for (p in patients) for (mv in movings) {
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
          patient_id = p, moving = paste0(p, "_", mv), reference = paste0(p, "_cycle1"),
          stages_separable = TRUE, stage_order = stages,
          matching = list(anchor_stage = "rigid", n_pairs = 1000,
                          pair_fraction = if (mode == "high") 0.9 else 0.62),
          stages = st, delta_vs_anchor = dv),
          file.path(d, paste0(p, "_", mv, "_seg_qc.json")), auto_unbox = TRUE)
        # VALIS's own error, where a real run publishes it (two levels below
        # registered/summary/, see test-run-qc.R) and in VALIS's own spelling:
        # `from` / `to`, no `img_name`. One csv per patient listing every moving
        # slide. `original` is the no-registration number and `rigid` the rigid-only
        # one; a pre-micro sibling exists only at depth 2, exactly as
        # bin/register.py writes it, and is what proves micro ran.
        if (valis && mv == movings[1]) {
          sdir <- file.path(root, arm, p, "registered", "summary", "preprocessed", "data")
          dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
          suffixes <- c("_summary.csv", if (dp == 2) "_summary_premicro.csv")
          for (suffix in suffixes)
            readr::write_csv(tibble::tibble(
              from = paste0(p, "_", movings), to = paste0(p, "_cycle1"),
              original_rTRE  = base * 0.05,
              rigid_rTRE     = rigid / 100,
              non_rigid_rTRE = if (suffix == "_summary.csv") vals[[if (dp == 2) "micro" else "non_rigid"]] / 100
                               else vals[["non_rigid"]] / 100,
              n_matches = 120L),
              file.path(sdir, paste0(p, suffix)))
        }
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

test_that("a three-preset sweep builds every figure, each preset with its own colour", {
  # The sweep gained `medium`; a two-colour positional scale stopped the knit there.
  expect_equal(.parse_arm_dir("valis_medium_micro1")$memory_mode, "medium")
  man  <- arm_manifest(arms_tree(modes = c("low", "medium", "high"), tiled = TRUE))
  expect_setequal(stats::na.omit(man$memory_mode), c("low", "medium", "high"))
  figs <- build_arm_figs(read_arms_seg_qc(man), read_arms_valis(man), man)
  for (p in figs) expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  expect_identical(levels(figs[["05_knob_effects"]]$data$memory_mode),
                   c("low", "medium", "high"))
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

# --- the channel-pair unit: main panels pooled, supplementary panels split -----------

test_that("a QC record's channel pair is patient-free, so it is shared across patients", {
  seg <- read_arms_seg_qc(arm_manifest(arms_tree(movings = c("cycle2", "cycle3"))))
  expect_setequal(unique(seg$pair), c("cycle1 vs cycle2", "cycle1 vs cycle3"))
  # The same pair in both patients is ONE label; a patient prefix left in would make it two.
  expect_equal(dplyr::n_distinct(seg$patient_id[seg$pair == "cycle1 vs cycle2"]), 2)
})

test_that("a record with no `reference` key still gets a pair label", {
  seg <- read_arms_seg_qc(arm_manifest(arms_tree(modes = "high", depths = 2, tiled = TRUE)))
  tl  <- dplyr::filter(seg, backend == "tiled")      # the tiled fixture writes no reference
  expect_true(all(is.na(tl$reference)))
  expect_true(all(tl$pair == "cycle2"))
})

test_that("the main panels hold one point per channel pair, pooled over patients", {
  man  <- arm_manifest(arms_tree(movings = c("cycle2", "cycle3"), valis = TRUE))
  figs <- build_arm_figs(read_arms_seg_qc(man), read_arms_valis(man), man)
  for (nm in c("01_final_residual_um_by_arm", "02_final_dice_by_arm",
               "02b_final_valis_error_by_arm")) {
    d <- figs[[nm]]$data
    expect_equal(nrow(d), 6 * 2 * 2, info = nm)      # arms x patients x pairs
    expect_equal(as.character(rlang::quo_get_expr(figs[[nm]]$mapping$x)), "arm", info = nm)
    expect_match(figs[[nm]]$labels$subtitle, "channel pair", info = nm)
  }
})

test_that("VALIS's `from`/`to` spelling is read as a slide and a pair", {
  man  <- arm_manifest(arms_tree(modes = "high", depths = 0, valis = TRUE,
                                 movings = c("cycle2", "cycle3")))
  long <- valis_error_long(read_arms_valis(man))
  # Not the csv's file name: without `from` every slide of a patient shared one label.
  expect_setequal(unique(long$slide_token), c("cycle2", "cycle3"))
  expect_setequal(unique(long$pair), c("cycle1 vs cycle2", "cycle1 vs cycle3"))
})

test_that("each main panel has a by-patient and a by-channel-pair supplement", {
  man  <- arm_manifest(arms_tree(movings = c("cycle2", "cycle3"), valis = TRUE))
  figs <- build_arm_figs(read_arms_seg_qc(man), read_arms_valis(man), man)
  supp <- grep("^S[0-9]+_", names(figs), value = TRUE)
  expect_setequal(supp, c("S1_residual_um_by_patient", "S2_residual_um_by_channel_pair",
                          "S3_dice_by_patient", "S4_dice_by_channel_pair",
                          "S5_valis_error_by_patient", "S6_valis_error_by_channel_pair"))
  for (nm in supp) {
    # Per arm, never pooled across arms: the split is WITHIN a configuration.
    expect_false(inherits(figs[[nm]]$facet, "FacetNull"), info = nm)
    # The same points as the main panel, regrouped — nothing re-summarised.
    main <- figs[[c(residual_um = "01_final_residual_um_by_arm", dice = "02_final_dice_by_arm",
                    valis_error = "02b_final_valis_error_by_arm")[[
                      sub("^S[0-9]+_(.*)_by_.*$", "\\1", nm)]]]]
    expect_equal(nrow(figs[[nm]]$data), nrow(main$data), info = nm)
    expect_match(figs[[nm]]$labels$subtitle, "n = 2 patients, n = 2 channel pairs", info = nm)
  }
})

test_that("with no VALIS summary the VALIS panels are skipped, not drawn empty", {
  man  <- arm_manifest(arms_tree())
  figs <- build_arm_figs(read_arms_seg_qc(man), tibble::tibble(), man)
  expect_false(any(grepl("valis_error", names(figs))))
  expect_true(all(c("S1_residual_um_by_patient", "S4_dice_by_channel_pair") %in% names(figs)))
})

test_that("an ashlar arm is keyed as ashlar, not as plain valis", {
  k <- .arm_kind(c("ashlar", "tiled", "valis", "valis"), c(NA, NA, NA, 2L))
  expect_equal(as.character(k), c("ashlar", "tiled (STARE)", "valis", "valis \u00b7 micro 2"))
  expect_false(anyNA(k))
})

test_that("pair labels that no two patients share are flagged, not drawn silently", {
  man <- arm_manifest(arms_tree(modes = "high", depths = 0))
  seg <- read_arms_seg_qc(man)
  expect_no_warning(build_arm_figs(seg, tibble::tibble(), man))
  # A stem that embeds the patient id mid-name survives .slide_token()'s prefix strip.
  seg$pair <- paste0("scan_", seg$patient_id, "_cycle2")
  expect_warning(build_arm_figs(seg, tibble::tibble(), man), "shared between patients")
})

# --- The manuscript panels: Fig 4(b)/(c) and Additional file 2 ----------------
test_that("the paper TRE panel ranks VALIS arms at their final stage and adds both baselines", {
  man  <- arm_manifest(arms_tree(tiled = TRUE, valis = TRUE))
  seg  <- read_arms_seg_qc(man)
  val  <- read_arms_valis(man)
  figs <- build_arm_paper_figs(seg, val, man)
  expect_setequal(names(figs), c("tre_by_arm", "dice_by_arm"))

  d <- figs$tre_by_arm$data
  expect_true(all(c(BASELINE_NONE, BASELINE_RIGID) %in% levels(d$arm)))
  # The baselines are the FIRST levels so coord_flip() puts them at the bottom.
  expect_equal(levels(d$arm)[1:2], c(BASELINE_NONE, BASELINE_RIGID))
  # Six VALIS arms + two baselines; the tiled arm has no VALIS error and is absent.
  expect_equal(nlevels(d$arm), 8)
  expect_false(any(grepl("tiled", levels(d$arm))))
  # "rigid only" is honest only where rigid means affine alone: depth 0.
  expect_setequal(unique(d$micro_reg[d$arm == BASELINE_RIGID]), 0L)
  expect_setequal(unique(d$micro_reg[d$arm == BASELINE_NONE]),  0L)
  # Depth-2 arms are ranked on the post-micro number, not the pre-micro one.
  d2 <- d[d$micro_reg == 2 & !d$arm %in% c(BASELINE_NONE, BASELINE_RIGID), ]
  expect_equal(unique(as.character(d2$final_stage)), "micro")
  # No fit line and no coefficient — the legend's claim, enforced on the object.
  expect_false(any(vapply(figs$tre_by_arm$layers,
                          function(l) inherits(l$stat, "StatSmooth"), logical(1))))
})

test_that("the paper Dice panel holds the same arms as the TRE panel, and the tiled comparator only on request", {
  man  <- arm_manifest(arms_tree(tiled = TRUE, valis = TRUE))
  seg  <- read_arms_seg_qc(man)
  figs <- build_arm_paper_figs(seg, read_arms_valis(man), man)
  # The legend: "plotted for the same arms as (b)". (b) cannot hold the tiled arm,
  # so by default (c) does not either — the two row sets are identical.
  expect_equal(levels(figs$dice_by_arm$data$arm), levels(figs$tre_by_arm$data$arm))
  expect_false(any(grepl("tiled", levels(figs$dice_by_arm$data$arm))))

  figs <- build_arm_paper_figs(seg, read_arms_valis(man), man,
                               backends = c("valis", "tiled"))
  d <- figs$dice_by_arm$data
  expect_equal(levels(d$arm)[1:2], c(BASELINE_NONE, BASELINE_RIGID))
  expect_equal(nlevels(d$arm), 9)                       # 6 VALIS + tiled + 2 baselines
  expect_true(any(grepl("tiled", levels(d$arm))))
  expect_setequal(unique(d$micro_reg[d$arm == BASELINE_RIGID]), 0L)
  # Baselines come from the native / rigid STAGES of the depth-0 runs.
  expect_equal(unique(as.character(d$final_stage[d$arm == BASELINE_NONE])),  "native")
  expect_equal(unique(as.character(d$final_stage[d$arm == BASELINE_RIGID])), "rigid")
  # The fixture makes native the worst overlap, so the baseline must sit below every arm.
  med <- tapply(d$dice_matched, as.character(d$arm), stats::median)
  expect_true(med[[BASELINE_NONE]] < min(med[names(med) != BASELINE_NONE]))
})

test_that("the paper panels build without VALIS summaries, dropping only the TRE panel", {
  man  <- arm_manifest(arms_tree(valis = FALSE))
  figs <- build_arm_paper_figs(read_arms_seg_qc(man), read_arms_valis(man), man)
  expect_null(figs$tre_by_arm)
  expect_s3_class(figs$dice_by_arm, "ggplot")
})

test_that("Additional file 2 is one row per arm plus the two baselines, with both metrics", {
  man <- arm_manifest(arms_tree(tiled = TRUE, valis = TRUE))
  seg <- read_arms_seg_qc(man); val <- read_arms_valis(man)
  tab <- arm_paper_table(seg, val, man)
  expect_true(all(c("arm", "backend", "memory_mode", "micro_reg", "n_slides", "stage",
                    "valis_tre", "disp_um_p50", "dice_matched", "pair_fraction")
                  %in% names(tab)))
  expect_equal(nrow(tab), 9)                            # 6 + tiled + 2 baselines
  expect_true(all(c(BASELINE_NONE, BASELINE_RIGID) %in% tab$arm))
  expect_true(is.na(tab$valis_tre[grepl("tiled", tab$arm)]))
  expect_equal(tab$stage[tab$arm == BASELINE_NONE],  "native")
  expect_equal(tab$stage[tab$arm == BASELINE_RIGID], "rigid")
  # n counts SLIDES, which for the baselines is slides x depth-0 presets.
  expect_equal(tab$n_slides[tab$arm == BASELINE_NONE], 4L)
  expect_equal(tab$n_slides[tab$arm == "high / micro 2"], 2L)
})

# --- placeholder mode (code/placeholders.R) ---------------------------------------
# arms.csv is the FULL plan. With the switch on, an arm it lists with no tree yet —
# and a patient an arm has not reached — are synthesised, flagged, and drawn as such.
# With the switch off the same tree reads exactly as it always did.
ph_tree <- function() {
  root <- arms_tree(manifest = TRUE, valis = TRUE)
  unlink(file.path(root, "valis_low_micro2"), recursive = TRUE)        # not run yet
  unlink(file.path(root, "valis_high_micro1", "052"), recursive = TRUE) # half run
  root
}

test_that("placeholder mode OFF: a planned-but-unrun arm stays absent", {
  withr::local_options(ihc.placeholder_missing = FALSE)
  root <- ph_tree()
  man  <- suppressMessages(arm_manifest(root))
  expect_equal(nrow(man), 5)
  expect_false("on_disk" %in% names(man))
  seg <- read_arms_seg_qc(man)
  expect_false("is_placeholder" %in% names(seg))
  expect_false("valis_low_micro2" %in% seg$arm_dir)
})

test_that("placeholder mode ON: the plan's missing slides are synthesised, flagged, and real rows are untouched", {
  root <- ph_tree()
  withr::local_options(ihc.placeholder_missing = FALSE)
  seg0 <- read_arms_seg_qc(suppressMessages(arm_manifest(root)))
  val0 <- read_arms_valis(suppressMessages(arm_manifest(root)))

  withr::local_options(ihc.placeholder_missing = TRUE)
  man <- suppressMessages(arm_manifest(root))
  expect_equal(nrow(man), 6)
  expect_false(man$on_disk[man$arm_dir == "valis_low_micro2"])
  seg <- suppressMessages(read_arms_seg_qc(man))
  expect_identical(dplyr::select(seg[!seg$is_placeholder, ], -is_placeholder, -placeholder_rule),
                   seg0)

  syn <- seg[seg$is_placeholder, ]
  lm2 <- dplyr::filter(syn, arm_dir == "valis_low_micro2")
  # 2 patients x 1 moving slide x the depth-2 VALIS ladder, micro included.
  expect_equal(nrow(lm2), 8)
  expect_setequal(as.character(lm2$stage), c("native", "rigid", "non_rigid", "micro"))
  expect_true(all(grepl("^no row", lm2$placeholder_rule)))
  # A half-run arm is topped up for the missing patient only, on ITS OWN stage list:
  # depth 1 reported no micro stage, so none is invented.
  hm1 <- dplyr::filter(syn, arm_dir == "valis_high_micro1")
  expect_equal(nrow(hm1), 3)
  expect_true(all(hm1$patient_id == "052"))
  expect_false("micro" %in% as.character(hm1$stage))
  # delta-vs-rigid is undefined at the rigid stage itself and stays NA.
  expect_true(all(is.na(syn$d_disp_um_vs_rigid[syn$stage == "rigid"])))
  # One pairing per slide -> one pair_fraction per synthetic slide.
  expect_equal(dplyr::n_distinct(lm2$pair_fraction[lm2$patient_id == "046"]), 1)
  expect_s3_class(seg$stage, "factor")

  # VALIS's own summary: the pre-micro file exists only at depth 2.
  val <- suppressMessages(read_arms_valis(man))
  expect_identical(dplyr::select(val[!val$is_placeholder, ], -is_placeholder, -placeholder_rule),
                   val0)
  vsyn <- val[val$is_placeholder, ]
  expect_setequal(vsyn$stage_scope[vsyn$arm_dir == "valis_low_micro2"], c("final", "pre-micro"))
  expect_setequal(vsyn$stage_scope[vsyn$arm_dir == "valis_high_micro1"], "final")
  vl <- valis_error_long(val)
  expect_true("micro" %in% vl$stage[vl$is_placeholder %in% TRUE & vl$arm_dir == "valis_low_micro2"])

  # Every figure that draws a synthetic point says so; the tables count them.
  figs <- suppressWarnings(build_arm_figs(seg, val, man))
  expect_match(figs[["01_final_residual_um_by_arm"]]$labels$subtitle, "^PLACEHOLDER — ")
  expect_match(figs[["07_valis_intrinsic_by_arm"]]$labels$subtitle, "^PLACEHOLDER — ")
  rt <- arm_ranking_table(seg)
  expect_equal(sum(rt$n_placeholder), nrow(arm_final_stage(seg)[arm_final_stage(seg)$is_placeholder, ]))
  pt <- arm_paper_table(seg, val, man)
  expect_true(all(c("n_placeholder", "n_placeholder_tre") %in% names(pt)))
  pf <- build_arm_paper_figs(seg, val, man)
  expect_match(pf$dice_by_arm$labels$subtitle, "^PLACEHOLDER — ")
})

# --- DRAPE arms ------------------------------------------------------------------
# The DRAPE branch names its tiled-backend arms tiled_<tier>_s<stride>[_<cross>], and
# its arms.csv still writes the pre-rename "tiled (STARE, ...)" label.
test_that("stride arms are labelled STARE with tier and stride, never 'tiled (STARE, defaults)'", {
  expect_identical(.drape_arm_label("tiled_high_s128"), "STARE (high, stride 128 px)")
  expect_identical(.drape_arm_label("tiled_low_s64_segstardist"),
                   "STARE (low, stride 64 px) [QC seg: stardist]")
  expect_identical(.drape_arm_label("tiled_medium_s256_segcellsam_pairmutual_nn"),
                   "STARE (medium, stride 256 px) [QC seg: cellsam] [QC pairing: mutual_nn]")
  # The manifest's own label is kept, with the method named STARE.
  expect_identical(.drape_arm_label("tiled_high_s128", "tiled (STARE, high, stride 128 px) [QC seg: stardist]"),
                   "STARE (high, stride 128 px) [QC seg: stardist]")
  # Not a DRAPE arm: untouched.
  expect_identical(.drape_arm_label(c("tiled_defaults", "tiled_high_gate05", "valis_high_micro2"),
                                    c(NA, "tiled (STARE, high, gate 0.5 px)", "high / micro 2")),
                   c(NA, "tiled (STARE, high, gate 0.5 px)", "high / micro 2"))

  root <- arms_tree(modes = "high", depths = 2, tiled = TRUE)
  for (a in c("tiled_high_s128", "tiled_low_s64_segstardist"))
    fs::dir_copy(file.path(root, "tiled_defaults"), file.path(root, a))
  man <- suppressMessages(arm_manifest(root))
  lab <- stats::setNames(man$arm, man$arm_dir)
  expect_identical(unname(lab["tiled_high_s128"]), "STARE (high, stride 128 px)")
  expect_identical(unname(lab["tiled_low_s64_segstardist"]),
                   "STARE (low, stride 64 px) [QC seg: stardist]")
  expect_true(all(man$backend[man$arm_dir %in% c("tiled_high_s128", "tiled_low_s64_segstardist")] == "tiled"))
  expect_false(any(grepl("DRAPE|tiled \\(", lab[startsWith(names(lab), "tiled_") & grepl("_s[0-9]", names(lab))])))
})

# --- stage figures: the shipped arm of each backend, with the after-transform record ---
# Add mirage's `full_transform` record to every seg_qc.json of a fixture tree, the way
# bin/warp_seg_qc.py writes it: the final stage's name, its own Dice and its own
# pair_fraction (a different pair set from the rigid-anchored ladder).
add_full_transform <- function(root, dice = 0.71, pair_fraction = 0.55) {
  for (f in fs::dir_ls(root, recurse = TRUE, glob = "*_seg_qc.json")) {
    d <- jsonlite::fromJSON(f, simplifyVector = FALSE)
    last <- d$stage_order[[length(d$stage_order)]]
    d$full_transform <- list(
      stage = last, n_pairs = 800, dice_matched = dice, iou_mean = 0.6,
      displacement_um_p50 = 0.4,
      matching = list(anchor_stage = last, n_pairs = 800, pair_fraction = pair_fraction))
    jsonlite::write_json(d, f, auto_unbox = TRUE)
  }
  root
}

stage_fixture <- function(full = TRUE) {
  root <- arms_tree(modes = "high", depths = c(0, 2), tiled = TRUE,
                    movings = c("CD3_P53", "CD4_CD8"))
  fs::dir_copy(file.path(root, "tiled_defaults"), file.path(root, "tiled_high_s128"))
  if (full) add_full_transform(root)
  man <- suppressMessages(arm_manifest(root))
  list(man = man, seg = read_arms_seg_qc(man), full = read_arms_seg_qc_full(man))
}

test_that("the after-transform record is read with its own Dice and pair fraction", {
  fx <- stage_fixture()
  v  <- dplyr::filter(fx$full, arm_dir == "valis_high_micro2")
  expect_equal(nrow(v), 4)                       # 2 patients x 2 moving slides
  expect_true(all(v$dice_matched == 0.71))
  expect_true(all(v$pair_fraction == 0.55))      # NOT the ladder's 0.9
  expect_true(all(v$final_stage == "micro"))
  expect_true(all(dplyr::filter(fx$full, arm_dir == "tiled_high_s128")$final_stage == "refined"))
  # A tree scored before the record existed yields no row, not an error.
  expect_equal(nrow(stage_fixture(full = FALSE)$full), 0)
})

test_that("stage figures draw the shipped arm of each backend on its own stage axis", {
  fx <- stage_fixture()
  d  <- arm_stage_frame(fx$seg, fx$full)
  expect_setequal(unique(d$arm_dir), c("valis_high_micro2", "tiled_high_s128"))
  st <- function(a) as.character(sort(unique(d$stage[d$arm_dir == a])))
  # Rigid and the re-paired whole transform ONLY: the rungs between are scored on the
  # rigid pairs and read, beside a re-paired box, as a progression they are not.
  expect_identical(st("valis_high_micro2"), c("rigid", "full_transform"))
  expect_identical(st("tiled_high_s128"), c("rigid", "full_transform"))
  # Both panels name their backend: a bare "high / micro 2" beside STARE does not.
  expect_setequal(levels(d$arm), c("VALIS (high / micro 2)", "STARE (high, stride 128 px)"))
  # The last box holds the re-paired score, not a copy of the ladder's last stage.
  expect_true(all(d$dice_matched[d$stage == "full_transform"] == 0.71))
  expect_false(any(d$dice_matched[d$stage == "rigid"] == 0.71))
  # The channel set is patient-free, so one colour collects a pair across patients.
  expect_setequal(unique(d$channels[d$backend == "valis"]), c("CD3_P53", "CD4_CD8"))
})

test_that("each stage figure exists plain, by patient and by channels", {
  fx   <- stage_fixture()
  figs <- build_arm_stage_figs(fx$seg, fx$full)
  expect_setequal(names(figs), c(
    "stage_dice", "stage_dice_by_patient", "stage_dice_by_channels",
    "stage_residual_um", "stage_residual_um_by_patient", "stage_residual_um_by_channels"))
  for (p in figs) expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  colour_of <- function(p) rlang::as_label(p$layers[[2]]$mapping$colour)
  expect_length(figs$stage_dice$layers, 1)                       # boxes only
  expect_match(colour_of(figs$stage_dice_by_patient), "patient_id")
  expect_match(colour_of(figs$stage_dice_by_channels), "channels")
  # More levels than the house palette still get one distinct colour each.
  cols <- .many_cols(sprintf("p%02d", 1:14))
  expect_length(unique(cols), 14)
})

test_that("old scores draw the rigid stage only, and an absent arm falls back out loud", {
  fx <- stage_fixture(full = FALSE)
  d  <- arm_stage_frame(fx$seg, fx$full)
  expect_false("full_transform" %in% as.character(d$stage))
  expect_match(build_arm_stage_figs(fx$seg, fx$full)$stage_dice$labels$subtitle,
               "rigid stage only")
  seg <- dplyr::filter(fx$seg, arm_dir != "tiled_high_s128")
  expect_message(d2 <- arm_stage_frame(seg, fx$full), "tiled_high_s128")
  expect_true("tiled_defaults" %in% d2$arm_dir)
})

# --- the headline is scored after the whole transform ---------------------------
test_that("the ranking takes the re-paired record, with its own pair fraction, and says so", {
  fx  <- stage_fixture()
  fin <- arm_final_score(fx$seg, fx$full)
  expect_equal(nrow(fin), nrow(arm_final_stage(fx$seg)))      # same rows, other values
  expect_true(all(fin$scored_on == SCORED_FULL))
  expect_true(all(fin$dice_matched == 0.71))
  expect_true(all(fin$pair_fraction == 0.55))                 # NOT the anchor's
  expect_true(all(is.na(fin$d_disp_um_vs_rigid)))             # no delta across pair sets
  tbl <- arm_ranking_table(fx$seg, fx$full)
  expect_true(all(tbl$dice_matched == 0.71))
  expect_true(all(tbl$scored_on == SCORED_FULL))
  # The pair-fraction and Dice panels draw the same re-paired values.
  figs <- build_arm_figs(fx$seg, tibble::tibble(), fx$man, fx$full)
  expect_true(all(figs[["04_pair_fraction_by_arm"]]$data$pair_fraction == 0.55))
  expect_true(all(figs[["02_final_dice_by_arm"]]$data$dice_matched == 0.71))
  expect_match(figs[["02_final_dice_by_arm"]]$labels$subtitle, "after the whole transform")
})

test_that("an arm scored before the record existed falls back to the ladder, labelled", {
  fx  <- stage_fixture(full = FALSE)
  fin <- arm_final_score(fx$seg, fx$full)
  expect_true(all(fin$scored_on == SCORED_LADDER))
  expect_equal(fin$dice_matched, arm_final_stage(fx$seg)$dice_matched)
  # One arm re-scored, the rest not: both kinds survive and the note says MIXED.
  fx2 <- stage_fixture()
  one <- dplyr::filter(fx2$full, arm_dir == "valis_high_micro2")
  mix <- arm_final_score(fx2$seg, one)
  expect_setequal(unique(mix$scored_on), c(SCORED_FULL, SCORED_LADDER))
  expect_match(scored_on_note(mix), "MIXED")
})

test_that("the ranking table carries a mean beside each median", {
  seg <- read_arms_seg_qc(arm_manifest(arms_tree()))
  fin <- arm_final_stage(seg)
  tbl <- arm_ranking_table(seg)
  one <- tbl[1, ]
  v   <- fin$disp_um_p50[fin$arm == one$arm]
  expect_equal(one$disp_um_p50, stats::median(v))
  expect_equal(one$disp_um_mean, mean(v))
  expect_true(all(c("dice_mean", "pair_fraction_mean") %in% names(tbl)))
  # Skewed slides: the two summaries must differ, so the mean is not a relabelled median.
  seg$disp_um_p50[seg$arm == one$arm & seg$patient_id == "046"] <- 50
  sk <- dplyr::filter(arm_ranking_table(seg), arm == one$arm)
  fs <- arm_final_stage(seg); w <- fs$disp_um_p50[fs$arm == one$arm]
  expect_equal(sk$disp_um_mean, mean(w))
  expect_equal(sk$disp_um_p50, stats::median(w))
})

# --- STARE's microns: converted from pixels with the slide's own pixel size -------
# Rewrite a fixture tree the way the real runs are written: VALIS records carry
# params.pixel_size_um; tiled records carry pixels and no micron value at all.
strip_tiled_um <- function(root, ps = 0.325) {
  for (f in fs::dir_ls(root, recurse = TRUE, glob = "*_seg_qc.json")) {
    d <- jsonlite::fromJSON(f, simplifyVector = FALSE)
    tiled <- grepl("/tiled_", f)
    d$params <- list(pixel_size_um = if (tiled) NULL else ps)
    if (tiled) {
      for (st in names(d$stages)) d$stages[[st]]$displacement_um_p50 <- NULL
      d$delta_vs_anchor <- NULL
    }
    jsonlite::write_json(d, f, auto_unbox = TRUE, null = "null")
  }
  root
}

test_that("a pixel-only STARE record gets microns from the same slide's pixel size", {
  root <- strip_tiled_um(arms_tree(modes = "high", depths = 2, tiled = TRUE))
  seg  <- suppressMessages(read_arms_seg_qc(arm_manifest(root)))
  tl   <- dplyr::filter(seg, backend == "tiled")
  expect_true(all(tl$um_from_px))
  expect_false(any(dplyr::filter(seg, backend == "valis")$um_from_px))
  # px x pixel size: the fixture wrote px = um / 0.325 (to write_json's 4 digits).
  expect_equal(tl$disp_um_p50[tl$stage == "refined"], rep(1.15, 2), tolerance = 1e-3)
  expect_message(read_arms_seg_qc(arm_manifest(root)), "converted")
  # STARE now has a place in the micron ranking.
  tbl <- arm_ranking_table(seg)
  expect_true(is.finite(tbl$disp_um_p50[tbl$backend == "tiled"]))
})

test_that("with no pixel size anywhere, nothing is invented", {
  root <- strip_tiled_um(arms_tree(modes = "high", depths = 2, tiled = TRUE))
  for (f in fs::dir_ls(root, recurse = TRUE, glob = "*_seg_qc.json")) {
    d <- jsonlite::fromJSON(f, simplifyVector = FALSE); d$params <- NULL
    jsonlite::write_json(d, f, auto_unbox = TRUE, null = "null")
  }
  expect_warning(seg <- read_arms_seg_qc(arm_manifest(root)), "no pixel size")
  expect_true(all(is.na(dplyr::filter(seg, backend == "tiled")$disp_um_p50)))
})
