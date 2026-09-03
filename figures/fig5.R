#!/usr/bin/env Rscript
# =============================================================================
# figures/fig5.R  —  Figure 5. Proof-of-concept concordance of MIRAGE + FlowPath
#                    imaging immune quantification with two orthogonal
#                    transcriptomic proxies, six head-and-neck cases.
#
#   (a) Registered, phenotyped cases: whole slide, tumour vs rest, beside
#       the tumour-richest region coloured by compartment              [computed]
#   (b) mIF CD45+/all-cells, three hot vs three cold                 [computed]
#   (c) quanTIseq deconvolution vs imaging fraction, by population   [computed]
#
# As with fig4.R, nothing here is new science: (a)-(c) are paper_phenotype_map(),
# paper_immune_fraction_hotcold() and paper_deconv_scatter() out of
# code/paper_figures.R, the same objects analysis/paper_figures.Rmd prints. This
# script selects, strips, lays out and exports.
#
# WHICH CASES PANEL (a) SHOWS is not hard-coded, and that is deliberate. The legend
# marks it "[author to supply]", and the honest answer is "one representative hot and
# one representative cold" — which means the SAME hot/cold axis panel (b) groups on.
# Deriving it here rather than pasting two patient IDs means (a) and (b) cannot drift:
# change the clinical table and both move together. Override with FIG5A_CASES if a
# specific pair is wanted for image quality reasons.
#
# NO TEST AND NO COEFFICIENT ON (b) OR (c). That is the legend's claim and it is
# enforced in the builders, not here — paper_immune_fraction_hotcold() deliberately
# never calls paired_spearman(), and paper_deconv_scatter() draws no fit line. Do
# not add either at assembly time: (c)'s axes use different denominators (imaging
# counts cells, deconvolution estimates a bulk-mixture fraction), so a regression
# line invites exactly the absolute-agreement reading the legend disclaims. Ranking
# is the claim.
#
# Run:  Rscript figures/fig5.R
# =============================================================================

source(file.path(tryCatch(here::here(), error = function(e) normalizePath(".")),
                 "figures", "_common.R"))
suppressPackageStartupMessages({
  library(dplyr); library(tibble)
})

root <- here_root
source(file.path(root, "code", "validation_helpers.R"))
source(file.path(root, "code", "paper_figures.R"))
source(file.path(root, "code", "arm_cells.R"))

# Set to a character vector of patient ids to pin panel (a); NULL derives them.
FIG5A_CASES    <- NULL
# "hot_score" ranks the continuous score and takes the top/bottom k. "immuno_phe"
# uses the clinical category instead. Same switch, same meaning, as the Rmd.
HOTCOLD_SOURCE <- "hot_score"

# --- Load --------------------------------------------------------------------
# ONE arm, the same one analysis/paper_figures.Rmd draws from. Regions are arm-local
# (see CLAUDE.md), so this cannot pool arms.
ARM     <- "massimo2"
as_spec <- arm_spec(ARM)
if (!dir.exists(as_spec$region_csv$path))
  stop("fig5: no cells under ", as_spec$region_csv$path, ". Expected\n",
       "  data/", ARM, "/csv/<patient>/<patient>_<A|B|C>.csv\n",
       "  data/", ARM, "/annotation/<patient>/<patient>_<A|B|C>.geojson")
as_cells  <- arm_cells(as_spec)
as_ucells <- arm_union_tier_cells(as_spec)
if (!nrow(as_cells)) stop("fig5: arm_cells() read nothing for ", ARM)

# arm_annotations() errors without sf on purpose — a caller that asked for an
# outline should hear that it cannot be drawn. The map is still correct without one,
# so this degrades to NULL rather than aborting the figure.
pids     <- unique(c(as_cells$patient_id, as_ucells$patient_id))
as_polys <- tryCatch(arm_annotations(as_spec, "region", patient_ids = pids),
  error = function(e) { warning("fig5: no annotation outlines (", conditionMessage(e),
                                ")", call. = FALSE); NULL })
as_upoly <- tryCatch(arm_annotations(as_spec, "union", patient_ids = pids),
                     error = function(e) NULL)
.prom    <- arm_promote_unregioned(as_spec, as_cells, as_polys, as_ucells, as_upoly)
as_cells <- .prom$cells; as_polys <- .prom$polys
as_per   <- arm_metrics(as_spec, as_cells, as_polys, "per_annotation")
as_union <- arm_metrics(as_spec, as_cells, as_polys, "union",
                        union_cells = as_ucells, union_polys = as_upoly)

# --- The hot/cold axis, shared by (a) and (b) --------------------------------
groups    <- NULL
clin_path <- file.path(root, "data", "clinical_data.xlsx")
if (file.exists(clin_path)) {
  # select(any_of(...)) so a clinical table missing either column loses the column
  # rather than aborting — same idiom as molecular_hot_cold.Rmd.
  clin <- readxl::read_excel(clin_path) |>
    filter(!is.na(`ID PATIENT`)) |>
    mutate(patient_id = slide_key(`ID PATIENT`)) |>
    select(patient_id,
           hot_score  = any_of("HOT score"),
           immuno_phe = any_of("Immuno-phenotype")) |>
    filter(patient_id %in% slide_key(as_union$patient_id))
  if ("hot_score" %in% names(clin))
    clin$hot_score <- suppressWarnings(as.numeric(clin$hot_score))

  if (HOTCOLD_SOURCE == "hot_score" && "hot_score" %in% names(clin) &&
      sum(is.finite(clin$hot_score)) >= 2) {
    k <- min(3, floor(sum(is.finite(clin$hot_score)) / 2))
    r <- rank(-clin$hot_score, ties.method = "first")
    groups <- clin |>
      mutate(group = case_when(r <= k ~ "hot",
                               r > n() - k ~ "cold",
                               TRUE ~ NA_character_)) |>
      filter(!is.na(group)) |> select(patient_id, group)
  } else if ("immuno_phe" %in% names(clin)) {
    groups <- clin |> transmute(patient_id, group = immuno_phe) |> filter(!is.na(group))
  }
}
if (is.null(groups))
  stop("fig5: panel (b) needs data/clinical_data.xlsx for the hot/cold axis.")

# --- (a) One representative case per group -----------------------------------
cases <- FIG5A_CASES %||% {
  avail <- slide_key(unique(as_cells$patient_id))
  g     <- filter(groups, slide_key(patient_id) %in% avail)
  # First of each group under the ranking already applied above, so "representative"
  # means "the most extreme case that has imaging", not an arbitrary pick.
  c(head(g$patient_id[g$group == "hot"],  1),
    head(g$patient_id[g$group == "cold"], 1))
}
cases <- cases[!is.na(cases)]
if (!length(cases)) stop("fig5: no case is both grouped and present in the imaging.")

# One ROW per case: the whole slide (tumour vs everything else, every annotation
# outlined) beside the tumour-richest annotation cut to its box and coloured by
# compartment. paper_phenotype_map_pair() returns the two panels; the layout is here.
rows <- lapply(cases, function(pid) {
  pr <- paper_phenotype_map_pair(as_cells, pid, annots = as_polys, per = as_per)
  if (is.null(pr$overview)) return(NULL)
  ov <- for_panel(pr$overview)
  if (is.null(pr$inset)) return(wrap_elements(full = ov))
  ov | for_panel(pr$inset)
})
rows <- rows[!vapply(rows, is.null, logical(1))]
if (!length(rows)) stop("fig5: paper_phenotype_map_pair() returned nothing for ",
                        paste(cases, collapse = ", "))

# wrap_elements() so the case rows read as ONE tagged panel. Without it, tag_levels
# tags each map separately and the real (b) becomes (c).
# The two colour keys (tumour / non-tumour on the left, the four compartments on the
# right) are pinned UNDER the maps rather than left to default. The default put
# them top-left, where they collided with the "(a)" tag — and they belong to panel
# (a) alone in any case: (b) and (c) do not use these colours, so hoisting them into
# the figure-wide strip would file a key under a figure two thirds of which never
# refers to it.
p5a <- wrap_elements(full = Reduce(`/`, rows) +
                       plot_layout(guides = "collect") &
                       theme(legend.position = "bottom"))

# --- (b) and (c) -------------------------------------------------------------
p5b <- for_panel(paper_immune_fraction_hotcold(as_union, groups))
if (is.null(p5b)) stop("fig5: panel (b) empty after joining groups to metrics.")

paired_path <- file.path(root, "output", "paired_deconv.rds")
if (!file.exists(paired_path))
  stop("fig5: panel (c) needs output/paired_deconv.rds. Knit ",
       "analysis/molecular_hot_cold.Rmd first — its last chunk caches the paired ",
       "frame so this panel does not re-run immunedeconv.")
p5c <- for_panel(paper_deconv_scatter(readRDS(paired_path), method = "quantiseq"))
if (is.null(p5c)) stop("fig5: no quantiseq rows in ", paired_path)

save_panel(p5a, "p5a"); save_panel(p5b, "p5b"); save_panel(p5c, "p5c")

# --- Assemble ----------------------------------------------------------------
# (a) is one row per case, two maps wide, so it takes the larger share of the
# height. (b) is one axis with six points and (c) one shared-axis cloud, so the
# bottom row is split 1:2 rather than evenly — an equal split gives (b) whitespace
# it does not use and squeezes (c) below the width where the ranking is readable.
fig5 <- p5a / (p5b | p5c) +
  plot_layout(heights = c(1.6, 1), widths = c(1, 2), guides = "collect") +
  plot_annotation(tag_levels = TAG$tag_levels,
                  tag_prefix = TAG$tag_prefix, tag_suffix = TAG$tag_suffix) &
  theme(plot.tag = element_text(face = "bold"), legend.position = "bottom")

export_figure(fig5, "Fig5", width_mm = MM[["two_col"]], height_mm = 225)

message("fig5: panel (a) cases = ", paste(cases, collapse = ", "),
        " | groups = ", nrow(groups))
