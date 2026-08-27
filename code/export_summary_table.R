# =============================================================================
# export_summary_table.R  —  the cross-source SUMMARY EXPORT, for analysis
#                            outside this repo (SPSS, Prism, a colleague's script)
#
#   Rscript code/export_summary_table.R                 # -> output/
#   Rscript code/export_summary_table.R --outdir /path  # somewhere else
#
# Four quantities that live in four places and at THREE DIFFERENT GRAINS:
#
#   bulk RNA immune markers          per bulk-RNA sample, keyed in CRF space
#   neoplastic cellularity           per (arm, patient, ANNOTATION_k) — arm-local
#   in-annotation cell fractions     per (arm, patient, ANNOTATION_k) — arm-local
#   IHC immune score + response      per patient
#
# So it writes TWO files rather than one, and which grain a number belongs to is
# the whole point of the split:
#
#   summary_patient.csv   one row per patient. Everything collapsed to the
#                         patient, cross-arm columns side by side.
#   summary_region.csv    one row per (arm, patient, annotation). LONG ON ARM,
#                         never wide — see below.
#
# WHY THE REGION FILE IS LONG ON ARM AND THE PATIENT FILE IS WIDE. massimo1's
# ANNOTATION_2 and massimo2's ANNOTATION_2 are DIFFERENT POLYGONS: the arms were
# annotated in separate sessions over the same tissue (code/arms.R). A wide region
# file would put two unrelated polygons on one row and invite exactly the join this
# project forbids. At patient level the arms ARE comparable — a patient is a
# patient — so there the columns sit side by side and carry the arm in their name.
#
# The two files join on `patient_id`. RNA and the immune score are deliberately
# NOT repeated down the region file: a patient with three regions would then vote
# three times in any regression fitted on the region grain.
#
# NOTHING HERE COMPUTES A QUANTITY OF ITS OWN. Every number comes from the same
# function the website's pages call — region_ratios() via arm_metrics(),
# ihc_marker_fraction(), asin_sqrt(), the marker_gene_map averages — so a value in
# this csv and the corresponding point in a figure cannot drift apart. This file
# only selects, renames and joins.
#
# SCOPES ARE MIXED ON PURPOSE, and the column names say so rather than leaving it
# to be inferred: `*_over_inside` is inside the tumour annotation, `*_wholeslide`
# is every imaged cell with no polygon consulted. Same protein, different
# denominator; the two are not interchangeable.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(here)
})

# --- Arguments ---------------------------------------------------------------
.args  <- commandArgs(trailingOnly = TRUE)
.opt   <- function(flag, default) {
  i <- match(flag, .args)
  if (is.na(i) || i == length(.args)) default else .args[i + 1L]
}
OUTDIR <- .opt("--outdir", here::here("output"))

# --- Preflight ---------------------------------------------------------------
# NAME THE MISSING INPUT, AND STOP. The alternative — carrying on and writing the
# columns that happen to be reachable — produces a csv whose NA columns are
# indistinguishable from a genuine measurement failure, and it is exactly the
# shape a collaborator would analyse without noticing. An absent tree is a
# legitimate state for a *page* (it renders empty panels and says so); it is not a
# legitimate state for an export whose whole purpose is to leave this repo.
.required <- c(
  "counts.RData"       = here::here("data", "counts.RData"),
  "clinical_data.xlsx" = here::here("data", "clinical_data.xlsx")
)
.missing <- .required[!file.exists(.required)]
# The arm trees are matched case-insensitively (the producer ships `Massimo1`,
# the registry spells it lowercase) — arm_spec() does that resolution, so ask it
# rather than testing a hardcoded path that is right on macOS and wrong on Linux.
source(here::here("code", "arms.R"))
for (a in c("massimo1", "massimo2")) {
  p <- arm_spec(a)$region_csv$path
  if (!dir.exists(p)) .missing[[paste0("arm ", a)]] <- p
}
if (length(.missing))
  stop("export_summary_table.R: missing input(s)\n",
       paste0("  ", names(.missing), "  ->  ", .missing, collapse = "\n"),
       "\n  symlink the cluster trees in:  ln -s <share>/Massimo1 data/massimo1",
       call. = FALSE)

# --- Load --------------------------------------------------------------------
# Same order, and the same capture.output() wrapper, as every analysis page: the
# loaders emit a LOADED/whole-slide-check message per arm that belongs in the
# cluster log, not interleaved with a progress bar.
message("loading...")
invisible(capture.output(suppressMessages(source(here::here("code", "load_data.R")))))
source(here::here("code", "validation_helpers.R"))
source(here::here("code", "membership.R"))

ARMS <- c("massimo1", "massimo2")

# --- Membership, once per arm ------------------------------------------------
# membership_data() reads the arm's OWN polygons — passing an annotation set is an
# error, not an ignored argument, because it would score one arm's cells against
# another arm's regions.
mem <- lapply(ARMS, function(a) {
  message("resolving membership: ", a)
  membership_data(a, ihc_for(a), neoplastic_for(a))
})
names(mem) <- ARMS

for (a in ARMS)
  if (nrow(mem[[a]]$cells) == 0)
    stop("arm ", a, " loaded no cells — every fraction would be NA. ",
         "Check ", arm_spec(a)$root_path, call. = FALSE)

# =============================================================================
# The five in-annotation fractions
# =============================================================================
# The requested set, in the order the immune lineage reads: total immune (CD45+),
# marker-gated T cells (CD3+CD45+), then the three phenotype lineages. All five are
# `n_<population> / n_inside` over the cells inside the tumour annotation, so all
# five are true 0..1 proportions and asin_sqrt() applies to every one of them.
#
# These are region_ratios() column names, NOT new definitions — `cd45_over_inside`
# and `frac_CD4T` are computed once, in validation_helpers.R, and merely selected
# here. Renaming them in this file would create a second vocabulary for the same
# quantity, which is how the project ended up with two copies of the filename
# parser it now has one of.
FRACTION_COLS <- c("cd45_over_inside", "cd3cd45_over_inside",
                   "frac_CD4T", "frac_CD8T", "frac_Treg")

# Each fraction plus its arcsine square-root twin. asin_sqrt() (validation_helpers.R)
# clamps to [0,1] and returns NA outside it, so a column that is not a proportion
# comes back empty rather than silently wrong.
add_asin <- function(df, cols = FRACTION_COLS) {
  cols <- intersect(cols, names(df))
  for (cc in cols) df[[paste0(cc, "_asin")]] <- asin_sqrt(df[[cc]])
  df
}

# =============================================================================
# FILE 2 — one row per (arm, patient, annotation)
# =============================================================================
# Built first because the patient file reuses its pathologist aggregation.
#
# `annotation` is ANNOTATION_<k> for a scored region and `whole_slide` for a
# patient with no polygon drawn at all (arm 2's 24086). Those are not the same
# statement and the label keeps them apart; `source` says how membership was
# actually decided — `sf` from the polygon, `flag` from the export's own
# out-of-annotation column, `whole_slide` when the arm's convention is that every
# cell counts.
region_tbl <- purrr::map_dfr(ARMS, function(a) {
  path_long <- neoplastic_for(a) |>
    mutate(patient_id = slide_key(SAMPLE)) |>
    filter(!is.na(path_pct)) |>
    select(patient_id, annotation, path_pct)

  mem[[a]]$per_annotation |>
    # LEFT join, not inner: a region the pathologist did not score still has real
    # cell fractions, and dropping it would quietly shrink the denominator of any
    # count taken off this file. The empty path_pct says which.
    left_join(path_long, by = c("patient_id", "annotation")) |>
    transmute(arm = a, patient_id, annotation, source,
              path_pct,
              n_inside, tumor_over_inside,
              across(all_of(intersect(FRACTION_COLS, names(mem[[a]]$per_annotation))))) |>
    add_asin() |>
    arrange(patient_id, annotation)
})

# =============================================================================
# FILE 1 — one row per patient
# =============================================================================

# --- Bulk RNA ----------------------------------------------------------------
# Two readings of the same matrix, because they answer different questions:
#   rna_<marker>   the marker_gene_map average — CD3 is the mean of CD3D/CD3E/CD3G
#                  and CD8 of CD8A/CD8B. This is the value molecular_*.Rmd plots.
#   gene_<SYMBOL>  the raw per-gene value, so an average hiding one discordant gene
#                  is visible rather than assumed away.
# The `gene_` prefix is not decoration: the marker CD4 and the gene CD4 are the
# same string, and two columns cannot share a name.
RNA_MARKERS <- c("CD45", "CD3", "CD4", "CD8", "FOXP3")

crf2slide <- crf_to_slide_map(clinical_data)

expr_long <- counts_data |>
  mutate(crf_id = norm_id(Sample)) |>
  pivot_longer(-c(Sample, crf_id), names_to = "gene", values_to = "expr")

# marker -> its genes, one row per (marker, gene). Explicit Map()/strsplit rather
# than separate_rows(), matching molecular_body.Rmd.
gene_rows <- do.call(rbind, Map(function(m, g)
  data.frame(marker = m, gene = trimws(strsplit(g, ",")[[1]]), stringsAsFactors = FALSE),
  marker_gene_map$marker, marker_gene_map$genes))
gene_rows <- gene_rows[gene_rows$marker %in% RNA_MARKERS, , drop = FALSE]

# A marker's genes that are actually columns of the RNA matrix. A gene in the map
# but absent from the matrix is REPORTED and excluded from the mean — averaging
# over a gene that was never measured would pull the marker toward zero.
absent <- setdiff(gene_rows$gene, unique(expr_long$gene))
if (length(absent))
  message("genes absent from the RNA matrix (excluded from the marker means): ",
          paste(absent, collapse = ", "))
marker_genes <- gene_rows[gene_rows$gene %in% unique(expr_long$gene), , drop = FALSE]

# One value per (marker, patient): mean over that marker's AVAILABLE genes, then
# CRF space -> slide space. crf_to_slide_map() is a named vector rather than a
# join because dplyr join column-propagation is unreliable on this project's build
# (see validation_helpers.R).
rna_marker <- expr_long |>
  inner_join(marker_genes, by = "gene", relationship = "many-to-many") |>
  group_by(marker, crf_id) |>
  summarise(rna = mean(expr), .groups = "drop") |>
  mutate(patient_id = unname(crf2slide[crf_id])) |>
  filter(!is.na(patient_id)) |>
  mutate(marker = paste0("rna_", marker)) |>
  pivot_wider(id_cols = c(patient_id, crf_id), names_from = marker, values_from = rna)

rna_gene <- expr_long |>
  filter(gene %in% marker_genes$gene) |>
  mutate(patient_id = unname(crf2slide[crf_id])) |>
  filter(!is.na(patient_id)) |>
  group_by(patient_id, gene) |>
  # A patient with two RNA samples is averaged, and the sample count is carried
  # out on the patient table so an averaged row is never mistaken for a single one.
  summarise(expr = mean(expr), .groups = "drop") |>
  mutate(gene = paste0("gene_", gene)) |>
  pivot_wider(names_from = gene, values_from = expr)

# The CRF ids a patient's RNA came from, and how many. Both are carried out of
# expr_long HERE because the patient-level mean below groups by `patient_id` and
# would otherwise drop the CRF column entirely — leaving the file with no record
# of which bulk samples produced its RNA values, and a `n_rna_samples` of 2
# indistinguishable from a 1 unless it is stated.
rna_n <- expr_long |>
  mutate(patient_id = unname(crf2slide[crf_id])) |>
  filter(!is.na(patient_id)) |>
  distinct(patient_id, Sample) |>
  group_by(patient_id) |>
  summarise(crf_samples   = paste(sort(Sample), collapse = ";"),
            n_rna_samples = dplyr::n(), .groups = "drop")

rna_tbl <- rna_marker |>
  group_by(patient_id) |>
  summarise(across(starts_with("rna_"), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  left_join(rna_gene, by = "patient_id") |>
  left_join(rna_n,    by = "patient_id")

# --- Pathologist cellularity, per arm, collapsed to the patient ---------------
# UNWEIGHTED mean of the patient's scored regions — the same collapse the site's
# union panels use (`path_union` in clinical_body.Rmd). It is a choice, not the
# only one: a plain mean lets a small region and a large one vote equally, which
# is why `n_regions_<arm>` ships beside every value. A patient scored whole_slide
# has exactly one row and the mean is that row.
path_patient <- function(a) {
  neoplastic_for(a) |>
    mutate(patient_id = slide_key(SAMPLE)) |>
    filter(!is.na(path_pct)) |>
    group_by(patient_id) |>
    summarise(path_pct = mean(path_pct), n_regions = dplyr::n(), .groups = "drop") |>
    rename_with(~ paste0(.x, "_", a), c(path_pct, n_regions))
}

# --- The union scope, per arm ------------------------------------------------
# One row per patient: the cells inside the DISSOLVED union of that patient's
# polygons, de-duplicated so a cell falling in two overlapping regions counts once.
#
# `mem$union` IS TAKEN WHOLE — deliberately not filtered to `annotation == "union"`.
# A patient whose arm has no annotation drawn for it (massimo2's 24086) carries a
# union row labelled `whole_slide` instead, and under that arm's own convention
# (`bare_region_is`) the whole slide IS its annotated region — so that row is its
# annotation_all value and belongs in these columns. Filtering on the label would
# empty every massimo2_* cell for that patient while the rest of the row stayed
# populated, which reads as a failed measurement rather than as a patient nobody
# annotated. code/scope_compare.R's .scope_union_labels() states the same rule for
# the site's three-scope figures.
union_tbl <- function(a, prefix, cols) {
  cols <- intersect(cols, names(mem[[a]]$union))
  out  <- mem[[a]]$union |>
    select(patient_id, n_inside, tumor_over_inside, all_of(cols)) |>
    # `tumor_over_inside` is asin'd on EVERY arm, not merely on whichever call
    # happened to list it: leaving it to the caller's column list gave arm 1 a
    # transformed tumour fraction and arm 2 none, an asymmetry with no meaning
    # that a reader would reasonably take for one.
    add_asin(unique(c(cols, "tumor_over_inside")))
  # Prefix everything but the key, so two arms' columns can sit on one row and
  # each still says which arm it came from.
  rename_with(out, ~ paste0(prefix, .x), -patient_id)
}

# --- The IHC immune score (arm 1), both scopes -------------------------------
# `_wholeslide` reproduces the site's "IHC immune score by induction response"
# figure exactly: CD45_posfrac over EVERY imaged cell of the patient, no polygon
# consulted (clinical_body.Rmd). `_inside` is the same protein over the cells
# inside the tumour union. They are different denominators and will not agree;
# both are here so the difference is visible instead of being a choice made
# silently upstream.
#
# marker_pos() counts a cell no gate evaluated as NOT positive, so a sparsely
# gated marker reads low rather than dropping out — check marker_gated() before
# reading a near-zero score as biology.
immune_ws <- ihc_marker_fraction(mem$massimo1$cells) |>
  select(patient_id, immune_score_massimo1_wholeslide = CD45_posfrac)

immune_in <- mem$massimo1$union |>
  select(patient_id, immune_score_massimo1_inside = cd45_over_inside)

# --- Clinical ----------------------------------------------------------------
# RESPONSE is the induction-treatment response the immune score is grouped by.
clin_tbl <- clinical_data |>
  mutate(patient_id = slide_key(`ID PATIENT`)) |>
  select(patient_id, response = any_of("RESPONSE")) |>
  distinct()

# --- Assemble ----------------------------------------------------------------
# FULL joins across the MEASUREMENTS, so the spine is the union of every measured
# source's patients: a patient with cells but no RNA (or the reverse) is a fact
# about the cohort and belongs in the file as a row with holes. An inner join
# would delete it and the file would describe a cohort that was never assembled.
#
# THE CLINICAL TABLE IS JOINED ON, NEVER JOINED WITH. It covers the whole study,
# not the handful of patients that were imaged and sequenced, so full-joining it
# would turn a 6-patient measurement file into a ~24-row one whose extra rows
# carry a `response` and nothing else. Those rows are not partial observations —
# nothing was measured for them here — and they would inflate every n taken off
# this file while each row still looked well-formed. A left join keeps `response`
# an ATTRIBUTE of a measured patient, which is what it is.
patient_tbl <- list(
  rna_tbl,
  path_patient("massimo1"), path_patient("massimo2"),
  union_tbl("massimo1", "massimo1_", "tumor_over_inside"),
  union_tbl("massimo2", "massimo2_", FRACTION_COLS),
  immune_ws, immune_in
) |>
  purrr::reduce(full_join, by = "patient_id") |>
  filter(!is.na(patient_id)) |>
  left_join(clin_tbl, by = "patient_id") |>
  arrange(patient_id)

# Column order is READING order, not the order the joins happened to produce:
# identity, then RNA, then the pathologist's two arms, then arm 2's fractions,
# then the immune score and the clinical outcome.
patient_tbl <- patient_tbl |>
  select(patient_id, any_of(c("crf_samples", "n_rna_samples")),
         starts_with("rna_"), starts_with("gene_"),
         starts_with("path_pct_"), starts_with("n_regions_"),
         starts_with("massimo1_"), starts_with("massimo2_"),
         starts_with("immune_score_"), any_of("response"),
         everything())

# --- Write -------------------------------------------------------------------
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
p_path <- file.path(OUTDIR, "summary_patient.csv")
r_path <- file.path(OUTDIR, "summary_region.csv")
readr::write_csv(patient_tbl, p_path)
readr::write_csv(region_tbl,  r_path)

# --- Say what was written ----------------------------------------------------
# The cluster log is the only thing a reader has when the job is over, so it
# states the n of each file and, per column, how many values are actually present.
# A column that joined to nothing is a silent all-NA column in the csv and an
# obvious 0/6 here.
cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("%-24s %d rows x %d cols  ->  %s\n",
            "summary_patient.csv", nrow(patient_tbl), ncol(patient_tbl), p_path))
cat(sprintf("%-24s %d rows x %d cols  ->  %s\n",
            "summary_region.csv",  nrow(region_tbl),  ncol(region_tbl),  r_path))
cat(strrep("=", 70), "\n\n", sep = "")

cat("patient file — non-missing values per column (of ", nrow(patient_tbl), " patients):\n", sep = "")
filled <- vapply(patient_tbl, function(x) sum(!is.na(x)), integer(1))
for (nm in names(filled))
  cat(sprintf("  %-42s %d\n", nm, filled[[nm]]))

cat("\nregion file — rows per arm:\n")
print(as.data.frame(dplyr::count(region_tbl, arm, source, name = "rows")))
cat("\nregions with no pathologist score:",
    sum(is.na(region_tbl$path_pct)), "of", nrow(region_tbl), "\n")
