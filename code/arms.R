# =============================================================================
# arms.R  —  the ARM REGISTRY: which files on disk belong to which arm, and how
# a filename names the region it holds.
#
# Three arms phenotype THE SAME SLIDES three different ways, so the diff between
# any two of them isolates a method rather than a cohort. Each publishes its cells
# and its pathologist polygons in its own tree, with its own file-naming
# convention, and this file is the ONE place that knows the mapping. Nothing
# downstream parses a path.
#
#   arm                 cells                                polygons
#   ------------------  -----------------------------------  ----------------------
#   massimo1            FlowPath_csv_selected/<pid>_a<k>.csv  annotation_selected/
#                       FlowPath_csv_all/<pid>/<pid>.csv        <pid>/<pid>_a<k>.geojson
#                                                            annotation_all/
#                                                              <pid>/<pid>.geojson
#   massimo1_inverted   csv_inverted-classification_          the SAME polygons as
#                         modified-thrPANCK/<pid>_a<k>.csv    massimo1 — it re-classifies
#                                                             arm 1's regions, it does
#                                                             not re-draw them
#   massimo2            csv/<pid>/<pid>_<A|B|C>.csv          annotation/
#                                                              <pid>/<pid>_<A|B|C>.geojson
#
# THE ARMS' REGIONS ARE INDEPENDENTLY DRAWN. massimo1's `a1` is NOT massimo2's `A`.
# They are separate annotation sessions over the same tissue, so ANNOTATION_<k> is
# an ARM-LOCAL label, not a cohort-wide identity. Two consequences, both enforced
# rather than documented:
#   - each arm carries its OWN pathologist table (neoplastic_massimo1 /
#     neoplastic_massimo2 in load_data.R), because a score keyed ANNOTATION_2
#     means a different polygon in each;
#   - any cross-arm figure joins on PATIENT, never on (patient, annotation).
#     Pairing arm 1's a2 with arm 2's B would produce a plausible correlation
#     between two regions that do not overlap.
# 24086 is the case that makes this concrete: three regions in massimo1, none at
# all in massimo2 (where it is a whole-slide export). Both are correct.
#
# WHY massimo1_inverted LIVES INSIDE massimo1's ROOT. It is a re-classification of
# arm 1's own regions at a modified PANCK threshold — it ships only the cells and
# borrows both polygon trees. Keeping it under `massimo1/` means one symlink
# supplies all three trees and the borrowed paths need no cross-root plumbing. It
# is still a full arm with its own mode string; only its storage is shared.
#
# TWO TIERS, TWO SCOPES. massimo1 publishes an `_all` tier (one whole-slide csv and
# one union polygon per patient) alongside its `_selected` regions, and the two ARE
# the two scopes `membership_data()` returns: `_selected` answers per_annotation,
# `_all` answers union. Arms with no `_all` tier derive the union the way
# arm_metrics() always has — pool the region files, de-duplicate, dissolve the
# polygons. `union_csv`/`union_poly` being NULL is what says "derive it".
#
# Depends on nothing but fs + base R, so it sits beside cell_tables.R at the bottom
# of the stack and can be tested without sf, DESeq2 or any data on disk.
# =============================================================================
suppressPackageStartupMessages({
  library(fs)
})

# --- Filename -> (patient, region) -------------------------------------------
# Four naming conventions across the three arms. Each parser answers the same
# question — which patient and which region is this file? — and every one of them
# lets the DIRECTORY NAME WIN over the filename stem where a directory exists. A
# file renamed by hand would otherwise be keyed to the wrong patient silently, and
# the only symptom is an implausible in-annotation count.
#
#   nested_letter  <root>/<pid>/<pid>_<A|B|C>.ext   letter, by alphabet POSITION
#   nested_digit   <root>/<pid>/<pid>_a<k>.ext      digit k, taken literally
#   nested_bare    <root>/<pid>/<pid>.ext           the patient's single file
#   flat_digit     <root>/<pid>_a<k>.ext            no directory: stem is the key
#   flat_letter    <root>/<pid>_<A|B|C>.ext         the same, letter-suffixed
#
# No arm uses `flat_letter` today. It is here because the axis is nested/flat x
# letter/digit/bare and a registry that describes five of six cases invites the sixth
# to be re-implemented somewhere else — which is exactly how validation_helpers.R
# ended up with a second copy of this parser. .annotation_key() delegates here.
#
# `annotation` comes back as ANNOTATION_<k>, or NA for a bare file — the CALLER
# decides what a bare file means, because it differs by tier: a bare csv in
# massimo2's region tier is a whole-slide export, while a bare geojson in
# massimo1's `_all` tier is that patient's union polygon. Deciding here would
# force one of those two readings onto the other.
ARM_NAME_PATTERNS <- c("nested_letter", "nested_digit", "nested_bare",
                       "flat_digit", "flat_letter")

# "A" -> 1, by alphabet position, not by file order: "C" is region 3 whether or
# not "B" was ever exported.
arm_letter_index <- function(letter) match(toupper(letter), LETTERS)

arm_parse_name <- function(path, pattern = ARM_NAME_PATTERNS) {
  pattern <- match.arg(pattern)
  stem    <- fs::path_ext_remove(fs::path_file(path))
  parent  <- fs::path_file(fs::path_dir(path))

  none <- function(pid) list(patient = pid, annotation = NA_character_)

  if (pattern == "flat_digit") {
    m <- regmatches(stem, regexec("^(.*)_a([0-9]+)$", stem))[[1]]
    if (length(m) == 3)
      return(list(patient = m[2], annotation = paste0("ANNOTATION_", as.integer(m[3]))))
    return(none(stem))
  }

  if (pattern == "flat_letter") {
    m <- regmatches(stem, regexec("^(.*)_([A-Za-z])$", stem))[[1]]
    if (length(m) == 3) {
      idx <- arm_letter_index(m[3])
      if (!is.na(idx)) return(list(patient = m[2], annotation = paste0("ANNOTATION_", idx)))
    }
    return(none(stem))
  }

  # Every remaining pattern is nested, so the parent directory is the patient.
  if (pattern == "nested_bare") return(none(parent))

  if (pattern == "nested_digit") {
    m <- regmatches(stem, regexec("^(.*)_a([0-9]+)$", stem))[[1]]
    if (length(m) == 3)
      return(list(patient = parent, annotation = paste0("ANNOTATION_", as.integer(m[3]))))
    return(none(parent))
  }

  # nested_letter. The trailing character must be a LETTER; "046_a1" ends in a
  # digit and correctly falls through to the bare reading rather than being
  # mis-parsed, which is what keeps one parser from silently eating another arm's
  # files if a tree is ever pointed at the wrong spec.
  m <- regmatches(stem, regexec("^(.*)_([A-Za-z])$", stem))[[1]]
  if (length(m) == 3) {
    idx <- arm_letter_index(m[3])
    if (!is.na(idx)) return(list(patient = parent, annotation = paste0("ANNOTATION_", idx)))
  }
  none(parent)
}

# --- The registry ------------------------------------------------------------
# One entry per arm. `root` is relative to data/; a tier is a (dir, pattern) pair
# or NULL, and NULL on a union tier means "derive it from the region tier".
#
# `bare_region_is` is the one genuinely arm-specific rule and it is spelled out
# rather than inferred: it says what a region-tier file with no region suffix
# means. In massimo2 that is a patient with no annotation drawn at all, whose
# whole slide counts as one region ("whole_slide" — the layout's own stated
# convention). No other arm has such a file, so the field is NA for them and a
# stray bare csv is a warning instead of a silent extra region.
ARM_SPECS <- list(
  massimo1 = list(
    arm            = "massimo1",
    root           = "massimo1",
    label          = "FlowPath, selected regions",
    region_csv     = list(dir = "FlowPath_csv_selected", pattern = "flat_digit"),
    region_poly    = list(dir = "annotation_selected",   pattern = "nested_digit"),
    union_csv      = list(dir = "FlowPath_csv_all",      pattern = "nested_bare"),
    union_poly     = list(dir = "annotation_all",        pattern = "nested_bare"),
    bare_region_is = NA_character_
  ),
  massimo1_inverted = list(
    arm            = "massimo1_inverted",
    root           = "massimo1",
    label          = "FlowPath, inverted classification at a modified PANCK threshold",
    region_csv     = list(dir = "csv_inverted-classification_modified-thrPANCK",
                          pattern = "flat_digit"),
    # Borrowed from massimo1: same regions, re-classified cells.
    region_poly    = list(dir = "annotation_selected", pattern = "nested_digit"),
    union_csv      = NULL,
    union_poly     = list(dir = "annotation_all",     pattern = "nested_bare"),
    bare_region_is = NA_character_
  ),
  massimo2 = list(
    arm            = "massimo2",
    root           = "massimo2",
    label          = "the all-slide export, letter-suffixed regions",
    region_csv     = list(dir = "csv",        pattern = "nested_letter"),
    region_poly    = list(dir = "annotation", pattern = "nested_letter"),
    union_csv      = NULL,
    union_poly     = NULL,
    bare_region_is = "whole_slide"
  )
)

ARM_MODES <- names(ARM_SPECS)

# DIRECTORY LOOKUP IS CASE-INSENSITIVE, DELIBERATELY.
#
# The producer ships `Massimo1`/`Massimo2`; the registry spells them lowercase to
# match the mode strings. macOS resolves that mismatch silently because its
# filesystem is case-insensitive, so a tree copied as `data/Massimo1` loads fine on
# a laptop and finds NOTHING on the Linux cluster — where the whole site then
# renders empty with no error, because an absent tree is a legitimate state.
#
# That failure is invisible exactly where it matters, so the name is resolved
# against what is actually on disk rather than assumed. An exact match always wins;
# a unique case-insensitive match is accepted and is what makes `Massimo1` and
# `massimo1` the same tree on both platforms. An AMBIGUOUS match (both spellings
# present, which only a case-sensitive filesystem can even represent) falls back to
# the exact name rather than guessing which the user meant.
.arm_resolve_dir <- function(parent, name) {
  direct <- file.path(parent, name)
  if (dir.exists(direct) || !dir.exists(parent)) return(direct)
  kids <- list.dirs(parent, full.names = FALSE, recursive = FALSE)
  hit  <- kids[tolower(kids) == tolower(name)]
  if (length(hit) == 1) file.path(parent, hit) else direct
}

# The spec, with its tier directories resolved to absolute paths. Every reader
# takes THIS rather than a root string, so no function downstream joins a path.
arm_spec <- function(arm = ARM_MODES, data_dir = here::here("data")) {
  arm  <- match.arg(arm)
  spec <- ARM_SPECS[[arm]]
  root <- .arm_resolve_dir(data_dir, spec$root)
  spec$root_path <- root
  for (tier in c("region_csv", "region_poly", "union_csv", "union_poly"))
    if (!is.null(spec[[tier]]))
      spec[[tier]]$path <- .arm_resolve_dir(root, spec[[tier]]$dir)
  spec
}

# Does this arm ship a dedicated whole-slide tier, or must the union be derived by
# pooling and de-duplicating the region files? Asked in several places, so it is
# named once rather than being an `is.null()` test scattered around.
arm_has_union_tier <- function(spec) !is.null(spec$union_csv)

# Which tiers are actually on disk. Reported by every page's provenance table:
# an arm whose region tier is missing yields empty panels, and the reader should
# be able to say which directory it looked in rather than "no cells".
arm_tier_status <- function(spec) {
  tiers <- c("region_csv", "region_poly", "union_csv", "union_poly")
  tibble::tibble(
    arm    = spec$arm,
    tier   = tiers,
    dir    = vapply(tiers, function(t) if (is.null(spec[[t]])) NA_character_ else spec[[t]]$dir,
                    character(1)),
    exists = vapply(tiers, function(t) !is.null(spec[[t]]) && dir.exists(spec[[t]]$path),
                    logical(1)))
}
