# Code

Shared R sourced by the analyses in `analysis/`, plus standalone scripts.

## What is where

| file | owns |
|---|---|
| `load_data.R` | the raw loaders (`dds`, `clinical_data`, `counts_data`, and the **per-arm** `neoplastic_massimo1/2` + `ihc_massimo1/2/1_inverted`) |
| `cell_tables.R` | the **single-cell export schema** — one vocabulary over three upstream formats |
| `validation_helpers.R` | the derived quantities (region ratios, composition, marker/lineage fractions, invasive margin, agreement stats) |
| `mirage_cells.R` | the mirage cell source — `phenotypes.csv` + `morphology.csv` joined per patient |
| `arms.R` | the **arm registry** — which files belong to which arm, and how a filename names its region |
| `arm_cells.R` | the **arm cell source** — one reader for all three arms: region tier, whole-slide tier, metrics, provenance |
| `membership.R` | **where the cells come from and which are inside a tumour annotation** — `membership_data(mode)`, the one knob each clinical page turns |
| `aggregation_compare.R` | the annotation-aggregation sensitivity grid |
| `plot_theme.R` | the house figure style (see below) |
| `pdf_export.R` | `export_pdf_figures(slug)` — collect a page's PDFs into `output/figures/<slug>/` |
| `benchmark_plots.R` | the benchmark sweep figures (vendored fork of mirage's `plots.R`) |
| `registration_accuracy_plots.R` | the **sweep** registration-accuracy figures — the only place they are built |
| `registration_arms.R` | the **arm sweep** on the real slides — one run per configuration (both registration backends), ranked, with the cross-arm comparability guard |
| `run_qc.R` | the **run's own** QC on the study samples: readers for mirage's per-patient QC artifacts, plus `build_run_qc_figs()` |
| `paper_figures.R` | the **manuscript panels** — re-cuts of existing quantities in the shape each figure legend asks for |

The dependency order is `cell_tables.R` + `arms.R` → `validation_helpers.R` → `membership.R`
(→ `mirage_cells.R`, `arm_cells.R`); sourcing `validation_helpers.R` pulls in the first two and
`plot_theme.R`, so an analysis that sources it is loaded and styled with nothing
further to call. `cell_tables.R` and `arms.R` are the bottom of the stack — base R + `tibble` / `fs`
respectively — so both can be sourced and tested on their own, with no data on disk.
`arms.R` in particular is pure path logic, which is why `tests/testthat/test-arms.R`
runs anywhere and is the first thing to fail if a producer renames a directory. `sf` is a **lazy** dependency
of `validation_helpers.R`: only the six geojson functions require it, so the
flag-membership reports source and run on a machine without it.

## Reading a cell table

Three upstream tools export "one row per cell" and spell the same four facts three
different ways: FlowPath's `PhenotypeCsvExporter`, mirage's `bin/phenotype_cells.py`,
and mirage's `bin/join_flowpath.py --out-table`. `cell_tables.R` is where that is
reconciled, and **no analysis should name an export's column directly**:

| instead of | use |
|---|---|
| `str_extract(phenotype, "(?<=\\().*?(?=\\))")` | `clean_phenotype()` / `cell_phenotype(cells)` |
| `cells$Out_of_annotation` | `cell_outside(cells)`, `has_outside_flag(cells)` |
| `cells$CD3_sign` | `marker_pos(cells, "CD3")` |
| `cells$centroid_x / 0.325` | `cell_centroids_px(cells, um_per_px)` |
| `"cell_id"` for a dedup | `cell_key_cols(cells)` |

Each accessor probes the columns it finds, so a marker the export never gated
contributes no cells rather than aborting the report — and a mirage table whose
coordinates are already in pixels is not rescaled a second time. `read_cell_csv()`
applies `normalise_cell_flags()`, which forces the boolean columns to real logicals so
two exports of the same cohort can be bound without a type clash. The file header
documents all three schemas.

## Three phenotyping arms

Three tools phenotype **the same slides** three different ways, so the diff between
any two of them isolates a *method* rather than a cohort. `arms.R` is the one place
that knows which files belong to which arm:

```
data/massimo1/                                   arm 1  (+ arm 3, which shares its root)
  FlowPath_csv_selected/<pid>_a<k>.csv                   region cells   (FLAT)
  annotation_selected/<pid>/<pid>_a<k>.geojson           region polygons
  FlowPath_csv_all/<pid>/<pid>.csv                       whole-slide cells
  annotation_all/<pid>/<pid>.geojson                     union polygon
  csv_inverted-classification_modified-thrPANCK/         arm 3 region cells (FLAT)
    <pid>_a<k>.csv
data/massimo2/                                   arm 2
  csv/<pid>/<pid>_<A|B|C>.csv                            region cells
  csv/<pid>/<pid>.csv                                    whole slide, no regions
  annotation/<pid>/<pid>_<A|B|C>.geojson                 region polygons
```

Symlink the cluster trees in — one link per root, and arm 3 comes with arm 1:

```sh
ln -s <share>/Massimo1 data/massimo1
ln -s <share>/Massimo2 data/massimo2
```

Four things about this differ from every other source, and all four are decisions
rather than accidents:

- **The arms' regions are drawn INDEPENDENTLY.** `massimo1`'s `ANNOTATION_2` is not
  `massimo2`'s `ANNOTATION_2` — they are separate annotation sessions over the same
  tissue. So `ANNOTATION_<k>` is an **arm-local** label, each arm carries its own
  pathologist table (`neoplastic_massimo1` / `neoplastic_massimo2`), and any cross-arm
  figure joins on **patient**, never on `(patient, annotation)`. 24086 is the case that
  makes it concrete: three regions in arm 1, none at all in arm 2, where it is a bare
  whole-slide export. Both are correct.
- **Region numbering follows the arm's own suffix.** Arm 2 suffixes with a LETTER read
  by **alphabet position** — `C` is region 3 whether or not `B` was exported. Arm 1
  suffixes with `_a<k>` and `k` is taken literally. `arm_parse_name()` is the only
  parser; `.annotation_key()` in `validation_helpers.R` delegates to it rather than
  keeping a second copy.
- **The two tiers ARE the two scopes.** An arm shipping an `_all` tier answers `union`
  from it directly — a real whole-slide export against a real union polygon, no
  inference. An arm without one dissolves its region polygons and de-duplicates its
  pooled region files. A patient with a whole-slide export but **no** region files
  (arm 1's 10338 and 15897) has its union polygon promoted to `ANNOTATION_1`, so it
  keeps a per-region row instead of vanishing from every per-annotation panel.
- **No annotation directory means everything is inside — per arm.** That is arm 2's
  own stated convention (`bare_region_is = "whole_slide"` in the registry), and it is
  spelled out per-arm rather than applied globally: everywhere else a missing polygon
  is a reason to *drop* a patient, and reinterpreting it as "all in" everywhere would
  turn a data problem into a silently 100 %-inside patient. It is recorded in the
  metrics frame's `source` column as `whole_slide`.

**The export shape is REPORTED, not assumed.** Whether a patient's region files repeat
the same cells (each holding the whole slide with `Out_of_annotation` computed for that
region) or partition it is a property of the producer, and it sets every cohort-level
denominator while leaving every number plausible either way. The code is correct under
both — per-region metrics intersect each file with its own polygon, and
`arm_cohort_cells()` de-duplicates on `cell_key_cols()` — and `arm_overlap_report()`
says which shape the data on disk actually is on every knit.

**Arm 1 is the ground truth for that de-duplication.** It is the only arm publishing
both a whole-slide export *and* region files, so it is the only place the procedure
arms 2 and 3 are *forced* to use can be checked rather than trusted:
`arm_wholeslide_reconciliation()` runs the dedup on arm 1's region files and compares
it against arm 1's own export. The deduped union is expected to be a **subset** (the
selected regions need not cover the slide); deduped **>** whole-slide is the one
direction coverage cannot explain, and it means `cell_key_cols()` has stopped
identifying cells — inflating every arm without an `_all` tier by the same factor.

Membership inside a region prefers, in order: the region's own geojson (`sf` — the
only source that knows the region's **area**, hence the only one that yields
densities), then the export's `Out_of_annotation` flag (`flag`), then every row
(`whole_slide`). Whichever was used is in the `source` column, never inferred.

Reach any arm through `membership_data("<arm>")`. Unlike `"mirage"` the arms read their
own polygons, so passing `annots` is an **error** rather than being ignored — an
outside set would score one arm's cells against another arm's regions.

**`"all_slide"` was renamed to `"massimo2"` on 2026-08-26**, when the second and third
arms arrived and one export stopped being *the* export. Same tree, same rules. The old
name now fails `match.arg()` loudly rather than partial-matching onto something
plausible. **Three modes were removed on 2026-08-11** (`geojson`, `flag`, `flag_old`)
because their layouts are no longer produced.

**There is no bare `ihc_data` or `neoplastic_data` any more.** Every page declares
`ARM <- "..."` next to its `SLUG` and takes `ihc_for(ARM)` / `neoplastic_for(ARM)`. A
default would let a page use one arm without saying so, and the resulting figure would
be indistinguishable from a figure about a different arm.

## VALIS's own error: three stages, two columns, two files

**There is no micro column.** VALIS's `error_df` schema is `from`/`filename`,
`rigid_D`, `non_rigid_D` — that is the whole of it. Micro-registration is not a fourth
stage with a column of its own; it is an **update to the non-rigid field**
(`register_micro()` re-runs `measure_error()` and overwrites `<name>_summary.csv`,
composing the micro residual into that same field). So the micro number does not live
in a column: **it lives in the difference between the two files.**

`*_summary_premicro.csv` is a copy taken just before that overwrite, and it is written
**only at `reg_micro_reg = 2`**. The stage axis is therefore recovered from *which file*
a value came from, not from which column — the same reconstruction mirage's own report
performs in `bin/generate_qc_report.py:_RECONCILE_TRE_SOURCE`:

| stage | column | file |
|---|---|---|
| `rigid` | `rigid_D` | final — unchanged by micro |
| `non_rigid` | `non_rigid_D` | **pre-micro** — the only file that isolates it |
| `micro` | `non_rigid_D` | **final** — same column, other file |

`valis_error_long()` does exactly this. Faceting by source file — what it used to do —
drew two panels in which `rigid_D` was the same number twice and the one differing point
carried the same label in both, so they looked identical because they mostly were.

**A blank `micro` stage is a claim, not a gap.** With no premicro sibling the final
`non_rigid_D` *is* the pre-micro value, so it becomes `non_rigid` and no `micro` level is
emitted. Emitting a byte-for-byte duplicate would read as "micro bought nothing", which
is not the same statement as "micro did not run". mirage makes the same choice on the
segmentation-overlap side (`docs/registration_qc.md`).

`original` stays in `VALIS_STAGE_LEVELS` but is normally **absent**: VALIS's per-patient
summary has no such column. The benchmark sweep's `registration_valis_rtre.csv` does — a
different artifact, from `make_tables.py` — so the level survives and the column is
detected rather than named, and a build that emits one plots instead of being dropped.

STARE never reaches this axis: the tiled backend writes no `registered/summary/*.csv`, so
it has no micro stage by construction. Its intrinsic TRE is read from `*_tre.json` in
**pixels** and plotted separately. mirage's report folds both backends into one slide dict
by reusing the `rigid_D` / `non_rigid_D` keys; this project keeps them apart because the
units differ and a shared axis invites a comparison that is not one.

`valis_error_long()` also **passes a caller's extra grouping columns through**
(`arm`, `backend`, `memory_mode`, `micro_reg`). It returning a narrower frame than it was
handed is how the arm sweep's figure lost its grouping and quietly stopped being built.
## Two registration backends, two stage vocabularies

`lib/WarpBackends.groovy` in mirage declares one `reg_qc=2` stage list per
registration method, and they are **not** the same vocabulary:

| `registration_method` | stages |
|---|---|
| `valis` (default) | `native → rigid → non_rigid → micro` |
| `tiled` (STARE) | `native → rigid → refined` |

`QC_STAGE_LEVELS` in `run_qc.R` is the union, and it has to be: `read_seg_qc()`
factors `stage` against it and filters out non-matches, so a VALIS-only list made
every tiled run's `refined` rows **disappear with no error and no warning** — and the
arm's final stage then came out as `rigid`, which reads as the backend performing far
worse than it does.

Two rules follow, both enforced in `registration_arms.R` rather than documented:

- **"The last stage" is read from each run's own `stage_order`** (`stage_index`), never
  from `QC_STAGE_LEVELS`' ordering. The two vocabularies interleave under any single
  ordering, so a shared ordering compares one backend's stage position against the
  other's.
- **`arm_comparable_stages()` returns only `native` for a mixed-backend sweep.** Only
  `native` is shared as both a spelling and a meaning; `rigid` is a shared *word* and
  not a shared *operation* (VALIS: affine, composed with micro-rigid at
  `reg_micro_reg ≥ 1`; STARE: the coarse global anchor before mesh refinement).

What *is* comparable across backends is the metric itself: mirage's segmentation-overlap
scorer takes `--method` and builds its warper from either a VALIS registrar pickle or a
STARE transform manifest, so it is the same measurement either way. Each backend's
**own** reported error is not: VALIS writes rTRE (a fraction of the image diagonal) to
`registered/summary/*.csv`, STARE writes TRE in **pixels** to
`qc/registration/*_tre.json`. Separate readers, separate figures, never one axis.

## Child documents

**There are two, and they are what keeps the three arms from drifting apart.**
`analysis/_children/clinical_body.Rmd` and `molecular_body.Rmd` hold the shared
analysis; six thin parents (`clinical_massimo1/2/1_inverted`, `molecular_...`)
supply only the YAML title, `SLUG`, `ARM`, the arm's objects and the
`export_pdf_figures(SLUG)` call. A parent is ~70 lines; the bodies are ~900 and
~1200. Six full copies would drift the first time one was edited, and the drift
would be invisible because every copy still knits.

Three rules, the first two enforced by `tests/testthat/test-child-documents.R`:

- Children live in **`analysis/_children/`**, not beside their parents. An underscore
  on the *file* only hides it from `render_site()`, which skips `^[_.]` resources;
  it does **not** hide it from `Sys.glob()`, which is what `wflow_build("analysis/*.Rmd")`
  expands — glob's `*` refuses a leading dot, not a leading underscore. That glob would
  hand the body to the builder as a page and the build dies on the child's contract
  check. The underscore on the *directory* hides it from `render_site()`, and a
  subdirectory is out of reach of a non-recursive glob, so both routes are closed.
- The `child=` path must be **absolute**, via `here::here("analysis", "_children", ...)`.
  knitr resolves it against the knit working directory, and `_workflowr.yml` sets
  `knit_root_dir: "."` — the project root — so a bare filename is looked for beside
  `_workflowr.yml` rather than beside its parent. The failure is
  `Error in file(con, "r") : cannot open the connection` partway through the parent,
  naming neither the child nor the path, because that is `readLines()`'s internal
  call inside `knitr:::call_block`.
- A parent must define everything the body reads before including it: `SLUG`, `ARM`,
  `ihc_data`, `neoplastic`, and the sourced helpers. `tests/testthat/test-chunk-gates.R`
  splices children into their parents exactly as knitr does, so a gate defined in a
  parent's setup and used in a child's chunk is checked as one page.

**The `child =` chunk option is BANNED; bodies are spliced with
`knit_child(text = readLines(here::here(...)))`.** workflowr's `fig.path` option hook
does `options$fig.path <- create_figure_path(knitr::current_input())`, and inside a
`child =` chunk `current_input()` returns the **child** — so all three clinical
parents would write their plots to `figure/clinical_body/`, one shared directory, and
whichever page knitted last would silently overwrite the other two.
`export_pdf_figures(SLUG)` then copies those same files into all three
`output/figures/<slug>/`, so a panel labelled massimo1 would show massimo2's cells
with nothing on the page to say so. `knit_child(text = ...)` has no input file of its
own, so `current_input()` stays the PARENT: figures land in `figure/<parent>/`, and
workflowr's *"custom fig.path was ignored"* warning stops firing as a side effect —
that warning was the visible symptom of the collision, never the problem itself.

`tests/testthat/test-child-documents.R` pins the idiom and fails on any `child =`.

## Sweep QC vs run QC

Two different questions, two different trees, and they must not be conflated:

| | reads | measures |
|---|---|---|
| `benchmark_plots.R`, `registration_accuracy_plots.R` | `data/benchmark/` — mirage's `benchmarks/` sweep | cost, scaling and accuracy on **synthetic** images with a known injected offset |
| `registration_arms.R` | `data/registration_arms/<arm>/` — one mirage run per configuration | which configuration to ship, measured on the **real** slides |
| `run_qc.R` | `data/mirage/<patient>/` — an ordinary pipeline run | how well **the study samples** were registered |

`data/mirage` must **be** the mirage `--outdir`, i.e. contain the per-patient
directories directly — symlink the run rather than copying a subtree.

`run_qc.R` deliberately does not source `load_data.R` or `mirage_cells.R`: the QC
page reads the run's QC artifacts, not its cells, so it renders from a pipeline run
alone with no `counts.RData`, no clinical table and no FlowPath gating.

**Search the artifact directories recursively.** Nextflow's `publishDir` preserves a
process's producer subdirectory, so `REGISTER`'s VALIS summaries — declared as
`path("preprocessed/data/*.csv")` — land at
`<patient>/registered/summary/preprocessed/data/*.csv`, not in `registered/summary/`
itself. A one-level listing finds them on a hand-flattened tree and never on a real
run. It also breaks *patient detection*, which keys on finding files: that directory
holds only a directory, so a VALIS run with `reg_qc < 2` used to
render "No mirage QC found" with its rTRE sitting right there. Join moving slides on
`slide_token`, not on `moving`: the tiled path names the artifact
`<patient_id>_<channels>_tre.json` while `seg_qc` carries the native image stem.

Its one rule: **VALIS rTRE and STARE TRE are intrinsic** — each method scoring itself
on the features it registered on — so neither is evidence on its own, and they are
not comparable to each other (only one path runs). The matched-nucleus Dice from
`*_seg_qc.json` is computed from DAPI segmentation instead, which is what makes it
the independent check; §4 of the page plots the two against each other.

## The two phenotype vocabularies

FlowPath and mirage name the **same taxonomy** differently — FlowPath's gate tree
writes `PANCK+Tumor`, `T helper`, `CD8+ T reg`; mirage's `panel.yaml` writes
`PANCK_Tumor`, `T_helper`, `CD8_Treg`. An unmapped label does not error: it joins to
`lineage = NA` and quietly empties the composition panels. So `phenotype_lineage`
(in `cell_tables.R`) lists both spellings and is joined on `pheno_join_key()`, which
strips punctuation and aliases the two names that genuinely differ (`NK_cell`,
`Activated_NK`). Use `cell_lineage(phenotype_clean)`, never a direct join.

Two differences are real rather than cosmetic:

- **`Myeloid` / `Macrophage_M2`** are mirage-only leaves; FlowPath's tree dead-ends
  at plain `Immune` on that branch, so both fold into `Immune_other`.
- **mirage never writes `Unknown`.** It emits four reserved outcomes —
  `Unclassified`, `Ambiguous`, `Conflict`, `Artefact` — for cells its constraint
  solver could not commit. `is_unresolved_phenotype()` covers all of them plus
  FlowPath's `Unknown`; matching only `"Unknown"` would keep mirage's unresolved
  cells in the QC-filtered denominators and drop FlowPath's, which is exactly the
  asymmetry that makes the two tools look different when they are not.

`tests/testthat/test-cell-tables.R` asserts every phenotype `panel.yaml` can emit is
mapped — if mirage adds a leaf, that test fails rather than the panels going quiet.

## Figure style

`plot_theme.R` is the single house style for every figure on the site. It is
sourced by `validation_helpers.R` and by `benchmark_plots.R`, so any analysis
that loads either one is styled automatically — there is nothing to call.

It exports:

| | |
|---|---|
| `theme_paper(base_size, grid, axis_lines)` | the theme; applied via `theme_set()` on source |
| `theme_paper_tile()` | heatmap variant (no axis line, no ticks) |
| `theme_paper_panels()` | dense facet grid: hairline panel border + wider gutter, for `scales = "free"` grids where adjacent panels do NOT share an axis |
| `oi`, `oi_ext` | Okabe-Ito categorical palette (8) and its 16-colour extension |
| `scale_*_oi()`, `scale_*_ordinal()` | categorical / ordered-discrete scales |
| `scale_*_div()`, `scale_*_seq()` | diverging (blue-white-red) and sequential ramps |
| `guide_cbar()` | short horizontal colourbar |
| `hotcold_cols()`, `hotcold_order()` | immune-phenotype colours and level order |
| `REF_LINE`, `FIT_LINE` | colours for reference lines and fitted trends |
| `pt_line()`, `pt_text()` | pt -> ggplot2 `linewidth` / geom text `size` |

**The one rule:** never write `theme_classic()` / `theme_bw()` / `theme_minimal()`
in an analysis. Those *replace* the active theme, which is how the figures drifted
apart in the first place. A bare `theme(...)` layers on top and is fine — use it
for genuine per-plot deviations (rotated tick labels, horizontal facet strips).

Sourcing `plot_theme.R` also sets ggplot2's `discrete.colour`/`discrete.fill`
options, so a scale you *don't* specify falls back to the house palette instead of
ggplot2's hue rainbow. New plots are publication-ready by default.

**Do not name a font family.** `theme_paper()` uses `""`, which means "the device's
own default" — Helvetica on `pdf()`, a Helvetica-metric sans on cairo/ragg. A named
family has to be registered in the *device's* database (`pdfFonts()`), not merely
installed on the machine, so an OS face like `Arial` or `DejaVu Sans` aborts every
knit that renders PDFs with `invalid font type` — and the PNG pass gives no warning,
because cairo resolves OS fonts happily. If a journal demands a specific face,
register it with `grDevices::pdfFonts()` and set `options(ihc.plot.family = "...")`
before sourcing; `paper_family()` validates it and degrades with a warning.

For a final multi-panel figure, drop the per-panel `title`/`subtitle`
(`labs(title = NULL, subtitle = NULL)`) and let `patchwork::plot_annotation(tag_levels = "a")`
place the panel letters — a tag and a title share the same slot and will overlap.
