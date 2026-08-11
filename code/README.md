# Code

Shared R sourced by the analyses in `analysis/`, plus standalone scripts.

## What is where

| file | owns |
|---|---|
| `load_data.R` | the raw loaders (`dds`, `clinical_data`, `neoplastic_data`, `counts_data`, `ihc_data`) |
| `cell_tables.R` | the **single-cell export schema** — one vocabulary over three upstream formats |
| `validation_helpers.R` | the derived quantities (region ratios, composition, marker/lineage fractions, invasive margin, agreement stats) |
| `mirage_cells.R` | the mirage cell source — `phenotypes.csv` + `morphology.csv` joined per patient |
| `all_slide.R` | the **all-slide cell source** — one CSV *and* one geojson per annotation region, nested per patient |
| `membership.R` | **where the cells come from and which are inside a tumour annotation** — `membership_data(mode)`, the one knob the four `clinical_data` pages turn |
| `aggregation_compare.R` | the annotation-aggregation sensitivity grid |
| `plot_theme.R` | the house figure style (see below) |
| `pdf_export.R` | `export_pdf_figures(slug)` — collect a page's PDFs into `output/figures/<slug>/` |
| `benchmark_plots.R` | the benchmark sweep figures (vendored fork of mirage's `plots.R`) |
| `registration_accuracy_plots.R` | the **sweep** registration-accuracy figures — the only place they are built |
| `registration_arms.R` | the **arm sweep** on the real slides — one run per configuration (both registration backends), ranked, with the cross-arm comparability guard |
| `run_qc.R` | the **run's own** QC on the study samples: readers for mirage's per-patient QC artifacts, plus `build_run_qc_figs()` |
| `paper_figures.R` | the **manuscript panels** — re-cuts of existing quantities in the shape each figure legend asks for |

The dependency order is `cell_tables.R` → `validation_helpers.R` → `membership.R`
(→ `mirage_cells.R`, `all_slide.R`); sourcing `validation_helpers.R` pulls in the first and
`plot_theme.R`, so an analysis that sources it is loaded and styled with nothing
further to call. `cell_tables.R` is the bottom of the stack and stays base-R +
`tibble`, so it can be sourced and tested on its own. `sf` is a **lazy** dependency
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

## The all-slide layout

`all_slide.R` reads the fifth cell source, and the only one where the
pathologist's polygon and the cell export are matched **per region**:

```
data/all_slide/csv/<patient>/<patient>_<A|B|C>.csv
data/all_slide/csv/<patient>/<patient>.csv              (no regions)
data/all_slide/annotation/<patient>/<patient>_<A|B|C>.geojson
```

On the cluster this tree is
`/hpcnfs/techunits/imaging/work/ATTEND/Mirage/all-slide_new`; symlink it in with
`ln -s <that path> data/all_slide`.

Three things about it differ from every other source, and all three are decisions
rather than accidents:

- **Region letters are positions, not names.** `A` → `ANNOTATION_1`, `B` →
  `ANNOTATION_2`, and so on by alphabet index, so the regions line up with
  `neoplastic_data`'s `ANNOTATION_1..3` columns. `C` stays region 3 whether or not
  `B` was exported.
- **No annotation directory means everything is inside.** That is this layout's
  own stated convention, and it is the reason `all_slide.R` exists instead of the
  older loaders being widened: everywhere else a missing polygon is a reason to
  *drop* a patient, and reinterpreting it as "all in" globally would turn a data
  problem into a silently 100 %-inside patient. Here it is deliberate and is
  recorded in the metrics frame's `source` column as `whole_slide`.
- **The union de-duplicates, and the export shape is REPORTED rather than assumed.**
  Whether a patient's region files repeat the same cells (each holding the whole slide
  with `Out_of_annotation` computed for that region) or partition it (each holding only
  its own region) is a property of the producer, not of the layout — and it sets every
  cohort-level denominator while leaving every number plausible either way. The code is
  correct under both: the per-region metrics intersect each file with its own polygon,
  and `all_slide_union_cells()` keys on `cell_key_cols()` (or the rounded centroid)
  before pooling. `all_slide_overlap_report()` says which shape the data on disk is, and
  `clinical_data.Rmd` prints it on every knit.

Membership inside a region prefers, in order: the region's own geojson (`sf` —
the only source that knows the region's **area**, hence the only one that yields
densities), then the export's `Out_of_annotation` flag (`flag`), then every row
(`whole_slide`). Whichever was used is in the `source` column, never inferred.

Reach it through `membership_data("all_slide")` like any other mode — it reads its
own polygons, so unlike `"mirage"` it takes no `annots` argument.

**Three membership modes were removed on 2026-08-11** (`geojson`, `flag`, `flag_old`)
because their layouts are no longer produced: `geojson` wanted a whole-slide
`data/flowpath/<patient>.csv` plus a FLAT `data/annotation/<patient>_a<k>.geojson`, and
the flag modes wanted `data/flowpath/per_annotation/` with an `old/` overlay. Two pages
went with them (`clinical_data_per_annotation`, `..._per_annotation_old`). What survives
is `all_slide` + `mirage`, whose diff isolates the phenotyping method.

`load_annotations()` reads **both** the nested tree and the flat legacy one, keyed by
`.annotation_key()` — the same parser `annotation_membership_qc()` uses on the cell
csvs, which is what guarantees a csv and the polygon it is compared against agree on
which region they are.

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

The four `clinical_data` pages are thin parents over `analysis/_clinical_data_body.Rmd`.
Two rules, both enforced by `tests/testthat/test-child-documents.R`:

- The **leading underscore** is what keeps `render_site()` (and so workflowr) from
  building the body as a page of its own.
- The `child=` path must be **absolute**, via `here::here("analysis", ...)`. knitr
  resolves it against the knit working directory, and `_workflowr.yml` sets
  `knit_root_dir: "."` — the project root — so a bare filename is looked for beside
  `_workflowr.yml` rather than beside its parent. The failure is
  `Error in file(con, "r") : cannot open the connection` partway through the parent,
  naming neither the child nor the path, because that is `readLines()`'s internal
  call inside `knitr:::call_block`.

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
