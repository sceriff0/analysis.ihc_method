# Code

Shared R sourced by the analyses in `analysis/`, plus standalone scripts.

## What is where

| file | owns |
|---|---|
| `load_data.R` | the raw loaders (`dds`, `clinical_data`, `neoplastic_data`, `counts_data`, `ihc_data`) |
| `cell_tables.R` | the **single-cell export schema** — one vocabulary over three upstream formats |
| `validation_helpers.R` | the derived quantities (region ratios, composition, marker/lineage fractions, invasive margin, agreement stats) |
| `mirage_cells.R` | the mirage cell source — `phenotypes.csv` + `morphology.csv` joined per patient |
| `membership.R` | **where the cells come from and which are inside a tumour annotation** — `membership_data(mode)`, the one knob the four `clinical_data` pages turn |
| `aggregation_compare.R` | the annotation-aggregation sensitivity grid |
| `plot_theme.R` | the house figure style (see below) |
| `pdf_export.R` | `export_pdf_figures(slug)` — collect a page's PDFs into `output/figures/<slug>/` |
| `benchmark_plots.R` | the benchmark sweep figures (vendored fork of mirage's `plots.R`) |
| `registration_accuracy_plots.R` | the **sweep** registration-accuracy figures — the only place they are built |
| `run_qc.R` | the **run's own** QC on the study samples: readers for mirage's per-patient QC artifacts, plus `build_run_qc_figs()` |

The dependency order is `cell_tables.R` → `validation_helpers.R` → `membership.R`
(→ `mirage_cells.R`); sourcing `validation_helpers.R` pulls in the first and
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

## Sweep QC vs run QC

Two different questions, two different trees, and they must not be conflated:

| | reads | measures |
|---|---|---|
| `benchmark_plots.R`, `registration_accuracy_plots.R` | `data/benchmark/` — mirage's `benchmarks/` sweep | cost, scaling and accuracy on **synthetic** images with a known injected offset |
| `run_qc.R` | `data/mirage/<patient>/` — an ordinary pipeline run | how well **the study samples** were registered and phenotyped |

`run_qc.R` deliberately does not source `load_data.R` or `mirage_cells.R`: the QC
page reads the run's QC artifacts, not its cells, so it renders from a pipeline run
alone with no `counts.RData`, no clinical table and no FlowPath gating.

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
