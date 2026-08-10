# Code

Shared R sourced by the analyses in `analysis/`, plus standalone scripts.

## What is where

| file | owns |
|---|---|
| `load_data.R` | the raw loaders (`dds`, `clinical_data`, `neoplastic_data`, `counts_data`, `ihc_data`) |
| `cell_tables.R` | the **single-cell export schema** — one vocabulary over three upstream formats |
| `validation_helpers.R` | the derived quantities (region ratios, composition, marker/lineage fractions, invasive margin, agreement stats) |
| `membership.R` | **which cells are inside a tumour annotation** — `membership_data(mode)`, the one knob the three `clinical_data` pages turn |
| `aggregation_compare.R` | the annotation-aggregation sensitivity grid |
| `plot_theme.R` | the house figure style (see below) |
| `pdf_export.R` | `export_pdf_figures(slug)` — collect a page's PDFs into `output/figures/<slug>/` |
| `benchmark_plots.R` | the benchmark sweep figures (vendored fork of mirage's `plots.R`) |
| `registration_accuracy_plots.R` | the registration-accuracy figures — **the only place** they are built |

The dependency order is `cell_tables.R` → `validation_helpers.R` → `membership.R`;
sourcing `validation_helpers.R` pulls in the first and `plot_theme.R`, so an analysis
that sources it is loaded and styled with nothing further to call.

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
