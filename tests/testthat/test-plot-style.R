# The house style, as assertions.
#
# plot_theme.R can only keep the figures looking like one figure set if the call
# sites actually go through it, and nothing about ggplot forces them to: a plot
# with a hand-rolled palette or a hand-picked fig.width renders perfectly well and
# simply looks different from its neighbours. Those are the drifts this file
# catches, because they are invisible in a code review and obvious in a printed
# manuscript.
#
# The three properties under test:
#   1. SIZE      one shrink factor for every figure, so type is one size on paper
#   2. COLOUR    one named palette, so a population is one colour across the site
#   3. n         every categorical axis states how many observations back each level
source(here::here("code", "plot_theme.R"))

RMDS <- list.files(here::here("analysis"), pattern = "[.]Rmd$", full.names = TRUE)

# ---- 1. size ---------------------------------------------------------------

test_that("every column width lands the base type at the same size on paper", {
  # This is the whole point of fig_width(): NOT that widths are equal, but that
  # rendered_width * FIG_SCALE == the target column, for every column.
  pt <- vapply(c("single", "oneandhalf", "double"), function(col) {
    10 * FIG_COLUMN_MM[[col]] / (fig_width(col) * 25.4)
  }, numeric(1))
  expect_equal(unname(diff(range(pt))), 0, tolerance = 0.05)
  # ... and inside the 5-8pt range journals ask for.
  expect_true(all(pt > 5 & pt < 8.5))
})

test_that("panel_width converts placed millimetres at the same scale", {
  expect_equal(panel_width(180), fig_width("double"))
  expect_equal(panel_width(85),  fig_width("single"))
  # The mm/20 shorthand the chunk headers are read by.
  expect_equal(panel_width(140), 7.0, tolerance = 0.05)
})

test_that("fig_height never exceeds the printed-page ceiling", {
  expect_lte(fig_height("double", aspect = 10) * 25.4 * FIG_SCALE, 230)
})

test_that("no analysis chunk hard-codes a figure width", {
  # A literal fig.width is how the set drifted to eight different widths (3.4in to
  # 11in) and therefore eight different printed type sizes. Widths come from
  # fig_width() or panel_width() so the shrink factor is stated, not implied.
  offenders <- Filter(function(f) any(grepl("fig[.]width\\s*=\\s*[0-9]", readLines(f, warn = FALSE))),
                      RMDS)
  expect_equal(basename(offenders), character(0))
})

test_that("every analysis sources the house style, unconditionally", {
  # Sourcing it inside `if (have_data)` means a page whose data is absent silently
  # renders on ggplot2's default grey theme.
  missing <- Filter(function(f) {
    txt <- readLines(f, warn = FALSE)
    plots <- any(grepl("fig[.]width|ggplot|build_.*figs", txt))
    plots && !any(grepl('source\\(here\\("code", "plot_theme[.]R"\\)\\)', txt))
  }, RMDS)
  expect_equal(basename(missing), character(0))
})

# ---- 2. colour -------------------------------------------------------------

test_that("the recurring palettes are named, not positional", {
  # A named palette maps level -> colour; an unnamed one maps POSITION -> colour and
  # silently recolours a population whenever the level set changes.
  for (pal in list(LINEAGE_COLS, LINEAGE_PALETTE, ARM_KIND_COLS, GATED_COLS)) {
    expect_false(is.null(names(pal)))
    expect_true(all(nzchar(names(pal))))
  }
})

test_that("the lineage palette covers every population any figure can draw", {
  expect_true(all(c(LEGIBLE_LINEAGES, "other") %in% names(LINEAGE_PALETTE)))
  # The marker-gated readouts share the lineage axis on the composition panels.
  expect_true(all(names(GATED_COLS) %in% names(LINEAGE_PALETTE)))
})

test_that("a population keeps its colour when other populations are absent", {
  # The regression this guards: clinical_flowpath.Rmd drew CD8T blue while the Fig
  # 5(a) map drew it vermillion, because the former used positional `oi`.
  hex <- function(levels) {
    p <- ggplot(data.frame(g = factor(levels, levels = levels), y = seq_along(levels)),
                aes(g, y, fill = g)) + geom_col() + scale_fill_lineage()
    d <- ggplot2::ggplot_build(p)$data[[1]]
    stats::setNames(d$fill, levels)
  }
  full  <- hex(c("Tumor", "CD8T", "CD4T", "Treg", "NK"))
  subset <- hex(c("CD8T", "NK"))
  expect_equal(unname(full[["CD8T"]]), unname(subset[["CD8T"]]))
  expect_equal(unname(full[["NK"]]),   unname(subset[["NK"]]))
  expect_equal(unname(full[["CD8T"]]), unname(LINEAGE_COLS[["CD8T"]]))
})

# Comment lines are stripped before any of the source scans below. Without that a
# file can satisfy the check with a comment SAYING to use the shared scale while the
# code beside it does the opposite — which is exactly what happened when this test
# was first written: the mutation "swap scale_colour_lineage() for a positional
# scale_colour_manual(values = oi)" passed, because the doc block above the function
# still mentioned the right name.
.code_only <- function(path) {
  txt <- readLines(path, warn = FALSE)
  paste(txt[!grepl("^\\s*#", txt)], collapse = "\n")
}

test_that("files that colour by lineage go through the shared scale", {
  srcs <- c(RMDS, list.files(here::here("code"), pattern = "[.]R$", full.names = TRUE))
  srcs <- setdiff(srcs, here::here("code", "plot_theme.R"))
  offenders <- Filter(function(f) {
    txt <- .code_only(f)
    maps_lineage <- grepl("(fill|colour|color)\\s*=\\s*lineage", txt) ||
                    grepl("aes\\(lineage", txt)
    maps_lineage && !grepl("scale_(fill|colour|color)_lineage", txt)
  }, srcs)
  expect_equal(basename(offenders), character(0))
})

test_that("legends show display names, not column codes", {
  expect_equal(unname(.lineage_labels("Immune_other")), "Immune (other)")
  expect_equal(unname(.lineage_labels("Tumor")), "Tumour")
  # An unmapped level passes through rather than becoming NA and vanishing.
  expect_equal(unname(.lineage_labels("CD45+")), "CD45+")
})

# ---- 3. n ------------------------------------------------------------------

test_that("label_n counts what was plotted and survives an absent level", {
  f <- label_n(c("a", "a", "b"))
  expect_equal(f(c("a", "b")), c("a\n(n = 2)", "b\n(n = 1)"))
  # A level present as a break but absent from the data must not read "(n = NA)".
  expect_equal(f("c"), "c")
  expect_equal(label_n(c("a", "a"), sep = " ")("a"), "a (n = 2)")
})

test_that("label_n ignores NA, which would otherwise be counted as a group", {
  expect_equal(label_n(c("a", NA, "a"))("a"), "a\n(n = 2)")
})

test_that("n_note names its unit, because n is ambiguous in this project", {
  # 24 patients, 24 slides and 24 runs are all plausible and mean different things.
  expect_equal(n_note(24, "patients"), "n = 24 patients")
  expect_equal(n_note(c("a", "b", "a"), "slides"), "n = 2 slides")
  expect_equal(with_n("subtitle", 3, "runs"), "subtitle · n = 3 runs")
  expect_equal(with_n(NULL, 3, "runs"), "n = 3 runs")
})
