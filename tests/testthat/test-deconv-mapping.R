# The deconvolution <-> IHC mapping. Each test pins one decision taken against the
# phenotyping gate tree (CD45 -> CD3 -> CD8/CD4 -> GZMB/FOXP3; CD3- -> CD56 -> GZMB;
# CD45- -> SMA -> PANCK -> VIMENTIN), because a wrong mapping still knits: it just
# correlates two different populations and reports a rho.
source(here::here("code", "cell_tables.R"))
source(here::here("code", "validation_helpers.R"))

.ct <- function(method, cell_type) tibble::tibble(method = method, cell_type = cell_type)

.lineages_of <- function(contrib, cell_type) sort(contrib$lineage[contrib$cell_type == cell_type])

# --- Deconvolution side: labels -> lineages ------------------------------------
test_that("follicular helper T cells are CD4T: the tree's T helper is CD3+CD8-CD4+FOXP3-", {
  expect_equal(deconv_to_lineage("T cells follicular helper"), "CD4T")
  expect_equal(deconv_to_lineage("T cell follicular helper"),  "CD4T")
})

test_that("gamma delta, MAIT and NKT have no single leaf and stay out of the subset lineages", {
  expect_true(all(is.na(deconv_to_lineage(c("T cells gamma delta", "T cell MAIT",
                                             "T cell gamma delta VD2", "T cell NK")))))
})

test_that("MCP-counter's pan-T and cytotoxic signatures get a counterpart", {
  expect_equal(deconv_to_lineage("T cell"),             "T_total")
  expect_equal(deconv_to_lineage("cytotoxicity score"), "Cytotoxic")
  expect_equal(deconv_to_lineage("Cancer associated fibroblast"), "Stroma")
})

test_that("a method with no Treg output compares its CD4 against CD4 incl. Treg", {
  # EPIC's tumour CD4 reference is 'CD4 but not CD8' — Tregs included.
  c <- deconv_contributions(.ct("epic", c("T cell CD4+", "T cell CD8+", "NK cell")))
  expect_equal(.lineages_of(c, "T cell CD4+"), c("CD4T_all", "Immune_total", "T_total"))
  expect_false("CD4T" %in% c$lineage)
})

test_that("an additive method with Tregs gets CD4T, Treg and their sum", {
  c <- deconv_contributions(.ct("quantiseq", c("T cell CD4+ (non-regulatory)",
                                               "T cell regulatory (Tregs)", "T cell CD8+",
                                               "NK cell", "uncharacterized cell")))
  expect_equal(.lineages_of(c, "T cell CD4+ (non-regulatory)"), c("CD4T", "CD4T_all", "Immune_total", "T_total"))
  expect_equal(.lineages_of(c, "T cell regulatory (Tregs)"),    c("CD4T_all", "Immune_total", "T_total", "Treg"))
  expect_equal(.lineages_of(c, "T cell CD8+"),                  c("CD8T", "Cytotoxic", "Immune_total", "T_total"))
  expect_equal(.lineages_of(c, "NK cell"),                      c("Cytotoxic", "Immune_total", "NK"))
  expect_equal(.lineages_of(c, "uncharacterized cell"),         "NonImmune")
})

test_that("CIBERSORTx's CD4 subsets and Tfh all sum into CD4T", {
  lm22 <- c("T cells CD4 naive", "T cells CD4 memory resting",
            "T cells CD4 memory activated", "T cells follicular helper",
            "T cells regulatory (Tregs)", "T cells gamma delta")
  c <- deconv_contributions(.ct("cibersortx", lm22))
  expect_setequal(c$cell_type[c$lineage == "CD4T"], lm22[1:4])
  # gamma delta joins only the pan-T sum
  expect_equal(.lineages_of(c, "T cells gamma delta"), "T_total")
})

test_that("a non-additive method uses its whole-lineage label and never sums", {
  # xCell reports both 'CD4+ T-cells' and its subsets; enrichment scores do not add.
  c <- deconv_contributions(.ct("xcell", c("T cell CD4+ (non-regulatory)", "T cell CD4+ Th1",
                                           "T cell CD4+ memory", "T cell regulatory (Tregs)")))
  expect_equal(c$cell_type[c$lineage == "CD4T"], "T cell CD4+ (non-regulatory)")
  expect_false(any(c$lineage %in% c("T_total", "CD4T_all", "Cytotoxic", "NonImmune")))
})

test_that("MCP-counter's own T-cell signature is used, not a sum of subsets", {
  c <- deconv_contributions(.ct("mcp_counter", c("T cell", "T cell CD8+", "cytotoxicity score")))
  expect_equal(c$cell_type[c$lineage == "T_total"],   "T cell")
  expect_equal(c$cell_type[c$lineage == "Cytotoxic"], "cytotoxicity score")
})

test_that("EPIC's non-immune share is uncharacterized + CAF + endothelial", {
  c <- deconv_contributions(.ct("epic", c("uncharacterized cell", "Cancer associated fibroblast",
                                          "Endothelial cell", "T cell CD8+")))
  expect_setequal(c$cell_type[c$lineage == "NonImmune"],
                  c("uncharacterized cell", "Cancer associated fibroblast", "Endothelial cell"))
  expect_equal(c$cell_type[c$lineage == "Stroma"], "Cancer associated fibroblast")
})

test_that("relative CIBERSORT is compared against CD45+ cells, everything else against all cells", {
  expect_equal(deconv_ihc_denominator(c("cibersortx", "cibersort", "quantiseq", "cibersort_abs",
                                        "not_a_method")),
               c("CD45", "CD45", "all", "all", "all"))
})

# --- IHC side: gate-tree leaves -> the compared populations ----------------------
.ihc <- function() {
  tibble::tibble(
    patient_id = "046",
    phenotype_clean = c("T cytotoxic", "Activated T cytotoxic", "CD8+ T reg",
                        "T helper", "CD4+ Treg", "Natural Killer", "Activated Natural Killer",
                        "Immune", "Immune", "PANCK+Tumor", "VIM+Tumor", "Stroma", "Unknown"),
    CD3_sign = c("+", "+", "+", "+", "+", "-", "-", "+", "-", "+", "-", "-", "-"))
}

test_that("CD8+ T regs count as CD8T, and Treg is CD4+ Treg only", {
  f <- ihc_comparison_fraction(.ihc())
  n <- stats::setNames(f$n, f$lineage)
  expect_equal(n[["CD8T"]], 3)
  expect_equal(n[["Treg"]], 1)
  expect_equal(n[["CD4T"]], 1)
  expect_equal(n[["CD4T_all"]], 2)
})

test_that("T_total adds CD3+ Immune cells (the CD3+CD8-CD4- branch) but not CD3- ones", {
  f <- ihc_comparison_fraction(.ihc())
  # 5 T leaves + the one CD3+ Immune cell; the CD3+ tumour cell is CD45- and ignored
  expect_equal(f$n[f$lineage == "T_total"], 6)
})

test_that("the CD45 denominator is the immune branch of the tree", {
  f <- ihc_comparison_fraction(.ihc())
  expect_true(all(f$n_all == 13))
  expect_true(all(f$n_cd45 == 9))
  expect_equal(f$frac_cd45[f$lineage == "NK"], 2 / 9)
  expect_equal(f$frac_all[f$lineage == "NonImmune"], 4 / 13)
})

test_that("mirage's spellings reach the same populations", {
  d <- .ihc()
  d$phenotype_clean <- c("T_cytotoxic", "Activated_T_cytotoxic", "CD8_Treg", "T_helper",
                         "CD4_Treg", "NK_cell", "Activated_NK", "Immune", "Immune",
                         "PANCK_Tumor", "VIM_Tumor", "Stroma", "Unknown")
  expect_equal(ihc_comparison_fraction(d)$n, ihc_comparison_fraction(.ihc())$n)
})

test_that("a lineage a patient lacks reads 0, not a missing row", {
  d <- .ihc()
  d <- rbind(d, transform(d[10, ], patient_id = "999"))
  f <- ihc_comparison_fraction(d)
  expect_equal(nrow(f), 2 * nrow(deconv_comparison_lineages))
  expect_equal(f$frac_all[f$patient_id == "999" & f$lineage == "CD8T"], 0)
})

test_that("pairing takes the IHC fraction on the method's denominator", {
  f <- ihc_comparison_fraction(.ihc())
  scores <- tibble::tibble(method = c("cibersortx", "quantiseq"), patient_id = "046",
                           lineage = "NK", score = c(.1, .2))
  p <- deconv_pair_with_ihc(scores, f)
  expect_equal(p$ihc_frac[p$method == "cibersortx"], 2 / 9)
  expect_equal(p$ihc_frac[p$method == "quantiseq"],  2 / 13)
  expect_equal(p$ihc_denominator, c("CD45", "all"))
})

test_that("every rule and every compared lineage is documented", {
  expect_true(all(nzchar(deconv_lineage_rules$why)))
  expect_true(all(stats::na.omit(deconv_lineage_rules$lineage) %in% deconv_comparison_lineages$lineage))
  expect_true(all(nzchar(deconv_comparison_lineages$ihc_gate)))
})

# --- Other immune cells and the macro categories -------------------------------
test_that("B, myeloid and granulocyte labels are the tree's CD3- CD56- Immune leaf", {
  lab <- c("B cell", "B cells memory", "Macrophages M2", "Macrophage/Monocyte",
           "Monocyte non-conventional", "Myeloid dendritic cell",
           "Plasmacytoid dendritic cell", "Dendritic cells resting",
           "Mast cells activated", "Eosinophils", "Neutrophil", "Basophil")
  expect_true(all(deconv_to_lineage(lab) == "Immune_other"))
})

test_that("plasma cells stay out: CD45-dim on IHC, so the tree may not call them immune", {
  expect_true(all(is.na(deconv_to_lineage(c("Plasma cells", "B cell plasma immature")))))
  # ...but 'plasmacytoid' is a dendritic cell, not a plasma cell
  expect_equal(deconv_to_lineage("Plasmacytoid dendritic cell"), "Immune_other")
})

test_that("Immune_other is a sum, so only an additive method gets one", {
  lm22 <- c("B cells naive", "Macrophages M2", "Neutrophils", "Plasma cells", "T cells CD8")
  c <- deconv_contributions(.ct("cibersortx", lm22))
  expect_setequal(c$cell_type[c$lineage == "Immune_other"], lm22[1:3])
  expect_true(all(c$route[c$lineage == "Immune_other"] == "sum"))
  expect_false("Plasma cells" %in% c$cell_type)

  m <- deconv_contributions(.ct("mcp_counter", c("B cell", "Monocyte", "T cell")))
  expect_false("Immune_other" %in% m$lineage)
})

test_that("Immune_total sums every immune label, only on the all-cells denominator", {
  q <- deconv_contributions(.ct("quantiseq", c("T cell CD8+", "NK cell", "B cell",
                                               "Macrophage M2", "uncharacterized cell")))
  expect_setequal(q$cell_type[q$lineage == "Immune_total"],
                  c("T cell CD8+", "NK cell", "B cell", "Macrophage M2"))
  # CIBERSORTx sums to 1 over leukocytes: its Immune_total is constant, so skipped
  x <- deconv_contributions(.ct("cibersortx", c("T cells CD8", "B cells naive")))
  expect_false("Immune_total" %in% x$lineage)
})

test_that("every compared population has a level from the gate tree", {
  expect_true(all(deconv_comparison_lineages$level %in%
                    c("macro", "intermediate", "leaf", "cross-cutting")))
  lv <- stats::setNames(deconv_comparison_lineages$level, deconv_comparison_lineages$lineage)
  expect_equal(unname(lv[c("Immune_total", "NonImmune")]), c("macro", "macro"))
  expect_equal(unname(lv[c("CD8T", "CD4T", "Treg")]), rep("leaf", 3))
})

test_that("Immune_other on IHC is the CD3- Immune cells plus mirage's myeloid calls", {
  d <- rbind(.ihc(), tibble::tibble(patient_id = "046",
                                    phenotype_clean = c("Myeloid", "Macrophage_M2"),
                                    CD3_sign = "-"))
  f <- ihc_comparison_fraction(d)
  n <- stats::setNames(f$n, f$lineage)
  expect_equal(n[["Immune_other"]], 3)          # one CD3- Immune + two myeloid
  expect_equal(n[["T_total"]], 6)               # unchanged: the CD3+ Immune only
  # Immune_total is the whole CD45+ branch, i.e. the CD45 denominator itself
  expect_equal(n[["Immune_total"]], f$n_cd45[1])
})

test_that("with CD3 never gated, every Immune cell is Immune_other and none is T", {
  d <- .ihc(); d$CD3_sign <- NULL
  f <- ihc_comparison_fraction(d)
  expect_equal(f$n[f$lineage == "Immune_other"], 2)
  expect_equal(f$n[f$lineage == "T_total"], 5)
})
