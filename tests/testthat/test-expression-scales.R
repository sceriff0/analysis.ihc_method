# expression_scales.R is base R over a gene x sample matrix, so it is testable with
# no DESeq2 object and no data on disk: the dds-reading branches are exercised only
# where DESeq2 is installed, the arithmetic and the file readers everywhere.
source(here::here("code", "expression_scales.R"))

counts <- matrix(c(10, 20, 30,  0, 5, 50), nrow = 3,
                 dimnames = list(c("A", "B", "C"), c("s1", "s2")))

test_that("TPM columns sum to a million and scale inversely with gene length", {
  tpm <- tpm_from_counts(counts, c(A = 1000, B = 2000, C = 500))
  expect_equal(unname(colSums(tpm)), c(1e6, 1e6))
  # Same count, twice the length -> half the rate: A (10 / 1 kb) vs B (20 / 2 kb).
  expect_equal(tpm["A", "s1"], tpm["B", "s1"])
  # A zero-count gene is zero, never NaN.
  expect_equal(tpm["A", "s2"], 0)
})

test_that("a per-sample length matrix (tximport's avgTxLength) is accepted", {
  len <- matrix(c(1000, 2000, 500, 1000, 4000, 500), nrow = 3,
                dimnames = dimnames(counts))
  tpm <- tpm_from_counts(counts, len)
  expect_equal(unname(colSums(tpm)), c(1e6, 1e6))
  # s2's B is twice as long as s1's, so B's rate halves in s2 and C's share rises
  # relative to a vector-length run — the matrix is used per sample, not collapsed.
  vec <- tpm_from_counts(counts, c(A = 1000, B = 2000, C = 500))
  expect_equal(tpm[, "s1"], vec[, "s1"])
  expect_gt(tpm["C", "s2"], vec["C", "s2"])
})

test_that("genes with no length are dropped with a message, not silently zeroed", {
  expect_message(tpm <- tpm_from_counts(counts, c(A = 1000, B = 2000)), "no length")
  expect_setequal(rownames(tpm), c("A", "B"))
})

test_that("a gene_lengths.tsv on disk is the last-resort length source", {
  f <- tempfile(fileext = ".tsv")
  writeLines(c("gene\tlength", "A\t1000", "B\t2000", "ZZZ\t9"), f)
  gl <- gene_lengths_from_dds(NULL, genes = rownames(counts), lengths_file = f)
  expect_equal(gl$source, "data/gene_lengths.tsv")
  expect_equal(gl$lengths[c("A", "B")], c(A = 1000, B = 2000))
  expect_true(is.na(gl$lengths[["C"]]))
  # Nothing anywhere -> NULL lengths and a source that says so.
  none <- gene_lengths_from_dds(NULL, genes = rownames(counts), lengths_file = tempfile())
  expect_null(none$lengths)
})

test_that("expression_scales() on a bare matrix gives raw and library sizes, and TPM only with lengths", {
  sc <- expression_scales(counts, samples = "s1", lengths = NULL)
  expect_equal(colnames(sc$raw), "s1")
  expect_equal(unname(sc$library_size), 60)
  expect_null(sc$tpm)
  sc2 <- expression_scales(counts, lengths = c(A = 1000, B = 2000, C = 500))
  expect_equal(names(sc2$library_size), c("s1", "s2"))
  expect_equal(unname(colSums(sc2$tpm)), c(1e6, 1e6))
  # A bare matrix has no size factors, so `normalized` is absent, not fabricated.
  expect_null(sc2$normalized)
})

test_that("scale_long() stacks the available scales and names them", {
  sc <- expression_scales(counts, lengths = c(A = 1000, B = 2000, C = 500))
  long <- scale_long(sc)
  expect_setequal(unique(long$scale), c("raw", "tpm"))
  expect_equal(nrow(long), 2 * length(counts))
  expect_true(all(c("scale", "sample", "gene", "value") %in% names(long)))
})

# --- CIBERSORTx: the mixture we hand it and the results it hands back ---------
test_that("the mixture file is tab-delimited, genes first, one column per sample", {
  f <- tempfile(fileext = ".txt")
  write_cibersortx_mixture(counts, f)
  L <- readLines(f)
  expect_equal(L[1], "GeneSymbol\ts1\ts2")
  expect_equal(L[2], "A\t10\t0")
  expect_false(any(grepl('"', L)))                    # no quoting: the container chokes
})

test_that("read_cibersortx_results() prefers the batch-corrected file and drops the QC columns", {
  d <- tempfile(); dir.create(d)
  hdr <- "Mixture\tT cells CD8\tT cells regulatory (Tregs)\tNK cells resting\tP-value\tCorrelation\tRMSE"
  writeLines(c(hdr, "s1\t0.5\t0.2\t0.3\t0.01\t0.9\t0.4", "s2\t0.1\t0.1\t0.8\t0.2\t0.5\t0.9"),
             file.path(d, "CIBERSORTx_Results.txt"))
  writeLines(c(hdr, "s1\t0.6\t0.1\t0.3\t0.01\t0.9\t0.4", "s2\t0.2\t0.1\t0.7\t0.2\t0.5\t0.9"),
             file.path(d, "CIBERSORTx_Adjusted.txt"))
  r <- read_cibersortx_results(d)
  expect_equal(unique(r$method), "cibersortx")
  expect_equal(unique(r$source_file), "CIBERSORTx_Adjusted.txt")
  expect_setequal(unique(r$cell_type), c("T cells CD8", "T cells regulatory (Tregs)", "NK cells resting"))
  expect_equal(r$score[r$sample == "s1" & r$cell_type == "T cells CD8"], 0.6)
  expect_false(any(c("P-value", "Correlation", "RMSE") %in% r$cell_type))
  # No output yet -> empty frame, not an error: the page must knit before the run.
  expect_equal(nrow(read_cibersortx_results(tempfile())), 0)
})

test_that("LM22's cell-type names reach the four comparable lineages through deconv_to_lineage()", {
  source(here::here("code", "validation_helpers.R"))
  lm22 <- c("T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting",
            "T cells CD4 memory activated", "T cells regulatory (Tregs)",
            "NK cells resting", "NK cells activated", "T cells follicular helper",
            "T cells gamma delta", "Macrophages M2")
  got <- deconv_to_lineage(lm22)
  expect_equal(got[1:7], c("CD8T", "CD4T", "CD4T", "CD4T", "Treg", "NK", "NK"))
  expect_true(all(is.na(got[8:10])))
})
