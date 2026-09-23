# =============================================================================
# expression_scales.R  —  the bulk-RNA matrix on its three scales
#
# The DESeq2 object holds RAW counts and, after DESeq(), size factors. Every page
# used to read ONE derived scale (size-factor-normalised counts) and nothing else,
# so nobody could see what the normalisation had done, and the deconvolution
# methods that ask for TPM (quanTIseq, EPIC, CIBERSORT) were fed something else
# with a caveat. This file derives all three from the same object:
#
#   raw         counts(dds)                     what the aligner counted
#   normalized  counts(dds, normalized = TRUE)  raw / size factor — library depth out
#   tpm         transcripts per million         raw / gene length, then per-million —
#                                               depth AND length out; the scale the
#                                               signature-based methods were fit on
#
# TPM NEEDS GENE LENGTHS AND A DESeqDataSet DOES NOT PROMISE THEM. Where they can
# live, in the order tried:
#   1. assays(dds)$avgTxLength     tximport's per-sample effective lengths
#   2. mcols(dds)$basepairs        what DESeq2::fpkm() reads
#   3. rowRanges(dds)              exon widths, when the object came from
#                                  summarizeOverlaps() or carried a GRangesList
#   4. data/gene_lengths.tsv       two columns, gene and length (bp), that YOU make
#                                  from the GTF the counts were made against
# Without any of these `tpm` is NULL and the page says so; it does not fall back to
# a made-up length, because a TPM on the wrong lengths is worse than none.
#
# Base R + tibble, like cell_tables.R: testable on a plain matrix with no DESeq2.
# DESeq2 / SummarizedExperiment are touched only when handed a real dds.
# =============================================================================

# TPM from a gene x sample count matrix. `lengths` is a named vector (bp) or a
# matrix with the counts' dimnames (per-sample lengths). Genes with no length are
# dropped, and said so: a silent zero would read as "not expressed".
tpm_from_counts <- function(counts, lengths) {
  stopifnot(is.matrix(counts), !is.null(rownames(counts)))
  if (is.matrix(lengths)) {
    keep <- rownames(counts)[rownames(counts) %in% rownames(lengths)]
    len  <- lengths[keep, colnames(counts), drop = FALSE]
    ok   <- keep[rowSums(!is.finite(len) | len <= 0) == 0]
  } else {
    keep <- rownames(counts)[rownames(counts) %in% names(lengths)]
    ok   <- keep[is.finite(lengths[keep]) & lengths[keep] > 0]
  }
  dropped <- setdiff(rownames(counts), ok)
  if (length(dropped))
    message("tpm_from_counts(): ", length(dropped), " of ", nrow(counts),
            " genes have no length and are dropped from the TPM matrix")
  if (!length(ok)) return(counts[integer(0), , drop = FALSE])
  x   <- counts[ok, , drop = FALSE]
  kb  <- if (is.matrix(lengths)) lengths[ok, colnames(counts), drop = FALSE] / 1000
         else matrix(lengths[ok] / 1000, nrow = length(ok), ncol = ncol(x),
                     dimnames = dimnames(x))
  rate  <- x / kb
  total <- colSums(rate)
  total[total == 0] <- 1                      # an all-zero library stays all zero
  sweep(rate, 2, total, "/") * 1e6
}

# Where the gene lengths are, if anywhere. Returns list(lengths, source);
# `lengths` is NULL when nothing was found. `genes` orders and completes the vector
# read from the tsv so a missing gene is an NA the caller can see.
gene_lengths_from_dds <- function(dds, genes = NULL,
                                  lengths_file = here::here("data", "gene_lengths.tsv")) {
  if (!is.null(dds) && requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    if (is.null(genes)) genes <- rownames(dds)
    an <- SummarizedExperiment::assayNames(dds)
    if ("avgTxLength" %in% an) {
      m <- SummarizedExperiment::assay(dds, "avgTxLength")
      return(list(lengths = as.matrix(m), source = "assays(dds)$avgTxLength (tximport)"))
    }
    mc <- SummarizedExperiment::mcols(dds)
    if ("basepairs" %in% names(mc))
      return(list(lengths = stats::setNames(as.numeric(mc$basepairs), rownames(dds)),
                  source = "mcols(dds)$basepairs"))
    rr <- tryCatch(SummarizedExperiment::rowRanges(dds), error = function(e) NULL)
    if (!is.null(rr) && requireNamespace("GenomicRanges", quietly = TRUE)) {
      w <- if (inherits(rr, "GRangesList"))
             sum(GenomicRanges::width(GenomicRanges::reduce(rr)))
           else if (inherits(rr, "GRanges")) GenomicRanges::width(rr)
           else NULL
      if (!is.null(w) && any(w > 0))
        return(list(lengths = stats::setNames(as.numeric(w), rownames(dds)),
                    source = "rowRanges(dds) exon widths"))
    }
  }
  if (file.exists(lengths_file)) {
    tab <- utils::read.delim(lengths_file, check.names = FALSE, stringsAsFactors = FALSE)
    stopifnot(ncol(tab) >= 2)
    len <- stats::setNames(as.numeric(tab[[2]]), as.character(tab[[1]]))
    if (!is.null(genes)) len <- stats::setNames(len[genes], genes)
    return(list(lengths = len, source = "data/gene_lengths.tsv"))
  }
  list(lengths = NULL, source = "none")
}

# The three scales from one object. `x` is a DESeqDataSet or a plain count
# matrix (then `normalized` is NULL: a bare matrix has no size factors, and this
# function will not invent them). `lengths` = "auto" looks them up as above; pass
# a vector / matrix / NULL to override.
expression_scales <- function(x, samples = NULL, lengths = "auto") {
  if (is.matrix(x)) {
    raw <- x; normalized <- NULL
    if (identical(lengths, "auto")) lengths <- NULL
    src <- if (is.null(lengths)) "none" else "supplied"
  } else {
    stopifnot(requireNamespace("DESeq2", quietly = TRUE))
    raw        <- DESeq2::counts(x, normalized = FALSE)
    normalized <- tryCatch(DESeq2::counts(x, normalized = TRUE), error = function(e) NULL)
    if (identical(lengths, "auto")) {
      gl <- gene_lengths_from_dds(x); lengths <- gl$lengths; src <- gl$source
    } else src <- if (is.null(lengths)) "none" else "supplied"
  }
  if (!is.null(samples)) {
    keep <- intersect(colnames(raw), samples)
    if (!length(keep)) stop("expression_scales(): none of the requested samples is in the matrix")
    raw <- raw[, keep, drop = FALSE]
    if (!is.null(normalized)) normalized <- normalized[, keep, drop = FALSE]
  }
  storage.mode(raw) <- "double"
  tpm <- if (!is.null(lengths)) tpm_from_counts(raw, lengths) else NULL
  if (!is.null(tpm) && nrow(tpm) == 0) tpm <- NULL
  list(raw = raw, normalized = normalized, tpm = tpm,
       library_size  = colSums(raw),
       length_source = src,
       samples       = colnames(raw))
}

# Long form over whichever scales exist: (scale, sample, gene, value).
scale_long <- function(scales) {
  parts <- lapply(c("raw", "normalized", "tpm"), function(nm) {
    m <- scales[[nm]]
    if (is.null(m)) return(NULL)
    tibble::tibble(scale  = nm,
                   sample = rep(colnames(m), each = nrow(m)),
                   gene   = rep(rownames(m), times = ncol(m)),
                   value  = as.vector(m))
  })
  out <- do.call(rbind, parts[!vapply(parts, is.null, logical(1))])
  out$scale <- factor(out$scale, levels = c("raw", "normalized", "tpm"))
  out
}

# --- CIBERSORTx --------------------------------------------------------------
# The container is not an R package: the page WRITES the mixture it should read,
# the job runs outside R (code/cibersortx_fractions.sh), and the page READS the
# result table back on the next knit. Both halves live here so the format has
# one owner.

# The mixture file: tab-delimited, unquoted, genes down the first column, one
# column per sample. The container rejects quoted headers and scientific
# notation is fine.
write_cibersortx_mixture <- function(mat, path) {
  stopifnot(is.matrix(mat), !is.null(rownames(mat)))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  df <- data.frame(GeneSymbol = rownames(mat), mat, check.names = FALSE,
                   stringsAsFactors = FALSE)
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = TRUE)
  invisible(path)
}

# The fractions table back: prefer the batch-corrected `*_Adjusted.txt` a B-mode
# run writes over the plain `*_Results.txt`, drop the per-sample QC columns
# (P-value, Correlation, RMSE, Absolute score), and return one row per
# (sample, cell_type) with method = "cibersortx". Empty when there is no output
# yet, so the page knits before the job has run.
read_cibersortx_results <- function(dir) {
  empty <- tibble::tibble(cell_type = character(), sample = character(),
                          score = numeric(), method = character(),
                          source_file = character())
  if (!dir.exists(dir)) return(empty)
  files <- list.files(dir, pattern = "_(Adjusted|Results)\\.txt$", full.names = TRUE)
  files <- files[!grepl("_Mixtures_Adjusted\\.txt$", files)]
  if (!length(files)) return(empty)
  f <- if (any(grepl("_Adjusted\\.txt$", files))) files[grepl("_Adjusted\\.txt$", files)][1]
       else files[grepl("_Results\\.txt$", files)][1]
  tab <- utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  qc  <- names(tab) %in% c("P-value", "Correlation", "RMSE") |
         grepl("^Absolute score", names(tab))
  cell_cols <- names(tab)[-1][!qc[-1]]
  out <- tibble::tibble(
    cell_type   = rep(cell_cols, each = nrow(tab)),
    sample      = rep(as.character(tab[[1]]), times = length(cell_cols)),
    score       = as.numeric(unlist(tab[cell_cols], use.names = FALSE)),
    method      = "cibersortx",
    source_file = basename(f))
  out
}
