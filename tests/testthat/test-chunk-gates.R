# GATED CHUNKS MUST NOT STRAND THEIR DEPENDENTS.
#
# Several pages skip their analysis sections when the chosen arm has no cells on
# disk (`eval = have_arm`, `have_qc`, `have_margin` — the pattern this project
# already used for the invasive-margin section). That is what turns a missing data
# tree from a crash into a page that explains itself.
#
# It also creates a hazard that is invisible on a machine WITH data: a chunk left
# ungated that reads an object only a gated chunk creates. The knit then dies with
# "object '<name>' not found" — and a defensive `if (!is.null(x))` in the dependent
# chunk does NOT save it, because the name was never bound at all. That is exactly
# how deconv-cache-paired broke after the deconvolution section was gated.
#
# This test is the static check: for every page, an ungated chunk may not read a
# name that only gated chunks bind at top level.
.rmd_chunks <- function(path) {
  L <- readLines(path)
  starts <- grep("^```\\{r", L)
  lapply(starts, function(s) {
    e   <- s + which(L[(s + 1):length(L)] == "```")[1]
    hdr <- L[s]
    list(label = trimws(sub("^```\\{r[ ,]*([^,}]*).*$", "\\1", hdr)),
         gate  = if (grepl("eval\\s*=\\s*have_", hdr))
                   sub(".*eval\\s*=\\s*(have_[a-z_]+).*", "\\1", hdr) else NA_character_,
         code  = if (e > s + 1) L[(s + 1):(e - 1)] else character(0))
  })
}

# TOP-LEVEL bindings only. A name appearing as an NSE column inside group_by() or
# data.frame() is not a binding, and counting it produces false positives — which
# is why this walks the parse tree instead of matching text.
.top_level_binds <- function(code) {
  p <- tryCatch(parse(text = paste(code, collapse = "\n")), error = function(e) NULL)
  if (is.null(p)) return(character(0))
  out <- character(0)
  for (e in as.list(p)) {
    if (!is.call(e)) next
    head1 <- as.character(e[[1]])[1]
    if (head1 %in% c("<-", "=", "<<-") && is.name(e[[2]]))
      out <- c(out, as.character(e[[2]]))
    if (head1 == "for") out <- c(out, as.character(e[[2]]))
  }
  unique(out)
}

.reads <- function(code) {
  p <- tryCatch(parse(text = paste(code, collapse = "\n")), error = function(e) NULL)
  if (is.null(p)) character(0) else unique(all.vars(p))
}

PAGES <- c("clinical_flowpath", "clinical_membership_qc", "molecular_hot_cold",
           "marker_qc", "paper_figures")

for (page in PAGES) local({
  pg   <- page
  path <- here::here("analysis", paste0(pg, ".Rmd"))

  test_that(paste0(pg, ": no ungated chunk reads a gated-only object"), {
    skip_if_not(file.exists(path))
    ch     <- .rmd_chunks(path)
    gated  <- Filter(function(c) !is.na(c$gate), ch)
    open   <- Filter(function(c)  is.na(c$gate), ch)
    skip_if(length(gated) == 0, "page has no gated chunks")

    gated_binds <- unique(unlist(lapply(gated, function(c) .top_level_binds(c$code))))
    open_binds  <- unique(unlist(lapply(open,  function(c) .top_level_binds(c$code))))
    only_gated  <- setdiff(gated_binds, open_binds)

    offenders <- character(0)
    for (c in open) {
      hit <- intersect(.reads(c$code), only_gated)
      if (length(hit))
        offenders <- c(offenders, sprintf("%s reads %s", c$label, paste(hit, collapse = ", ")))
    }
    expect_equal(offenders, character(0),
                 info = paste0("\nUngated chunk(s) depending on a gated chunk in ", pg,
                               ".Rmd:\n  ", paste(offenders, collapse = "\n  "),
                               "\nGate them too, or move the binding to an ungated chunk."))
  })

  test_that(paste0(pg, ": every gate is defined before the chunks that use it"), {
    skip_if_not(file.exists(path))
    ch <- .rmd_chunks(path)
    skip_if(!any(!is.na(vapply(ch, function(c) c$gate, character(1)))), "no gates")
    seen <- character(0)
    for (c in ch) {
      if (!is.na(c$gate))
        expect_true(c$gate %in% seen,
                    info = paste0(pg, ": chunk '", c$label, "' is gated on '", c$gate,
                                  "', which no earlier ungated chunk defines"))
      if (is.na(c$gate)) seen <- c(seen, .top_level_binds(c$code))
    }
  })
})
