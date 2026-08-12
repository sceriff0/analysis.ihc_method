# A `child=` chunk path is resolved by knitr against the KNIT WORKING DIRECTORY, not
# against the parent document. workflowr sets that to the project root
# (_workflowr.yml: knit_root_dir: "."), so a bare "_body.Rmd" written next to its
# parent in analysis/ is looked for at <root>/_body.Rmd and the build dies with
#   Error in file(con, "r") : cannot open the connection
# partway through the parent — readLines()'s internal call, with nothing naming the
# child. Every child reference must therefore be absolute via here::here().

# Parse the chunk header as R rather than regexing it: the child path is commonly a
# call with its own commas (here::here("analysis", "x.Rmd")), which a comma-splitting
# regex truncates into a syntax error.
parent_child_refs <- function() {
  refs <- list()
  for (f in list.files(here::here("analysis"), pattern = "[.]Rmd$", full.names = TRUE)) {
    for (line in grep("^```\\{r.*\\bchild\\s*=", readLines(f, warn = FALSE), value = TRUE)) {
      header <- sub("^```\\{r\\s*", "", sub("\\}\\s*$", "", line))   # "body, child = ..."
      opts <- tryCatch(eval(parse(text = paste0("alist(", header, ")"))),
                       error = function(e) NULL)
      if (is.null(opts) || !"child" %in% names(opts)) next
      refs[[length(refs) + 1]] <- list(parent = basename(f), expr = opts$child)
    }
  }
  refs
}

test_that("every child= path resolves from the knit root, not the parent's directory", {
  refs <- parent_child_refs()
  skip_if(length(refs) == 0, "no child chunks in analysis/")
  # Resolve exactly as knitr would: evaluate the option, then read it with the
  # working directory set to the knit root.
  withr::with_dir(here::here(), {
    for (r in refs) {
      path <- eval(r$expr)
      expect_true(file.exists(path),
                  info = paste0(r$parent, " -> ", deparse1(r$expr), " resolves to '", path,
                                "', which does not exist from the knit root"))
    }
  })
})

test_that("child paths are absolute, so the parent also knits from analysis/", {
  # RStudio's Knit button uses the document's own directory. A path that happens to
  # work from the project root ("analysis/_body.Rmd") would break there; here::here()
  # works from both.
  for (r in parent_child_refs())
    expect_match(deparse1(r$expr), "here::here|here\\(",
                 info = paste(r$parent, "child path is not anchored with here::here()"))
})

test_that("a bare child path really does fail under workflowr's knit root", {
  # The regression itself, reproduced end to end — this is what the build hit.
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  root <- file.path(tempdir(), paste0("childrepro-", sample(1e6, 1)))
  dir.create(file.path(root, "analysis"), recursive = TRUE)
  writeLines("child text", file.path(root, "analysis", "_body.Rmd"))
  writeLines(c("---", "title: p", "---", "", '```{r body, child = "_body.Rmd"}', "```"),
             file.path(root, "analysis", "bare.Rmd"))
  writeLines(c("---", "title: p", "---", "",
               sprintf('```{r body, child = "%s"}', file.path(root, "analysis", "_body.Rmd")),
               "```"),
             file.path(root, "analysis", "absolute.Rmd"))

  # The bare case is EXPECTED to warn on the missing file before it errors — that is
  # the bug being reproduced, not a problem with the test.
  render_from_root <- function(f) suppressWarnings(tryCatch({
    rmarkdown::render(file.path(root, "analysis", f), knit_root_dir = root,
                      output_dir = tempdir(), quiet = TRUE)
    "ok"
  }, error = function(e) conditionMessage(e)))

  expect_match(render_from_root("bare.Rmd"), "cannot open the connection")
  expect_equal(render_from_root("absolute.Rmd"), "ok")
})

# The regression the reorganisation was triggered by: `wflow_build("analysis/*.Rmd")`
# crashed because Sys.glob() handed it the shared clinical body as if it were a page.
# render_site()'s `^[_.]` rule had hidden it from the no-argument build, so the file
# looked correctly excluded right up until someone globbed explicitly.
test_that("no child document is reachable by the glob wflow_build expands", {
  globbed <- basename(Sys.glob(here::here("analysis", "*.Rmd")))
  expect_false(any(grepl("^_", globbed)),
               info = paste("underscored file(s) reachable by analysis/*.Rmd:",
                            paste(grep("^_", globbed, value = TRUE), collapse = ", "),
                            "- glob's `*` refuses a leading dot, not a leading",
                            "underscore. Move children to analysis/_children/."))
  # And every child that a parent references really does live outside that glob.
  for (r in parent_child_refs())
    expect_false(basename(dirname(eval(r$expr))) == "analysis",
                 info = paste(r$parent, "includes a child that sits directly in",
                              "analysis/, where wflow_build(\"analysis/*.Rmd\") will",
                              "try to build it as a page"))
})

# Every page writes its vector figures to output/figures/<SLUG>/. If SLUG drifts from
# the filename, a renamed page silently keeps exporting under its old name and the
# figure directory stops matching the site.
test_that("each page's SLUG equals its own filename", {
  for (f in Sys.glob(here::here("analysis", "*.Rmd"))) {
    slug <- grep("^SLUG\\s*<-", readLines(f, warn = FALSE), value = TRUE)
    if (length(slug) == 0) next            # about / license / index define none
    expect_equal(sub('^SLUG\\s*<-\\s*"([^"]+)".*$', "\\1", slug[1]),
                 sub("[.]Rmd$", "", basename(f)),
                 info = paste(basename(f), "SLUG does not match its filename"))
  }
})
