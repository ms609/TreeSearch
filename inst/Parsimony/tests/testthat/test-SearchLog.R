test_that("Search log workflow produces expected outputs", {
  app <- AppDriver$new(
    app_dir = "../../",
    seed = 0,
    load_timeout = 200000,
    shiny_args = list(test.mode = TRUE),
    name = "SearchLog"
  )
  on.exit(app$stop(), add = TRUE)

  # The legacy tests/shinytest/SearchLog.R snapshotted the saveZip download
  # (session R-script log) wholesale. That is NOT portable across a CI matrix:
  # the log embeds search-result specifics (e.g. `startTree <- trees[[23]]`, the
  # first-optimal index; `allTrees[1:150]`, the tree count) and `nThreads = 4L`
  # (the runner's core count), all of which vary by machine/run because the
  # parallel parsimony search is not byte-deterministic. So instead of a
  # brittle content snapshot we assert the DETERMINISTIC log-generation markers:
  # the app must map the chosen inputs to the correct generated code. This
  # ports v1's intent (regression-guard the session log) in a portable form.
  expect_log_contains <- function(markers) {
    path <- app$get_download("dl-saveZip")
    expect_true(file.exists(path) && file.size(path) > 0)
    log <- readLines(path)
    for (m in markers) {
      expect_true(any(grepl(m, log, fixed = TRUE)),
                  info = paste("session log missing marker:", m))
    }
  }
  # saveNwk / saveNex carry the search-result trees (same non-determinism):
  # assert structural validity, not content.
  expect_valid_trees <- function(output, reader) {
    path <- app$get_download(output)
    expect_true(file.exists(path) && file.size(path) > 0)
    trees <- reader(path)
    expect_true(inherits(trees, c("phylo", "multiPhylo")))
  }

  # --- EW search (implied weights off, sprint) ---
  app$set_inputs(`data-dataSource` = "Wills2012", timeout_ = 4000)
  app$click("search-searchConfig")
  app$wait_for_idle(timeout = 5000)
  app$set_inputs(`search-concavity` = 1.1)
  app$set_inputs(`search-epsilon` = 1)
  app$set_inputs(`search-implied.weights` = "off")
  app$set_inputs(`search-strategy` = "sprint")
  app$set_inputs(`search-maxReplicates` = 5)
  app$set_inputs(`search-targetHits` = 3)
  app$click("search-modalGo")
  # ExtendedTask returns immediately; poll exported counter
  app$wait_for_value(export = "searchCount",
                     ignore = list(NULL, 0L),
                     timeout = 120000)
  app$click("search-searchConfig")
  app$wait_for_idle(timeout = 5000)

  expect_log_contains(c(
    'library("TreeSearch")',
    'dataFile <- system.file("datasets/Wills2012.nex"',
    "Rogue::QuickRogue",
    "MaximizeParsimony(",
    'strategy = "sprint"',
    "concavity = Inf",       # implied weights off
    "targetHits = 3"
  ))
  expect_valid_trees("dl-saveNwk", ape::read.tree)

  # --- IW search (implied weights on, default) ---
  app$set_inputs(`search-implied.weights` = "on")
  app$set_inputs(`search-strategy` = "default")
  app$set_inputs(`search-maxReplicates` = 3)
  app$set_inputs(`search-targetHits` = 2)
  app$set_inputs(`search-epsilon` = 0)
  app$click("search-modalGo")
  app$wait_for_value(export = "searchCount",
                     ignore = list(NULL, 0L, 1L),
                     timeout = 200000)

  expect_log_contains(c(
    'strategy = "default"',
    "extended_iw = FALSE",   # implied weights on
    "targetHits = 2"
  ))
  expect_valid_trees("dl-saveNex", ape::read.nexus)
})
