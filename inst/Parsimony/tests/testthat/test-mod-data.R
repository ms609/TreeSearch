library(shiny)
library(shinyjs)

# Stub globals that the module references from global.R
# Must be in global env BEFORE source() so debounce() can find them
palettes <- list("#cc9966")
Notification <- function(...) invisible(NULL)
logging <- FALSE
aJiffy <- 42
typingJiffy <- 105
aFewTrees <- 48L
fastSampleAbove <- 1000L
LogMsg <- function(...) invisible(NULL)

# Source the module under test (relative to tests/testthat/)
source("../../server/mod_data.R", local = TRUE)

# Stub logging functions
noop <- function(...) invisible(NULL)
stub_log_fns <- list(
  LogMsg     = noop,
  LogComment = noop,
  LogCode    = noop,
  CacheInput = noop,
  LastFile   = function(type) "stub.txt"
)

# Stub callbacks
stub_callbacks <- list(
  DisplayTreeScores       = noop,
  UpdateKeepNTipsRange    = noop,
  UpdateDroppedTaxaDisplay = noop,
  UpdateOutgroupInput     = noop
)

test_that("data_server returns expected reactive list", {
  r <- reactiveValues(
    trees       = NULL,
    allTrees    = NULL,
    treeHash    = NULL,
    dataset     = NULL,
    dataHash    = NULL,
    nTree       = 1L,
    treeRange   = c(1L, 1L),
    oldNTree    = NULL,
    oldTreeRange = NULL,
    updatingTrees = FALSE,
    newTrees    = NULL,
    dataFileVisible = FALSE,
    chars       = NULL,
    charNotes   = NULL,
    readDataFile = NULL,
    sortTrees   = FALSE,
    bestSearchScore = NULL,
    excelFiles  = 0,
    searchTotalHits = 0L,
    searchTotalReps = 0L,
    searchWithout = character(0),
    visibleConfigs = character(0),
    outgroup    = NULL,
    keepNTips   = 0L
  )

  shiny::testServer(
    data_server,
    args = list(
      r              = r,
      parent_session = NULL,
      callbacks      = stub_callbacks,
      log_fns        = stub_log_fns
    ),
    {
      returned <- session$getReturned()
      expect_true(is.list(returned))
      expect_true(all(c("AnyTrees", "HaveData", "tipLabels", "nChars",
                         "TaxonOrder", "DatasetMatchesTrees",
                         "UpdateAllTrees", "UpdateActiveTrees",
                         "dataSource") %in% names(returned)))
    }
  )
})

make_data_state <- function(...) {
  reactiveValues(
    trees       = NULL,
    allTrees    = NULL,
    treeHash    = NULL,
    dataset     = NULL,
    dataHash    = NULL,
    nTree       = 1L,
    treeRange   = c(1L, 1L),
    oldNTree    = NULL,
    oldTreeRange = NULL,
    updatingTrees = FALSE,
    newTrees    = NULL,
    dataFileVisible = FALSE,
    chars       = NULL,
    charNotes   = NULL,
    readDataFile = NULL,
    sortTrees   = FALSE,
    bestSearchScore = NULL,
    excelFiles  = 0,
    searchTotalHits = 0L,
    searchTotalReps = 0L,
    searchWithout = character(0),
    visibleConfigs = character(0),
    outgroup    = NULL,
    keepNTips   = 0L,
    ...
  )
}

test_that("AnyTrees and HaveData respond to state", {
  r <- make_data_state()

  shiny::testServer(
    data_server,
    args = list(
      r              = r,
      parent_session = NULL,
      callbacks      = stub_callbacks,
      log_fns        = stub_log_fns
    ),
    {
      returned <- session$getReturned()
      expect_false(returned$AnyTrees())
      expect_false(returned$HaveData())

      r$trees <- ape::rmtree(3, 6)
      session$flushReact()
      expect_true(returned$AnyTrees())
    }
  )
})

test_that("tipLabels returns tree tips when trees present", {
  r <- make_data_state()
  trees <- ape::rmtree(3, 6)

  shiny::testServer(
    data_server,
    args = list(r = r, parent_session = NULL,
                callbacks = stub_callbacks, log_fns = stub_log_fns),
    {
      returned <- session$getReturned()
      # No trees => tipLabels is an empty character vector (mod_data.R guards
      # `if (!length(r$trees)) return(character(0L))` -- safer than NULL for the
      # downstream intersect()/choices= consumers).
      expect_equal(returned$tipLabels(), character(0L))

      r$trees <- trees
      session$flushReact()
      expect_equal(returned$tipLabels(), trees[[1]]$tip.label)
    }
  )
})

test_that("nChars returns 0 when no dataset", {
  r <- make_data_state()

  shiny::testServer(
    data_server,
    args = list(r = r, parent_session = NULL,
                callbacks = stub_callbacks, log_fns = stub_log_fns),
    {
      returned <- session$getReturned()
      expect_equal(returned$nChars(), 0L)
    }
  )
})

test_that("HaveData requires phyDat class", {
  r <- make_data_state()
  r$dataset <- list(a = 1)  # Not a phyDat object

  shiny::testServer(
    data_server,
    args = list(r = r, parent_session = NULL,
                callbacks = stub_callbacks, log_fns = stub_log_fns),
    {
      returned <- session$getReturned()
      expect_false(returned$HaveData())
    }
  )
})

# Helper: minimal fake phyDat (class only; enough for HaveData + DatasetMatchesTrees)
fake_phyDat <- function(tip_names) {
  structure(setNames(vector("list", length(tip_names)), tip_names),
            class = "phyDat")
}

test_that("dataset observer preserves compatible trees (T-151 regression)", {
  # Reproduce the blank-plot bug: UpdateData() sets r$dataset AND r$allTrees
  # in the same reactive flush; the observeEvent(r$dataset,...) must NOT clear
  # trees that are already compatible with the new dataset.
  r <- make_data_state()
  tips <- paste0("t", 1:6)

  shiny::testServer(
    data_server,
    args = list(r = r, parent_session = NULL,
                callbacks = stub_callbacks, log_fns = stub_log_fns),
    {
      returned <- session$getReturned()

      trees <- ape::rmtree(3, 6)
      for (i in seq_along(trees)) trees[[i]]$tip.label <- tips

      # Simulate UpdateData(): allTrees set first, then dataset triggers observer
      r$allTrees <- trees
      r$trees    <- trees
      r$dataset  <- fake_phyDat(tips)
      session$flushReact()

      # Trees must survive the observer — fix for T-151
      expect_equal(length(r$allTrees), 3L,
                   label = "compatible trees cleared by dataset observer (T-151)")
    }
  )
})

test_that("UpdateActiveTrees does not crash when allTrees cleared with stale nTree (T-292)", {
  # Regression: UpdateAllTrees(list()) clears r$allTrees but does NOT reset
  # r$nTree (the `if (nTrees > 0L)` guard skips it). UpdateActiveTrees() then
  # reaches the else branch, computes r$allTrees[1:old_nTree] = list of NULLs,
  # and WideSample() throws "trees must be a multiPhylo object".
  r <- make_data_state()
  tips <- paste0("t", 1:6)
  trees <- ape::rmtree(5, 6)
  for (i in seq_along(trees)) trees[[i]]$tip.label <- tips

  shiny::testServer(
    data_server,
    args = list(r = r, parent_session = NULL,
                callbacks = stub_callbacks, log_fns = stub_log_fns),
    {
      returned <- session$getReturned()

      # Simulate post-search state: 5 trees
      r$allTrees  <- trees
      r$trees     <- trees
      r$nTree     <- 5L
      r$treeRange <- c(1L, 5L)
      session$flushReact()

      # Simulate UpdateAllTrees(list()): clears allTrees but leaves nTree = 5
      r$allTrees <- list()
      session$flushReact()

      expect_no_error(returned$UpdateActiveTrees())
      expect_null(r$trees)
    }
  )
})

test_that("dataset observer clears incompatible trees on dataset switch", {
  # When old trees have different taxa than the new dataset, they must be cleared.
  r <- make_data_state()
  tips_new <- paste0("x", 1:5)  # Different from tree tips

  shiny::testServer(
    data_server,
    args = list(r = r, parent_session = NULL,
                callbacks = stub_callbacks, log_fns = stub_log_fns),
    {
      returned <- session$getReturned()

      old_trees <- ape::rmtree(2, 6)  # 6-taxon trees: tip labels t1..t6
      r$allTrees <- old_trees
      r$trees    <- old_trees

      # Switch to a 5-taxon dataset with completely different names
      r$dataset <- fake_phyDat(tips_new)
      session$flushReact()

      expect_null(r$allTrees,
                  label = "incompatible trees should be cleared on dataset switch")
    }
  )
})

test_that("PrethinIndex() and WideSampleCode() agree", {
  expect_null(PrethinIndex(1000, 48, maxTrees = 1000))
  expect_null(PrethinIndex(5000, 6000, maxTrees = 1000))
  expect_equal(PrethinIndex(10000, 48, maxTrees = 1000),
               round(seq(1, 10000, length.out = 1000)))
  expect_equal(anyDuplicated(PrethinIndex(1001, 48, maxTrees = 1000)), 0)
  # Asking for more trees than the cap keeps that many, not the cap
  expect_length(PrethinIndex(10000, 5000, maxTrees = 1000), 5000)

  # Small ranges log the unchanged pre-fix expression (download snapshots)
  expect_equal(WideSampleCode("trees", c(77, 125), 48, maxTrees = 1000),
               "trees <- WideSample(trees[77:125], 48)")
  code <- WideSampleCode("allTrees", c(11, 10010), 48, maxTrees = 1000)
  expect_equal(code, paste0("trees <- WideSample(allTrees[11:10010]",
                            "[round(seq(1, 10000, length.out = 1000))], 48)"))
  allTrees <- 1:20000
  WideSample <- function(x, n) x
  eval(parse(text = code))
  expect_equal(trees, (11:10010)[PrethinIndex(10000, 48, maxTrees = 1000)])
})

# Upload a tree file through input$treeFile, recording notifications and
# logged code.  Large tree sets froze EasyTrees() for minutes with no message:
# WideSample() built a 10 000 x 10 000 distance matrix (2026-11 report).
upload_trees <- function(trees, maxTrees, wideSample = NULL) {
  modEnv <- environment(data_server)
  oldMax <- get("fastSampleAbove", modEnv)
  oldNotification <- get("Notification", modEnv)
  on.exit({
    assign("fastSampleAbove", oldMax, envir = modEnv)
    assign("Notification", oldNotification, envir = modEnv)
    if (exists("WideSample", envir = modEnv, inherits = FALSE)) {
      rm("WideSample", envir = modEnv)
    }
  })
  rec <- new.env()
  rec$notes <- list()
  rec$code <- character(0)
  assign("fastSampleAbove", maxTrees, envir = modEnv)
  assign("Notification", function(ui, ..., type = "default") {
    rec$notes[[length(rec$notes) + 1L]] <- list(text = ui, type = type)
  }, envir = modEnv)
  if (!is.null(wideSample)) {
    assign("WideSample", wideSample, envir = modEnv)
  }
  log_fns <- stub_log_fns
  log_fns$LogCode <- function(...) rec$code <- c(rec$code, ...)

  treeFile <- tempfile(fileext = ".nex")
  on.exit(unlink(treeFile), add = TRUE)
  ape::write.nexus(trees, file = treeFile)

  r <- make_data_state(readTreeFile = NULL, thinningSeed = NULL)
  shiny::testServer(
    data_server,
    args = list(r = r, parent_session = shiny::MockShinySession$new(),
                callbacks = stub_callbacks, log_fns = log_fns),
    {
      session$setInputs(nTree = 48L)
      session$setInputs(treeFile = data.frame(name = "trees.nex",
                                              datapath = treeFile))
      rec$allTrees <- r$allTrees
      rec$trees <- r$trees
      rec$newTrees <- r$newTrees
      rec$readTreeFile <- r$readTreeFile
    }
  )
  as.list(rec)
}

dotted_trees <- function(nTree) {
  # Punctuated labels as in the reporting MorphoBank matrix
  tips <- c("Cf._Peytoia", "Cricocosmia_n._sp.", "Scolecofurca_",
            "MAOTIANSHANIA_cylindrica", paste0("t", 5:9))
  set.seed(0)
  trees <- ape::rmtree(nTree, length(tips), tip.label = tips)
  for (i in seq_along(trees)) trees[[i]][["edge.length"]] <- NULL
  trees
}

test_that("large uploaded tree sets are thinned from evenly spaced trees", {
  trees <- dotted_trees(60)
  out <- upload_trees(trees, maxTrees = 20L)

  types <- vapply(out$notes, `[[`, "", "type")
  texts <- vapply(out$notes, `[[`, "", "text")
  expect_false("error" %in% types)
  expect_equal(out$readTreeFile, "read.nexus(treeFile)")
  expect_length(out$allTrees, 60)
  expect_setequal(out$allTrees[[1]]$tip.label, trees[[1]]$tip.label)
  # nTree (48) exceeds maxTrees (20), so 48 evenly spaced trees are kept
  expect_length(out$trees, 48)
  keep <- round(seq(1, 60, length.out = 48))
  expect_equal(unname(out$trees), unname(out$allTrees[keep]))
  expect_true(any(grepl("Sampled 48 trees from 48 evenly spaced among 60",
                        texts, fixed = TRUE)))
  expect_true(any(grepl("Loaded 60 trees (48 sampled)", texts, fixed = TRUE)))

  # Logged code reproduces the selection
  logged <- grep("WideSample(", out$code, fixed = TRUE, value = TRUE)
  expect_equal(logged, paste0("trees <- WideSample(allTrees[1:60]",
                              "[round(seq(1, 60, length.out = 48))], 48)"))
  allTrees <- out$allTrees
  eval(parse(text = logged))
  expect_equal(unname(trees), unname(out$trees))
})

test_that("tree sets within the cap are thinned from every tree", {
  out <- upload_trees(dotted_trees(60), maxTrees = 1000L)
  expect_false("error" %in% vapply(out$notes, `[[`, "", "type"))
  expect_length(out$trees, 48)
  expect_equal(grep("WideSample(", out$code, fixed = TRUE, value = TRUE),
               "trees <- WideSample(allTrees[1:60], 48)")
})

test_that("errors processing uploaded trees are shown, not swallowed", {
  out <- upload_trees(dotted_trees(60), maxTrees = 1000L,
                      wideSample = function(...) stop("boom"))
  types <- vapply(out$notes, `[[`, "", "type")
  texts <- vapply(out$notes, `[[`, "", "text")
  expect_equal(types, "error")
  expect_match(texts, "Read 60 trees, but could not process them: boom",
               fixed = TRUE)
  # Forgotten, so that a re-upload is not skipped as "unchanged"
  expect_null(out$newTrees)
})
