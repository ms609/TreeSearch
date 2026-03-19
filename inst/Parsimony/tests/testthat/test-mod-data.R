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
      # No trees => tipLabels is NULL (NULL[[1]][["tip.label"]])
      expect_null(returned$tipLabels())

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
