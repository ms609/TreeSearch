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

test_that("AnyTrees and HaveData respond to state", {
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
      expect_false(returned$AnyTrees())
      expect_false(returned$HaveData())

      r$trees <- ape::rmtree(3, 6)
      session$flushReact()
      expect_true(returned$AnyTrees())
    }
  )
})
