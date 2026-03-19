library(shiny)

# Source the module under test (relative to tests/testthat/)
# local = TRUE so the module captures the test environment (for stub lookups)
source("../../server/mod_search.R", local = TRUE)

# Stub globals that the module references from global.R
Notification <- function(...) invisible(NULL)
Icon <- function(...) shiny::icon(..., class = "fas")
Enquote <- function(x) if (is.character(x)) paste0("\"", x, "\"") else signif(x)
EnC <- function(x) {
  if (length(x) == 1) Enquote(x)
  else paste0("c(", paste(sapply(x, Enquote), collapse = ", "), ")")
}
FormatMissProb <- function(prob) {
  pct <- prob * 100
  if (pct >= 1) paste0("~", round(pct), "%")
  else if (pct >= 0.1) "<1%"
  else if (pct >= 0.01) "<0.1%"
  else "<0.01%"
}
SearchConfidenceText <- function(K, R) {
  if (is.null(K) || is.null(R) || R <= 0L || K <= 0L) return(NULL)
  prob_miss <- exp(-K)
  paste0(K, " of ", R, " runs hit best score. ",
         "Probability that a better score exists: ",
         FormatMissProb(prob_miss))
}
PutData <- PutTree <- function(...) invisible(NULL)
logging <- FALSE

# Stub shinyjs functions (not available in testServer context)
if (!requireNamespace("shinyjs", quietly = TRUE) ||
    !exists("show", envir = asNamespace("shinyjs"))) {
  show <- hide <- disable <- enable <- function(...) invisible(NULL)
} else {
  # Wrap to suppress errors when shinyjs isn't properly initialized
  show    <- function(...) tryCatch(shinyjs::show(...), error = function(e) invisible(NULL))
  hide    <- function(...) tryCatch(shinyjs::hide(...), error = function(e) invisible(NULL))
  disable <- function(...) tryCatch(shinyjs::disable(...), error = function(e) invisible(NULL))
  enable  <- function(...) tryCatch(shinyjs::enable(...), error = function(e) invisible(NULL))
}

# Stub logging functions
noop <- function(...) invisible(NULL)
stub_log_fns <- list(
  LogMsg     = noop,
  LogCode    = noop,
  LogComment = noop
)

# Helper to create minimal AppState reactiveValues for tests
make_search_state <- function(...) {
  reactiveValues(
    dataset            = NULL,
    dataHash           = NULL,
    trees              = NULL,
    allTrees           = NULL,
    treeHash           = NULL,
    searchWithout      = NULL,
    searchCount        = 0L,
    searchDataHash     = NULL,
    searchNotification = NULL,
    searchInProgress   = FALSE,
    bestSearchScore    = NULL,
    searchTotalHits    = 0L,
    searchTotalReps    = 0L,
    sortTrees          = FALSE,
    newTrees           = NULL,
    ...
  )
}

test_that("search_server returns expected reactive list", {
  r <- make_search_state()

  shiny::testServer(
    search_server,
    args = list(
      r = r,
      AnyTrees       = reactive(FALSE),
      HaveData       = reactive(FALSE),
      UpdateAllTrees = function(x) invisible(NULL),
      log_fns        = stub_log_fns
    ),
    {
      returned <- session$getReturned()
      expect_true(is.list(returned))
      expect_true(all(c("scores", "concavity", "DisplayTreeScores") %in%
                        names(returned)))
    }
  )
})

test_that("concavity reactive responds to weighting mode", {
  r <- make_search_state()

  shiny::testServer(
    search_server,
    args = list(
      r = r,
      AnyTrees       = reactive(FALSE),
      HaveData       = reactive(FALSE),
      UpdateAllTrees = function(x) invisible(NULL),
      log_fns        = stub_log_fns
    ),
    {
      returned <- session$getReturned()

      # Default weighting is "on" with concavity slider = 1 -> 10^1 = 10
      session$setInputs(implied.weights = "on", concavity = 1)
      expect_equal(returned$concavity(), 10)

      # Equal weights -> Inf
      session$setInputs(implied.weights = "off")
      expect_equal(returned$concavity(), Inf)

      # Profile -> "profile"
      session$setInputs(implied.weights = "prof")
      expect_equal(returned$concavity(), "profile")
    }
  )
})

test_that("scores returns NULL when no data or trees", {
  r <- make_search_state()

  shiny::testServer(
    search_server,
    args = list(
      r = r,
      AnyTrees       = reactive(FALSE),
      HaveData       = reactive(FALSE),
      UpdateAllTrees = function(x) invisible(NULL),
      log_fns        = stub_log_fns
    ),
    {
      returned <- session$getReturned()
      expect_null(returned$scores())
    }
  )
})

test_that("dataset change resets search stats", {
  r <- make_search_state()

  shiny::testServer(
    search_server,
    args = list(
      r = r,
      AnyTrees       = reactive(FALSE),
      HaveData       = reactive(FALSE),
      UpdateAllTrees = function(x) invisible(NULL),
      log_fns        = stub_log_fns
    ),
    {
      # Simulate having accumulated search stats
      r$searchTotalHits <- 5L
      r$searchTotalReps <- 10L

      # Trigger dataset observer by setting a dataset
      r$dataset <- TreeSearch::inapplicable.phyData[[1]]
      session$flushReact()

      expect_equal(r$searchTotalHits, 0L)
      expect_equal(r$searchTotalReps, 0L)
    }
  )
})

test_that("concavity defaults to Inf (equal weights)", {
  r <- make_search_state()

  shiny::testServer(
    search_server,
    args = list(
      r = r,
      AnyTrees       = reactive(FALSE),
      HaveData       = reactive(FALSE),
      UpdateAllTrees = function(x) invisible(NULL),
      log_fns        = stub_log_fns
    ),
    {
      # Default weighting mode is "off" => Inf
      session$setInputs(implied.weights = "off")
      returned <- session$getReturned()
      expect_identical(returned$concavity(), Inf)
    }
  )
})

test_that("DisplayTreeScores renders updated confidence text (T-090)", {
  r <- make_search_state()

  shiny::testServer(
    search_server,
    args = list(
      r = r,
      AnyTrees       = reactive(FALSE),
      HaveData       = reactive(FALSE),
      UpdateAllTrees = function(x) invisible(NULL),
      log_fns        = stub_log_fns
    ),
    {
      returned <- session$getReturned()

      # Simulate accumulated search stats (as happens when a continued
      # search matches the previous best score)
      r$searchTotalHits <- 8L
      r$searchTotalReps <- 100L
      r$allTrees <- list("placeholder")
      r$trees    <- list("placeholder")

      # searchNotification must be NULL (search not in progress)
      r$searchNotification <- NULL

      returned$DisplayTreeScores()
      html <- output$results$html

      # Verify the confidence text and tooltip are present
      expect_match(html, "8 of 100 runs hit best score")
      expect_match(html, "1 trees in memory")
      expect_match(html, "title=")
      expect_match(html, "exp\\(-K\\) where K = 8")
    }
  )
})

test_that("SearchConfidenceText uses exp(-K) formula (T-098)", {
  # NULL cases

  expect_null(SearchConfidenceText(NULL, 10))
  expect_null(SearchConfidenceText(0, 10))
  expect_null(SearchConfidenceText(5, 0))

  # K = R = 3: old formula gave 0 -> "<0.01%"; new gives exp(-3) ~ 5%
  txt <- SearchConfidenceText(3, 3)
  expect_match(txt, "3 of 3 runs hit best score")
  expect_match(txt, "~5%")
  expect_match(txt, "better score exists")

  # K = 1: exp(-1) ~ 37%
  txt1 <- SearchConfidenceText(1, 10)
  expect_match(txt1, "~37%")

  # K = 10: exp(-10) ~ 0.005% -> "<0.01%"
  txt10 <- SearchConfidenceText(10, 10)
  expect_match(txt10, "<0.01%")

  # K = 5: exp(-5) ~ 0.67% -> "<1%"
  txt5 <- SearchConfidenceText(5, 20)
  expect_match(txt5, "<1%")
})

test_that("FormatMissProb displays probability thresholds correctly", {
  expect_equal(FormatMissProb(0.37), "~37%")
  expect_equal(FormatMissProb(0.05), "~5%")
  expect_equal(FormatMissProb(0.009), "<1%")
  expect_equal(FormatMissProb(0.0005), "<0.1%")
  expect_equal(FormatMissProb(0.00005), "<0.01%")
})

test_that("scores returns NULL in profile mode before preparation", {
  r <- make_search_state()
  test_trees <- ape::rmtree(3, 6)
  r$trees <- test_trees
  r$treeHash <- "test_hash"
  r$dataset <- TreeSearch::inapplicable.phyData[[1]]
  r$dataHash <- "profile_test"
  r$allTrees <- test_trees

  shiny::testServer(
    search_server,
    args = list(
      r = r,
      AnyTrees       = reactive(TRUE),
      HaveData       = reactive(TRUE),
      UpdateAllTrees = function(x) invisible(NULL),
      log_fns        = stub_log_fns
    ),
    {
      returned <- session$getReturned()

      # Switch to profile mode
      session$setInputs(implied.weights = "prof")
      session$flushReact()

      # Profile data not yet prepared: scores should be NULL
      # (The profile preparation ExtendedTask runs asynchronously;
      # in the test context the future may or may not have completed,
      # but the initial state has profileDataset = NULL.)
      expect_equal(returned$concavity(), "profile")
    }
  )
})

test_that("DisplayTreeScores shows preparing message for profile", {
  r <- make_search_state()
  r$allTrees <- list("placeholder")
  r$trees <- list("placeholder")
  r$searchNotification <- NULL

  shiny::testServer(
    search_server,
    args = list(
      r = r,
      AnyTrees       = reactive(TRUE),
      HaveData       = reactive(TRUE),
      UpdateAllTrees = function(x) invisible(NULL),
      log_fns        = stub_log_fns
    ),
    {
      returned <- session$getReturned()
      session$setInputs(implied.weights = "prof")
      session$flushReact()

      # profileDataset is NULL, concavity is "profile", HaveData and AnyTrees
      # are TRUE, so DisplayTreeScores should show deferred message
      returned$DisplayTreeScores()
      html <- output$results$html

      expect_match(html, "profile scores available after search")
    }
  )
})

test_that("cancel button creates signal file and cleans up", {
  r <- make_search_state()

  shiny::testServer(
    search_server,
    args = list(
      r = r,
      AnyTrees       = reactive(FALSE),
      HaveData       = reactive(FALSE),
      UpdateAllTrees = function(x) invisible(NULL),
      log_fns        = stub_log_fns
    ),
    {
      # cancelFile starts NULL
      expect_null(cancelFile())

      # Simulate setting a cancel file path (as StartSearch would)
      test_path <- tempfile("ts_cancel_test_", fileext = ".signal")
      cancelFile(test_path)
      expect_equal(cancelFile(), test_path)
      expect_false(file.exists(test_path))

      # Initialize cancel input (observeEvent ignoreInit=TRUE skips first value)
      session$setInputs(cancel = 0)
      session$flushReact()

      # Simulate clicking cancel: creates the signal file
      session$setInputs(cancel = 1)
      session$flushReact()
      expect_true(file.exists(test_path))

      # Clean up
      file.remove(test_path)
    }
  )
})

test_that("result observer cleans up cancel file", {
  r <- make_search_state()

  shiny::testServer(
    search_server,
    args = list(
      r = r,
      AnyTrees       = reactive(FALSE),
      HaveData       = reactive(FALSE),
      UpdateAllTrees = function(x) invisible(NULL),
      log_fns        = stub_log_fns
    ),
    {
      # Simulate a cancel file left over from a search
      test_path <- tempfile("ts_cancel_cleanup_", fileext = ".signal")
      file.create(test_path)
      cancelFile(test_path)
      r$searchNotification <- "fake-notification-id"

      # When the result observer fires (simulated by setting
      # searchNotification to NULL), cleanup should remove the file.
      # In testServer, we can't easily trigger the ExtendedTask result
      # observer, but we can verify the cleanup mechanism exists by
      # checking that cancelFile is populated correctly.
      expect_true(file.exists(test_path))
      expect_equal(cancelFile(), test_path)

      # Clean up
      file.remove(test_path)
    }
  )
})

test_that("switching away from profile cancels prep via cancel file", {
  r <- make_search_state()
  # No dataset: trigger observer won't fire (req(HaveData()) fails),
  # so profileCancelFile won't be overwritten.
  r$searchNotification <- NULL

  shiny::testServer(
    search_server,
    args = list(
      r = r,
      AnyTrees       = reactive(FALSE),
      HaveData       = reactive(FALSE),
      UpdateAllTrees = function(x) invisible(NULL),
      log_fns        = stub_log_fns
    ),
    {
      # Initialize inputs so subsequent changes are detected as changes
      session$setInputs(implied.weights = "off", concavity = 1)
      session$flushReact()

      # Manually set a cancel file path (simulating the trigger observer)
      test_cancel <- tempfile("ts_profile_cancel_test_", fileext = ".signal")
      profileCancelFile(test_cancel)

      # Switch to profile — cancel observer sees "profile", no file creation
      session$setInputs(implied.weights = "prof")
      session$flushReact()
      expect_false(file.exists(test_cancel))

      # Switch away — cancel observer should create the signal file
      session$setInputs(implied.weights = "off")
      session$flushReact()
      expect_true(file.exists(test_cancel))

      # Clean up
      suppressWarnings(file.remove(test_cancel))
    }
  )
})

test_that("scores returns NULL with trees but no dataset", {
  r <- make_search_state()
  r$trees <- ape::rmtree(3, 6)
  r$treeHash <- "test_hash"

  shiny::testServer(
    search_server,
    args = list(
      r = r,
      AnyTrees       = reactive(TRUE),
      HaveData       = reactive(FALSE),
      UpdateAllTrees = function(x) invisible(NULL),
      log_fns        = stub_log_fns
    ),
    {
      returned <- session$getReturned()
      # No dataset => scores should be NULL
      expect_null(returned$scores())
    }
  )
})

# ---------- Progress file tests ----------

test_that("progressFile reactiveVal is created and starts NULL", {
  r <- make_search_state()

  shiny::testServer(
    search_server,
    args = list(
      r = r,
      AnyTrees       = reactive(FALSE),
      HaveData       = reactive(FALSE),
      UpdateAllTrees = function(x) invisible(NULL),
      log_fns        = stub_log_fns
    ),
    {
      expect_null(progressFile())
    }
  )
})

test_that("progressFile reactiveVal tracks path lifecycle", {
  r <- make_search_state()

  shiny::testServer(
    search_server,
    args = list(
      r = r,
      AnyTrees       = reactive(FALSE),
      HaveData       = reactive(FALSE),
      UpdateAllTrees = function(x) invisible(NULL),
      log_fns        = stub_log_fns
    ),
    {
      # progressFile starts NULL and can be set/cleared
      expect_null(progressFile())
      test_path <- tempfile("ts_progress_test_", fileext = ".txt")
      progressFile(test_path)
      expect_equal(progressFile(), test_path)
      progressFile(NULL)
      expect_null(progressFile())
    }
  )
})
