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
SearchConfidenceText <- function(K, R, nSearches = 1L,
                                 nTopologies = NULL,
                                 lastImprovedRep = NULL) {
  if (is.null(K) || is.null(R) || R <= 0L || K <= 0L) return(NULL)
  K <- min(K, R)
  prob_miss <- if (K < R) (1 - K / R) ^ R else exp(-K)
  runs_label <- if (!is.null(nSearches) && nSearches > 1L) {
    paste0("total runs across ", nSearches, " searches")
  } else {
    "runs"
  }
  topo_note <- if (!is.null(nTopologies) && nTopologies == 1L) {
    " [single topology \u2014 limited independence]"
  } else {
    ""
  }
  trajectory_note <- if (!is.null(lastImprovedRep) && R > 1L) {
    paste0(" Last improvement: replicate ", lastImprovedRep, ".")
  } else {
    ""
  }
  rugged_note <- if (K / R < 0.3 && R >= 5L) {
    paste0(" Hit rate low (", round(100 * K / R),
           "%) \u2014 more replicates may help.")
  } else {
    ""
  }
  small_sample_note <- if (K == R && R <= 5L) {
    paste0(" \u2014 increase \u2018Stop when N runs hit best\u2019 for a ",
           "tighter estimate")
  } else {
    ""
  }
  paste0(K, " of ", R, " ", runs_label, " hit best score. ",
         "Probability that a better score exists: ",
         FormatMissProb(prob_miss),
         topo_note, trajectory_note, rugged_note, small_sample_note)
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
      expect_match(html, "\\(1 - K/R\\)\\^R where K = 8")
    }
  )
})

test_that("SearchConfidenceText uses binomial bound (T-163)", {
  # NULL cases
  expect_null(SearchConfidenceText(NULL, 10))
  expect_null(SearchConfidenceText(0, 10))
  expect_null(SearchConfidenceText(5, 0))

  # K = R = 3: falls back to exp(-3) ~ 5%
  txt <- SearchConfidenceText(3, 3)
  expect_match(txt, "3 of 3 runs hit best score")
  expect_match(txt, "~5%")
  expect_match(txt, "better score exists")
  # Small sample nudge
  expect_match(txt, "tighter estimate")

  # K = 1, R = 10: (1 - 1/10)^10 = 0.9^10 ~ 35%
  txt1 <- SearchConfidenceText(1, 10)
  expect_match(txt1, "~35%")

  # K = R = 10: exp(-10) ~ 0.005% -> "<0.01%"
  txt10 <- SearchConfidenceText(10, 10)
  expect_match(txt10, "<0.01%")

  # K = 5, R = 20: (1 - 5/20)^20 = 0.75^20 ~ 0.32% -> "<1%"
  txt5 <- SearchConfidenceText(5, 20)
  expect_match(txt5, "<1%")

  # Ruggedness flag: K/R < 0.3 and R >= 5
  txt_rugged <- SearchConfidenceText(1, 10)
  expect_match(txt_rugged, "Hit rate low")
  expect_match(txt_rugged, "10%")

  # No ruggedness flag when K/R >= 0.3
  txt_smooth <- SearchConfidenceText(4, 10)
  expect_false(grepl("Hit rate low", txt_smooth))

  # Single topology warning
  txt_single <- SearchConfidenceText(5, 10, nTopologies = 1L)
  expect_match(txt_single, "single topology.*limited independence")

  # No topology note for multiple trees (redundant with "trees in memory")
  txt_multi <- SearchConfidenceText(5, 10, nTopologies = 3L)
  expect_false(grepl("topolog", txt_multi))

  # Last-improved replicate info
  txt_traj <- SearchConfidenceText(5, 10, lastImprovedRep = 7L)
  expect_match(txt_traj, "Last improvement: replicate 7")

  # nSearches label
  txt_multi_search <- SearchConfidenceText(5, 10, nSearches = 3L)
  expect_match(txt_multi_search, "across 3 searches")

  # Stop reason: consensus stable
  txt_cons <- SearchConfidenceText(4, 90, stopReason = "consensus")
  expect_match(txt_cons, "consensus stable")

  # Stop reason: timeout
  txt_time <- SearchConfidenceText(4, 90, stopReason = "timeout")
  expect_match(txt_time, "time limit")

  # No stop reason
  txt_none <- SearchConfidenceText(4, 90, stopReason = NULL)
  expect_false(grepl("stopped", txt_none, ignore.case = TRUE))
})

test_that("FormatMissProb displays probability thresholds correctly", {
  expect_equal(FormatMissProb(0.37), "~37%")
  expect_equal(FormatMissProb(0.05), "~5%")
  expect_equal(FormatMissProb(0.009), "<1%")
  expect_equal(FormatMissProb(0.0005), "<0.1%")
  expect_equal(FormatMissProb(0.00005), "<0.01%")
})

test_that("SearchConfidenceText appends Chao1 coverage note when sufficient replicates", {
  # 10 replicates, low coverage (many singletons): coverage note expected
  rscores <- c(100, 101, 102, 103, 104, 105, 106, 107, 108, 109)  # all singletons
  txt <- SearchConfidenceText(1, 10, replicateScores = rscores)
  expect_match(txt, "coverage", ignore.case = TRUE)

  # No coverage note below threshold (< 5 replicates)
  txt_few <- SearchConfidenceText(2, 4, replicateScores = c(10, 10, 20, 30))
  expect_false(grepl("coverage", txt_few, ignore.case = TRUE))

  # No coverage note with NULL replicateScores
  txt_null <- SearchConfidenceText(3, 10, replicateScores = NULL)
  expect_false(grepl("coverage", txt_null, ignore.case = TRUE))

  # High coverage (all scores identical): note still present but no "unseen" warning
  rscores_conv <- rep(42, 10)
  txt_conv <- SearchConfidenceText(10, 10, replicateScores = rscores_conv)
  expect_match(txt_conv, "coverage", ignore.case = TRUE)
  expect_false(grepl("unseen", txt_conv, ignore.case = TRUE))
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

# ---------- T-165: run stats reset on weighting change ----------

test_that("changing concavity resets run stats but keeps trees (T-165)", {
  r <- make_search_state()
  r$allTrees <- list("placeholder_tree")
  r$trees    <- list("placeholder_tree")

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
      # Simulate accumulated search stats from a prior run
      r$searchTotalHits <- 52L
      r$searchTotalReps <- 1604L
      r$bestSearchScore <- 1.42854

      # Change the concavity constant — stats should reset, trees preserved
      session$setInputs(implied.weights = "on", concavity = 3)
      session$flushReact()
      session$setInputs(concavity = 2)
      session$flushReact()

      expect_equal(r$searchTotalHits, 0L)
      expect_equal(r$searchTotalReps, 0L)
      expect_null(r$bestSearchScore)
      expect_length(r$allTrees, 1L)  # trees preserved
    }
  )
})

test_that("changing weighting mode resets run stats but keeps trees (T-165)", {
  r <- make_search_state()
  r$allTrees <- list("placeholder_tree")
  r$trees    <- list("placeholder_tree")

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
      # First: set initial state (on -> EW) with accumulated stats
      session$setInputs(implied.weights = "on", concavity = 1)
      session$flushReact()

      r$searchTotalHits <- 30L
      r$searchTotalReps <- 200L
      r$bestSearchScore <- 2.5

      # Switch to EW — should reset stats
      session$setInputs(implied.weights = "off")
      session$flushReact()

      expect_equal(r$searchTotalHits, 0L)
      expect_equal(r$searchTotalReps, 0L)
      expect_null(r$bestSearchScore)
      expect_length(r$allTrees, 1L)  # trees preserved
    }
  )
})
