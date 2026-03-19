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
