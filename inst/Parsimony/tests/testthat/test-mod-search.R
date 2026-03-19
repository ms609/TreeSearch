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
SearchConfidenceText <- function(K, R) {
  if (is.null(K) || is.null(R) || R <= 0L || K <= 0L) return(NULL)
  paste0(K, "/", R, " reps")
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
