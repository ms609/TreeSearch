library(shiny)
library(shinyjs)  # mod_clustering's clThresh observer calls runjs()

# Source the module under test (relative to tests/testthat/)
source("../../server/mod_clustering.R")

# Stub globals that the module references from global.R
palettes <- list(
  "#cc9966",
  c("#cc9966", "#336699"),
  c("#cc9966", "#336699", "#669933")
)
Notification <- function(...) invisible(NULL)
EnC <- function(x) paste0("c(", paste(paste0("\"", x, "\""), collapse = ", "), ")")

# Stub logging functions
noop <- function(...) invisible(NULL)
stub_log_fns <- list(
  LogMsg      = noop,
  LogCommentP = noop,
  LogCodeP    = noop,
  LogIndent   = noop,
  BeginLogP   = noop,
  LogExprP    = noop
)

test_that("clustering_server returns expected reactive list", {
  r <- reactiveValues(
    trees    = ape::rmtree(5, 10),
    treeHash = "hash1"
  )

  shiny::testServer(
    clustering_server,
    args = list(
      r        = r,
      distMeth = reactive("rf"),
      log_fns  = stub_log_fns
    ),
    {
      returned <- session$getReturned()
      expect_true(is.list(returned))
      expect_true(all(c("distances", "LogDistances", "silThreshold",
                         "clusterings", "LogClusterings") %in%
                        names(returned)))
    }
  )
})

test_that("silThreshold tracks clThresh input", {
  r <- reactiveValues(
    trees    = ape::rmtree(3, 6),
    treeHash = "hash2"
  )

  shiny::testServer(
    clustering_server,
    args = list(
      r        = r,
      distMeth = reactive("rf"),
      log_fns  = stub_log_fns
    ),
    {
      session$setInputs(clThresh = 0.7)
      returned <- session$getReturned()
      expect_equal(returned$silThreshold(), 0.7)
    }
  )
})

test_that("distances returns matrix for single tree", {
  r <- reactiveValues(
    trees    = ape::rmtree(1, 6),
    treeHash = "hash3"
  )

  shiny::testServer(
    clustering_server,
    args = list(
      r        = r,
      distMeth = reactive("rf"),
      log_fns  = stub_log_fns
    ),
    {
      returned <- session$getReturned()
      d <- returned$distances()
      expect_equal(dim(d), c(0, 0))
    }
  )
})

test_that("clusterings returns 'none' for fewer than 3 trees", {
  r <- reactiveValues(
    trees    = ape::rmtree(2, 6),
    treeHash = "hash4"
  )

  shiny::testServer(
    clustering_server,
    args = list(
      r        = r,
      distMeth = reactive("rf"),
      log_fns  = stub_log_fns
    ),
    {
      session$setInputs(clThresh = 0.5)
      returned <- session$getReturned()
      cl <- returned$clusterings()
      expect_equal(cl$method, "no significant clustering")
      expect_equal(cl$n, 1)
    }
  )
})

test_that("distances returns valid dist object for multiple trees", {
  set.seed(4821)
  r <- reactiveValues(
    trees    = ape::rmtree(5, 8),
    treeHash = "hash5"
  )

  shiny::testServer(
    clustering_server,
    args = list(
      r        = r,
      distMeth = reactive("rf"),
      log_fns  = stub_log_fns
    ),
    {
      returned <- session$getReturned()
      d <- returned$distances()
      dm <- as.matrix(d)
      expect_equal(nrow(dm), 5)
      expect_equal(ncol(dm), 5)
      # Diagonal should be zero
      expect_equal(unname(diag(dm)), rep(0, 5))
    }
  )
})

test_that("clusterings with many trees returns valid structure", {
  set.seed(7293)
  r <- reactiveValues(
    trees    = ape::rmtree(10, 8),
    treeHash = "hash6"
  )

  shiny::testServer(
    clustering_server,
    args = list(
      r        = r,
      distMeth = reactive("rf"),
      log_fns  = stub_log_fns
    ),
    {
      session$setInputs(clThresh = 0.5)
      returned <- session$getReturned()
      cl <- returned$clusterings()
      expect_true(cl$n >= 1)
      expect_true(is.numeric(cl$sil))
      # Cluster vector length matches tree count when clustering is found;
      # may be scalar 1 when no significant clustering
      expect_true(length(cl$cluster) %in% c(1, 10))
    }
  )
})
