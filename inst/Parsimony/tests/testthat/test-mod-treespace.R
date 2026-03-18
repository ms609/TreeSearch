library(shiny)

# Source the module under test (relative to tests/testthat/)
source("../../server/mod_treespace.R")

# Stub globals that the module references from global.R
palettes <- list(
  "#cc9966",
  c("#cc9966", "#336699"),
  c("#cc9966", "#336699", "#669933")
)
badToGood <- hcl.colors(108, "Temps")
Notification <- function(...) invisible(NULL)
ErrorPlot <- function(msg) {
  plot.new()
  text(0.5, 0.5, msg)
}
Enquote <- function(x) paste0("\"", x, "\"")
EnC <- function(x) paste0("c(", paste(Enquote(x), collapse = ", "), ")")
LogMsg <- function(...) invisible(NULL)

# Stub logging functions
noop <- function(...) invisible(NULL)
stub_log_fns <- list(
  BeginLogP      = noop,
  LogCommentP    = noop,
  LogCodeP       = noop,
  LogIndent      = noop,
  LogClusterings = noop
)

# Minimal clustering stub
stub_clustering <- list(sil = -1, n = 1, cluster = rep(1, 3), method = "none")

# Stub distances (moved to clustering module; treespace now receives as arg)
stub_distances <- reactive({
  matrix(0, 0, 0)
})
stub_LogDistances <- function() invisible(NULL)

test_that("treespace_server returns expected reactive list", {
  r <- reactiveValues(
    trees    = ape::rmtree(5, 10),
    treeHash = "hash1"
  )

  shiny::testServer(
    treespace_server,
    args = list(
      r            = r,
      clusterings  = reactive(stub_clustering),
      silThreshold = reactive(0.5),
      scores       = reactive(rep(10, 5)),
      concavity    = reactive(Inf),
      distMeth     = reactive("cid"),
      plotFormat    = reactive("cons"),
      distances    = stub_distances,
      LogDistances = stub_LogDistances,
      log_fns      = stub_log_fns
    ),
    {
      # Module should return a named list
      returned <- session$getReturned()
      expect_true(is.list(returned))
      expect_true(all(c("mapping", "dims", "nProjDim",
                         "TreeCols", "treePch", "saveDetails",
                         "TreespacePlot", "LogTreespacePlot") %in%
                        names(returned)))
    }
  )
})

test_that("saveDetails returns correct metadata per plot format", {
  r <- reactiveValues(
    trees    = ape::rmtree(3, 6),
    treeHash = "hash2"
  )

  shiny::testServer(
    treespace_server,
    args = list(
      r            = r,
      clusterings  = reactive(stub_clustering),
      silThreshold = reactive(0.5),
      scores       = reactive(NULL),
      concavity    = reactive(Inf),
      distMeth     = reactive("rf"),
      plotFormat    = reactive("space"),
      distances    = stub_distances,
      LogDistances = stub_LogDistances,
      log_fns      = stub_log_fns
    ),
    {
      returned <- session$getReturned()
      sd <- returned$saveDetails()
      expect_equal(sd$fileName, "TreeSpace")
      expect_equal(sd$asp, 1L)
    }
  )
})
