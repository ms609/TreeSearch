library(shiny)

source("../../server/mod_downloads.R")

# ---------------------------------------------------------------------------
# Minimal stubs
# ---------------------------------------------------------------------------
make_args <- function() {
  stub_state <- reactiveValues(
    allTrees    = NULL,
    trees       = NULL,
    plottedTree = NULL,
    dataFiles   = 0L,
    excelFiles  = 0L,
    treeFiles   = 0L
  )
  list(
    state         = stub_state,
    dataSource    = reactive("Agnarsson2004"),
    plotSize      = reactive(600L),
    cmdLogFile    = tempfile(fileext = ".R"),
    stashTrees    = function(trees) invisible(NULL),
    dataFileName  = function(n) paste0("dataFile-", formatC(n, width = 2, flag = "0"), ".txt"),
    excelFileName = function(n) paste0("excelFile-", formatC(n, width = 2, flag = "0"), ".xlsx"),
    treeFileName  = function(n) paste0("treeFile-", formatC(n, width = 2, flag = "0"), ".txt"),
    lastFile      = function(type) NULL,
    mainPlot      = function() invisible(NULL),
    # Nine lines so the testmode sub() for lines 5:9 doesn't go out of bounds
    rCode         = reactive(c("# placeholder", "", "# line3", "# line4",
                               "# - pkg ver1",  "# - pkg ver2", "# - pkg ver3",
                               "# - pkg ver4",  "# - pkg ver5")),
    saveDetails   = reactive(list(fileName = "TestFile", title = "Test", asp = 1))
  )
}

# ---------------------------------------------------------------------------
# Tests
# testServer() executes both the filename and content functions when
# output$xxx is accessed, returning the path of the written temp file.
# ---------------------------------------------------------------------------

test_that("saveZip testmode: copies cmdLogFile and filename is correct", {
  args <- make_args()
  writeLines("# sentinel log", con = args$cmdLogFile)

  with_options(list(shiny.testmode = TRUE), {
    testServer(downloads_server, args = args, {
      path <- output$saveZip
      expect_true(file.exists(path))
      expect_equal(basename(path), "TreeSearch-session.zip")
      expect_equal(readLines(path)[1], "# sentinel log")
    })
  })
})

test_that("savePlotZip testmode: filename uses saveDetails and content is written", {
  args <- make_args()

  with_options(list(shiny.testmode = TRUE), {
    testServer(downloads_server, args = args, {
      path <- output$savePlotZip
      expect_true(file.exists(path))
      expect_equal(basename(path), "TestFile.zip")
      # In testmode the content function calls writeLines() on the file
      content <- readLines(path)
      expect_true(length(content) > 0)
    })
  })
})

test_that("savePlotNwk and savePlotNex use state$plottedTree", {
  library(ape)
  args <- make_args()
  args$state$plottedTree <- rtree(5)

  testServer(downloads_server, args = args, {
    nwk_path <- output$savePlotNwk
    expect_true(file.exists(nwk_path))
    expect_equal(basename(nwk_path), "TreeSearch-consensus.nwk")
    recovered <- read.tree(nwk_path)
    expect_equal(Ntip(recovered), 5L)
  })
})

test_that("saveNwk and saveNex use state$trees", {
  library(ape)
  args <- make_args()
  args$state$trees <- c(rtree(4), rtree(4))

  testServer(downloads_server, args = args, {
    nwk_path <- output$saveNwk
    expect_true(file.exists(nwk_path))
    recovered <- read.tree(nwk_path)
    expect_equal(length(recovered), 2L)
  })
})
