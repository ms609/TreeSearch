# Module: Downloads
#
# Owns all 8 downloadHandler outputs:
#   saveZip, savePlotZip, savePng, savePdf,
#   savePlotNwk, savePlotNex, saveNwk, saveNex
#
# Reactive args (passed from server.R top-level input):
#   dataSource  reactive(input$dataSource)
#   plotSize    reactive(input$plotSize)
#
# Callback args (functions/values from sourced server files):
#   cmdLogFile    character — path to session R-script log
#   stashTrees    function(trees)  — writes trees to temp file (logging.R)
#   dataFileName  function(n)      — (logging.R)
#   excelFileName function(n)      — (logging.R)
#   treeFileName  function(n)      — (logging.R)
#   lastFile      function(type)   — (logging.R)
#   mainPlot      function()       — renders current plot (consensus.R)
#   rCode         reactive         — plot R-script lines (consensus.R)
#   saveDetails   reactive         — list(fileName, title, asp) (treespace.R)

# ---------------------------------------------------------------------------
# UI helpers — returns a named list so scattered buttons can be placed
# individually in ui.R without duplicating ns() logic.
# ---------------------------------------------------------------------------
downloads_ui <- function(id) {
  ns <- NS(id)
  list(
    save_zip      = downloadButton(ns("saveZip"),     "Save log",  icon = Icon("download")),
    save_nwk      = downloadButton(ns("saveNwk"),     "Newick",    icon = Icon("download")),
    save_nex      = downloadButton(ns("saveNex"),     "Nexus",     icon = Icon("download")),
    save_plot_zip = downloadButton(ns("savePlotZip"), "R script",  icon = Icon("download")),
    save_pdf      = downloadButton(ns("savePdf"),     "PDF",       icon = Icon("download")),
    save_png      = downloadButton(ns("savePng"),     "PNG",       icon = Icon("download")),
    save_plot_nwk = downloadButton(ns("savePlotNwk"), "Newick",    icon = Icon("download")),
    save_plot_nex = downloadButton(ns("savePlotNex"), "Nexus",     icon = Icon("download"))
  )
}

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
downloads_server <- function(id, state, dataSource, plotSize,
                              cmdLogFile, stashTrees,
                              dataFileName, excelFileName, treeFileName, lastFile,
                              mainPlot, rCode, saveDetails) {
  moduleServer(id, function(input, output, session) {

    output$saveZip <- downloadHandler(
      filename = function() "TreeSearch-session.zip",
      content = function(file) {
        if (isTRUE(getOption("shiny.testmode"))) {
          file.copy(cmdLogFile, file)
        } else {
          zipDir <- tempfile("zip-")
          dir.create(zipDir)
          on.exit(unlink(zipDir))
          rFile <- paste0(zipDir, "/TreeSearch-session.R")
          file.copy(cmdLogFile, rFile, overwrite = TRUE)
          zip::zip(file, c(
            rFile,
            if (state$dataFiles)
              paste0(tempdir(), "/", dataFileName(seq_len(state$dataFiles))),
            if (state$excelFiles)
              paste0(tempdir(), "/", excelFileName(seq_len(state$excelFiles))),
            if (state$treeFiles)
              paste0(tempdir(), "/", treeFileName(seq_len(state$treeFiles)))
          ), compression_level = 9, mode = "cherry-pick")
        }
      }
    )

    output$savePlotZip <- downloadHandler(
      filename = function() paste0(saveDetails()$fileName, ".zip"),
      content = function(file) {
        stashTrees(state$allTrees)

        if (isTRUE(getOption("shiny.testmode"))) {
          rCode_val <- rCode()
          rCode_val <- sub("TreeSearch plot log: 2[\\d\\-]{9} [012][\\d:]{7}",
                           "TreeSearch plot log: <DATE-AND-TIME>",
                           rCode_val, perl = TRUE)
          rCode_val[4] <- "# System: <SYS-INFO>"
          rCode_val[5:9] <- sub("^(# \\- \\w+ ).*$", "\\1<VERSION>",
                                rCode_val[5:9], perl = TRUE)
          rCode_val <- sub("dataFile <- .*$",
                           paste0("dataFile <- system.file(\"datasets/",
                                  dataSource(),
                                  ".nex\", package = \"TreeSearch\") # FALSE CODE for TEST MODE"),
                           rCode_val,
                           perl = TRUE)
          rCode_val <- sub("treeFile <- .*$",
                           "treeFile <- dataFile # Test mode",
                           rCode_val,
                           perl = TRUE)
          writeLines(rCode_val, con = file)
        } else {
          tempDir <- tempfile("plot-zip-")
          dir.create(tempDir)
          on.exit(unlink(tempDir))
          rFile <- paste0(tempDir, "/", saveDetails()$fileName, ".R")
          writeLines(rCode(), con = rFile)

          zip::zip(file, c(
            rFile,
            paste0(tempdir(), "/", lastFile("data")),
            paste0(tempdir(), "/", lastFile("excel")),
            paste0(tempdir(), "/", lastFile("tree"))
          ), compression_level = 9, mode = "cherry-pick")
        }
      }
    )

    output$savePng <- downloadHandler(
      filename = function() paste0(saveDetails()$fileName, ".png"),
      content = function(file) {
        png(file, width = plotSize(), height = plotSize())
        mainPlot()
        dev.off()
      }
    )

    output$savePdf <- downloadHandler(
      filename = function() paste0(saveDetails()$fileName, ".pdf"),
      content = function(file) {
        width <- 8
        pdf(
          file,
          title = saveDetails()$title,
          width = width,
          height = saveDetails()$asp * width
        )
        mainPlot()
        dev.off()
      }
    )

    output$savePlotNwk <- downloadHandler(
      filename = "TreeSearch-consensus.nwk",
      content = function(file) {
        write.tree(state$plottedTree, file = file)
      }
    )

    output$savePlotNex <- downloadHandler(
      filename = "TreeSearch-consensus.nex",
      content = function(file) {
        write.nexus(state$plottedTree, file = file)
      }
    )

    output$saveNwk <- downloadHandler(
      filename = "TreeSearch.nwk",
      content = function(file) {
        write.tree(state$trees, file = file, tree.names = TRUE)
      }
    )

    output$saveNex <- downloadHandler(
      filename = "TreeSearch.nex",
      content = function(file) {
        write.nexus(state$trees, file = file)
      }
    )
  })
}
