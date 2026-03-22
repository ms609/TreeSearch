# Module: Data loading and tree management
#
# Absorbs data.R + trees.R + data/tree event bindings from events.R.
# Owns inputs: dataSource, dataFile, readxl.sheet, readxlSkip, readxlSkipCols,
#   treeFile, nTree, treeRange.
# Writes most data/tree state fields in AppState.
#
# Returns a list of reactives/functions consumed by other modules/source'd files.

data_ui <- function(id) {
  ns <- NS(id)
  list(
    data_source = selectInput(
      ns("dataSource"), "Dataset",
      c("< Load from file below >" = "file",
        "Agnarsson 2004" = "Agnarsson2004",
        "Sun et al. 2018" = "Sun2018",
        "Wills et al. 2012" = "Wills2012",
        if (logging) setNames(names(inapplicable.datasets),
                              names(inapplicable.datasets)))
    ),
    data_file = fileInput(
      ns("dataFile"),
      tags$span(
        tags$i(class = "fas fa-solid fa-table"),
        tags$span("Load data from file")
      ),
      placeholder = "No data file selected"
    ),
    readxl_options = hidden(tags$span(
      id = ns("readxl_options"),
      selectInput(ns("readxl_sheet"), "Excel sheet to read:",
                  "Sheet 1", "Sheet 1"),
      tags$span("First character row & column:"),
      numericInput(ns("readxlSkip"), label = NULL,
                   min = 2L, value = 2L, step = 1L),
      numericInput(ns("readxlSkipCols"), label = NULL,
                   min = 2L, value = 2L, step = 1L),
      htmlOutput(ns("readxl_chars"), style = "clear: both;"),
      htmlOutput(ns("readxl_taxa"), style = "clear: both; margin-bottom: 1em;")
    )),
    tree_file = fileInput(
      ns("treeFile"),
      label = tags$span(
        tags$i(class = "fas fa-solid fa-tree"),
        tags$span("Load trees")
      ),
      placeholder = "No tree file selected"
    ),
    nTree_input = numericInput(ns("nTree"),
                   label = HTML("Sample <i>n</i> trees from range:"),
                   min = 1L, value = 1L, step = 1L),
    treeRange_input = sliderInput(ns("treeRange"), label = "",
                  min = 1L, max = 1L, step = 1L, value = c(1, 1))
  )
}

#' @param id Module namespace id.
#' @param r AppState reactiveValues.
#' @param parent_session The top-level Shiny session (for cross-module
#'   \code{updateXxxInput} calls targeting non-namespaced inputs).
#' @param callbacks Named list of callback functions from events.R / consensus.R
#'   that the module triggers on tree updates:
#'   \code{DisplayTreeScores}, \code{UpdateKeepNTipsRange},
#'   \code{UpdateDroppedTaxaDisplay}, \code{UpdateOutgroupInput},
#'   \code{KeptTips}.
#' @param log_fns Named list of logging functions from logging.R:
#'   \code{LogMsg}, \code{LogComment}, \code{LogCode}, \code{CacheInput},
#'   \code{LastFile}.
data_server <- function(id, r, parent_session, callbacks, log_fns) {
  moduleServer(id, function(input, output, session) {

    # Unpack logging
    LogMsg      <- log_fns$LogMsg
    LogComment  <- log_fns$LogComment
    LogCode     <- log_fns$LogCode
    CacheInput  <- log_fns$CacheInput
    LastFile    <- log_fns$LastFile

    # Unpack callbacks (from events.R / consensus.R — use isolate-safe pattern)
    DisplayTreeScores      <- callbacks$DisplayTreeScores
    UpdateKeepNTipsRange   <- callbacks$UpdateKeepNTipsRange
    UpdateDroppedTaxaDisplay <- callbacks$UpdateDroppedTaxaDisplay
    UpdateOutgroupInput    <- callbacks$UpdateOutgroupInput

    # Cross-module shinyjs helpers (target top-level DOM ids, not namespaced)
    parentShow <- function(id) {
      runjs(paste0("$('#", id, "').removeClass('shinyjs-hide').show()"))
    }
    parentHide <- function(id) {
      runjs(paste0("$('#", id, "').hide()"))
    }

    ############################################################################
    # Helper reactives (from data.R)
    ############################################################################

    AnyTrees <- reactive({
      !is.null(r$trees) && length(r$trees) > 0
    })

    HaveData <- reactive({
      !is.null(r$dataset) && length(r$dataset) > 0 &&
        inherits(r$dataset, "phyDat")
    })

    tipLabels <- reactive({
      if (!length(r$trees)) return(character(0L))
      r$trees[[1]][["tip.label"]]
    })

    nChars <- reactive({
      if (HaveData()) {
        as.integer(length(attr(r$dataset, "index")))
      } else {
        0L
      }
    })

    TaxonOrder <- reactive({
      if (HaveData()) {
        names(r$dataset)
      } else {
        tipLabels()
      }
    })

    DatasetMatchesTrees <- reactive({
      length(intersect(names(r$dataset), tipLabels())) == length(r$dataset)
    })

    ############################################################################
    # Tree management (from trees.R)
    ############################################################################

    UpdateNTree <- function(n) {
      if (is.null(n) || length(n) == 0) return(FALSE)
      if (n > length(r$allTrees)) {
        r$oldNTree <- n
        n <- length(r$allTrees)
      }
      if (r$nTree == n) {
        FALSE
      } else {
        LogMsg("UpdateNTree(", r$nTree, " -> ", n, ")")
        r$nTree <- n
        if (input$nTree != n) {
          updateNumericInput(session, "nTree", value = n)
        }
        TRUE
      }
    }

    UpdateTreeRange <- function(range) {
      if (is.null(range) || length(range) == 0) return(FALSE)
      if (identical(range, r$treeRange)) {
        FALSE
      } else {
        LogMsg("UpdateTreeRange([", paste(r$treeRange, collapse = ", "),
               "] -> [", paste(range, collapse = ", "), "])")
        r$treeRange <- range
        span <- r$treeRange[2] - r$treeRange[1]
        if (r$nTree > span + 1L) {
          UpdateNTree(span + 1L)
        }
        TRUE
      }
    }

    UpdateActiveTrees <- reactive({
      if (r$updatingTrees) {
        LogMsg("   Skipping UpdateActiveTrees()")
        return()
      }
      r$updatingTrees <- TRUE
      on.exit(r$updatingTrees <- FALSE)
      LogMsg("UpdateActiveTrees()")

      nTrees <- length(r$allTrees)
      if (r$nTree == nTrees &&
          r$treeRange[1] == 1L && r$treeRange[2] == nTrees) {
        thinnedTrees <- r$allTrees
        if (!is.null(r$allTrees) && !identical(r$trees, thinnedTrees)) {
          LogCode("trees <- allTrees")
        }
      } else {
        rangedTrees <- r$allTrees[r$treeRange[1]:r$treeRange[2]]
        thinnedTrees <- WideSample(rangedTrees, r$nTree)

        if (!is.null(r$allTrees) && !identical(r$trees, thinnedTrees)) {
          LogCode(paste0(
            "trees <- WideSample(allTrees[",
            r$treeRange[1], ":", r$treeRange[2],
            "], ", r$nTree, ")"))
        }
      }

      r$trees <- thinnedTrees
      r$treeHash <- rlang::hash(r$trees)

      DisplayTreeScores()

      # Consensus module observes r$treeHash for whichTree, keepNTips,
      # neverDrop, outgroup, droppedTips updates (T-063).

      updateSelectizeInput(session = parent_session,
                           inputId = "treespace-relators",
                           choices = tipLabels(),
                           selected = parent_session$input[["treespace-relators"]])
    })

    UpdateAllTrees <- function(newTrees) {
      LogMsg("UpdateAllTrees()")
      on.exit(LogMsg("/UpdateAllTrees()"), add = TRUE)

      newTrees <- c(newTrees)
      if (length(newTrees) > 1L) {
        newTrees <- RenumberTips(newTrees, newTrees[[1]]$tip.label)
      }
      if (identical(newTrees, r$newTrees)) {
        LogMsg("   <Trees unchanged; returning>")
        return()
      }
      r$newTrees <- newTrees

      oldNTrees <- length(r$allTrees)

      if (!identical(r$allTrees, newTrees)) {
        LogCode("allTrees <- newTrees")
        r$allTrees <- newTrees
      }
      nTrees <- length(newTrees)

      if (nTrees != oldNTrees) {
        if (nTrees > 0L) {
          if (!identical(input$treeRange, c(1L, nTrees))) {
            r$oldTreeRange <- input$treeRange
          }
          UpdateTreeRange(c(1L, nTrees))
          updateSliderInput(session, "treeRange",
                            min = 1L, max = nTrees,
                            value = r$treeRange)

          r$oldNTree <- input$nTree
          UpdateNTree(min(max(input$nTree, aFewTrees), nTrees))
          updateNumericInput(session, "nTree", max = nTrees,
                             value = r$nTree)
        }
        # When nTrees == 0, skip slider updates — the tree manipulation panel
        # is hidden by the parentHide("manipulateTreeset") call below, so no
        # visible element needs updating and we avoid min > max warnings.
      }

      UpdateActiveTrees()
      if (AnyTrees()) {
        parentShow("manipulateTreeset")
      } else {
        parentHide("manipulateTreeset")
      }
    }

    # Debounced nTree / treeRange watchers
    FetchNTree <- debounce(reactive({
      if (!is.null(r$oldNTree)) {
        if (!identical(input$nTree, r$oldNTree)) {
          r$oldNTree <- NULL
        }
      } else {
        if (UpdateNTree(input$nTree)) {
          UpdateActiveTrees()
        }
      }
    }), typingJiffy)

    FetchTreeRange <- debounce(reactive({
      if (!is.null(r$oldTreeRange)) {
        if (!identical(input$treeRange, r$oldTreeRange)) {
          r$oldTreeRange <- NULL
        }
      } else {
        if (UpdateTreeRange(input$treeRange)) {
          UpdateActiveTrees()
        }
      }
    }), aJiffy)

    # Force evaluation of the debounced reactives
    observe(FetchNTree())
    observe(FetchTreeRange())

    ############################################################################
    # Data loading (from data.R + events.R bindings)
    ############################################################################

    UpdateData <- reactive({
      source <- input$dataSource
      if (source == "file") {
        if (!r$dataFileVisible) {
          showElement("dataFile")
          r$dataFileVisible <- TRUE
          dfId <- session$ns("dataFile")
          runjs(paste0("console.log($('#", dfId, "-label'));"))
          runjs(paste0(
            "$('#", dfId, "-label').parent()",
            ".css({'outline': 'dashed #428bca 20px', ",
            "'width': '100%'})",
            ".animate({'outline-width': '0px'}, 'slow');"))
          return()
        }

        fileInput <- input$dataFile
        r$dataset <- NULL
        r$chars <- NULL
        if (is.null(fileInput)) {
          Notification(type = "error", "No data file selected")
          return("No data file selected.")
        }
        dataFile <- fileInput$datapath
        if (is.null(dataFile)) {
          Notification(type = "error", "No data file found.")
          return("No data file specified.")
        }

        LogMsg("UpdateData(): from file")
        r$sortTrees <- FALSE
        r$readDataFile <- NULL
        r$bestSearchScore <- NULL

        if (length(grep("\\.xlsx?$", dataFile))) {
          if (!requireNamespace("readxl", quietly = TRUE)) {
            install.packages("readxl")
          }
          showElement("readxl_options")

          r$dataset <- tryCatch({
            sheets <- readxl::excel_sheets(dataFile)
            updateSelectInput(session,
                              inputId = "readxl_sheet",
                              choices = setNames(sheets, sheets),
                              selected = if (input$readxl_sheet %in% sheets) {
                                input$readxl_sheet
                              } else {
                                sheets[1]
                              })

            tibble <- readxl::read_excel(
              path = dataFile,
              sheet = match(input$readxl_sheet, sheets, nomatch = 1L),
              skip = max(0L, input$readxlSkip - 2L),
              .name_repair = "minimal",
              col_types = "text"
            )

            firstCol <- input$readxlSkipCols - 1L
            chars <- colnames(tibble)[-seq_len(firstCol)]
            taxNames <- gsub(" ", "_", trimws(unlist(tibble[, firstCol])))
            output$readxl_taxa <- renderUI(HTML(paste(
              "<em>Taxon names</em>:",
              paste(head(taxNames, 3), collapse = ", "), "...\n")))
            output$readxl_chars <- renderUI(HTML(paste(
              "<em>Character names</em>:",
              paste(head(chars, 3), collapse = ", "), "...")))
            r$chars <- chars

            dat <- as.matrix(tibble[, -seq_len(firstCol)])
            rownames(dat) <- taxNames
            dat <- MatrixToPhyDat(dat)
            if (attr(dat, "nr") == 0) {
              stop("No characters loaded; throw error")
            }

            LogComment("Load data from spreadsheet", 2)
            if (r$excelFiles == 0 ||
                tools::md5sum(dataFile) !=
                tools::md5sum(paste0(tempdir(), "/", LastFile("excel")))) {
              CacheInput("excel", dataFile)
            }
            LogCode(c(
              paste0("dataFile <- \"", LastFile("excel"), "\""),
              "excelSheet <- readxl::read_excel(",
              "  path = dataFile,",
              paste0("  sheet = ",
                     match(input$readxl_sheet, sheets, 1L), ","),
              paste0("  skip = ", max(0L, input$readxlSkip - 2L), ","),
              "  .name_repair = \"minimal\",",
              "  col_types = \"text\"",
              ")",
              paste0("dat <- as.matrix(excelSheet[, -seq_len(",
                     firstCol, ")])"),
              paste0("rownames(dat) <- unlist(excelSheet[, ",
                     firstCol, "])"),
              "dataset <- MatrixToPhyDat(dat)"
            ))

            dat
          }, error = function(e) NULL)
        } else {
          hideElement("readxl_options")
        }

        if (is.null(r$dataset)) suppressWarnings({
          r$dataset <- tryCatch({
            r$readDataFile <- "ReadTntAsPhyDat(dataFile)"
            ReadTntAsPhyDat(dataFile)
          }, error = function(e) tryCatch({
            r$chars <- tryCatch(
              ReadCharacters(dataFile),
              error = function(e) {
                Notification(type = "error",
                             "Error reading characters from file")
                NULL
              })

            r$charNotes <- tryCatch(
              ReadNotes(dataFile),
              error = function(e) {
                Notification(type = "error",
                             "Error reading character notes")
                NULL
              })

            r$readDataFile <- "ReadAsPhyDat(dataFile)"
            ReadAsPhyDat(dataFile)
          }, error = function(e) {
            r$readDataFile <- NULL
            NULL
          }))

          if (!is.null(r$dataset)) {
            LogComment("Load data from file", 2)
            CacheInput("data", dataFile)
            LogCode(c(
              paste0("dataFile <- \"", LastFile("data"), "\""),
              paste0("dataset <- ", r$readDataFile)
            ))
          }
        })
      } else {
        LogMsg("UpdateData(): from package")

        r$sortTrees <- TRUE
        r$bestSearchScore <- NULL

        r$dataFileVisible <- FALSE
        hideElement("dataFile")

        dataFile <- system.file(paste0("datasets/", source, ".nex"),
                                package = "TreeSearch")
        CacheInput("data", dataFile)
        r$chars <- ReadCharacters(dataFile)
        r$charNotes <- ReadNotes(dataFile)
        r$readDataFile <- "ReadAsPhyDat(dataFile)"
        r$dataset <- ReadAsPhyDat(dataFile)
        LogComment("Load dataset file from TreeSearch package")
        LogCode(c(
          paste0("dataFile <- system.file(\"datasets/", source,
                 ".nex\", package = \"TreeSearch\")"),
          "dataset <- ReadAsPhyDat(dataFile)"
        ))
      }

      if (is.null(r$dataset)) {
        Notification(type = "error", "Could not read data from file")
        # Consensus module observes nChars() for plottedChar/searchChar (T-063)
        return("Could not read data from file")
      } else {
        Notification(type = "message",
                     paste("Loaded", nChars(), "characters and",
                           length(r$dataset), "taxa"))
        # Consensus module observes nChars() for plottedChar/searchChar (T-063)
      }

      tryCatch({
        # suppressWarnings: ape::read.nexus emits a spurious recycling warning
        # when a NEXUS file has unequal counts of [ and ] comment brackets on
        # a line (upstream ape bug; does not affect parsing correctness).
        dataFileTrees <- suppressWarnings(read.nexus(dataFile))
        LogComment("Read trees from dataset file")
        LogCode("newTrees <- read.nexus(dataFile)")
        UpdateAllTrees(dataFileTrees)
        CacheInput("tree", dataFile)
        r$readTreeFile <- "read.nexus(treeFile)"
      }, error = function(e) {
        # Data file has no trees — clear stale trees only if they don't
        # match the new dataset (prevents blank plot from incompatible tips).
        # Keep trees if they match (e.g., re-selecting same dataset after search).
        if (AnyTrees() && !DatasetMatchesTrees()) {
          UpdateAllTrees(list())
        }
      })
      if (AnyTrees() && DatasetMatchesTrees()) {
        parentShow("displayConfig")
      }
      # Button labels reactively managed by mod_search.R

      DisplayTreeScores()
    })

    ############################################################################
    # Tree file loading (from events.R)
    ############################################################################

    observeEvent(input$treeFile, {
      tmpFile <- input$treeFile$datapath
      newTrees <- tryCatch({
        r$readTreeFile <- "read.tree(treeFile)"
        LogMsg("Trying read.tree()")
        read.tree(tmpFile)
      },
      error = function(x) tryCatch({
        r$readTreeFile <- "read.nexus(treeFile)"
        LogMsg("Trying read.nexus()")
        suppressWarnings(read.nexus(tmpFile))
      },
      error = function(err) tryCatch({
        if (grepl("NA/NaN argument", err)) {
          LogMsg("Terminating tree block")
          withEnd <- tempfile()
          on.exit(unlink(withEnd))
          # suppressWarnings: readLines emits "incomplete final line" warning
          # for files without trailing newline; benign, does not affect parsing.
          writeLines(c(suppressWarnings(readLines(tmpFile)), "\nEND;"), withEnd)
          read.nexus(withEnd)
        } else {
          stop("Next handler, please")
        }
      },
      error = function(x) tryCatch(
        # withCallingHandlers muffles the benign readLines "incomplete final
        # line" warning from ReadTntTree before it reaches the outer warning
        # handler (which is for genuine TNT tip-label warnings only).
        withCallingHandlers(
          {
            r$readTreeFile <- "ReadTntTree(treeFile)"
            ReadTntTree(tmpFile)
          },
          warning = function(w) {
            if (grepl("incomplete final line", conditionMessage(w),
                      ignore.case = TRUE)) {
              invokeRestart("muffleWarning")
            }
          }
        ),
        warning = function(x) tryCatch({
          Notification(as.character(x), type = "warning")
          tryLabels <- TipLabels(r$dataset)
          if (length(tryLabels) > 2) {
            Notification("Inferring tip labels from dataset",
                         type = "warning")
            r$readTreeFile <-
              "ReadTntTree(treeFile, tipLabels = TipLabels(dataset))"
            ReadTntTree(tmpFile, tipLabels = tryLabels)
          } else {
            NULL
          }
        }, error = function(e) NULL),
        error = function(e) NULL))))

      if (is.null(newTrees)) {
        # No trees found: check whether the file is a data file uploaded to
        # the wrong input.  Mirror the data-loader fallback chain
        # (ReadTntAsPhyDat → ReadAsPhyDat).
        autoData <- tryCatch(
          suppressWarnings(ReadTntAsPhyDat(tmpFile)),
          error = function(e) tryCatch(
            suppressWarnings(ReadAsPhyDat(tmpFile)),
            error = function(e) NULL
          )
        )
        if (!is.null(autoData)) {
          # Treat as a data file: load it as the active dataset.
          # observeEvent(r$dataset) handles tree-clearing + hash update.
          r$dataset         <- autoData
          r$chars           <- tryCatch(suppressWarnings(ReadCharacters(tmpFile)),
                                        error = function(e) NULL)
          r$charNotes       <- tryCatch(suppressWarnings(ReadNotes(tmpFile)),
                                        error = function(e) NULL)
          r$readDataFile    <- "ReadAsPhyDat(dataFile)"
          r$sortTrees       <- FALSE
          r$bestSearchScore <- NULL
          Notification(
            paste0("No trees found \u2014 loaded ",
                   length(autoData), " taxa and ",
                   length(attr(autoData, "index")), " characters as dataset"),
            type = "message"
          )
        } else {
          Notification("Trees not in a recognized format", type = "error")
        }
      } else {
        LogComment("Load tree from file", 2)
        CacheInput("tree", tmpFile)
        LogCode(paste0("treeFile <- \"", LastFile("tree"), "\""))
        LogCode(paste0("newTrees <- ", r$readTreeFile))

        UpdateAllTrees(newTrees)

        removeModal()
        Notification(paste("Loaded", length(r$trees), "trees"),
                     type = "message")
        # Button labels reactively managed by mod_search.R
        parentShow("displayConfig")
      }
    })

    ############################################################################
    # Data event bindings (from events.R)
    ############################################################################

    observeEvent(input$dataSource, UpdateData(), ignoreInit = TRUE)
    observeEvent(input$dataFile, UpdateData(), ignoreInit = TRUE)
    observeEvent(input$readxl_sheet, UpdateData(), ignoreInit = TRUE)
    observeEvent(input$readxlSkip, UpdateData(), ignoreInit = TRUE)
    observeEvent(input$readxlSkipCols, UpdateData(), ignoreInit = TRUE)

    observeEvent(r$dataset, {
      r$dataHash <- rlang::hash(r$dataset)
      # Clear stale trees only when they are incompatible with the new dataset.
      # UpdateData() may call UpdateAllTrees() *before* this observer fires, so
      # trees from the same .nex file are already in r$allTrees and are
      # compatible. Unconditionally clearing them blanks the plot and resets
      # the tree count to 0 for all 31 bundled example datasets (T-151).
      if (!HaveData() || !DatasetMatchesTrees()) {
        r$allTrees <- NULL
        r$trees <- NULL
        r$treeHash <- NULL
        r$newTrees <- NULL
        parentHide("manipulateTreeset")
      }
      # Search stat reset + timeout default handled by mod_search.R
    })

    ############################################################################
    # Return reactives/functions for other modules
    ############################################################################

    list(
      AnyTrees            = AnyTrees,
      HaveData            = HaveData,
      tipLabels           = tipLabels,
      nChars              = nChars,
      TaxonOrder          = TaxonOrder,
      DatasetMatchesTrees = DatasetMatchesTrees,
      UpdateAllTrees      = UpdateAllTrees,
      UpdateActiveTrees   = UpdateActiveTrees,
      dataSource          = reactive(input$dataSource)
    )
  })
}
