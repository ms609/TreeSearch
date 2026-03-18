  ##############################################################################
  # Load data
  ##############################################################################
  
  tipLabels <- reactive({r$trees[[1]][["tip.label"]]})
  
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
  
  UpdateData <- reactive({
    source <- input$dataSource
    if (source == "file") {
      if (!r$dataFileVisible) {
        showElement("dataFile")
        r$dataFileVisible <- TRUE
        runjs("console.log($('#dataFile-label'))")
        runjs(paste0(
          "$('#dataFile-label').parent()",
          ".css({'outline': 'dashed #428bca 20px', ", 
          "'width': '100%'})", 
          ".animate({'outline-width': '0px'}, 'slow');"))
        return()
      }
      
      fileInput <- input$dataFile
      r$dataset <- NULL
      r$chars <- NULL
      if (is.null(fileInput)) {
        # How can this be?
        Notification(type = "error", "No data file selected")
        return("No data file selected.")
      }
      dataFile <- fileInput$datapath
      if (is.null(dataFile)) {
        Notification(type = "error", "No data file found.")
        return ("No data file specified.")
      }
      
      LogMsg("UpdateData(): from file")
      r$sortTrees <- FALSE # Trees loaded from dataset may be in sequence
      r$readDataFile <- NULL
      r$bestSearchScore <- NULL  # Reset accumulated score on new dataset
      
      if (length(grep("\\.xlsx?$", dataFile))) {
        if (!requireNamespace("readxl", quietly = TRUE)) {
          install.packages("readxl")
        }
        showElement("readxl.options", anim = TRUE)
        
        r$dataset <- tryCatch({
          sheets <- readxl::excel_sheets(dataFile)
          updateSelectInput(session,
                            inputId = "readxl.sheet",
                            choices = setNames(sheets, sheets),
                            selected = if (input$readxl.sheet %in% sheets) {
                                input$readxl.sheet
                              } else {
                                sheets[1]
                              })

          tibble <- readxl::read_excel(
            path = dataFile,
            sheet = match(input$readxl.sheet, sheets, nomatch = 1L),
            skip = max(0L, input$readxlSkip - 2L),
            .name_repair = "minimal",
            col_types = "text"
          )
          
          firstCol <- input$readxlSkipCols - 1L
          chars <- colnames(tibble)[-seq_len(firstCol)]
          taxNames <- gsub(" ", "_", trimws(unlist(tibble[, firstCol])))
          output$readxl.taxa <- renderUI(HTML(paste(
            "<em>Taxon names</em>:",
            paste(head(taxNames, 3), collapse = ", "),
            "...\n")))
          output$readxl.chars <- renderUI(HTML(paste(
            "<em>Character names</em>:",
            # not r$chars, which may be modified before output updated
            paste(head(chars, 3), collapse = ", "),
            "..."
          )))
          r$chars <- chars
          
          dat <- as.matrix(tibble[, -seq_len(firstCol)])
          rownames(dat) <- taxNames
          dat <- MatrixToPhyDat(dat)
          if (attr(dat, "nr") == 0) {
            stop("No characters loaded; throw error")
          }
          
          # Lines that could cause an error must come before log
          
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
            paste0("  sheet = ", match(input$readxl.sheet, sheets, 1L), ","),
            paste0("  skip = ", max(0L, input$readxlSkip - 2L), ","),
            "  .name_repair = \"minimal\",",
            "  col_types = \"text\"",
            ")",
            paste0("dat <- as.matrix(excelSheet[, -seq_len(", firstCol, ")])"),
            paste0("rownames(dat) <- unlist(excelSheet[, ", firstCol, "])"),
            "dataset <- MatrixToPhyDat(dat)"
          ))
          
          # Return:
          dat
        }, error = function(e) {
          NULL
        })
      } else {
        hideElement("readxl.options")
      }
      
      if (is.null(r$dataset)) suppressWarnings({
        r$dataset <- tryCatch({
          r$readDataFile <- "ReadTntAsPhyDat(dataFile)"
          
          # Return:
          ReadTntAsPhyDat(dataFile)
        }, error = function(e) tryCatch({
          r$chars <- tryCatch(
            ReadCharacters(dataFile),
            error = function(e) {
              Notification(type = "error", "Error reading characters from file")
              # Return:
              NULL
            })
          
          r$charNotes <- tryCatch(
            ReadNotes(dataFile),
            error = function(e) {
              Notification(type = "error", "Error reading character notes")
              # Return:
              NULL
            })
          
          r$readDataFile <- "ReadAsPhyDat(dataFile)"
          
          # Return:
          ReadAsPhyDat(dataFile)
        }, error = function(e) {
          r$readDataFile <- NULL
          # Return:
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
      
      r$sortTrees <- TRUE # Nicer plots
      r$bestSearchScore <- NULL  # Reset accumulated score on new dataset
      
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
      
      updateNumericInput(session, "plottedChar", min = 0L,
                         max = 0L, value = 0L)
      updateSelectizeInput(session, "searchChar", choices = NULL)
      return ("Could not read data from file")
    } else {
      Notification(type = "message", 
                       paste("Loaded", nChars(), "characters and",
                             length(r$dataset), "taxa"))
      
      updateNumericInput(session, "plottedChar", min = 0L,
                         max = nChars(), value = 1L)
      updateSelectizeInput(session, "searchChar",
                           choices = paste0(seq_len(nChars()), ": ", 
                                            colnames(r$chars)),
                           selected = "",
                           server = TRUE)
    }
    
    tryCatch({
      dataFileTrees <- read.nexus(dataFile)
      LogComment("Read trees from dataset file")
      LogCode("newTrees <- read.nexus(dataFile)")
      UpdateAllTrees(dataFileTrees)
      CacheInput("tree", dataFile)
      r$readTreeFile <- "read.nexus(treeFile)"
    }, error = function (e) NULL)
    if (!AnyTrees() || !DatasetMatchesTrees()) {
      updateActionButton(session, "go", "New search")
    } else {
      show("displayConfig")
    }
    
    DisplayTreeScores()
  })
  
