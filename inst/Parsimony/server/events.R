  ##############################################################################
  # Event listeners
  ##############################################################################
  
  observeEvent(input$dataSource, UpdateData(), ignoreInit = TRUE)
  observeEvent(input$dataFile, UpdateData(), ignoreInit = TRUE)
  observeEvent(input$readxl.sheet, UpdateData(), ignoreInit = TRUE)
  observeEvent(input$readxlSkip, UpdateData(), ignoreInit = TRUE)
  observeEvent(input$readxlSkipCols, UpdateData(), ignoreInit = TRUE)
  
  observeEvent(r$dataset, {
    r$dataHash <- rlang::hash(r$dataset)
  })
  observeEvent(input$plotSize, {
    px <- paste0("'", input$plotSize, "px'")
    runjs(paste0("$('#treePlot').css({height: ", px, ", width: ", px, "});"))
  })
  
  observeEvent(input$searchConfig, {
    updateSelectInput(session, "implied.weights",
                      selected = input$implied.weights)
    updateSliderInput(session, "concavity", value = input$concavity)
    updateNumericInput(session, "epsilon", value = input$epsilon)
    updateSelectInput(session, "strategy", selected = input$strategy)
    updateSliderInput(session, "maxReplicates", value = input$maxReplicates)
    updateSliderInput(session, "targetHits", value = input$targetHits)
    updateSliderInput(session, "timeout", value = input$timeout)
    showModal(modalDialog(
      easyClose = TRUE,
      fluidPage(column(6,
      tagList(
              selectInput("implied.weights", "Step weighting",
                         list("Implied" = "on", "Profile" = "prof",
                              "Equal" = "off"), "on"),
              sliderInput("concavity", "Step weight concavity constant", min = 0L,
                         max = 3L, pre = "10^", value = 1L),
              numericInput("epsilon", "Keep if suboptimal by \u2264", min = 0,
                          value = 0),
              selectInput("strategy", "Search strategy",
                         list("Auto" = "auto", "Sprint" = "sprint",
                              "Default" = "default", "Thorough" = "thorough"),
                         "auto"),
              sliderInput("timeout", "Maximum run duration", min = 1,
                          max = 600, value = 30, post = "min", step = 1),
      )), column(6,
             tagList(
              sliderInput("maxReplicates", "Maximum replicates", min = 1L,
                          max = 500L, value = 100L, step = 1L),
              sliderInput("targetHits", "Stop after N hits to best score",
                          min = 1L, max = 50L, value = 10L, step = 1L),
              selectizeInput("searchWithout", "Exclude taxa", DatasetTips(),
                             r$searchWithout, multiple = TRUE)
             ))
      ),
      title = "Tree search settings",
      footer = tagList(modalButton("Close", icon = Icon("rectangle-xmark")),
                       actionButton("modalGo", icon = Icon("magnifying-glass"),
                                    if(length(r$trees)) {
                                      "Continue search"
                                    } else {
                                      "Start search"
                                    }))
    ))
    show("go")
  })
  
  observeEvent(input$treeFile, {
    tmpFile <- input$treeFile$datapath
    newTrees <- tryCatch({
        r$readTreeFile <- "read.tree(treeFile)"
        LogMsg("Trying read.tree()")
        read.tree(tmpFile)
      },
      error = function (x) tryCatch({
        r$readTreeFile <- "read.nexus(treeFile)"
          LogMsg("Trying read.nexus()")
          read.nexus(tmpFile)
        },
        error = function (err) tryCatch(
          {
            if (grepl("NA/NaN argument", err)) {
              LogMsg("Terminating tree block")
              # Unterminated tree block, perhaps because a search is ongoing
              withEnd <- tempfile()
              on.exit(unlink(withEnd))
              writeLines(c(readLines(tmpFile), "\nEND;"), withEnd)
              read.nexus(withEnd)
            } else {
              stop("Next handler, please")
            }
          },
          error = function (x) tryCatch({
              r$readTreeFile <- "ReadTntTree(treeFile)"
              ReadTntTree(tmpFile)
            }, warning = function (x) tryCatch({
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
            }, error = NULL
            )
          )
        )
      )
    )
    if (is.null(newTrees)) {
      Notification("Trees not in a recognized format", type = "error")
    } else {
      LogComment("Load tree from file", 2)
      CacheInput("tree", tmpFile)
      LogCode(paste0("treeFile <- \"", LastFile("tree"), "\""))
      LogCode(paste0("newTrees <- ", r$readTreeFile))
      
      UpdateAllTrees(newTrees) # updates r$trees
      
      removeModal()
      Notification(paste("Loaded", length(r$trees), "trees"), type = "message")
      updateActionButton(session, "modalGo", "Continue search")
      updateActionButton(session, "go", "Continue")
      show("displayConfig")
    }
    
  })
  
  observeEvent(input$implied.weights, {
    switch(input$implied.weights,
           "on" = show("concavity"),
           hide("concavity")
    )
    DisplayTreeScores()
  })
  
  weighting <- reactive(
    if (length(input$implied.weights) > 0) {
      input$implied.weights
    } else {
      "on"
    }
  )
  wtType <- reactive(switch(weighting(),
                            "on" = paste0("k = ", signif(concavity(), 3)),
                            "off" = "EW",
                            "prof" = "PP"))
  
  scores <- bindCache(reactive({
    if (!HaveData() || !AnyTrees()) {
      return(NULL)
    }
    PutTree(r$trees)
    PutData(r$dataset)
    LogMsg("scores(): Recalculating scores with k = ", concavity())
    withProgress(tryCatch(
      signif(TreeLength(
        RootTree(r$trees, 1),
        r$dataset,
        concavity = concavity()
      )),
      error = function (x) {
        if (HaveData() && AnyTrees()) {
          cli::cli_alert(x[[2]])
          cli::cli_alert_danger(x[[1]])
          Notification(type = "error",
                       "Could not score all trees with dataset")
        }
        NULL
     }),
     value = 0.85, message = "Scoring trees")
  }), r$treeHash, r$dataHash, concavity())
  
  DisplayTreeScores <- function () {
    LogMsg("DisplayTreeScores()")
    treeScores <- scores()
    score <- if (is.null(treeScores)) {
      "; could not be scored from dataset"
    } else if (length(unique(treeScores)) == 1) {
      paste0(", each with score ", treeScores[1], " (", wtType(), ")")
    } else {
      paste0(" with scores ", min(treeScores), " to ", max(treeScores),
             " (", wtType(), ")")
    }
    
    msg <- paste0(
      length(r$allTrees), " trees in memory: ",
      length(r$trees), " sampled", 
      score
    )
    output$results <- renderText(msg)
    msg
  }
  
  observeEvent(input$concavity, {
    DisplayTreeScores()
  }, ignoreInit = TRUE)
  
  TipsInTree <- reactive({
    if (AnyTrees()) {
      length(r$trees[[1]]$tip.label)
    } else {
      0L
    }
  })
  
  UpdateKeepNTipsRange <- reactive({
    if (AnyTrees() && "consConfig" %in% r$visibleConfigs) {
      nTip <- TipsInTree()
      LogMsg("UpdateKeepNTipsRange(", input$keepNTips, " -> ", nTip, ")")
      r$keepNTips <- nNonRogues()
      if (r$keepNTips != input$keepNTips) {
        r$oldkeepNTips <- input$keepNTips
      }
      updateNumericInput(inputId = "keepNTips",
                         label = paste0("Tips to show (/", nTip, "):"),
                         min = max(3L, length(input$neverDrop)),
                         max = nTip,
                         value = nNonRogues())
    }
  })
  
  UpdateExcludedTipsInput <- reactive({
    if (AnyTrees() && "consConfig" %in% r$visibleConfigs) {
      LogMsg("UpdateExcludedTipsInput()")
      dropList <- dropSeq()[seq_along(DroppedTips())]
      updateSelectInput(inputId = "excludedTip",
                        choices = dropList,
                        selected = if(input$excludedTip %in% DroppedTips())
                          input$excludedTip else dropSeq()[1])
      html("droppedList",
           paste0("<label class=\"control-label\">Dropped tips:</label>", 
                  "<ul>", 
                  paste0("<li style=\"color: ", TipCols()[dropList], "\">",
                         dropList, "</li>", collapse = "\r\n"),
                  "</ul>"))
    }
  })
  
  UpdateDroppedTaxaDisplay <- reactive({
    LogMsg("UpdateDroppedTaxaDisplay()")
    if ("consConfig" %in% r$visibleConfigs) {
      if (length(DroppedTips())) {
        UpdateExcludedTipsInput()
        if ("droppedTips" %in% r$visibleConfigs) {
          show("droppedTips")
        }
        if ("droppedList" %in% r$visibleConfigs) {
          show("droppedList")
        }
      } else {
        hide("droppedTips")
        hide("droppedList")
      }
    }
  })
  
  observeEvent(r$visibleConfigs, {
    UpdateDroppedTaxaDisplay()
  })
  
  UpdateOutgroupInput <- reactive({
    if (AnyTrees() && "treePlotConfig" %in% r$visibleConfigs) {
      LogMsg("UpdateOutgroupInput()")
      r$outgroup <- intersect(r$outgroup, KeptTips())
      if (length(r$outgroup) == 0) {
        r$outgroup <- if (HaveData()) {
          intersect(names(r$dataset), KeptTips())[1]
        } else {
          KeptTips()[1]
        }
      }
      
      if (!identical(sort(r$outgroup), sort(input$outgroup))) {
        r$oldOutgroup <- if (is.null(input$outgroup)) {
          NO_OUTGROUP
        } else {
          input$outgroup
        }
      }
      
      updateSelectizeInput(
        inputId = "outgroup",
        selected = r$outgroup,
        choices = KeptTips()
        )
    }
  })
  
  observeEvent(input$implied.weights, {
    switch(input$implied.weights,
           "on" = show("concavity"),
           hide("concavity")
    )
  })
  
  ShowConfigs <- function (visible = character(0)) {
    allConfigs <- c("whichTree", "charChooser",
                    "consConfig", "clusConfig",
                    "clusLegend", "branchLegend",
                    "spaceConfig", "treePlotConfig",
                    "mapConfig", "savePlottedTrees",
                    "droppedTips", "droppedList")
    r$visibleConfigs <- visible
    lapply(visible, show)
    lapply(setdiff(allConfigs, visible), hide)
  }
  
  observeEvent(input$plotFormat, {
    ShowConfigs(switch(input$plotFormat,
                       "ind" = c("whichTree", "charChooser",
                                 "treePlotConfig"),
                       "cons" = c("consConfig", "droppedTips",
                                  "savePlottedTrees",
                                  "treePlotConfig", "branchLegend"),
                       "clus" = c("clusConfig", "clusLegend",
                                  "savePlottedTrees",
                                  "consConfig", "droppedList",
                                  "treePlotConfig"),
                       "space" = c("clusConfig", "clusLegend",
                                   "spaceConfig", "mapConfig"),
                       ""))
  })
  
  
  output$branchLegend <- renderUI({
    if (!AnyTrees()) {
      return()
    }
    LogMsg("renderUI(branchLegend)")
    on.exit(LogMsg("/renderUI(branchLegend)"))
    kept <- KeptTips()
    dropped <- DroppedTips()
    
    if (length(dropped) &&
        length(input$excludedTip) &&
        nchar(input$excludedTip) &&
        input$excludedTip %in% tipLabels()) {
      consTrees <- lapply(r$trees, DropTip, setdiff(dropped, input$excludedTip))
      plotted <- TreeTools::RoguePlot(
        trees = consTrees,
        tip = input$excludedTip,
        p = consP(),
        plot = FALSE
      )
      tagList(
        tags$span(class = "legendLeft", "1 tree"),
        tags$span(id = "blackToGreen", class = "legendBar", "\ua0"),
        tags$span(class = "legendRight",
                  paste(max(c(plotted$onEdge, plotted$atNode)), "trees")),
      )
    }
  })
  
