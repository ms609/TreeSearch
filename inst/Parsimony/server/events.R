  ##############################################################################
  # Event listeners
  ##############################################################################
  
  # Data input observers + r$dataset + treeFile now in mod_data.R
  
  observeEvent(input$plotSize, {
    px <- paste0("'", input$plotSize, "px'")
    runjs(paste0("$('#treePlot').css({height: ", px, ", width: ", px, "});"))
  })
  
  # searchConfig modal now handled by mod_search.R
  # treeFile handler now in mod_data.R
  
  # weighting, wtType, scores, DisplayTreeScores, concavity handler
  # now handled by mod_search.R
  
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
  
  # implied.weights show/hide now handled by mod_search.R
  
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
  
