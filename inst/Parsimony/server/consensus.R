  UserRoot <- function(tree) {
    outgroupTips <- intersect(r$outgroup, tree$tip.label)
    if (length(outgroupTips)) {
      # DELETE? tr <- deparse(substitute(tree))
      RootTree(tree, outgroupTips)
    } else {
      tree
    }
  }
  
  LogUserRoot <- function(tree = "cons", dropped = character(0)) {
    outgroupTips <- setdiff(r$outgroup, dropped)
    if (length(outgroupTips)) {
      LogCommentP("Root tree")
      LogCodeP(paste0(tree, " <- RootTree(", tree, ", ", EnC(outgroupTips), ")"))
    }
  }
  
  PlottedChar <- debounce(reactive({
    typed <- max(0L, as.integer(input$plottedChar), na.rm = TRUE)
    if (nChars() > 0 && typed > nChars()) {
      Notification(type = "warning",
                   paste("Dataset contains", nChars(), "characters.")
      )
      updateNumericInput(session, "plottedChar", value = nChars())
    }
    min(typed, nChars())
  }), aJiffy)
  
  observeEvent(PlottedChar(), {
    if (PlottedChar() > 0) {
      showElement("mapDisplay")
    } else {
      hideElement("mapDisplay")
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$searchChar, {
    searchResult <- as.numeric(strsplit(input$searchChar, ": ")[[1]][1])
    if (!is.na(searchResult)) {
      updateNumericInput(session, "plottedChar", value = searchResult)
    }
  })
  
  whichTree <- debounce(reactive(input$whichTree), aJiffy)
  
  PlottedTree <- reactive({
    if (length(r$trees) > 0L) {
      plottedTree <- if (whichTree() > 0) {
        r$trees[[whichTree()]]
      } else {
        Consensus(r$trees, p = 1)
      }
      plottedTree <- UserRoot(plottedTree)
      plottedTree <- SortEdges(plottedTree)
      if (!("tipsRight" %in% input$mapDisplay)) {
        plottedTree$edge.length <- rep_len(2, dim(plottedTree[["edge"]])[[1]])
      }
      plottedTree
    }
  })
  LogPlottedTree <- function() {
    if (whichTree() > 0) {
      LogCodeP(paste0("plottedTree <- trees[[", whichTree(), "]]"))
    } else {
      LogCodeP("plottedTree <- Consensus(trees, p = 1)")
    }
    LogUserRoot("plottedTree")
    if (!("tipsRight" %in% input$mapDisplay)) {
      LogCommentP("Set uniform edge length", 0)
      LogCodeP(
        "plottedTree$edge.length <- rep.int(2, nrow(plottedTree$edge))"
      )
    }
    LogSortEdges("plottedTree")
  }
  
  Instab <- reactive({
    TipInstability(r$trees)
  })
  
  dropSeq <- reactive({
    LogMsg("dropSeq()")
    Rogues()$taxon[-1]
  })
  
  stableCol <- reactive({
    Rogue::ColByStability(r$trees)
  })
  
  Rogues <- bindCache(reactive({
    if (AnyTrees() && inherits(r$trees, "multiPhylo")) {
      LogComment("Check for rogue taxa", 2)
      LogComment(paste0(
        "Use RogueTaxa() in place of QuickRogue() for a more complete ",
        "analysis"))
      LogCode(c(
        "rogues <- Rogue::QuickRogue(",
        "  trees,",
        if (length(input$neverDrop)) paste0(
          "  neverDrop = ", EnC(input$neverDrop), ","
        ),
        "  fullSeq = TRUE,",
        paste0("  p = ", Enquote(consP())),
        ")",
        "print(rogues) # Detailed results of rogue analysis",
        "print(rogues$taxon[-1]) # Sequence of taxa to drop"
      ))
      withProgress(
        message = "Identifying rogues", value = 0.99,
        rogues <- Rogue::QuickRogue(r$trees, neverDrop = input$neverDrop,
                          fullSeq = TRUE, p = consP())
      )
      # TODO delete once Rogue 2.1.2 released -- return QuickRogue above.
      rogues[!rogues$taxon %in% input$neverDrop, ]
    } else {
      data.frame(num = 0, taxNum = NA_integer_, taxon = NA_character_,
                 rawImprovement = NA_real_, IC = 0)
    }
  }), r$treeHash, input$neverDrop, consP())
  
  unitEdge <- reactive({
    TRUE
  })
  
  SortEdges <- function (tr, force = FALSE) {
    if (force || r$sortTrees) {
      # Return:
      SortTree(tr, order = TaxonOrder())
    } else {
      # Return:
      tr
    }
  }
  LogSortEdges <- function(tr) (
    if (r$sortTrees) {
      LogCommentP("Rotate nodes, to display clades in order of size", 0)
      LogCodeP(paste0(
        tr, " <- SortTree(", tr, ", order = ", 
        if (HaveData()) {
          "names(dataset)"
        } else {
          "trees[[1]]$tip.label"
        },
        ")"
      ))
    }
  )
  
  nNonRogues <- reactive({
    LogMsg("nNonRogues()")
    on.exit(LogMsg("nNonRogues: ", nrow(Rogues()) - which.max(Rogues()$IC)))
    nrow(Rogues()) - which.max(Rogues()$IC)
  })
  
  TipCols <- reactive(stableCol()) # TODO allow user to choose how to colour
  
  TipColLegend <- function() {
    PlotTools::SpectrumLegend(
      "bottomleft", horiz = TRUE, inset = 0.01, bty = "n", xpd = NA,
      palette = hcl.colors(131, "inferno")[1:101],
      legend = c("Stable", "Unstable"),
      title = "Leaf stability",
      title.font = 2
    )
  }
  
  consP <- debounce(reactive(signif(input$consP)), 50)
  observeEvent(consP(), {
    if (AnyTrees()) {
      LogMsg("Observed consP()")
      UpdateKeepNTipsRange()
      UpdateDroppedTaxaDisplay()
      r$concordance <- list()
    }
  }, ignoreInit = TRUE)
  
  concordance <- bindCache(reactive({
    LogMsg("concordance()")
    # Return:
    switch(input$concordance,
          "p" = SplitFrequency(r$plottedTree, r$trees) / length(r$trees),
          "qc" = QuartetConcordance(r$plottedTree, r$dataset),
          "mcc" = MutualClusteringConcordance(r$plottedTree, r$dataset),
          "spc" = SharedPhylogeneticConcordance(r$plottedTree, r$dataset),
          "clc" = ClusteringConcordance(r$plottedTree, r$dataset),
          "phc" = PhylogeneticConcordance(r$plottedTree, r$dataset),
          NULL
    )
  }), r$plottedTree, r$treeHash, r$dataHash, input$concordance)
  
  LabelConcordance <- \() {
    LogMsg("LabelConcordance()")
    if (input$concordance != "none" &&
        !is.null(r$plottedTree)) {
      LabelSplits(r$plottedTree, signif(concordance(), 3),
                  col = SupportColor(concordance()),
                  frame = "none", pos = 3L)
    }
  }
  
  LogConcordance <- function(plottedTree = "plottedTree") {
    if (input$concordance != "none") {
      LogCommentP("Calculate split concordance", 1)
      concCode <- switch(
        input$concordance,
        "p"   = paste0("SplitFrequency(", plottedTree,
                       ", trees) / length(trees)"),
        "qc"  = paste0("QuartetConcordance(", plottedTree, ", dataset)"),
        "clc" = paste0("ClusteringConcordance(", plottedTree, ", dataset)"),
        "phc" = paste0("PhylogeneticConcordance(", plottedTree, ", dataset)"),
        "mcc" = paste0("MutualClusteringConcordance(", plottedTree,
                       ", dataset)"),
        "spc" = paste0("SharedPhylogeneticConcordance(", plottedTree,
                       ", dataset)"),
        NULL
      )
      LogCodeP(paste0("concordance <- ", concCode))
      LogCommentP("Annotate splits by concordance", 1)
      LogCodeP("LabelSplits(",
              paste0("  tree = ", plottedTree, ","),
              "  labels = signif(concordance, 3),",
              "  col = SupportColor(concordance),",
              "  frame = \"none\",",
              "  pos = 3",
              ")")
    }
  }
  
  observeEvent(input$keepNTips, {
    if (!is.null(r$oldkeepNTips)) {
      if (!identical(input$keepNTips, r$oldkeepNTips)) {
        r$oldkeepNTips <- NULL
      }
    } else {
      LogMsg("Observed input$keepNTips -> ", EnC(input$keepNTips))
      r$keepNTips <- max(length(input$neverDrop), 3L,
                         min(input$keepNTips, TipsInTree()))
      UpdateOutgroupInput()
      UpdateDroppedTaxaDisplay()
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$neverDrop, {
    LogMsg("Observed input$neverDrop -> ", EnC(input$neverDrop))
    UpdateKeepNTipsRange()
    UpdateOutgroupInput()
    UpdateDroppedTaxaDisplay()
  }, ignoreInit = TRUE)
  
  observeEvent(input$outgroup, {
    if (!is.null(r$oldOutgroup)) {
      if (!identical(input$outgroup, r$oldOutgroup)) {
        r$oldOutgroup <- NULL
      }
    } else {
      LogMsg("Observed input$outgroup -> ", EnC(input$outgroup))
      r$outgroup <- input$outgroup
    }
  }, ignoreInit = TRUE)
  
  DatasetTips <- reactive(names(r$dataset))
  SearchTips <- reactive(setdiff(DatasetTips(), r$searchWithout))
  
  KeptTips <- reactive({
    LogMsg("KeptTips()")
    n <- r$keepNTips
    maxN <- length(tipLabels())
    if (is.na(n) || is.null(n)) {
      n <- maxN
    }
    if (n < 3L) {
      n <- 3L 
    }
    nNeverDrop <- length(input$neverDrop)
    if (n < nNeverDrop) {
      n <- nNeverDrop
    }
    nFromDropSeq <- n - nNeverDrop
    
    # Return:
    if (nFromDropSeq > length(dropSeq())) {
      c(input$neverDrop, dropSeq())
    } else {
      c(input$neverDrop, rev(dropSeq())[seq_len(nFromDropSeq)])
    }
  })
  
  DroppedTips <- reactive({
    LogMsg("DroppedTips()")
    if (length(KeptTips()) > 1) {
      setdiff(tipLabels(), KeptTips())
    } else {
      character(0)
    }
  })
  
  ConsensusPlot <- function() {
    LogMsg("ConsensusPlot()")
    on.exit(LogMsg("/ConsensusPlot()"))
    
    par(mar = rep(0, 4), cex = 0.9)
    kept <- KeptTips()
    dropped <- DroppedTips()
    
    if (length(dropped) &&
        length(input$excludedTip) &&
        nchar(input$excludedTip) &&
        input$excludedTip %in% tipLabels()) {
      
      if (length(setdiff(dropped, input$excludedTip))) {
        consTrees <- lapply(r$trees, DropTip,
                            setdiff(dropped, input$excludedTip))
      } else {
        consTrees <- r$trees
      }
      
      plotted <- TreeTools::RoguePlot(
        consTrees,
        input$excludedTip,
        p = consP(),
        edgeLength = 1,
        outgroupTips = r$outgroup,
        tip.color = TipCols()[intersect(consTrees[[1]]$tip.label, kept)]
      )
      r$plottedTree <- plotted$cons
      
      LabelConcordance()
    } else {
      without <- intersect(dropped, tipLabels()) # `dropped` might be outdated
      if (length(without)) {
      } else {
      }
      cons <- ConsensusWithout(r$trees, without, p = consP())
      cons <- UserRoot(cons)
      
      if (unitEdge()) {
        cons$edge.length <- rep.int(1, dim(cons$edge)[1])
      }
      cons <- SortEdges(cons)
      
      r$plottedTree <- cons
      plot(r$plottedTree, tip.color = TipCols()[intersect(cons$tip.label, kept)])
      LabelConcordance()
    }
  }
  
  LogConsensusPlot <- function() {
    BeginLogP()
    LogPar()
    dropped <- DroppedTips()
    
    if (length(dropped) &&
        length(input$excludedTip) &&
        nchar(input$excludedTip) &&
        input$excludedTip %in% tipLabels()) {
      
      LogCommentP("Prepare reduced consensus tree", 1)
      if (length(setdiff(dropped, input$excludedTip))) {
        LogCodeP(paste0("exclude <- ",
                       EnC(setdiff(dropped, input$excludedTip))))
        LogCodeP("consTrees <- lapply(trees, DropTip, exclude)")
        LogCodeP("labels <- setdiff(consTrees[[1]]$tip.label, exclude)")
      } else {
        LogCodeP("consTrees <- trees",
                "labels <- consTrees[[1]]$tip.label")
      }
      
      LogCommentP(paste0(
        "Colour tip labels according to their original 'instability' ",
        "(Smith 2022)")
      )
      LogCodeP(
        "tipCols <- Rogue::ColByStability(trees)",
        paste0(
          "tipCols <- tipCols[setdiff(labels, ",
          Enquote(input$excludedTip), ")]"
        )
      )
      LogCommentP(paste0(
        "Plot the reduced consensus tree, showing position of ",
        gsub("_", " ", input$excludedTip, fixed = TRUE))
      )
      LogCodeP("plotted <- RoguePlot(",
              "  trees = consTrees,",
              paste0("  tip = ", Enquote(input$excludedTip), ","),
              paste0("  p = ", consP(), ","),
              "  edgeLength = 1,",
              if(length(r$outgroup)) {
                  paste0("  outgroupTips = ", EnC(r$outgroup), ",")
              },
              "  tip.color = tipCols",
              ")")
      
      LogCommentP("Store tree to plot concordance")
      LogCodeP("plottedTree <- plotted$cons")
      
      LogConcordance()
    } else {
      without <- intersect(dropped, tipLabels()) # `dropped` might be outdated
      LogCommentP("Calculate consensus tree")
      if (length(without)) {
        LogCodeP(
          "cons <- ConsensusWithout(",
          "  trees,",
          paste0("  ", EnC(without), ","),
          paste0("  p = ", consP()),
          ")")
      } else {
        LogCodeP(paste0(
          "cons <- Consensus(trees, p = ", consP(), ")"
        ))
      }
      LogUserRoot(dropped = without)
      if (unitEdge()) {
        LogCodeP("cons$edge.length <- rep.int(1L, nrow(cons$edge))")
      }
      LogSortEdges("cons")
      LogCommentP("Plot consensus tree")
      LogCodeP(
        "tipCols <- Rogue::ColByStability(trees)[cons$tip.label]",
        "plot(cons, tip.color = tipCols)")
      LogConcordance("cons")
    }
  }
  
  PolEscVal <- reactive({
    LengthAdded(r$trees,
                r$dataset[tipLabels(), PlottedChar()],
                concavity())
  })
  
  CharacterwisePlot <- function() {
    par(mar = rep(0, 4), cex = 0.9)
    n <- PlottedChar()
    if (whichTree() > 0) {
      LogMsg("Plotting PlottedTree(", whichTree(), ", ", n, ")")
    }
    r$plottedTree <- PlottedTree()
    if (length(n) && n > 0L) {
      pc <- tryCatch({
        extraLen <- PolEscVal()
        roguishness <- if (max(extraLen) == 0) {
          "black"
        } else {
          hcl.colors(256, "inferno")[
            (192 * extraLen[r$plottedTree$tip.label] / max(extraLen)) + 1
          ]
        }
        PlotCharacter(
          if (whichTree() > 0) {
            MakeTreeBinary(r$plottedTree)
          } else {
            lapply(r$trees, function(t) MakeTreeBinary(UserRoot(t)))
          },
          r$dataset,
          n,
          edge.width = 2.5,
          updateTips = "updateTips" %in% input$mapDisplay,
          tip.color = roguishness,
          Display = function(tr) {
            tr <- UserRoot(tr)
            if (unitEdge()) {
              tr$edge.length <- rep.int(1, dim(tr$edge)[[1]])
            }
            SortEdges(tr)
          }
        )
        if (max(extraLen) > 0) {
          PlotTools::SpectrumLegend(
            "bottomleft", bty = "n",
            palette = hcl.colors(256, "inferno")[1:193],
            title = "Mean tree score\nimpact",
            title.font = 2,
            y.intersp = 1.42,
            legend = c(signif(4:1 * max(extraLen) / 4, 3), "No impact")
          )
        }
      },
      error = function (cond) {
        cli::cli_alert_danger(cond)
        Notification(type = "error",
                     "Could not match dataset to taxa in trees")
        ErrorPlot("Load dataset with\n", "character codings\n",
                  "for taxa on tree")
        return()
      }
      )
      
      LabelConcordance()
    } else {
      plot(r$plottedTree, tip.color = TipCols()[r$plottedTree$tip.label])
      TipColLegend()
    }
  }
  
  LogPar <- function() {
    LogCommentP("Set up plotting area")
    LogCodeP(c(
      "par(",
      "  mar = c(0, 0, 0, 0), # Zero margins",
      "  cex = 0.9            # Smaller font size",
      ")"
    ))
  }
  
  LogCharacterwisePlot <- function() {
    BeginLogP()
    LogPar()
    n <- PlottedChar()
    if (whichTree() > 0) {
      LogComment(paste("Select tree", whichTree(), "from tree set"))
    }
    LogPlottedTree()
    if (length(n) && n > 0L) {
      if (whichTree() > 0) {
        LogCommentP(paste("Map character", n, "onto tree", whichTree()))
      } else {
        LogCommentP(paste("Map character", n, "onto consensus tree"))
      }
      LogCodeP(
        "PlotCharacter(",
        if (whichTree() > 0) "  tree = MakeTreeBinary(plottedTree)," else 
          paste0("  tree = lapply(RootTree(trees, ", EnC(r$outgroup),
                 "), MakeTreeBinary),"),
        "  dataset = dataset,",
        paste0("  char = ", n, ","),
        paste0("  updateTips = ", "updateTips" %in% input$mapDisplay, ","),
        "  Display = function(tr) {",
        paste0("    tr <- RootTree(tr, ", EnC(r$outgroup), ")"),
        "    tr$edge.length <- rep.int(2, nrow(tr$edge))",
        "    SortTree(tr)",
        "  },",
        "  edge.width = 2.5",
        ")"
      )
      LogConcordance()
    } else {
      LogCommentP("Plot single tree")
      LogCodeP(
        "tipCols <- Rogue::ColByStability(trees)[plottedTree$tip.label]",
        "plot(plottedTree, tip.color = tipCols)"
      )
    }
  }
  
  MainPlot <- function() {
    if (AnyTrees()) {
      LogMsg("MainPlot()")
      switch(
        input$plotFormat,
        "cons" = {
          ConsensusPlot()
        },
        "clus" = {
          PlotClusterCons()
        },
        "ind" = {
          CharacterwisePlot()
        },
        "space" = {
          TreespacePlot()
        }
      ) # end switch
    }
  }
  ReactiveMainPlot <- reactive({MainPlot()})
  
  output$treePlot <- renderCachedPlot(
    ReactiveMainPlot(),
    cacheKeyExpr = { # Must be identical to RCode below
      switch(
        input$plotFormat,
        
        "clus" = list(r$treeHash, input$plotFormat,
                      r$keepNTips, input$excludedTip,
                      consP(),
                      input$neverDrop, r$outgroup,
                      input$distMeth,
                      input$concordance,
                      silThreshold()),
        "cons" = list(r$treeHash, input$plotFormat,
                      r$keepNTips, input$excludedTip,
                      consP(),
                      input$neverDrop, r$outgroup,
                      input$concordance),
        "ind" = list(PlottedChar(),
                     whichTree(),
                     input$concordance,
                     r$outgroup,
                     concavity(),
                     input$mapDisplay,
                     r$dataHash, r$treeHash), 
        "space" = list(r$treeHash, input$plotFormat,
                       min(dims(), nProjDim()),
                       TreeCols(),
                       treePch(),
                       input$distMeth,
                       ts$spaceCol(),
                       ts$mapLines(),
                       concavity(),
                       ts$spacePch(),
                       if (ts$spacePch() == "relat") ts$relators(),
                       silThreshold())
      )
    },
    sizePolicy = function(x) rep(input$plotSize, 2)
  )
  
  RCode <- bindCache(reactive({
    switch(
      input$plotFormat,
      "cons" = {
        LogConsensusPlot()
      },
      "clus" = {
        LogPlotClusterCons()
      },
      "ind" = {
        LogCharacterwisePlot()
      },
      "space" = {
        LogTreespacePlot()
      }
    )
    
    # Return:
    r$plotLog
  }),  # Must be identical to output$treePlot above
    switch(
      input$plotFormat,
      
      "clus" = list(r$treeHash, input$plotFormat,
                    r$keepNTips, input$excludedTip,
                    consP(),
                    input$neverDrop, r$outgroup,
                    input$distMeth,
                    input$concordance,
                    silThreshold()),
      "cons" = list(r$treeHash, input$plotFormat,
                    r$keepNTips, input$excludedTip,
                    consP(),
                    input$neverDrop, r$outgroup,
                    input$concordance),
      "ind" = list(PlottedChar(),
                   whichTree(),
                   input$concordance,
                   r$outgroup,
                   concavity(),
                   input$mapDisplay,
                   r$dataHash, r$treeHash), 
      "space" = list(r$treeHash, input$plotFormat,
                     min(dims(), nProjDim()),
                     TreeCols(),
                     treePch(),
                     input$distMeth,
                     ts$spaceCol(),
                     ts$mapLines(),
                     concavity(),
                     ts$spacePch(),
                     if (ts$spacePch() == "relat") ts$relators(),
                     silThreshold())
    )
  )
  
  UCFirst <- function (str) {
    paste0(toupper(substr(str, 1, 1)),
           substr(str, 2, nchar(str)))
  }
  
  nonAmbigContrast <- reactive({
    cont <- attr(r$dataset, "contrast")
    applic <- cont[, setdiff(colnames(cont), "-")]
    cont[rowSums(applic) == dim(applic)[[2]], ] <- 0
    
    # Return:
    cont
  })
  
  plottedTokens <- reactive({
    n <- PlottedChar()
    # `phyDat[,]` returns a new phyDat object with a different "contrast"
    # Hence we manually extract the compressed character tokens:
    phyColumn <- vapply(r$dataset, `[[`, integer(1),
                        attr(r$dataset, "index")[[n]], USE.NAMES = FALSE)
    tokens <- colSums(nonAmbigContrast()[phyColumn, ]) > 0L
    names(tokens[tokens])
  })
  
  output$charMapLegend <- bindCache(
    renderUI({
      n <- PlottedChar()
      if (length(n) && n > 0L && !is.null(r$chars)) {
      
        pal <- c("#00bfc6", "#ffd46f", "#ffbcc5", "#c8a500",
                 "#ffcaf5", "#d5fb8d", "#e082b4", "#25ffd3",
                 "#a6aaff", "#e6f3cc", "#67c4ff", "#9ba75c",
                 "#60b17f")
        
        states <- attr(r$chars, "state.labels")[[n]]
        tokens <- plottedTokens()
        appTokens <- setdiff(tokens, "-")
        datApp <- setdiff(attr(r$dataset, "levels"), "-")
        .State <- function (glyph, text = "Error?", col = "red") {
          if (is.numeric(glyph)) {
            if (glyph > length(appTokens)) {
              return(NULL)
            }
            level <- match(appTokens[[glyph]], datApp)
            text <- states[[level]]
            col <- pal[[level]]
            glyph <- appTokens[[glyph]]
          }
          
          tags$li(style = "margin-bottom: 2px;",
                  tags$span(glyph,
                            style = paste("display: inline-block;",
                                          "border: 1px solid;",
                                          "width: 1em;",
                                          "text-align: center;",
                                          "line-height: 1em;",
                                          "margin-right: 0.5em;",
                                          "background-color:", col, ";")
                  ),
                  tags$span(UCFirst(text)))
        }
        
        tagList(
          tags$h3(colnames(r$chars)[n]),
          tags$ul(style = "list-style: none;",
                  .State(1), .State(2), .State(3), .State(4), .State(5),
                  .State(6), .State(7), .State(8), .State(9),
                  .State(10), .State(11), .State(12), .State(13),
                  if ("-" %in% tokens) 
                    .State("-", "Inapplicable", "lightgrey"),
                  .State("?", "Ambiguous", "grey")
          )
        )
      }
    }),
    PlottedChar(),
    r$chars,
    r$dataset
  )
    
  
  output$charNotes <- bindCache(
    renderUI({
      n <- PlottedChar()
      if (length(n) && n > 0L
          && is.list(r$charNotes) && is.list(r$charNotes[[1]])
          && length(r$charNotes) >= n) {
      
        charNotes <- r$charNotes[[n]]
        description <- charNotes[[1]]
        notes <- charNotes[[2]]
        states <- attr(r$chars, "state.labels")[[n]]
        tokens <- plottedTokens()
        
        tagList(
          if (length(description) > 0) {
            tags$div(id = "char-description",
                     lapply(strsplit(description, "\n")[[1]], tags$p))
          },
          if (!is.null(notes)) tags$ul(class = "state-notes", {
            PrintNote <- function(note) {
              taxa <- names(note)[note]
              tags$li(class = "state-note",
                      tags$span(class = "state-note-label",
                                paste(gsub("_", " ", fixed = TRUE,
                                           taxa), collapse = ", ")),
                      tags$span(class = "state-note-detail",
                                notes[taxa[1]]))
            }
            
            DuplicateOf <- function(x) {
              duplicates <- duplicated(x)
              masters <- x[!duplicates]
              vapply(masters, function(d) x == d, logical(length(x)))
            }
            if (length(notes) == 1) {
              onlyOne <- TRUE
              names(onlyOne) <- names(notes)
              PrintNote(onlyOne)
            } else {
              notes <- notes[order(names(notes))]
              duplicates <- DuplicateOf(toupper(notes))
              apply(duplicates, 2, PrintNote)
            }
          }),
          if (!states[[1]] %in% c("", "''")
              && any(tokens == "-")) {
            tags$p(tags$em("Brazeau et al. (2019) advise that neomorphic (0/1) characters should not contain inapplicable tokens (-)."))
          }
        )
      }
    }),
    PlottedChar(),
    r$dataset,
    r$chars,
    r$charNotes
  )
  
  # QualityPlot, LogScore, output$pcQuality moved to mod_treespace.R
  
  observeEvent(input$clThresh, {
    classes <- c("meaningless", "weak", "good", "strong")
    liveClass <- classes[as.integer(cut(input$clThresh, c(0, 0.25, 0.5, 0.7, 1),
                                        include.lowest = TRUE, right = FALSE))]
    addClass("clThresh-label", liveClass)
    removeClass("clThresh-label", setdiff(classes, liveClass))
  })
