# Module: Consensus & Main Plot
#
# Absorbs consensus.R + residual clustering.R + consensus-related bindings
# from events.R. Owns the main plot dispatch, consensus tree plotting,
# character mapping, stability / rogue analysis, concordance, cluster
# consensus plotting, plot code logging, and associated UI updates.
#
# Owns inputs: consP, keepNTips, neverDrop, outgroup, concordance,
#   plottedChar, searchChar, mapDisplay, whichTree, excludedTip.
#
# Owns outputs: treePlot, charMapLegend, charNotes, branchLegend.
#
# Reactive args:
#   r                AppState reactiveValues
#   AnyTrees         reactive logical (data module)
#   HaveData         reactive logical (data module)
#   tipLabels        reactive character (data module)
#   nChars           reactive integer (data module)
#   TaxonOrder       reactive character (data module)
#   concavity        reactive (search module)
#   clusterings      reactive list (clustering module)
#   silThreshold     reactive numeric (clustering module)
#   LogClusterings   function (clustering module)
#   TreespacePlot    function (treespace module)
#   LogTreespacePlot function (treespace module)
#   dims             reactive integer (treespace module)
#   nProjDim         reactive integer (treespace module)
#   TreeCols         reactive character (treespace module)
#   treePch          reactive (treespace module)
#   ts_spaceCol      reactive character (treespace module)
#   ts_mapLines      reactive character (treespace module)
#   ts_spacePch      reactive character (treespace module)
#   ts_relators      reactive character (treespace module)
#   plotFormat       reactive character (top-level input)
#   plotSize         reactive integer (top-level input)
#   distMeth         reactive character (top-level input)
#   log_fns          named list of logging functions
#
# Returns:
#   MainPlot, RCode, UpdateKeepNTipsRange,
#   UpdateDroppedTaxaDisplay, UpdateOutgroupInput

# ---------------------------------------------------------------------------
# UI — returns named list for scattered placement in ui.R
# ---------------------------------------------------------------------------
consensus_ui <- function(id) {
  ns <- NS(id)
  list(
    tree_plot = plotOutput(ns("treePlot"), height = "600px"),

    which_tree = tagList(
      sliderInput(ns("whichTree"), "Tree to plot", value = 0L,
                  min = 0L, max = 1L, step = 1L),
      htmlOutput(ns("clusterLabel"), inline = TRUE)
    ),

    tree_plot_config = tagList(
      selectizeInput(ns("outgroup"), "Root on:", multiple = TRUE,
                     choices = list()),
      selectizeInput(
        ns("concordance"), "Split support:",
        choices = list(
          "None" = "none",
          "% trees containing" = "p",
          "Quartet concordance" = "qc",
          "Clustering concordance" = "clc",
          "Phylogenetic concordance" = "phc",
          "Mutual Clustering conc." = "mcc",
          "Shared Phylog. conc." = "spc"
        ))
    ),

    char_chooser = tagList(
      tags$div(
        numericInput(ns("plottedChar"), "Character to map:", value = 1L,
                     min = 0L, max = 1L, step = 1L, width = 200),
        selectizeInput(ns("searchChar"), "Search characters:",
                       multiple = FALSE, choices = list()),
        checkboxGroupInput(ns("mapDisplay"), "", list(
          "Align tips" = "tipsRight",
          "Infer tips" = "updateTips"
        )),
        style = "float: right; width: 200px; margin-left: 2em;"
      ),
      htmlOutput(ns("charMapLegend")),
      htmlOutput(ns("charNotes"))
    ),

    cons_config = tagList(
      tags$div(style = "float: right; width: 200px; margin-left: 2em;",
        sliderInput(ns("consP"), "Majority:", value = 1,
                    min = 0.5, max = 1, width = 200),
        numericInput(ns("keepNTips"), "Tips to show:", value = 0L,
                     min = 3L, max = 2L, step = 1L, width = 200),
        selectizeInput(ns("neverDrop"), "Never drop:", multiple = TRUE,
                       choices = c())
      ),
      tags$div(id = "consLegend",
        tags$span(id = "instabLegend",
          tagList(
            tags$span(class = "legendLeft", "Stable"),
            tags$span(class = "infernoScale legendBar", "\ua0"),
            tags$span(class = "legendRight", "Unstable")
          )
        ),
        # Wrapper keeps top-level id for ShowConfigs show/hide
        tags$span(id = "branchLegend",
          htmlOutput(ns("branchLegend"), inline = TRUE)
        )
      ),
      tags$div(id = "droppedTips",
        selectInput(ns("excludedTip"), "Show excluded tip", choices = list())
      ),
      tags$div(id = "droppedList", style = "float: left;")
    )
  )
}

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
consensus_server <- function(id, r,
                             AnyTrees, HaveData, tipLabels, nChars, TaxonOrder,
                             concavity,
                             clusterings, silThreshold, LogClusterings,
                             TreespacePlot, LogTreespacePlot,
                             dims, nProjDim, TreeCols, treePch,
                             ts_spaceCol, ts_mapLines, ts_spacePch, ts_relators,
                             plotFormat, plotSize, distMeth,
                             log_fns) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Unpack logging
    LogMsg      <- log_fns$LogMsg
    LogComment  <- log_fns$LogComment
    LogCode     <- log_fns$LogCode
    LogCommentP <- log_fns$LogCommentP
    LogCodeP    <- log_fns$LogCodeP
    LogIndent   <- log_fns$LogIndent
    BeginLogP   <- log_fns$BeginLogP
    LogExprP    <- log_fns$LogExprP

    ############################################################################
    # Cross-module shinyjs helpers (target top-level DOM ids)
    ############################################################################

    parentShow <- function(id) {
      runjs(paste0("$('#", id, "').removeClass('shinyjs-hide').show()"))
    }
    parentHide <- function(id) {
      runjs(paste0("$('#", id, "').hide()"))
    }
    parentHtml <- function(id, html) {
      escaped <- gsub("'", "\\'", html, fixed = TRUE)
      runjs(paste0("$('#", id, "').html('", escaped, "')"))
    }

    ############################################################################
    # Core helpers
    ############################################################################

    UserRoot <- function(tree) {
      outgroupTips <- intersect(r$outgroup, tree$tip.label)
      if (length(outgroupTips)) {
        RootTree(tree, outgroupTips)
      } else {
        tree
      }
    }

    LogUserRoot <- function(tree = "cons", dropped = character(0)) {
      outgroupTips <- setdiff(r$outgroup, dropped)
      if (length(outgroupTips)) {
        LogCommentP("Root tree")
        LogCodeP(paste0(tree, " <- RootTree(", tree, ", ",
                        EnC(outgroupTips), ")"))
      }
    }

    unitEdge <- reactive(TRUE)

    SortEdges <- function(tr, force = FALSE) {
      if (force || r$sortTrees) {
        SortTree(tr, order = TaxonOrder())
      } else {
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

    LogPar <- function() {
      LogCommentP("Set up plotting area")
      LogCodeP(c(
        "par(",
        "  mar = c(0, 0, 0, 0), # Zero margins",
        "  cex = 0.9            # Smaller font size",
        ")"
      ))
    }

    UCFirst <- function(str) {
      paste0(toupper(substr(str, 1, 1)),
             substr(str, 2, nchar(str)))
    }

    TipsInTree <- reactive({
      if (AnyTrees()) {
        length(r$trees[[1]]$tip.label)
      } else {
        0L
      }
    })

    ############################################################################
    # Debounced input reactives
    ############################################################################

    PlottedChar <- debounce(reactive({
      typed <- max(0L, as.integer(input$plottedChar), na.rm = TRUE)
      if (nChars() > 0 && typed > nChars()) {
        Notification(type = "warning",
                     paste("Dataset contains", nChars(), "characters."))
        updateNumericInput(session, "plottedChar", value = nChars())
      }
      min(typed, nChars())
    }), aJiffy)

    whichTree <- debounce(reactive(input$whichTree), aJiffy)

    output$clusterLabel <- renderUI({
      wt <- whichTree()
      if (is.null(wt) || wt < 1L) return(NULL)
      cl <- clusterings()
      if (cl$n < 2L) return(NULL)
      clId <- cl$cluster[wt]
      col <- palettes[[min(length(palettes), cl$n)]][clId]
      tags$span(
        paste0("Cluster ", clId),
        style = paste0("color:", col, ";font-weight:bold;margin-left:4px;")
      )
    })

    consP <- debounce(reactive(signif(input$consP)), 50)

    ############################################################################
    # Stability / rogue analysis
    ############################################################################

    Instab <- reactive({
      TipInstability(r$trees)
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
        rogues[!rogues$taxon %in% input$neverDrop, ]
      } else {
        data.frame(num = 0, taxNum = NA_integer_, taxon = NA_character_,
                   rawImprovement = NA_real_, IC = 0)
      }
    }), r$treeHash, input$neverDrop, consP())

    dropSeq <- reactive({
      LogMsg("dropSeq()")
      Rogues()$taxon[-1]
    })

    nNonRogues <- reactive({
      LogMsg("nNonRogues()")
      on.exit(LogMsg("nNonRogues: ", nrow(Rogues()) - which.max(Rogues()$IC)))
      nrow(Rogues()) - which.max(Rogues()$IC)
    })

    TipCols <- reactive(stableCol())

    TipColLegend <- function() {
      PlotTools::SpectrumLegend(
        "bottomleft", horiz = TRUE, inset = 0.01, bty = "n", xpd = NA,
        palette = hcl.colors(131, "inferno")[1:101],
        legend = c("Stable", "Unstable"),
        title = "Leaf stability",
        title.font = 2
      )
    }

    ############################################################################
    # Tip subsetting
    ############################################################################

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

    ############################################################################
    # Concordance
    ############################################################################

    concordance <- bindCache(reactive({
      LogMsg("concordance()")
      if (input$concordance %in% c("qc", "mcc", "spc", "clc", "phc") &&
          !setequal(TipLabels(r$plottedTree), names(r$dataset))) {
        return(NULL)
      }
      switch(input$concordance,
             "p"   = SplitFrequency(r$plottedTree, r$trees) / length(r$trees),
             "qc"  = QuartetConcordance(r$plottedTree, r$dataset),
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
          inherits(r$plottedTree, "phylo")) {
        conc <- concordance()
        if (is.null(conc)) return(invisible())
        LabelSplits(r$plottedTree, signif(conc, 3),
                    col = SupportColor(conc),
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

    ############################################################################
    # Tree plotting
    ############################################################################

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
          plottedTree$edge.length <-
            rep_len(2, dim(plottedTree[["edge"]])[[1]])
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

    ############################################################################
    # Consensus plot
    ############################################################################

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
        without <- intersect(dropped, tipLabels())
        cons <- ConsensusWithout(r$trees, without, p = consP())
        cons <- UserRoot(cons)

        if (unitEdge()) {
          cons$edge.length <- rep.int(1, dim(cons$edge)[1])
        }
        cons <- SortEdges(cons)

        r$plottedTree <- cons
        plot(r$plottedTree,
             tip.color = TipCols()[intersect(cons$tip.label, kept)])
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
        without <- intersect(dropped, tipLabels())
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

    ############################################################################
    # Character-wise plot
    ############################################################################

    PolEscVal <- reactive({
      tl <- tipLabels()
      dl <- names(r$dataset)
      # Skip if taxa don't match exactly: tipLabels() may include taxa absent
      # from the dataset (e.g. trees loaded from a superset dataset), causing
      # a matrix-dimension mismatch inside LengthAdded / TreeLength.
      if (!setequal(tl, dl)) return(NULL)
      LengthAdded(r$trees,
                  r$dataset[tl, PlottedChar()],
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
              if ("tipsRight" %in% input$mapDisplay) {
                # Cladogram: tips aligned to the right
                tr$edge.length <- NULL
              } else {
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
        error = function(cond) {
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

    ############################################################################
    # Cluster consensus plot (absorbed from clustering.R)
    ############################################################################

    # Per-edge colors for cluster consensus: unique splits get the full
    # cluster color; splits shared by other clusters fade towards grey.
    ClusterEdgeCols <- function(tree, cluster_col, all_splits, cluster_idx) {
      n_tip <- Ntip(tree)
      n_edge <- nrow(tree$edge)
      edge_col <- rep(cluster_col, n_edge)

      my_splits <- all_splits[[cluster_idx]]
      n_clusters <- length(all_splits)
      if (length(my_splits) == 0 || n_clusters < 2) return(edge_col)

      other_idx <- setdiff(seq_len(n_clusters), cluster_idx)
      split_nodes <- as.integer(names(my_splits))

      shared <- integer(length(my_splits))
      for (j in other_idx) {
        if (length(all_splits[[j]]) > 0) {
          shared <- shared + as.integer(my_splits %in% all_splits[[j]])
        }
      }
      uniqueness <- 1 - shared / length(other_idx)

      grey_rgb <- col2rgb("grey70")[, 1]
      col_rgb <- col2rgb(cluster_col)[, 1]
      edge_child <- tree$edge[, 2]
      for (e in seq_len(n_edge)) {
        child <- edge_child[e]
        if (child > n_tip) {
          sidx <- match(child, split_nodes)
          if (!is.na(sidx)) {
            u <- uniqueness[sidx]
            bl <- grey_rgb + (col_rgb - grey_rgb) * u
            edge_col[e] <- rgb(bl[1], bl[2], bl[3], maxColorValue = 255)
          }
        }
      }
      edge_col
    }

    PlotClusterCons <- function() {
      LogMsg("PlotClusterCons()")
      on.exit(LogMsg("/PlotClusterCons()"))

      cl <- clusterings()

      kept <- KeptTips()
      dropped <- if (length(kept) > 1) {
        setdiff(TipLabels(r$trees[[1]]), kept)
      } else {
        character(0)
      }
      par(mar = c(0.2, 0, 0.2, 0), xpd = NA)
      if (cl$sil > silThreshold()) {
        nRow <- ceiling(cl$n / 3)
        r$plottedTree <- vector("list", cl$n)
        par(mfrow = c(nRow, ceiling(cl$n / nRow)))

        # Phase 1: compute all cluster consensus trees
        all_cons <- vector("list", cl$n)
        for (i in seq_len(cl$n)) {
          cons <- ConsensusWithout(r$trees[cl$cluster == i], dropped,
                                   p = consP())
          cons <- UserRoot(cons)
          if (unitEdge()) {
            cons$edge.length <- rep.int(1, dim(cons$edge)[1])
          }
          all_cons[[i]] <- SortEdges(cons)
        }
        all_splits <- lapply(all_cons, as.Splits)

        # Phase 2: plot with uniqueness-based edge coloring
        for (i in seq_len(cl$n)) {
          col <- palettes[[min(length(palettes), cl$n)]][i]
          PutTree(r$trees)
          PutData(cl$cluster)

          cons <- all_cons[[i]]
          r$plottedTree[[i]] <- cons
          edge_col <- ClusterEdgeCols(cons, col, all_splits, i)
          plot(cons, edge.width = 2, font = 3, cex = 0.83,
               edge.color = edge_col, tip.color = TipCols()[cons$tip.label])
          legend("topright", paste0("Cluster ", i), pch = 15, col = col,
                 pt.cex = 1.5, bty = "n")
          LabelConcordance()
        }
      } else {
        PutTree(r$trees)
        cons <- ConsensusWithout(r$trees, dropped, p = consP())
        cons <- UserRoot(cons)
        if (unitEdge()) {
          cons$edge.length <- rep.int(1, dim(cons$edge)[1])
        }
        cons <- SortEdges(cons)
        r$plottedTree <- cons
        plot(cons, edge.width = 2, font = 3, cex = 0.83,
             edge.color = palettes[[1]],
             tip.color = TipCols()[cons$tip.label])
        LabelConcordance()
        legend("topright", "No clustering", pch = 16, col = palettes[[1]],
               bty = "n")
      }
    }

    LogPlotClusterCons <- function() {
      LogMsg("PlotClusterCons()")
      on.exit(LogMsg("/PlotClusterCons()"))

      BeginLogP()

      cl <- clusterings()
      LogClusterings()

      kept <- KeptTips()
      dropped <- if (length(kept) > 1) {
        setdiff(TipLabels(r$trees[[1]]), kept)
      } else {
        character(0)
      }
      if (cl$sil > silThreshold()) {
        nRow <- ceiling(cl$n / 3)
        LogCommentP("Plot consensus of each tree cluster", 2)
        LogCodeP(paste0(
          "par(mfrow = c(", nRow, ", ",
          ceiling(cl$n / nRow), "))",
          " # Plotting area layout"
        ))
        LogCodeP(
          paste0(
            "tipCols <- Rogue::ColByStability(trees)",
            " # Colour tips by stability"
          )
        )
        LogCommentP("Compute all cluster consensus trees:", 1)
        LogCodeP(
          paste0("allCons <- lapply(seq_len(", cl$n, "), function(i) {"),
          "  clusterTrees <- trees[clustering == i]",
          "  cons <- ConsensusWithout(",
          "    trees = clusterTrees,",
          paste0("    tip = ", EnC(dropped), ","),
          paste0("    p = ", consP()),
          "  )"
        )
        LogUserRoot(dropped = dropped)
        if (unitEdge()) {
          LogExprP("  cons$edge.length <- rep.int(1, nrow(cons$edge))")
        }
        LogCodeP("  TreeTools::SortTree(cons)", "})")
        LogCommentP(paste0(
          "Compare splits across clusters to highlight unique edges"
        ))
        LogCodeP("allSplits <- lapply(allCons, TreeTools::as.Splits)")
        LogCommentP("Plot each consensus tree in turn:", 1)
        LogCodeP(paste0("for (i in seq_len(", cl$n, ")) {"))
        LogIndent(+2)
        LogCodeP(
          "cons <- allCons[[i]]",
          "nTip <- ape::Ntip(cons)",
          "mySplits <- allSplits[[i]]",
          paste0("otherIdx <- setdiff(seq_len(", cl$n, "), i)"),
          "shared <- integer(length(mySplits))",
          "for (j in otherIdx) {",
          "  if (length(allSplits[[j]]) > 0)",
          "    shared <- shared + (mySplits %in% allSplits[[j]])",
          "}",
          "uniqueness <- 1 - shared / length(otherIdx)",
          "greyRgb <- col2rgb(\"grey70\")[, 1]",
          "colRgb <- col2rgb(clusterCol[i])[, 1]",
          "edgeCol <- rep(clusterCol[i], nrow(cons$edge))",
          "splitNodes <- as.integer(names(mySplits))",
          "for (e in seq_len(nrow(cons$edge))) {",
          "  child <- cons$edge[e, 2]",
          "  if (child > nTip) {",
          "    si <- match(child, splitNodes)",
          "    if (!is.na(si)) {",
          "      bl <- greyRgb + (colRgb - greyRgb) * uniqueness[si]",
          "      edgeCol[e] <- rgb(bl[1], bl[2], bl[3], maxColorValue = 255)",
          "    }",
          "  }",
          "}"
        )
        LogCodeP("plot(",
                 "  cons,",
                 "  edge.width = 2,",
                 "  font = 3,",
                 "  cex = 0.83,",
                 "  edge.color = edgeCol,",
                 "  tip.color = tipCols[cons$tip.label]",
                 ")")
        LogCodeP("legend(",
                 "  \"bottomright\",",
                 "  paste(\"Cluster\", i),",
                 "  pch = 15,",
                 "  pt.cex = 1.5,",
                 "  col = clusterCol[i],",
                 "  bty = \"n\"",
                 ")")
        LogConcordance("cons")
        LogIndent(-2)
        LogCodeP("}")
      } else {
        LogCommentP("No clustering structure: Plot consensus tree")
        LogCodeP(
          if (length(dropped)) {
            c("cons <- ConsensusWithout(",
              "  trees = trees,",
              paste0("  tip = ", EnC(dropped), ","),
              paste0("  p = ", consP()),
              ")"
            )
          } else {
            paste0("cons <- Consensus(trees, p = ", consP(), ")")
          }
        )
        LogUserRoot("cons", dropped = dropped)
        if (unitEdge()) {
          LogCommentP("Set unit edge length", 0)
          LogCodeP("cons$edge.length <- rep.int(1, nrow(cons$edge))")
        }
        LogSortEdges("cons")
        LogCodeP("plottedTree <- cons # Store for future reference")

        LogCodeP("tipCols <- Rogue::ColByStability(trees)[cons$tip.label]")
        LogCommentP("Plot consensus tree")
        LogCodeP(
          "plot(",
          "  cons,",
          "  edge.width = 2, # Widen lines",
          "  font = 3,       # Italicize labels",
          "  cex = 0.83,     # Shrink tip font size",
          "  tip.color = tipCols",
          ")"
        )
        LogConcordance()
      }
    }

    ############################################################################
    # Main plot dispatch
    ############################################################################

    MainPlot <- function() {
      if (AnyTrees()) {
        LogMsg("MainPlot()")
        switch(
          plotFormat(),
          "cons" = ConsensusPlot(),
          "clus" = PlotClusterCons(),
          "ind"  = CharacterwisePlot(),
          "space" = TreespacePlot()
        )
      }
    }
    ReactiveMainPlot <- reactive({ MainPlot() })

    output$treePlot <- renderCachedPlot(
      ReactiveMainPlot(),
      cacheKeyExpr = {
        switch(
          plotFormat(),

          "clus" = list(r$treeHash, plotFormat(),
                        r$keepNTips, input$excludedTip,
                        consP(),
                        input$neverDrop, r$outgroup,
                        distMeth(),
                        input$concordance,
                        silThreshold()),
          "cons" = list(r$treeHash, plotFormat(),
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
          "space" = list(r$treeHash, plotFormat(),
                         min(dims(), nProjDim()),
                         TreeCols(),
                         treePch(),
                         distMeth(),
                         ts_spaceCol(),
                         ts_mapLines(),
                         concavity(),
                         ts_spacePch(),
                         if (ts_spacePch() == "relat") ts_relators(),
                         silThreshold())
        )
      },
      sizePolicy = function(x) rep(plotSize(), 2)
    )

    ############################################################################
    # R code logging for plots (for downloads)
    ############################################################################

    RCode <- bindCache(reactive({
      switch(
        plotFormat(),
        "cons" = LogConsensusPlot(),
        "clus" = LogPlotClusterCons(),
        "ind"  = LogCharacterwisePlot(),
        "space" = LogTreespacePlot()
      )
      r$plotLog
    }),
      switch(
        plotFormat(),

        "clus" = list(r$treeHash, plotFormat(),
                      r$keepNTips, input$excludedTip,
                      consP(),
                      input$neverDrop, r$outgroup,
                      distMeth(),
                      input$concordance,
                      silThreshold()),
        "cons" = list(r$treeHash, plotFormat(),
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
        "space" = list(r$treeHash, plotFormat(),
                       min(dims(), nProjDim()),
                       TreeCols(),
                       treePch(),
                       distMeth(),
                       ts_spaceCol(),
                       ts_mapLines(),
                       concavity(),
                       ts_spacePch(),
                       if (ts_spacePch() == "relat") ts_relators(),
                       silThreshold())
      )
    )

    ############################################################################
    # Character map legend + notes (htmlOutput)
    ############################################################################

    nonAmbigContrast <- reactive({
      cont <- attr(r$dataset, "contrast")
      applic <- cont[, setdiff(colnames(cont), "-")]
      cont[rowSums(applic) == dim(applic)[[2]], ] <- 0
      cont
    })

    plottedTokens <- reactive({
      n <- PlottedChar()
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
          .State <- function(glyph, text = "Error?", col = "red") {
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
              tags$p(tags$em(paste0(
                "Brazeau et al. (2019) advise that neomorphic (0/1) ",
                "characters should not contain inapplicable tokens (-)."
              )))
            }
          )
        }
      }),
      PlottedChar(),
      r$dataset,
      r$chars,
      r$charNotes
    )

    ############################################################################
    # Branch legend (from events.R)
    ############################################################################

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
        consTrees <- lapply(r$trees, DropTip,
                            setdiff(dropped, input$excludedTip))
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

    ############################################################################
    # Update functions (from events.R) — used by data module via callbacks
    ############################################################################

    UpdateKeepNTipsRange <- reactive({
      if (AnyTrees() && "consConfig" %in% r$visibleConfigs) {
        nTip <- TipsInTree()
        # isolate() prevents re-triggering when user manually edits keepNTips
        currentInput <- isolate(input$keepNTips)
        LogMsg("UpdateKeepNTipsRange(", currentInput, " -> ", nTip, ")")
        r$keepNTips <- nNonRogues()
        if (r$keepNTips != currentInput) {
          r$oldkeepNTips <- currentInput
        }
        updateNumericInput(session, inputId = "keepNTips",
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
        updateSelectInput(session, inputId = "excludedTip",
                          choices = dropList,
                          selected = if (input$excludedTip %in% DroppedTips())
                            input$excludedTip else dropSeq()[1])
        # droppedList is a top-level div — use runjs
        droppedHtml <- paste0(
          "<label class=\"control-label\">Dropped tips:</label>",
          "<ul>",
          paste0("<li style=\"color: ", TipCols()[dropList], "\">",
                 dropList, "</li>", collapse = "\r\n"),
          "</ul>")
        parentHtml("droppedList", droppedHtml)
      }
    })

    UpdateDroppedTaxaDisplay <- reactive({
      LogMsg("UpdateDroppedTaxaDisplay()")
      if ("consConfig" %in% r$visibleConfigs) {
        if (length(DroppedTips())) {
          UpdateExcludedTipsInput()
          if ("droppedTips" %in% r$visibleConfigs) {
            parentShow("droppedTips")
          }
          if ("droppedList" %in% r$visibleConfigs) {
            parentShow("droppedList")
          }
        } else {
          parentHide("droppedTips")
          parentHide("droppedList")
        }
      }
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
          updateSelectizeInput(
            session,
            inputId = "outgroup",
            selected = r$outgroup,
            choices = KeptTips()
          )
        } else {
          # Update available choices without changing selection, to avoid
          # stealing focus from the widget after user adds/removes a taxon.
          updateSelectizeInput(
            session,
            inputId = "outgroup",
            choices = KeptTips()
          )
        }
      }
    })

    # Force reactive UI-update functions to run whenever their dependencies
    # change. Without these observers, the reactives are never consumed on
    # initial load, leaving inputs with their placeholder values.
    observe(UpdateKeepNTipsRange())
    observe(UpdateOutgroupInput())

    ############################################################################
    # Input observers
    ############################################################################

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

    observeEvent(consP(), {
      if (AnyTrees()) {
        LogMsg("Observed consP()")
        UpdateKeepNTipsRange()
        UpdateDroppedTaxaDisplay()
        r$concordance <- list()
      }
    }, ignoreInit = TRUE)

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

    observeEvent(r$visibleConfigs, {
      UpdateDroppedTaxaDisplay()
    })

    ############################################################################
    # Cross-module reactivity: observe state changes -> update module inputs
    # Replaces parent_session updateXxxInput calls from mod_data.R
    ############################################################################

    # When dataset changes: update plottedChar range + searchChar choices
    observeEvent(r$dataHash, {
      if (HaveData()) {
        n <- nChars()
        updateNumericInput(session, "plottedChar",
                           min = 0L, max = n, value = 1L)
        updateSelectizeInput(session, "searchChar",
                             choices = paste0(seq_len(n), ": ",
                                              colnames(r$chars)),
                             selected = "",
                             server = TRUE)
      } else {
        updateNumericInput(session, "plottedChar",
                           min = 0L, max = 0L, value = 0L)
        updateSelectizeInput(session, "searchChar", choices = NULL)
      }
    }, ignoreInit = TRUE)

    # When trees change: update whichTree slider range + neverDrop choices
    observeEvent(r$treeHash, {
      if (AnyTrees()) {
        nTr <- length(r$trees)
        updateSliderInput(session, "whichTree",
                          min = 0L, max = nTr, value = 0L)
        updateSelectizeInput(session, "neverDrop",
                             choices = tipLabels(),
                             selected = input$neverDrop)
        showElement("keepNTips")
        showElement("neverDrop")
      } else {
        hideElement("keepNTips")
        hideElement("neverDrop")
      }
    }, ignoreInit = TRUE)

    # Resize plot via CSS when plotSize changes
    observe({
      px <- paste0("'", plotSize(), "px'")
      runjs(paste0("$('#", ns("treePlot"), "').css({height: ",
                   px, ", width: ", px, "});"))
    })

    ############################################################################
    # Return values for other modules / server.R
    ############################################################################

    list(
      MainPlot                 = MainPlot,
      RCode                    = RCode,
      UpdateKeepNTipsRange     = UpdateKeepNTipsRange,
      UpdateDroppedTaxaDisplay = UpdateDroppedTaxaDisplay,
      UpdateOutgroupInput      = UpdateOutgroupInput
    )
  })
}
