# Module: Tree space visualization
#
# Absorbs treespace.R + plotsettings.R. Owns inputs: spaceDim, spaceCol,
# spacePch, relators, mapLines. Reads: r$trees, r$treeHash, clusterings(),
# silThreshold(), scores(), concavity(). Receives top-level distMeth and
# plotFormat as reactive args.
#
# Returns a list of reactives consumed by other source'd server files:
#   distances, mapping, dims, nProjDim, TreeCols, treePch,
#   saveDetails, TreespacePlot, LogTreespacePlot, mstEnds

treespace_ui <- function(id) {
  ns <- NS(id)
  tags$div(
    id = "spaceConfig",
    tags$div(id = "spaceLegend",
             style = "float: left;",
             plotOutput(outputId = ns("pcQuality"),
                        height = "72px", width = "240px"),
             htmlOutput(ns("stressLegend"), inline = TRUE)
    ),
    tags$div(
      style = "float: right; width: 200px; margin-left: 2em;",
      sliderInput(ns("spaceDim"), "Dimensions:", value = 5,
                  min = 1, max = 12, step = 1, width = 200),
      selectInput(ns("spaceCol"), "Colour trees by:",
                  list("Cluster membership" = "clust",
                       "Parsimony score" = "score",
                       "When first found" = "firstHit")),
      selectInput(ns("spacePch"), "Plotting symbols:",
                  selected = "relat",
                  list("Cluster membership" = "clust",
                       "Relationships" = "relat",
                       "Tree name" = "name")),
      selectizeInput(ns("relators"), "Show relationship between:",
                     choices = list(), multiple = TRUE),
    ),
  )
}

#' @param id Module namespace id.
#' @param r AppState reactiveValues.
#' @param clusterings Reactive returning clustering result list.
#' @param silThreshold Reactive returning silhouette threshold.
#' @param scores Reactive returning tree scores.
#' @param concavity Reactive returning concavity value.
#' @param distMeth Reactive wrapping top-level \code{input$distMeth}.
#' @param plotFormat Reactive wrapping top-level \code{input$plotFormat}.
#' @param log_fns Named list of logging functions from logging.R:
#'   BeginLogP, LogCommentP, LogCodeP, LogIndent, LogClusterings.
treespace_server <- function(id, r, clusterings, silThreshold, scores,
                             concavity, distMeth, plotFormat, log_fns) {
  moduleServer(id, function(input, output, session) {

    # Unpack logging functions
    BeginLogP      <- log_fns$BeginLogP
    LogCommentP    <- log_fns$LogCommentP
    LogCodeP       <- log_fns$LogCodeP
    LogIndent      <- log_fns$LogIndent
    LogClusterings <- log_fns$LogClusterings

    ############################################################################
    # Plot settings (from plotsettings.R)
    ############################################################################

    spaceCex <- reactive(1.7)
    spaceLwd <- reactive(2)

    FirstHit <- reactive({
      r$trees <- WhenFirstHit(r$trees)
      attr(r$trees, "firstHit")
    })

    LogFirstHit <- function() {
      LogCodeP("whenHit <- gsub(\"(seed|start|ratch\\\\d+|final)_\\\\d+\", \"\\\\1\",
              names(trees), perl = TRUE)")
      LogCodeP("attr(trees, \"firstHit\") <- table(whenHit)[unique(whenHit)]")
    }

    FirstHitCols <- reactive({
      if (is.null(FirstHit())) {
        palettes[[1]]
      } else {
        hcl.colors(length(FirstHit()), "viridis")
      }
    })

    LogFirstHitCols <- reactive({
      if (is.null(FirstHit())) {
        paste0(palettes[[1]], " # Arbitrarily")
      } else {
        "hcl.colors(length(firstHit), \"viridis\")"
      }
    })

    TreeCols <- reactive({
      switch(
        input$spaceCol,
        "clust" = {
          cl <- clusterings()
          if (cl$sil > silThreshold()) {
            palettes[[min(length(palettes), cl$n)]][cl$cluster]
          } else {
            palettes[[1]]
          }
        }, "score" = {
          if (is.null(scores()) || length(unique(scores())) == 1L) {
            palettes[[1]]
          } else {
            norm <- scores() - min(scores())
            norm <- (length(badToGood) - 1L) * norm / max(norm)
            rev(badToGood)[1 + norm]
          }
        }, "firstHit" = {
          if (is.null(FirstHit())) {
            Notification("Data not available; were trees loaded from file?",
                         type = "warning")
            palettes[[1]]
          } else {
            rep(FirstHitCols(), FirstHit())
          }
        },
        "black"
      )
    })

    LogTreeCols <- reactive({
      beige <- paste0("treeCols <- ", Enquote(palettes[[1]]), " # Arbitrarily")
      switch(
        input$spaceCol,
        "clust" = {
          cl <- clusterings()
          if (cl$sil > silThreshold()) {
            paste0("treeCols <- ",
                   EnC(palettes[[min(length(palettes), cl$n)]]),
                   "[clustering]")
          } else {
            beige
          }
        }, "score" = {
          if (is.null(scores()) || length(unique(scores())) == 1L) {
            beige
          } else {
            c(paste0("scores <- TreeLength(trees, dataset, concavity = ",
                     Enquote(concavity()), ")"),
              "normalized <- scores - min(scores)",
              "normalized <- 107 * normalized / max(normalized)",
              "goodToBad <- hcl.colors(108, \"Temps\")",
              "treeCols <- goodToBad[1 + normalized]"
            )
          }
        }, "firstHit" = {
          if (is.null(FirstHit())) {
            beige
          } else {
            c("trees <- WhenFirstHit(trees)",
              "firstHit <- attr(trees, \"firstHit\")",
              paste0("treeCols <- rep(", LogFirstHitCols(), ", firstHit))")
            )
          }
        },
        "treeCols <- black"
      )
    })

    treeNameClustering <- reactive({
      ClusterStrings(names(r$trees))
    })

    treePch <- reactive({
      switch(
        input$spacePch,
        "clust" = {
          cl <- clusterings()
          if (cl$sil > silThreshold()) {
            cl$cluster - 1
          } else {
            16
          }
        }, "relat" = {
          quartet <- input$relators
          if (length(quartet) == 4) {
            QuartetResolution(r$trees, input$relators)
          } else {
            Notification("Select four taxa to show relationships")
            0
          }
        }, "name" = {
          if (is.null(names(r$trees))) {
            Notification("Trees lack names", type = "warning")
            16
          } else {
            indices <- treeNameClustering()
            c(1, 3, 4, 2, seq_len(max(indices))[-(1:4)])[indices]
          }
        }, 0)
    })

    LogTreePch <- function() {
      switch(
        input$spacePch,
        "clust" = {
          cl <- clusterings()
          if (cl$sil > silThreshold()) {
            "cl$cluster - 1"
          } else {
            "16 # No clustering structure: Use filled circle"
          }
        }, "relat" = {
          quartet <- input$relators
          if (length(quartet) == 4) {
            paste0("QuartetResolution(trees, ", EnC(input$relators), ")")
          } else {
            "0 # Square"
          }
        }, "name" = {
          if (is.null(names(r$trees))) {
            "16 # Filled circle"
          } else {
            "ClusterStrings(names(trees))"
          }
        }, "0 # Square")
    }

    maxProjDim <- reactive({
      min(12, max(0L, length(r$trees) - 1L))
    })

    nProjDim <- reactive({
      dim(mapping())[2]
    })

    dims <- debounce(reactive({
      min(input$spaceDim, maxProjDim())
    }), 400)

    Quartet <- function(...) {
      if (!requireNamespace("Quartet", quietly = TRUE)) {
        Notification("Installing required package \"Quartet\"",
                     type = "warning", duration = 20)
        install.packages("Quartet")
      }
      as.dist(Quartet::QuartetDivergence(
        Quartet::ManyToManyQuartetAgreement(...), similarity = FALSE))
    }

    distances <- bindCache(reactive({
      LogMsg("distances(): ", distMeth())
      if (length(r$trees) > 1L) {
        Dist <- switch(distMeth(),
                       "cid" = TreeDist::ClusteringInfoDistance,
                       "pid" = TreeDist::PhylogeneticInfoDistance,
                       "msid" = TreeDist::MatchingSplitInfoDistance,
                       "rf" = TreeDist::RobinsonFoulds,
                       "qd" = Quartet)
        withProgress(
          message = "Initializing distances...", value = 0.99,
          Dist(r$trees)
        )
      } else {
        matrix(0, 0, 0)
      }
    }), distMeth(), r$treeHash)

    LogDistances <- function() {
      LogCommentP("Compute tree distances")
      LogCodeP(switch(
        distMeth(),
        "cid" = "dists <- TreeDist::ClusteringInfoDistance(trees)",
        "pid" = "dists <- TreeDist::PhylogeneticInfoDistance(trees)",
        "msid" = "dists <- TreeDist::MatchingSplitInfoDistance(trees)",
        "rf" = "dists <- TreeDist::RobinsonFoulds(trees)",
        "qd" = c("dists <- as.dist(Quartet::QuartetDivergence(",
                 "  Quartet::ManyToManyQuartetAgreement(trees),",
                 "  similarity = FALSE)", ")")
      ))
    }

    mapping <- bindCache(reactive({
      LogMsg("mapping()")
      if (maxProjDim() > 1L) {
        withProgress(
          message = "Mapping trees",
          value = 0.99,
          tryCatch(cmdscale(distances(), k = maxProjDim()),
                   warning = function(e) {
                     nDim <- as.integer(substr(e$message, 6, 7))
                     updateSliderInput(inputId = "spaceDim",
                                       value = min(nDim, input$spaceDim),
                                       max = nDim)
                     message("Max dimensions available for mapping: ", nDim, ".")
                     cmdscale(distances(), k = nDim)
                   })
        )
      } else {
        matrix(0, 0, 0)
      }
    }), r$treeHash, distMeth(), maxProjDim())

    LogMapping <- function() {
      k <- dim(mapping())[2]
      if (!is.null(k) && k > 0) {
        LogCommentP(paste0(
          "Generate first ", k, " dimensions of tree space using PCoA"
        ))
        LogCodeP(paste0("map <- cmdscale(dists, k = ", k, ")"))
      }
    }

    mstEnds <- bindCache(reactive({
      dist <- as.matrix(distances())
      withProgress(message = "Calculating MST", {
        edges <- MSTEdges(dist)
      })
      edges
    }), distMeth(), r$treeHash)

    ############################################################################
    # Tree space plot (from treespace.R)
    ############################################################################

    TreespacePlot <- function() {
      if (length(r$trees) < 3) {
        return(ErrorPlot("Need at least\nthree trees to\nmap tree space"))
      }

      cl <- clusterings()
      map <- mapping()

      nDim <- min(dims(), nProjDim())
      if (nDim < 2) {
        if (dim(map)[2] == 1L) {
          map <- cbind(map, 0)
        } else {
          map[, 2] <- 0
        }
        nDim <- 2L
        nPanels <- 1L
      } else {
        plotSeq <- matrix(0, nDim, nDim)
        nPanels <- nDim * (nDim - 1L) / 2L
        plotSeq[upper.tri(plotSeq)] <- seq_len(nPanels)
        if (nDim > 2) {
          plotSeq[nDim - 1, 2] <- max(plotSeq) + 1L
        }
        layout(t(plotSeq[-nDim, -1]))
      }

      par(mar = rep(0.2, 4))
      withProgress(message = "Drawing plot", {
        for (i in 2:nDim) for (j in seq_len(i - 1)) {
          incProgress(1 / nPanels)
          plot(map[, j], map[, i], ann = FALSE, axes = FALSE,
               frame.plot = nDim > 2L,
               type = "n", asp = 1, xlim = range(map), ylim = range(map))

          if ("seq" %in% input$mapLines) {
            lines(map[, j], map[, i], col = "#ffcc33", lty = 2)
          }

          if ("mst" %in% input$mapLines) {
            segments(map[mstEnds()[, 1], j], map[mstEnds()[, 1], i],
                     map[mstEnds()[, 2], j], map[mstEnds()[, 2], i],
                     col = "#bbbbbb", lty = 1)
          }

          points(map[, j], map[, i], pch = treePch(),
                 col = paste0(TreeCols(), as.hexmode(200)),
                 cex = spaceCex(),
                 lwd = spaceLwd()
          )

          if (cl$sil > silThreshold() && "hull" %in% input$mapLines) {
            for (clI in seq_len(cl$n)) {
              inCluster <- cl$cluster == clI
              clusterX <- map[inCluster, j]
              clusterY <- map[inCluster, i]
              hull <- chull(clusterX, clusterY)
              polygon(clusterX[hull], clusterY[hull], lty = 1, lwd = 2,
                      border = palettes[[min(length(palettes), cl$n)]][clI])
            }
          }
        }
        if (nDim > 2) {
          plot.new()
        }
        if (input$spacePch == "relat") {
          if (length(input$relators) == 4L) {
            legend(
              "topright",
              bty = "n",
              pch = 1:3,
              xpd = NA,
              pt.cex = spaceCex(),
              pt.lwd = spaceLwd(),
              gsub("_", " ", fixed = TRUE,
                   paste(input$relators[2:4], "&", input$relators[[1]]))
            )
          }
        } else if (input$spacePch == "name") {
          clstr <- treeNameClustering()
          clusters <- unique(clstr)
          if (length(clusters) > 1L) {
            legend(bty = "n", "topright", xpd = NA,
                   pch = c(1, 3, 4, 2,
                           seq_len(max(clstr))[-(1:4)])[clusters],
                   paste0("~ ", attr(clstr, "med"), " (", table(clstr), ")"))
          }
        }
        if (input$spaceCol == "firstHit" && length(FirstHit())) {
          legend(bty = "n", "topleft", pch = 16, col = FirstHitCols(),
                 pt.cex = spaceCex(),
                 names(FirstHit()), title = "Iteration first hit")
        } else if (input$spaceCol == "score") {
          legendRes <- length(badToGood)
          leg <- rep(NA, legendRes)
          leg[c(legendRes, 1)] <- signif(range(scores()))
          legend("bottomright", bty = "n", border = NA,
                 legend = leg, fill = rev(badToGood),
                 y.intersp = 0.04, cex = 1.1)
        }
      })
    }

    LogTreespacePlot <- function() {
      BeginLogP()

      LogClusterings()
      LogMapping()

      map <- mapping()
      nDim <- min(dims(), nProjDim())
      if (nDim < 2) {
        LogCommentP("Prepare 1D map", 0)
        if (dim(map)[2] == 1L) {
          LogCodeP("map <- cbind(map, 0)")
        } else {
          LogCodeP("map[, 2] <- 0")
        }
        nDim <- 2L
        nPanels <- 1L
      } else {
        LogCommentP("Prepare plot layout")

        LogCodeP(c(
          paste0("nDim <- ", nDim, " # Number of dimensions to plot"),
          "nPanels <- nDim * (nDim - 1L) / 2L # Lower-left triangle",
          "plotSeq <- matrix(0, nDim, nDim)",
          "plotSeq[upper.tri(plotSeq)] <- seq_len(nPanels)",
          if (nDim > 2) {
            "plotSeq[nDim - 1, 2] <- max(plotSeq) + 1L"
          },
          "layout(t(plotSeq[-nDim, -1]))"
        ))
      }

      LogCommentP("Set plot margins", 0)
      LogCodeP("par(mar = rep(0.2, 4))")

      LogCommentP("Set up tree plotting symbols")
      LogCodeP(paste0("treePch <- ", LogTreePch()),
               LogTreeCols(),
               "treeCols <- paste0(treeCols, as.hexmode(200)) # Semitransparent"
      )

      LogCodeP("for (i in 2:nDim) for (j in seq_len(i - 1)) {")
      LogIndent(+2)
      LogCommentP("Set up blank plot")
      LogCodeP("plot(",
               "  x = map[, j],",
               "  y = map[, i],",
               "  ann = FALSE,        # No annotations",
               "  axes = FALSE,       # No axes",
               paste0("  frame.plot = ",
                      if (nDim > 2L) {
                        "TRUE,  # Border around plot"
                      } else {
                        "FALSE, # No border around plot"
                      }),
               "  type = \"n\",         # Don't plot any points yet",
               "  asp = 1,            # Fix aspect ratio to avoid distortion",
               "  xlim = range(map),  # Constant X range for all dimensions",
               "  ylim = range(map)   # Constant Y range for all dimensions",
               ")")

      if ("seq" %in% input$mapLines) {
        LogCommentP("Connect trees in sequence")
        LogCodeP("lines(",
                 "  x = map[, j],",
                 "  y = map[, i],",
                 "  col = \"#ffcc33\", # Orange",
                 "  lty = 2 # dashed",
                 ")")
      }

      if ("mst" %in% input$mapLines) {
        LogCommentP("Plot minimum spanning tree (Gower 1969)")
        LogCodeP(
          "mst <- MSTEdges(as.matrix(dists))",
          "segments(",
          "  x0 = map[mst[, 1], j],",
          "  y0 = map[mst[, 1], i],",
          "  x1 = map[mst[, 2], j],",
          "  y1 = map[mst[, 2], i],",
          "  col = \"#bbbbbb\", # Light grey",
          "  lty = 1          # Solid lines",
          ")"
        )
      }

      LogCommentP("Add points")
      LogCodeP(
        "points(",
        "  x = map[, j],",
        "  y = map[, i],",
        "  pch = treePch,",
        "  col = treeCols,",
        paste0("  cex = ", spaceCex(), ", # Point size"),
        paste0("  lwd = ", spaceLwd(), " # Line width"),
        ")"
      )

      cl <- clusterings()
      if (cl$sil > silThreshold() && "hull" %in% input$mapLines) {
        LogCommentP("Mark clusters")
        LogCodeP("for (clI in seq_len(nClusters)) {")
        LogIndent(+2)
        LogCodeP(
          "inCluster <- clustering == clI",
          "clusterX <- map[inCluster, j]",
          "clusterY <- map[inCluster, i]",
          "hull <- chull(clusterX, clusterY)",
          "polygon(",
          "  x = clusterX[hull],",
          "  y = clusterY[hull],",
          "  lty = 1, # Solid line style",
          "  lwd = 2, # Wider line width",
          "  border = clusterCol[clI]",
          ")")
        LogIndent(-2)
        LogCodeP("}")
      }

      LogIndent(-2)
      LogCodeP("}")

      if (nDim > 2) {
        LogCodeP("plot.new() # Use new panel to plot legends")
      }

      if (input$spacePch == "relat") {
        if (length(input$relators) == 4L) {
          LogCommentP("Add legend for plotting symbols")
          LogCodeP(
            "legend(",
            "  \"topright\",",
            "  bty = \"n\", # No legend border box",
            "  pch = 1:3, # Legend symbols",
            "  xpd = NA, # Display overflowing text",
            paste0("  pt.cex = ", spaceCex(), ", # Point size"),
            paste0("  pt.lwd = ", spaceLwd(), ", # Line width"),
            paste0("  ",
                   EnC(gsub("_", " ", fixed = TRUE,
                            paste(input$relators[2:4], "&",
                                  input$relators[[1]])))
            ), ")"
          )
        }
      } else if (input$spacePch == "name") {
        clstr <- treeNameClustering()
        clusters <- unique(clstr)
        if (length(clusters) > 1L) {
          LogCommentP("Add legend for plotting symbols")
          LogCodeP(
            "nameClusters <- ClusterStrings(names(trees))",
            "uniqueClusters <- unique(nameClusters)",
            "legend(",
            "  \"topright\",",
            "  bty = \"n\", # No legend border box",
            "  xpd = NA, # Display overflowing text",
            paste0(
              "  pch = ",
              EnC(c(1, 3, 4, 2,
                    seq_len(max(clstr))[-(1:4)])[clusters]),
              ", # Legend symbols"
            ), paste0("  ",
                      EnC(paste0("~ ", attr(clstr, "med"),
                                 " (", table(clstr), ")"))
            ),
            ")")
        }
      }
      if (input$spaceCol == "firstHit" && length(FirstHit())) {
        LogCommentP("Record when trees first hit")
        LogFirstHit()

        LogCommentP("Add legend for symbol colours")
        LogCodeP(
          "legend(",
          "  \"topleft\",",
          "  bty = \"n\", # No legend border box",
          "  pch = 16, # Circle symbol",
          "  xpd = NA, # Display overflowing text",
          paste0("  col = ", LogFirstHitCols(), ","),
          paste0("  pt.cex = ", spaceCex(), ", # Point size"),
          paste0("  ", EnC(names(FirstHit())), ","),
          "  title = \"Iteration first hit\"",
          ")"
        )
      } else if (input$spaceCol == "score") {
        LogCommentP("Add legend for symbol colours")
        LogCodeP(
          "goodToBad <- hcl.colors(108, \"Temps\")",
          "leg <- rep_len(NA, 108)",
          paste0("leg[c(1, 108)] <- ",
                 EnC(rev(signif(range(scores()))))),
          "legend(",
          "  \"bottomright\",",
          "  legend = leg,",
          "  bty = \"n\", # No legend border box",
          "  border = NA, # No border around plot icons",
          "  xpd = NA, # Display overflowing text",
          "  fill = goodToBad,",
          "  y.intersp = 0.04, # Compress squares to make gradient scale",
          "  cex = 1.1 # Increase font and icon size slightly",
          ")"
        )
      }
    }

    ############################################################################
    # saveDetails (shared with downloads)
    ############################################################################

    saveDetails <- reactive({
      switch(plotFormat(),
             "cons" = list(
               fileName = "ConsensusTrees",
               title = "Consensus tree - TreeSearch",
               asp = 2L
             ),
             "clus" = list(
               fileName = "ClusterCons",
               title = "Cluster Consensus trees - TreeSearch",
               asp = 1.6
             ),
             "ind" = list(
               fileName = "OptimalTree",
               title = "Optimal tree - TreeSearch",
               asp = 2L
             ),
             "space" = list(
               fileName = "TreeSpace",
               title = "Tree space - TreeSearch",
               asp = 1L
             ))
    })

    ############################################################################
    # Mapping quality (moved from consensus.R)
    ############################################################################

    LogScore <- function(x) {
      (-(log10(1 - pmin(1, x) + 1e-2))) / 2
    }

    QualityPlot <- function(quality) {
      par(mar = c(2, 0, 0, 0))
      nStop <- length(badToGood) + 1L

      plot(NULL, xlim = c(0, 1), ylim = c(-1.5, 2.5),
           ann = FALSE, axes = FALSE)
      x <- seq.int(from = 0, to = 1, length.out = nStop)
      segments(x[-nStop], numeric(nStop), x[-1], lwd = 5, col = badToGood)

      trust <- quality[["Trustworthiness"]]
      cont  <- quality[["Continuity"]]
      txc   <- quality[["sqrtTxC"]]

      if (trust > 1) {
        LogMsg("Preternaturally high Trustworthiness: ", trust)
      }
      if (cont > 1) {
        LogMsg("Preternaturally high Continuity: ", cont)
      }
      LogMsg(trust * nStop)
      segments(LogScore(txc), -1, y1 = 1, lty = 3)
      text(LogScore(trust), 1, "T",
           col = badToGood[LogScore(trust) * nStop])
      text(LogScore(cont), -1, "C",
           col = badToGood[LogScore(cont) * nStop])

      tickPos <- c(0, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0)
      ticks <- LogScore(tickPos)

      axis(1, at = ticks, labels = NA, line = 0)
      axis(1, tick = FALSE, at = ticks, labels = tickPos, line = 0)
      axis(1, line = -1, tick = FALSE,
           at = ticks[-1] - ((ticks[-1] - ticks[-length(ticks)]) / 2),
           labels = c("", "dire", "", "ok", "gd", "excellent"))
      axis(3, at = 0.5, tick = FALSE, line = -2,
           paste0(dims(), "D mapping quality (trustw. / contin.):"))
    }

    output$pcQuality <- renderCachedPlot({
      if (length(r$trees) < 3) {
        return()
      }
      dstnc <- distances()
      mppng <- mapping()
      mppng <- mapping()[, seq_len(min(dim(mppng)[2], dims()))]
      neighbs <- min(10L, length(r$trees) / 2)
      future_promise(
        TreeDist::MappingQuality(dstnc, dist(mppng), neighbs),
        seed = TRUE) %...>% QualityPlot
    }, cacheKeyExpr = {
      list(r$treeHash, distMeth(), dims())
    },
      sizePolicy = function(dims) dims
    )

    ############################################################################
    # Return reactives for other modules
    ############################################################################

    list(
      distances        = distances,
      mapping          = mapping,
      dims             = dims,
      nProjDim         = nProjDim,
      TreeCols         = TreeCols,
      treePch          = treePch,
      mstEnds          = mstEnds,
      saveDetails      = saveDetails,
      TreespacePlot    = TreespacePlot,
      LogTreespacePlot = LogTreespacePlot,
      LogDistances     = LogDistances,
      # Expose input values for cache keys in consensus.R
      spaceCol  = reactive(input$spaceCol),
      spacePch  = reactive(input$spacePch),
      mapLines  = reactive(input$mapLines),
      relators  = reactive(input$relators)
    )
  })
}
