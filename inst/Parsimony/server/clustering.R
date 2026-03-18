  # Cluster consensus plotting functions
  #
  # These depend heavily on consensus.R functions (KeptTips, ConsensusWithout,
  # UserRoot, unitEdge, consP, TipCols, LabelConcordance, PutTree, PutData)
  # and will migrate to mod_consensus.R (T-063). For now they stay source'd.
  #
  # Uses clusterings() and silThreshold() exposed from the clustering module.

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

      for (i in seq_len(cl$n)) {
        col <- palettes[[min(length(palettes), cl$n)]][i]
        PutTree(r$trees)
        PutData(cl$cluster)

        cons <- ConsensusWithout(r$trees[cl$cluster == i], dropped, p = consP())
        cons <- UserRoot(cons)
        if (unitEdge()) {
          cons$edge.length <- rep.int(1, dim(cons$edge)[1])
        }
        cons <- SortEdges(cons)
        r$plottedTree[[i]] <- cons
        plot(cons, edge.width = 2, font = 3, cex = 0.83,
             edge.color = col, tip.color = TipCols()[cons$tip.label])
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
           edge.color = palettes[[1]], tip.color = TipCols()[cons$tip.label])
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
      LogCommentP("Plot each consensus tree in turn:", 1)
      LogCodeP(paste0("for (i in seq_len(", cl$n, ")) {"))
      LogIndent(+2)
      LogCodeP(
        "clusterTrees <- trees[clustering == i]",
        "cons <- ConsensusWithout(",
        "  trees = clusterTrees,",
        paste0("  tip = ", EnC(dropped), ","),
        paste0("  p = ", consP()),
        ")"
      )
      LogUserRoot(dropped = dropped)
      if (unitEdge()) {
        LogExprP("cons$edge.length <- rep.int(1, nrow(cons$edge))")
      }
      LogSortEdges("cons")
      LogCodeP("plot(",
              "  cons,",
              "  edge.width = 2,             # Widen lines",
              "  font = 3,                   # Italicize labels",
              "  cex = 0.83,                 # Shrink tip font size",
              "  edge.color = clusterCol[i], # Colour tree",
              "  tip.color = tipCols[cons$tip.label]",
              ")")
      LogCodeP("legend(",
              "  \"bottomright\",",
              "  paste(\"Cluster\", i),",
              "  pch = 15,            # Filled circle icon",
              "  pt.cex = 1.5,        # Increase icon size",
              "  col = clusterCol[i],",
              "  bty = \"n\"            # Don't plot legend in box",
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
