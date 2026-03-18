  silThreshold <- debounce(reactive({
    input$clThresh
  }), 50)
  
  ##############################################################################
  # Clusterings
  ##############################################################################
  clusterings <- bindCache(reactive({
    ## CAUTION: Update LogClusterings() to reflect any changes made
    ## to this function 
    LogMsg("clusterings()")
    maxCluster <- min(15L, length(r$trees) - 1L)
    if (maxCluster > 1L) {
      possibleClusters <- 2:maxCluster
      
      hSil <- pamSil <- -99
      dists <- distances()
      
      nMethodsChecked <- 3L
      cli::cli_progress_bar("Computing clusterings", "K-means",
                            total = nMethodsChecked)
      
      nK <- length(possibleClusters)
    
      kClusters <- lapply(possibleClusters,
                          function (k) TreeDist::KMeansPP(dists, k))
      kSils <- vapply(kClusters, function (kCluster) {
        mean(cluster::silhouette(kCluster$cluster, dists)[, 3])
      }, double(1))
      bestK <- which.max(kSils)
      kSil <- kSils[bestK]
      kCluster <- kClusters[[bestK]]$cluster
      
      cli::cli_progress_update(1, status = "PAM")
      pamClusters <- lapply(possibleClusters, function (k) {
        cluster::pam(dists, k = k)
      })
      pamSils <- vapply(pamClusters, function (pamCluster) {
        mean(cluster::silhouette(pamCluster)[, 3])
      }, double(1))
      bestPam <- which.max(pamSils)
      pamSil <- pamSils[bestPam]
      pamCluster <- pamClusters[[bestPam]]$cluster
      
      cli::cli_progress_update(1, status = "Hierarchical")
      hTree <- protoclust::protoclust(dists)
      hClusters <- lapply(possibleClusters, function (k) cutree(hTree, k = k))
      hSils <- vapply(hClusters, function (hCluster) {
        mean(cluster::silhouette(hCluster, dists)[, 3])
      }, double(1))
      bestH <- which.max(hSils)
      hSil <- hSils[bestH]
      hCluster <- hClusters[[bestH]]
      cli::cli_progress_update(1, status = "Done")
      
      bestCluster <- c("none", "pam", "hmm", "kmn")[
        which.max(c(silThreshold(), pamSil, hSil, kSil))]
    } else {
      bestCluster <- "none"
    }
     
    LogMsg("Best clustering: ", bestCluster, 
        "; sil: ", signif(switch(bestCluster, pam = pamSil, hmm = hSil, kmn = kSil, 0)))
    # Return:
    list(method = switch(bestCluster, pam = "part. around medoids",
                                      hmm = "minimax linkage",
                                      kmn = "k-means",
                                      none = "no significant clustering"),
         n = 1 + switch(bestCluster, pam = bestPam, hmm = bestH, kmn = bestK, 0),
         sil = switch(bestCluster, pam = pamSil, hmm = hSil, kmn = kSil, 0), 
         cluster = switch(bestCluster, pam = pamCluster, hmm = hCluster, kmn = kCluster, 1)
    )

  }), r$treeHash, silThreshold(), input$distMeth)
  
  LogClusterings <- function() {
    maxCluster <- min(15L, length(r$trees) - 1L)
    if (maxCluster > 1L) {
      possibleClusters <- paste(2, maxCluster, sep = ":")
      
      hSil <- pamSil <- -99
      LogDistances()
      dists <- distances()
      
      LogCommentP("Compute clusters of trees", 2)
      nK <- length(possibleClusters)
      LogCommentP("Try K-means++ clustering (Arthur & Vassilvitskii 2007):")
      LogCodeP(
        paste0(
          "kClusters <- lapply(", possibleClusters, ", ",
          "function (k) KMeansPP(dists, k)", ")"
        ),
        "kSils <- vapply(kClusters, function (kCluster) {",
        "  mean(cluster::silhouette(kCluster$cluster, dists)[, 3])",
        "}, double(1))",
        "bestK <- which.max(kSils)",
        "kSil <- kSils[bestK] # Best silhouette coefficient",
        "kCluster <- kClusters[[bestK]]$cluster # Best solution"
      )
      
      LogCommentP("Try partitioning around medoids (Maechler et al. 2019):")
      LogCodeP(
        paste0(
          "pamClusters <- lapply(", possibleClusters, ", ",
          "function (k) cluster::pam(dists, k = k)", ")"
        ),
        "pamSils <- vapply(pamClusters, function (pamCluster) {",
        "  mean(cluster::silhouette(pamCluster)[, 3])",
        "}, double(1))",
        "bestPam <- which.max(pamSils)",
        "pamSil <- pamSils[bestPam] # Best silhouette coefficient",
        "pamCluster <- pamClusters[[bestPam]]$cluster # Best solution"
      )
      
      
      LogCommentP(
        paste("Try hierarchical clustering with minimax linkage",
              "(Bien & Tibshirani 2011):")
      )
      LogCodeP(
        "hTree <- protoclust::protoclust(dists)",
        paste0(
          "hClusters <- lapply(", possibleClusters, ", ", 
          "function (k) cutree(hTree, k = k)", ")"
        ),
        "hSils <- vapply(hClusters, function (hCluster) {",
        "  mean(cluster::silhouette(hCluster, dists)[, 3])",
        "}, double(1))",
        "bestH <- which.max(hSils)",
        "hSil <- hSils[bestH] # Best silhouette coefficient",
        "hCluster <- hClusters[[bestH]] # Best solution"
      )
      
      LogCommentP("Set threshold for recognizing meaningful clustering")
      LogCommentP("no support < 0.25 < weak < 0.5 < good < 0.7 < strong", 0)
      LogCodeP(paste0("threshold <- ", silThreshold()))
      
      LogCommentP("Compare silhouette coefficients of each method")
      LogCodeP(
        "bestMethodId <- which.max(c(threshold, pamSil, hSil, kSil))",
        "bestCluster <- c(\"none\", \"pam\", \"hmm\", \"kmn\")[bestMethodId]"
      )
      if (clusterings()$n == 1) {
        LogCommentP("No significant clustering was found.")
        LogCodeP("clustering <- 1 # Assign all trees to single cluster")
      } else {
        LogCommentP(paste0("Best clustering was ", clusterings()$method, ":"))
        LogCommentP(paste0("Silhouette coefficient = ",
                          signif(clusterings()$sil)), 0)
        LogCommentP(paste0("Store the cluster to which each tree is ",
                          "optimally assigned:"))
        LogCodeP(paste0(
          "clustering <- switch(bestCluster, pam = pamCluster, hmm = hCluster,",
          " kmn = kCluster, 1)"),
          paste0("nClusters <- length(unique(clustering))"),
          paste0(
          "clusterCol <- ",
          EnC(palettes[[min(length(palettes), clusterings()$n)]]),
          " # Arbitrarily"
          )
        )
      }
    } else {
      LogCommentP("Not enough trees for clustering analysis")
      LogCodeP("bestCluster <- \"none\"")
      LogCodeP("nClusters <- 1")
    }
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
  
