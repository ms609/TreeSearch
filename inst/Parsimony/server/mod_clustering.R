# Module: Clustering analysis
#
# Owns inputs: clThresh. Owns distances computation (shared with treespace).
# Reads: r$trees, r$treeHash.
# Receives top-level distMeth as reactive arg.
#
# Returns a list of reactives:
#   distances, LogDistances, silThreshold, clusterings, LogClusterings

clustering_ui <- function(id) {
  ns <- NS(id)
  sliderInput(ns("clThresh"), "Cluster threshold:", value = 0.5,
              min = 0, max = 1, width = 200)
}

#' @param id Module namespace id.
#' @param r AppState reactiveValues.
#' @param distMeth Reactive wrapping top-level \code{input$distMeth}.
#' @param log_fns Named list of logging functions from logging.R:
#'   LogMsg, LogCommentP, LogCodeP, LogIndent, BeginLogP, LogExprP.
clustering_server <- function(id, r, distMeth, log_fns) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Unpack logging functions
    LogMsg      <- log_fns$LogMsg
    LogCommentP <- log_fns$LogCommentP
    LogCodeP    <- log_fns$LogCodeP
    LogIndent   <- log_fns$LogIndent
    BeginLogP   <- log_fns$BeginLogP
    LogExprP    <- log_fns$LogExprP

    ############################################################################
    # Silhouette threshold (debounced clThresh input)
    ############################################################################

    silThreshold <- debounce(reactive({
      input$clThresh
    }), 50)

    ############################################################################
    # Tree distances (moved from treespace module)
    ############################################################################

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

    ############################################################################
    # Clusterings
    ############################################################################

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
                            function(k) TreeDist::KMeansPP(dists, k))
        kSils <- vapply(kClusters, function(kCluster) {
          mean(cluster::silhouette(kCluster$cluster, dists)[, 3])
        }, double(1))
        bestK <- which.max(kSils)
        kSil <- kSils[bestK]
        kCluster <- kClusters[[bestK]]$cluster

        cli::cli_progress_update(1, status = "PAM")
        pamClusters <- lapply(possibleClusters, function(k) {
          cluster::pam(dists, k = k)
        })
        pamSils <- vapply(pamClusters, function(pamCluster) {
          mean(cluster::silhouette(pamCluster)[, 3])
        }, double(1))
        bestPam <- which.max(pamSils)
        pamSil <- pamSils[bestPam]
        pamCluster <- pamClusters[[bestPam]]$cluster

        cli::cli_progress_update(1, status = "Hierarchical")
        hTree <- protoclust::protoclust(dists)
        hClusters <- lapply(possibleClusters, function(k) cutree(hTree, k = k))
        hSils <- vapply(hClusters, function(hCluster) {
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
             "; sil: ", signif(switch(bestCluster,
                                      pam = pamSil, hmm = hSil,
                                      kmn = kSil, 0)))
      # Return:
      list(method = switch(bestCluster,
                           pam = "part. around medoids",
                           hmm = "minimax linkage",
                           kmn = "k-means",
                           none = "no significant clustering"),
           n = 1 + switch(bestCluster, pam = bestPam, hmm = bestH,
                          kmn = bestK, 0),
           sil = switch(bestCluster, pam = pamSil, hmm = hSil,
                        kmn = kSil, 0),
           cluster = switch(bestCluster, pam = pamCluster, hmm = hCluster,
                            kmn = kCluster, 1)
      )

    }), r$treeHash, silThreshold(), distMeth())

    ############################################################################
    # LogClusterings
    ############################################################################

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
        LogCommentP(
          "no support < 0.25 < weak < 0.5 < good < 0.7 < strong", 0)
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
          LogCommentP(paste0("Best clustering was ",
                            clusterings()$method, ":"))
          LogCommentP(paste0("Silhouette coefficient = ",
                            signif(clusterings()$sil)), 0)
          LogCommentP(paste0("Store the cluster to which each tree is ",
                            "optimally assigned:"))
          LogCodeP(paste0(
            "clustering <- switch(bestCluster, pam = pamCluster,",
            " hmm = hCluster, kmn = kCluster, 1)"),
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

    ############################################################################
    # clThresh label CSS class (color-codes threshold strength)
    ############################################################################

    observeEvent(input$clThresh, {
      classes <- c("meaningless", "weak", "good", "strong")
      liveClass <- classes[as.integer(cut(
        input$clThresh, c(0, 0.25, 0.5, 0.7, 1),
        include.lowest = TRUE, right = FALSE
      ))]
      labelId <- ns("clThresh-label")
      runjs(paste0(
        "$('#", labelId, "').removeClass('", paste(classes, collapse = " "),
        "').addClass('", liveClass, "');"
      ))
    })

    ############################################################################
    # Return reactives for server.R and other modules
    ############################################################################

    list(
      distances      = distances,
      LogDistances   = LogDistances,
      silThreshold   = silThreshold,
      clusterings    = clusterings,
      LogClusterings = LogClusterings
    )
  })
}
