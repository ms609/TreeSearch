#' Cluster similar strings
#' 
#' Calculate string similarity using the Levenshtein distance and return
#' clusters of similar strings. 
#' 
#' @param x Character vector.
#' @param maxCluster Integer specifying maximum number of clusters to consider.
#' @return `NameClusters()` returns an integer assigning each element of `x`
#' to a cluster, with an attribute `med` specifying the median string in each 
#' cluster, and `silhouette` reporting the silhouette coefficient of the optimal
#' clustering.  Coefficients < 0.5 indicate weak structure, and no clusters are
#' returned.  If the number of unique elements of `x` is less than `maxCluster`,
#' all occurrences of each entry are assigned to an individual cluster.
#' 
#' @examples
#' ClusterStrings(c(paste0("FirstCluster ", 1:5),
#'                  paste0("SecondCluster.", 8:12),
#'                  paste0("AnotherCluster_", letters[1:6])))
#' @template MRS
#' @importFrom utils adist
#' @importFrom cluster pam silhouette
#' @importFrom protoclust protoclust
#' @importFrom stats as.dist cutree
#' @family utility functions
#' @export
ClusterStrings <- function (x, maxCluster = 12) {
  if (maxCluster < 2L) {
    stop("`maxCluster` must be at least two.")
  }
  
  if (length(unique(x)) < maxCluster) {
    nom <- unique(x)
    structure(match(x, nom), "med" = nom)
  } else {
    possibleClusters <- 2:maxCluster
    hSil <- pamSil <- -99
    dists <- adist(x) # approximate string distance
    
    nMethodsChecked <- 2
    methInc <- 1 / nMethodsChecked
    nK <- length(possibleClusters)
    kInc <- 1 / (nMethodsChecked * nK)
    
    # pamClusters <- lapply(possibleClusters, function (k) {
    #   pam(dists, k = k)
    # })
    # pamSils <- vapply(pamClusters, function (pamCluster) {
    #   mean(silhouette(pamCluster)[, 3])
    # }, double(1))
    # bestPam <- which.max(pamSils)
    # pamSil <- pamSils[bestPam]
    # pamCluster <- pamClusters[[bestPam]]$cluster
    
    hTree <- protoclust(as.dist(dists))
    hClusters <- lapply(possibleClusters, function (k) cutree(hTree, k = k))
    hSils <- vapply(hClusters, function (hCluster) {
      mean(silhouette(hCluster, dists)[, 3])
    }, double(1))
    bestH <- which.max(hSils)
    hSil <- hSils[bestH]
    hCluster <- hClusters[[bestH]]
    pamSil <- hSil
    pamCluster <- hCluster
    
    bestCluster <- c("none", "pam", "hmm")[which.max(c(0.5, pamSil, hSil))]
    
    clustering <- switch(bestCluster, pam = pamCluster, hmm = hCluster, 1)
    
    medians <- vapply(seq_len(max(clustering)),
                      function (i) {
                        these <- clustering == i
                        x[these][which.min(colSums(dists[these, these]))]
                      }, character(1))
    
    structure(clustering,
              silhouette = switch(bestCluster, pam = pamSil, hmm = hSil,
                                  max(pamSil, hSil)),
              med = medians)
  }
}

