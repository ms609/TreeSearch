#' Tree search using successive approximations
#'
#' Searches for a tree that is optimal under the Successive Approximations 
#' criterion \insertCite{Farris1969}{TreeSearch}.
#'
#' @template treeParam
#' @template datasetParam
#' @param outgroup if not NULL, taxa on which the tree should be rooted
#' @param k Constant for successive approximations, see Farris 1969 p. 379
#' @param maxSuccIter maximum iterations of successive approximation
#' @param ratchetHits maximum hits for parsimony ratchet 
#' @param searchHits maximum hits in tree search
#' @param searchIter maximum iterations in tree search
#' @param ratchetIter maximum iterations of parsimony ratchet
#' @template verbosityParam
#' @param suboptimal retain trees that are this proportion less optimal than the optimal tree
#' 
#' @return `SuccessiveApproximations()` returns a list of class `multiPhylo`
#' containing optimal (and slightly suboptimal, if suboptimal > 0) trees.
#'
#' @references 
#' \insertAllCited{}
#'
#' @importFrom ape consensus root
#' @family custom search functions
#' @export
SuccessiveApproximations <- function (tree, dataset, outgroup = NULL, k = 3,
                                      maxSuccIter = 20, ratchetHits = 100,
                                      searchHits = 50, searchIter = 500,
                                      ratchetIter = 5000, verbosity = 0,
                                      suboptimal = 0.1) {
  
  if (k < 1) stop ("k should be at least 1, see Farris 1969 p.379")
  attr(dataset, "sa.weights") <- rep.int(1, length(attr(dataset, "weight")))
  collectSuboptimal <- suboptimal > 0
  
  max.node <- max(tree[["edge"]][, 1])
  n.tip <- length(tree[["tip.label"]])
  n.node <- max.node - n.tip
  bests <- vector("list", maxSuccIter + 1L)
  bestsConsensus <- vector("list", maxSuccIter + 1L)
  best <- bests[[1]] <- bestsConsensus[[1]] <- root(tree, outgroup, resolve.root=TRUE)
  for (i in seq_len(maxSuccIter) + 1L) {
    if (verbosity > 0) message("\nSuccessive Approximations Iteration ", i - 1L)
    attr(best, "score") <- NULL
    if (suboptimal > 0) {
      suboptimalSearch <- suboptimal * sum(attr(dataset, "sa.weights") *
                                             attr(dataset, "weight"))
    }
    trees <- Ratchet(best, dataset, TreeScorer = SuccessiveWeights,
                     all = collectSuboptimal, 
                     suboptimal = suboptimalSearch,
                     rearrangements = "NNI",
                     ratchetHits=ratchetHits, searchHits = searchHits,
                     searchIter = searchIter, ratchetIter = ratchetIter,
                     outgroup = outgroup, verbosity = verbosity - 1)
    trees <- unique(trees)
    bests[[i]] <- trees
    suboptimality <- Suboptimality(trees)
    bestsConsensus[[i]] <- consensus(trees[suboptimality == 0])
    if (all.equal(bestsConsensus[[i]], bestsConsensus[[i - 1]])) {
      return(bests[2:i])
    }
    best <- trees[suboptimality == 0][[1]]
    l.i <- CharacterLength(best, dataset, compress = TRUE)
    p.i <- l.i / (n.node - 1)
    w.i <- ((p.i)^-k) - 1
    attr(dataset, "sa.weights") <- w.i
  }
  message("Stability not reached.")
  
  # Return:
  structure(bests, class = "multiPhylo")
}

#' Tree suboptimality
#'
#' How suboptimal is a tree?
#'
#' @param trees list of trees, to include an optimal tree
#' @param proportional logical stating whether to normalise results to lowest
#' score
#' @return `Suboptimality()` returns a vector listing, for each tree, how much
#' its score differs from the optimal (lowest) score.
#' @keywords internal
#' @export
Suboptimality <- function (trees, proportional = FALSE) {
  scores <- vapply(trees, attr, double(1), "score")
  
  # Return:
  if (proportional) {
    (scores - min(scores)) / min(scores)
  } else {
    scores - min(scores)
  }
}

#' @rdname SuccessiveApproximations
#' @return `SuccessiveWeights()` returns the score of a tree, given the
#' weighting instructions specified in the attributes of the dataset.
#' 
#' @keywords internal
#' @export
SuccessiveWeights <- function(tree, dataset) {
  # Data
  if (inherits(dataset, "phyDat")) dataset <- PrepareDataSA(dataset)
  if (!inherits(dataset, "saDat")) {
    stop("Invalid data type; prepare data with PhyDat() or PrepareDataSA().")
  }
  at <- attributes(dataset)
  weight <- at[["weight"]]
  sa.weights <- at[["sa.weights"]]
  if (is.null(sa.weights)) sa.weights <- rep.int(1, length(weight))
  steps <- CharacterLength(tree, dataset, compress = TRUE)
  
  # Return:
  sum(steps * sa.weights * weight)
}

PrepareDataSA <- function (dataset) {
# Written with reference to phangorn:::prepareDataFitch
  at <- attributes(dataset)
  nam <- at[["names"]]
  nLevel <- length(at[["levels"]])
  nChar <- at[["nr"]]
  cont <- attr(dataset, "contrast")
  nTip <- length(dataset)
  
  at[["names"]] <- NULL
  powers.of.2 <- 2L ^ c(0L:(nLevel - 1L))
  tmp <- cont %*% powers.of.2
  tmp <- as.integer(tmp)
  dataset <- unlist(dataset, FALSE, FALSE)
  ret <- tmp[dataset] 
  ret <- as.integer(ret)
  attributes(ret) <- at
  inappLevel <- which(at[["levels"]] == "-")
  attr(ret, "dim") <- c(nChar, nTip)  
  applicableTokens <- setdiff(powers.of.2, 2 ^ (inappLevel - 1))
  attr(ret, "split.sizes") <- t(apply(ret, 1, function(x) vapply(applicableTokens, function (y) sum(x == y), integer(1))))
  dimnames(ret) <- list(NULL, nam)
  attr(ret, "bootstrap") <- list("split.sizes", "sa.weights")
  class(ret) <- "saDat"
  ret
}
