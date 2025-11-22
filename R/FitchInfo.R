#' Information required to encode a character using a Fitch-like algorithm
#' 
#' Given a tree, the amount of information required to encode a character is
#' determined by the number of extra steps the character implies on a tree.
#' 
#' If there are \eqn{n} ways to produce the observed distribution of tokens with
#' \eqn{s} extra steps, then the minimum encoding length of the observed
#' distribution (given the number of states observed, the number of steps, and
#' the tree) is \eqn{n\log{n}}.
#' 
#' This encoding captures how closely a character constrains a tree, which is
#' not quite what we are interested in: there are many more trees on which
#' a character fits rather poorly than on which it fits maximally poorly.
#' 
#' Hence we use a cumulative information measure, asking instead how many trees
#' require _at most_ as many steps as the number that are observed.
#' This fits with the possibility that a parsimonious distribution may arise
#' through multiple steps, even if such steps would not be distinguished
#' by a parsimonious reconstruction.
#' 
#' This value is most readily interpreted when normalized against the expected
#' entropy of a random shuffling of states.
#' 
#' 
#' @section Quirks:
#' Characters that are not parsimony informative (i.e. singletons) contain no
#' phylogenetic information, so contribute zero to both numerator
#' and denominator.  If no character is parsimony informative, we define the
#' return value as 1.
#' 
#' @examples
#' dataset <- inapplicable.phyData[["Vinther2008"]]
#' tree <- TreeTools::NJTree(dataset)
#' FitchInfo(tree, dataset)
#' @importFrom TreeTools LnUnrooted
#' @importFrom TreeDist Entropy
#' @template MRS
#' @export
FitchInfo <- function(tree, dataset) {
  # Check inputs
  if (is.null(dataset)) {
    warning("Cannot calculate concordance without `dataset`.")
    return(NULL)
  }
  if (is.null(tree)) {
    warning("Cannot calculate concordance without `dataset`.")
    return(NULL)
  }
  
  keep <- MatchStrings(TipLabels(tree), names(dataset), warning)
  if (length(keep) == 0) {
    return(NULL)
  }
  dataset <- dataset[keep]
  
  at <- attributes(dataset)
  cont <- at[["contrast"]]
  if ("-" %in% colnames(cont)) {
    cont[cont[, "-"] > 0, ] <- 1
  }
  ambiguous <- rowSums(cont) != 1
  
  mat <- matrix(as.integer(unlist(dataset)), length(dataset), byrow = TRUE)
  mat[mat %in% which(ambiguous)] <- NA_integer_
  maxToken <- max(mat, na.rm = TRUE)
  tokens <- as.character(seq_len(maxToken))
  mat <- apply(mat, 2, function (x) {
    uniques <- tabulate(x, maxToken) == 1
    x[x %in% tokens[uniques]] <- NA_integer_
    tokens <- sort(unique(x[!is.na(x)]))
    newLabel <- tabulate(tokens)
    newLabel[newLabel > 0] <- seq_len(sum(newLabel > 0)) - 1
    newLabel[x]
  }) # retains phyDat character compression
  
  .TwoStateH <- function(char, tree) {
    # Prune tree to fit, and process
    tr <- KeepTip(tree, !is.na(char))
    char <- char[!is.na(char)]
    tab <- table(char, deparse.level = 0)
    if (length(tab) != 2) stop("More than two states")
    logP <- vapply(seq_len(min(tab)), LogCarter1, 1.0, tab[[1]], tab[[2]]) - 
      LnUnrooted(sum(tab))
    cumH <- .LogCumSumExp(logP) / -log(2)
    cumH[abs(cumH) < sqrt(.Machine$double.eps)] <- 0
    steps <- CharacterLength(tr, StringToPhyDat(paste0(char), tr))
    h <- cumH[[steps]]
    expH <- sum(cumH * exp(logP))
    c(norm = (h - expH) / (cumH[[1]] - expH),
      h = h,
      hMax = cumH[[1]],
      expH = expH)
  }
  
  charH <- apply(mat, 2, function(char) {
    obs <- table(char, deparse.level = 0)
    pursue <- names(obs[obs > 1])
    if (length(pursue) < 2) {
      c(h = 0, hMax = 0)
    } else {
      # Prune tree to fit, and process
      tr <- KeepTip(tree, !is.na(char))
      char <- char[!is.na(char)]
      tab <- table(char, deparse.level = 0)
      steps <- CharacterLength(tr, StringToPhyDat(paste0(char), tr))
      
      nLevels <- length(tab)
      logP <- if (nLevels > 2) {
        bTab <- double(2 ^ nLevels - 1)
        bTab[2 ^ (seq_along(tab) - 1)] <- tab
        MaddisonSlatkin((nLevels - 1):steps, bTab)
      } else {
        LogCarter1(seq_len(steps), rep(tab[[1]], steps), rep(tab[[2]], steps)) -
          LnUnrooted(sum(tab))
      }
      cumH <- LogCumSumExp(logP) / -log(2)
      h <- cumH[[steps - (nLevels - 2)]]
      expH <- sum(cumH * exp(logP))
      c(norm = (h - expH) / (cumH[[1]] - expH),
        h = h,
        hMax = cumH[[1]],
        expH = expH)
    }
  })[, attr(dataset, "index"), drop = FALSE]
  
  totalHMax <- sum(charH["hMax", ])
  # Return:
  structure(
    if(totalHMax == 0) 1 else
      sum(charH["h", , drop = FALSE] * charH["hMax", , drop = FALSE] / totalHMax),
    byChar = charH
    )
}

#' Here we take an approach that is perhaps closer to our original intent.
#' 
#' Here we take as input a tree topology, a distribution of states, and a
#' number of steps.
#' 
#' We ask "how many ways are there to produce the observed state frequencies by
#' placing at most the observed number of steps on the observed tree".
#' 
#' Well, actually, we are asking "what is the probability of seeing at most
#' the observed number of steps, given the entropy of the observed state
#' frequencies".
#' 
#' @examples
#' tree <- as.phylo(0, 6)
#' dataset <- MatrixToPhyDat(`rownames<-`(rbind(1, 1, 2, 2, 3, 3), TipLabels(6)))
#' FitchInfo2(tree, dataset)
#' @export
FitchInfo2 <- function(tree, dataset) {
  # Check inputs
  if (is.null(dataset)) {
    warning("Cannot calculate concordance without `dataset`.")
    return(NULL)
  }
  if (is.null(tree)) {
    warning("Cannot calculate concordance without `dataset`.")
    return(NULL)
  }
  
  keep <- MatchStrings(TipLabels(tree), names(dataset), warning)
  if (length(keep) == 0) {
    return(NULL)
  }
  dataset <- dataset[keep]
  
  at <- attributes(dataset)
  cont <- at[["contrast"]]
  if ("-" %in% colnames(cont)) {
    cont[cont[, "-"] > 0, ] <- 1
  }
  ambiguous <- rowSums(cont) != 1
  
  mat <- matrix(as.integer(unlist(dataset)), length(dataset), byrow = TRUE)
  mat[mat %in% which(ambiguous)] <- NA_integer_
  maxToken <- max(mat, na.rm = TRUE)
  tokens <- as.character(seq_len(maxToken))
  mat <- apply(mat, 2, function (x) {
    uniques <- tabulate(x, maxToken) == 1
    x[x %in% tokens[uniques]] <- NA_integer_
    tokens <- sort(unique(x[!is.na(x)]))
    newLabel <- tabulate(tokens)
    newLabel[newLabel > 0] <- seq_len(sum(newLabel > 0)) - 1
    newLabel[x]
  }) # retains phyDat character compression
  
  # Work down the tree counting, for each node,
  # for each reconstruction (i.e. 0, 1, 02, ...)
  # (i) the number of leaves above it in each state; (ii) the number of steps
  # encountered thus far.  Then the root will enumerate all combinations
  # and number of steps.
  # This has something in common with the recursive formula of
  # Maddison & Slater 1991

  charH <- apply(mat, 2, function(char) {
    # Prune tree to fit, and process
    tr <- KeepTip(tree, !is.na(char))
    char <- char[!is.na(char)]
    
    edge <- tr[["edge"]]
    nodes <- unique(edge[, 1]) # in preorder
    desc <- lapply(seq_len(max(nodes)), function (node) edge[edge[, 1] == node, 2])
    parent <- lapply(seq_len(max(nodes)), function (node) edge[edge[, 2] == node, 1])
    parent[lengths(parent) == 0] <- which(lengths(parent) == 0)
    parent <- unlist(parent)
    nTip <- NTip(tr)
    nEdge <- nTip + nTip - 2
    nVert <- length(nodes) + nTip
    rootNode <- nTip + 1
    
    # Calculate entropy of character
    stateFreq <- tabulate(char + 1)
    stateP <- stateFreq / nTip
    rawInfo <- Entropy(stateP) * nTip
    nState <- length(stateFreq)
    if (nState < 2) {
      return(c("h" = 1, "hObs" = 0, "hMax" = 0))
    }
    dp <- `mode<-`(.DownpassOutcome(nState), "integer") # for tapply
    dpStep <- attr(dp, "step")
    dpUp <- dpFlat <- dp
    dpUp[!dpStep] <- 0
    dpFlat[dpStep] <- 0
    
    dpState <- array(0, dim = c( # 3D array with dimensions:
      2 ^ nState - 1,            # 1. Binary encoding of possible tokens
      nVert,                     # 2. Vertices of tree
      nTip - 1                   # 3. Number of steps observed
      ))
    dpState[2 ^ (seq_len(nState) - 1), 1:nTip, 1] <- stateP
    levels <- as.character(seq_len(2 ^ nState - 1))
    obsSteps <- 0
    maxSteps <- rep(0, nTip)
    
    solution <- 2 ^ char
    
  
    for (node in rev(nodes)) { # Postorder traversal
      maxSteps[node] <- sum(1, maxSteps[desc[[node]]])
      for (priorSteps in ((maxSteps[node] - 1):0)) {
        for (i in 0:priorSteps) {
          j <- priorSteps - i
          childPI <- dpState[, desc[[node]][[1]], i + 1]
          childPJ <- dpState[, desc[[node]][[2]], j + 1]
          nodeP <- outer(childPI, childPJ)
          noNewStep <- tapply(outer(childPI, childPJ), dpFlat, sum)[levels]
          dpState[, node, priorSteps + 1] <- dpState[, node, priorSteps + 1] + noNewStep
          if (priorSteps < (nTip - 2)) {
            newStep <- tapply(outer(childPI, childPJ), dpUp, sum)[levels]
            newStep[is.na(newStep)] <- 0
            dpState[, node, priorSteps + 1 + 1] <- dpState[, node, priorSteps + 1 + 1] + newStep
          }
        }
      }
      
      childSol <- rbind(solution[desc[[node]]])
      obsSteps <- obsSteps + dpStep[childSol]
      solution[node] <- dp[childSol]
    }
    stopifnot(isTRUE(all.equal(sum(dpState[, rootNode, ]), 1)))
    # + 1 because dpState[, , 1] corresponds to zero steps
    finalP <- colSums(dpState[, rootNode, ])
    c("h" = log2(sum(finalP[1:(obsSteps + 1)])) / log2(sum(finalP[1:nState])),
      "hObs" = -log2(sum(finalP[1:(obsSteps + 1)])),
      "hMax" = -log2(sum(finalP[1:nState])))
  })[, attr(dataset, "index"), drop = FALSE]
  
  totalHMax <- sum(charH["hMax", ])
  # Return:
  structure(
    if(totalHMax == 0) 1 else
      sum(charH["h", , drop = FALSE] * charH["hMax", , drop = FALSE] / totalHMax),
    byChar = charH
    )
}

.DownpassOutcome <- function(nStates) {
  # Assume we have 3 states.  Number them from zero.
  # Encode each state as 2^(0, 1, 2)
  # Each node must exhibit a state: cannot be all FALSE
  stateSpace <- as.raw(seq_len(2 ^ nStates - 1))
  intersect <- outer(stateSpace, stateSpace, `&`)
  union <- outer(stateSpace, stateSpace, `|`)
  ret <- intersect
  newStep <- intersect == as.raw(0)
  ret[newStep] <- union[newStep]
  structure(ret, step = newStep)
}


### Pseudocode from Fitch (1971)
# 
# if intersect(ancestor, node) == ancestor: # I: node contains all states in ancestor
#   node <- ancestor # II: Eliminate states not in ancestor
# else: # III
#   if descendants have no nodes in common:
#   node <- union(node, ancestor) # IV
# else:
#   node <- union(node, intersect(ancestor, union(descendants)) # V
.UppassOutcome <- function(nStates) {
  spaceSize <- 2 ^ nStates - 1
  stateSpace <- as.raw(seq_len(spaceSize))
  node1 <- .DownpassOutcome(nStates)
  node <- vapply(stateSpace, function(x) node1, node1)
  ancestor <- vapply(stateSpace, matrix, node1, spaceSize, spaceSize)
  descUnion <- vapply(stateSpace, function(x) outer(stateSpace, stateSpace, `|`),
                      node1)
  descIntersect <- vapply(stateSpace,
                          function(x) outer(stateSpace, stateSpace, `&`), node1)
  descCommon <- descIntersect != as.raw(0)
  
  # Unfamiliar logic here: we start with the 'else', then work our way through
  # the 'if' cases to overwrite any 'elses' that shouldn't have been triggered.
  
  ret <- node | (ancestor & descUnion)
  
  caseIII <- !descCommon
  ret[caseIII] <- (node | ancestor)[caseIII]
  
  caseI <- (ancestor & node) == ancestor
  ret[caseI] <- ancestor[caseI]
  ret
}
