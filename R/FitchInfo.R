#' Information required to encode a character using a Fitch-like algorithm
#' 
#' @examples
#' dataset <- inapplicable.phyData[["Vinther2008"]]
#' tree <- TreeTools::NJTree(dataset)
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
  }) # retains compression
  
  apply(mat, 2, function(char) {
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
    
    # Calculate entropy of character
    stateFreq <- tabulate(char + 1)
    stateP <- stateFreq / nTip
    rawInfo <- Entropy(stateP) * nTip
    nState <- length(stateFreq)
    dp <- `mode<-`(.DownpassOutcome(nState), "integer") # for tapply
    up <- `mode<-`(.UppassOutcome(nState), "integer")
    
    dpState <- matrix(0, 2 ^ nState - 1, nVert)
    dpState[2 ^ (seq_len(nState) - 1), 1:nTip] <- stateP
    solution <- 2 ^ char
    
    for (node in rev(nodes)) { # Postorder traversal
      childP <- dpState[, desc[[node]]]
      nodeP <- outer(childP[, 1], childP[, 2])
      dpState[, node] <- tapply(outer(childP[, 1], childP[, 2]), dp, sum)
      
      childSol <- solution[desc[[node]]]
      solution[node] <- dp[childSol[[1]], childSol[[2]]]
    }
    
    for (node in nodes) { # Preorder traversal
      # Populate solution with standard Fitch algorithm
      childSol <- solution[desc[[node]]]
      ancSol <- solution[[parent[[node]]]]
      solution[[node]] <- up[childSol[[1]], childSol[[2]], ancSol]
    }
    
    upState <- dpState # Can probably avoid a copy here and overwrite dpState
    info <- `length<-`(double(0), nVert)
    
    for (node in nodes) { # Preorder traversal
      # Prior probabilities at node
      descs <- desc[[node]]
      childP <- dpState[, descs]
      childPMat <- outer(childP[, 1], childP[, 2])
      ancSol <- solution[[parent[[node]]]]
      nodeSol <- solution[[node]]
      
      # Now we can provide some information
      info[[node]] <- -log2(upState[nodeSol, node])
      
      condP <- childPMat * (up[, , ancSol] == nodeSol)
      condP <- condP / sum(condP)
      
      # Update the probabilities of the left child conditional on this node's
      # state, which we know now:
      left <- descs[[1]]
      upState[, left] <- rowSums(condP)
      
      cp <- childP[, 2] * (up[solution[[left]], , ancSol] == nodeSol)
      
      # By the time we get to the right child, we will also know the left
      # child's state, so we must condition on that too:
      condP <- condP[solution[[left]], ]
      condP <- condP / sum(condP)
      stopifnot(abs(cp / sum(cp) - condP) < sqrt(.Machine$double.eps))
      
      right <- descs[[2]]
      upState[, right] <- condP
    }
    
    
    for (tip in seq_len(nTip)) {
      info[[tip]] <- -log2(upState[2 ^ char[tip], tip])
    }
    info
    
  })
  
}

.DownpassOutcome <- function(nStates) {
  # Assume we have 3 states.  Number them from zero.
  # Encode each state as 2^(0, 1, 2)
  # Each node must exhibit a state: cannot be all FALSE
  stateSpace <- as.raw(seq_len(2 ^ nStates - 1))
  intersect <- outer(stateSpace, stateSpace, `&`)
  union <- outer(stateSpace, stateSpace, `|`)
  ret <- intersect
  ret[intersect == as.raw(0)] <- union[intersect == as.raw(0)]
  ret
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
