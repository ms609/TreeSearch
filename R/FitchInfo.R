#' Information required to encode a character using a Fitch-like algorithm
#' 
#' @examples
#' dataset <- inapplicable.phyData[["Vinther2008"]]
#' tree <- TreeTools::NJTree(dataset)
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
    cont <- cont[, colnames(cont) != "-"]
  }
  states <- vapply(dataset, function (taxon) {
    cont[taxon, ] == 1
  }, matrix(logical(0), attr(dataset, "nr"), dim(cont)[[2]]))
  minimal <- apply(states, 1, function(state) {
    essential <- rowSums(state[, colSums(state) == 1]) > 0
    while(any(colSums(state[essential, , drop = FALSE]) == 0)) {
      stop("Not yet implemented")
    }
    state[essential, , drop = FALSE]
  }, simplify = FALSE)
  charStates <- vapply(minimal, nrow, 1)
  
  edge <- Postorder(tree)[["edge"]]
  nodes <- unique(edge[, 1]) # in postorder
  desc <- lapply(seq_len(max(nodes)), function (node) edge[edge[, 1] == node, 2])
  parent <- lapply(seq_len(max(nodes)), function (node) edge[edge[, 2] == node, 1])
  parent[lengths(parent) == 0] <- which(lengths(parent) == 0)
  parent <- unlist(parent)
  counts <- CharacterLength(tree, dataset)
  nTip <- NTip(tree)
  nEdge <- nTip + nTip - 2
  ic <- double(length(charStates))
  
  for (nStates in sort(unique(nStates[nStates > 1]))) {
    nStateChars <- which(charStates == nStates)
    i <- seq_len(nStates - 1)
    overl2 <- 1 / log(2)
    # TODO vectorize
    for (idx in nStateChars) {
      
      char <- minimal[[idx]]
      stateCounts <- rowSums(char)
      
      dp <- cbind(char, matrix(NA, nrow = nStates, ncol = length(nodes)))
      colnames(dp)[colnames(dp) == ""] <- which(colnames(dp) == "")
      cost <- double(dim(char)[[2]])
      for (node in nodes) {
        kids <- desc[[node]]
        common <- apply(dp[, kids], 1, all)
        if (any(common)) {
          dp[, node] <- common
        } else {
          dp[, node] <- apply(dp[, desc[[node]]], 1, any)
        }
      }
      
      up <- dp
      root <- nodes[[length(nodes)]]
      # TODO test assertion that this will lead to the lowest overall entropy
      up[, root] <- 1:nStates == which.max(stateCounts * up[, root])
      info <- -log2(stateCounts[up[, root]] / sum(stateCounts))
      
      for (node in c(rev(nodes), seq_len(dim(char)[[2]]))) {
        prnt <- parent[[node]]
        if (identical(up[, prnt], up[, node])) {
          next
        }
        parentState <- up[, prnt]
        up[, node] <- if (any(apply(up[, desc[[node]]], 1, all))) {
          apply(up[, c(desc[[node]], node), drop = FALSE], 1, all)
        } else {
          apply(up[, c(prnt, node)], 1, all)
        }
        if (identical(up[, prnt], up[, node])) {
          next
        }
        message(node, " from ", which(parentState), " to ", which(up[, node]))
        message(stateCounts[up[, node]], " / ", sum(stateCounts[!parentState]))
        info <- info - log2(stateCounts[up[, node]] / sum(stateCounts[!parentState]))
      }
      
      up
      info
      
      TreeDist::Ntropy(tabulate(apply(char, 2, which.max))) * dim(char)[[2]]
      
      transitions <- counts[[idx]]
      iTransition <- lchoose(nEdge, transitions) * overl2
      
      ic[[idx]] <- iTransition + info
    }
  }
  
}
