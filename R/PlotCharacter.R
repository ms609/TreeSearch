#' Plot the distribution of a character on a tree
#' 
#' Reconstructs the distribution of a character on a tree topology using the
#' modified Fitch algorithm presented in Brazeau _et al._ (2019).
#' 
#' @template treeParam
#' @template datasetParam
#' @char Index of character to plot.
#' @param \dots Further arguments to pass to `plot.phylo()`.
#' @references 
#' - \insertRef{Brazeau2019}{TreeSearch}
#' @examples
#' tree <- TreeTools::BalancedTree(12)
#' ## A character with inapplicable data
#' dataset <- StringToPhyDat("23--1??--032", tips = tree)
#' PlotCharacter(tree, dataset)
#' 
#' 
#' 
#' data("Lobo", package = "TreeTools")
#' dataset <- Lobo.phy
#' tree <- TreeTools::NJTree(dataset)
#' oPar <- par(mar = rep(0, 4))
#' PlotCharacter(tree, dataset, 1)
#' par(oPar)
#' @export
PlotCharacter <- function (tree, dataset, char = 1L, ...) {
  
  # Read tree
  tree <- Postorder(tree)
  nNode <- tree$Nnode
  nTip <- NTip(tree)
  edge <- tree$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  left <- integer(nNode + nTip)
  right <- left
  parentOf <- integer(nNode + nTip)
  for (e in seq_len(dim(edge)[1])) {
    pa <- parent[e]
    ch <- child[e]
    parentOf[ch] <- pa
    if (right[pa]) {
      left[pa] <- ch
    } else {
      right[pa] <- ch
    }
  }
  postOrderNodes <- unique(parent)
  preOrderNodes <- rev(postOrderNodes)
  rootNode <- preOrderNodes[1]
  parentOf[rootNode] <- rootNode
  tips <- seq_len(nTip)
  
  # Read states
  if (!is.phylo(dataset)) {
    dataset <- MatrixToPhyDat(dataset)
  }
  character <- dataset[, char]
  contrast <- attr(character, 'contrast') == 1
  levels <- colnames(contrast)
  state <- rbind(contrast[as.integer(character), ],
                 matrix(NA, nNode, dim(contrast)[2]))
  
  if (is.na(match('-', levels))) {
    # Standard Fitch
    for (n in postOrderNodes) {
      lState <- state[left[n], ]
      rState <- state[right[n], ]
      common <- lState & rState
      if (any(common)) {
        state[n, ] <- common
      } else {
        state[n, ] <- lState | rState
        # Also add to score
      }
    }
    
    for (n in preOrderNodes) {
      nState <- state[n, ]
      aState <- state[parentOf[n], ]
      lState <- state[left[n], ]
      rState <- state[right[n], ]
      inherited <- nState & aState
      if (all(inherited == aState)) {
        state[n, ] <- inherited
      } else if (any(lState & rState)) {
        state[n, ] <- nState | (aState & (lState | rState))
      } else {
        state[n, ] <- aState | nState
      }
    }
    
  } else {
    
  }
    
  plot(tree)
  
  nodelabels(apply(state, 1, function (n) {
    paste0(levels[n], collapse = '')
  }), seq_len(nTip + nNode))
  }
