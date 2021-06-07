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
        #state[n, ] <- nState | (aState & (lState | rState))
        state[n, ] <- aState & lState & rState
      } else {
        state[n, ] <- aState | nState
      }
    }
    
    for (n in tips) {
      nState <- state[n, ]
      aState <- state[parentOf[n], ]
      common <- aState & nState
      if (any(common)) {
        state[n, ] <- common
      }
    }
    
  } else {
    # Inapplicable Fitch, Brazeau, Guillerme & Smith 2019
    inappLevel <- levels == '-'
    appLevels <- !inappLevel
    
    # First downpass
    for (n in postOrderNodes) {
      lState <- state[left[n], ]
      rState <- state[right[n], ]
      common <- lState & rState
      if (any(common)) { # 2
        # If the token in common is only the inapplicable token, 
        # and both descendants have an applicable token
        if (all(common == inappLevel) && 
            any(lState[appLevels]) && 
            any(rState[appLevels])
            ) {
          # Set the node’s state to be the union of the descendants’ states
          state[n, ] <- lState | rState
        } else {
          # set the node’s state to be the token in common between both descendants
          state[n, ] <- common
        }
      } else { # 3
        # If both descendants have an applicable token
        if (any(lState[appLevels]) && any(rState[appLevels])) {
          # set the node’s state to be the union of both descendants’ states
          # without the inapplicable token
          state[n, ] <- (lState | rState) & appLevels
        } else {
          # set the node’s state to be the union of its descendants’ states
          state[n, ] <- lState | rState
        }
      }
    }
    
    # First uppass
    for (n in preOrderNodes) {
      nState <- state[n, ]
      aState <- state[parentOf[n], ]
      lState <- state[left[n], ]
      rState <- state[right[n], ]
      # 1. If the node has the inapplicable token
      if (any(nState[inappLevel])) {
        # 2. If the node also has an applicable token
        if (any(nState[appLevels])) {
          # 3. If the node’s ancestor has the inapplicable token
          if (any(aState[inappLevel])) {
            # set the node’s state to be the inapplicable token only
            state[n, ] <- inappLevel
          } else {
            # remove the inapplicable token from the current node’s state
            state[n, ] <- nState & appLevels
          }
        } else {
          # 4. If the node’s ancestor has the inapplicable token
          if (any(aState[inappLevel])) {
            # set the node’s state to be the inapplicable token only
          } else {
            # 5. If any of the descendants have an applicable token
            if (any(lState[appLevels]) || any(rState[appLevels])) {
              # set the node’s state to be the union of the applicable states
              # of its descendants
              state[n, ] <- (lState | rState) & appLevels
            } else {
              # set the node’s state to be the inapplicable token only
              state[n, ] <- inappLevel
            }
          }
        }
      }
    }
    for (n in tips) {
      nState <- state[n, ]
      aState <- state[parentOf[n], ]
      # 6. If the unvisited tip includes both inapplicable and applicable tokens
      if (any(nState[inappLevel]) && any(nState[appLevels])) {
        # 7. If the current node has only the inapplicable token
        if (all(aState == inappLevel)) {
          # set the tip’s state to the inapplicable token only
          state[n, ] <- inappLevel
        } else {
          # remove the inapplicable token from the tip’s state
          state[n, ] <- nState & appLevels
        }
      }
    }
  }
    
  plot(tree)
  
  nodelabels(apply(state, 1, function (n) {
    paste0(levels[n], collapse = '')
  }), seq_len(nTip + nNode))
  }
