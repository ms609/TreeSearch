#' Plot the distribution of a character on a tree
#' 
#' Reconstructs the distribution of a character on a tree topology using the
#' modified Fitch algorithm presented in Brazeau _et al._ (2019).
#' 
#' @template treeParam
#' @template datasetParam
#' @param char Index of character to plot.
#' @param \dots Further arguments to pass to `plot.phylo()`.
#' @return `PlotCharacter()` returns a matrix in which each row corresponds
#' to a numbered tip or node of `tree`, and each column corresponds to a 
#' token; the tokens that might parsimoniously be present at each point
#' on a tree are denoted with `TRUE`.
#' 
#' @references 
#' - \insertRef{Brazeau2019}{TreeSearch}
#' @examples
#' tree <- TreeTools::BalancedTree(12)
#' tree <- ape::read.tree(text = 
#'   "((((((a, b), c), d), e), f), (g, (h, (i, (j, (k, l))))));")
#' ## A character with inapplicable data
#' dataset <- StringToPhyDat("23--1??--032", tips = tree)
#' PlotCharacter(tree, dataset)
#' 
#' 
#' data("Lobo", package = "TreeTools")
#' dataset <- Lobo.phy
#' tree <- TreeTools::NJTree(dataset)
#' oPar <- par(mar = rep(0, 4))
#' PlotCharacter(tree, dataset, 1)
#' par(oPar)
#' @template MRS
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
  if (!inherits(dataset, 'phylo')) {
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
      message ("DP1: Set node ", n, " to: ", paste0(levels[state[n, ]], collapse = ''))
    }
    
    # First uppass
    for (n in preOrderNodes) {
      nState <- state[n, ]
      aState <- if (n == rootNode && !all(state[n, ] == inappLevel)) {
        state[n, ] & appLevels
      } else {
        state[parentOf[n], ]
      }
      
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
      message ("UP1: Set node ", n, " to: ", paste0(levels[state[n, ]], collapse = ''))
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
      message ("UP1: Set tip  ", n, " to: ", paste0(levels[state[n, ]], collapse = ''))
    }
    
    # Second downpass
    for (n in postOrderNodes) {
      nState <- state[n, ]
      lState <- state[left[n], ]
      rState <- state[right[n], ]
      # If the node had an applicable token in the first uppass
      if (any(nState[appLevels])) {
        # 3. If there is any token in common between both descendants
        common <- lState & rState
        if (any(common)) {
          # 4. If the tokens in common are applicable
          if (any(common[appLevels])) {
            # set the node’s state to be the tokens held in common,
            #  without the inapplicable token
            state[n, ] <- common & appLevels
          } else {
            # set the node’s state to be the inapplicable token
            state[n, ] <- inappLevel
          }
        } else {
          # 5. Set the node’s state to be the union of the states of both
          # descendants (if present) without the inapplicable token
          state[n, ] <- (lState | rState) & appLevels 
        }
      }
      message ("DP2: Set node ", n, " to: ", paste0(levels[state[n, ]], collapse = ''))
    }
    
    # Second uppass
    for (n in preOrderNodes) {
      nState <- state[n, ]
      aState <- state[parentOf[n], ]
      lState <- state[left[n], ]
      rState <- state[right[n], ]
      # 1. If the node has any applicable token 
      if (any(nState[appLevels])) {
        # 2. If the node’s ancestor has any applicable token
        if (any(aState[appLevels])) {
          # 3. If the node’s state is NOT the same as its ancestor’s
          if (any(nState != aState)) {
            # 4. If there is any token in common between the node’s descendants
            common <- lState & rState
            if (any(common)) {
              # 5. Add to the current node’s state any token in common between
              #  its ancestor and its descendants
              state[n, ] <- nState | (aState & common)
            } else {
              # 6. If the states of the node’s descendants both contain the
              #  inapplicable token
              if (any(lState[inappLevel]) && any(rState[inappLevel])) {
                # 7. If there is any token in common between either of the 
                # node’s descendants and its ancestor
                if (any(lState & aState) || any(rState & aState)) {
                  # set the node’s state to be its ancestor’s state
                  state[n, ] <- aState
                } else { 
                  # set the current node’s state to be all applicable tokens 
                  # common to both its descendants and ancestor
                  state[n, ] <- appLevels & common & aState
                }
              } else {
                # 8. Add to the node’s state the tokens of its ancestor
                state[n, ] <- nState | aState
              }
            }
          }
        }
      }
      message ("UP2: Set node ", n, " to: ", paste0(levels[state[n, ]], collapse = ''))
    }
    
    for (n in tips) {
      nState <- state[n, ]
      aState <- state[parentOf[n], ]
      common <- aState & nState
      if (any(common)) {
        state[n, ] <- common
        message ("UP2: Set tip ", n, " to: ", paste0(levels[state[n, ]], collapse = ''))
      }
    }
  }
    
  plot(tree)
  
  nodelabels(apply(state, 1, function (n) paste0(levels[n], collapse = '')),
             seq_len(nTip + nNode))
  nodelabels(adj = 3, bg = 'yellow')
  
  # Return:
  state
}
