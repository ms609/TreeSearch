#' Plot the distribution of a character on a tree
#' 
#' Reconstructs the distribution of a character on a tree topology using the
#' modified Fitch algorithm presented in 
#' \insertCite{Brazeau2019;textual}{TreeSearch}.
#' 
#' @param tree A bifurcating tree of class `phylo`, or a list or `multiPhylo`
#' object containing such trees.
#' @inheritParams MaximizeParsimony
#' @param char Index of character to plot.
#' @param updateTips Logical; if `FALSE`, tips will be labelled with their
#' original state in `dataset`.
#' @param plot Logical specifying whether to plot the output.
#' @param tokenCol Palette specifying colours to associate with each token in
#' turn, in the sequence listed in `attr(dataset, "levels")`.
#' @param ambigCol,ambigLty,inappCol,inappLty,plainLty Colours and line types
#' to apply to ambiguous, inapplicable and applicable tokens.  See the `lty` 
#' [graphical parameter] for details of line styles.  Overrides `tokenCol`.
#' @param tipOffset Numeric: how much to offset tips from their labels.
#' @param unitEdge Logical: Should all edges be plotted with a unit length?
#' @param Display Function that takes argument `tree` and returns a tree
#' of class `phylo`, formatted as it will be plotted.
#' @param \dots Further arguments to pass to `plot.phylo()`.
#' 
#' @return `PlotCharacter()` invisibly returns a matrix in which each row
#' corresponds to a numbered tip or node of `tree`, and each column corresponds
#' to a token; the tokens that might parsimoniously be present at each point
#' on a tree are denoted with `TRUE`.
#' If multiple trees are supplied, the strict consensus of all trees and
#' reconstructions will be returned; i.e. if a node is reconstructed as $0$
#' in one tree, and $2$ in another, it will be labelled $(02)$.
#' 
#' @references 
#' \insertAllCited{}
#' @examples
#' # Set up plotting area
#' oPar <- par(mar = rep(0, 4))
#'
#' tree <- ape::read.tree(text = 
#'   "((((((a, b), c), d), e), f), (g, (h, (i, (j, (k, l))))));")
#' ## A character with inapplicable data
#' dataset <- TreeTools::StringToPhyDat("23--1??--032", tips = tree)
#' plotted <- PlotCharacter(tree, dataset)
#' plotted
#' 
#' # Character from a real dataset 
#' data("Lobo", package = "TreeTools")
#' dataset <- Lobo.phy
#' tree <- TreeTools::NJTree(dataset)
#' PlotCharacter(tree, dataset, 14)
#' par(oPar)
#' @template MRS
#' @importFrom ape plot.phylo nodelabels 
#' @importFrom graphics par
#' @importFrom TreeTools PostorderOrder
#' @export
PlotCharacter <- function(tree, dataset, char = 1L,
                          updateTips = FALSE,
                          plot = TRUE,
                          
                          tokenCol = NULL,
                          ambigCol = "grey",
                          inappCol = "lightgrey",
                          
                          ambigLty = "dotted",
                          inappLty = "dashed",
                          plainLty = par("lty"),
                          
                          tipOffset = 1,
                          unitEdge = FALSE,
                          Display = function(tree) tree,
                          ...
) {
  UseMethod("PlotCharacter")
}

#' @rdname PlotCharacter
#' @export
PlotCharacter.phylo <- function(tree, dataset, char = 1L,
                                updateTips = FALSE,
                                plot = TRUE,
                                
                                tokenCol = NULL,
                                ambigCol = "grey",
                                inappCol = "lightgrey",
                                
                                ambigLty = "dotted",
                                inappLty = "dashed",
                                plainLty = par("lty"),
                                
                                tipOffset = 1,
                                unitEdge = FALSE,
                                Display = function(tree) tree,
                                ...
) {
  
  # Reconcile labels
  datasetTaxa <- names(dataset)
  tree <- Display(tree)
  treeTaxa <- tree[["tip.label"]]
  if(!all(treeTaxa %fin% datasetTaxa)) {
    stop("Taxa in tree missing from dataset:\n  ",
         paste0(setdiff(treeTaxa, datasetTaxa), collapse = ", "))
  }
  dataset <- dataset[treeTaxa]
  
  # Read tree
  postorder <- PostorderOrder(tree)
  edgeLength <- tree[["edge.length"]][postorder]
  if (!is.null(edgeLength) && length(unique(edgeLength)) == 1) {
    tree[["edge.length"]] <- edgeLength
  }
  nNode <- tree[["Nnode"]]
  nTip <- NTip(tree)
  if (nNode != nTip - 1) {
    stop("`tree` must be bifurcating. Try TreeTools::MakeTreeBinary(tree).")
  }
  edge <- tree[["edge"]][postorder, ]
  parent <- edge[, 1]
  child <- edge[, 2]
  left <- integer(nNode + nTip)
  right <- left
  parentOf <- integer(nNode + nTip)
  for (e in seq_len(dim(edge)[1])) {
    pa <- parent[[e]]
    ch <- child[[e]]
    parentOf[[ch]] <- pa
    if (right[[pa]]) {
      left[[pa]] <- ch
    } else {
      right[[pa]] <- ch
    }
  }
  preOrderNodes <- unique(rev(parent)) # Root guaranteed first
  postOrderNodes <- rev(preOrderNodes)
  rootNode <- preOrderNodes[[1]]
  parentOf[[rootNode]] <- rootNode
  tips <- seq_len(nTip)
  
  # Read states
  if (!inherits(dataset, "phyDat")) {
    dataset <- MatrixToPhyDat(dataset)
  }
  character <- dataset[, char]
  contrast <- attr(character, "contrast") == 1
  levels <- colnames(contrast)
  inputState <- contrast[as.integer(character), , drop = FALSE]
  state <- rbind(inputState, matrix(NA, nNode, dim(contrast)[2]))
  
  if (is.na(match("-", levels))) {
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
      aState <- state[parentOf[[n]], ]
      lState <- state[left[[n]], ]
      rState <- state[right[[n]], ]
      inherited <- nState & aState
      if (all(inherited == aState)) {
        state[n, ] <- inherited
      } else if (any(lState & rState)) {
        state[n, ] <- nState | (aState & (lState | rState))
      } else {
        state[n, ] <- aState | nState
      }
    }
    
    for (n in tips) {
      nState <- state[n, ]
      aState <- state[parentOf[[n]], ]
      common <- aState & nState
      if (any(common)) {
        state[n, ] <- common
      }
    }
    
  } else {
    # Inapplicable Fitch, Brazeau, Guillerme & Smith 2019
    inappLevel <- levels == "-"
    appLevels <- !inappLevel
    
    # First downpass
    for (n in postOrderNodes) {
      lState <- state[left[[n]], ]
      rState <- state[right[[n]], ]
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
      # message ("DP1: Set node ", n, " to: ", paste0(levels[state[n, ]], collapse = ""))
    }
    
    # First uppass
    for (n in preOrderNodes) {
      nState <- state[n, ]
      aState <- if (n == rootNode && !all(state[n, ] == inappLevel)) {
        state[n, ] & appLevels
      } else {
        state[parentOf[[n]], ]
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
      # message ("UP1: Set node ", n, " to: ", paste0(levels[state[n, ]], collapse = ""))
    }
    for (n in tips) {
      nState <- state[n, ]
      aState <- state[parentOf[[n]], ]
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
      # message ("UP1: Set tip  ", n, " to: ", paste0(levels[state[n, ]], collapse = ""))
    }
    
    # Second downpass
    for (n in postOrderNodes) {
      nState <- state[n, ]
      lState <- state[left[[n]], ]
      rState <- state[right[[n]], ]
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
      # message ("DP2: Set node ", n, " to: ", paste0(levels[state[n, ]], collapse = ""))
    }
    
    # Second uppass
    for (n in preOrderNodes) {
      nState <- state[n, ]
      aState <- state[parentOf[n], ]
      lState <- state[left[[n]], ]
      rState <- state[right[[n]], ]
      # 1. If the node has any applicable token 
      if (any(nState[appLevels])) {
        # 2. If the node’s ancestor has any applicable token
        if (any(aState[appLevels])) {
          #2A [ADDED IN ERRATUM?]
          common <- aState & nState
          if (any(common) && all(common == aState)) {
            state[n, ] <- aState
          } else 
          # 3. If the node’s state is NOT the same as its ancestor’s
          # if (any(nState != aState))
            {
            # 4. If there is any token in common between the node’s descendants
            common <- lState & rState
            if (any(common)) {
              # 5. Add to the current node’s state any token in common between
              #  its ancestor and *either of* its descendants
              state[n, ] <- nState | (aState & (lState | rState))
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
      # message ("UP2: Set node ", n, " to: ", paste0(levels[state[n, ]], collapse = ""))
    }
    
    for (n in tips) {
      nState <- state[n, ]
      aState <- state[parentOf[[n]], ]
      common <- aState & nState
      if (any(common)) {
        state[n, ] <- common
        # message ("UP2: Set tip ", n, " to: ", paste0(levels[state[n, ]], collapse = ""))
      }
    }
  }
  
  if (!updateTips) {
    state[seq_len(nTip), ] <- inputState
  }
  
  hasToken <- if (length(setdiff(colnames(state), "-")) > 1L) {
    as.logical(rowSums(!state[, colnames(state) != "-", drop = FALSE]))
  } else {
    !logical(nrow(state))
  }
  anywhere <- as.logical(colSums(state[hasToken, , drop = FALSE]))
  slimState <- state[, anywhere, drop = FALSE]
  
  if (isTRUE(plot)) {
    .PlotCharacter(tree, nTip, state, levels, tokenCol, ambigCol, inappCol,
                   ambigLty, inappLty, plainLty, tipOffset, unitEdge, ...)
  }
  
  # Return:
  invisible(slimState)
}

.PlotCharacter <- function(tree, nTip, state, tokens,
                           tokenCol, ambigCol, inappCol,
                           ambigLty, inappLty, plainLty,
                           tipOffset, unitEdge, ...) {
  tokens <- colnames(state)
  
  hasToken <- if (length(setdiff(colnames(state), "-")) > 1L) {
    as.logical(rowSums(!state[, colnames(state) != "-", drop = FALSE]))
  } else {
    !logical(nrow(state))
  }
  anywhere <- as.logical(colSums(state[hasToken, , drop = FALSE]))
  slimState <- state[, anywhere, drop = FALSE]
  
  if (is.null(tokenCol)) {
    tokenCol <- tokens
    tokenCol[tokens != "-"] <- c("#00bfc6",
                                 "#ffd46f",
                                 "#ffbcc5",
                                 "#c8a500",
                                 "#ffcaf5",
                                 "#d5fb8d",
                                 "#e082b4",
                                 "#25ffd3",
                                 "#a6aaff",
                                 "#e6f3cc",
                                 "#67c4ff",
                                 "#9ba75c",
                                 "#60b17f")[seq_along(setdiff(tokens, "-"))]
    tokenCol[tokens == "-"] <- inappCol
  }
  nodeStyle <- apply(state, 1, function (tkn) {
    if (length(tkn) == 0) {
      c(col = ambigCol, lty = ambigLty)
    } else if (sum(tkn) > 1L) {
      c(col = ambigCol, lty = ambigLty)
    } else {
      c(col = tokenCol[tkn],
        lty = ifelse(tokens[tkn] == "-", inappLty, plainLty))
    }
  })
  if (unitEdge) {
    tree[["edge.length"]] <- rep_len(1, dim(tree[["edge"]])[1])
  }
  plot.phylo(tree,
             node.color = nodeStyle["col", , drop = FALSE],
             node.lty = nodeStyle["lty", , drop = FALSE],
             label.offset = tipOffset,
             ...)
  
  .NodeText <- function (n) {
    if (length(n) == 0 || (
      sum(n) > 1L && all(n[anywhere & names(n) != "-"]))) {
      "?"
    } else {
      paste0(tokens[n], collapse = "")
    }
  }
  nodelabels(apply(state, 1, .NodeText),
             seq_len(nTip + tree[["Nnode"]]),
             bg = nodeStyle["col", , drop = FALSE])
}

#' @rdname PlotCharacter
#' @importFrom TreeTools as.Splits Consensus DescendantTips TipLabels
#' @export
PlotCharacter.multiPhylo <- function(tree, dataset, char = 1L,
                                     updateTips = FALSE,
                                     plot = TRUE,
                                     
                                     tokenCol = NULL,
                                     ambigCol = "grey",
                                     inappCol = "lightgrey",
                                     
                                     ambigLty = "dotted",
                                     inappLty = "dashed",
                                     plainLty = par("lty"),
                                     
                                     tipOffset = 1,
                                     unitEdge = FALSE,
                                     Display = function(tree) tree,
                                     ...) {
  
  if (length(tree) == 1) {
    return(PlotCharacter(tree[[1]], dataset, char, updateTips, plot,
                         tokenCol, ambigCol, inappCol,
                         ambigLty, inappLty, plainLty,
                         tipOffset, unitEdge, Display, ...))
  }
  
  tipLabels <- unique(lapply(lapply(tree, TipLabels), sort))
  if (length(tipLabels) != 1) {
    stop("All trees must have the same tip labels")
  }
  tipLabels <- tipLabels[[1]]
  nTip <- length(tipLabels)
  tokens <- attr(dataset, "levels")
  reconstructions <- lapply(tree, PlotCharacter,
                            dataset = dataset, char = char,
                            updateTips = updateTips, plot = FALSE,
                            Display = function(tree) tree, ...)
  # Check labels: definitely identical, possibly in different sequence
  consTree <- Display(Consensus(tree, p = 1, check.labels = TRUE))
  .TreeClades <- function(tr) {
    ed <- tr[["edge"]]
    lab <- TipLabels(tr)
    apply(DescendantTips(ed[, 1], ed[, 2],
                         node = seq_len(nTip + tr[["Nnode"]])),
          1, function (tips) {
      paste0(sort(lab[tips]), collapse = " @||@ ")
    })
  }
  consClades <- .TreeClades(consTree)
  .Recon <- function(i) {
    reconstructions[[i]][
      match(consClades, .TreeClades(tree[[i]])), , drop = FALSE]
  }
  recon <- matrix(FALSE, nrow = length(consClades), ncol = length(tokens),
                  dimnames = list(NULL, tokens))
  for (i in seq_along(tree)) {
    ri <- .Recon(i)
    recon[, colnames(ri)] <- recon[, colnames(ri)] | ri
  }
  
  if (isTRUE(plot)) {
    .PlotCharacter(consTree, nTip, recon, tokens, tokenCol, ambigCol, inappCol,
                   ambigLty, inappLty, plainLty, tipOffset, unitEdge, ...)
  }
  
  invisible(recon)
}

#' @rdname PlotCharacter
#' @export
PlotCharacter.list <- function(tree, dataset, char = 1L,
                               updateTips = FALSE,
                               plot = TRUE,
                               
                               tokenCol = NULL,
                               ambigCol = "grey",
                               inappCol = "lightgrey",
                               
                               ambigLty = "dotted",
                               inappLty = "dashed",
                               plainLty = par("lty"),
                               
                               tipOffset = 1,
                               unitEdge = FALSE,
                               Display = function(tree) tree,
                               ...
) {
  if (all(vapply(tree, inherits, logical(1), "phylo"))) {
    PlotCharacter.multiPhylo(tree, dataset, char, updateTips, plot,
                             tokenCol, ambigCol, inappCol,
                             ambigLty, inappLty, plainLty,
                             tipOffset, unitEdge, Display, ...)
  } else {
    stop("Elements of `tree` must be of class `phylo`")
  }
}
