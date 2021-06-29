#' Calculate the parsimony score of a tree given a dataset
#'
#' `TreeLength()` uses the Morphy library \insertCite{Brazeau2017}{TreeSearch}
#' to calculate a parsimony score for a tree, handling inapplicable data 
#' according to the algorithm of \insertCite{Brazeau2019;textual}{TreeSearch}.
#' Tree scoring can employ implied weights \insertCite{Goloboff1993}{TreeSearch}
#' or profile parsimony \insertCite{Faith2001}{TreeSearch}.
#'
#' @param tree A tree of class `phylo`, a list thereof (optionally of class
#' `multiPhylo`), or an integer -- in which case `tree` random trees will be 
#' uniformly sampled.
#' @template datasetParam
#' @template concavityParam
#' 
#' @return `TreeLength()` returns a numeric vector containing the score for
#' each tree.
#' 
#' @examples
#' data("inapplicable.datasets")
#' tree <- TreeTools::BalancedTree(inapplicable.phyData[[1]])
#' TreeLength(tree, inapplicable.phyData[[1]])
#' TreeLength(tree, inapplicable.phyData[[1]], concavity = 10)
#' TreeLength(tree, inapplicable.phyData[[1]], concavity = 'profile')
#' TreeLength(5, inapplicable.phyData[[1]])
#' @seealso 
#' - Conduct tree search using [`MaximizeParsimony()`] (command line), 
#' [`EasyTrees()`] (graphical user interface), or [`TreeSearch()`]
#' (custom optimality criteria).
#' 
#' - See score for each character: [`CharacterLength()`].
#' @family tree scoring 
#' 
#' @references
#' \insertAllCited{}
#' @author Martin R. Smith (using Morphy C library, by Martin Brazeau)
#' @importFrom phangorn phyDat
#' @importFrom TreeTools Renumber RenumberTips TreeIsRooted
#' @export
TreeLength <- function (tree, dataset, concavity = Inf) UseMethod('TreeLength')

#' @rdname TreeLength
#' @export
TreeLength.phylo <- function (tree, dataset, concavity = Inf) {
  if (is.finite(concavity)) {
    if (!('min.length' %in% names(attributes(dataset)))) {
      dataset <- PrepareDataIW(dataset)
    }
    at <- attributes(dataset)
    nChar  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
    weight <- at$weight
    steps <- CharacterLength(tree, dataset, compress = TRUE)
    minLength <- at$min.length
    homoplasies <- steps - minLength
    
    # This check was once triggered - possibly fixed but remains
    # under investigation...
    if (any(homoplasies < 0)) { #nocov start
      stop("Minimum steps have been miscalculated.\n", 
           "       Please report this bug at:\n", 
           "       https://github.com/ms609/TreeSearch/issues/new\n\n",
           "       See above for full tree: ", dput(tree))
    } #nocov end
    fit <- homoplasies / (homoplasies + concavity)
    # Return:
    sum(fit * weight)
    
  } else if (.UseProfile(concavity)) {
    dataset <- PrepareDataProfile(dataset)
    steps <- CharacterLength(tree, dataset, compress = TRUE)
    info <- attr(dataset, 'info.amounts')
    
    # Return:
    sum(vapply(which(steps > 0), function (i) info[steps[i], i],
               double(1)) * attr(dataset, 'weight')[steps > 0])
  } else {
    tree <- RenumberTips(Renumber(tree), names(dataset))
    if (!TreeIsRooted(tree)) stop("`tree` must be rooted; try RootTree(tree)")
    morphyObj <- PhyDat2Morphy(dataset)
    on.exit(morphyObj <- UnloadMorphy(morphyObj))
    MorphyTreeLength(tree, morphyObj)
  }
}


#' @rdname TreeLength
#' @importFrom TreeTools RandomTree
#' @export
#TODO could be cleverer still and allow TreeLength.edge
TreeLength.numeric <- function (tree, dataset, concavity = Inf) {
  TreeLength(lapply(!logical(tree), RandomTree, tips = dataset), 
             dataset = dataset, concavity = concavity)
}

#' @rdname TreeLength
#' @export
TreeLength.list <- function (tree, dataset, concavity = Inf) {
  # Define constants
  iw <- is.finite(concavity)
  profile <- .UseProfile(concavity)
  
  tree[] <- RenumberTips(tree, dataset)
  edges <- vapply(tree, `[[`, tree[[1]]$edge, 'edge')
  
  # Initialize data
  if (profile) {
    dataset <- PrepareDataProfile(dataset)
    profiles <- attr(dataset, 'info.amounts')
  }
  if (iw || profile) {
    at <- attributes(dataset)
    characters <- PhyToString(dataset, ps = '', useIndex = FALSE,
                              byTaxon = FALSE, concatenate = FALSE)
    weight <- at$weight
    informative <- at$informative
    charSeq <- seq_along(characters) - 1L
    
    # Save time by dropping uninformative characters
    if (!is.null(informative)) charSeq <- charSeq[informative]
    morphyObjects <- lapply(characters, SingleCharMorphy)
    on.exit(morphyObjects <- vapply(morphyObjects, UnloadMorphy, integer(1)),
            add = TRUE)
  } else {
    morphyObj <- PhyDat2Morphy(dataset)
    on.exit(morphyObj <- UnloadMorphy(morphyObj), add = TRUE)
    weight <- unlist(MorphyWeights(morphyObj)[1, ]) # exact == approx
  }
  
  # Return:
  if (iw) {
    minLength <- at$min.length
    if (is.null(minLength)) {
      minLength <- attr(PrepareDataIW(dataset), 'min.length')
    }
    apply(edges, 3, morphy_iw, morphyObjects, weight, minLength, charSeq,
          concavity, Inf)
  } else if (profile) {
    apply(edges, 3, morphy_profile, morphyObjects, weight, charSeq, profiles,
          Inf)
  } else {
    apply(edges, 3, preorder_morphy, morphyObj)
  }
  
}


#' @rdname TreeLength
#' @export
TreeLength.multiPhylo <- TreeLength.list

#' @rdname TreeLength
#' @export
Fitch <- function (tree, dataset) {
  .Deprecated('TreeLength')
  TreeLength(tree, dataset, Inf)
}



#' Character length
#' 
#' Homoplasy length of each character in a dataset on a specified tree.
#' 
#' @template treeParam
#' @template datasetParam
#' @template compressParam
#'
#' @return `CharacterLength()` returns a vector listing the contribution of each
#' character to tree score, according to the algorithm of
#' \insertCite{Brazeau2018;textual}{TreeTools}.
#'
#' @examples
#' data('inapplicable.datasets')
#' dataset <- inapplicable.phyData[[12]]
#' tree <- TreeTools::NJTree(dataset)
#' CharacterLength(tree, dataset)
#' CharacterLength(tree, dataset, compress = TRUE)
#' @template MRS
#' @family tree scoring
#' @references
#' \insertAllCited{}
#' @importFrom TreeTools Renumber RenumberTips
#' @export
CharacterLength <- function (tree, dataset, compress = FALSE) {
  if (!inherits(dataset, 'phyDat')) {
    stop("Dataset must be of class phyDat, not ", class(dataset), '.')
  }
  if (!inherits(tree, 'phylo')) {
    stop("Tree must be of class phylo, not ", class(tree), '.')
  }
  if (is.null(tree$tip.label))
  if (length(tree$tip.label) < length(dataset)) {
    if (all(tree$tip.label %in% names(dataset))) {
      dataset[!(names(dataset)%in% tree$tip.label)] <- NULL
    } else {
      stop ("Tree tips ", 
            paste(tree$tip.label[!(tree$tip.label %in% names(dataset))],
                  collapse = ', '), 
            " not found in dataset.")
    }
  }
  
  tree <- RenumberTips(Renumber(tree), names(dataset))  
  
  ret <- FastCharacterLength(tree, dataset)
  # Return:
  if (compress) {
    ret
  } else {
    ret[attr(dataset, 'index')]
  }
  
}

#' @rdname CharacterLength
#' @export
FitchSteps <- function (tree, dataset) {
  .Deprecated("CharacterLength")
  CharacterLength(tree, dataset, compress = TRUE)
}

#' @describeIn CharacterLength Do not perform checks.  Use with care: may cause
#' erroneous results or  software crash if variables are in the incorrect format.
FastCharacterLength <- function (tree, dataset) {
  characters <- PhyToString(dataset, ps = '', useIndex = FALSE, byTaxon = FALSE,
                            concatenate = FALSE)
  morphyObjects <- lapply(characters, SingleCharMorphy)
  on.exit(morphyObjects <- vapply(morphyObjects, UnloadMorphy, integer(1)))
  
  # Return:
  vapply(morphyObjects, MorphyTreeLength, tree = tree, integer(1))
}

#' Calculate parsimony score from Morphy object
#' 
#' This function must be passed a valid Morphy object, or R may crash.
#' For most users, the function [`TreeLength()`] will be more appropriate.
#' 
#' @template labelledTreeParam
#' @template morphyObjParam
#'
#' @return `MorphyTreeLength()` returns the length of the tree,
#' after applying weighting.
#'
#' @seealso PhyDat2Morphy
#'
#' @family tree scoring
#' @author Martin R. Smith
#' @keywords internal
#' @export
MorphyTreeLength <- function (tree, morphyObj) {
  if (!is.morphyPtr(morphyObj)) {
    stop("`morphyObj` must be a valid morphy pointer")
  }
  nTaxa <- mpl_get_numtaxa(morphyObj)
  if (nTaxa != length(tree$tip.label)) {
    stop ("Number of taxa in morphy object (", nTaxa,
          ") not equal to number of tips in tree")
  }
  treeOrder <- attr(tree, 'order')
  inPostorder <- (!is.null(treeOrder) && treeOrder == "postorder")
  treeEdge <- tree$edge

  # Return:
  MorphyLength(treeEdge[, 1], treeEdge[, 2], morphyObj, inPostorder, nTaxa)
}

#' @describeIn MorphyTreeLength Faster function that requires internal tree
#'   parameters. Node numbering must increase monotonically away from root.
#' @template treeParent
#' @template treeChild
#' @author Martin R. Smith
#' @keywords internal
#' @importFrom TreeTools Postorder Preorder
#' @export
MorphyLength <- function (parent, child, morphyObj, inPostorder = FALSE,
                          nTaxa = mpl_get_numtaxa(morphyObj)) {
  if (!inPostorder) {
    edgeList <- Postorder(Preorder(cbind(parent, child)))
    parent <- edgeList[, 1]
    child <- edgeList[, 2]
  }
  if (!inherits(morphyObj, 'morphyPtr')) {
    stop("morphyObj must be a morphy pointer. See ?LoadMorphy().")
  }
  if (nTaxa < 1L) {
    # Run this test after we're sure that morphyObj is a morphyPtr, or lazy
    # evaluation of nTaxa will cause a crash.
    stop("Error: ", mpl_translate_error(nTaxa))
  }
  
  maxNode <- nTaxa + mpl_get_num_internal_nodes(morphyObj)
  rootNode <- nTaxa + 1L
  allNodes <- rootNode:maxNode
  
  parentOf <- parent[match(seq_len(maxNode), child)]
  parentOf[rootNode] <- rootNode # Root node's parent is a dummy node
  leftChild <- child[length(parent) + 1L - match(allNodes, rev(parent))]
  rightChild <- child[match(allNodes, parent)]
  
  # Return:
  .Call('MORPHYLENGTH', as.integer(parentOf - 1L), as.integer(leftChild - 1L), 
               as.integer(rightChild - 1L), morphyObj)
}

#' @describeIn MorphyTreeLength Fastest function that requires internal tree parameters
#' @template parentOfParam
#' @template leftChildParam
#' @template rightChildParam
#' @family tree scoring
#' @author Martin R. Smith
#' @keywords internal
#' @export
GetMorphyLength <- function (parentOf, leftChild, rightChild, morphyObj) {
  # Return:
  .Call('MORPHYLENGTH', as.integer(parentOf), as.integer(leftChild), 
               as.integer(rightChild), morphyObj)
}

#' @describeIn MorphyTreeLength Direct call to C function. Use with caution.
#' @param parentOf For each node, numbered in postorder, the number of its parent node.
#' @param leftChild  For each internal node, numbered in postorder, the number of its left 
#'                   child node or tip.
#' @param rightChild For each internal node, numbered in postorder, the number of its right
#'                   child node or tip.
#' @keywords internal
#' @export
C_MorphyLength <- function (parentOf, leftChild, rightChild, morphyObj) {
  .Call('MORPHYLENGTH', as.integer(parentOf - 1L), as.integer(leftChild - 1L), 
               as.integer(rightChild - 1L), morphyObj)
}
