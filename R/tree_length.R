#' @title Calculate parsimony score with inapplicable data
#'
#' @description Uses code modified from the Morphy library to calculate a 
#' parsimony score in datasets that contain inapplicable data.
#' 
#' 
#' 
#' If `concavity` is non-infinite, using implied weights
#' (Goloboff 1997).
#'
#' @param tree A tree of class `phylo`, a list thereof (optionally of class
#' `multiPhylo`, or an integer -- in which case `tree` random trees will be 
#' uniformly sampled.
#' @template datasetParam
#' @template concavityParam
#' 
#' @return `TreeLength()` returns the elements from a list containing:
#'    \itemize{
#' \item     The total parsimony score
#' \item     The parsimony score associated with each character 
#' \item     A matrix comprising character reconstructions for each node
#'           after the final pass
#'   }
#'   
#' The elements to return are specified by the parameter `detail`.  
#' If a single element is requested (default) then just that element will be returned
#' If multiple elements are requested then these will be returned in a list.
#' 
#' 
#' @return #TODO if concavity, then The 'fit', `h / h + k`, where `h` is the amount of homoplasy ('extra steps') 
#'         and `k` is a constant (the 'concavity constant')
#' 
#' 
#' @examples
#' data("inapplicable.datasets")
#' tree <- TreeTools::BalancedTree(inapplicable.phyData[[1]])
#' TreeLength(tree, inapplicable.phyData[[1]])
#' TreeLength(tree, inapplicable.phyData[[1]], 10)
#' @seealso 
#' - [`TreeSearch()`]
#' @family tree scoring
#' 
#' 
#' @references
#'  - \insertRef{Goloboff1997}{TreeSearch}
#'  
#'  - \insertRef{Smith2019}{TreeSearch}
#'
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
    steps <- CharacterLength(tree, dataset)
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
    
  } else if (tolower(concavity) == 'profile') {
    if (!('info.amounts' %in% names(attributes(dataset)))) {
      dataset <- PrepareDataProfile(dataset)
    }
    steps <- CharacterLength(tree, dataset)
    info <- attr(dataset, 'info.amounts')
    # Return:
    sum(vapply(seq_along(steps), function (i) info[steps[i], i], double(1)) *
          attr(dataset, 'weight'))
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
  profile <- tolower(concavity) == 'profile'
  
  edges <- vapply(tree, `[[`, tree[[1]]$edge, 'edge')
  
  # Initialize data
  if (profile) {
    dataset <- PrepareDataProfile(dataset)
    originalLevels <- attr(dataset, 'levels')
    if ('-' %in% originalLevels) {
      #TODO Fixing this will require updating the counts table cleverly
      # Or we could use approximate info amounts, e.g. by treating '-' as 
      # an extra token
      message("Inapplicable tokens '-' treated as ambiguous '?' for profile parsimony")
      cont <- attr(dataset, 'contrast')
      cont[cont[, '-'] != 0, ] <- 1
      attr(dataset, 'contrast') <- cont[, colnames(cont) != '-']
      attr(dataset, 'levels') <- originalLevels[originalLevels != '-']
    }
    profiles <- attr(dataset, 'info.amounts')
  }
  if (iw || profile) {
    at <- attributes(dataset)
    characters <- PhyToString(dataset, ps = '', useIndex = FALSE,
                              byTaxon = FALSE, concatenate = FALSE)
    weight <- at$weight
    charSeq <- seq_along(characters) - 1L
    morphyObjects <- lapply(characters, SingleCharMorphy)
    on.exit(morphyObjects <- vapply(morphyObjects, UnloadMorphy, integer(1)),
            add = TRUE)
  } else {
    morphyObj <- PhyDat2Morphy(dataset)
    on.exit(morphyObj <- UnloadMorphy(morphyObj), add = TRUE)
    weight <- unlist(MorphyWeights(morphyObj)[1, ]) # exact == approx
  }
  
  
  if (iw) {
    nLevel <- length(at$level)
    nChar <- at$nr
    nTip <- length(dataset)
    cont <- at$contrast
    if (is.null(colnames(cont))) colnames(cont) <- as.character(at$levels)
    
    inappLevel <- at$levels == '-'
    
    if (any(inappLevel)) {
      # TODO this is a workaround until MinimumLength can handle {-, 1}
      cont[cont[, inappLevel] > 0, ] <- 0
      ambiguousToken <- at$allLevels == '?'
      cont[ambiguousToken, ] <- colSums(cont[!ambiguousToken, ]) > 0
    }
    
    # Perhaps replace with previous code:
    # inappLevel <- which(at$levels == "-")
    # cont[, inappLevel] <- 0
    
    powersOf2 <- 2L ^ c(0L, seq_len(nLevel - 1L))
    tmp <- as.integer(cont %*% powersOf2)
    unlisted <- unlist(dataset, use.names = FALSE)
    binaryMatrix <- matrix(tmp[unlisted], nChar, nTip, byrow = FALSE)
    minLength <- apply(binaryMatrix, 1, MinimumLength)
  }
  
  # Return:
  if (iw) {
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
#'
#' @return `CharacterLength()` returns a vector listing the contribution of each
#' character to tree score, according to the algorithm of Brazeau, Guillerme 
#' and Smith (2019).
#'
#' @examples
#' data('inapplicable.datasets')
#' dataset <- inapplicable.phyData[[12]]
#' tree <- TreeTools::NJTree(dataset)
#' CharacterLength(tree, dataset)
#'
#' @family tree scoring
#' @references
#'  \insertRef{Brazeau2018}{TreeTools}
#' @importFrom TreeTools Renumber RenumberTips
#' @export
CharacterLength <- function (tree, dataset) {
  if (!inherits(dataset, 'phyDat')) {
    stop ("Dataset must be of class phyDat, not ", class(dataset))
  }
  if (length(tree$tip.label) < length(dataset)) {
    if (all(tree$tip.label %in% names(dataset))) {
      dataset[!(names(dataset)%in% tree$tip.label)] <- NULL
    } else {
      stop ("Tree tips", 
            paste(tree$tip.label[!(tree$tip.label %in% names(dataset))], sep = ', '), 
            "not found in dataset.")
    }
  }
  
  tree <- RenumberTips(Renumber(tree), names(dataset))  
  
  # Return:
  FastCharacterLength(tree, dataset)
}

#' @rdname CharacterLength
FitchSteps <- function (tree, dataset) {
  .Deprecated(CharacterLength)
  CharacterLength(tree, dataset)
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

#' Calculate parsimony score with inapplicable data
#' 
#' @template labelledTreeParam
#' @template morphyObjParam
#'
#' @return The length of the tree (after weighting)
#'
#' @seealso PhyDat2Morphy
#'
#' @family tree scoring
#' @author Martin R. Smith
#' @keywords internal
#' @export
MorphyTreeLength <- function (tree, morphyObj) {
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
#'   parameters
#' @template treeParent
#' @template treeChild
#' @author Martin R. Smith
#' @keywords internal
#' @importFrom TreeTools Postorder
#' @export
MorphyLength <- function (parent, child, morphyObj, inPostorder = FALSE,
                          nTaxa = mpl_get_numtaxa(morphyObj)) {
  if (!inPostorder) {
    edgeList <- Postorder(cbind(parent, child))
    parent <- edgeList[, 1]
    child <- edgeList[, 2]
  }
  if (nTaxa < 1L) stop("Error: ", mpl_translate_error(nTaxa))
  if (!inherits(morphyObj, 'morphyPtr')) {
    stop("morphyObj must be a morphy pointer. See ?LoadMorphy().")
  }
  
  maxNode <- nTaxa + mpl_get_num_internal_nodes(morphyObj)
  rootNode <- nTaxa + 1L
  allNodes <- rootNode:maxNode
  
  parentOf <- parent[match(1:maxNode, child)]
  parentOf[rootNode] <- rootNode # Root node's parent is a dummy node
  leftChild <- child[length(parent) + 1L - match(allNodes, rev(parent))]
  rightChild <- child[match(allNodes, parent)]
  
  # Return:
  .Call('MORPHYLENGTH', as.integer(parentOf -1L), as.integer(leftChild -1L), 
               as.integer(rightChild -1L), morphyObj)
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
  .Call('MORPHYLENGTH', as.integer(parentOf -1L), as.integer(leftChild -1L), 
               as.integer(rightChild -1L), morphyObj)
}

#' Extract character data from dataset
#'
#' Specifies how characters are stored in a dataset object
#' 
#' @param dataset a matrix or list containing symbols associated with each tip;
#'  for example, a dataset of class \code{phyDat}
#' @param tips vector detailing the tips to be selected, whether as their
#'        names or numbers corresponding to their rows/columns
#' @return an integer vector, listing the tokens associated with each character for each tip in turn
#'         - ready to send to FITCH or equivalent C routine
#' @keywords internal
#' @export
TipsAreNames <- function(dataset, tips) as.integer(unlist(dataset[tips]))

#TODO Github issue #2
###   #' @describeIn TipsAreNames use if each row in a matrix corresponds to a tip
###   #' @keywords internal
###   #' @export
###   TipsAreRows <- function(dataset, tips) as.integer(dataset[tips, ])

#' @describeIn TipsAreNames use if each column in a matrix corresponds to a tip
#' @keywords internal
#' @export
TipsAreColumns <- function(dataset, tips) as.integer(dataset[, tips])
