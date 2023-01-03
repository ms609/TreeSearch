#' Reverse neighbour breeding tree
#' 
#' Generates a starting tree by pairing taxa and computing ancestors.
#' (#TODO explain properly.)
#' 
#' @inheritParams MaximizeParsimony
#' @param sequence Character or numeric vector listing sequence in which to add
#' taxa. Randomized if not provided.
#' @examples
#' data("Lobo", package = "TreeTools")
#' PedigreeTree(Lobo.phy)
#' @template MRS
#' @return `PedigreeTree()` returns a tree of class `phylo`, rooted on
#' `sequence[1]`.
#' @importFrom TreeTools AddUnconstrained Hamming RenumberEdges
#' PectinateTree
#' @importFrom cli cli_progress_bar cli_progress_update
#' @family tree generation functions
#' @seealso
#' 
#' Impose a constraint: [`TreeTools::ImposeConstraint()`](https://ms609.github.io/TreeTools/reference/ImposeConstraint)
#' 
#' Neighbour-joining trees: [`TreeTools::NJTree()`](https://ms609.github.io/TreeTools/reference/NJTree.html);
#' [`TreeTools::ConstrainedNJ()`](https://ms609.github.io/TreeTools/reference/ConstrainedNJ)
#' @export
PedigreeTree <- function (dataset, concavity = Inf, constraint, sequence) {
  # Initialize missing parameters
  taxa <- names(dataset)
  if (missing(sequence)) {
    sequence <- taxa[1]
  } else if (is.numeric(sequence)) {
    sequence <- taxa[sequence]
  }
  nTaxa <- length(taxa)
  if (length(taxa) < 4) {
    return(PectinateTree(taxa))
  }
  unlisted <- setdiff(taxa, sequence)
  if (length(unlisted) > 0) {
    sequence <- c(sequence, sample(unlisted))
  }
  if (!missing(constraint)) {
    constraint <- AddUnconstrained(constraint, taxa)
  }
  
  nEdge <- nTaxa + nTaxa - 2
  nNode <- nTaxa - 1
  nVisits <- nNode - 1
  parent <- integer(nEdge)
  child <- integer(nEdge)
  dat <- dataset
  names(dat) <- seq_along(dat)
  
  cli_progress_bar("Pedigree tree", total = nVisits)
  for (node in nTaxa + seq_len(nVisits)) {
    cli_progress_update(1)
    # TODO prohibit merges that violate constraint
    distances <- Hamming(dat, ratio = FALSE)
    minima <- which(as.matrix(distances) == min(distances), arr.ind = TRUE)
    unbreed <- minima[sample.int(nrow(minima), 1) ,]
    
    edgeI <- (node - nTaxa) * 2 - 1:0
    parent[edgeI] <- node
    child[edgeI] <- as.numeric(names(dat)[unbreed])
    
    newNode <- PhyDatToMatrix(.FitchAncestor(dat[unbreed, ]))
    rownames(newNode) <- node
    dat <- MatrixToPhyDat(rbind(
      PhyDatToMatrix(dat[-unbreed, ]),
      newNode
    ))
  }
  edgeI <- nTaxa * 2 - 3:2
  parent[edgeI] <- nTaxa + nNode
  child[edgeI] <- as.numeric(names(dat))
  
  tree <- PectinateTree(taxa)
  tree$edge <- do.call(cbind, RenumberEdges(parent, child))
  tree
}

.FitchAncestor <- function (children) {
  cont <- data.frame(t(attr(children, "contrast") > 0))
  token <- children[[1]]
  ret <- children[1, ]
  for (i in seq_len(attr(children, "nr"))) {
    pair <- children[, i]
    union <- cont[, pair[[1]]] & cont[, pair[[2]]]
    if (any(union)) {
      token[i] <- match(data.frame(union), cont)
    } else {
      intersect <- data.frame(cont[, pair[[1]]] | cont[, pair[[2]]],
                              fix.empty.names = FALSE)
      token[i] <- match(intersect, cont)
      if(is.na(token[i])) {
        cont <- cbind(cont, intersect)
        token[i] <- length(cont)
        attr(ret, "allLevels") <- .AddLevel(attr(ret, "allLevels"))
      }
    }
  }
  ret[[1]] <- token
  attr(ret, "contrast") <- t(as.matrix(cont) * 1)
  rownames(attr(ret, "contrast")) <- NULL
  ret
}

.AddLevel <- function(x) {
  c(x, setdiff(c(letters, LETTERS, 0:9), x)[1])
}