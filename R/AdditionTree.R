#' Addition tree
#' 
#' Generates a starting tree by adding each taxon in turn to the most
#' parsimonious location.
#' 
#' @inheritParams MaximizeParsimony
#' @param sequence Character or numeric vector listing sequence in which to add
#' taxa. Randomized if not provided.
#' @examples
#' data("inapplicable.phyData", package = "TreeSearch")
#' AdditionTree(inapplicable.phyData[["Longrich2010"]], concavity = 10)
#' @template MRS
#' @return `AdditionTree()` returns a tree of class `phylo`, rooted on
#' `sequence[1]`.
#' @importFrom TreeTools PectinateTree Renumber
#' @family tree generation functions
#' @seealso 
#' 
#' Impose a constraint: [`TreeTools::ImposeConstraint()`](
#' https://ms609.github.io/TreeTools/reference/ImposeConstraint)
#' 
#' Neighbour-joining trees: [`TreeTools::NJTree()`](
#' https://ms609.github.io/TreeTools/reference/NJTree.html);
#' [`TreeTools::ConstrainedNJ()`](
#' https://ms609.github.io/TreeTools/reference/ConstrainedNJ)
#' @export
AdditionTree <- function(dataset, concavity = Inf, constraint, sequence) {

  if (!inherits(dataset, "phyDat")) {
    stop("`dataset` must be a `phyDat` object")
  }
  taxa <- names(dataset)
  nTaxa <- length(taxa)

  if (nTaxa < 4L) {
    return(PectinateTree(taxa))
  }

  # Build addition order
  if (missing(sequence)) {
    sequence <- taxa[[1]]
  } else if (is.numeric(sequence)) {
    sequence <- taxa[sequence]
  }
  if (anyNA(sequence) || !all(sequence %in% taxa)) {
    stop("`sequence` must list only taxa present in `dataset` ",
         "(by name, or by valid index)")
  }
  unlisted <- setdiff(taxa, sequence)
  if (length(unlisted) > 0L) {
    sequence <- c(sequence, sample(unlisted))
  }
  addition_order <- match(sequence, taxa)

  # Profile parsimony: simplify data and extract info_amounts
  useProfile <- !missing(concavity) && identical(concavity, "profile")
  profileArgs <- list()
  if (useProfile) {
    dataset <- PrepareDataProfile(dataset)
    infoAmounts <- attr(dataset, "info.amounts")
    if (!is.null(infoAmounts) && length(infoAmounts) > 0L) {
      profileArgs$infoAmounts <- infoAmounts
    }
    concavity <- Inf
  }
  # NaN/NA slip past `is.finite() && <= 0` and would reach the kernel as a
  # non-finite double, silently selecting equal weights; reject them explicitly.
  if (!is.numeric(concavity) || length(concavity) != 1L || is.na(concavity)) {
    stop("`concavity` must be a single number (or Inf for equal weights, ",
         "or \"profile\" for profile parsimony).")
  }
  if (is.finite(concavity) && concavity <= 0) {
    stop("`concavity` must be positive (or Inf for equal weights, ",
         "or \"profile\" for profile parsimony).")
  }

  # Extract data matrices
  at <- attributes(dataset)
  contrast <- at$contrast
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = nTaxa, byrow = TRUE)
  weight <- .ScaleWeight(at$weight)
  levels <- at$levels

  # Constraint
  consArgs <- list()
  if (!missing(constraint)) {
    consArgs <- .PrepareConstraint(constraint, dataset)
  }

  # Call C++ Wagner tree
  searchArgs <- list(
    contrast = contrast,
    tip_data = tip_data,
    weight = weight,
    levels = levels,
    addition_order = addition_order,
    concavity = as.double(concavity)
  )
  result <- do.call(ts_wagner_tree, c(searchArgs, consArgs, profileArgs))

  # Reconstruct phylo from edge matrix
  tree <- list(
    edge = result$edge,
    tip.label = taxa,
    Nnode = nTaxa - 1L
  )
  class(tree) <- "phylo"
  Renumber(tree)
}


.ConstraintConstrains <- function(constraint) {
  if (is.null(constraint) || length(constraint) == 0L) return(FALSE)
  if (length(constraint[[1]]) < 1) {
    FALSE
  } else {
    contrast <- attr(constraint, "contrast")
    if (is.null(contrast) || dim(contrast)[[2]] < 2) {
      FALSE
    } else {
      cont <- `mode<-`(contrast, "logical")
      nLevel <- dim(contrast)[[1]]
      exclude <- vapply(seq_len(nLevel), function(i) {
        colSums(apply(cont, 1, `&`, cont[i, ])) == 0
      }, logical(nLevel))
      splits <- exclude * tabulate(unlist(constraint), nLevel)
      any(splits[lower.tri(splits)] > 1 & t(splits)[lower.tri(splits)] > 1)
    }
  }
}


.Recompress <- function(dataset) {
  TreeTools::MatrixToPhyDat(TreeTools::PhyDatToMatrix(dataset))
}
