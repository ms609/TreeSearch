#' Addition tree
#' 
#' Generates a starting tree by adding each taxon in turn to the most
#' parsimonious location.
#' 
#' @inheritParams MaximizeParsimony
#' @param sequence Character or numeric vector listing sequence in which to add
#' taxa. Randomized if not provided.
#' @param concavity Determines the degree to which extra steps beyond the first
#' are penalized.  Specify a numeric value to use implied weighting
#' \insertCite{Goloboff1993}{TreeSearch}; `concavity` specifies _k_ in
#'  _k_ / _e_ + _k_. A value of 10 is recommended;
#' TNT sets a default of 3, but this is too low in some circumstances
#' \insertCite{Goloboff2018,Smith2019}{TreeSearch}.
#' Better still explore the sensitivity of results under a range of
#' concavity values, e.g. `k = 2 ^ (1:7)`.
#' Specify `Inf` to weight each additional step equally,
#' (which underperforms step weighting approaches
#' \insertCite{Goloboff2008,Goloboff2018,Goloboff2019,Smith2019}{TreeSearch}).
#' Specify `"profile"` to employ an approximation of profile parsimony
#' \insertCite{Faith2001}{TreeSearch}.
#' Note that tips are always placed using an equal-weights proxy, so a
#' numeric `concavity` value has no effect on the tree topology returned
#' by `AdditionTree()`: the topology has been observed to be identical
#' whatever numeric value of `concavity` is specified (including `Inf`,
#' i.e. equal weights). `AdditionTree()` does not return a score, so this
#' has no user-visible effect at all.
#' Specifying `concavity = "profile"` _does_ affect the returned topology,
#' because the underlying character data are recoded before tree
#' construction begins.
#' @examples
#' data("inapplicable.phyData", package = "TreeSearch")
#' # concavity = 10 has (empirically) no effect on the tree topology
#' # returned: placement always uses an equal-weights proxy.
#' AdditionTree(inapplicable.phyData[["Longrich2010"]], concavity = 10)
#' @template MRS
#' @return `AdditionTree()` returns a tree of class `phylo`. The tree carries a
#' degree-two root, but its root position is an arbitrary by-product of the
#' order in which taxa were added and is not necessarily `sequence[1]`;
#' parsimony scores are unaffected by rooting, so root the result yourself with
#' [`TreeTools::RootTree()`](https://ms609.github.io/TreeTools/reference/RootTree)
#' if the position of the root matters to you. With fewer than four
#' taxa there is nothing to optimise, and a pectinate tree of the dataset's taxa
#' is returned without consulting `sequence` or `constraint`.
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
    # Reject non-positive, fractional, out-of-range or duplicated indices before
    # subsetting: R's `taxa[i]` would otherwise silently drop (`i <= 0`),
    # truncate (fractional) or recycle, yielding a tree that ignores the
    # requested order rather than throwing an error.
    if (anyNA(sequence) || any(sequence != round(sequence)) ||
        any(sequence < 1L) || any(sequence > nTaxa) ||
        anyDuplicated(sequence)) {
      stop("numeric `sequence` must be distinct whole-number indices ",
           "between 1 and ", nTaxa, " (the number of taxa in `dataset`)")
    }
    sequence <- taxa[sequence]
  }
  if (anyNA(sequence) || !all(sequence %in% taxa)) {
    stop("`sequence` must list only taxa present in `dataset` ",
         "(by name, or by valid index)")
  }
  # A duplicated taxon poisons the C++ kernel's addition order: the repeated
  # tip is inserted twice and a different tip is never added, so the returned
  # tree silently contains one taxon twice and drops another (the numeric path
  # already rejects duplicates; mirror that here for character `sequence`).
  if (anyDuplicated(sequence)) {
    stop("`sequence` must not list any taxon more than once")
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

  # IW: minimum step counts per character, needed so `result$score` (an IW
  # score when `concavity` is finite) isn't computed against min_steps = 0.
  # Placement itself ignores this: see @param concavity above.
  minSteps <- if (is.finite(concavity)) {
    as.integer(MinimumLength(dataset, compress = TRUE))
  } else {
    integer(0)
  }

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
    min_steps = minSteps,
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
