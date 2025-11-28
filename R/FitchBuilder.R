#' Fitch skeleton builder
#' 
#' @param steps Integer; Number of steps to build to
#' @param tokenCount Integer vector: Number of leaves displaying each state.
#'  
#' `FitchBuilder()` is a first step to computing how many phylogenetic trees, with
#' leaves labelled unambiguously by character states, exhibit _steps_ steps
#' under Fitch parsimony.
#' 
#' The approach we take inverts the standard Fitch parsimony procedure.
#' 
#' We first generate 'skeleton' trees with `FitchBuilder()`.  These trees
#' comprise only the minimum number of leaves necessary to produce _steps_
#' steps.  
#' Each leaf in the returned skeletons corresponds to one of the Fitch regions
#' that must be occupied by at least one observed leaf with a given label.
#' Once we have a complete set of such leaves, we work out how many
#' ways there are to add the remaining leaves, **without introducing any more
#' steps**.
#' 
#' To generate the complete list of skeleton trees, we simply work our way
#' up from a specified root.  For convenience, as the position of the root
#' is immaterial to the Fitch parsimony score, we will root our tree next to
#' the first clade containing solely token 1. This will allow us to compute
#' the number of unrooted trees with the desired score.
#' 
#' Rooting in this position, there will be no steps incurred on the left subtree.
#' We will therefore add our first step on the right subtree.
#' The root of the right subtree must thus be one
#' of the states that will trigger a step when paired with token `1`: perhaps
#' `2` or `[2, 4]`.
#' 
#' Now let's consider that node. We need to consider each unique pair of states
#' that could yield the observed node label; and we can observe whether it
#' introduces an additional step.  For each permutation we can move up another
#' node.  We continue until we've used up all our steps, or we've run out of
#' tokens.
#' 
#' We will end up generating some skeletons that are impossible, because the
#' recursive algorithm doesn't know what the right branch is doing whilst it's
#' traversing the left branch in preorder.  These will get filtered out when it
#' comes to calculating how many ways each skeleton can be populated with the
#' observed leaf labels.
#' 
#' @examples
#' # Number of trees with 2 steps for character 0011122
#' FitchBuilder(2, c(2, 3, 2)) |> FitchPopulateCount()
#' 
#' 
#' @export
FitchBuilder <- function(steps, tokenCount) {
  
  # Pay attention now, this is subtle.
  # We don't really care about the position of the root: (a, (b, (c, d))) is 
  # not distinct from (b, (a, (c, d))) or ((a, b), (c, d)).
  # What we'll do therefore is arbitrarily define the position of the root as
  # alongside our first clade of token 1s.
  
  nTaxa <- sum(tokenCount)
  nLevels <- length(tokenCount)
  nStates <- 2 ^ nLevels - 1
  states <- integer(nStates)
  states[2 ^ (seq_len(nLevels) - 1)] <- tokenCount
  dp <- `mode<-`(.DownpassOutcome(nLevels), "integer")
  dpStep <- attr(dp, "step")
  diag(dp) <- NA_integer_
  dp[lower.tri(dp, diag = TRUE)] <- NA_integer_
  
  .BR <- function(x) .BinaryRepresentation(x, digits = nLevels)
  
  .Recurse <- function(rootState, tokensAvailable, stepsAvailable) {
    
    if (stepsAvailable < 1) {
      # As we have no more steps to expend, all our descendant leaves must be
      # identical.  As no leaves are ambiguous, we must be in a region with a
      # single state.
      stopifnot(sum(.BR(rootState)) == 1)
      return(rootState)
    }
    
    # Because our leaves are unambiguous, a compound state requires at least
    # one of each of its leaves 'spare' to deploy later. For example,
    # [2, 3] can only be present at this node if there is at least
    # one 2 and at least one 3 among its children
    if (any(tokensAvailable[.BR(rootState)] < 1)) return(rootState)
    
    ret <- apply(which(dp == rootState, arr.ind = TRUE), 1, function(childStates) {
      lState <- childStates[[1]]
      rState <- childStates[[2]]
      
      if (lState == 0 && stepsAvailable == 0) {
        # We've hit the end of the line
        return(rootState)
      }
      
      addsStep <- dpStep[rbind(childStates)]
      lTokens <- .BR(lState)
      rTokens <- .BR(rState)
      
      tokensNeeded <- lTokens + rTokens
      if (any(tokensNeeded > tokensAvailable)) return(c('steps0' = NA_integer_))
      
      lStepsNeeded <- sum(lTokens) - 1L
      rStepsNeeded <- sum(rTokens) - 1L
      if (addsStep) {
        stepsAvailable <- stepsAvailable - 1L
      }
      if (lStepsNeeded + rStepsNeeded > stepsAvailable) return(c('steps0' = NA_integer_))
      
      result <- lapply(lStepsNeeded:(stepsAvailable - rStepsNeeded), function(lSteps) {
        rSteps <- stepsAvailable - lSteps
        list("L" = .Recurse(lState, tokensAvailable - rTokens, lSteps),
             "R" = .Recurse(rState, tokensAvailable - lTokens, rSteps))
      })
      `names<-`(result, paste0("steps", seq_along(result)))
    })
    `names<-`(ret, paste0("rt", seq_along(ret))) |>
      unlist(recursive = FALSE)
  }
  
  # We root next to, not within, the 1-clade.
  skeleta <- lapply(which(dpStep[1, ]), function(rState) {
    list("L" = 1,
         "R" = .Recurse(
           rState,
           c(tokensAvailable[1] - 1, tokensAvailable[-1]),
           steps - 1)
    )
  }) |>
    .ExpandChoices()
  
  message(paste(names(skeleta), sapply(skeleta, .AsNewick), collapse = "\n"))
  vapply(skeleta, function(x) tabulate(unlist(x, use.names = FALSE), nLevels),
    integer(nLevels)
  )
}

BuildCountuh <- function(...) {
  apply(FitchBuilder(...), 2, AssignLeavesToRegions, tokenCount) |> sum()
}

#' Count the number of labelled rooted binary trees with leaves in specified regions
#'
#' `AssignLeavesToRegions()` computes the total number of **labelled rooted binary trees**
#' in which a specified number of leaves are assigned to distinct regions (edges) of a
#' starting tree. Each region receives at least one leaf, and leaves are considered
#' distinct (labelled).
#'
#' For a single region with `e` chosen edges and `n` labelled leaves, the number of
#' trees is computed as:
#'
#' \deqn{
#' \text{LeafAssignment}(e, n) = \sum_{k_1 + \dots + k_e = n, k_i \ge 1}
#' \binom{n}{k_1, \dots, k_e} \prod_{i=1}^e NRooted(k_i),
#' }
#'
#' where
#'
#' \eqn{NRooted(k) = (2k - 3)!!} 
#' is the number of **labelled rooted binary trees** with `k` leaves, and
#' \eqn{\binom{n}{k_1, \dots, k_e}} is the multinomial coefficient giving the number
#' of ways to assign `n` distinct leaves to the `e` edges according to the
#' composition `(k_1, ..., k_e)`.
#'
#' The total number of trees across multiple regions is obtained by taking the product
#' over each region:
#'
#' \deqn{
#' \text{TotalTrees} = \prod_{i} \text{LeafAssignment}(regions[[i]], leaves[[i]])
#' }
#'
#' @param regions A vector of integers, where `regions[[i]]` is the number of chosen edges
#'   in region `i` that must receive at least one leaf.
#' @param leaves A vector of integers, where `leaves[[i]]` is the total number of labelled
#'   leaves to be added to region `i`.
#'
#' @return `AssignLeavesToRegions()` returns a numeric scalar giving the total number
#'   of labelled rooted binary trees consistent with the allocation of leaves to regions.
#'
#' @examples
#' # One region with 1 chosen edge, 3 leaves
#' .AssignLeavesToRegions(list(1), list(3))
#'
#' # Two regions: 1 edge with 2 leaves, 2 edges with 4 leaves
#' .AssignLeavesToRegions(list(1, 2), list(2, 4))
#'
#' @seealso [NRooted()] for the double factorial formula counting labelled rooted trees.
#' @export
AssignLeavesToRegions <- function(regions, leaves) {
  vapply(seq_along(regions), function(i) {
    .LeafAssignment(regions[[i]], leaves[[i]])
    }, double(1)) |> prod()
}
  
.LeafAssignment <- function(e, n) {
  if (n < e) return(0L)   # impossible to assign ≥1 leaf to each edge
  comps <- .CompositionsPos(n, e)
  total <- 0
  for (kvec in comps) {
    # Multinomial coefficient: ways to assign labels to the chosen edges
    multinom <- factorial(n) / prod(factorial(kvec))
    # Product of NRooted for each edge
    total <- total + multinom * prod(vapply(kvec, NRooted, double(1)))
  }
  total
}

# Compositions of n into e positive integers
.CompositionsPos <- function(n, e) {
  if (e == 1) return(list(c(n)))
  res <- list()
  k <- 1L
  for (first in 1:(n - e + 1)) {
    tails <- .CompositionsPos(n - first, e - 1)
    for (t in tails) {
      res[[k]] <- c(first, t)
      k <- k + 1L
    }
  }
  res
}


.ExpandChoices <- function(x) {
  
  is_leaf <- function(x) !is.list(x)
  
  is_binary_subtree <- function(x) {
    is.list(x) &&
      length(x) == 2 &&
      identical(sort(names(x)), c("L", "R"))
  }
  
  is_choice_set <- function(x) {
    is.list(x) &&
      !is_binary_subtree(x)
  }
  
  # 1. Leaf
  if (is_leaf(x)) return(list(x))
  
  # 2. Binary subtree: recurse on L and R
  if (is_binary_subtree(x)) {
    left_exp  <- .ExpandChoices(x$L)
    right_exp <- .ExpandChoices(x$R)
    
    out <- vector("list", length(left_exp) * length(right_exp))
    k <- 1L
    for (l in left_exp) for (r in right_exp) {
      out[[k]] <- list(L = l, R = r)
      k <- k + 1
    }
    return(out)
  }
  
  # 3. Choice set
  if (is_choice_set(x)) {
    # Union of all possibilities
    possibilities <- unlist(lapply(x, .ExpandChoices), recursive = FALSE)
    return(possibilities[!vapply(possibilities, function(x) any(is.na(x)), logical(1))])
  }
  
  stop("Unreachable: unexpected structure")
}


.AsNewick <- function(x) {
  if (!is.list(x)) {
    # atomic: just return as character
    return(as.character(x))
  }
    # recursive case
  children <- vapply(x, .AsNewick, FUN.VALUE = character(1))
  brackets <- if (any(c("L", "R") %in% names(children))) c("(", ")") else c("[", "]")
  paste0(brackets[[1]], paste(children, collapse = ", "), brackets[[2]])
}

.BinaryRepresentation <- function(x, digits = 32) {
  intToBits(x)[1:digits] != as.raw(0)
}

# TODO Delete: # For reference, here's the definition of:
# .DownpassOutcome <- function(nStates) {
#   # Assume we have 3 states.  Number them from zero.
#   # Encode each state as 2^(0, 1, 2)
#   # Each node must exhibit a state: cannot be all FALSE
#   stateSpace <- as.raw(seq_len(2 ^ nStates - 1))
#   intersect <- outer(stateSpace, stateSpace, `&`)
#   union <- outer(stateSpace, stateSpace, `|`)
#   ret <- intersect
#   newStep <- intersect == as.raw(0)
#   ret[newStep] <- union[newStep]
#   structure(ret, step = newStep)
# }

