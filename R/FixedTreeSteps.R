.LogCumSumExp <- function (x) { 
  n <- length(x)
  Lk <- c(x[1], double(n - 1L))
  for (k in 1L + seq_len(n - 1L)) {
    Lk[k] <- Lk[k - 1]
    Lk[k] <- max(x[k], Lk[k]) + log1p(exp(-abs(x[k] - Lk[k])))
  }
  
  # Return:
  Lk
}

# R is the probability that in a randomly selected tree of _n_ taxa, the
# smaller of the two basal subclades will have _m_ taxa.
# 
# The formula in Maddison & Slatkin 1991 overcounts by a factor of 2; rather
# than doubling the count where n != 2m (why?), we should halve the count
# where n == 2m: because choosing A B : C D gives the same trees as C D : A B,
# choose(2m, m) counts each tree twice.
# 
# Without this correction, we can see that R can exceed 1,
# e.g. R(1 | 3) = 2 * 4 * 1 * 3 / 15, using T = NRooted.
# (Note the misplaced brackets in Maddison & Slatkin 1991's definition of T(n))
# 
# TODO simplify by cancelling the choose(n,m) term with the one in .LogD
.LogR <- function(m, n) {
  log(if (n == m + m) 1 / 2 else 1) + 
    lchoose(n, m) +
    LnRooted(m) + LnRooted(n - m) - LnRooted(n)
}

# D is the probability that, in a randomly selected tree on `leaves`, the
# smaller subclade of taxa will receive taxa with labels `drawn`.
# If both subclades are the same size, the 'smaller' subclade is drawn at
# random
.LogD <- function(drawn, leaves) {
  sum(lchoose(leaves, drawn)) - lchoose(sum(leaves), sum(drawn))
}

.LogRD <- function(drawn, leaves) {
  n <- sum(leaves)
  m <- sum(drawn)
  
  log(if (n == m + m) 1 / 2 else 1) + 
    LnRooted(m) + LnRooted(n - m) - LnRooted(n) +
    sum(lchoose(leaves, drawn))
}

.LogB_cache <- new.env(parent = emptyenv())

#' @importFrom stringi stri_paste
.HashLeaves <- function(v) {
  stri_paste(v, collapse = ",")
}

# B(b | tokens) is the probability that `token` is reconstructed at the base
# of a clade with leaves labelled `leaves`
.LogB <- function(token, leaves, dp) {
  if (token < 1 || token > length(leaves)) {
    stop("`token` must be 1 .. length(leaves)")
  }
  
  leafHash <- .HashLeaves(leaves)
  sub <- .LogB_cache[[leafHash]]
  if (is.null(sub)) {
    sub <- new.env(parent = emptyenv())
    .LogB_cache[[leafHash]] <- sub
  }
  val <- sub[[as.character(token)]]
  if (!is.null(val)) {
    return(val)
  }
  
  if (missing(dp)) {
    nLevels <- floor(log2(length(leaves))) + 1
    nStates <- 2 ^ nLevels - 1
    dp <- `mode<-`(.DownpassOutcome(nLevels), "integer")
  }
  n <- sum(leaves)
  result <- if (n == 1) {
    # If the one leaf we're looking at bears the token, p = 1; else p = 0
    log(leaves[[token]] == 1)
  } else if (n == 2) {
    result <- if (any(leaves == 2)) {
      which.max(leaves)
    } else {
      dp[rbind(which(leaves > 0))]
    }
    log(if (token == result) 1 else 0)
  } else {
    apply(.ValidDraws(leaves), 1, function(drawn) {
      m <- sum(drawn)
      undrawn <- leaves - drawn
      # if both subtrees are equal size, D will give half the probability we want.
      balancedCorrection <- if (m + m == n) log(2) else log(1)
      if (all(drawn == undrawn)) {
        # Double counting in symmetrical case
        balancedCorrection <- balancedCorrection - log(2)
      }
      
      # Return:
      balancedCorrection + 
        .LogRD(drawn, leaves) +
        LogSumExp(apply(
          which(dp == token, arr.ind = TRUE), 1, function(pair) {
            LogProdExp(list(.LogB(pair[[1]], drawn, dp),
                            .LogB(pair[[2]], undrawn, dp)))
          }))
    }) |>
      LogSumExp()
  }
  sub[[as.character(token)]] <- result
  result
}


.ValidDraws <- function(leaves) {
  grid <- do.call(expand.grid, lapply(leaves, function(mc) seq.int(0, mc))) |>
    as.matrix()
  nDrawn <- rowSums(grid)
  n <- sum(leaves)
  evenSplits <- which(nDrawn == n / 2)
  duplicated <- evenSplits[
    apply(grid[evenSplits, , drop = FALSE], 1,
          function(drawn) {
            undrawn <- leaves - drawn
            cmp <- NA_integer_
            for (i in seq_along(drawn)) {
              if (drawn[[i]] < undrawn[[i]]) {cmp <- -1L; break }
              if (drawn[[i]] > undrawn[[i]]) {cmp <- +1L; break }
            }
            # Return:
            !is.na(cmp) && cmp == +1L
          })]
  validDraws <- setdiff(which(nDrawn > 0 & nDrawn <= floor(n / 2)), duplicated)
  grid[validDraws, , drop = FALSE]
}

.LogP_cache <- new.env(parent = emptyenv())


#' @rdname Carter1
#' @examples
#' # Number of trees with 2 steps for character 0011122
#' MaddisonSlatkin(2, c("0" = 2, "1" = 3, "01" = 0, "2" = 2)) * NUnrooted(7)
#' 
#' @importFrom TreeTools LnRooted
#' @export
MaddisonSlatkin <- function(steps, states) {
  nTaxa <- sum(states)
  nLevels <- floor(log2(length(states))) + 1
  nStates <- 2 ^ nLevels - 1
  length(states) <- nStates
  states[is.na(states)] <- 0
  dp <- `mode<-`(.DownpassOutcome(nLevels), "integer")
  
  # Probability of seeing s steps given tokens at leaves, and `token` at root
  .LogP <- function(s, leaves, token) {
    n <- sum(leaves)
    if (n == 1) {
      if(leaves[[token]] != 1) {
        # Impossible to have `token` at base given `leaves`
        return(log(0))
      }
      return(log(s == 0))
    }
    if (n == 2) {
      return(log(if (any(leaves == 2)) {
        s == 0
      } else {
        s == attr(dp, "step")[rbind(which(leaves > 0))]
      }))
    }
    
    
    leafHash <- .HashLeaves(leaves)
    sub <- .LogP_cache[[leafHash]]
    if (is.null(sub)) {
      sub <- new.env(parent = emptyenv())
      .LogP_cache[[leafHash]] <- sub
    }
    outHash <- paste(s, token, collapse = ",")
    val <- sub[[outHash]]
    if (!is.null(val)) {
      return(val)
    }
    
    denominator <- .LogB(token, leaves, dp)
    result <- if (is.finite(denominator)) {
      
      LogSumExp(
        apply(.ValidDraws(leaves), 1, function(drawn) {
          m <- sum(drawn)
          undrawn <- leaves - drawn
          
          # Return:
          .LogRD(drawn, leaves) +
            # If the two subtrees are the same size, we don't care whether the
            # "smaller" or "larger" corresponds to our draw, as we have filtered
            # the symmetrical duplicate out of .ValidDraws()
            # But if the draw is its own mirror, we don't want to count twice!
            log(if (sum(drawn) == sum(undrawn) && !all(drawn == undrawn)) 2 else 1) +
            LogSumExp(
              
              
              # 0:s where the conjunction doesn't add a step
              vapply(0:s, function(r) {
                apply(which(dp == token & !attr(dp, "step"), arr.ind = TRUE), 1, function(pair)
                  LogProdExp(list(
                    .LogP(r, drawn, pair[[1]]),
                    .LogB(pair[[1]], drawn, dp),
                    .LogP(s - r, undrawn, pair[[2]]),
                    .LogB(pair[[2]], undrawn, dp)
                  ))
                ) |> LogSumExp()
              }, double(1)),
              
              # 0:(s - 1) where the conjunction adds a step
              vapply(seq_len(s) - 1, function(r) {
                which(dp == token & attr(dp, "step"), arr.ind = TRUE) |>
                  # Empty call to which will deliver pair = c(0, 0)
                  apply(1, function(pair) if (pair[[1]] == 0) -Inf else {
                    LogProdExp(list(
                      .LogP(r, drawn, pair[[1]]),
                      .LogB(pair[[1]], drawn, dp),
                      # s - r - 1 because we're adding a step here
                      .LogP(s - r - 1, undrawn, pair[[2]]),
                      .LogB(pair[[2]], undrawn, dp)
                    ))
                  }) |> LogSumExp()
              }, double(1))
            )
        })
      ) - denominator
    } else {
      denominator
    }
    
    sub[[outHash]] <- result
    result
  }
  
  
  p <- log(0)
  for (state in seq_len(nStates)) {
    p <- LogSumExp(p, LogProdExp(list(
      .LogB(state, states, dp),
      .LogP(steps, states, state))
    ))
  }
  p
}

.DownpassOutcome <- function(nStates) {
  # Assume we have 3 states.  Number them from zero.
  # Encode each state as 2^(0, 1, 2)
  # Each node must exhibit a state: cannot be all FALSE
  stateSpace <- as.raw(seq_len(2 ^ nStates - 1))
  intersect <- outer(stateSpace, stateSpace, `&`)
  union <- outer(stateSpace, stateSpace, `|`)
  ret <- intersect
  newStep <- intersect == as.raw(0)
  ret[newStep] <- union[newStep]
  structure(ret, step = newStep)
}

#' @param tree A binary tree of class phylo
#' @param steps Integer vector: maximum number of steps to compute
#' @param tokens Integer vector: Occurrences of each token
#' @examples
#' tree <- TreeTools::BalancedTree(7)
#' tokens <- c(2, 3, 2) # e.g. 0 0 1 1 1 2 2
#' FixedTreeCount(tree, 2:4, tokens)
#' 
#' Note: setting `Inf` for steps will give all possible outcomes.
#' Setting a lower value will allow some recursions to terminate early,
#' potentially improving runtime - but probably not by much.
#' 
#' @importFrom TreeTools CladeSizes
#' @export
FixedTreeCount <- function(tree, tokens, steps = Inf) {
  
  # 1. Setup and Pre-processing
  nTip <- length(tree$tip.label)
  if (nTip != sum(tokens)) {
    stop(nTip, " leaves but ", sum(tokens), " tokens.")
  }
  maxSteps <- min(steps, nTip - length(tokens) + 1)
  
  nodeSizes <- TreeTools::CladeSizes(tree, internal = FALSE)
  edge <- tree$edge
  
  tokens <- tokens[tokens > 0]
  nStates <- length(tokens)
  nRows <- 2 ^ nStates - 1 # e.g., 7 rows for 3 states (Masks 1 through 7)
  
  # 2. Memoization Cache
  cache <- new.env(hash = TRUE, parent = emptyenv())
  
  # 3. Helper Functions 
  
  .ValidSplits <- function(availableTokens, leftSize) {
    # (Same as before)
    ranges <- lapply(availableTokens, function(x) 0:x)
    grid <- do.call(expand.grid, ranges)
    valid <- grid[rowSums(grid) == leftSize, , drop = FALSE]
    splitList <- split(as.matrix(valid), row(valid))
    names(splitList) <- NULL
    splitList
  }
  
  .FitchOp <- function(mask1, mask2) {
    # (Same as before, masks are 0-based integers)
    intersection <- bitwAnd(mask1, mask2)
    if (intersection > 0) {
      list(state = intersection, cost = 0)
    } else {
      list(state = bitwOr(mask1, mask2), cost = 1)
    }
  }
  
  # 4. The Recursive Core
  .Recurse <- function(node, currentTokens) {
    
    tokenKey <- paste(currentTokens, collapse = ",")
    cacheKey <- paste(node, tokenKey, sep = "|")
    
    if (exists(cacheKey, envir = cache)) {
      return(cache[[cacheKey]])
    }
    
    resMat <- matrix(0, nrow = nRows, ncol = maxSteps + 1,
                     dimnames = list(NULL, as.character(0:maxSteps)))
    
    # 🌟 BASE CASE: Leaf Node
    if (node <= nTip) {
      # The token corresponding to the leaf's label must be 1.
      # Example: If currentTokens is c(1, 0, 0), tokenIndex is 1.
      tokenIndex <- which(currentTokens == 1)
      
      # State mask is 2^(index - 1)
      stateMask <- 2^(tokenIndex - 1)
      
      # Matrix Row Index = State Mask + 1 (1-based index)
      # Column Index = 1 (for 0 steps)
      resMat[stateMask, 1] <- 1 
      
      assign(cacheKey, resMat, envir = cache)
      return(resMat)
    }
    
    # RECURSIVE STEP: Internal Node (Simplified indexing for clarity)
    children <- edge[edge[, 1] == node, 2]
    leftNode <- children[1]
    rightNode <- children[2]
    leftSize <- nodeSizes[[leftNode]]
    
    splits <- .ValidSplits(currentTokens, leftSize)
    
    for (tokensL in splits) {
      tokensR <- currentTokens - tokensL
      multiplier <- prod(choose(currentTokens, tokensL))
      if (multiplier == 0) next
      
      matL <- .Recurse(leftNode, tokensL)
      matR <- .Recurse(rightNode, tokensR)
      
      idxL <- which(matL > 0, arr.ind = TRUE)
      idxR <- which(matR > 0, arr.ind = TRUE)
      
      if (length(idxL) == 0 || length(idxR) == 0) next
      
      # Convolve Left and Right outcomes
      for (i in 1:nrow(idxL)) {
        sL  <- idxL[i, 1]       # R row index (1-based)
        stL <- idxL[i, 2] - 1      # Steps Left (0-based)
        countL <- matL[sL, stL + 1]
        
        for (j in 1:nrow(idxR)) {
          sR  <- idxR[j, 1]     # R row index (1-based)
          stR <- idxR[j, 2] - 1    # Steps Right (0-based)
          countR <- matR[sR, stR + 1]
          
          fitch <- .FitchOp(sL, sR)
          totalSteps <- stL + stR + fitch$cost
          
          if (totalSteps <= maxSteps) {
            rowIdx <- fitch$state
            # Map to 1-based row index
            colIdx <- totalSteps + 1 
            
            addedCount <- countL * countR * multiplier
            resMat[rowIdx, colIdx] <- resMat[rowIdx, colIdx] + addedCount
          }
        }
      }
    }
    
    assign(cacheKey, resMat, envir = cache)
    return(resMat)
  }
  
  # 5. Execute and Format Output
  rootNode <- nTip + 1
  #debug(.Recurse)
  finalMat <- .Recurse(rootNode, tokens)
  
  totalPerStep <- colSums(finalMat)
  step_labels <- 0:maxSteps
  names(totalPerStep) <- step_labels
  
  if (maxSteps < length(totalPerStep) - 1) {
    totalPerStep[seq_len(maxSteps + 1)]
  } else {
    totalPerStep
  }
}

#' @param tree A binary tree of class phylo
#' @param steps Integer vector: maximum number of steps to compute
#' @param tokens Integer vector: Occurrences of each token
#' @examples
#' tree <- TreeTools::BalancedTree(7)
#' tokens <- c(2, 3, 2) # e.g. 0 0 1 1 1 2 2
#' FixedTreeCount(tree, 2:4, tokens)
#' 
#' Note: setting `Inf` for steps will give all possible outcomes.
#' Setting a lower value will allow some recursions to terminate early,
#' potentially improving runtime - but probably not by much.
#' 
#' @importFrom TreeTools CladeSizes
#' @export
FixedTreeCountPrompt <- function(tree, steps = Inf, tokens) {
  edge <- Preorder(tree)[["edge"]]
  sizes <- CladeSizes(tree, internal = FALSE)
  nTip <- NTip(tree)
  
  tokens <- tokens[tokens > 0]
  
  states <- double(2 ^ length(tokens) - 1)
  states[2 ^ (seq_along(tokens) - 1)] <- tokens
  
  rootNode <- edge[[1]]
  
  # Recursive function
  .Ntrees <- function(rootNode, tips) {
    # Returns a matrix where each row corresponds to a state,
    # and each column corresponds to a number of steps (col 1 = 0 steps, 2 = 1...).
    # Each entry in the matrix denotes the number of trees with STATE at the 
    # root (under the Fitch algorithm) and STEPS steps.
    # 
    
    if (rootNode < nTip) {
      # If rootNode is a leaf, then this will be state = state at leaf; steps = 0.
      #return(<...>)
    }
    
    nodeChildren <- edge[edge[, 1] == rootNode, 2]
    
    # Get list of valid draws
    # This may be inspired by the existing .ValidDraws, but will differ because
    # the size of the subtree is known
    # tipsOnRight can be inferred from leaves - sum(tipsOnLeft)
    # Need to consider symmetry - not I think a problem because
    # ((t1 = 0, t2 = 0), (t3 = 1, t4 = 1)) is distinct from 
    # ((t1 = 1, t2 = 1), (t3 = 0, t4 = 0)) – ought to check this though!
    validDraws <- .ValidDrawsFixed(tips, tipsOnLeft = sizes[[nodeChildren[[1]]]])
    # .ValidDrawsFixed should return a list where entry 1 = drawn for left leaves;
    # element 2 = draw for leaves of right subclade
    # 
    ret <- matrix(0, ...)
    for (draw in validDraws) {
      leftSteps <- .NTrees(nodeChildren[[1]])
      rightSteps <- .NTrees(nodeChildren[[2]])
      # Entries where this node introduces a step need to be promoted
      newSteps <- attr(.DownpassOutcome(...), "step")
      destinationColumn <- .NewNSteps(leftSteps + rightSteps + newSteps)
      destinationRow <- .DownpassOutcome(leftSteps, rightSteps) # or equivalent
      
      # sum leftSteps + rightSteps count entries into appropriate row and column
      newEntries <- .Combine(leftSteps, rightSteps, destinationRow, destinationCol)
      
      ret <- ret + newEntries
    }
  }
  
  # Return:
  colSums(.NTrees(rootNode))
  
}