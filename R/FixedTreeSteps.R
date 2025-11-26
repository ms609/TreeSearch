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
            
            addedCount <- countL * countR
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
