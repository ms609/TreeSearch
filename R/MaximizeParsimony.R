#' Find most parsimonious trees
#' 
#' Search for most parsimonious trees using the parsimony ratchet and 
#' \acronym{TBR} rearrangements, treating inapplicable data as such using the
#' algorithm of \insertCite{Brazeau2019;textual}{TreeSearch}.
#'  
#' Tree search will be conducted from a specified or automatically-generated
#' starting tree in order to find a tree with an optimal parsimony score,
#' under implied or equal weights, treating inapplicable characters as such
#' in order to avoid the artefacts of the standard Fitch algorithm
#' \insertCite{@see @Maddison1993; @Brazeau2019}{TreeSearch}.
#' Tree length is calculated using the MorphyLib C library
#' \insertCite{Brazeau2017}{TreeSearch}.
#' 
#' Tree search commences with `ratchIter` iterations of the parsimony ratchet
#' \insertCite{Nixon1999}{TreeSearch}, which bootstraps the input dataset 
#' in order to escape local optima.
#' A final round of tree bisection and reconnection (\acronym{TBR})
#' is conducted to broaden the sampling of trees.
#' 
#' This function can be called using the R command line / terminal, or through
#' the 'shiny' graphical user interface app (type `EasyTrees()` to launch).
#' 
#'  
#' For detailed documentation of the 'TreeSearch' package, including full 
#' instructions for loading phylogenetic data into R and initiating and 
#' configuring tree search, see the 
#' [package documentation](https://ms609.github.io/TreeSearch/).
#'  
#' 
#' 
#' @template datasetParam
#' @param tree (optional) A bifurcating tree of class \code{\link{phylo}},
#' containing only the tips listed in `dataset`, from which the search
#' should begin.
#' If unspecified, an [addition tree][AdditionTree()] will be generated from
#'  `dataset`, respecting any supplied `constraint`.
#' Edge lengths are not supported and will be deleted.
#' @param ratchIter Numeric specifying number of iterations of the 
#' parsimony ratchet \insertCite{Nixon1999}{TreeSearch} to conduct.
#' @param tbrIter Numeric specifying the maximum number of \acronym{TBR}
#' break points to evaluate before concluding each search.
#' The counter is reset to zero each time tree score improves.
#' The counter is reset to zero each time tree score improves.
#' One 'iteration' comprises breaking a single branch and evaluating all 
#' possible reconnections.
#' @param startIter Numeric: an initial round of tree search with
#' `startIter` &times; `tbrIter` \acronym{TBR} break points is conducted in
#' order to locate a local optimum before beginning ratchet searches. 
#' @param finalIter Numeric: a final round of tree search will evaluate
#' `finalIter` &times; `tbrIter` \acronym{TBR} break points, in order to
#' sample the final optimal neighbourhood more intensely.
#' @param maxHits Numeric specifying the maximum times that an optimal
#' parsimony score may be hit before concluding a ratchet iteration or final 
#' search concluded.
#' @param maxTime Numeric: after `maxTime` minutes, stop tree search at the
#' next opportunity.
#' @param quickHits Numeric: iterations on subsampled datasets
#'  will retain `quickHits` &times; `maxHits` trees with the best score.
#' @param concavity Numeric specifying concavity constant for implied step 
#' weighting.
#' The most appropriate value will depend on the dataset, but values around
#' 10--15 often perform well \insertCite{Goloboff2018,Smith2019}{TreeSearch}.
#' The character string "profile" employs an approximation of profile parsimony
#' \insertCite{Faith2001}{TreeSearch}.
#' Set as `Inf` for equal step weights, which underperforms step weighting
#' approaches
#' \insertCite{Goloboff2008,Goloboff2018,Goloboff2019,Smith2019}{TreeSearch}.
#' @param tolerance Numeric specifying degree of suboptimality to tolerate
#' before rejecting a tree.  The default, `sqrt(.Machine$double.eps)`, retains
#' trees that may be equally parsimonious but for rounding errors.  
#' Setting to larger values will include trees suboptimal by up to `tolerance`
#' in search results, which may improve the accuracy of the consensus tree
#' (at the expense of resolution) \insertCite{Smith2019}{TreeSearch}.
#' @param constraint An object of class `phyDat`; returned trees will be
#' perfectly compatible with each character in `constraint`.
#' See [`ImposeConstraint()`] and 
#' [vignette](https://ms609.github.io/TreeSearch/articles/tree-search.html)
#' for further examples.
#' @param verbosity Integer specifying level of messaging; higher values give
#' more detailed commentary on search progress. Set to `0` to run silently.
#' @param \dots Additional parameters to `MaximizeParsimony()`.
#' 
#' @return `MaximizeParsimony()` returns a list of trees with class
#' `multiPhylo`. This lists all trees found during each search step that
#' are within `tolerance` of the optimal score, listed in the sequence that
#' they were first visited, and named according to the step in which they were
#' first found; it may contain more than `maxHits` elements.
#' Note that the default search parameters may need to be increased in order for
#' these trees to be the globally optimal trees; examine the messages printed
#' during tree search to evaluate whether the optimal score has stabilized.
#' 
#' The return value has the attribute `firstHit`, a named integer vector listing
#' the number of optimal trees visited for the first time in each stage of
#' the tree search. Stages are named:
#' - `seed`: starting trees;
#' - `start`: Initial TBR search;
#' - `ratchN`: Ratchet iteration `N`;
#' - `final`: Final TBR search.
#' The first tree hit for the first time in ratchet iteration three is named
#' `ratch3_1`.
#' 
#' @examples
#' ## Only run examples in interactive R sessions
#' if (interactive()) {
#'   # launch 'shiny' point-and-click interface
#'   EasyTrees()
#'   
#'   # Here too, use the "continue search" function to ensure that tree score
#'   # has stabilized and a global optimum has been found
#' }
#' 
#' 
#' # Load data for analysis in R
#' library('TreeTools')
#' data('congreveLamsdellMatrices', package = 'TreeSearch')
#' dataset <- congreveLamsdellMatrices[[42]]
#' 
#' # A very quick run for demonstration purposes
#' trees <- MaximizeParsimony(dataset, ratchIter = 0, startIter = 0,
#'                            tbrIter = 1, maxHits = 4, maxTime = 1/100,
#'                            concavity = 10, verbosity = 4)
#' names(trees)
#'
#' # In actual use, be sure to check that the score has converged on a global
#' # optimum, conducting additional iterations and runs as necessary.
#'  
#' if (interactive()) {
#' # Jackknife resampling
#' nReplicates <- 10
#' jackTrees <- replicate(nReplicates,
#'   #c() ensures that each replicate returns a list of trees
#'   c(Resample(dataset, trees, ratchIter = 0, tbrIter = 2, startIter = 1,
#'              maxHits = 5, maxTime = 1 / 10,
#'              concavity = 10, verbosity = 0))
#'  )
#' 
#' # In a serious analysis, more replicates would be conducted, and each
#' # search would undergo more iterations.
#' 
#' # Now we must decide what to do with the multiple optimal trees from
#' # each replicate.
#' 
#' # Treat each tree equally
#' JackLabels(ape::consensus(trees), unlist(jackTrees, recursive = FALSE))
#' 
#' # Take the strict consensus of all trees for each replicate
#' JackLabels(ape::consensus(trees), lapply(jackTrees, ape::consensus))
#' 
#' # Take a single tree from each replicate (the first; order's irrelevant)
#' JackLabels(ape::consensus(trees), lapply(jackTrees, `[[`, 1))
#' }
#' 
#' # Tree search with a constraint
#' constraint <- MatrixToPhyDat(c(a = 1, b = 1, c = 0, d = 0, e = 0, f = 0))
#' characters <- MatrixToPhyDat(matrix(
#'   c(0, 1, 1, 1, 0, 0,
#'     1, 1, 1, 0, 0, 0), ncol = 2,
#'   dimnames = list(letters[1:6], NULL)))
#' MaximizeParsimony(characters, constraint = constraint, verbosity = 0)
#' 
#' @template MRS
#' 
#' @importFrom cli cli_alert cli_alert_danger cli_alert_info cli_alert_success cli_alert_warning
#' cli_h1 
#' cli_progress_bar cli_progress_done cli_progress_update
#' @importFrom fastmatch fmatch
#' @importFrom stats runif
#' @importFrom TreeTools
#' AddUnconstrained 
#' AllDescendantEdges
#' CharacterInformation
#' ConstrainedNJ 
#' DropTip
#' ImposeConstraint
#' MakeTreeBinary
#' NTip
#' @references
#' \insertAllCited{}
#' @seealso
#' Tree search _via_ graphical user interface: [`EasyTrees()`]
#' 
#' @encoding UTF-8
#' @export
MaximizeParsimony <- function (dataset, tree,
                               ratchIter = 6L,
                               tbrIter = 2L,
                               startIter = 2L, finalIter = 1L,
                               maxHits = NTip(dataset) * 1.8,
                               maxTime = 60,
                               quickHits = 1 / 3,
                               concavity = Inf,
                               tolerance = sqrt(.Machine$double.eps),
                               constraint,
                               verbosity = 3L) {

  ### User messaging functions ###
  .Message <- function (level, ...) {
    if (level < verbosity) {
      cli_alert(paste0(...))
    }
  }
  .Heading <- function (text, ...) {
    if (0 < verbosity) {
      cli_h1(text)
      cli_alert(paste0(...))
    }
  }
  .Info <- function (level, ...) {
    if (level < verbosity) {
      cli_alert_info(paste0(...))
    }
  }
  .Success <- function (level, ...) {
    if (level < verbosity) {
      cli_alert_success(paste0(...))
    }
  }
  
  ### Tree score functions ###
  .EWScore <- function (edge, morphyObj, ...) {
    preorder_morphy(edge, morphyObj)
  }
  
  .IWScore <- function (edge, morphyObjs, weight, charSeq, concavity, 
                        minLength, target = Inf) {
    morphy_iw(edge, morphyObjs, weight, minLength, charSeq,
              concavity, target + epsilon)
  } 
  
  # Must have same order of parameters as .IWScore, even though minLength unused
  .ProfileScore <- function (edge, morphyObjs, weight, charSeq, profiles, 
                             minLength, target = Inf) {
    morphy_profile(edge, morphyObjs, weight, charSeq, profiles,
                   target + epsilon)
  }
  
  .Score <- function (edge) {
    if (length(dim(edge)) == 3L) {
      edge <- edge[, , 1]
    }
    if (profile) {
      .ProfileScore(edge, morphyObjects, startWeights, charSeq, profiles)
    } else if (iw) {
      .IWScore(edge, morphyObjects, startWeights, charSeq, concavity, minLength)
    } else {
      preorder_morphy(edge, morphyObj)
    }
  }
  
  ### Tree search functions ###
  .TBRSearch <- function (Score, name,
                          edge, morphyObjs, weight,
                          tbrIter, maxHits,
                          minLength = NULL, charSeq = NULL, concavity = NULL) {
  
    iter <- 0L
    nHits <- 1L
    hold <- array(NA, dim = c(dim(edge), max(maxHits * 1.1, maxHits + 10L)))
    maxHits <- ceiling(maxHits)
    hold[, , 1] <- edge
    bestScore <- Score(edge, morphyObjs, weight, charSeq, concavity, minLength)
    bestPlusEps <- bestScore + epsilon
    cli_progress_bar(name, total = maxHits, 
                     auto_terminate = FALSE,
                     clear = verbosity < 3L,
                     format_done = paste0("  - TBR rearrangement at depth {iter}",
                                          " found score {signif(bestScore)}",
                                          " {nHits} time{?s}."))
    
    while (iter < tbrIter) {
      iter <- iter + 1L
      optTbr <- sample(3:(nTip * 2 - 2))
      .Message(4L, "New TBR iteration (depth ", iter, 
               ", score ", signif(bestScore), ")")
      cli_progress_update(set = 0, total = length(optTbr))
      
      for (brk in optTbr) {
        cli_progress_update(1, status = paste0('D', iter, ", score ",
                                               signif(bestScore), ", hit ",
                                               nHits, "."))
        .Message(7L, "Break ", brk)
        moves <- TBRMoves(edge, brk)
        improvedScore <- FALSE
        nMoves <- length(moves)
        moveList <- sample.int(nMoves)
        for (i in seq_along(moveList)) {
          move <- moves[[moveList[i]]]
          if (.Forbidden(move)) {
            .Message(10L, "Skipping prohibited topology")
            next
          }
          moveScore <- Score(move, morphyObjs, weight, charSeq, concavity, 
                             minLength, bestPlusEps)
          if (moveScore < bestPlusEps) {
            edge <- move
            if (moveScore < bestScore) {
              improvedScore <- TRUE
              iter <- 0L
              bestScore <- moveScore
              bestPlusEps <- bestScore + epsilon
              nHits <- 1L
              hold[, , 1] <- edge
              .Message(5L, "New best score ", signif(bestScore),
                       " at break ", fmatch(brk, optTbr), "/", length(optTbr))
              break
            } else {
              .Message(6L, "Best score ", signif(bestScore),
                       " hit again (", nHits, "/", ceiling(maxHits), ")")
              nHits <- nHits + 1L
              hold[, , nHits] <- edge
              if (nHits >= maxHits) break
            }
          }
          # If an early iteration improves the score, a later iteration will
          # probably improve it even more; we may as well keep working through
          # the list instead of calculating a new one (which takes time)
          if (improvedScore && runif(1) < (i / nMoves) ^ 2) break
        }
        if (nHits >= maxHits) break
        pNextTbr <- (fmatch(brk, optTbr) / length(optTbr)) ^ 2
        if (improvedScore && runif(1) < pNextTbr) break
      }
      if (nHits >= maxHits) break
    }
    cli_progress_done()
    
    # Return:
    unique(hold[, , seq_len(nHits), drop = FALSE], MARGIN = 3L)
  
  }

  
  .Search <- function (name = 'TBR search', .edge = edge, .hits = searchHits,
                       .weight = startWeights) {
    if (length(dim(.edge)) == 3L) {
      .edge <- .edge[, , 1]
    }
    if (profile) {
      .TBRSearch(.ProfileScore, name, edge = .edge, morphyObjects, 
                 tbrIter = searchIter, maxHits = .hits,
                 weight = .weight, minLength = minLength, charSeq = charSeq,
                 concavity = profiles)
  
    } else if (iw) {
      .TBRSearch(.IWScore, name, edge = .edge, morphyObjects, 
                 tbrIter = searchIter, maxHits = .hits,
                 weight = .weight, minLength = minLength, charSeq = charSeq,
                 concavity = concavity)
    } else {
      .TBRSearch(.EWScore, name, edge = .edge, morphyObj, 
                 tbrIter = searchIter, maxHits = .hits)
    }
  }
  
  .Timeout <- function () {
    if (Sys.time() > stopTime) {
      .Info(1L, "Stopping search at ", Sys.time(), ": ", maxTime,
            " minutes have elapsed.",
            "  Best score was ", signif(.Score(bestEdges[, , 1])), '.',
            if (maxTime == 60) "\nIncrease `maxTime` for longer runs.")
      return (TRUE)
    }
    
    FALSE
  }
  
  .ReturnValue <- function (bestEdges) {
    if (verbosity > 0L) {
      cli_alert_success(paste0(Sys.time(), 
                               ": Tree search terminated with score {.strong ",
                               "{signif(.Score(bestEdges[, , 1]))}}"))
    }
    firstHit <- attr(bestEdges, 'firstHit')
    structure(lapply(seq_len(dim(bestEdges)[3]), function (i) {
      tr <- tree
      tr$edge <- bestEdges[, , i]
      if (any(is.na(outgroup))) {
        tr
      } else {
        RootTree(tr, outgroup)
      }
    }),
    firstHit = firstHit,
    names = paste0(rep(names(firstHit), firstHit), '_', unlist(lapply(firstHit, seq_len))),
    class = 'multiPhylo')
  }
  
  
  # Define constants
  epsilon <- tolerance
  pNextTbr <- 0.33
  profile <- .UseProfile(concavity)
  iw <- is.finite(concavity)
  constrained <- !missing(constraint)
  startTime <- Sys.time()
  stopTime <- startTime + as.difftime(maxTime, units = 'mins')
  
  # Initialize tree
  startTrees <- NULL
  if (missing(tree)) {
    tree <- AdditionTree(dataset, constraint = constraint,
                         concavity = concavity)
  } else if (inherits(tree, 'multiPhylo')) {
    startTrees <- unique(tree)
    sampledTree <- sample.int(length(tree), 1)
    .Info(2L, paste0("Starting search from {.var tree[[", sampledTree, "]]}"))
    tree <- tree[[sampledTree]]
  } else if (inherits(tree, 'phylo')) {
    startTrees <- c(tree)
  }
  if (dim(tree$edge)[1] != 2 * tree$Nnode) {
    cli_alert_warning("`tree` is not bifurcating; collapsing polytomies at random")
    tree <- MakeTreeBinary(tree)
    if (dim(tree$edge)[1] != 2 * tree$Nnode) {
      cli_alert_warning("Rooting `tree` on first leaf")
      tree <- RootTree(tree, 1)
    }
    if (dim(tree$edge)[1] != 2 * tree$Nnode) {
      stop("Could not make `tree` binary.")
    }
  }
  
  # Check tree labels matches dataset
  leaves <- tree$tip.label
  taxa <- names(dataset)
  treeOnly <- setdiff(leaves, taxa) 
  datOnly <- setdiff(taxa, leaves) 
  if (length(treeOnly)) {
    cli_alert_warning(paste0("Ignoring taxa on tree missing in dataset:\n>   ",
                      paste0(treeOnly, collapse = ', ')))
    warning("Ignored taxa on tree missing in dataset:\n   ",
             paste0(treeOnly, collapse = ', '))
    tree <- DropTip(tree, treeOnly)
    startTrees <- DropTip(startTrees, treeOnly)
  }
  if (length(datOnly)) {
    cli_alert_warning(paste0("Ignoring taxa in dataset missing on tree:\n>   ",
                      paste0(datOnly, collapse = ', ')))
    warning("Ignored taxa in dataset missing on tree:\n>   ",
            paste0(datOnly, collapse = ', '))
    dataset <- dataset[-fmatch(datOnly, taxa)]
  }
  if (constrained) {
    consTaxa <- names(constraint)
    treeOnly <- setdiff(tree$tip.label, consTaxa)
    if (length(treeOnly)) {
      constraint <- AddUnconstrained(constraint, treeOnly)
    }
    consOnly <- setdiff(consTaxa, tree$tip.label)
    if (length(consOnly)) {
      cli_alert_warning(paste0("Ignoring taxa in constraint missing on tree:\n>   ", 
              paste0(consOnly, collapse = ', ')))
      warning("Ignored taxa in constraint missing on tree:\n   ", 
              paste0(consOnly, collapse = ', '))
      constraint <- constraint[-fmatch(consOnly, consTaxa)]
    }
    constraint <- constraint[names(dataset)]
  }
  
  
  tree <- Preorder(RenumberTips(tree, names(dataset)))
  nTip <- NTip(tree)
  edge <- tree$edge
  
  # Initialize constraints
  if (constrained) {
    morphyConstr <- PhyDat2Morphy(constraint)
    on.exit(morphyConstr <- UnloadMorphy(morphyConstr), add = TRUE)
    constraintWeight <- attr(constraint, 'weight')
    if (any(constraintWeight > 1)) {
      cli_alert_warning("Some constraints are exact duplicates.")
    }
    # Calculate constraint minimum score
    constraintLength <- sum(MinimumLength(constraint, compress = TRUE) *
                              constraintWeight)
    
    .Forbidden <- function (edges) {
      preorder_morphy(edges, morphyConstr) != constraintLength
    }
    
    # Check that starting tree is consistent with constraints 
    if (.Forbidden(edge)) {
      cli_alert_warning("Modifying `tree` to match `constraint`...")
      outgroup <- edge[AllDescendantEdges(edge[, 1], edge[, 2])[1, ], 2]
      outgroup <- outgroup[outgroup <= nTip]
      tree <- RootTree(ImposeConstraint(tree, constraint), outgroup)
      # RootTree leaves `tree` in preorder
      edge <- tree$edge
      if (.Forbidden(edge)) {
        stop("Could not reconcile starting tree with `constraint`. ",
             "Are all constraints compatible?")
      }
    }
    
    cli_alert_success(paste0("Initialized ", length(constraintWeight),
                             " distinct constraints."))
    
  } else {
    .Forbidden <- function (edges) FALSE
  }
  
  
  if (edge[1, 2] > nTip) {
    outgroup <- edge[AllDescendantEdges(edge[, 1], edge[, 2])[1, ], 2]
    outgroup <- outgroup[outgroup <= nTip]
    if (length(outgroup) > nTip / 2L) {
      outgroup <- seq_len(nTip)[-outgroup]
    }
    tree <- RootTree(tree, 1)
    edge <- tree$edge
  } else {
    outgroup <- NA
  }
  
  # Initialize data
  if (profile) {
    dataset <- PrepareDataProfile(dataset)
    originalLevels <- attr(dataset, 'levels')
    if ('-' %fin% originalLevels) {
      #TODO Fixing this will require updating the counts table cleverly
      # Or we could use approximate info amounts, e.g. by treating '-' as 
      # an extra token
      cli_alert_info(paste0("Inapplicable tokens '-' treated as ambiguous '?' ",
                            "for profile parsimony"))
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
    startWeights <- at$weight
    minLength <- MinimumLength(dataset, compress = TRUE)
    morphyObjects <- lapply(characters, SingleCharMorphy)
    on.exit(morphyObjects <- vapply(morphyObjects, UnloadMorphy, integer(1)),
            add = TRUE)
    
    nLevel <- length(at$level)
    nChar <- at$nr
    nTip <- length(dataset)
    cont <- at$contrast
    if (is.null(colnames(cont))) colnames(cont) <- as.character(at$levels)
    simpleCont <- ifelse(rowSums(cont) == 1,
                         apply(cont != 0, 1, function (x) colnames(cont)[x][1]),
                         '?')
  
    
    unlisted <- unlist(dataset, use.names = FALSE)
    tokenMatrix <- matrix(simpleCont[unlisted], nChar, nTip)
    charInfo <- apply(tokenMatrix, 1, CharacterInformation)
    needsInapp <- rowSums(tokenMatrix == '-') > 2
    inappSlowdown <- 3L # A guess
    # Crude estimate of score added per unit processing time
    rawPriority <- charInfo / ifelse(needsInapp, inappSlowdown, 1)
    priority <- startWeights * rawPriority
    informative <- needsInapp | charInfo > 0
    # Will work from end of sequence to start.
    charSeq <- seq_along(charInfo)[informative][order(priority[informative])] - 1L
  } else {
    morphyObj <- PhyDat2Morphy(dataset)
    on.exit(morphyObj <- UnloadMorphy(morphyObj), add = TRUE)
    startWeights <- unlist(MorphyWeights(morphyObj)[1, ]) # exact == approx
  }
  
  # Initialize variables and prepare search
  
  .Heading(paste0("BEGIN TREE SEARCH (k = ", concavity, ")"))
  nHits <- 1L
  tbrStart <- startIter > 0
  tbrEnd <- finalIter > 0
  if (is.null(startTrees)) {
    bestEdges <- edge
    dim(bestEdges) <- c(dim(bestEdges), 1)
    bestScore <- .Score(edge)
  } else {
    starters <- RenumberTips(startTrees, names(dataset))
    startEdges <- vapply(lapply(starters, Preorder), `[[`, startTrees[[1]]$edge,
                        'edge')
    startScores <- apply(startEdges, 3, .Score)
    bestScore <- min(startScores)
    bestEdges <- startEdges[, , startScores == bestScore, drop = FALSE]
  }
  nStages <- sum(tbrStart, ratchIter, tbrEnd)
  attr(bestEdges, 'firstHit') <- c('seed' = dim(bestEdges)[3],
    setNames(double(nStages),
             c(if(tbrStart) 'start',
               if(ratchIter > 0) paste0('ratch', seq_len(ratchIter)),
               if(tbrEnd) 'final')))
  
  
  
  
  
  # Find a local optimum
  
  if (tbrStart) {
    searchIter <- tbrIter * startIter
    searchHits <- maxHits
    
    .Heading("Find local optimum",
             " TBR depth ", as.integer(searchIter),
             "; keeping max ", as.integer(searchHits),
             " trees; k = ", concavity, ".")
    .Info(1L, Sys.time(), ": Score to beat: ", signif(bestScore))

    newEdges <- .Search('TBR search 1')
    
    newBestScore <- .Score(newEdges)
    scoreImproved <- newBestScore + epsilon < bestScore
    bestEdges <- if (scoreImproved) {
      .ReplaceResults(bestEdges, newEdges, 2)
    } else {
      .CombineResults(bestEdges, newEdges, 2)
    }
    if (.Timeout()) {
      return(.ReturnValue(bestEdges))                                           # nocov
    }
    edge <- bestEdges[, , 1L]
  }
  
  searchIter <- tbrIter
  searchHits <- maxHits * quickHits
  bestScore <- .Score(edge)
  bestPlusEps <- bestScore + epsilon
  
  if (ratchIter > 0L) {
    
    .Heading("Escape local optimum", "{ratchIter} ratchet iterations; ", 
             "TBR depth {ceiling(searchIter)}; ",
             "max. {ceiling(searchHits)} hits; ",
             "k = {concavity}.")
    cli_alert("{Sys.time()}: Score to beat: {.strong {signif(bestScore)}}")
    
    iter <- 0L
    while (iter < ratchIter) {
      iter <- iter + 1L
      .Message(1L, "Ratchet iteration {iter} @ ",
               "{format(Sys.time(), '%H:%M:%S')}",
               "; score to beat: {.strong {signif(bestScore)} }")
      verbosity <- verbosity - 1L
      eachChar <- seq_along(startWeights)
      deindexedChars <- rep.int(eachChar, startWeights)
      resampling <- tabulate(sample(deindexedChars, replace = TRUE),
                             length(startWeights))
      if (profile || iw) {
        priority <- resampling * rawPriority
        sampled <- informative & resampling > 0
        ratchSeq <- seq_along(charInfo)[sampled][order(priority[sampled])] - 1L
        ratchetTrees <- .Search('Bootstrapped search', .weight = resampling)
      } else {
        errors <- vapply(eachChar, function (i) 
          mpl_set_charac_weight(i, resampling[i], morphyObj), integer(1))
        if (any(errors)) {                                                      # nocov start
          stop ("Error resampling morphy object: ",
                mpl_translate_error(unique(errors[errors < 0L])))
        }
        if (mpl_apply_tipdata(morphyObj) -> error) {
          stop("Error applying tip data: ", mpl_translate_error(error))
        }                                                                       # nocov end
        
        ratchetTrees <- .Search('Bootstrapped search')
        
        errors <- vapply(eachChar, function (i) 
          mpl_set_charac_weight(i, startWeights[i], morphyObj), integer(1))
        if (any(errors)) stop ("Error resampling morphy object: ",
                               mpl_translate_error(unique(errors[errors < 0L])))
        if (mpl_apply_tipdata(morphyObj) -> error) {
          stop("Error applying tip data: ", mpl_translate_error(error))
        }
      }
      if (.Timeout()) {
        return(.ReturnValue(bestEdges))                                         # nocov
      }
      
      verbosity <- verbosity + 1L
      ratchetStart <- ratchetTrees[, , sample.int(dim(ratchetTrees)[3], 1)]
      
      ratchetImproved <- .Search('TBR search', .edge = ratchetStart,
                                 .hits = maxHits)
      ratchetScore <- .Score(ratchetImproved[, , 1])
      
      if (ratchetScore < bestPlusEps) {
        if (ratchetScore + epsilon < bestScore) {
          .Success(2L, "{.strong New best score}: {signif(ratchetScore)}")
          bestScore <- ratchetScore
          bestPlusEps <- bestScore + epsilon
          bestEdges <- .ReplaceResults(bestEdges, ratchetImproved,
                                       1 + tbrStart + iter)
          edge <- ratchetImproved[, , sample.int(dim(ratchetImproved)[3], 1)]
        } else {
          .Info(3L, "Hit best score {.strong {signif(bestScore)}} again")

          edge <- ratchetImproved[, , sample.int(dim(ratchetImproved)[3], 1)]
          bestEdges <- .CombineResults(bestEdges, ratchetImproved,
                                       1 + tbrStart + iter)
        }
      } else {
        if (3L < verbosity) {
          cli_alert_danger("Did not hit best score {signif(bestScore)}")
        }
      }
      if (.Timeout()) {
        return(.ReturnValue(bestEdges))                                         # nocov
      }
    }
  }
  
  # Branch breaking
  if (tbrEnd) {
    searchIter <- tbrIter * finalIter
    searchHits <- maxHits
    
    .Heading("Sample local optimum",
             "TBR depth {searchIter}; keeping {searchHits}",
             " trees; k = {concavity}")
    .Info(1L, Sys.time(), ": Score: ", signif(bestScore))
    finalEdges <- .Search('Final search')
    newBestScore <- .Score(finalEdges[, , 1])
    improved <- newBestScore + epsilon < bestScore
    bestEdges <- if (improved) {
      .ReplaceResults(bestEdges, finalEdges, 1 + tbrStart + ratchIter + 1)
    } else {
      .CombineResults(bestEdges, finalEdges, 1 + tbrStart + ratchIter + 1)
    }
  }
  
  # Return:
  .ReturnValue(bestEdges)
}

#' Combine two edge matrices
#' 
#' @param x,y 3D arrays, each slice containing an edge matrix from a tree
#' of class `phylo`.  `x` should not contain duplicates.
#' @return A single 3D array containing each unique edge matrix from (`x` and)
#' `y`, with a `firstHit` attribute as documented in [`MaximizeParsimony()`].
#' @template MRS
#' @keywords internal
.CombineResults <- function (x, y, stage) {
  xDim <- dim(x)
  if (length(xDim) == 2L) {
    xDim <- c(xDim, 1L)
  }
  if (any(duplicated(x, MARGIN = 3L))) {
    warning(".CombineResults(x) should not contain duplicates.")
  }
  
  res <- unique(array(c(x, y), dim = xDim + c(0, 0, dim(y)[3])), MARGIN = 3L)
  firstHit <- attr(x, 'firstHit')
  firstHit[stage] <- dim(res)[3] - xDim[3]
  attr(res, 'firstHit') <- firstHit
  
  # Return:
  res
}

#' @rdname dot-CombineResults
#' @param old old array of edge matrices with `firstHit` attribute.
#' @param new new array of edge matrices.
#' @param stage Integer specifying element of `firstHit` in which new hits
#' should be recorded.
#' @keywords internal
.ReplaceResults <- function (old, new, stage) {
  hit <- attr(old, 'firstHit')
  hit[] <- 0
  hit[stage] <- dim(new)[3]
  structure(new, 'firstHit' = hit)
}

#' @rdname MaximizeParsimony
#' 
#' @param method Unambiguous abbreviation of `jackknife` or `bootstrap` 
#' specifying how to resample characters.  Note that jackknife is considered
#' to give more meaningful results.
#' 
#' @param proportion Numeric between 0 and 1 specifying what proportion of 
#' characters to retain under jackknife resampling.
#' 
#' @section Resampling:
#' Note that bootstrap support is a measure of the amount of data supporting
#' a split, rather than the amount of confidence that should be afforded the
#' grouping.
#' "Bootstrap support of 100% is not enough, the tree must also be correct" 
#' \insertCite{Phillips2004}{TreeSearch}.
#' See discussion in \insertCite{Egan2006;textual}{TreeSearch};
#' \insertCite{Wagele2009;textual}{TreeSearch};
#' \insertCite{Simmons2011}{TreeSearch};
#' \insertCite{Kumar2012;textual}{TreeSearch}.
#' 
#' For a discussion of suitable search parameters in resampling estimates, see
#' \insertCite{Muller2005;textual}{TreeSearch}.
#' The user should decide whether to start each resampling
#' from the optimal tree (which may be quicker, but result in overestimated 
#' support values as searches get stuck in local optima close to the 
#' optimal tree) or a random tree (which may take longer as more rearrangements
#' are necessary to find an optimal tree on each iteration).
#' 
#' For other ways to estimate clade concordance, see [`SiteConcordance()`].
#' 
#' @return `Resample()` returns a `multiPhylo` object containing a list of
#' trees obtained by tree search using a resampled version of `dataset`.
#' @family split support functions
#' @encoding UTF-8
#' @export
Resample <- function (dataset, tree, method = 'jack',
                      proportion = 2/3,
                      ratchIter = 1L, tbrIter = 8L, finalIter = 3L,
                      maxHits = 12L, concavity = Inf,
                      tolerance = sqrt(.Machine$double.eps),
                      constraint,
                      verbosity = 2L,
                      ...) {
  if (!inherits(dataset, 'phyDat')) {
    stop("`dataset` must be of class `phyDat`.")
  }
  index <- attr(dataset, 'index')
  kept <- switch(pmatch(tolower(method), c('jackknife', 'bootstrap')),
         {
           nKept <- ceiling(proportion * length(index))
           if (nKept < 1L) {
             stop("No characters retained. `proportion` must be positive.")
           }
           if (nKept == length(index)) {
             stop("`proportion` too high; no characters deleted.")
           }
           sample(index, nKept)
         }, {
           sample(index, length(index), replace = TRUE)
         })
  if (is.null(kept)) {
    stop("`method` must be either 'jackknife' or 'bootstrap'.")
  }
  
  attr(dataset, 'index') <- kept
  attr(dataset, 'weight') <- vapply(seq_len(attr(dataset, 'nr')),
                                    function (x) sum(kept == x),
                                    integer(1))
  
  MaximizeParsimony(dataset, tree = tree,
                    ratchIter = ratchIter, tbrIter = tbrIter,
                    finalIter = finalIter,
                    maxHits = maxHits,
                    concavity = concavity,
                    tolerance = tolerance, constraint = constraint,
                    verbosity = verbosity, ...) 
}

#' Launch tree search graphical user interface
#' 
#' @rdname MaximizeParsimony
#' @importFrom cluster pam silhouette
#' @importFrom future future
#' @importFrom promises future_promise
#' @importFrom protoclust protoclust
#' @importFrom Rogue ColByStability
#' @importFrom shiny runApp
#' @importFrom shinyjs useShinyjs
#' @importFrom TreeDist ClusteringInfoDistance
#' @export
EasyTrees <- function () {#nocov start
  shiny::runApp(system.file('Parsimony', package = 'TreeSearch'))
}

#' @rdname MaximizeParsimony
#' @export
EasyTreesy <- EasyTrees
#nocov end

.UseProfile <- function (concavity) {
  pmatch(tolower(concavity), "profile", -1L) == 1L
}
