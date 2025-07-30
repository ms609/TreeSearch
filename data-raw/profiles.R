library("TreeTools", quietly = TRUE)
library("TreeSearch")
nTip <- 9

Profiles <- function (nTip) {
  if (nTip < 4) return (0)
  message(Sys.time(), ": Profiling ", nTip)
  smallest1 <- nTip
  smallest2 <- 2:(nTip / 2)
  two <- cbind(smallest = smallest2, largest = nTip - smallest2)
  nTrees <- NUnrooted(nTip)
  
  allEdges <- vapply(seq_len(nTrees) - 1L, function (x) 
    do.call(cbind,
            RenumberEdges(TreeTools:::num_to_parent(TreeTools:::.Int64.to.C(x), 
                                                    nTip),
                          seq_len(nTip + nTip - 2L))),
    matrix(0L, nTip + nTip - 2L, 2))
  
  MakeChar <- function (x) {
    SingleCharMorphy(rep.int(seq_along(x) - 1L, x))
  }
  Subtract <- function (x) x - log2(nTrees)
  
  twoChars <- apply(two, 1, MakeChar)
  on.exit(vapply(twoChars, UnloadMorphy, 0))
  
  twoScores <- apply(allEdges, 3, TreeSearch:::preorder_morphy_by_char,
                     twoChars)
  tab <- if (is.null(dim(twoScores))) {
    list(table(twoScores))
  } else {
    apply(twoScores, 1, table)
  }
  
  list(
    0,
    lapply(lapply(lapply(tab, cumsum), log2), Subtract)
  )
}

profiles <- lapply(1:10, Profiles) # Error w/ 11 tips: cannot allocate vector of size 5.1 Gb
usethis::use_data(profiles, overwrite = TRUE)
