# Extracted from test-Concordance.R:159

# prequel ----------------------------------------------------------------------
library("TreeTools", quietly = TRUE)

# test -------------------------------------------------------------------------
tree <- BalancedTree(8)
splits <- as.Splits(tree)
mataset <- matrix(c(0, 0, 0, 0, 0, 0, 0, 1,
                      rep("?", 8)), 8,
                    dimnames = list(paste0("t", 1:8), NULL))
dat <- MatrixToPhyDat(mataset)
expect_equal(unname(ClusteringConcordance(tree, dat)), rep(NA_real_, 5))
tree <- ape::read.tree(text = "((a, b, c, d, e), (f, g, h));")
split <- as.Splits(tree)
mataset <- matrix(c(0, 0, 0, 0, 0, 0, 0, 1,
                      0, 0, 0, 0, 0, 1, 1, 1, # Matches split
                      0, 0, 0, 0, 1, 1, 1, 1, # Consistent but not identical
                      0, 0, 0, 1, 1, 1, 1, 1, # Consistent, more different
                      0, 0, 0, 0, 0, 0, 1, 1, # Consistent other way
                      0, 1, 0, 1, 0, 1, 0, 1, # Worst possible
                      0, 0, 0, 0, rep("?", 4), # No information
                      0, 0, 1, 1, rep("?", 4), # No relevant information
                      rep("?", 8)), 8,
                    dimnames = list(letters[1:8], NULL))
dat <- MatrixToPhyDat(mataset)
cc <- ClusteringConcordance(tree, dat, return = "all")[, "10", ]
.Entropy <- function(...) {
    TreeDist::Entropy(c(...) / sum(...))
  }
.NormExp <- function(a, b, ab) {
    .Rezero(
      (.Entropy(a) + .Entropy(b) - .Entropy(ab)) / .Entropy(a),
      .ExpectedMI(a, b) / .Entropy(a)
    )
  }
