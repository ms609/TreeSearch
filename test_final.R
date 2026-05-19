devtools::load_all(quiet = TRUE, compile = FALSE)
library(TreeTools)

# Exact replica of the failing test
true_tree <- ape::read.tree(text = "(((((1,2),3),4),5),6);")
dataset <- TreeTools::StringToPhyDat("110000 111000 111100", 1:6, byTaxon = FALSE)
start_tree <- TreeTools::RenumberTips(ape::read.tree(
  text = "(((1, 6), 3), (2, (4, 5)));"), true_tree$tip.label)

cat("Test assertions:\n")
cat("1. TreeLength(start_tree, dataset) == 6:",
    TreeLength(start_tree, dataset) == 6, "\n")
cat("2. TreeLength(true_tree, dataset) == 3:",
    TreeLength(true_tree, dataset) == 3, "\n")

# The failing test
cat("\nRunning Ratchet...\n")
ratchetScore <- attr(Ratchet(start_tree, dataset,
                swappers = list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap),
                ratchIter = 3, searchHits = 5, verbosity = 0), "score")

cat("Ratchet score:", ratchetScore, "\n")
cat("Expected (TreeLength of true_tree):", TreeLength(true_tree, dataset), "\n")
cat("Test passes:", TreeLength(true_tree, dataset) == ratchetScore, "\n")
