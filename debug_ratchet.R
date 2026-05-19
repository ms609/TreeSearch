devtools::load_all(quiet = TRUE)
library(TreeTools)

# Replicate the failing test
true_tree <- ape::read.tree(text = "(((((1,2),3),4),5),6);")
dataset <- TreeTools::StringToPhyDat("110000 111000 111100", 1:6, byTaxon = FALSE)
start_tree <- TreeTools::RenumberTips(ape::read.tree(
  text = "(((1, 6), 3), (2, (4, 5)));"), true_tree$tip.label)

cat("Start tree score:", TreeLength(start_tree, dataset), "\n")
cat("True tree score:", TreeLength(true_tree, dataset), "\n")

# Run Ratchet
set.seed(42)  # For reproducibility
ratchetResult <- Ratchet(start_tree, dataset,
                swappers = list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap),
                ratchIter = 3, searchHits = 5, verbosity = 2)

ratchetScore <- attr(ratchetResult, "score")
cat("Ratchet score:", ratchetScore, "\n")
cat("Expected score:", TreeLength(true_tree, dataset), "\n")
cat("Test passes:", isTRUE(all.equal(TreeLength(true_tree, dataset), ratchetScore)), "\n")
