devtools::load_all(quiet = TRUE, compile = FALSE)
library(TreeTools)

# Simple Ratchet test
true_tree <- ape::read.tree(text = "(((((1,2),3),4),5),6);")
dataset <- TreeTools::StringToPhyDat("110000 111000 111100", 1:6, byTaxon = FALSE)
start_tree <- TreeTools::RenumberTips(ape::read.tree(
  text = "(((1, 6), 3), (2, (4, 5)));"), true_tree$tip.label)

cat("TreeLength of start_tree:", TreeLength(start_tree, dataset), "\n")
cat("TreeLength of true_tree:", TreeLength(true_tree, dataset), "\n")

cat("\nCalling Ratchet...\n")
flush.console()

ratchetResult <- Ratchet(start_tree, dataset,
                         swappers = list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap),
                         ratchIter = 3, searchHits = 5, verbosity = 2)

cat("Ratchet completed!\n")
ratchetScore <- attr(ratchetResult, "score")
cat("Ratchet score:", ratchetScore, "\n")
