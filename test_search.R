devtools::load_all(quiet = TRUE, compile = FALSE)
library(TreeTools)

# Replicate the failing test
true_tree <- ape::read.tree(text = "(((((1,2),3),4),5),6);")
dataset <- TreeTools::StringToPhyDat("110000 111000 111100", 1:6, byTaxon = FALSE)
start_tree <- TreeTools::RenumberTips(ape::read.tree(
  text = "(((1, 6), 3), (2, (4, 5)));"), true_tree$tip.label)

cat("Start tree score:", TreeLength(start_tree, dataset), "\n")
cat("True tree score:", TreeLength(true_tree, dataset), "\n")

# Try just a simple TreeSearch call
cat("\nTrying TreeSearch with SPRSwap...\n")
result <- TreeSearch(start_tree, dataset,
                    EdgeSwapper = SPRSwap,
                    maxIter = 100,
                    maxHits = 10,
                    verbosity = 1)
cat("TreeSearch result score:", attr(result, "score"), "\n")
