cat("Starting test...\n")
devtools::load_all(quiet = TRUE, compile = FALSE)
cat("Package loaded\n")
library(TreeTools)
cat("TreeTools loaded\n")

# Replicate the failing test
true_tree <- ape::read.tree(text = "(((((1,2),3),4),5),6);")
cat("true_tree created\n")
dataset <- TreeTools::StringToPhyDat("110000 111000 111100", 1:6, byTaxon = FALSE)
cat("dataset created\n")
start_tree <- TreeTools::RenumberTips(ape::read.tree(
  text = "(((1, 6), 3), (2, (4, 5)));"), true_tree$tip.label)
cat("start_tree created\n")

cat("Start tree score:", TreeLength(start_tree, dataset), "\n")
cat("Done\n")
