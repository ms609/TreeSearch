library(TreeTools)
library(TreeSearch)

true_tree <- ape::read.tree(text = "(((((1,2),3),4),5),6);")
dataset <- TreeTools::StringToPhyDat("110000 111000 111100", 1:6, byTaxon = FALSE)
start_tree <- TreeTools::RenumberTips(ape::read.tree(
  text = "(((1, 6), 3), (2, (4, 5)));"), true_tree$tip.label)

cat("First TreeLength call...\n")
flush.console()
result1 <- TreeLength(start_tree, dataset)
cat("Result 1:", result1, "\n")

cat("Second TreeLength call...\n")
flush.console()
result2 <- TreeLength(true_tree, dataset)
cat("Result 2:", result2, "\n")

cat("Done!\n")
