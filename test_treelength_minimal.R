cat("Testing TreeLength directly...\n")
library(TreeTools)
library(TreeSearch)

dataset <- TreeTools::StringToPhyDat("110000 111000 111100", 1:6, byTaxon = FALSE)
start_tree <- TreeTools::RenumberTips(ape::read.tree(
  text = "(((1, 6), 3), (2, (4, 5)));"), c("1", "2", "3", "4", "5", "6"))

cat("Calling TreeLength...\n")
flush.console()
result <- TreeLength(start_tree, dataset)
cat("Result:", result, "\n")
