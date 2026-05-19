cat("Loading...\n")
devtools::load_all(quiet = TRUE, compile = FALSE)
library(TreeTools)

cat("Creating test data...\n")
true_tree <- ape::read.tree(text = "(((((1,2),3),4),5),6);")
dataset <- TreeTools::StringToPhyDat("110000 111000 111100", 1:6, byTaxon = FALSE)
start_tree <- TreeTools::RenumberTips(ape::read.tree(
  text = "(((1, 6), 3), (2, (4, 5)));"), true_tree$tip.label)

cat("Dataset tip labels:", names(dataset), "\n")
cat("Start tree tip labels:", start_tree$tip.label, "\n")

cat("Calculating TreeLength for start_tree...\n")
flush.console()
result <- TreeLength(start_tree, dataset)
cat("Result:", result, "\n")
