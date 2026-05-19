cat("Loading...\n")
devtools::load_all(quiet = TRUE, compile = FALSE)
library(TreeTools)

cat("Creating test data...\n")
true_tree <- ape::read.tree(text = "(((((1,2),3),4),5),6);")
dataset <- TreeTools::StringToPhyDat("110000 111000 111100", 1:6, byTaxon = FALSE)
cat("Dataset created successfully\n")

cat("Creating start tree...\n")
start_tree <- TreeTools::RenumberTips(ape::read.tree(
  text = "(((1, 6), 3), (2, (4, 5)));"), true_tree$tip.label)
cat("Start tree created successfully\n")

cat("Checking if tree is rooted...\n")
is_rooted <- ape::is.rooted(start_tree)
cat("Tree is rooted:", is_rooted, "\n")

cat("Checking tree structure...\n")
cat("Number of tips:", length(start_tree$tip.label), "\n")
cat("Number of nodes:", start_tree$Nnode, "\n")
cat("Edge matrix:\n")
print(start_tree$edge)

cat("Done!\n")
