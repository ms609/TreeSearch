cat("Loading library...\n")
library(TreeSearch)

cat("Testing simple C++ call...\n")
# Try calling a simple exported C++ function
#result <- .Call("_TreeSearch_ts_fitch_score")
#cat("Result:", result, "\n")

cat("Calling loaded package check...\n")
cat("Package functions:\n")
ls(getNamespace("TreeSearch"))[grep("ts_", ls(getNamespace("TreeSearch")))] |> head(10) |> print()

cat("Done!\n")
