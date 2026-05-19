cat("Testing with 3 taxa...\n")
library(TreeSearch)

# 3-taxon tree: ((1,2),3)
edge <- matrix(c(4, 4, 5, 5, 1, 2, 4, 3), nrow = 4, ncol = 2)
cat("Edge matrix:\n")
print(edge)

contrast <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)
tip_data <- matrix(c(1, 1, 2, 1, 1, 1), nrow = 3, ncol = 2)
weight <- c(1, 1)
levels <- c("0", "1")

cat("Calling ts_fitch_score...\n")
flush.console()
result <- TreeSearch:::ts_fitch_score(edge, contrast, tip_data, weight, levels)
cat("Result:", result, "\n")
