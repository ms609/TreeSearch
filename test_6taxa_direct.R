cat("Testing with 6 taxa (direct call)...\n")
library(TreeSearch)

# Use the same edge matrix as the failing test
edge <- matrix(c(7, 8, 8, 9, 9, 8, 7, 10, 10, 11,
                 8, 9, 9, 1, 6, 3, 10, 2, 11, 4), nrow = 10, ncol = 2)
cat("Edge matrix:\n")
print(edge)

contrast <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)
tip_data <- matrix(c(1,1,1,1,1,1,1,1,2,2,1,1,1,2,2,2,1,2,2,2), nrow = 6, ncol = 2, byrow = FALSE)

cat("Tip data:\n")
print(tip_data)

weight <- c(1, 1, 1)  # 3 characters total
levels <- c("0", "1")

cat("Calling ts_fitch_score...\n")
flush.console()
result <- TreeSearch:::ts_fitch_score(edge, contrast, tip_data, weight, levels)
cat("Result:", result, "\n")
