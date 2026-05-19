cat("Testing with 6 taxa and 3 characters...\n")
library(TreeSearch)

# Use the same edge matrix as the failing test
edge <- matrix(c(7, 8, 8, 9, 9, 8, 7, 10, 10, 11,
                 8, 9, 9, 1, 6, 3, 10, 2, 11, 4), nrow = 10, ncol = 2)

contrast <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)

# 6 taxa x 3 characters (matching the test data: "110000 111000 111100")
#  char1: 1 1 1 1 1 1
#  char2: 1 1 1 2 2 2
#  char3: 0 0 1 1 2 2 (but scaled to 1-based: 1 1 2 2 3 3)
# Actually let me make it binary matching the input
# "110000" = [1,1,0,0,0,0] -> [1,1,1,1,1,1] (0->1, 1->2)
# "111000" = [1,1,1,0,0,0] -> [2,2,2,1,1,1]
# "111100" = [1,1,1,1,0,0] -> [2,2,2,2,1,1]

tip_data <- matrix(c(1, 1, 1, 1, 1, 1,    # char 1
                     2, 2, 2, 1, 1, 1,    # char 2
                     2, 2, 2, 2, 1, 1),   # char 3
                   nrow = 6, ncol = 3, byrow = FALSE)

cat("Tip data (6x3):\n")
print(tip_data)

weight <- c(1, 1, 1)
levels <- c("0", "1")

cat("Calling ts_fitch_score...\n")
flush.console()
result <- TreeSearch:::ts_fitch_score(edge, contrast, tip_data, weight, levels)
cat("Result:", result, "\n")
