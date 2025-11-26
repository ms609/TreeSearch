library("Rcpp", quietly = TRUE)
library("testthat", quietly = TRUE)
library("TreeTools", quietly = TRUE)
library("gtools", quietly = TRUE) # Used to generate all permutations

# --- Enumeration Helper Function (CORRECTED) ---

# Note: Keeping the helper function logic simple and explicit for testing
calculate_parsimony <- function(tree, character_vector) {
  # (Fitch algorithm implementation from previous step, assumed correct)
  names(character_vector) <- tree$tip.label
  nodes_postorder <- rev(tree$edge[, 1])
  n_tips <- length(tree$tip.label)
  n_all <- tree$Nnode + n_tips
  node_states <- vector("list", n_all)
  
  tip_states <- character_vector[as.character(1:n_tips)]
  for (i in 1:n_tips) {
    node_states[[i]] <- as.integer(tip_states[i]) 
  }
  
  score <- 0
  
  for (parent_node in unique(nodes_postorder)) {
    children <- tree$edge[tree$edge[, 1] == parent_node, 2]
    set_L <- node_states[[children[1]]]
    set_R <- node_states[[children[2]]]
    set_P <- intersect(set_L, set_R)
    
    if (length(set_P) == 0) {
      node_states[[parent_node]] <- union(set_L, set_R)
      score <- score + 1
    } else {
      node_states[[parent_node]] <- set_P
    }
  }
  return(score)
}


enumerate_fixed_tree_scores <- function(tree, states) {
  n_tips <- length(tree$tip.label)
  
  state_indices <- which(states > 0)
  state_tokens <- state_indices - 1 
  
  base_labels <- integer(0)
  for (i in seq_along(state_tokens)) {
    base_labels <- c(base_labels, rep(state_tokens[i], states[state_tokens[i] + 1]))
  }
  
  # CRITICAL FIX: Use 'unique()' to filter the 720 non-unique permutations down to 90.
  labelings_matrix <- unique(gtools::permutations(n_tips, n_tips, v = base_labels))
  
  scores <- apply(labelings_matrix, 1, function(labeling) {
    calculate_parsimony(tree, labeling)
  })
  
  score_counts <- table(scores)
  
  return(score_counts)
}

# -------------------------------------------------------------------------
# Test Case Setup: Balanced Tree (6 taxa)
# -------------------------------------------------------------------------
nTip <- 6
tree <- TreeTools::BalancedTree(nTip)
tree$tip.label <- as.character(1:nTip) 
edge <- tree$edge

tree_ptr <- FixedTree_Preprocess(edge, nTip)

# -------------------------------------------------------------------------
# Test 1: Core 3-Token Character (n0=2, n1=2, n2=2). N=6.
# -------------------------------------------------------------------------
states3_6 <- c(2, 2, 2) 
log_total3_6 <- LogMultinomial(states3_6) # Log(90)

# Calculate expected counts via enumeration (Now correct)
expected_counts3_6 <- enumerate_fixed_tree_scores(tree, states3_6)
# Expected values for BalancedTree(6): k=4: 40, k=5: 50.

test_that("FixedTree_Count 3-Token (2, 2, 2) matches Enumeration", {
  
  exp_k4 <- as.numeric(expected_counts3_6["4"])
  exp_k5 <- as.numeric(expected_counts3_6["5"])
  
  # Test k=4 (Min score)
  log_count_k4 <- FixedTree_Count(tree_ptr, steps=4, states=states3_6)
  expect_equal(exp(log_count_k4[1]), exp_k4, tolerance = 1e-6)
  
  # Test k=5
  log_count_k5 <- FixedTree_Count(tree_ptr, steps=5, states=states3_6)
  expect_equal(exp(log_count_k5[1]), exp_k5, tolerance = 1e-6)
  
  # Check that the sum of the counts equals the total number of labelings
  total_calculated <- log(exp(log_count_k4) + exp(log_count_k5))
  expect_equal(total_calculated, log_total3_6, tolerance = 1e-6)
})

# -------------------------------------------------------------------------
# Test 2: 2-Token Character (n0=3, n1=3). N=6. Matches Enumeration
# -------------------------------------------------------------------------
states2_6 <- c(3, 3) # n0=3, n1=3 (only two states used)
log_total2_6 <- LogMultinomial(states2_6) # Log(20)

# Calculate expected counts via enumeration (Now correct)
expected_counts2_6 <- enumerate_fixed_tree_scores(tree, states2_6)
# Expected values for BalancedTree(6): k=2 (20 ways)

test_that("FixedTree_Count 2-Token (3, 3) matches Enumeration", {
  
  exp_k2 <- as.numeric(expected_counts2_6["2"])
  
  # Test k=1 (Should be 0 labelings, min score is 2)
  log_count_k1 <- FixedTree_Count(tree_ptr, steps=1, states=states2_6)
  expect_true(is.infinite(log_count_k1[1]) && log_count_k1[1] < 0)
  
  # Test k=2 (Min/Max score)
  log_count_k2 <- FixedTree_Count(tree_ptr, steps=2, states=states2_6)
  expect_equal(exp(log_count_k2[1]), exp_k2, tolerance = 1e-6)
  
  # Check total matches Multinomial
  total_calculated <- log(exp(log_count_k2))
  expect_equal(total_calculated, log_total2_6, tolerance = 1e-6)
})
