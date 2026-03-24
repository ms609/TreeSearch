context("Strategy tracker (Thompson sampling bandit)")

# Exercises StrategyTracker via the ts_test_strategy_tracker() bridge.
# The bridge creates a fresh tracker, draws `n_draws` selections, then
# applies 5 successes to arm 0 (WAGNER_RANDOM) and 5 failures to arm 1
# (WAGNER_GOLOBOFF), decays by 0.5, and draws again. It also produces
# a round-robin sequence.

test_that("Strategy tracker has 4 fresh-start arms with correct names", {
  res <- TreeSearch:::ts_test_strategy_tracker(42L, 100L)
  expect_equal(res$n_strategies, 4L)
  expect_equal(
    res$strategy_names,
    c("wag_rand", "wag_golob", "wag_entropy", "rand_tree")
  )
})

test_that("Initial priors: Beta(1,1) for most, Beta(1,2) for RANDOM_TREE", {
  res <- TreeSearch:::ts_test_strategy_tracker(1L, 1L)
  expect_equal(res$alpha_init, c(1, 1, 1, 1))
  # RANDOM_TREE (index 4) has beta=2
  expect_equal(res$beta_init, c(1, 1, 1, 2))
})

test_that("All 4 arms selected over many draws", {
  res <- TreeSearch:::ts_test_strategy_tracker(8371L, 10000L)
  expect_true(all(res$initial_counts > 0))
  expect_equal(sum(res$initial_counts), 10000L)
})

test_that("RANDOM_TREE selected less often due to pessimistic prior", {
  # Over many draws from fresh priors, RANDOM_TREE (Beta(1,2)) should
  # average ~33% vs ~50% for the Beta(1,1) arms
  res <- TreeSearch:::ts_test_strategy_tracker(2847L, 50000L)
  # 4 arms: 3 with Beta(1,1) and 1 with Beta(1,2)
  # RANDOM_TREE should be selected notably less often than the others.
  rand_tree_frac <- res$initial_counts[4] / 50000
  other_mean_frac <- mean(res$initial_counts[1:3]) / 50000
  expect_true(rand_tree_frac < other_mean_frac)
})

test_that("Update modifies alpha/beta correctly", {
  res <- TreeSearch:::ts_test_strategy_tracker(1L, 1L)
  # arm 0 (WAGNER_RANDOM): 5 successes -> alpha += 5
  expect_equal(res$alpha_after_update[1], 1 + 5)
  expect_equal(res$beta_after_update[1], 1)  # no failures added
  # arm 1 (WAGNER_GOLOBOFF): 5 failures -> beta += 5
  expect_equal(res$alpha_after_update[2], 1)  # no successes added
  expect_equal(res$beta_after_update[2], 1 + 5)
  # arm 2 (WAGNER_ENTROPY): unchanged
  expect_equal(res$alpha_after_update[3], 1)
  expect_equal(res$beta_after_update[3], 1)
})

test_that("Decay halves excess over prior", {
  res <- TreeSearch:::ts_test_strategy_tracker(1L, 1L)
  # After update: arm 0 alpha=6, beta=1. Decay 0.5:
  #   alpha = max(1, 1 + (6-1)*0.5) = 1 + 2.5 = 3.5
  #   beta  = max(1, 1 + (1-1)*0.5) = 1
  expect_equal(res$alpha_after_decay[1], 3.5)
  expect_equal(res$beta_after_decay[1], 1.0)
  # arm 1: alpha=1, beta=6 -> beta = 1 + (6-1)*0.5 = 3.5
  expect_equal(res$alpha_after_decay[2], 1.0)
  expect_equal(res$beta_after_decay[2], 3.5)
  # RANDOM_TREE (arm 4): alpha=1 beta=2 -> beta = 1 + (2-1)*0.5 = 1.5
  expect_equal(res$alpha_after_decay[4], 1.0)
  expect_equal(res$beta_after_decay[4], 1.5)
})

test_that("Biased counts favour arm 0 after 5 successes (post-decay)", {
  # After arm 0 gets 5 successes + decay -> alpha=3.5, beta=1
  # Mean = 3.5/(3.5+1) = 0.78 — should dominate Thompson sampling
  res <- TreeSearch:::ts_test_strategy_tracker(5319L, 10000L)
  # arm 0 should be most-selected
  expect_equal(which.max(res$biased_counts), 1L)
  # arm 0 should get majority of draws
  expect_true(res$biased_counts[1] > 3000)
})

test_that("Round-robin sequence cycles through 4 arms", {
  res <- TreeSearch:::ts_test_strategy_tracker(1L, 1L)
  rr <- res$round_robin
  expect_length(rr, 12)
  # Cycles 0,1,2,3, 0,1,2,3, 0,1,2,3
  expect_equal(rr, rep(0:3, 3))
})

test_that("Adaptive search returns strategy_diagnostics attribute", {
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]
  set.seed(7021)
  res <- MaximizeParsimony(
    ds, maxReplicates = 4L, targetHits = 2L,
    adaptiveStart = TRUE, verbosity = 0L
  )
  diag <- attr(res, "strategy_diagnostics")
  expect_type(diag, "list")
  expect_named(diag, c("attempts", "successes"))
  expect_length(diag$attempts, 4L)
  # Total attempts should equal replicates completed
  expect_equal(sum(diag$attempts), attr(res, "replicates"))
  # At least one strategy should have been attempted
  expect_true(any(diag$attempts > 0))
})
