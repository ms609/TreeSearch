# `stopPatience` is a flat replicate patience: stop after N consecutive replicates fail to
# improve, with no reference to the hit count.  It exists because both other no-improvement
# rules are indexed on replicates via hits, so a change that makes a replicate individually
# better but slower delays the stop instead of improving the answer.  It must be inert at its
# default of 0, must never stop the search LATER than the other rules, and -- since it is
# shipped by `sprint`/`default` under implied weights only -- must leave equal weights alone.

testDataset <- function() {
  TreeTools::MatrixToPhyDat(rbind(
    a = c(0, 0, 0, 0, 1, 0), b = c(0, 0, 0, 1, 1, 1), c = c(1, 1, 0, 0, 0, 1),
    d = c(1, 1, 1, 0, 1, 0), e = c(1, 0, 1, 1, 0, 0), f = c(0, 1, 1, 1, 0, 1),
    g = c(1, 0, 0, 1, 1, 0), h = c(0, 1, 0, 0, 1, 1)
  ))
}

test_that("stopPatience is inert at its default", {
  dataset <- testDataset()
  reps <- 20L
  # Both count rules off, so only the replicate cap can end the search.
  r <- MaximizeParsimony(dataset, targetHits = 99999L, perturbStopFactor = 0L,
                         consensusStableReps = 0L, maxReplicates = reps, verbosity = 0L)
  expect_equal(attr(r, "replicates"), reps)
  expect_equal(SearchControl()$stopPatience, 0L)
})

test_that("stopPatience stops the search early and keeps the best score", {
  dataset <- testDataset()
  reps <- 200L
  args <- list(dataset = dataset, targetHits = 99999L, perturbStopFactor = 0L,
               consensusStableReps = 0L, maxReplicates = reps, verbosity = 0L)
  full <- do.call(MaximizeParsimony, args)
  short <- do.call(MaximizeParsimony, c(args, list(stopPatience = 2L)))
  expect_lt(attr(short, "replicates"), reps)
  # This dataset is tiny, so the optimum is found immediately: an early stop must not cost
  # score.  (A patience rule can only truncate; it cannot alter earlier replicates.)
  expect_equal(min(attr(short, "score")), min(attr(full, "score")))
  expect_true(attr(short, "perturb_stop"))
})

test_that("stopPatience fires at lastImprovement + patience", {
  # The identity the campaign's terminator audit relies on, and the reason a pre-registered
  # check of `replicates == patience + 1` was wrong: the counter RESETS on every improvement,
  # so the stop replicate is measured from the last improvement, not from the start.
  dataset <- testDataset()
  patience <- 4L
  r <- MaximizeParsimony(dataset, targetHits = 99999L, perturbStopFactor = 0L,
                         consensusStableReps = 0L, maxReplicates = 200L,
                         stopPatience = patience, verbosity = 0L)
  expect_equal(attr(r, "replicates"),
               attr(r, "last_improved_rep") + patience)
})

test_that("A larger stopPatience runs at least as long as a smaller one", {
  dataset <- testDataset()
  args <- list(dataset = dataset, targetHits = 99999L, perturbStopFactor = 0L,
               consensusStableReps = 0L, maxReplicates = 200L, verbosity = 0L)
  tight <- do.call(MaximizeParsimony, c(args, list(stopPatience = 2L)))
  loose <- do.call(MaximizeParsimony, c(args, list(stopPatience = 25L)))
  expect_lte(attr(tight, "replicates"), attr(loose, "replicates"))
})

test_that("stopPatience = 0 disables the rule, and a negative value errors", {
  dataset <- testDataset()
  reps <- 15L
  r <- MaximizeParsimony(dataset, targetHits = 99999L, perturbStopFactor = 0L,
                         consensusStableReps = 0L, maxReplicates = reps,
                         stopPatience = 0L, verbosity = 0L)
  expect_equal(attr(r, "replicates"), reps)
  # Silently treating -20 as "off" would swallow an obvious typo.
  expect_error(SearchControl(stopPatience = -5L), "non-negative")
  expect_error(SearchControl(stopPatience = c(1L, 2L)), "single")
  expect_error(SearchControl(stopPatience = NA_integer_), "non-negative")
  # A non-numeric value reaches the guard as NA via as.integer(), with a warning.
  expect_error(suppressWarnings(SearchControl(stopPatience = "twenty")),
               "non-negative")
})

test_that("stopPatience survives the SearchControl round trip", {
  # The C++ side reads this field only `if (ctrl.containsElementNamed(...))`, so a field added
  # to SearchControl() but dropped anywhere in the plumbing fails SILENTLY as "patience off"
  # rather than erroring.  Assert the value arrives, not merely that it is accepted.
  ctrl <- SearchControl(stopPatience = 7L)
  expect_equal(ctrl$stopPatience, 7L)
  expect_true("stopPatience" %in% attr(ctrl, "explicit"))
  dataset <- testDataset()
  r <- MaximizeParsimony(dataset, targetHits = 99999L, perturbStopFactor = 0L,
                         consensusStableReps = 0L, maxReplicates = 200L,
                         control = ctrl, verbosity = 0L)
  expect_equal(attr(r, "replicates"), attr(r, "last_improved_rep") + 7L)
})

test_that("stopPatience also stops the parallel search", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  # The parallel path has its OWN implementation of this rule: it counts a dry spell over
  # replicates completed into the shared pool, not the serial `unsuccessful_reps`, so the
  # `lastImprovement + patience` identity above does NOT transfer -- the firing replicate
  # varies with `nThreads` (patience 5 on this matrix: 6 replicates serial, 21 on two
  # threads, 65 on four).  Assert only what holds on both paths.
  #
  # `testDataset()` is deliberately NOT used here.  Its replicates are so cheap that all 120
  # finish inside the parallel monitor's first 200 ms sleep, so the loop's next act is to
  # observe `replicates_done >= max_replicates` and break -- the rule is never evaluated and
  # the search runs to the cap.  That polling granularity is shared by every stopping rule on
  # this path, not specific to `stopPatience`; it needs replicates slow enough for a poll to
  # observe the dry spell.  Hence a 40-tip matrix under implied weights.
  #
  # Worth exercising for a second reason: nThreads >= 2 has a history of MinGW heap
  # corruption here (PR #258), and this change restructured a guard in the replicate loop.
  set.seed(2)
  dataset <- TreeTools::MatrixToPhyDat(
    matrix(sample(0:1, 40 * 30, TRUE), nrow = 40,
           dimnames = list(paste0("t", seq_len(40)), NULL))
  )
  cap <- 120L
  # `.rung = "none"` is load-bearing: 40 tips and 30 characters make `auto` resolve to
  # `default`, which under implied weights now ships stopPatience 15 -- so a control arm left
  # on `auto` stops at ~21 replicates and the comparison measures nothing.
  args <- list(dataset = dataset, .rung = "none", concavity = 10,
               targetHits = 99999L, perturbStopFactor = 0L,
               consensusStableReps = 0L, maxReplicates = cap, nThreads = 2L,
               verbosity = 0L)
  full <- do.call(MaximizeParsimony, args)
  short <- do.call(MaximizeParsimony, c(args, list(stopPatience = 3L)))
  expect_equal(attr(full, "replicates"), cap)      # nothing else ends the search
  expect_lt(attr(short, "replicates"), cap)
  expect_true(attr(short, "perturb_stop"))
  # Deliberately NOT `expect_equal(short score, full score)`.  That assertion held on
  # this matrix but is not a property the rule guarantees, and it broke under covr:
  # instrumentation slows every replicate, the parallel path evaluates the dry spell on a
  # 200 ms monitor poll over replicates completed into the shared pool, so far fewer
  # replicates land per poll and patience 3 fires genuinely earlier in the search
  # (13.063 against 12.988).  Stopping sooner is *allowed* to cost score -- that is the
  # trade-off the parameter exists to offer -- so equality was asserting a coincidence of
  # this machine's timing.  What the rule does guarantee is asserted above: the search
  # stops before the cap, and it stops for the no-improvement reason.  The score is only
  # checked for being a valid, finite improvement over a random start, which holds at any
  # speed.  See [[parallel-stop-rules-poll-granularity]].
  # What the rule actually
  # guarantees is that stopping early can't do BETTER than letting the search run on.
  expect_true(is.finite(min(attr(short, "score"))))
  expect_lt(min(attr(short, "score")), min(attr(full, "score")) * 1.5)
  expect_gte(min(attr(short, "score")), min(attr(full, "score")))
})

# ---- the shipped implied-weights operating point --------------------------------------------
# `sprint`/`default` take a deeper ratchet paid for by this patience, under implied weights
# ONLY (2026-07-28, 4624 cells).  Equal weights and profile parsimony were never measured, so
# they must come through untouched -- these tests pin the SCOPE, which is the part a later
# refactor is most likely to break silently.

test_that(".IwStopPackage applies to sprint and default under implied weights", {
  expect_equal(.IwStopPackage("sprint", 10),
               list(ratchetCycles = 12L, ratchetPerturbProb = 0.25,
                    stopPatience = 20L))
  expect_equal(.IwStopPackage("default", 10),
               list(ratchetCycles = 20L, stopPatience = 15L))
})

test_that(".IwStopPackage leaves equal weights and profile parsimony alone", {
  expect_null(.IwStopPackage("sprint", Inf))
  expect_null(.IwStopPackage("default", Inf))
  expect_null(.IwStopPackage("sprint", "profile"))
  expect_null(.IwStopPackage("default", "profile"))
})

test_that(".IwStopPackage does not fire for the strategies it never measured", {
  # Disjoint from .IwRatchetDepth() by strategy, so the two can never both set ratchetCycles.
  for (strategy in c("thorough", "large", "intensive", "none")) {
    expect_null(.IwStopPackage(strategy, 10))
  }
  expect_null(.IwStopPackage(character(0), 10))
})

test_that(".IwStopPackage never overrides a field the caller set", {
  expect_equal(.IwStopPackage("sprint", 10, userSet = "ratchetCycles"),
               list(ratchetPerturbProb = 0.25, stopPatience = 20L))
  expect_equal(.IwStopPackage("default", 10, userSet = "stopPatience"),
               list(ratchetCycles = 20L))
  expect_null(.IwStopPackage("default", 10,
                             userSet = c("ratchetCycles", "stopPatience")))
})

test_that("the implied-weights package reaches a real search, including via auto", {
  dataset <- testDataset()                      # 8 tips -> auto resolves to sprint
  expect_equal(.AutoRung(8L, 6L), 1L)
  # A user-set value must win over the package even when the strategy would impose one.
  r <- MaximizeParsimony(dataset, effort = 0L, concavity = 10,
                         maxReplicates = 200L, targetHits = 99999L,
                         perturbStopFactor = 0L, consensusStableReps = 0L,
                         stopPatience = 3L, verbosity = 0L)
  expect_equal(attr(r, "replicates"), attr(r, "last_improved_rep") + 3L)
  # Left to itself, `auto` under implied weights takes sprint's patience of 20.
  auto <- MaximizeParsimony(dataset, effort = 0L, concavity = 10,
                            maxReplicates = 500L, targetHits = 99999L,
                            perturbStopFactor = 0L, consensusStableReps = 0L,
                            verbosity = 0L)
  expect_equal(attr(auto, "replicates"), attr(auto, "last_improved_rep") + 20L)
  # Equal weights is out of scope, so nothing stops the search but the cap.
  ew <- MaximizeParsimony(dataset, effort = 0L, maxReplicates = 30L,
                          targetHits = 99999L, perturbStopFactor = 0L,
                          consensusStableReps = 0L, verbosity = 0L)
  expect_equal(attr(ew, "replicates"), 30L)
})
