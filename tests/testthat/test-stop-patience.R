# TS_STOP_PATIENCE is an experimental, opt-in flat replicate patience: stop after N
# consecutive replicates fail to improve, with no reference to the hit count.  It exists
# because both shipped rules are indexed on replicates via hits, so a change that makes a
# replicate individually better but slower delays the stop instead of improving the answer.
# It must be inert unless set, and must never stop the search LATER than the shipped rules.

withPatience <- function(value, code) {
  old <- Sys.getenv("TS_STOP_PATIENCE", unset = NA)
  if (is.na(old)) {
    on.exit(Sys.unsetenv("TS_STOP_PATIENCE"), add = TRUE)
  } else {
    on.exit(Sys.setenv(TS_STOP_PATIENCE = old), add = TRUE)
  }
  Sys.setenv(TS_STOP_PATIENCE = value)
  force(code)
}

testDataset <- function() {
  TreeTools::MatrixToPhyDat(rbind(
    a = c(0, 0, 0, 0, 1, 0), b = c(0, 0, 0, 1, 1, 1), c = c(1, 1, 0, 0, 0, 1),
    d = c(1, 1, 1, 0, 1, 0), e = c(1, 0, 1, 1, 0, 0), f = c(0, 1, 1, 1, 0, 1),
    g = c(1, 0, 0, 1, 1, 0), h = c(0, 1, 0, 0, 1, 1)
  ))
}

test_that("TS_STOP_PATIENCE is inert when unset", {
  skip_if(nzchar(Sys.getenv("TS_STOP_PATIENCE")),
          "TS_STOP_PATIENCE is set in this environment")
  dataset <- testDataset()
  reps <- 20L
  # Both count rules off, so only the replicate cap can end the search.
  r <- MaximizeParsimony(dataset, targetHits = 99999L, perturbStopFactor = 0L,
                         consensusStableReps = 0L, maxReplicates = reps, verbosity = 0L)
  expect_equal(attr(r, "replicates"), reps)
})

test_that("A flat patience stops the search early and keeps the best score", {
  dataset <- testDataset()
  reps <- 200L
  args <- list(dataset = dataset, targetHits = 99999L, perturbStopFactor = 0L,
               consensusStableReps = 0L, maxReplicates = reps, verbosity = 0L)
  full <- do.call(MaximizeParsimony, args)
  short <- withPatience("2", do.call(MaximizeParsimony, args))
  expect_lt(attr(short, "replicates"), reps)
  # This dataset is tiny, so the optimum is found immediately: an early stop must not cost
  # score.  (A patience rule can only truncate; it cannot alter earlier replicates.)
  expect_equal(min(attr(short, "score")), min(attr(full, "score")))
})

test_that("A larger patience runs at least as long as a smaller one", {
  dataset <- testDataset()
  args <- list(dataset = dataset, targetHits = 99999L, perturbStopFactor = 0L,
               consensusStableReps = 0L, maxReplicates = 200L, verbosity = 0L)
  tight <- withPatience("2", do.call(MaximizeParsimony, args))
  loose <- withPatience("25", do.call(MaximizeParsimony, args))
  expect_lte(attr(tight, "replicates"), attr(loose, "replicates"))
})

test_that("A non-positive or unparseable patience is ignored", {
  dataset <- testDataset()
  reps <- 15L
  args <- list(dataset = dataset, targetHits = 99999L, perturbStopFactor = 0L,
               consensusStableReps = 0L, maxReplicates = reps, verbosity = 0L)
  for (bad in c("0", "-5", "not-a-number", "")) {
    r <- withPatience(bad, do.call(MaximizeParsimony, args))
    expect_equal(attr(r, "replicates"), reps)
  }
})
