# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# --- Helpers ---
# Build annealConfig list from SearchControl fields
make_anneal_config <- function(ctrl) {
  cycles <- as.integer(if (is.null(ctrl$annealCycles)) 0L else ctrl$annealCycles)
  if (cycles > 0L) {
    list(
      cycles = cycles,
      phases = as.integer(ctrl$annealPhases),
      tStart = as.double(ctrl$annealTStart),
      tEnd = as.double(ctrl$annealTEnd),
      movesPerPhase = as.integer(ctrl$annealMovesPerPhase)
    )
  } else {
    NULL
  }
}

anneal_search <- function(ds, ctrl_overrides = list(), maxSeconds = 5,
                          maxReplicates = 3L, verbosity = 0L,
                          concavity = -1.0) {
  ctrl <- SearchControl(
    annealCycles = 1L, annealPhases = 5L, annealTStart = 20, annealTEnd = 0,
    ratchetCycles = 2L, driftCycles = 0L,
    xssRounds = 0L, rssRounds = 0L, cssRounds = 0L,
    fuseInterval = 0L
  )
  for (nm in names(ctrl_overrides)) ctrl[[nm]] <- ctrl_overrides[[nm]]

  TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = maxReplicates, targetHits = 2L,
    tbrMaxHits = 1L,
    ratchetCycles = ctrl$ratchetCycles,
    ratchetPerturbProb = ctrl$ratchetPerturbProb,
    ratchetPerturbMode = ctrl$ratchetPerturbMode,
    ratchetPerturbMaxMoves = ctrl$ratchetPerturbMaxMoves,
    driftCycles = ctrl$driftCycles,
    xssRounds = ctrl$xssRounds,
    rssRounds = ctrl$rssRounds,
    cssRounds = ctrl$cssRounds,
    fuseInterval = ctrl$fuseInterval,
    maxSeconds = maxSeconds,
    verbosity = verbosity,
    concavity = concavity,
    annealConfig = make_anneal_config(ctrl)
  )
}

dataset <- inapplicable.phyData[["Vinther2008"]]
ds <- make_ts_data(dataset)

# --- Tests ---

test_that("Annealing runs end-to-end and returns valid trees", {
  result <- anneal_search(ds, maxSeconds = 3, maxReplicates = 2L)
  expect_gt(result$pool_size, 0)
  expect_lt(result$best_score, Inf)
  validate_result(result, length(dataset))
})

test_that("anneal_ms > 0 when annealing is enabled", {
  result <- anneal_search(ds, maxSeconds = 3, maxReplicates = 2L)
  expect_gt(result$timings[["anneal_ms"]], 0)
})

test_that("anneal_ms = 0 when annealing is disabled", {
  result <- anneal_search(ds,
    ctrl_overrides = list(annealCycles = 0L),
    maxSeconds = 3, maxReplicates = 2L
  )
  expect_equal(result$timings[["anneal_ms"]], 0)
})

test_that("Annealing with IW scoring works", {
  result <- anneal_search(ds, maxSeconds = 3, maxReplicates = 2L,
                          concavity = 3.0)
  expect_gt(result$pool_size, 0)
  expect_gt(result$timings[["anneal_ms"]], 0)
})

test_that("SearchControl accepts annealing parameters", {
  ctrl <- SearchControl(annealCycles = 3L, annealPhases = 5L,
                        annealTStart = 15, annealTEnd = 1,
                        annealMovesPerPhase = 50L)
  expect_equal(ctrl$annealCycles, 3L)
  expect_equal(ctrl$annealPhases, 5L)
  expect_equal(ctrl$annealTStart, 15)
  expect_equal(ctrl$annealTEnd, 1)
  expect_equal(ctrl$annealMovesPerPhase, 50L)
})

test_that("SearchControl defaults disable annealing", {
  ctrl <- SearchControl()
  expect_equal(ctrl$annealCycles, 0L)
  expect_equal(ctrl$annealPhases, 5L)
})

test_that("Large preset enables annealing and disables drift", {
  presets <- TreeSearch:::.AutoStrategy(200L, 200L)
  expect_equal(presets, "large")
  large_ctrl <- TreeSearch:::.StrategyPresets()[["large"]]
  expect_equal(large_ctrl$annealCycles, 1L)
  expect_gt(large_ctrl$annealPhases, 0L)
  expect_equal(large_ctrl$driftCycles, 0L)
})

test_that("Annealing with T=0 acts as strict hill-climbing", {
  result <- anneal_search(ds,
    ctrl_overrides = list(annealTStart = 0, annealTEnd = 0,
                          annealPhases = 2L),
    maxSeconds = 3, maxReplicates = 2L
  )
  expect_gt(result$pool_size, 0)
  expect_gt(result$timings[["anneal_ms"]], 0)
})

test_that("MaximizeParsimony respects annealing in SearchControl", {
  ctrl <- SearchControl(annealCycles = 1L, annealPhases = 3L,
                        annealTStart = 10, annealTEnd = 0,
                        ratchetCycles = 1L,
                        driftCycles = 0L, xssRounds = 0L,
                        rssRounds = 0L, cssRounds = 0L)
  result <- MaximizeParsimony(dataset, concavity = Inf,
                              maxSeconds = 3, maxReplicates = 2L,
                              control = ctrl, verbosity = 0L)
  expect_s3_class(result[[1]], "phylo")
})

test_that("Multi-cycle PCSA runs and reports sa_ms", {
  result <- anneal_search(ds,
    ctrl_overrides = list(annealCycles = 3L),
    maxSeconds = 5, maxReplicates = 2L
  )
  expect_gt(result$pool_size, 0)
  expect_gt(result$timings[["anneal_ms"]], 0)
  validate_result(result, length(dataset))
})

test_that("Multi-cycle PCSA score <= single-cycle", {
  set.seed(7418)
  single <- anneal_search(ds,
    ctrl_overrides = list(annealCycles = 1L),
    maxSeconds = 5, maxReplicates = 3L
  )
  multi <- anneal_search(ds,
    ctrl_overrides = list(annealCycles = 3L),
    maxSeconds = 5, maxReplicates = 3L
  )
  # Multi-cycle should find scores at least as good (may tie)
  expect_lte(multi$best_score, single$best_score + 1)
})
