# Deep-search escalation: a CALLER-SET `targetHits` at 2x its default deepens
# per-replicate perturbation.  Caller-set is half the contract: the effort ladder
# raises `targetHits` itself from rung 5, and that must NOT engage the bundle.
library("TreeTools", quietly = TRUE)
data("inapplicable.phyData", package = "TreeSearch")
ds <- inapplicable.phyData[["Vinther2008"]]  # 23 tips

test_that(".TargetHitsEscalation reports the raise ratio, floored at 1", {
  TE <- TreeSearch:::.TargetHitsEscalation
  expect_equal(TE(10L, 10L), 1)      # at the default
  expect_equal(TE(20L, 10L), 2)      # doubled
  expect_equal(TE(96L, 96L), 1)
  expect_equal(TE(192L, 96L), 2)
  expect_equal(TE(48L, 96L), 1)      # LOWERED -> never below 1
  expect_equal(TE(4L, 10L), 1)       # the documented "one tree, quickly" idiom
  # Degenerate inputs fall back to "not raised" rather than erroring.
  expect_equal(TE(NA_integer_, 10L), 1)
  expect_equal(TE(10L, 0L), 1)
  expect_equal(TE(Inf, 10L), 1)
})

test_that(".ReachEscalationDeltas is the exact expected lever set", {
  expect_identical(
    TreeSearch:::.ReachEscalationDeltas(),
    list(
      ratchetPerturbMaxMoves = 0L,
      driftCycles = 25L,
      postRatchetSectorial = TRUE,
      stallEscalateFactor = 1.5,
      intraFuse = TRUE,
      poolSuboptimal = 3
    )
  )
})

test_that("every escalation delta is a real SearchControl field", {
  expect_true(all(names(TreeSearch:::.ReachEscalationDeltas()) %in%
                    names(SearchControl())))
})

test_that("escalation does NOT touch ratchetCycles", {
  # Ratchet depth belongs to .IwRatchetDepth, which scales it continuously
  # against a 36-matrix calibration (48 cycles, up to 115). A flat value here
  # would silently clobber that under implied weights, so the bundle must never
  # carry ratchetCycles.
  expect_false("ratchetCycles" %in% names(TreeSearch:::.ReachEscalationDeltas()))
  ctrl <- TreeSearch:::.ApplyReachEscalation(SearchControl(ratchetCycles = 48L),
                                            "thorough", escalation = 10,
                                            userSetHits = TRUE)
  expect_identical(ctrl[["ratchetCycles"]], 48L)
})

test_that("the two escalations compose without fighting over ratchetCycles", {
  # Mirrors the call-site order: .IwRatchetDepth sets the depth, then the bundle
  # is applied on top. The depth must survive, at its calibrated value.
  ctrl <- TreeSearch:::.ApplyStrategyPreset(
    SearchControl(), TreeSearch:::.StrategyPresets()[["thorough"]]
  )
  iw <- TreeSearch:::.IwRatchetDepth("thorough", concavity = 10,
                                     targetHits = 20L, defaultHits = 10L)
  expect_equal(iw, 96L)                    # 48 * escalation(2), under the 115 cap
  ctrl[["ratchetCycles"]] <- iw
  ctrl <- TreeSearch:::.ApplyReachEscalation(ctrl, "thorough", escalation = 2,
                                             userSetHits = TRUE)
  expect_identical(ctrl[["ratchetCycles"]], 96L)   # NOT clobbered by the bundle
  expect_identical(ctrl[["driftCycles"]], 25L)     # bundle still applied
})

test_that(".ApplyReachEscalation applies all deltas at or above the ratio", {
  deltas <- TreeSearch:::.ReachEscalationDeltas()
  for (strat in c("thorough", "large")) {
    for (esc in c(2, 2.5, 20)) {
      ctrl <- TreeSearch:::.ApplyReachEscalation(SearchControl(), strat,
                                                escalation = esc,
                                                userSetHits = TRUE)
      for (nm in names(deltas)) {
        expect_identical(ctrl[[nm]], deltas[[nm]], info = paste(strat, nm))
      }
    }
  }
})

test_that(".ApplyReachEscalation is inert below the ratio", {
  stock <- SearchControl()
  for (esc in c(1, 1.5, 1.99)) {
    expect_identical(TreeSearch:::.ApplyReachEscalation(stock, "thorough",
                                                        escalation = esc,
                                                        userSetHits = TRUE),
                     stock)
  }
  # Degenerate escalation must not escalate.
  for (esc in list(NA_real_, numeric(0), Inf)) {
    expect_identical(TreeSearch:::.ApplyReachEscalation(stock, "thorough",
                                                        escalation = esc,
                                                        userSetHits = TRUE),
                     stock)
  }
})

test_that(".ApplyReachEscalation requires the CALLER to have set targetHits", {
  # The ratio alone cannot distinguish "search harder" from a rung change:
  # .RungSpec()'s hitMultiplier doubles `targetHits` at rung 5 exactly when the
  # user did NOT set it, landing the ratio on 2.0 and tripping the `>=` gate.
  # Two independent lines measure that axis flat (see .ApplyReachEscalation), so
  # the bundle must stay shut unless the caller named the number themselves.
  stock <- SearchControl()
  for (notSet in list(FALSE, NA, NULL, logical(0))) {
    expect_identical(
      TreeSearch:::.ApplyReachEscalation(stock, "thorough", escalation = 10,
                                         userSetHits = notSet),
      stock, info = paste("userSetHits", format(notSet))
    )
  }
  # ... and open when they did, at the same ratio.
  expect_identical(
    TreeSearch:::.ApplyReachEscalation(stock, "thorough", escalation = 10,
                                       userSetHits = TRUE)[["driftCycles"]],
    25L
  )
})

test_that(".ApplyReachEscalation is scoped to thorough/large", {
  # `sprint` and `default` document themselves as fast/shallow ("3 ratchet
  # cycles, no drift"), and the A/B measured 0 better / 0 worse across their
  # 70 small- and medium-tier cells -- pure wall cost. They must not escalate.
  stock <- SearchControl()
  for (strat in c("sprint", "default", "none", NA_character_, character(0))) {
    expect_identical(
      TreeSearch:::.ApplyReachEscalation(stock, strat, escalation = 10,
                                         userSetHits = TRUE),
      stock, info = paste("strategy", strat)
    )
  }
})

test_that(".ApplyReachEscalation preserves caller-set fields", {
  ctrl <- TreeSearch:::.ApplyReachEscalation(
    SearchControl(), "thorough", escalation = 4, userSetHits = TRUE,
    userSet = c("driftCycles", "intraFuse")
  )
  expect_identical(ctrl[["driftCycles"]], SearchControl()[["driftCycles"]])
  expect_identical(ctrl[["intraFuse"]], SearchControl()[["intraFuse"]])
  # a field the caller did not set still escalates
  expect_identical(ctrl[["postRatchetSectorial"]], TRUE)
})

test_that("the default targetHits never escalates", {
  # Off-by-default guarantee: the size-scaled default is its own reference, so an
  # un-tuned search sits at ratio 1 -- below the trigger -- for any tip count.
  TE <- TreeSearch:::.TargetHitsEscalation
  for (n in c(5L, 23L, 50L, 88L, 200L, 482L, 4062L)) {
    defHits <- max(10L, as.integer(n / 5))
    expect_lt(TE(defHits, defHits), TreeSearch:::.reachEscalationMinRatio)
  }
})

test_that("escalation measurably deepens the search end to end", {
  # Not a smoke test: asserts a signal that DISAPPEARS if the feature is removed.
  # The bundle multiplies per-replicate work (drift 2 -> 25, deep kick, an extra
  # sectorial pass), so at matched replicates the escalated run must evaluate
  # substantially more candidates. Vinther2008 (23 tips): default targetHits = 10,
  # so 20 is exactly 2x.
  # Vinther2008 is 23 tips, so .AutoRung() gives rung 1 (`sprint`); `effort = 2`
  # is rung 3 (`thorough`), which is in scope and below the rung-5 hitMultiplier.
  runCand <- function(hits) {
    set.seed(4242)
    r <- MaximizeParsimony(ds, effort = 2L, maxReplicates = 2L,
                           targetHits = hits, maxSeconds = 0, verbosity = 0L)
    list(cand = as.double(attr(r, "candidates_evaluated")), res = r)
  }
  ordinary <- runCand(10L)
  escalated <- runCand(20L)
  expect_s3_class(escalated$res, "multiPhylo")
  expect_true(is.finite(attr(escalated$res, "score")))
  expect_true(is.finite(ordinary$cand) && ordinary$cand > 0)
  # Same seed and same replicate count, so with the feature removed the two runs
  # would be identical and the ratio exactly 1 (as the `sprint` test below shows).
  # Measured here at ~1.49; 1.2 leaves room for stochastic drift while still
  # failing outright if the escalation stops engaging.  NB the ratio is much
  # larger against a shallower preset -- `thorough` already drifts, so the
  # increment over it is smaller than over `sprint`.
  expect_gt(escalated$cand, 1.2 * ordinary$cand)
})

test_that("sprint is NOT escalated end to end", {
  # The scope guard, observably. Under a preset outside the escalation's scope,
  # raising targetHits cannot change per-replicate work at all: with the replicate
  # count fixed the two runs are identical, so the candidate counts must match
  # EXACTLY. This fails the moment the strategy gate is loosened.
  runCand <- function(hits) {
    set.seed(99L)
    r <- MaximizeParsimony(ds, effort = 0L, maxReplicates = 2L,
                           targetHits = hits, maxSeconds = 0, verbosity = 0L)
    as.double(attr(r, "candidates_evaluated"))
  }
  expect_identical(runCand(20L), runCand(10L))
})

test_that("an escalated search still returns only best-score trees", {
  # The bundle raises `poolSuboptimal` internally (intraFuse needs recipients).
  # That must not leak into the result: with collapse = FALSE the pool is returned
  # verbatim, so without the guard the caller would silently get trees up to 3
  # steps worse than attr(, "score") from a result documented as the best found.
  set.seed(31L)
  r <- MaximizeParsimony(ds, effort = 2L, maxReplicates = 3L,
                         targetHits = 20L, collapse = FALSE, verbosity = 0L)
  best <- attr(r, "score")
  expect_true(is.finite(best))
  scores <- vapply(r, function(t) TreeLength(t, ds, concavity = Inf), double(1))
  expect_true(all(scores == best))
})

test_that("the effort ladder's own targetHits rise does NOT deepen the search", {
  # The provenance test, with every VALUE held equal.  At rung 5 the ladder
  # doubles the 23-tip default of 10 to 20 by itself; the second run names 20,
  # so the ladder skips its multiplier and leaves it at 20.  Both runs therefore
  # search with targetHits = 20, the same preset (`large`), the same replicate
  # cap and the same seed -- the ONLY difference is who set the number.  Only the
  # caller's version may deepen the perturbation.  With the gate reading the
  # ratio alone (as it first did) both runs escalate and the counts match, which
  # would fire this bundle on every dataset over 120 tips at `effort = 1`.
  runCand <- function(...) {
    set.seed(808L)
    r <- MaximizeParsimony(ds, effort = 4L, maxReplicates = 2L,
                           maxSeconds = 0, verbosity = 0L, ...)
    as.double(attr(r, "candidates_evaluated"))
  }
  ladder <- runCand()
  asked <- runCand(targetHits = 20L)
  expect_true(is.finite(ladder) && ladder > 0)
  expect_gt(asked, 1.2 * ladder)
  # And the ladder run is indistinguishable from the rung below it, whose ratio
  # is 1: rung 5 with a user-set `maxReplicates` differs only in the hit target.
  set.seed(808L)
  rung4 <- as.double(attr(
    MaximizeParsimony(ds, effort = 3L, maxReplicates = 2L, maxSeconds = 0,
                      verbosity = 0L),
    "candidates_evaluated"
  ))
  expect_identical(ladder, rung4)
})

test_that("a caller's own poolSuboptimal is still honoured", {
  # The guard above must not steal the documented behaviour from someone who
  # asked for suboptimal trees themselves.
  set.seed(31L)
  r <- MaximizeParsimony(ds, effort = 2L, maxReplicates = 3L,
                         targetHits = 20L, poolSuboptimal = 3,
                         collapse = FALSE, verbosity = 0L)
  expect_s3_class(r, "multiPhylo")
  scores <- vapply(r, function(t) TreeLength(t, ds, concavity = Inf), double(1))
  expect_true(all(scores <= attr(r, "score") + 3))
})
