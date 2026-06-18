# rasStarts > 1 triggers Goloboff-1999 RSS re-solve (RAS + TBR restarts) inside the
# sector search; rasStarts = 1 (default) is the prior single-TBR polish. These
# guard the R -> DrivenParams -> SectorParams plumbing (the kernel itself lives in
# ts_sector.cpp). Serial path (nThreads = 1L), where the sector params are wired.

data("inapplicable.phyData", package = "TreeSearch")
ds <- inapplicable.phyData[["Wortley2006"]]

test_that("rasStarts is exposed and defaulted by SearchControl", {
  expect_equal(SearchControl()$rasStarts, 1L)
  expect_equal(SearchControl(rasStarts = 3L)$rasStarts, 3L)
})

test_that("rasStarts > 1 makes the sector search do strictly more work", {
  set.seed(42)
  r1 <- MaximizeParsimony(ds, rasStarts = 1L, maxReplicates = 2L, targetHits = 99L,
                          nThreads = 1L, verbosity = 0L)
  set.seed(42)
  r3 <- MaximizeParsimony(ds, rasStarts = 3L, maxReplicates = 2L, targetHits = 99L,
                          nThreads = 1L, verbosity = 0L)
  expect_s3_class(r3, "multiPhylo")
  # Each extra RAS restart rebuilds + TBRs the sector from scratch, so re-solve
  # evaluates strictly more candidate rearrangements than the single-TBR polish.
  expect_gt(attr(r3, "candidates_evaluated"), attr(r1, "candidates_evaluated"))
})
