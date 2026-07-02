# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()
## Tests for impose_constraint() — post-hoc topology repair (T-213)

library("TreeTools")

# Reuse helpers from test-ts-constraint-small.R
make_ds5 <- function() {
  phangorn::phyDat(
    matrix(c("0", "0", "0", "1", "1",
             "0", "1", "0", "1", "0"),
           nrow = 5, dimnames = list(paste0("t", 1:5), NULL)),
    type = "USER", levels = c("0", "1")
  )
}

check_constraint <- function(tree, constraint) {
  tips <- sort(constraint$tip.label)
  tree_sp <- as.Splits(tree, tipLabels = tips)
  cons_sp <- as.Splits(constraint, tipLabels = tips)
  tm <- as.logical(tree_sp)
  cm <- as.logical(cons_sp)
  if (!is.matrix(tm)) tm <- matrix(tm, nrow = 1)
  if (!is.matrix(cm)) cm <- matrix(cm, nrow = 1)
  all(apply(cm, 1, function(c_row) {
    any(apply(tm, 1, function(t_row) {
      all(c_row == t_row) || all(c_row == !t_row)
    }))
  }))
}

# Larger dataset for meaningful NNI perturbation
make_ds12 <- function() {
  set.seed(8113)
  m <- matrix(sample(c("0", "1"), 12 * 6, replace = TRUE),
              nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  phangorn::phyDat(m, type = "USER", levels = c("0", "1"))
}

# ----- NNI perturbation + constraint repair -----

test_that("T-213: NNI perturbation works under constraints (5 tips)", {
  ds5 <- make_ds5()
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  set.seed(4529)
  result <- MaximizeParsimony(
    ds5, constraint = cons,
    maxReplicates = 2L, verbosity = 0L,
    control = SearchControl(nniPerturbCycles = 3L,
                            nniPerturbFraction = 0.5)
  )
  expect_s3_class(result, "multiPhylo")
  for (i in seq_along(result)) {
    expect_true(check_constraint(result[[i]], cons),
                info = paste("tree", i))
  }
})

test_that("T-213: NNI perturbation under single constraint split", {
  ds5 <- make_ds5()
  cons1 <- ape::read.tree(text = "((t1,t2),t3,t4,t5);")

  set.seed(6719)
  result <- MaximizeParsimony(
    ds5, constraint = cons1,
    maxReplicates = 2L, verbosity = 0L,
    control = SearchControl(nniPerturbCycles = 5L,
                            nniPerturbFraction = 0.8)
  )
  for (i in seq_along(result)) {
    expect_true(check_constraint(result[[i]], cons1))
  }
})

test_that("T-213: NNI perturbation under constraints (12 tips)", {
  ds12 <- make_ds12()
  cons12 <- ape::read.tree(
    text = "((t1,t2,t3,t4),(t5,t6,(t7,(t8,t9,t10,t11,t12))));"
  )

  set.seed(3146)
  result <- MaximizeParsimony(
    ds12, constraint = cons12,
    maxReplicates = 3L, verbosity = 0L,
    control = SearchControl(nniPerturbCycles = 5L,
                            nniPerturbFraction = 0.5)
  )
  expect_s3_class(result, "multiPhylo")
  for (i in seq_along(result)) {
    expect_true(check_constraint(result[[i]], cons12),
                info = paste("tree", i))
  }
})

test_that("T-213: nested constraints with NNI perturbation", {
  ds12 <- make_ds12()
  # Two nested constraint splits: {t1,t2,t3,t4} and {t1,t2}
  cons_nested <- ape::read.tree(
    text = "(((t1,t2),t3,t4),(t5,t6,t7,t8,t9,t10,t11,t12));"
  )

  set.seed(5612)
  result <- MaximizeParsimony(
    ds12, constraint = cons_nested,
    maxReplicates = 3L, verbosity = 0L,
    control = SearchControl(nniPerturbCycles = 4L,
                            nniPerturbFraction = 0.6)
  )
  for (i in seq_along(result)) {
    expect_true(check_constraint(result[[i]], cons_nested),
                info = paste("tree", i))
  }
})

# ----- Fuse + constraint repair -----

test_that("T-213: fuse under constraints preserves constraint", {
  ds12 <- make_ds12()
  cons12 <- ape::read.tree(
    text = "((t1,t2,t3,t4),(t5,t6,(t7,(t8,t9,t10,t11,t12))));"
  )

  set.seed(2754)
  result <- MaximizeParsimony(
    ds12, constraint = cons12,
    maxReplicates = 4L, verbosity = 0L,
    control = SearchControl(fuseInterval = 2L)
  )
  for (i in seq_along(result)) {
    expect_true(check_constraint(result[[i]], cons12),
                info = paste("tree", i))
  }
})

# ----- IW scoring + NNI perturbation + constraints -----

test_that("T-213: IW scoring + NNI perturbation + constraints", {
  ds5 <- make_ds5()
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  set.seed(9833)
  result <- MaximizeParsimony(
    ds5, constraint = cons, concavity = 10,
    maxReplicates = 2L, verbosity = 0L,
    control = SearchControl(nniPerturbCycles = 3L)
  )
  for (i in seq_along(result)) {
    expect_true(check_constraint(result[[i]], cons),
                info = paste("tree", i))
  }
})

# ----- Root-child move in impose_constraint -----
# When fuse produces a tree where constraint tips span the root, fixing the
# violation requires moving a root child. This tests the topology_spr helper
# that handles the root-child case (previously skipped by spr_clip guard).

test_that("impose_constraint repairs root-child violations (8 tips)", {
  ds8 <- phangorn::phyDat(
    matrix(c("0","0","0","0","1","1","1","1",
             "0","1","0","1","0","1","0","1",
             "0","0","1","1","0","0","1","1"),
           nrow = 8, dimnames = list(paste0("t", 1:8), NULL)),
    type = "USER", levels = c("0", "1")
  )
  # Constraint requires {t1,t2,t3,t4} on one side — violations likely
  # when fuse produces trees splitting this group across the root.
  cons8 <- ape::read.tree(text = "((t1,t2,t3,t4),(t5,t6,t7,t8));")

  n_ok <- 0L
  n_total <- 0L
  for (s in c(1147L, 2258L, 3369L, 4470L, 5581L)) {
    set.seed(s)
    result <- MaximizeParsimony(
      ds8, constraint = cons8,
      maxReplicates = 6L, verbosity = 0L,
      control = SearchControl(adaptiveStart = TRUE)
    )
    for (i in seq_along(result)) {
      n_total <- n_total + 1L
      if (check_constraint(result[[i]], cons8)) n_ok <- n_ok + 1L
    }
  }
  expect_equal(n_ok, n_total,
               info = paste(n_ok, "/", n_total, "satisfy constraint"))
})

test_that("impose_constraint repairs root-child violations (12 tips, nested)", {
  set.seed(6293)
  ds12 <- phangorn::phyDat(
    matrix(sample(0:1, 12 * 6, replace = TRUE),
           nrow = 12, dimnames = list(paste0("t", 1:12), NULL)),
    type = "USER", levels = c("0", "1")
  )
  # Nested constraint: {t1..t6} and within it {t1,t2,t3}
  cons12 <- ape::read.tree(text = "((t1,t2,t3),(t4,t5,t6),(t7,t8,t9,t10,t11,t12));")

  n_ok <- 0L
  n_total <- 0L
  for (s in c(7104L, 8215L, 9326L)) {
    set.seed(s)
    result <- MaximizeParsimony(
      ds12, constraint = cons12,
      maxReplicates = 4L, verbosity = 0L, nThreads = 2L,
      control = SearchControl(adaptiveStart = TRUE)
    )
    for (i in seq_along(result)) {
      n_total <- n_total + 1L
      if (check_constraint(result[[i]], cons12)) n_ok <- n_ok + 1L
    }
  }
  expect_equal(n_ok, n_total,
               info = paste(n_ok, "/", n_total, "satisfy constraint"))
})

# ----- T-327: stale best_node in impose_one_pass -> cyclic tree -----
# impose_one_pass() captured `best_node` once and reused it across a batch of
# topology_spr() moves. When a move relocated that node id (its parent WAS
# best_node), the stale reference let a later graft target land adjacent to or
# inside the moved subtree, producing a cyclic / double-parented tree; the next
# DFS (collect_edges_*, build_postorder) then walked the cycle and allocated
# without bound -> std::bad_alloc. Trigger: a topological constraint plus
# nniPerturbCycles > 0 (opt-in) on a matrix with both misplaced-in and
# misplaced-out tips relative to best_node. These configs crashed the process
# before the snapshot-validate-revert guard; they must now complete and honour
# the constraint.

test_that("T-327: impose_one_pass survives stale best_node (no bad_alloc)", {
  cons_a <- ape::read.tree(text = "((t1,t2),(t3,t4),(t5,t6,t7,t8,t9,t10,t11,t12));")
  cons_b <- ape::read.tree(text = "((t1,t2,t3,t4,t5,t6),(t7,t8),(t9,t10,t11,t12));")

  # nThreads = 1L is deliberate: serial mode makes set.seed() drive the C++
  # search RNG deterministically (make_rng seeds from R's RNG only when serial),
  # so these seed->crash mappings reproduce across machines/CI. Every (seed,
  # cons) pair below crashed the process with std::bad_alloc on the pre-fix
  # build (confirmed by scratchpad/repro_sweep_serial.R).
  n_ok <- 0L
  n_total <- 0L
  for (seed in c(22L, 23L, 24L)) {
    set.seed(seed)
    m <- matrix(sample(c("0", "1"), 12 * 20, replace = TRUE),
                nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
    pd <- phangorn::phyDat(m, type = "USER", levels = c("0", "1"))
    for (cons in list(cons_a, cons_b)) {
      set.seed(seed)
      expect_no_error(
        result <- MaximizeParsimony(
          pd, constraint = cons,
          maxReplicates = 5L, verbosity = 0L, nThreads = 1L,
          control = SearchControl(nniPerturbCycles = 8L,
                                  nniPerturbFraction = 0.7,
                                  ratchetCycles = 0L)
        )
      )
      for (i in seq_along(result)) {
        n_total <- n_total + 1L
        if (check_constraint(result[[i]], cons)) n_ok <- n_ok + 1L
      }
    }
  }
  # Every returned tree must honour its constraint (the verify-before-capture
  # guard, T-213, discards any tree impose_one_pass could not fully repair).
  expect_equal(n_ok, n_total,
               info = paste(n_ok, "/", n_total, "satisfy constraint"))
})
