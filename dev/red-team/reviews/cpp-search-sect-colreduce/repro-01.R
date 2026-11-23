# repro-01.R  --  F1: subtree.total_words / n_blocks desync after column reduction
#
# NOT RUN by the reviewer (Windows build; needs the flagged build). This script
# WOULD demonstrate the crash / wrong score if executed with TS_SECT_COLREDUCE=1
# against the uncommitted ts_sector.cpp.
#
# Mechanism: build_reduced_dataset() sets rd.subtree.total_words = ds.total_words
# and rd.subtree.n_blocks = ds.n_blocks (ts_sector.cpp:519-520). reduce_sector_
# columns_ew() then shrinks rd.data.total_words/n_blocks (ts_sector.cpp:401-402)
# but NEVER re-syncs rd.subtree.*. The subtree state arrays are allocated with
# the NEW rd.data.total_words (651) yet load_tip_states + every fitch_* scorer
# index by the STALE (large) rd.subtree.total_words -> heap overflow on both the
# source (rd.data.tip_states, now sized n_tip*ntw) and the destination (prelim,
# sized n_node*ntw).  Fires on ANY sector that has >=1 constant-within-sector
# character (i.e. essentially every sector).
#
# Need: a dataset big enough to trigger sectorial search (n_tips >= 2*sector_min
# ~ 2*30 = 60 by default), EW, weight 1, no NA, with many invariant columns so
# the reduction definitely fires.

library(TreeSearch)

set.seed(1)
n_tip  <- 80
n_char <- 200
# Make ~80% of columns constant (single state shared by all tips) so reduction
# fires hard in every sector; the rest random binary.
m <- matrix("0", nrow = n_tip, ncol = n_char)
informative_cols <- sample(seq_len(n_char), size = round(0.2 * n_char))
for (j in informative_cols) {
  m[, j] <- sample(c("0", "1"), n_tip, replace = TRUE)
}
rownames(m) <- paste0("t", seq_len(n_tip))
dat <- TreeTools::MatrixToPhyDat(m)

start <- TreeTools::RandomTree(n_tip, root = TRUE)

run <- function(flag) {
  Sys.setenv(TS_SECT_COLREDUCE = flag)
  set.seed(42)
  tr <- MaximizeParsimony(dat, tree = start, ratchIter = 2,
                          tbrIter = 2, maxHits = 5, verbosity = 0)
  TreeLength(tr[[1]], dat)
}

len_off <- run("0")
len_on  <- run("1")   # EXPECT: crash (ASan) or len_on != len_off

cat(sprintf("OFF=%s  ON=%s  dScore=%s\n", len_off, len_on, len_on - len_off))
stopifnot(identical(len_off, len_on))  # FAILS under the bug
