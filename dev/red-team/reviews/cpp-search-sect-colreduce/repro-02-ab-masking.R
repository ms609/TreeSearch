# repro-02-ab-masking.R -- WHY the off-vs-on dScore=0 A/B can FALSE-PASS.
#
# The sector accept gate (rss_search/xss_search, ts_sector.cpp:1328-1375) is:
#   improved   = sector_best < sector_current        # both from the BROKEN reduced scorer
#   if (accept && sector_best <= sector_current) {
#       reinsert_sector(tree, rd);
#       new_score = score_tree(tree, ds);             # EXACT full-ds rescore (AUTHORITATIVE)
#       if (new_score < result.best_score) keep
#       else if (== && accept_equal)      keep
#       else                              restore_clade(...)  # REVERT
#   }
# => the exact full rescore + revert already tolerates a WRONG reduced score
#    (the HTU is itself an approximation by design). So when TS_SECT_COLREDUCE=1
#    corrupts the reduced search (stride bug F1), the worst observable effect on
#    the FINAL score is: garbage candidates get reverted, search does less useful
#    work, and the FINAL MP score can still equal the OFF run => dScore = 0.
#
# CONSEQUENCE: a 9-run (3 datasets x 3 seeds) "final dScore must be 0" A/B is
# INSUFFICIENT. It can report PASS while:
#   (a) the feature is doing OOB reads/writes on every sector pick (UB; may also
#       silently corrupt the heap and crash only sometimes / only on some libc), and
#   (b) the "byte-identical sector trajectory" claim is FALSE (the reduced search
#       is corrupted, not faster-but-equal).
#
# The A/B is only trustworthy if it ALSO:
#   1. runs under ASan (branch *asan*) -- the source-side OOB read in
#      load_tip_states (n_sector_tips*OLD_tw words memcpy'd from a
#      n_sector_tips*NEW_tw buffer) is GUARANTEED to fault whenever the
#      reduction fires; and
#   2. PROVES the reduction actually fired on the test data (else it validates
#      nothing). Use the TS_AUDIT_PROBE counters (g_sect_inf_chars/tot_chars,
#      fp_blocks < tot_blocks) and assert >=1 column was dropped on >=1 sector.
#      With small/clean matrices or sectors where surv.size()>=tot_active, NO
#      column is dropped, NEW==OLD, the bug never bites, and dScore=0 is vacuous.
#
# This script just makes (2) observable at R level: build with -DTS_AUDIT_PROBE,
# run a sectorial search, and confirm the "SECT_REDUCE ... fp/tot_blocks" line
# shows fp < tot (i.e. drops happen) on the chosen dataset before trusting any
# colreduce A/B.
library(TreeSearch)
data("Lobo", package = "TreeSearch", envir = environment())
dat <- Lobo.phy
# (build the package with PKG_CPPFLAGS += -DTS_AUDIT_PROBE first; the
#  SECT_REDUCE diagnostic prints to stderr every 2000 sector picks)
Sys.setenv(TS_SECT_COLREDUCE = "1")
invisible(MaximizeParsimony(dat, ratchIter = 0, tbrIter = 2,
                            maxHits = 20, maxTime = 5/60, verbosity = 0))
# Inspect stderr: a fp/tot_blocks < 1.0 confirms columns are dropped => the
# stride mismatch IS exercised => the ASan A/B is meaningful.
