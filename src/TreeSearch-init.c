#define USE_RINTERNALS

#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <stdlib.h> /* for NULL */
#include <R_ext/Rdynload.h>

#include "build_postorder.h"

extern SEXP _TreeSearch_nni(SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_spr(SEXP, SEXP);
extern SEXP _TreeSearch_spr_moves(SEXP);
extern SEXP _TreeSearch_tbr(SEXP, SEXP);
// extern SEXP _TreeSearch_tbr_moves(SEXP);
extern SEXP _TreeSearch_all_spr(SEXP, SEXP);
extern SEXP _TreeSearch_all_tbr(SEXP, SEXP);

extern SEXP _TreeSearch_expected_mi(SEXP, SEXP);
extern SEXP _TreeSearch_mi_key(SEXP, SEXP);

// extern SEXP _TreeSearch_astar_search_r(SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_quartet_concordance(SEXP, SEXP);
extern SEXP _TreeSearch_ts_fitch_score(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_na_char_steps(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_char_steps(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_na_debug_char(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_nni_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_debug_clip(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_spr_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_tbr_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_drift_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_test_indirect(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_ratchet_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_compute_splits(SEXP, SEXP);
extern SEXP _TreeSearch_ts_trees_equal(SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_pool_test(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_wagner_tree(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_random_wagner_tree(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_tree_fuse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_sector_diag(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_rss_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_xss_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_driven_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_resample_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_parallel_resample(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_successive_approx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_bench_tbr_phases(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_simplify_diag(SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_MaddisonSlatkin(SEXP, SEXP);
extern SEXP _TreeSearch_MaddisonSlatkin_clear_cache();
extern SEXP _TreeSearch_ts_hsj_score(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_sankoff_test(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_mc_fitch_scores(SEXP, SEXP);
extern SEXP _TreeSearch_ts_wagner_bias_bench(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
/* ts_stochastic_tbr and ts_parallel_temper removed — on feature/parallel-temper */
extern SEXP _TreeSearch_ts_test_strategy_tracker(SEXP, SEXP);
extern SEXP _TreeSearch_ts_tbr_diagnostics(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_ev_cache_key_probe(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_ls_fit(SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_ls_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_collapsed_flags_debug(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_collapse_pool(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef callMethods[] = {
  {"_TreeSearch_nni",               (DL_FUNC) &_TreeSearch_nni, 3},
  {"_TreeSearch_spr",               (DL_FUNC) &_TreeSearch_spr, 2},
  {"_TreeSearch_all_spr",           (DL_FUNC) &_TreeSearch_all_spr, 2},
  {"_TreeSearch_spr_moves",         (DL_FUNC) &_TreeSearch_spr_moves, 1},
  {"_TreeSearch_tbr",               (DL_FUNC) &_TreeSearch_tbr, 2},
  {"_TreeSearch_all_tbr",           (DL_FUNC) &_TreeSearch_all_tbr, 2},
//  {"_TreeSearch_tbr_moves",         (DL_FUNC) &_TreeSearch_tbr_moves, 1},
  {"_TreeSearch_expected_mi",       (DL_FUNC) &_TreeSearch_expected_mi, 2},
  {"_TreeSearch_mi_key",            (DL_FUNC) &_TreeSearch_mi_key, 2},

  // {"_TreeSearch_astar_search_r",    (DL_FUNC) &_TreeSearch_astar_search_r, 3},
  {"_TreeSearch_quartet_concordance",(DL_FUNC) &_TreeSearch_quartet_concordance, 2},
  {"_TreeSearch_ts_fitch_score",    (DL_FUNC) &_TreeSearch_ts_fitch_score, 12},
  {"_TreeSearch_ts_na_char_steps", (DL_FUNC) &_TreeSearch_ts_na_char_steps, 5},
  {"_TreeSearch_ts_char_steps", (DL_FUNC) &_TreeSearch_ts_char_steps, 5},
  {"_TreeSearch_ts_na_debug_char", (DL_FUNC) &_TreeSearch_ts_na_debug_char, 6},
  {"_TreeSearch_ts_nni_search",    (DL_FUNC) &_TreeSearch_ts_nni_search, 8},
  {"_TreeSearch_ts_debug_clip",    (DL_FUNC) &_TreeSearch_ts_debug_clip, 6},
  {"_TreeSearch_ts_spr_search",    (DL_FUNC) &_TreeSearch_ts_spr_search, 8},
  {"_TreeSearch_ts_tbr_search",    (DL_FUNC) &_TreeSearch_ts_tbr_search, 10},
  {"_TreeSearch_ts_drift_search",  (DL_FUNC) &_TreeSearch_ts_drift_search, 11},
  {"_TreeSearch_ts_test_indirect", (DL_FUNC) &_TreeSearch_ts_test_indirect, 8},
  {"_TreeSearch_ts_ratchet_search", (DL_FUNC) &_TreeSearch_ts_ratchet_search, 14},
  {"_TreeSearch_ts_compute_splits", (DL_FUNC) &_TreeSearch_ts_compute_splits, 2},
  {"_TreeSearch_ts_trees_equal",   (DL_FUNC) &_TreeSearch_ts_trees_equal, 3},
  {"_TreeSearch_ts_pool_test",     (DL_FUNC) &_TreeSearch_ts_pool_test, 5},
  {"_TreeSearch_ts_wagner_tree",   (DL_FUNC) &_TreeSearch_ts_wagner_tree, 14},
  {"_TreeSearch_ts_random_wagner_tree", (DL_FUNC) &_TreeSearch_ts_random_wagner_tree, 13},
  {"_TreeSearch_ts_sector_diag", (DL_FUNC) &_TreeSearch_ts_sector_diag, 6},
  {"_TreeSearch_ts_rss_search",  (DL_FUNC) &_TreeSearch_ts_rss_search, 13},
  {"_TreeSearch_ts_xss_search",  (DL_FUNC) &_TreeSearch_ts_xss_search, 12},
  {"_TreeSearch_ts_tree_fuse",   (DL_FUNC) &_TreeSearch_ts_tree_fuse, 11},
  {"_TreeSearch_ts_driven_search", (DL_FUNC) &_TreeSearch_ts_driven_search, 10},
  {"_TreeSearch_ts_resample_search", (DL_FUNC) &_TreeSearch_ts_resample_search, 25},
{"_TreeSearch_ts_parallel_resample", (DL_FUNC) &_TreeSearch_ts_parallel_resample, 27},
  {"_TreeSearch_ts_successive_approx", (DL_FUNC) &_TreeSearch_ts_successive_approx, 25},
  {"_TreeSearch_ts_bench_tbr_phases", (DL_FUNC) &_TreeSearch_ts_bench_tbr_phases, 7},
  {"_TreeSearch_ts_simplify_diag", (DL_FUNC) &_TreeSearch_ts_simplify_diag, 4},
  {"_TreeSearch_MaddisonSlatkin", (DL_FUNC) &_TreeSearch_MaddisonSlatkin, 2},
  {"_TreeSearch_MaddisonSlatkin_clear_cache", (DL_FUNC) &_TreeSearch_MaddisonSlatkin_clear_cache, 0},
  {"_TreeSearch_ts_hsj_score", (DL_FUNC) &_TreeSearch_ts_hsj_score, 9},
  {"_TreeSearch_ts_sankoff_test", (DL_FUNC) &_TreeSearch_ts_sankoff_test, 5},
  {"_TreeSearch_ts_wagner_bias_bench", (DL_FUNC) &_TreeSearch_ts_wagner_bias_bench, 10},
  /* ts_stochastic_tbr (9) and ts_parallel_temper (10) removed */

  {"_TreeSearch_mc_fitch_scores",    (DL_FUNC) &_TreeSearch_mc_fitch_scores, 2},
  {"RANDOM_TREE",                   (DL_FUNC) &RANDOM_TREE, 1},
  {"_TreeSearch_ts_test_strategy_tracker", (DL_FUNC) &_TreeSearch_ts_test_strategy_tracker, 2},
  {"_TreeSearch_ts_tbr_diagnostics", (DL_FUNC) &_TreeSearch_ts_tbr_diagnostics, 12},
  {"_TreeSearch_ts_ev_cache_key_probe", (DL_FUNC) &_TreeSearch_ts_ev_cache_key_probe, 9},
  {"_TreeSearch_ts_ls_fit", (DL_FUNC) &_TreeSearch_ts_ls_fit, 4},
  {"_TreeSearch_ts_ls_search", (DL_FUNC) &_TreeSearch_ts_ls_search, 6},
  {"_TreeSearch_ts_collapsed_flags_debug", (DL_FUNC) &_TreeSearch_ts_collapsed_flags_debug, 6},
  {"_TreeSearch_ts_collapse_pool", (DL_FUNC) &_TreeSearch_ts_collapse_pool, 9},
  {NULL, NULL, 0}
};

void R_init_TreeSearch(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
