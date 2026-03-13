#define USE_RINTERNALS

#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <stdlib.h> /* for NULL */
#include <R_ext/Rdynload.h>

#include "mpl.h"
#include "RMorphyUtils.h"
#include "RMorphy.h"
#include "build_postorder.h"

extern SEXP _TreeSearch_nni(SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_spr(SEXP, SEXP);
extern SEXP _TreeSearch_spr_moves(SEXP);
extern SEXP _TreeSearch_tbr(SEXP, SEXP);
// extern SEXP _TreeSearch_tbr_moves(SEXP);
extern SEXP _TreeSearch_all_spr(SEXP, SEXP);
extern SEXP _TreeSearch_all_tbr(SEXP, SEXP);
extern SEXP _TreeSearch_preorder_morphy(SEXP, SEXP);
extern SEXP _TreeSearch_preorder_morphy_by_char(SEXP, SEXP);
extern SEXP _TreeSearch_morphy_iw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_morphy_profile(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP _TreeSearch_expected_mi(SEXP, SEXP);
extern SEXP _TreeSearch_mi_key(SEXP, SEXP);

// extern SEXP _TreeSearch_astar_search_r(SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_quartet_concordance(SEXP, SEXP);
extern SEXP _TreeSearch_ts_fitch_score(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_na_char_steps(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_na_debug_char(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_nni_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_debug_clip(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_spr_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_tbr_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_drift_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_test_indirect(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_ratchet_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_compute_splits(SEXP, SEXP);
extern SEXP _TreeSearch_ts_trees_equal(SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_pool_test(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_wagner_tree(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_random_wagner_tree(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_tree_fuse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_sector_diag(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_rss_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_xss_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeSearch_ts_driven_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef callMethods[] = {
  {"_R_wrap_mpl_new_Morphy",        (DL_FUNC) &_R_wrap_mpl_new_Morphy, 0},
  {"_R_wrap_mpl_delete_Morphy",     (DL_FUNC) &_R_wrap_mpl_delete_Morphy, 1},
  {"_R_wrap_mpl_init_Morphy",       (DL_FUNC) &_R_wrap_mpl_init_Morphy, 3},
  {"_R_wrap_mpl_get_numtaxa",       (DL_FUNC) &_R_wrap_mpl_get_numtaxa, 1},
  {"_R_wrap_mpl_get_num_charac",    (DL_FUNC) &_R_wrap_mpl_get_num_charac, 1},
  {"_R_wrap_mpl_attach_symbols",    (DL_FUNC) &_R_wrap_mpl_attach_symbols, 2},
  {"_R_wrap_mpl_get_symbols",       (DL_FUNC) &_R_wrap_mpl_get_symbols, 1},
  {"_R_wrap_mpl_attach_rawdata",    (DL_FUNC) &_R_wrap_mpl_attach_rawdata, 2},
  {"_R_wrap_mpl_set_parsim_t",      (DL_FUNC) &_R_wrap_mpl_set_parsim_t, 3},
  {"_R_wrap_mpl_get_gaphandl",      (DL_FUNC) &_R_wrap_mpl_get_gaphandl, 1},
  {"_R_wrap_mpl_set_gaphandl",      (DL_FUNC) &_R_wrap_mpl_set_gaphandl, 2},
  {"_R_wrap_mpl_set_num_internal_nodes", (DL_FUNC) &_R_wrap_mpl_set_num_internal_nodes, 2},
  {"_R_wrap_mpl_get_num_internal_nodes", (DL_FUNC) &_R_wrap_mpl_get_num_internal_nodes, 1},
  {"_R_wrap_mpl_apply_tipdata",     (DL_FUNC) &_R_wrap_mpl_apply_tipdata, 1},
  {"_R_wrap_mpl_set_charac_weight", (DL_FUNC) &_R_wrap_mpl_set_charac_weight, 3},
  {"_R_wrap_mpl_get_charac_weight", (DL_FUNC) &_R_wrap_mpl_get_charac_weight, 2},
  {"_R_wrap_mpl_first_down_recon",  (DL_FUNC) &_R_wrap_mpl_first_down_recon, 4},
  {"_R_wrap_mpl_first_up_recon",    (DL_FUNC) &_R_wrap_mpl_first_up_recon, 5},
  {"_R_wrap_mpl_second_down_recon", (DL_FUNC) &_R_wrap_mpl_second_down_recon, 4},
  {"_R_wrap_mpl_second_up_recon",   (DL_FUNC) &_R_wrap_mpl_second_up_recon, 5},
  {"_R_wrap_mpl_update_tip",        (DL_FUNC) &_R_wrap_mpl_update_tip, 3},
  {"_R_wrap_mpl_update_lower_root", (DL_FUNC) &_R_wrap_mpl_update_lower_root, 3},
  {"_TreeSearch_nni",               (DL_FUNC) &_TreeSearch_nni, 3},
  {"_TreeSearch_spr",               (DL_FUNC) &_TreeSearch_spr, 2},
  {"_TreeSearch_all_spr",           (DL_FUNC) &_TreeSearch_all_spr, 2},
  {"_TreeSearch_spr_moves",         (DL_FUNC) &_TreeSearch_spr_moves, 1},
  {"_TreeSearch_tbr",               (DL_FUNC) &_TreeSearch_tbr, 2},
  {"_TreeSearch_all_tbr",           (DL_FUNC) &_TreeSearch_all_tbr, 2},
//  {"_TreeSearch_tbr_moves",         (DL_FUNC) &_TreeSearch_tbr_moves, 1},
  {"_TreeSearch_preorder_morphy",   (DL_FUNC) &_TreeSearch_preorder_morphy, 2},
  {"_TreeSearch_preorder_morphy_by_char",   (DL_FUNC) &_TreeSearch_preorder_morphy_by_char, 2},
  
  {"_TreeSearch_morphy_iw",         (DL_FUNC) &_TreeSearch_morphy_iw, 7},
  {"_TreeSearch_morphy_profile",    (DL_FUNC) &_TreeSearch_morphy_profile, 6},
  {"_TreeSearch_expected_mi",       (DL_FUNC) &_TreeSearch_expected_mi, 2},
  {"_TreeSearch_mi_key",            (DL_FUNC) &_TreeSearch_mi_key, 2},

  // {"_TreeSearch_astar_search_r",    (DL_FUNC) &_TreeSearch_astar_search_r, 3},
  {"_TreeSearch_quartet_concordance",(DL_FUNC) &_TreeSearch_quartet_concordance, 2},
  {"_TreeSearch_ts_fitch_score",    (DL_FUNC) &_TreeSearch_ts_fitch_score, 7},
  {"_TreeSearch_ts_na_char_steps", (DL_FUNC) &_TreeSearch_ts_na_char_steps, 5},
  {"_TreeSearch_ts_na_debug_char", (DL_FUNC) &_TreeSearch_ts_na_debug_char, 6},
  {"_TreeSearch_ts_nni_search",    (DL_FUNC) &_TreeSearch_ts_nni_search, 8},
  {"_TreeSearch_ts_debug_clip",    (DL_FUNC) &_TreeSearch_ts_debug_clip, 6},
  {"_TreeSearch_ts_spr_search",    (DL_FUNC) &_TreeSearch_ts_spr_search, 8},
  {"_TreeSearch_ts_tbr_search",    (DL_FUNC) &_TreeSearch_ts_tbr_search, 10},
  {"_TreeSearch_ts_drift_search",  (DL_FUNC) &_TreeSearch_ts_drift_search, 11},
  {"_TreeSearch_ts_test_indirect", (DL_FUNC) &_TreeSearch_ts_test_indirect, 8},
  {"_TreeSearch_ts_ratchet_search", (DL_FUNC) &_TreeSearch_ts_ratchet_search, 10},
  {"_TreeSearch_ts_compute_splits", (DL_FUNC) &_TreeSearch_ts_compute_splits, 2},
  {"_TreeSearch_ts_trees_equal",   (DL_FUNC) &_TreeSearch_ts_trees_equal, 3},
  {"_TreeSearch_ts_pool_test",     (DL_FUNC) &_TreeSearch_ts_pool_test, 5},
  {"_TreeSearch_ts_wagner_tree",   (DL_FUNC) &_TreeSearch_ts_wagner_tree, 7},
  {"_TreeSearch_ts_random_wagner_tree", (DL_FUNC) &_TreeSearch_ts_random_wagner_tree, 6},
  {"_TreeSearch_ts_sector_diag", (DL_FUNC) &_TreeSearch_ts_sector_diag, 6},
  {"_TreeSearch_ts_rss_search",  (DL_FUNC) &_TreeSearch_ts_rss_search, 13},
  {"_TreeSearch_ts_xss_search",  (DL_FUNC) &_TreeSearch_ts_xss_search, 12},
  {"_TreeSearch_ts_tree_fuse",   (DL_FUNC) &_TreeSearch_ts_tree_fuse, 11},
  {"_TreeSearch_ts_driven_search", (DL_FUNC) &_TreeSearch_ts_driven_search, 17},

  {"MORPHYLENGTH",                  (DL_FUNC) &MORPHYLENGTH, 4},
  {"RANDOM_TREE",                   (DL_FUNC) &RANDOM_TREE, 1},
  {"RANDOM_TREE_SCORE",             (DL_FUNC) &RANDOM_TREE_SCORE, 2},
  {NULL, NULL, 0}
};

void R_init_TreeSearch(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
