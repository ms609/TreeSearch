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

static const R_CallMethodDef callMethods[] = {
  {"_R_wrap_mpl_new_Morphy",        (DL_FUNC) &_R_wrap_mpl_new_Morphy, 0},
  {"_R_wrap_mpl_delete_Morphy",     (DL_FUNC) &_R_wrap_mpl_delete_Morphy, 1},
  {"_R_wrap_mpl_init_Morphy",       (DL_FUNC) &_R_wrap_mpl_init_Morphy, 3},
  {"_R_wrap_mpl_get_numtaxa",       (DL_FUNC) &_R_wrap_mpl_get_numtaxa, 1},
  {"_R_wrap_mpl_get_num_charac",    (DL_FUNC) &_R_wrap_mpl_get_num_charac, 1},
  {"_R_wrap_mpl_attach_symbols",    (DL_FUNC) &_R_wrap_mpl_attach_symbols, 2},
  {"_R_wrap_mpl_get_symbols",       (DL_FUNC) &_R_wrap_mpl_get_symbols, 1},
  {"_R_wrap_mpl_attach_rawdata",    (DL_FUNC) &_R_wrap_mpl_attach_rawdata, 2},
  {"_R_wrap_mpl_delete_rawdata",    (DL_FUNC) &_R_wrap_mpl_delete_rawdata, 1},
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
