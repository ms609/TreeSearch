// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// preorder_morphy
int preorder_morphy(IntegerMatrix edge, SEXP MorphyHandl);
RcppExport SEXP _TreeSearch_preorder_morphy(SEXP edgeSEXP, SEXP MorphyHandlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type MorphyHandl(MorphyHandlSEXP);
    rcpp_result_gen = Rcpp::wrap(preorder_morphy(edge, MorphyHandl));
    return rcpp_result_gen;
END_RCPP
}
// preorder_morphy_by_char
IntegerVector preorder_morphy_by_char(IntegerMatrix edge, List MorphyHandls);
RcppExport SEXP _TreeSearch_preorder_morphy_by_char(SEXP edgeSEXP, SEXP MorphyHandlsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< List >::type MorphyHandls(MorphyHandlsSEXP);
    rcpp_result_gen = Rcpp::wrap(preorder_morphy_by_char(edge, MorphyHandls));
    return rcpp_result_gen;
END_RCPP
}
// morphy_iw
double morphy_iw(IntegerMatrix edge, List MorphyHandls, NumericVector weight, IntegerVector minScore, IntegerVector sequence, NumericVector concavity, NumericVector target);
RcppExport SEXP _TreeSearch_morphy_iw(SEXP edgeSEXP, SEXP MorphyHandlsSEXP, SEXP weightSEXP, SEXP minScoreSEXP, SEXP sequenceSEXP, SEXP concavitySEXP, SEXP targetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< List >::type MorphyHandls(MorphyHandlsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type minScore(minScoreSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sequence(sequenceSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type concavity(concavitySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type target(targetSEXP);
    rcpp_result_gen = Rcpp::wrap(morphy_iw(edge, MorphyHandls, weight, minScore, sequence, concavity, target));
    return rcpp_result_gen;
END_RCPP
}
// morphy_profile
double morphy_profile(IntegerMatrix edge, List MorphyHandls, NumericVector weight, IntegerVector sequence, NumericMatrix profiles, NumericVector target);
RcppExport SEXP _TreeSearch_morphy_profile(SEXP edgeSEXP, SEXP MorphyHandlsSEXP, SEXP weightSEXP, SEXP sequenceSEXP, SEXP profilesSEXP, SEXP targetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< List >::type MorphyHandls(MorphyHandlsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sequence(sequenceSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type profiles(profilesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type target(targetSEXP);
    rcpp_result_gen = Rcpp::wrap(morphy_profile(edge, MorphyHandls, weight, sequence, profiles, target));
    return rcpp_result_gen;
END_RCPP
}
// nni
IntegerMatrix nni(const IntegerMatrix edge, const IntegerVector randomEdge, const IntegerVector whichSwitch);
RcppExport SEXP _TreeSearch_nni(SEXP edgeSEXP, SEXP randomEdgeSEXP, SEXP whichSwitchSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type randomEdge(randomEdgeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type whichSwitch(whichSwitchSEXP);
    rcpp_result_gen = Rcpp::wrap(nni(edge, randomEdge, whichSwitch));
    return rcpp_result_gen;
END_RCPP
}
// spr_moves
IntegerMatrix spr_moves(const IntegerMatrix edge);
RcppExport SEXP _TreeSearch_spr_moves(SEXP edgeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type edge(edgeSEXP);
    rcpp_result_gen = Rcpp::wrap(spr_moves(edge));
    return rcpp_result_gen;
END_RCPP
}
// spr
IntegerMatrix spr(const IntegerMatrix edge, const IntegerVector move);
RcppExport SEXP _TreeSearch_spr(SEXP edgeSEXP, SEXP moveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type move(moveSEXP);
    rcpp_result_gen = Rcpp::wrap(spr(edge, move));
    return rcpp_result_gen;
END_RCPP
}
// tbr
IntegerMatrix tbr(const IntegerMatrix edge, const IntegerVector move);
RcppExport SEXP _TreeSearch_tbr(SEXP edgeSEXP, SEXP moveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type move(moveSEXP);
    rcpp_result_gen = Rcpp::wrap(tbr(edge, move));
    return rcpp_result_gen;
END_RCPP
}
// asan_error
List asan_error(const IntegerMatrix x);
RcppExport SEXP _TreeSearch_asan_error(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(asan_error(x));
    return rcpp_result_gen;
END_RCPP
}
// all_spr
List all_spr(const IntegerMatrix edge, const IntegerVector break_order);
RcppExport SEXP _TreeSearch_all_spr(SEXP edgeSEXP, SEXP break_orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type break_order(break_orderSEXP);
    rcpp_result_gen = Rcpp::wrap(all_spr(edge, break_order));
    return rcpp_result_gen;
END_RCPP
}
// all_tbr
List all_tbr(const IntegerMatrix edge, const IntegerVector break_order);
RcppExport SEXP _TreeSearch_all_tbr(SEXP edgeSEXP, SEXP break_orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type break_order(break_orderSEXP);
    rcpp_result_gen = Rcpp::wrap(all_tbr(edge, break_order));
    return rcpp_result_gen;
END_RCPP
}
