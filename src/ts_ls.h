#ifndef TS_LS_H
#define TS_LS_H

// Least-squares (LS) branch-length fitting and topology search.
//
// Given a target symmetric dissimilarity matrix D over the tips and a fixed
// topology, fit branch lengths v that minimise the (optionally weighted)
// residual sum of squares
//
//   RSS = sum_{i<j} w_ij * ( d_tree(i,j) - D_ij )^2
//
// where d_tree(i,j) is the patristic (path-length) distance between tips i and
// j under the fitted lengths.  Two fitting modes are offered:
//
//   OLS  — ordinary least squares (closed form via the path design matrix and
//          the normal equations).  Branch lengths may be negative.
//   NNLS — non-negative least squares (Lawson & Hanson active-set), matching
//          phangorn::nnls.tree().  Branch lengths are constrained >= 0.
//
// The design matrix is built over the *unrooted* branches of the tree
// (2*n_tip - 3 of them: n_tip pendant edges + n_tip-3 internal splits).  The
// two edges incident to the rooted TreeState's root describe a single unrooted
// branch and are merged into one design column, so the fit matches the
// unrooted convention used by phangorn.
//
// This path is entirely independent of the parsimony/Fitch machinery: it only
// reads the tree topology (parent/left/right + postorder) and never touches
// the per-node character state arrays, so it operates on a TreeState built with
// total_words == 0.

#include "ts_tree.h"
#include <vector>
#include <cstdint>

namespace ts {

enum class LSMethod { OLS, NNLS };

// Target dissimilarities, flattened to the canonical pair order
//   for i in 0..n-1: for j in i+1..n-1   ->   pair index p
struct LSData {
  int n_tip = 0;
  int n_pairs = 0;                 // n_tip*(n_tip-1)/2
  std::vector<double> target;      // [n_pairs]  D_ij
  std::vector<double> sqrt_weight; // [n_pairs]  sqrt(w_ij); all 1 if unweighted
  bool weighted = false;
};

// Build LSData from a full n*n symmetric matrix (column-major == row-major
// because D is symmetric).  `wmat` may be null for unit weights; when given it
// is the per-pair weight matrix (e.g. Fitch-Margoliash 1/D^2), also symmetric.
LSData build_ls_data(const double* dmat, int n_tip, const double* wmat);

// Result of fitting branch lengths on a fixed topology.
struct LSFit {
  bool ok = false;                     // false if the linear solve failed
  double rss = 0.0;                    // residual sum of squares (weighted)
  int n_branch = 0;                    // 2*n_tip - 3
  std::vector<double> branch_length;   // [n_branch] fitted length per column
  std::vector<int> branch_node;        // [n_branch] child node id of each branch
};

// Fit branch lengths on the fixed topology `tree` (postorder must be current).
LSFit ls_fit(const TreeState& tree, const LSData& ls, LSMethod method);

// RSS only — used in the search hot loop.  Returns +Inf on solve failure.
double ls_score(const TreeState& tree, const LSData& ls, LSMethod method);

// NNI hill-climbing search minimising LS RSS.  Modifies `tree` in place to the
// best topology found.  First-improvement over a randomised edge order; repeats
// passes until no NNI improves the score (beyond `max_hits` equal-score moves).
struct LSSearchResult {
  double rss = 0.0;
  int n_moves = 0;        // improving (or accepted equal) moves applied
  int n_iterations = 0;   // candidate evaluations
};
LSSearchResult ls_nni_search(TreeState& tree, const LSData& ls,
                             LSMethod method, int max_hits);

// SPR hill-climbing search minimising LS RSS (full rescore per candidate).
LSSearchResult ls_spr_search(TreeState& tree, const LSData& ls,
                             LSMethod method, int max_hits);

} // namespace ts

#endif // TS_LS_H
