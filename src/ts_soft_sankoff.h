#ifndef TS_SOFT_SANKOFF_H
#define TS_SOFT_SANKOFF_H

// Soft-Sankoff scoring: Sankoff's dynamic program with `min` replaced by a
// soft-minimum at temperature T.
//
//   softmin_T(x) = -T * log( sum_j exp(-x_j / T) )
//   S_v(i)       = sum over children c of  softmin_T over j of
//                    [ cost_c(i, j) + S_c(j) ]
//
// T -> 0 recovers hard weighted parsimony; T = 1 with cost = -log P_ij(t)
// recovers Felsenstein's pruning algorithm exactly.  See
// dev/plans/2026-08-01-soft-sankoff-temperature-dial.md.
//
// EXPLORATION CODE.  This is deliberately a separate kernel from
// ts_sankoff.{h,cpp}: that struct is on the live x-transformation path, which
// carries open rooting defects (T-374, T-385), and nothing here should be
// built on top of an unsettled objective.  Nothing in this file is reachable
// from MaximizeParsimony() or any default scoring path.
//
// Like ts_sankoff.h, the implementation is FULL-RESCORE ONLY: there is no
// incremental variant.  It therefore prices the *constant factor* of soft
// scoring, not the loss of incremental subtree rescoring, which is the term
// that changes inner-loop complexity rather than a constant.

#include <cmath>
#include <limits>
#include <vector>

namespace ts {

// ---------------------------------------------------------------------------
// Soft minimum
// ---------------------------------------------------------------------------

// Numerically stable soft-minimum of x[0 .. n-1], shifted by the hard minimum
// so every exponent lies in [-inf, 0] and the sum lies in [1, n].
//
// `temperature <= 0` returns the hard minimum.  This is an internal branch
// rather than a delegation to ts::sankoff_score_char(), so that a T -> 0 test
// against the hard kernel compares two independent implementations rather than
// a function to itself.
//
// +infinity entries contribute exp(-inf) = 0 and need no filtering.  The one
// case needing a guard is *every* entry being non-finite, where the shifted
// form would evaluate inf - inf = NaN; there we return the minimum unchanged.
// Infinite tip costs are the normal case (0 for an observed state, +inf
// otherwise), not an edge case.
inline double softmin(const double* x, int n, double temperature) {
  double m = std::numeric_limits<double>::infinity();
  for (int i = 0; i < n; ++i) {
    if (x[i] < m) m = x[i];
  }
  if (temperature <= 0.0) return m;
  if (!std::isfinite(m)) return m;

  double acc = 0.0;
  for (int i = 0; i < n; ++i) {
    acc += std::exp(-(x[i] - m) / temperature);
  }
  return m - temperature * std::log(acc);
}

// ---------------------------------------------------------------------------
// Data
// ---------------------------------------------------------------------------

// One soft-Sankoff character.
struct SoftSankoffChar {
  int n_states;

  // Cost matrix, row-major, cost_matrix[from * n_states + to].  Applied to
  // every edge when `branch_costs` is empty.
  std::vector<double> cost_matrix;

  // Optional per-branch cost matrices, indexed by CHILD node (each edge is
  // uniquely named by the node below it), matching the R reference's
  // `cost[[child]]` convention:
  //   branch_costs[child * n_states * n_states + from * n_states + to]
  // Length n_node * n_states^2 when used; the root's slot is never read.
  // Needed for the T = 1 Felsenstein identity, where each branch carries its
  // own -log P_ij(t).
  std::vector<double> branch_costs;

  // Optional additive per-state cost at the root, length n_states.  Supply
  // -log(pi) to reproduce a likelihood with root frequencies pi; leave empty
  // for the parsimony convention of no root cost.
  std::vector<double> root_cost;
};

// Characters + tip costs.  Tip cost layout matches ts::SankoffData:
//   tip_costs[tip * stride() + ch * max_states + state]
// = 0 if the state is observed at that tip, +infinity if not.
struct SoftSankoffData {
  int n_tips;
  int n_chars;
  int max_states;
  std::vector<SoftSankoffChar> chars;
  std::vector<double> tip_costs;

  int stride() const { return n_chars * max_states; }
};

// ---------------------------------------------------------------------------
// Scoring (downpass only)
// ---------------------------------------------------------------------------

// Score one character at `temperature`.
//
// Topology arrays follow ts::sankoff_score_char():
//   left[i], right[i]: children of internal node (n_tip + i)
//   postorder[i]: internal nodes, leaves-to-root; postorder.back() is the root
//
// tip_costs_ch[tip * tip_stride + state], tip_stride >= sc.n_states.
//
// node_costs: if non-null, filled with per-node per-state conditional costs,
//   node_costs[node * sc.n_states + state]; tip rows are copies of the tip
//   costs.  Sized n_node * n_states by the caller.
//
// Returns softmin_T over s of (root_costs[s] + root_cost[s]).
double soft_sankoff_score_char(
    const int* left, const int* right,
    const int* postorder, int n_internal,
    int n_tip,
    const SoftSankoffChar& sc,
    const double* tip_costs_ch,
    int tip_stride,
    double temperature,
    double* node_costs = nullptr);

// Total soft score, summed over characters.
double soft_sankoff_score(
    const int* left, const int* right,
    const int* postorder, int n_internal,
    int n_tip,
    const SoftSankoffData& sd,
    double temperature);

} // namespace ts

#endif // TS_SOFT_SANKOFF_H
