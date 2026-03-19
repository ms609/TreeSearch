#ifndef TS_SANKOFF_H
#define TS_SANKOFF_H

// Sankoff parsimony optimization for step-matrix characters.
//
// Used for Goloboff et al. (2021) x-transformation recoding of hierarchical
// (inapplicable) characters.  Each hierarchy block is recoded into a single
// Sankoff character with up to 2^n+1 states and an asymmetric cost matrix
// (gain:loss = n+1 : 1).
//
// The implementation is full-rescore only (no incremental variant).  When
// integrated with the search pipeline, Fitch scores non-hierarchy characters
// and Sankoff scores only the recoded hierarchy characters.

#include <cmath>
#include <limits>
#include <vector>

namespace ts {

// Cost matrix and root-forcing info for one Sankoff character.
struct SankoffChar {
  int n_states;                        // number of distinct states
  std::vector<double> cost_matrix;     // [n_states x n_states], row-major
                                       // cost_matrix[from * n_states + to]
  int forced_root_state;               // -1 = unconstrained, 0..n_states-1
};

// Complete data for Sankoff scoring: characters + tip costs.
struct SankoffData {
  int n_tips;
  int n_chars;                         // number of Sankoff characters
  int max_states;                      // max(chars[i].n_states)
  std::vector<SankoffChar> chars;

  // Per-tip per-character per-state cost, flat array:
  //   tip_costs[tip * stride + ch * max_states + state]
  // where stride = n_chars * max_states.
  // = 0.0 if state is observed at that tip, +infinity if not.
  std::vector<double> tip_costs;

  int stride() const { return n_chars * max_states; }
};

// ---------------------------------------------------------------------------
// Scoring (downpass)
// ---------------------------------------------------------------------------

// Score a tree under the Sankoff criterion for all characters in sd.
//
// Topology arrays:
//   left[i], right[i]: children of internal node (n_tip + i), i = 0..n_internal-1
//   postorder[i]: internal nodes in leaves-to-root order, length n_internal
//   n_tip: number of tips (0..n_tip-1)
//
// Returns total parsimony score (sum over all characters).
double sankoff_score(
    const int* left, const int* right,
    const int* postorder, int n_internal,
    int n_tip,
    const SankoffData& sd);

// Score a single Sankoff character.
//
// tip_costs_ch: per-tip costs for this character.  Layout:
//   tip_costs_ch[tip * tip_stride + state], 0..n_tip-1, 0..sc.n_states-1
//   where tip_stride >= sc.n_states (allows strided access into SankoffData).
//
// node_costs: if non-null, filled with per-node per-state downpass costs.
//   Layout: node_costs[node * sc.n_states + state], node = 0..n_node-1.
//   For tips, these are copies of tip_costs_ch.
//
// Returns root score for this character.
double sankoff_score_char(
    const int* left, const int* right,
    const int* postorder, int n_internal,
    int n_tip,
    const SankoffChar& sc,
    const double* tip_costs_ch,
    int tip_stride,
    double* node_costs = nullptr);

// ---------------------------------------------------------------------------
// Uppass (state reconstruction)
// ---------------------------------------------------------------------------

// Assign an optimal state to each node (root-to-leaves), given downpass costs.
//
// node_costs: per-node per-state costs from sankoff_score_char().
//   Layout: node_costs[node * sc.n_states + state]
//
// optimal_states: output, length n_node.  optimal_states[node] = chosen state.
//
// The uppass assigns root = argmin(node_costs[root][s]) (or forced_root_state),
// then propagates to children choosing the state that minimises
// cost_matrix[parent_state][child_state] + node_costs[child][child_state].
void sankoff_uppass(
    const int* left, const int* right,
    const int* postorder, int n_internal,
    int n_tip,
    const SankoffChar& sc,
    const double* node_costs,
    int* optimal_states);

} // namespace ts

#endif // TS_SANKOFF_H
