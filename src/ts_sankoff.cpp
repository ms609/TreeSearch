#include "ts_sankoff.h"
#include <algorithm>
#include <cstring>

namespace ts {

// =========================================================================
// Single-character downpass scoring
// =========================================================================

double sankoff_score_char(
    const int* left, const int* right,
    const int* postorder, int n_internal,
    int n_tip,
    const SankoffChar& sc,
    const double* tip_costs_ch,
    int tip_stride,
    double* node_costs_out)
{
  const int ns = sc.n_states;
  const int n_node = n_tip + n_internal;
  const double INF = std::numeric_limits<double>::infinity();

  // Allocate per-node cost array (n_node * ns).
  // Use caller's buffer if provided, otherwise allocate locally.
  std::vector<double> local_buf;
  double* costs;
  if (node_costs_out) {
    costs = node_costs_out;
  } else {
    local_buf.resize(static_cast<size_t>(n_node) * ns, INF);
    costs = local_buf.data();
  }

  // Initialize tip costs
  for (int t = 0; t < n_tip; ++t) {
    double* dst = costs + static_cast<size_t>(t) * ns;
    const double* src = tip_costs_ch + static_cast<size_t>(t) * tip_stride;
    for (int s = 0; s < ns; ++s) {
      dst[s] = src[s];
    }
  }

  // Postorder traversal: compute costs at each internal node
  for (int i = 0; i < n_internal; ++i) {
    int node = postorder[i];
    int ni = node - n_tip;
    int c1 = left[ni];
    int c2 = right[ni];

    double* nc = costs + static_cast<size_t>(node) * ns;
    const double* cc1 = costs + static_cast<size_t>(c1) * ns;
    const double* cc2 = costs + static_cast<size_t>(c2) * ns;

    for (int s = 0; s < ns; ++s) {
      // min_t(cost_matrix[s][t] + cost[child][t]) for each child
      double best_c1 = INF;
      double best_c2 = INF;
      const double* cm_row = sc.cost_matrix.data() + static_cast<size_t>(s) * ns;

      for (int t = 0; t < ns; ++t) {
        double val1 = cm_row[t] + cc1[t];
        if (val1 < best_c1) best_c1 = val1;

        double val2 = cm_row[t] + cc2[t];
        if (val2 < best_c2) best_c2 = val2;
      }

      nc[s] = best_c1 + best_c2;
    }
  }

  // Root score
  int root = postorder[n_internal - 1];
  const double* root_costs = costs + static_cast<size_t>(root) * ns;

  if (sc.forced_root_state >= 0 && sc.forced_root_state < ns) {
    return root_costs[sc.forced_root_state];
  }

  double best = INF;
  for (int s = 0; s < ns; ++s) {
    if (root_costs[s] < best) best = root_costs[s];
  }
  return best;
}

// =========================================================================
// Multi-character scoring
// =========================================================================

double sankoff_score(
    const int* left, const int* right,
    const int* postorder, int n_internal,
    int n_tip,
    const SankoffData& sd)
{
  double total = 0.0;
  const int stride = sd.stride();

  for (int ch = 0; ch < sd.n_chars; ++ch) {
    // Pointer to this character's tip costs with appropriate stride
    // tip_costs layout: tip_costs[tip * stride + ch * max_states + state]
    // We pass a pointer offset to the first tip's data for this char,
    // with tip_stride = stride so the function skips correctly between tips.
    const double* ch_tip = sd.tip_costs.data() + ch * sd.max_states;

    total += sankoff_score_char(
        left, right, postorder, n_internal, n_tip,
        sd.chars[ch], ch_tip, stride);
  }

  return total;
}

// =========================================================================
// Uppass: optimal state reconstruction
// =========================================================================

void sankoff_uppass(
    const int* left, const int* right,
    const int* postorder, int n_internal,
    int n_tip,
    const SankoffChar& sc,
    const double* node_costs,
    int* optimal_states)
{
  const int ns = sc.n_states;
  const int n_node = n_tip + n_internal;
  const double INF = std::numeric_limits<double>::infinity();

  // --- Root assignment ---
  int root = postorder[n_internal - 1];
  const double* root_c = node_costs + static_cast<size_t>(root) * ns;

  if (sc.forced_root_state >= 0 && sc.forced_root_state < ns) {
    optimal_states[root] = sc.forced_root_state;
  } else {
    int best_s = 0;
    double best_v = root_c[0];
    for (int s = 1; s < ns; ++s) {
      if (root_c[s] < best_v) {
        best_v = root_c[s];
        best_s = s;
      }
    }
    optimal_states[root] = best_s;
  }

  // --- Preorder traversal (reverse postorder) ---
  for (int i = n_internal - 1; i >= 0; --i) {
    int node = postorder[i];
    int ni = node - n_tip;
    int c1 = left[ni];
    int c2 = right[ni];
    int parent_state = optimal_states[node];

    // Assign optimal state to each child
    const int children[2] = {c1, c2};
    for (int ci = 0; ci < 2; ++ci) {
      int child = children[ci];
      const double* child_c = node_costs + static_cast<size_t>(child) * ns;
      const double* cm_row =
          sc.cost_matrix.data() + static_cast<size_t>(parent_state) * ns;

      int best_s = 0;
      double best_v = cm_row[0] + child_c[0];
      for (int s = 1; s < ns; ++s) {
        double val = cm_row[s] + child_c[s];
        if (val < best_v) {
          best_v = val;
          best_s = s;
        }
      }
      optimal_states[child] = best_s;
    }
  }
}

} // namespace ts
