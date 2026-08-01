#include "ts_soft_sankoff.h"

namespace ts {

// =========================================================================
// Single-character downpass
// =========================================================================

double soft_sankoff_score_char(
    const int* left, const int* right,
    const int* postorder, int n_internal,
    int n_tip,
    const SoftSankoffChar& sc,
    const double* tip_costs_ch,
    int tip_stride,
    double temperature,
    double* node_costs_out)
{
  const int ns = sc.n_states;
  const int n_node = n_tip + n_internal;
  const double INF = std::numeric_limits<double>::infinity();

  // Degenerate 1-tip "tree": no root to read from postorder, and a single node
  // always admits some state.  Mirrors the guard in sankoff_score_char().
  if (n_internal <= 0) return 0.0;

  std::vector<double> local_buf;
  double* costs;
  if (node_costs_out) {
    costs = node_costs_out;
  } else {
    local_buf.resize(static_cast<size_t>(n_node) * ns, INF);
    costs = local_buf.data();
  }

  // Tip costs
  for (int t = 0; t < n_tip; ++t) {
    double* dst = costs + static_cast<size_t>(t) * ns;
    const double* src = tip_costs_ch + static_cast<size_t>(t) * tip_stride;
    for (int s = 0; s < ns; ++s) {
      dst[s] = src[s];
    }
  }

  const bool per_branch = !sc.branch_costs.empty();
  const size_t cm_size = static_cast<size_t>(ns) * ns;

  // Scratch for the inner reduction: cost_matrix[s][t] + child_cost[t] over t.
  std::vector<double> work(ns);

  for (int i = 0; i < n_internal; ++i) {
    const int node = postorder[i];
    const int ni = node - n_tip;
    const int children[2] = {left[ni], right[ni]};

    double* nc = costs + static_cast<size_t>(node) * ns;
    for (int s = 0; s < ns; ++s) nc[s] = 0.0;

    for (int ci = 0; ci < 2; ++ci) {
      const int child = children[ci];
      const double* cc = costs + static_cast<size_t>(child) * ns;
      // Per-branch matrices are indexed by the child node, so each edge picks
      // up its own costs.  An off-by-one here does not error, it yields inf.
      const double* cm = per_branch
          ? sc.branch_costs.data() + static_cast<size_t>(child) * cm_size
          : sc.cost_matrix.data();

      for (int s = 0; s < ns; ++s) {
        const double* cm_row = cm + static_cast<size_t>(s) * ns;
        for (int t = 0; t < ns; ++t) {
          work[t] = cm_row[t] + cc[t];
        }
        nc[s] += softmin(work.data(), ns, temperature);
      }
    }
  }

  // Root reduction.  Unlike the hard kernel this is a *soft* reduction, and it
  // carries the optional additive root cost -log(pi).
  const int root = postorder[n_internal - 1];
  const double* root_costs = costs + static_cast<size_t>(root) * ns;

  if (sc.root_cost.empty()) {
    return softmin(root_costs, ns, temperature);
  }

  std::vector<double> root_total(ns);
  for (int s = 0; s < ns; ++s) {
    root_total[s] = root_costs[s] + sc.root_cost[s];
  }
  return softmin(root_total.data(), ns, temperature);
}

// =========================================================================
// Multi-character scoring
// =========================================================================

double soft_sankoff_score(
    const int* left, const int* right,
    const int* postorder, int n_internal,
    int n_tip,
    const SoftSankoffData& sd,
    double temperature)
{
  double total = 0.0;
  const int stride = sd.stride();

  for (int ch = 0; ch < sd.n_chars; ++ch) {
    const double* ch_tip = sd.tip_costs.data() + ch * sd.max_states;
    total += soft_sankoff_score_char(
        left, right, postorder, n_internal, n_tip,
        sd.chars[ch], ch_tip, stride, temperature);
  }

  return total;
}

} // namespace ts
