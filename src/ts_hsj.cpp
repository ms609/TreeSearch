#include "ts_hsj.h"
#include "ts_fitch.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <limits>

namespace ts {

std::vector<int> partition_weights(
    const int* index_r, int n_orig_chars,
    const int* weight_r, int n_patterns,
    const std::vector<int>& hierarchy_chars)
{
  std::vector<int> adjusted(weight_r, weight_r + n_patterns);

  for (int c : hierarchy_chars) {
    if (c < 0 || c >= n_orig_chars) continue;
    int pat = index_r[c];
    if (pat >= 0 && pat < n_patterns && adjusted[pat] > 0) {
      --adjusted[pat];
    }
  }

  return adjusted;
}

// Fitch downpass + uppass for a single character represented as integer state
// labels.  After the downpass, state sets at internal nodes are ambiguous
// (intersection or union of children).  The uppass resolves each node to a
// single state so that parent–child mismatches can be detected for HSJ
// secondary dissimilarity.
// Returns number of Fitch steps (union operations in the downpass).
static int fitch_label_char(
    const TreeState& tree,
    const std::vector<int>& tip_labels,
    int char_idx,
    int n_orig_chars,
    int inapp_state,
    std::vector<uint32_t>& state_sets)
{
  int n_tip = tree.n_tip;

  // Initialize tips
  for (int t = 0; t < n_tip; ++t) {
    int label = tip_labels[t * n_orig_chars + char_idx];
    if (label < 0 || label > 30) {
      // Ambiguous: all states
      state_sets[t] = 0xFFFFFFFFu;
    } else {
      state_sets[t] = 1u << label;
    }
  }

  // --- Downpass ---
  int steps = 0;
  for (int i = 0; i < static_cast<int>(tree.postorder.size()); ++i) {
    int node = tree.postorder[i];
    int ni = node - n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];

    uint32_t inter = state_sets[lc] & state_sets[rc];
    if (inter != 0) {
      state_sets[node] = inter;
    } else {
      state_sets[node] = state_sets[lc] | state_sets[rc];
      ++steps;
    }
  }

  // --- Uppass: resolve each node to a single state ---
  // Root: pick lowest bit from its state set
  int root = tree.postorder.back();
  state_sets[root] = state_sets[root] & (~state_sets[root] + 1); // isolate lowest bit

  // Traverse preorder (reverse postorder) to resolve internal nodes and tips
  for (int i = static_cast<int>(tree.postorder.size()) - 1; i >= 0; --i) {
    int node = tree.postorder[i];
    int ni = node - n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];

    // Resolve each child: prefer parent's state if compatible
    for (int child : {lc, rc}) {
      uint32_t overlap = state_sets[child] & state_sets[node];
      if (overlap != 0) {
        state_sets[child] = state_sets[node]; // inherit parent's state
      } else {
        // Pick lowest bit from child's state set
        state_sets[child] = state_sets[child] & (~state_sets[child] + 1);
      }
    }
  }

  return steps;
}

// Check if two state sets have disjoint states (= mismatch for HSJ).
static inline bool states_mismatch(uint32_t a, uint32_t b) {
  return (a & b) == 0;
}

// Score one hierarchy block via the HSJ a(n)/p(n) dynamic programming.
//
// For each internal node n with children c1, c2:
//   a(n) = min cost assuming primary ABSENT at n
//   p(n) = min cost assuming primary PRESENT at n
//
// Branch cost(parent_state, child_state):
//   absent → absent:   0
//   absent → present:  1  (gain)
//   present → absent:  1  (loss)
//   present → present: alpha * d(parent, child) / m
//     where d = number of secondaries with disjoint state sets
static double score_hierarchy_block(
    const TreeState& tree,
    const HierarchyBlock& block,
    double alpha,
    const std::vector<int>& tip_labels,
    int n_orig_chars,
    int inapp_state)
{
  const int n_tip = tree.n_tip;
  const int n_node = tree.n_node;
  const int m = block.n_secondaries;
  const double INF = std::numeric_limits<double>::infinity();

  // Step 1: Run Fitch downpass on each secondary character
  // sec_states[j * n_node + node] = bitmask of possible states
  std::vector<uint32_t> sec_states(m * n_node, 0);
  {
    std::vector<uint32_t> buf(n_node);
    for (int j = 0; j < m; ++j) {
      fitch_label_char(tree, tip_labels, block.secondary_chars[j],
                       n_orig_chars, inapp_state, buf);
      for (int nd = 0; nd < n_node; ++nd) {
        sec_states[j * n_node + nd] = buf[nd];
      }
    }
  }

  // Step 2: Determine primary state at each tip
  // primary_present[tip] = true if primary is present (not absent_state)
  std::vector<bool> primary_present(n_tip, false);
  for (int t = 0; t < n_tip; ++t) {
    int label = tip_labels[t * n_orig_chars + block.primary_char];
    primary_present[t] = (label != block.absent_state);
  }

  // Step 3: a(n)/p(n) DP
  // a[node], p[node]
  std::vector<double> a(n_node, 0.0);
  std::vector<double> p(n_node, 0.0);

  // Initialize leaves
  for (int t = 0; t < n_tip; ++t) {
    if (primary_present[t]) {
      a[t] = INF;
      p[t] = 0.0;
    } else {
      a[t] = 0.0;
      p[t] = INF;
    }
  }

  // Helper: count secondary mismatches between two nodes
  auto count_mismatches = [&](int node1, int node2) -> int {
    int d = 0;
    for (int j = 0; j < m; ++j) {
      if (states_mismatch(sec_states[j * n_node + node1],
                          sec_states[j * n_node + node2])) {
        ++d;
      }
    }
    return d;
  };

  // Postorder traversal
  for (int i = 0; i < static_cast<int>(tree.postorder.size()); ++i) {
    int node = tree.postorder[i];
    int ni = node - n_tip;
    int c1 = tree.left[ni];
    int c2 = tree.right[ni];

    // Compute branch costs for each combination of parent/child states
    // parent absent:
    double bc_aa_c1 = 0.0;  // absent→absent
    double bc_ap_c1 = 1.0;  // absent→present (gain)
    double bc_aa_c2 = 0.0;
    double bc_ap_c2 = 1.0;

    // parent present:
    double bc_pa_c1 = 1.0;  // present→absent (loss)
    double bc_pa_c2 = 1.0;
    double bc_pp_c1 = 0.0;  // present→present (secondary dissimilarity)
    double bc_pp_c2 = 0.0;

    if (m > 0 && alpha > 0.0) {
      int d1 = count_mismatches(node, c1);
      int d2 = count_mismatches(node, c2);
      bc_pp_c1 = alpha * d1 / m;
      bc_pp_c2 = alpha * d2 / m;
    }

    // a(n): parent is absent
    double best_a = INF;
    // Try all 4 combinations of child states
    double cost_aa = (a[c1] + bc_aa_c1) + (a[c2] + bc_aa_c2);  // both absent
    double cost_ap = (a[c1] + bc_aa_c1) + (p[c2] + bc_ap_c2);  // c1 abs, c2 pres
    double cost_pa = (p[c1] + bc_ap_c1) + (a[c2] + bc_aa_c2);  // c1 pres, c2 abs
    double cost_pp = (p[c1] + bc_ap_c1) + (p[c2] + bc_ap_c2);  // both present
    best_a = std::min({cost_aa, cost_ap, cost_pa, cost_pp});
    a[node] = best_a;

    // p(n): parent is present
    double best_p = INF;
    cost_aa = (a[c1] + bc_pa_c1) + (a[c2] + bc_pa_c2);  // both absent
    cost_ap = (a[c1] + bc_pa_c1) + (p[c2] + bc_pp_c2);  // c1 abs, c2 pres
    cost_pa = (p[c1] + bc_pp_c1) + (a[c2] + bc_pa_c2);  // c1 pres, c2 abs
    cost_pp = (p[c1] + bc_pp_c1) + (p[c2] + bc_pp_c2);  // both present
    best_p = std::min({cost_aa, cost_ap, cost_pa, cost_pp});
    p[node] = best_p;
  }

  // Root score = min(a[root], p[root])
  int root = tree.postorder.back();
  return std::min(a[root], p[root]);
}

double hsj_score(
    TreeState& tree,
    const DataSet& ds,
    const std::vector<HierarchyBlock>& hierarchy_blocks,
    double alpha,
    const std::vector<int>& tip_labels,
    int n_orig_chars)
{
  // Score non-hierarchy characters via standard Fitch.
  // Use fitch_score_ew() directly to avoid infinite recursion when
  // score_tree() dispatches to hsj_score() for HSJ mode.
  double fitch_total = fitch_score_ew(tree, ds);

  // Score each hierarchy block via HSJ DP
  double hsj_total = 0.0;
  for (const auto& block : hierarchy_blocks) {
    hsj_total += score_hierarchy_block(
        tree, block, alpha, tip_labels, n_orig_chars, ds.inapp_state);
  }

  return fitch_total + hsj_total;
}

double hsj_score(TreeState& tree, const DataSet& ds)
{
  return hsj_score(tree, ds, ds.hierarchy_blocks, ds.hsj_alpha,
                   ds.tip_labels, ds.n_orig_chars);
}

} // namespace ts
