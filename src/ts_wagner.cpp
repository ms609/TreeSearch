#include "ts_wagner.h"
#include "ts_fitch.h"
#include <algorithm>
#include <numeric>
#include <climits>
#include <vector>
#include <R.h>  // GetRNGstate, PutRNGstate, unif_rand

namespace ts {

namespace {

// Initialize a full-sized TreeState for Wagner construction.
// Allocates all arrays but does not set topology (caller does that).
void init_wagner_state(TreeState& tree, const DataSet& ds) {
  tree.n_tip = ds.n_tips;
  tree.n_internal = ds.n_tips - 1;
  tree.n_node = 2 * ds.n_tips - 1;
  tree.total_words = ds.total_words;
  tree.n_blocks = ds.n_blocks;

  tree.parent.assign(tree.n_node, -1);
  tree.left.assign(tree.n_internal, -1);
  tree.right.assign(tree.n_internal, -1);

  size_t state_size = static_cast<size_t>(tree.n_node) * tree.total_words;
  tree.prelim.assign(state_size, 0ULL);
  tree.final_.assign(state_size, 0ULL);
  tree.down2.assign(state_size, 0ULL);
  tree.subtree_actives.assign(state_size, 0ULL);
  tree.local_cost.assign(
      static_cast<size_t>(tree.n_node) * tree.n_blocks, 0ULL);

  tree.load_tip_states(ds);
}

// Build the initial 3-taxon tree:
//
//        root (= n_tip)
//       /    \
//    int1     t2
//  (n_tip+1)
//   /    \
// t0      t1
void build_three_taxon_tree(TreeState& tree, int t0, int t1, int t2) {
  int root = tree.n_tip;
  int int1 = tree.n_tip + 1;

  tree.parent[root] = root;
  tree.parent[int1] = root;
  tree.parent[t2] = root;
  tree.parent[t0] = int1;
  tree.parent[t1] = int1;

  tree.left[0] = int1;  // root's left
  tree.right[0] = t2;   // root's right
  tree.left[1] = t0;    // int1's left
  tree.right[1] = t1;   // int1's right
}

// Insert tip at edge (above, below), creating new_internal between them.
//
// Before:     above          After:     above
//               |                        |
//             below                   new_internal
//                                     /        \
//                                   tip       below
void insert_tip_at_edge(TreeState& tree, int tip, int new_internal,
                        int above, int below) {
  int ni = new_internal - tree.n_tip;

  tree.parent[new_internal] = above;
  tree.left[ni] = tip;
  tree.right[ni] = below;
  tree.parent[tip] = new_internal;
  tree.parent[below] = new_internal;

  // Update above's child pointer
  int ai = above - tree.n_tip;
  if (tree.left[ai] == below) {
    tree.left[ai] = new_internal;
  } else {
    tree.right[ai] = new_internal;
  }
}

} // anonymous namespace

WagnerResult wagner_tree(TreeState& tree, const DataSet& ds,
                         const std::vector<int>& addition_order) {
  int n_tip = ds.n_tips;

  // Validate or create addition order
  std::vector<int> order;
  if (addition_order.empty()) {
    order.resize(n_tip);
    std::iota(order.begin(), order.end(), 0);
  } else {
    order = addition_order;
  }

  // Initialize full-sized TreeState
  init_wagner_state(tree, ds);

  // Build initial 3-taxon tree
  build_three_taxon_tree(tree, order[0], order[1], order[2]);
  tree.build_postorder();
  double score = score_tree(tree, ds);

  // Add remaining taxa one at a time
  for (int i = 3; i < n_tip; ++i) {
    int tip = order[i];
    int new_internal = n_tip + i - 1;

    const uint64_t* tip_prelim =
        &ds.tip_states[static_cast<size_t>(tip) * ds.total_words];

    // Find best insertion edge via DFS from root
    int best_above = -1, best_below = -1;
    int best_extra = INT_MAX;

    std::vector<int> stack;
    stack.push_back(n_tip);  // root

    while (!stack.empty()) {
      int node = stack.back();
      stack.pop_back();

      if (node < n_tip) continue;  // tip — no children to enumerate

      int ni = node - n_tip;
      int lc = tree.left[ni];
      int rc = tree.right[ni];

      // Evaluate edge (node, lc)
      int extra = fitch_indirect_length(tip_prelim, tree, ds, node, lc);
      if (extra < best_extra) {
        best_extra = extra;
        best_above = node;
        best_below = lc;
      }

      // Evaluate edge (node, rc)
      extra = fitch_indirect_length(tip_prelim, tree, ds, node, rc);
      if (extra < best_extra) {
        best_extra = extra;
        best_above = node;
        best_below = rc;
      }

      stack.push_back(lc);
      stack.push_back(rc);
    }

    // Insert tip at the best edge
    insert_tip_at_edge(tree, tip, new_internal, best_above, best_below);

    // Rebuild postorder and rescore
    tree.build_postorder();
    score = score_tree(tree, ds);
  }

  WagnerResult result;
  result.score = score;
  return result;
}

WagnerResult random_wagner_tree(TreeState& tree, const DataSet& ds) {
  int n_tips = ds.n_tips;
  std::vector<int> order(n_tips);
  std::iota(order.begin(), order.end(), 0);

  // Fisher-Yates shuffle using R's RNG
  GetRNGstate();
  for (int i = n_tips - 1; i > 0; --i) {
    int j = static_cast<int>(unif_rand() * (i + 1));
    if (j > i) j = i;  // guard against floating-point edge case
    std::swap(order[i], order[j]);
  }
  PutRNGstate();

  return wagner_tree(tree, ds, order);
}

} // namespace ts
