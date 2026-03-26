#include "ts_collapsed.h"
#include "ts_simd.h"
#include <cstring>

namespace ts {

void compute_collapsed_flags(
    const TreeState& tree,
    const DataSet& ds,
    std::vector<uint8_t>& collapsed) {

  collapsed.assign(tree.n_node, 0);

  // If all characters were simplified away (total_words == 0), there's
  // nothing to evaluate — leave all edges uncollapsed.
  if (tree.total_words == 0) return;

  // Detect whether any block has inapplicable characters.
  bool has_na = false;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) { has_na = true; break; }
  }

  const int nb = ds.n_blocks;
  const int tw = tree.total_words;
  const size_t word_bytes = static_cast<size_t>(tw) * sizeof(uint64_t);

  for (int c = 0; c < tree.n_node; ++c) {
    // Skip root (no parent edge).
    if (c == tree.n_tip) continue;
    int p = tree.parent[c];
    // Root's children can be clipped but removing one changes the root
    // structure; skip for safety.
    if (p == tree.n_tip) continue;

    int pi = p - tree.n_tip;  // internal index for left[]/right[]
    int s = (tree.left[pi] == c) ? tree.right[pi] : tree.left[pi];

    // --- Condition 1: zero standard-block cost at parent ---
    bool zero_std = true;
    for (int b = 0; b < nb && zero_std; ++b) {
      if (ds.blocks[b].has_inapplicable) continue;
      if (tree.local_cost[static_cast<size_t>(p) * nb + b])
        zero_std = false;
    }
    if (!zero_std) continue;

    // --- Condition 2: zero NA-block cost at parent ---
    if (has_na) {
      bool zero_na = true;
      int left = tree.left[pi];
      int right = tree.right[pi];
      size_t lb = static_cast<size_t>(left) * tw;
      size_t rb = static_cast<size_t>(right) * tw;
      size_t pb = static_cast<size_t>(p) * tw;

      for (int b = 0; b < nb && zero_na; ++b) {
        const CharBlock& blk = ds.blocks[b];
        if (!blk.has_inapplicable) continue;

        int k = blk.n_states;
        int off = ds.block_word_offset[b];

        // l_act = OR of applicable state words (states 1..k-1) for left
        uint64_t l_act = simd::or_reduce(
            &tree.subtree_actives[lb + off], k, 1);
        uint64_t r_act = simd::or_reduce(
            &tree.subtree_actives[rb + off], k, 1);

        if ((l_act & r_act) == 0) continue;  // auto-zero

        // Both subtrees have applicable tips — compute full condition.
        // ss_app = OR of applicable states at node p
        uint64_t ss_app = 0;
        for (int st = 1; st < k; ++st)
          ss_app |= tree.final_[pb + off + st];

        // any_isect = any D2 state intersection between children
        uint64_t any_isect = simd::any_hit_reduce(
            &tree.down2[lb + off], &tree.down2[rb + off], k);

        uint64_t needs_step =
            l_act & r_act & ~(ss_app & any_isect) & blk.active_mask;
        if (needs_step) zero_na = false;
      }
      if (!zero_na) continue;
    }

    // --- Condition 3: prelim[sibling] == prelim[parent] ---
    size_t sb = static_cast<size_t>(s) * tw;
    size_t pb = static_cast<size_t>(p) * tw;
    if (std::memcmp(&tree.prelim[sb], &tree.prelim[pb], word_bytes) != 0)
      continue;

    // --- Conditions 4–5 (NA only): down2 and subtree_actives preservation ---
    if (has_na) {
      if (std::memcmp(&tree.down2[sb], &tree.down2[pb], word_bytes) != 0)
        continue;
      if (std::memcmp(&tree.subtree_actives[sb],
                      &tree.subtree_actives[pb], word_bytes) != 0)
        continue;
    }

    collapsed[c] = 1;
  }
}

void compute_collapsed_regions(
    const TreeState& tree,
    const DataSet& ds,
    CollapsedRegions& info) {

  // Step 1: compute per-node collapsed flags.
  compute_collapsed_flags(tree, ds, info.collapsed);

  const int n_node = tree.n_node;
  info.region_id.assign(n_node, -1);
  info.n_collapsed = 0;
  info.n_regions = 0;

  // Count collapsed edges.
  for (int c = 0; c < n_node; ++c) {
    if (info.collapsed[c]) ++info.n_collapsed;
  }
  if (info.n_collapsed == 0) return;

  // Step 2: assign region IDs.
  //
  // A collapsed edge connects child c to parent p (when collapsed[c] == 1).
  // Two nodes share a region if they are connected by a collapsed edge.
  //
  // Process internal nodes in REVERSE postorder (= preorder for trees):
  // parents visited before children. When a parent creates or joins a
  // region, its children inherit the same region_id.
  //
  // The root itself is never collapsed (no parent edge) and root's children
  // are excluded by compute_collapsed_flags(), so root always has
  // region_id == -1.
  const auto& po = tree.postorder;
  for (int idx = static_cast<int>(po.size()) - 1; idx >= 0; --idx) {
    int node = po[idx];
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];

    // Process left child
    if (info.collapsed[lc]) {
      if (info.region_id[node] >= 0) {
        // Parent already has a region — child joins it.
        info.region_id[lc] = info.region_id[node];
      } else {
        // Start a new region for both parent and child.
        int rid = info.n_regions++;
        info.region_id[node] = rid;
        info.region_id[lc] = rid;
      }
    }

    // Process right child
    if (info.collapsed[rc]) {
      if (info.region_id[node] >= 0) {
        info.region_id[rc] = info.region_id[node];
      } else {
        int rid = info.n_regions++;
        info.region_id[node] = rid;
        info.region_id[rc] = rid;
      }
    }
  }
}

} // namespace ts
