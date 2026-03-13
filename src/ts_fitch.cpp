#include "ts_fitch.h"
#include <vector>

namespace ts {

// --- Node-level helpers ---

int fitch_downpass_node(
    const uint64_t* left_state,
    const uint64_t* right_state,
    uint64_t* node_state,
    int n_states,
    uint64_t active_mask)
{
  uint64_t any_intersect = 0;
  for (int s = 0; s < n_states; ++s) {
    uint64_t isect = left_state[s] & right_state[s];
    any_intersect |= isect;
  }

  uint64_t needs_union = ~any_intersect & active_mask;
  int steps = popcount64(needs_union);

  for (int s = 0; s < n_states; ++s) {
    uint64_t isect = left_state[s] & right_state[s];
    uint64_t uni = left_state[s] | right_state[s];
    node_state[s] = (isect & any_intersect) | (uni & needs_union);
  }

  return steps;
}

// Compute final states at node from ancestor's final and node's prelim.
// Returns true if any final state word changed.
static bool uppass_node(TreeState& tree, const DataSet& ds, int node) {
  int anc = tree.parent[node];
  size_t node_base = static_cast<size_t>(node) * tree.total_words;
  size_t anc_base = static_cast<size_t>(anc) * tree.total_words;
  bool changed = false;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    int offset = ds.block_word_offset[b];

    uint64_t any_intersect = 0;
    for (int s = 0; s < blk.n_states; ++s) {
      any_intersect |= (tree.final_[anc_base + offset + s]
                      & tree.prelim[node_base + offset + s]);
    }
    uint64_t has_isect = any_intersect;
    uint64_t no_isect = ~any_intersect & blk.active_mask;

    for (int s = 0; s < blk.n_states; ++s) {
      uint64_t isect = tree.final_[anc_base + offset + s]
                     & tree.prelim[node_base + offset + s];
      uint64_t new_val = (isect & has_isect)
                       | (tree.prelim[node_base + offset + s] & no_isect);
      if (new_val != tree.final_[node_base + offset + s]) changed = true;
      tree.final_[node_base + offset + s] = new_val;
    }
  }

  return changed;
}

// --- Full passes ---

int fitch_downpass(TreeState& tree, const DataSet& ds) {
  int total_steps = 0;

  for (int node : tree.postorder) {
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];

    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      int offset = ds.block_word_offset[b];

      const uint64_t* left_state =
          &tree.prelim[static_cast<size_t>(lc) * tree.total_words + offset];
      const uint64_t* right_state =
          &tree.prelim[static_cast<size_t>(rc) * tree.total_words + offset];
      uint64_t* node_state =
          &tree.prelim[static_cast<size_t>(node) * tree.total_words + offset];

      uint64_t any_intersect = 0;
      for (int s = 0; s < blk.n_states; ++s) {
        uint64_t isect = left_state[s] & right_state[s];
        any_intersect |= isect;
      }

      uint64_t needs_union = ~any_intersect & blk.active_mask;
      total_steps += blk.weight * popcount64(needs_union);

      // Store local cost
      tree.local_cost[static_cast<size_t>(node) * tree.n_blocks + b] =
          needs_union;

      for (int s = 0; s < blk.n_states; ++s) {
        uint64_t isect = left_state[s] & right_state[s];
        uint64_t uni = left_state[s] | right_state[s];
        node_state[s] = (isect & any_intersect) | (uni & needs_union);
      }
    }
  }

  return total_steps;
}

void fitch_uppass(TreeState& tree, const DataSet& ds) {
  // Root: final = prelim
  int root = tree.n_tip;
  size_t root_base = static_cast<size_t>(root) * tree.total_words;
  for (int w = 0; w < tree.total_words; ++w) {
    tree.final_[root_base + w] = tree.prelim[root_base + w];
  }

  // Reverse postorder (root to leaves) — skip root itself
  for (int i = static_cast<int>(tree.postorder.size()) - 1; i >= 0; --i) {
    int node = tree.postorder[i];
    if (node == root) continue;
    uppass_node(tree, ds, node);
  }
}

int fitch_score(TreeState& tree, const DataSet& ds) {
  int score = fitch_downpass(tree, ds);
  fitch_uppass(tree, ds);
  return score;
}

// --- Incremental passes for SPR clipping ---

int fitch_incremental_downpass(TreeState& tree, const DataSet& ds,
                               int start_node) {
  int length_delta = 0;
  int node = start_node;

  while (true) {
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];

    // Save old state before overwriting
    tree.save_node_state(node);

    bool changed = false;

    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      int offset = ds.block_word_offset[b];

      const uint64_t* left_state =
          &tree.prelim[static_cast<size_t>(lc) * tree.total_words + offset];
      const uint64_t* right_state =
          &tree.prelim[static_cast<size_t>(rc) * tree.total_words + offset];
      uint64_t* node_state =
          &tree.prelim[static_cast<size_t>(node) * tree.total_words + offset];

      // Subtract old local cost
      size_t cost_idx = static_cast<size_t>(node) * tree.n_blocks + b;
      length_delta -= blk.weight * popcount64(tree.local_cost[cost_idx]);

      // Compute new prelim
      uint64_t any_intersect = 0;
      for (int s = 0; s < blk.n_states; ++s) {
        uint64_t isect = left_state[s] & right_state[s];
        any_intersect |= isect;
      }

      uint64_t needs_union = ~any_intersect & blk.active_mask;
      length_delta += blk.weight * popcount64(needs_union);

      // Store new local cost
      tree.local_cost[cost_idx] = needs_union;

      for (int s = 0; s < blk.n_states; ++s) {
        uint64_t isect = left_state[s] & right_state[s];
        uint64_t uni = left_state[s] | right_state[s];
        uint64_t new_val = (isect & any_intersect) | (uni & needs_union);
        if (new_val != node_state[s]) changed = true;
        node_state[s] = new_val;
      }
    }

    // Stop if at root or prelim didn't change
    if (!changed || node == tree.n_tip) break;
    int p = tree.parent[node];
    if (p == node) break;
    node = p;
  }

  return length_delta;
}

void fitch_incremental_uppass(TreeState& tree, const DataSet& ds,
                              int clip_ancestor) {
  // Recompute final_ for nodes affected by the clip.
  // Start at root, propagate downward through changed nodes.
  // Only nodes whose ancestor's final changed need recomputation.

  // Root: final = prelim (may have changed during incremental downpass)
  int root = tree.n_tip;
  size_t root_base = static_cast<size_t>(root) * tree.total_words;
  bool root_changed = false;
  for (int w = 0; w < tree.total_words; ++w) {
    if (tree.final_[root_base + w] != tree.prelim[root_base + w]) {
      root_changed = true;
    }
    tree.final_[root_base + w] = tree.prelim[root_base + w];
  }

  if (!root_changed) {
    // Check if any node on the path had its prelim changed but root didn't.
    // The incremental downpass may have stopped before reaching root.
    // In that case, final_ only needs updating for the sister of the clip
    // point (whose new ancestor may have different final states).
  }

  // Use reverse postorder, but only visit nodes whose ancestor's final
  // may have changed. We track this with a "dirty" flag per node.
  std::vector<bool> dirty(tree.n_node, false);

  // Mark root as dirty (we just updated it; its children need checking)
  dirty[root] = root_changed;

  // Also mark the clip ancestor — its children definitely need checking
  // because the topology changed around it.
  if (clip_ancestor >= tree.n_tip) {
    dirty[clip_ancestor] = true;
  }

  // Reverse postorder traversal
  for (int i = static_cast<int>(tree.postorder.size()) - 1; i >= 0; --i) {
    int node = tree.postorder[i];
    if (node == root) continue;

    int anc = tree.parent[node];
    if (!dirty[anc]) continue;

    // Save state before modifying (if not already saved during downpass)
    // The downpass only saved nodes on the rootward path; the uppass may
    // touch additional nodes (siblings and their descendants).
    tree.save_node_state(node);

    bool changed = uppass_node(tree, ds, node);

    if (changed && node >= tree.n_tip) {
      dirty[node] = true;
    }
  }
}

// --- Indirect tree length calculation (Goloboff 1996) ---

int fitch_indirect_length(const uint64_t* clip_prelim,
                          const TreeState& tree,
                          const DataSet& ds,
                          int node_a, int node_d) {
  // Virtual root Y at edge (A, D) is the union of final states:
  //   Y = final(A) | final(D)
  // All states present at either endpoint are available on that edge.
  // Extra steps = count of characters where clip_prelim & Y = 0.
  //
  // This is exact for non-additive characters (Goloboff 1996).
  // Using Fitch intersection-then-union here would narrow Y and overcount.

  int extra_steps = 0;

  size_t a_base = static_cast<size_t>(node_a) * tree.total_words;
  size_t d_base = static_cast<size_t>(node_d) * tree.total_words;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    int offset = ds.block_word_offset[b];

    uint64_t any_hit = 0;
    for (int s = 0; s < blk.n_states; ++s) {
      uint64_t vroot = tree.final_[a_base + offset + s]
                     | tree.final_[d_base + offset + s];
      any_hit |= (clip_prelim[offset + s] & vroot);
    }

    uint64_t needs_step = ~any_hit & blk.active_mask;
    extra_steps += blk.weight * popcount64(needs_step);
  }

  return extra_steps;
}


// --- Per-character step extraction ---

void extract_char_steps(const TreeState& tree, const DataSet& ds,
                        std::vector<int>& char_steps) {
  // Standard blocks: count from local_cost
  for (int node : tree.postorder) {
    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      if (blk.has_inapplicable) continue;
      uint64_t mask =
          tree.local_cost[static_cast<size_t>(node) * tree.n_blocks + b];
      while (mask) {
        int c = ctz64(mask);
        char_steps[blk.pattern_index[c]] += 1;
        mask &= mask - 1;
      }
    }
  }

  // NA blocks: use same formula as Pass 3 of three-pass algorithm.
  // Requires subtree_actives (computed by fitch_na_score).
  for (int node : tree.postorder) {
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    size_t nb = static_cast<size_t>(node) * tree.total_words;
    size_t lb = static_cast<size_t>(lc) * tree.total_words;
    size_t rb = static_cast<size_t>(rc) * tree.total_words;

    for (int b = 0; b < ds.n_blocks; ++b) {
      if (!ds.blocks[b].has_inapplicable) continue;
      const CharBlock& blk = ds.blocks[b];
      int off = ds.block_word_offset[b];
      int k = blk.n_states;

      const uint64_t* L2 = &tree.down2[lb + off];
      const uint64_t* R2 = &tree.down2[rb + off];
      const uint64_t* D2 = &tree.down2[nb + off];

      // Node's applicable second-downpass states
      uint64_t ss_app = 0;
      for (int s = 1; s < k; ++s) ss_app |= D2[s];

      // D2 intersection (including NA state)
      uint64_t any_d2_isect = 0;
      for (int s = 0; s < k; ++s) any_d2_isect |= (L2[s] & R2[s]);

      // subtree_actives: applicable tips anywhere in children's subtrees
      const uint64_t* la = &tree.subtree_actives[lb + off];
      const uint64_t* ra = &tree.subtree_actives[rb + off];
      uint64_t l_act = 0, r_act = 0;
      for (int s = 1; s < k; ++s) { l_act |= la[s]; r_act |= ra[s]; }

      // Same formula as fitch_na_score Pass 3
      uint64_t needs_step = l_act & r_act
                          & ~(ss_app & any_d2_isect) & blk.active_mask;

      while (needs_step) {
        int c = ctz64(needs_step);
        char_steps[blk.pattern_index[c]] += 1;
        needs_step &= needs_step - 1;
      }
    }
  }
}

// --- IW computation ---

double compute_iw(const DataSet& ds, const std::vector<int>& char_steps) {
  double score = 0.0;
  double k = ds.concavity;
  for (int p = 0; p < ds.n_patterns; ++p) {
    int extra = char_steps[p] - ds.min_steps[p];
    if (extra > 0) {
      score += ds.pattern_freq[p] * (static_cast<double>(extra) / (k + extra));
    }
  }
  return score;
}

void precompute_iw_delta(const DataSet& ds,
                         const std::vector<int>& divided_steps,
                         std::vector<double>& iw_delta) {
  double k = ds.concavity;
  for (int p = 0; p < ds.n_patterns; ++p) {
    int e = divided_steps[p] - ds.min_steps[p];
    if (e < 0) e = 0;
    double old_cost = static_cast<double>(e) / (k + e);
    double new_cost = static_cast<double>(e + 1) / (k + e + 1);
    iw_delta[p] = ds.pattern_freq[p] * (new_cost - old_cost);
  }
}

double indirect_iw_length(
    const uint64_t* clip_prelim,
    const TreeState& tree, const DataSet& ds,
    int node_a, int node_d,
    double base_iw,
    const std::vector<double>& iw_delta) {

  double candidate_iw = base_iw;

  size_t a_base = static_cast<size_t>(node_a) * tree.total_words;
  size_t d_base = static_cast<size_t>(node_d) * tree.total_words;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    int offset = ds.block_word_offset[b];

    uint64_t any_hit = 0;
    for (int s = 0; s < blk.n_states; ++s) {
      uint64_t vroot = tree.final_[a_base + offset + s]
                     | tree.final_[d_base + offset + s];
      any_hit |= (clip_prelim[offset + s] & vroot);
    }

    uint64_t needs_step = ~any_hit & blk.active_mask;
    while (needs_step) {
      int c = ctz64(needs_step);
      candidate_iw += iw_delta[blk.pattern_index[c]];
      needs_step &= needs_step - 1;
    }
  }

  return candidate_iw;
}

// --- Unified scoring ---

double score_tree(TreeState& tree, const DataSet& ds) {
  // Check if any block has inapplicable characters
  bool has_na = false;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) { has_na = true; break; }
  }

  if (std::isinf(ds.concavity)) {
    // Equal weights
    if (has_na) {
      return static_cast<double>(fitch_na_score(tree, ds));
    } else {
      return static_cast<double>(fitch_score(tree, ds));
    }
  }

  // Implied weights: run appropriate scoring pass, then extract per-char steps
  if (has_na) {
    fitch_na_score(tree, ds);
  } else {
    fitch_score(tree, ds);
  }

  std::vector<int> char_steps(ds.n_patterns, 0);
  extract_char_steps(tree, ds, char_steps);
  return compute_iw(ds, char_steps);
}


// =========================================================================
// Inapplicable (NA) three-pass scoring (Brazeau et al. 2019)
// =========================================================================

#include "ts_fitch_na.inc"

} // namespace ts
