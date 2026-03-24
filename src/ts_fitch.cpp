#include "ts_fitch.h"
#include "ts_hsj.h"
#include "ts_sankoff.h"
#include <vector>
#include <R.h>

namespace ts {

// --- Node-level helpers ---

int fitch_downpass_node(
    const uint64_t* left_state,
    const uint64_t* right_state,
    uint64_t* node_state,
    int n_states,
    uint64_t active_mask)
{
  if (active_mask == 0) return 0;

  // Pass 1: compute any_intersect = OR( left[s] & right[s] )
  uint64_t any_intersect = simd::any_hit_reduce(left_state, right_state,
                                                 n_states);

  uint64_t needs_union = ~any_intersect & active_mask;
  int steps = popcount64(needs_union);

  // Pass 2: compute output states with broadcast masks
#if defined(TS_SIMD_SSE2) || defined(TS_SIMD_NEON)
  simd::v128 ai = simd::set1_64(any_intersect);
  simd::v128 nu = simd::set1_64(needs_union);
  int s = 0;
  for (; s + 2 <= n_states; s += 2) {
    simd::v128 l = simd::loadu128(&left_state[s]);
    simd::v128 r = simd::loadu128(&right_state[s]);
    simd::v128 isect = simd::and128(l, r);
    simd::v128 uni = simd::or128(l, r);
    simd::storeu128(&node_state[s],
        simd::or128(simd::and128(isect, ai), simd::and128(uni, nu)));
  }
  for (; s < n_states; ++s) {
    uint64_t isect = left_state[s] & right_state[s];
    uint64_t uni = left_state[s] | right_state[s];
    node_state[s] = (isect & any_intersect) | (uni & needs_union);
  }
#else
  for (int s = 0; s < n_states; ++s) {
    uint64_t isect = left_state[s] & right_state[s];
    uint64_t uni = left_state[s] | right_state[s];
    node_state[s] = (isect & any_intersect) | (uni & needs_union);
  }
#endif

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
    if (blk.active_mask == 0) continue;
    int offset = ds.block_word_offset[b];

    uint64_t any_intersect = simd::any_hit_reduce(
        &tree.final_[anc_base + offset],
        &tree.prelim[node_base + offset],
        blk.n_states);
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
      if (blk.active_mask == 0) continue;
      int offset = ds.block_word_offset[b];

      const uint64_t* left_state =
          &tree.prelim[static_cast<size_t>(lc) * tree.total_words + offset];
      const uint64_t* right_state =
          &tree.prelim[static_cast<size_t>(rc) * tree.total_words + offset];
      uint64_t* node_state =
          &tree.prelim[static_cast<size_t>(node) * tree.total_words + offset];

      uint64_t any_intersect = simd::any_hit_reduce(
          left_state, right_state, blk.n_states);

      uint64_t needs_union = ~any_intersect & blk.active_mask;
      int nu = popcount64(needs_union);
      if (blk.upweight_mask) nu += popcount64(needs_union & blk.upweight_mask);
      total_steps += blk.weight * nu;

      // Store local cost
      tree.local_cost[static_cast<size_t>(node) * tree.n_blocks + b] =
          needs_union;

      // Compute output states with broadcast masks
#if defined(TS_SIMD_SSE2) || defined(TS_SIMD_NEON)
      {
        simd::v128 ai = simd::set1_64(any_intersect);
        simd::v128 nuvec = simd::set1_64(needs_union);
        int s = 0;
        for (; s + 2 <= blk.n_states; s += 2) {
          simd::v128 l = simd::loadu128(&left_state[s]);
          simd::v128 r = simd::loadu128(&right_state[s]);
          simd::v128 isect = simd::and128(l, r);
          simd::v128 uni = simd::or128(l, r);
          simd::storeu128(&node_state[s],
              simd::or128(simd::and128(isect, ai), simd::and128(uni, nuvec)));
        }
        for (; s < blk.n_states; ++s) {
          uint64_t isect = left_state[s] & right_state[s];
          uint64_t uni = left_state[s] | right_state[s];
          node_state[s] = (isect & any_intersect) | (uni & needs_union);
        }
      }
#else
      for (int s = 0; s < blk.n_states; ++s) {
        uint64_t isect = left_state[s] & right_state[s];
        uint64_t uni = left_state[s] | right_state[s];
        node_state[s] = (isect & any_intersect) | (uni & needs_union);
      }
#endif
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
      if (blk.active_mask == 0) continue;
      int offset = ds.block_word_offset[b];

      const uint64_t* left_state =
          &tree.prelim[static_cast<size_t>(lc) * tree.total_words + offset];
      const uint64_t* right_state =
          &tree.prelim[static_cast<size_t>(rc) * tree.total_words + offset];
      uint64_t* node_state =
          &tree.prelim[static_cast<size_t>(node) * tree.total_words + offset];

      // Subtract old local cost
      size_t cost_idx = static_cast<size_t>(node) * tree.n_blocks + b;
      uint64_t old_cost = tree.local_cost[cost_idx];
      int old_nu = popcount64(old_cost);
      if (blk.upweight_mask) old_nu += popcount64(old_cost & blk.upweight_mask);
      length_delta -= blk.weight * old_nu;

      // Compute new prelim
      uint64_t any_intersect = simd::any_hit_reduce(
          left_state, right_state, blk.n_states);

      uint64_t needs_union = ~any_intersect & blk.active_mask;
      int new_nu = popcount64(needs_union);
      if (blk.upweight_mask) new_nu += popcount64(needs_union & blk.upweight_mask);
      length_delta += blk.weight * new_nu;

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
    if (blk.active_mask == 0) continue;
    int offset = ds.block_word_offset[b];

    uint64_t any_hit = simd::any_hit_reduce3(
        &clip_prelim[offset],
        &tree.final_[a_base + offset],
        &tree.final_[d_base + offset],
        blk.n_states);

    uint64_t needs_step = ~any_hit & blk.active_mask;
    int ns = popcount64(needs_step);
    if (blk.upweight_mask) ns += popcount64(needs_step & blk.upweight_mask);
    extra_steps += blk.weight * ns;
  }

  return extra_steps;
}


// Early-termination variant of fitch_indirect_length.
int fitch_indirect_length_bounded(const uint64_t* clip_prelim,
                                  const TreeState& tree,
                                  const DataSet& ds,
                                  int node_a, int node_d,
                                  int cutoff) {
  int extra_steps = 0;

  size_t a_base = static_cast<size_t>(node_a) * tree.total_words;
  size_t d_base = static_cast<size_t>(node_d) * tree.total_words;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int offset = ds.block_word_offset[b];

    uint64_t any_hit = simd::any_hit_reduce3(
        &clip_prelim[offset],
        &tree.final_[a_base + offset],
        &tree.final_[d_base + offset],
        blk.n_states);

    uint64_t needs_step = ~any_hit & blk.active_mask;
    int ns = popcount64(needs_step);
    if (blk.upweight_mask) ns += popcount64(needs_step & blk.upweight_mask);
    extra_steps += blk.weight * ns;
    if (extra_steps >= cutoff) return extra_steps;
  }

  return extra_steps;
}

// Precomputed-vroot variant with early termination.
int fitch_indirect_length_cached(const uint64_t* clip_prelim,
                                 const uint64_t* vroot,
                                 const DataSet& ds,
                                 int cutoff) {
  int extra_steps = 0;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int offset = ds.block_word_offset[b];

    uint64_t any_hit = simd::any_hit_reduce(
        &clip_prelim[offset], &vroot[offset], blk.n_states);

    uint64_t needs_step = ~any_hit & blk.active_mask;
    int ns = popcount64(needs_step);
    if (blk.upweight_mask) ns += popcount64(needs_step & blk.upweight_mask);
    extra_steps += blk.weight * ns;
    if (extra_steps >= cutoff) return extra_steps;
  }

  return extra_steps;
}


// --- Per-character step extraction ---

void extract_char_steps(const TreeState& tree, const DataSet& ds,
                        std::vector<int>& char_steps) {
  std::fill(char_steps.begin(), char_steps.end(), 0);
  // Standard blocks: count from local_cost
  for (int node : tree.postorder) {
    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      if (blk.has_inapplicable || blk.active_mask == 0) continue;
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
      if (!ds.blocks[b].has_inapplicable || ds.blocks[b].active_mask == 0)
        continue;
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
  for (int p = 0; p < ds.n_patterns; ++p) {
    int extra = char_steps[p] - ds.min_steps[p];
    if (extra > 0) {
      score += ds.pattern_freq[p] * ds.phi[p] *
               (static_cast<double>(extra) / (ds.eff_k[p] + extra));
    }
  }
  return score;
}

void precompute_iw_delta(const DataSet& ds,
                         const std::vector<int>& divided_steps,
                         std::vector<double>& iw_delta) {
  for (int p = 0; p < ds.n_patterns; ++p) {
    int e = divided_steps[p] - ds.min_steps[p];
    if (e < 0) {
      // Below minimum in divided tree: adding one step still leaves at or
      // below minimum, so reconnection incurs no IW cost increase.
      iw_delta[p] = 0.0;
      continue;
    }
    double k = ds.eff_k[p];
    // Guard e == 0: old_cost = 0 for any k >= 0, avoiding 0/0 NaN when k == 0
    double old_cost = (e == 0) ? 0.0
                                : static_cast<double>(e) / (k + e);
    double new_cost = static_cast<double>(e + 1) / (k + e + 1);
    iw_delta[p] = ds.pattern_freq[p] * ds.phi[p] * (new_cost - old_cost);
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
    if (blk.active_mask == 0) continue;
    int offset = ds.block_word_offset[b];

    uint64_t any_hit = simd::any_hit_reduce3(
        &clip_prelim[offset],
        &tree.final_[a_base + offset],
        &tree.final_[d_base + offset],
        blk.n_states);

    uint64_t needs_step = ~any_hit & blk.active_mask;
    while (needs_step) {
      int c = ctz64(needs_step);
      candidate_iw += iw_delta[blk.pattern_index[c]];
      needs_step &= needs_step - 1;
    }
  }

  return candidate_iw;
}

// Early-termination IW variant.
double indirect_iw_length_bounded(
    const uint64_t* clip_prelim,
    const TreeState& tree, const DataSet& ds,
    int node_a, int node_d,
    double base_iw,
    const std::vector<double>& iw_delta,
    double cutoff) {

  double candidate_iw = base_iw;

  size_t a_base = static_cast<size_t>(node_a) * tree.total_words;
  size_t d_base = static_cast<size_t>(node_d) * tree.total_words;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int offset = ds.block_word_offset[b];

    uint64_t any_hit = simd::any_hit_reduce3(
        &clip_prelim[offset],
        &tree.final_[a_base + offset],
        &tree.final_[d_base + offset],
        blk.n_states);

    uint64_t needs_step = ~any_hit & blk.active_mask;
    while (needs_step) {
      int c = ctz64(needs_step);
      candidate_iw += iw_delta[blk.pattern_index[c]];
      needs_step &= needs_step - 1;
    }
    if (candidate_iw >= cutoff) return candidate_iw;
  }

  return candidate_iw;
}

// Precomputed-vroot IW variant with early termination.
double indirect_iw_length_cached(
    const uint64_t* clip_prelim,
    const uint64_t* vroot,
    const DataSet& ds,
    double base_iw,
    const std::vector<double>& iw_delta,
    double cutoff) {

  double candidate_iw = base_iw;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int offset = ds.block_word_offset[b];

    uint64_t any_hit = simd::any_hit_reduce(
        &clip_prelim[offset], &vroot[offset], blk.n_states);

    uint64_t needs_step = ~any_hit & blk.active_mask;
    while (needs_step) {
      int c = ctz64(needs_step);
      candidate_iw += iw_delta[blk.pattern_index[c]];
      needs_step &= needs_step - 1;
    }
    if (candidate_iw >= cutoff) return candidate_iw;
  }

  return candidate_iw;
}

// --- Profile parsimony scoring ---

double compute_profile(const DataSet& ds, const std::vector<int>& char_steps) {
  double score = 0.0;
  for (int p = 0; p < ds.n_patterns; ++p) {
    // Add back precomputed_steps to get the original total step count,
    // since info_amounts is indexed by the full (unsimplified) step count.
    int s = char_steps[p];
    if (!ds.precomputed_steps.empty()) s += ds.precomputed_steps[p];
    // info_amounts is column-major: [(s-1) + info_max_steps * p]
    // s is total step count (1-based); row 0 = 1 total step
    int idx = s - 1;
    if (idx >= 0 && idx < ds.info_max_steps) {
      score += ds.pattern_freq[p] * ds.info_amounts[idx + ds.info_max_steps * p];
    }
    // s == 0: invariant character → 0 cost
    // s > info_max_steps: beyond table → treat as max cost
    else if (s > ds.info_max_steps && ds.info_max_steps > 0) {
      score += ds.pattern_freq[p] *
               ds.info_amounts[(ds.info_max_steps - 1) + ds.info_max_steps * p];
    }
  }
  return score;
}

void precompute_profile_delta(const DataSet& ds,
                               const std::vector<int>& divided_steps,
                               std::vector<double>& delta) {
  for (int p = 0; p < ds.n_patterns; ++p) {
    int s = divided_steps[p];  // reduced steps in divided tree
    // Add back precomputed_steps to get the original total step count,
    // since info_amounts is indexed by the full (unsimplified) step count.
    if (!ds.precomputed_steps.empty()) s += ds.precomputed_steps[p];
    int idx_old = s - 1;       // 0-based row for current step count
    int idx_new = s;           // 0-based row for step count + 1

    double old_cost = (idx_old >= 0 && idx_old < ds.info_max_steps)
        ? ds.info_amounts[idx_old + ds.info_max_steps * p] : 0.0;

    double new_cost;
    if (idx_new >= 0 && idx_new < ds.info_max_steps) {
      new_cost = ds.info_amounts[idx_new + ds.info_max_steps * p];
    } else if (idx_new >= ds.info_max_steps && ds.info_max_steps > 0) {
      // Beyond table: cap at maximum cost (delta = 0 from this point)
      new_cost = ds.info_amounts[(ds.info_max_steps - 1) + ds.info_max_steps * p];
    } else {
      new_cost = 0.0;
    }

    delta[p] = ds.pattern_freq[p] * (new_cost - old_cost);
  }
}

// --- Weighted scoring dispatch (IW or profile) ---

double compute_weighted_score(const DataSet& ds,
                               const std::vector<int>& char_steps) {
  if (ds.scoring_mode == ScoringMode::PROFILE)
    return compute_profile(ds, char_steps);
  return compute_iw(ds, char_steps);
}

void precompute_weighted_delta(const DataSet& ds,
                                const std::vector<int>& divided_steps,
                                std::vector<double>& delta) {
  if (ds.scoring_mode == ScoringMode::PROFILE)
    precompute_profile_delta(ds, divided_steps, delta);
  else
    precompute_iw_delta(ds, divided_steps, delta);
}

// --- Unified scoring ---

double fitch_score_ew(TreeState& tree, const DataSet& ds) {
  // Check if any block has inapplicable characters
  bool has_na = false;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) { has_na = true; break; }
  }

  if (ds.scoring_mode == ScoringMode::EW ||
      ds.scoring_mode == ScoringMode::HSJ) {
    // Equal weights — add back precomputed topology-independent steps
    // (HSJ non-hierarchy chars also use EW scoring)
    if (has_na) {
      return static_cast<double>(fitch_na_score(tree, ds)) + ds.ew_offset;
    } else {
      return static_cast<double>(fitch_score(tree, ds)) + ds.ew_offset;
    }
  }

  // Weighted scoring (IW or profile): run Fitch, extract steps, transform
  if (has_na) {
    fitch_na_score(tree, ds);
  } else {
    fitch_score(tree, ds);
  }

  std::vector<int> char_steps(ds.n_patterns, 0);
  extract_char_steps(tree, ds, char_steps);
  return compute_weighted_score(ds, char_steps);
}

double score_tree(TreeState& tree, const DataSet& ds) {
  if (ds.scoring_mode == ScoringMode::HSJ) {
    // HSJ: Fitch on non-hierarchy chars + HSJ DP on hierarchy blocks.
    // hsj_score() calls fitch_score_ew() internally, avoiding recursion.
    return hsj_score(tree, ds);
  }
  if (ds.scoring_mode == ScoringMode::XFORM) {
    // Xform: Fitch on non-hierarchy chars + Sankoff on recoded hierarchy chars.
    double fitch_part = fitch_score_ew(tree, ds);
    if (ds.sankoff_n_chars > 0) {
      // Build SankoffData on the fly from DataSet fields.
      // Reuse the Sankoff engine's multi-char scoring.
      SankoffData sd;
      sd.n_tips = tree.n_tip;
      sd.n_chars = ds.sankoff_n_chars;
      sd.max_states = ds.sankoff_max_states;
      sd.chars.resize(sd.n_chars);
      for (int ch = 0; ch < sd.n_chars; ++ch) {
        sd.chars[ch].n_states = ds.sankoff_n_states[ch];
        sd.chars[ch].forced_root_state = ds.sankoff_forced_root[ch];
        int ns = ds.sankoff_max_states;
        sd.chars[ch].cost_matrix.resize(ns * ns);
        const double* src = ds.sankoff_cost_matrices.data() +
            ch * ns * ns;
        std::copy(src, src + ns * ns, sd.chars[ch].cost_matrix.begin());
      }
      sd.tip_costs = ds.sankoff_tip_costs;

      fitch_part += sankoff_score(
          tree.left.data(), tree.right.data(),
          tree.postorder.data(), tree.n_internal,
          tree.n_tip, sd);
    }
    return fitch_part;
  }
  return fitch_score_ew(tree, ds);
}


// =========================================================================
// Inapplicable (NA) three-pass scoring (Brazeau et al. 2019)
// =========================================================================

#include "ts_fitch_na.h"

// =========================================================================
// Incremental NA-aware scoring for SPR/TBR candidate evaluation
// =========================================================================

#include "ts_fitch_na_incr.h"

} // namespace ts
