// This file is #included at the end of ts_fitch.cpp, after ts_fitch_na.inc.
// Incremental NA-aware scoring for SPR/TBR candidate evaluation.
//
// Provides:
//   fitch_na_incremental_downpass()  — NA-aware incremental Pass 1
//   fitch_na_incremental_uppass()    — NA-aware incremental Pass 2 + tips
//   fitch_na_pass3_score()           — Full Pass 3 on divided tree
//   fitch_na_indirect_length()       — NA-aware indirect length
//   indirect_na_iw_length()          — NA-aware indirect IW length

// =========================================================================
// Step 1: NA-aware incremental first downpass
// =========================================================================
//
// Walks rootward from start_node, recomputing prelim using NA-aware logic
// for inapplicable blocks and standard Fitch for standard blocks.
// Also maintains subtree_actives for NA blocks.
//
// Returns the length_delta for standard blocks (same accounting as
// fitch_incremental_downpass). NA block steps require Pass 3 for exact count.

int fitch_na_incremental_downpass(TreeState& tree, const DataSet& ds,
                                   int start_node) {
  int length_delta = 0;
  int node = start_node;

  while (true) {
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];

    tree.save_node_state(node);

    bool changed = false;

    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      if (blk.active_mask == 0) continue;
      int off = ds.block_word_offset[b];
      int k = blk.n_states;

      size_t nb = static_cast<size_t>(node) * tree.total_words;
      size_t lb = static_cast<size_t>(lc) * tree.total_words;
      size_t rb = static_cast<size_t>(rc) * tree.total_words;

      if (!blk.has_inapplicable) {
        // Standard Fitch — same as fitch_incremental_downpass
        const uint64_t* left_state = &tree.prelim[lb + off];
        const uint64_t* right_state = &tree.prelim[rb + off];
        uint64_t* node_state = &tree.prelim[nb + off];

        size_t cost_idx = static_cast<size_t>(node) * tree.n_blocks + b;
        uint64_t old_cost = tree.local_cost[cost_idx];
        int old_nu = popcount64(old_cost);
        if (blk.upweight_mask) old_nu += popcount64(old_cost & blk.upweight_mask);
        length_delta -= blk.weight * old_nu;

        uint64_t any_isect = 0;
        for (int s = 0; s < k; ++s) {
          any_isect |= (left_state[s] & right_state[s]);
        }
        uint64_t needs_union = ~any_isect & blk.active_mask;
        int new_nu = popcount64(needs_union);
        if (blk.upweight_mask) new_nu += popcount64(needs_union & blk.upweight_mask);
        length_delta += blk.weight * new_nu;

        tree.local_cost[cost_idx] = needs_union;

        for (int s = 0; s < k; ++s) {
          uint64_t isect = left_state[s] & right_state[s];
          uint64_t uni = left_state[s] | right_state[s];
          uint64_t new_val = (isect & any_isect) | (uni & needs_union);
          if (new_val != node_state[s]) changed = true;
          node_state[s] = new_val;
        }
      } else {
        // NA-aware first downpass (same logic as fitch_na_score Pass 1)
        const uint64_t* L = &tree.prelim[lb + off];
        const uint64_t* R = &tree.prelim[rb + off];
        uint64_t* N = &tree.prelim[nb + off];

        uint64_t I_app = 0, L_app = 0, R_app = 0;
        for (int s = 1; s < k; ++s) {
          I_app |= (L[s] & R[s]);
          L_app |= L[s];
          R_app |= R[s];
        }
        uint64_t I0 = L[0] & R[0];
        uint64_t both_app = L_app & R_app;
        uint64_t case_keep = I_app | (I0 & ~I_app & ~both_app);
        uint64_t case_strip = ~I0 & ~I_app & both_app;

        for (int s = 1; s < k; ++s) {
          uint64_t isect = L[s] & R[s];
          uint64_t new_val = (isect & case_keep) | ((L[s] | R[s]) & ~case_keep);
          if (new_val != N[s]) changed = true;
          N[s] = new_val;
        }
        uint64_t new_n0 = (I0 & case_keep)
                        | ((L[0] | R[0]) & ~case_keep & ~case_strip);
        if (new_n0 != N[0]) changed = true;
        N[0] = new_n0;

        // Subtree actives
        const uint64_t* la = &tree.subtree_actives[lb + off];
        const uint64_t* ra = &tree.subtree_actives[rb + off];
        uint64_t* na = &tree.subtree_actives[nb + off];
        na[0] = 0;
        for (int s = 1; s < k; ++s) {
          uint64_t new_sa = la[s] | ra[s];
          if (new_sa != na[s]) changed = true;
          na[s] = new_sa;
        }
      }
    }

    if (!changed || node == tree.n_tip) break;
    int p = tree.parent[node];
    if (p == node) break;
    node = p;
  }

  return length_delta;
}


// =========================================================================
// Step 2: NA-aware incremental first uppass + tip processing
// =========================================================================

void fitch_na_incremental_uppass(TreeState& tree, const DataSet& ds,
                                  int clip_ancestor) {
  int root = tree.n_tip;
  size_t root_base = static_cast<size_t>(root) * tree.total_words;
  bool root_changed = false;
  for (int w = 0; w < tree.total_words; ++w) {
    if (tree.final_[root_base + w] != tree.prelim[root_base + w]) {
      root_changed = true;
    }
    tree.final_[root_base + w] = tree.prelim[root_base + w];
  }

  std::vector<bool> dirty(tree.n_node, false);
  dirty[root] = root_changed;
  if (clip_ancestor >= tree.n_tip) {
    dirty[clip_ancestor] = true;
  }

  // Reverse postorder: internal nodes
  for (int i = static_cast<int>(tree.postorder.size()) - 1; i >= 0; --i) {
    int node = tree.postorder[i];
    if (node == root) continue;

    int anc = tree.parent[node];
    if (!dirty[anc]) continue;

    tree.save_node_state(node);

    bool any_changed = false;
    size_t nb = static_cast<size_t>(node) * tree.total_words;
    size_t ab = static_cast<size_t>(anc) * tree.total_words;

    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      if (blk.active_mask == 0) continue;
      int off = ds.block_word_offset[b];
      int k = blk.n_states;

      if (!blk.has_inapplicable) {
        // Standard uppass (same as uppass_node)
        uint64_t any_isect = 0;
        for (int s = 0; s < k; ++s) {
          any_isect |= (tree.final_[ab + off + s]
                      & tree.prelim[nb + off + s]);
        }
        uint64_t no_isect = ~any_isect & blk.active_mask;
        for (int s = 0; s < k; ++s) {
          uint64_t isect = tree.final_[ab + off + s]
                         & tree.prelim[nb + off + s];
          uint64_t new_val = (isect & any_isect)
                           | (tree.prelim[nb + off + s] & no_isect);
          if (new_val != tree.final_[nb + off + s]) any_changed = true;
          tree.final_[nb + off + s] = new_val;
        }
      } else {
        // NA-aware first uppass (same logic as fitch_na_score Pass 2)
        const uint64_t* Np = &tree.prelim[nb + off];
        const uint64_t* A = &tree.final_[ab + off];
        uint64_t* F = &tree.final_[nb + off];

        uint64_t npre_has_NA = Np[0];
        uint64_t npre_has_app = 0;
        for (int s = 1; s < k; ++s) npre_has_app |= Np[s];
        uint64_t anc_app = 0;
        for (int s = 1; s < k; ++s) anc_app |= A[s];
        uint64_t anc_is_NA = A[0] & ~anc_app;

        // node is internal (postorder only has internals)
        int ni2 = node - tree.n_tip;
        size_t cl = static_cast<size_t>(tree.left[ni2]) * tree.total_words + off;
        size_t cr = static_cast<size_t>(tree.right[ni2]) * tree.total_words + off;
        uint64_t children_app = 0;
        for (int s = 1; s < k; ++s) {
          children_app |= (tree.prelim[cl + s] | tree.prelim[cr + s]);
        }

        uint64_t case_pass = ~npre_has_NA & blk.active_mask;
        uint64_t case_strip = npre_has_NA & npre_has_app & ~anc_is_NA;
        uint64_t case_children = npre_has_NA & ~npre_has_app
                               & ~anc_is_NA & children_app;
        uint64_t case_force = blk.active_mask
                            & ~case_pass & ~case_strip & ~case_children;

        uint64_t new_f0 = (Np[0] & case_pass) | case_force;
        if (new_f0 != F[0]) any_changed = true;
        F[0] = new_f0;

        for (int s = 1; s < k; ++s) {
          uint64_t child_union =
              tree.prelim[static_cast<size_t>(tree.left[ni2]) * tree.total_words + off + s]
            | tree.prelim[static_cast<size_t>(tree.right[ni2]) * tree.total_words + off + s];
          uint64_t new_val = (Np[s] & (case_pass | case_strip))
                           | (child_union & case_children);
          if (new_val != F[s]) any_changed = true;
          F[s] = new_val;
        }
      }
    }

    if (any_changed && node >= tree.n_tip) {
      dirty[node] = true;
    }
  }

  // Process tips whose ancestor's final_ changed
  for (int tip = 0; tip < tree.n_tip; ++tip) {
    int anc = tree.parent[tip];
    if (!dirty[anc]) continue;

    tree.save_node_state(tip);

    size_t tb = static_cast<size_t>(tip) * tree.total_words;
    size_t ab = static_cast<size_t>(anc) * tree.total_words;

    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      if (blk.active_mask == 0) continue;
      int off = ds.block_word_offset[b];
      int k = blk.n_states;

      if (!blk.has_inapplicable) {
        // Standard tip uppass
        uint64_t any_isect = 0;
        for (int s = 0; s < k; ++s) {
          any_isect |= (tree.final_[ab + off + s]
                      & tree.prelim[tb + off + s]);
        }
        uint64_t no_isect = ~any_isect & blk.active_mask;
        for (int s = 0; s < k; ++s) {
          uint64_t isect = tree.final_[ab + off + s]
                         & tree.prelim[tb + off + s];
          tree.final_[tb + off + s] = (isect & any_isect)
                                    | (tree.prelim[tb + off + s] & no_isect);
        }
      } else {
        // NA-aware tip update (matches morphy's mpl_fitch_NA_tip_update)
        const uint64_t* T = &tree.prelim[tb + off];
        const uint64_t* A = &tree.final_[ab + off];
        uint64_t* F = &tree.final_[tb + off];

        uint64_t any_isect = 0;
        for (int s = 0; s < k; ++s) any_isect |= (T[s] & A[s]);

        uint64_t anc_app = 0;
        for (int s = 1; s < k; ++s) anc_app |= A[s];

        uint64_t strip_na = any_isect & anc_app;

        F[0] = T[0] & ~strip_na;
        for (int s = 1; s < k; ++s) F[s] = T[s];

        // down2 = final for tips
        uint64_t* D2 = &tree.down2[tb + off];
        for (int s = 0; s < k; ++s) D2[s] = F[s];

        // Update tip subtree_actives
        uint64_t* sa = &tree.subtree_actives[tb + off];
        for (int s = 1; s < k; ++s) {
          sa[s] = (T[s] & A[s] & any_isect) | (T[s] & ~any_isect);
        }
        sa[0] = 0;
      }
    }
  }
}


// =========================================================================
// Step 3: Full Pass 3 (second downpass) on divided tree
// =========================================================================
//
// Runs the full Pass 3 over tree.postorder. Computes down2 for all
// internal nodes and counts steps for both standard and NA blocks.
// Requires that Passes 1+2 (prelim, final_, tip down2, subtree_actives)
// are already current (either from full scoring or incremental passes).
//
// Returns the total EW score of the tree.

int fitch_na_pass3_score(TreeState& tree, const DataSet& ds) {
  int total_steps = 0;

  // Standard block steps from local_cost (set during Pass 1)
  for (int node : tree.postorder) {
    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      if (blk.has_inapplicable || blk.active_mask == 0) continue;
      size_t cost_idx = static_cast<size_t>(node) * tree.n_blocks + b;
      int nu = popcount64(tree.local_cost[cost_idx]);
      if (blk.upweight_mask) nu += popcount64(tree.local_cost[cost_idx] & blk.upweight_mask);
      total_steps += blk.weight * nu;
    }
  }

  // NA blocks: Pass 3 second downpass
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

      const uint64_t* F = &tree.final_[nb + off];
      const uint64_t* L2 = &tree.down2[lb + off];
      const uint64_t* R2 = &tree.down2[rb + off];
      uint64_t* D2 = &tree.down2[nb + off];

      uint64_t ss_app = 0;
      for (int s = 1; s < k; ++s) ss_app |= F[s];

      uint64_t any_isect = 0;
      for (int s = 0; s < k; ++s) any_isect |= (L2[s] & R2[s]);

      uint64_t I_app = 0;
      for (int s = 1; s < k; ++s) I_app |= (L2[s] & R2[s]);

      const uint64_t* la = &tree.subtree_actives[lb + off];
      const uint64_t* ra = &tree.subtree_actives[rb + off];
      uint64_t l_act = 0, r_act = 0;
      for (int s = 1; s < k; ++s) { l_act |= la[s]; r_act |= ra[s]; }

      uint64_t needs_step = l_act & r_act
                          & ~(ss_app & any_isect) & blk.active_mask;
      int ns_p3 = popcount64(needs_step);
      if (blk.upweight_mask) ns_p3 += popcount64(needs_step & blk.upweight_mask);
      total_steps += blk.weight * ns_p3;

      // Compute down2
      uint64_t na_only_isect = any_isect & ~I_app;
      for (int s = 1; s < k; ++s) {
        uint64_t isect = L2[s] & R2[s];
        uint64_t uni = L2[s] | R2[s];
        D2[s] = ss_app & ((isect & any_isect) | (uni & ~any_isect));
      }
      D2[0] = (~ss_app | na_only_isect) & blk.active_mask;

      // Update subtree_actives at this node (same as Pass 3 in fitch_na_score)
      uint64_t* na_out = &tree.subtree_actives[nb + off];
      na_out[0] = 0;
      for (int s = 1; s < k; ++s) na_out[s] = (la[s] | ra[s]);
    }
  }

  return total_steps;
}


// =========================================================================
// Step 4: NA-aware indirect length calculation
// =========================================================================
//
// For standard blocks: identical to fitch_indirect_length.
// For NA blocks: uses subtree_actives to suppress steps where either
// the clip subtree or the edge-below subtree has no applicable tips.
//
// clip_actives: subtree_actives for the clip node (total_words entries).
//   For SPR, points into tree.subtree_actives at clip_node offset.
//   For TBR rerootings, the same pointer works because rerooting doesn't
//   change which tips are in the clip subtree.

int fitch_na_indirect_length(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
    const TreeState& tree,
    const DataSet& ds,
    int node_a, int node_d) {
  int extra_steps = 0;

  size_t a_base = static_cast<size_t>(node_a) * tree.total_words;
  size_t d_base = static_cast<size_t>(node_d) * tree.total_words;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int off = ds.block_word_offset[b];
    int k = blk.n_states;

    if (!blk.has_inapplicable) {
      uint64_t any_hit = simd::any_hit_reduce3(
          &clip_prelim[off],
          &tree.final_[a_base + off],
          &tree.final_[d_base + off], k);
      uint64_t needs_step = ~any_hit & blk.active_mask;
      int ns = popcount64(needs_step);
      if (blk.upweight_mask) ns += popcount64(needs_step & blk.upweight_mask);
      extra_steps += blk.weight * ns;
    } else {
      uint64_t any_hit = simd::any_hit_reduce3_from1(
          &clip_prelim[off],
          &tree.final_[a_base + off],
          &tree.final_[d_base + off], k);

      uint64_t clip_has_active = simd::or_reduce(&clip_actives[off], k, 1);

      const uint64_t* d_act = &tree.subtree_actives[d_base + off];
      uint64_t below_has_active = simd::or_reduce(d_act, k, 1);

      uint64_t needs_step = ~any_hit
                          & clip_has_active
                          & below_has_active
                          & blk.active_mask;
      int ns = popcount64(needs_step);
      if (blk.upweight_mask) ns += popcount64(needs_step & blk.upweight_mask);
      extra_steps += blk.weight * ns;
    }
  }

  return extra_steps;
}


// NA-aware bounded indirect length (early termination)
int fitch_na_indirect_length_bounded(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
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
    int off = ds.block_word_offset[b];
    int k = blk.n_states;
    uint64_t needs_step;

    if (!blk.has_inapplicable) {
      uint64_t any_hit = simd::any_hit_reduce3(
          &clip_prelim[off],
          &tree.final_[a_base + off],
          &tree.final_[d_base + off], k);
      needs_step = ~any_hit & blk.active_mask;
    } else {
      uint64_t any_hit = simd::any_hit_reduce3_from1(
          &clip_prelim[off],
          &tree.final_[a_base + off],
          &tree.final_[d_base + off], k);
      uint64_t clip_has_active = simd::or_reduce(&clip_actives[off], k, 1);
      const uint64_t* d_act = &tree.subtree_actives[d_base + off];
      uint64_t below_has_active = simd::or_reduce(d_act, k, 1);
      needs_step = ~any_hit & clip_has_active & below_has_active
                 & blk.active_mask;
    }
    int ns = popcount64(needs_step);
    if (blk.upweight_mask) ns += popcount64(needs_step & blk.upweight_mask);
    extra_steps += blk.weight * ns;
    if (extra_steps >= cutoff) return extra_steps;
  }

  return extra_steps;
}

// NA-aware cached indirect length (pre-computed vroot + below_actives)
// below_actives: pre-computed per-edge OR of subtree_actives[D] applicable
// words (one uint64_t per edge, only for NA filtering).
int fitch_na_indirect_length_cached(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
    const uint64_t* vroot,
    const uint64_t* below_actives,
    const DataSet& ds,
    int cutoff) {
  int extra_steps = 0;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int off = ds.block_word_offset[b];
    int k = blk.n_states;
    uint64_t needs_step;

    if (!blk.has_inapplicable) {
      uint64_t any_hit = simd::any_hit_reduce(
          &clip_prelim[off], &vroot[off], k);
      needs_step = ~any_hit & blk.active_mask;
    } else {
      uint64_t any_hit = simd::any_hit_reduce_from1(
          &clip_prelim[off], &vroot[off], k);
      uint64_t clip_has_active = simd::or_reduce(&clip_actives[off], k, 1);
      needs_step = ~any_hit & clip_has_active & below_actives[b]
                 & blk.active_mask;
    }
    int ns = popcount64(needs_step);
    if (blk.upweight_mask) ns += popcount64(needs_step & blk.upweight_mask);
    extra_steps += blk.weight * ns;
    if (extra_steps >= cutoff) return extra_steps;
  }

  return extra_steps;
}

// NA-aware bounded indirect IW length
double indirect_na_iw_length_bounded(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
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
    int off = ds.block_word_offset[b];
    int k = blk.n_states;
    uint64_t needs_step;

    if (!blk.has_inapplicable) {
      uint64_t any_hit = simd::any_hit_reduce3(
          &clip_prelim[off],
          &tree.final_[a_base + off],
          &tree.final_[d_base + off], k);
      needs_step = ~any_hit & blk.active_mask;
    } else {
      uint64_t any_hit = simd::any_hit_reduce3_from1(
          &clip_prelim[off],
          &tree.final_[a_base + off],
          &tree.final_[d_base + off], k);
      uint64_t clip_has_active = simd::or_reduce(&clip_actives[off], k, 1);
      const uint64_t* d_act = &tree.subtree_actives[d_base + off];
      uint64_t below_has_active = simd::or_reduce(d_act, k, 1);
      needs_step = ~any_hit & clip_has_active & below_has_active
                 & blk.active_mask;
    }

    while (needs_step) {
      int c = ctz64(needs_step);
      candidate_iw += iw_delta[blk.pattern_index[c]];
      needs_step &= needs_step - 1;
    }
    if (candidate_iw >= cutoff) return candidate_iw;
  }

  return candidate_iw;
}

// NA-aware cached indirect IW length
double indirect_na_iw_length_cached(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
    const uint64_t* vroot,
    const uint64_t* below_actives,
    const DataSet& ds,
    double base_iw,
    const std::vector<double>& iw_delta,
    double cutoff) {

  double candidate_iw = base_iw;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int off = ds.block_word_offset[b];
    int k = blk.n_states;
    uint64_t needs_step;

    if (!blk.has_inapplicable) {
      uint64_t any_hit = simd::any_hit_reduce(
          &clip_prelim[off], &vroot[off], k);
      needs_step = ~any_hit & blk.active_mask;
    } else {
      uint64_t any_hit = simd::any_hit_reduce_from1(
          &clip_prelim[off], &vroot[off], k);
      uint64_t clip_has_active = simd::or_reduce(&clip_actives[off], k, 1);
      needs_step = ~any_hit & clip_has_active & below_actives[b]
                 & blk.active_mask;
    }

    while (needs_step) {
      int c = ctz64(needs_step);
      candidate_iw += iw_delta[blk.pattern_index[c]];
      needs_step &= needs_step - 1;
    }
    if (candidate_iw >= cutoff) return candidate_iw;
  }

  return candidate_iw;
}


// NA-aware indirect IW length
double indirect_na_iw_length(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
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
    int off = ds.block_word_offset[b];
    int k = blk.n_states;

    uint64_t needs_step;

    if (!blk.has_inapplicable) {
      uint64_t any_hit = simd::any_hit_reduce3(
          &clip_prelim[off],
          &tree.final_[a_base + off],
          &tree.final_[d_base + off], k);
      needs_step = ~any_hit & blk.active_mask;
    } else {
      uint64_t any_hit = simd::any_hit_reduce3_from1(
          &clip_prelim[off],
          &tree.final_[a_base + off],
          &tree.final_[d_base + off], k);
      uint64_t clip_has_active = simd::or_reduce(&clip_actives[off], k, 1);
      const uint64_t* d_act = &tree.subtree_actives[d_base + off];
      uint64_t below_has_active = simd::or_reduce(d_act, k, 1);

      needs_step = ~any_hit & clip_has_active & below_has_active
                 & blk.active_mask;
    }

    while (needs_step) {
      int c = ctz64(needs_step);
      candidate_iw += iw_delta[blk.pattern_index[c]];
      needs_step &= needs_step - 1;
    }
  }

  return candidate_iw;
}
