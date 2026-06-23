// This file is #included at the end of ts_fitch.cpp, after ts_fitch_na_incr.h.
// NA-aware dirty-set incremental rescore for the SPR accept path.
//
// Mirrors the EW dirty-set approach in fitch_dirty_downpass / fitch_dirty_uppass
// but with the NA-aware Pass 1 / Pass 2 logic from fitch_na_incremental_*.
// Visits each affected node exactly once in postorder — avoids the
// double-chain composition that produced the unresolved −3 in the prior
// T-300 attempt on the EW path.
//
// Provides:
//   fitch_na_dirty_downpass()  — NA-aware dirty-set Pass 1
//   fitch_na_dirty_uppass()    — NA-aware dirty-set Pass 2 + tip update
//
// Both functions take two seed nodes (start_a, start_b), mark the rootward
// paths from each, and update prelim/final_/local_cost/subtree_actives/down2
// for nodes on the union of paths.  Off-path nodes retain their state from
// the prior full rescore, which remains valid because their children are
// unchanged by the SPR move.

// =========================================================================
// NA-aware dirty-set Pass 1 (first downpass)
// =========================================================================

int fitch_na_dirty_downpass(TreeState& tree, const DataSet& ds,
                             int start_a, int start_b, int start_c) {
  std::vector<char> dirty(tree.n_node, 0);

  auto mark_path = [&](int node) {
    while (node >= tree.n_tip && !dirty[node]) {
      dirty[node] = 1;
      int p = tree.parent[node];
      if (p == node) break;  // root
      node = p;
    }
  };
  mark_path(start_a);
  mark_path(start_b);
  mark_path(start_c);   // third seed (-1 = no-op): TBR-reroot dirty region

  int length_delta = 0;

  for (int node : tree.postorder) {
    if (node < tree.n_tip) continue;
    if (!dirty[node]) continue;

    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];

    tree.save_node_state(node);

    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      if (blk.active_mask == 0) continue;
      int off = ds.block_word_offset[b];
      int k = blk.n_states;

      size_t nb = static_cast<size_t>(node) * tree.total_words;
      size_t lb = static_cast<size_t>(lc) * tree.total_words;
      size_t rb = static_cast<size_t>(rc) * tree.total_words;

      if (!blk.has_inapplicable) {
        // Standard Fitch — mirrors fitch_dirty_downpass.
        const uint64_t* left_state = &tree.prelim[lb + off];
        const uint64_t* right_state = &tree.prelim[rb + off];
        uint64_t* node_state = &tree.prelim[nb + off];

        size_t cost_idx = static_cast<size_t>(node) * tree.n_blocks + b;
        uint64_t old_cost = tree.local_cost[cost_idx];
        int old_nu = popcount64(old_cost);
        if (blk.upweight_mask) old_nu += popcount64(old_cost & blk.upweight_mask);
        length_delta -= blk.weight * old_nu;

        uint64_t any_intersect = simd::any_hit_reduce(
            left_state, right_state, k);
        uint64_t needs_union = ~any_intersect & blk.active_mask;
        int new_nu = popcount64(needs_union);
        if (blk.upweight_mask) new_nu += popcount64(needs_union & blk.upweight_mask);
        length_delta += blk.weight * new_nu;

        tree.local_cost[cost_idx] = needs_union;

        simd::fitch_combine(left_state, right_state, node_state,
                            k, any_intersect, needs_union);
      } else {
        // NA-aware first downpass (same logic as fitch_na_incremental_downpass).
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
          N[s] = (isect & case_keep) | ((L[s] | R[s]) & ~case_keep);
        }
        N[0] = (I0 & case_keep)
             | ((L[0] | R[0]) & ~case_keep & ~case_strip);

        // Subtree actives
        const uint64_t* la = &tree.subtree_actives[lb + off];
        const uint64_t* ra = &tree.subtree_actives[rb + off];
        uint64_t* na = &tree.subtree_actives[nb + off];
        na[0] = 0;
        for (int s = 1; s < k; ++s) {
          na[s] = la[s] | ra[s];
        }
      }
    }
  }

  return length_delta;
}


// =========================================================================
// NA-aware dirty-set Pass 2 (first uppass) + tip update
// =========================================================================

void fitch_na_dirty_uppass(TreeState& tree, const DataSet& ds,
                            int start_a, int start_b, int start_c) {
  // Step 1: root final_ = root prelim.  Root is on every rootward path so
  // its prelim may have changed in the downpass.
  int root = tree.n_tip;
  size_t root_base = static_cast<size_t>(root) * tree.total_words;
  for (int w = 0; w < tree.total_words; ++w) {
    tree.final_[root_base + w] = tree.prelim[root_base + w];
  }

  // Step 2: mark dirty_up = same rootward paths as the downpass.  Every
  // node on these paths had its prelim refreshed, so its final_ may shift
  // and its children must be re-checked.
  std::vector<char> dirty_up(tree.n_node, 0);
  auto mark_path = [&](int node) {
    while (node >= tree.n_tip && !dirty_up[node]) {
      dirty_up[node] = 1;
      int p = tree.parent[node];
      if (p == node) break;
      node = p;
    }
  };
  mark_path(start_a);
  mark_path(start_b);
  mark_path(start_c);   // third seed (-1 = no-op): TBR-reroot dirty region

  // Step 3: reverse postorder over internal nodes.  Visit any node whose
  // parent is dirty.  If that node's final_ changes, propagate the flag.
  for (int i = static_cast<int>(tree.postorder.size()) - 1; i >= 0; --i) {
    int node = tree.postorder[i];
    if (node == root) continue;
    if (node < tree.n_tip) continue;  // tips handled below

    int anc = tree.parent[node];
    if (!dirty_up[anc]) continue;

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
        // Standard uppass
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
        // NA-aware first uppass (same logic as fitch_na_incremental_uppass)
        const uint64_t* Np = &tree.prelim[nb + off];
        const uint64_t* A = &tree.final_[ab + off];
        uint64_t* F = &tree.final_[nb + off];

        uint64_t npre_has_NA = Np[0];
        uint64_t npre_has_app = 0;
        for (int s = 1; s < k; ++s) npre_has_app |= Np[s];
        uint64_t anc_app = 0;
        for (int s = 1; s < k; ++s) anc_app |= A[s];
        uint64_t anc_is_NA = A[0] & ~anc_app;

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
          uint64_t child_union = tree.prelim[cl + s] | tree.prelim[cr + s];
          uint64_t new_val = (Np[s] & (case_pass | case_strip))
                           | (child_union & case_children);
          if (new_val != F[s]) any_changed = true;
          F[s] = new_val;
        }
      }
    }

    if (any_changed) {
      dirty_up[node] = 1;
    }
  }

  // Step 4: process tips whose ancestor is dirty_up.  Mirrors the tip loop
  // in fitch_na_incremental_uppass; updates tip final_, down2, and
  // subtree_actives.
  for (int tip = 0; tip < tree.n_tip; ++tip) {
    int anc = tree.parent[tip];
    if (!dirty_up[anc]) continue;

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
